#!/usr/bin/env bash
# vim: set ft=sh ts=2 sw=0 tw=100 et fdm=syntax:
# shellcheck disable=SC2086

#===============================================================================
#        File:  match_control_seqs.sh
# Description:  Match sample's DNA sequences of target spikes using MMseqs2
#        Date:  2019-08-27 11:15:11-0400
#      Author:  Benjamin Leopold (cometsong)
#    Requires:  bash v4+, mmseqs
#  MMses2 ref:  https://github.com/soedinglab/mmseqs2
#===============================================================================
the_script=$(basename ${BASH_SOURCE[0]})

########## initialize default var settings: ##########
declare -i DEBUG # set to empty (DEBUG= ) to disable debug outputs, using '-v' to enable
TmpDir="$TMPDIR"
declare -i KeepDbs=0
#declare -i CreateDbIdx=0  #TODO: add opt to create query db index
declare -i Compressed=1
declare -i Threads=8  #TODO: allow '0' threads = all available processors (/proc/cpuinfo, or /proc/$$/status)

USAGE=\
"Match all sample's sequences of target (spike) seqs using MMseqs2
Usage: ${the_script} [options] <sample.fastq> <target-file>
  <sample-fastq>  File of sample sequences in fast[aq][.gz] format
  <target-file>   target (e.g. spike or control) fast[aq] file
                  or target.db if pre-existing
  -t <dir>  specific dir to hold dbs. Default '$TmpDir'
  -k        keep mmseqs db files
  -v        verbose (log level: DEBUG)
  -c <0|1>  db compression: off(0) or on(1). Default=$Compressed
  -p <int>  set number of processes (threads). Default=$Threads
  -h        show this help message
"

########## sequence file vars
declare -x    \
  target_file \
  sample_file \
  sample_name

########## resulting fields in mmseqs alignment
default_format_output='query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'
format_output="query,target,pident,alnlen,qlen,tlen,qcov,tcov" #N.B. This field order is used in filtering!
format_fields="${format_output//,/	}" #s/,/\t/g


########## Utility Functions ##########
#------- datestamps --------
get_datestamp() { local FORMAT=${*:-"%Y-%m-%dT%H:%M:%S%z"}; date +"$FORMAT"; }
get_date_ymd() { local FORMAT="%Y-%m-%d"; get_datestamp $FORMAT; }
#------- Output util funcs --------
print() { printf '%b\n' "$@"; }
err()   { print "$*" >&2; }
die()   { err "$@" && exit 11; }
usage() { print "$USAGE" && exit; }

#------- Log Output Funcs --------
log()       { err "[$(log.date)] $*"; return 0; }
log.date()  { get_datestamp '%Y%m%d-%H%M%S'; }
log.info()  { log "INFO:  $*"; return 0; }
log.error() { log "ERROR: $*"; return 0; }
log.debug() { [ $DEBUG ] && \
              log "DEBUG: $*"; return 0; }
log.die()   { log.error "$*"; exit 11; }

log_exec( ) {
  local logfile=${1?:-"$(log.date)-${the_script}.log"}
  exec 1> >(tee -a $logfile)
  exec 2>&1
  log.debug "Logging to '$logfile'"
}

# trap: clearing up exec redirects
clearing_up( ) {
  exec 1>&-;
  exec 2>&-;
}
trap clearing_up INT EXIT # HUP QUIT ABRT SEGV 

#=====================================================
# Func: verify arg to opt does not start with undesired char.
# Defaulting to '-' (hyphen)
# Args: <option name> <option arg> [first-char]
check_optarg( ) {
  local opt="$1";
  local arg="$2";
  local chr="${3:--}" # default '-'
  #log.debug "opt: $opt, arg: $arg"
  [ "${arg:0:1}" == "$chr" ] && \
    log.die "Option -'$opt' requires an argument not starting with '$chr'.";
  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}


#################### check cli opts/args ####################
log.debug "Next: Parsing CLI opts"
[ $# -eq 0 ] && usage

OPTIND=1
while getopts :hvkt:p:c: opt; do
  case $opt in
    h) usage; ;;
    v) DEBUG=1;   log.debug "All Debug messages will be shown and logged."; ;;
    k) KeepDbs=1; log.info  "Intermediate mmseqs alignment db files will be kept."; ;;

    t) check_optarg "$opt" "$OPTARG";
       TmpDir="${OPTARG}";
       log.info  "The tmp directory will be '$TmpDir'."; ;;
    p) check_optarg $opt $OPTARG;
       { [[ "$OPTARG" -ge 1 ]] && Threads=$OPTARG; } \
       || log.error "Processes must be >0";
       log.info  "Number of process threads is '$Threads'"; ;;
    c) check_optarg $opt $OPTARG;
       { [[ "$OPTARG" =~ ^[01]$ ]] && Compressed=$OPTARG; } \
       || log.error "Compression must be 0 or 1.";
       log.info  "Compression value is '$Compressed'"; ;;

    \?) log.error "Unknown option: -$OPTARG" $'\n'; usage; ;;
    :)  log.die   "Option -"$OPTARG" requires an argument."; ;;
  esac
done; shift $((OPTIND-1))

[[ $(type -p mmseqs) ]] || \
  log.die "ERROR: Required 'mmseqs' command is not found!"

[[ $# -ge 2 ]] || die "\n$USAGE\nERROR: Missing sample file or target file args."
[ -r $1 ] && sample_file="$1" || log.die "File '$sample_file' is not readable."
[ -r $2 ] && target_file="$2" || log.die "File '$target_file' is not readable."
[ -d $TmpDir ] || mkdir -p $TmpDir || log.die "TmpDir '$TmpDir' is having trouble."


########## Loggin It ##########
sample_base=$(basename ${sample_file})
sample_name=${sample_base%%.*}
logfile="${sample_name}-mmseqs_match.log"
log_exec $logfile


########## Functional! ##########

#=====================================================
# Func: count reads in .fastq or .fasta file; can be .gz
#   NB: presumes fasta = 2 lines and fastq 4 per read
# Args: <seqfile-name>
# Returns: sets num_reada = number of reads in file
count_seq_reads( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  local seqfile="${1}";
  exts="${seqfile#*.}"
  if [ "${exts: -3}" == ".gz" ]; then
    CAT='zcat';
    exts="${exts%.gz}"
  else
    CAT='cat';
  fi
  if   [ "${exts: -1}" == "q" ]; then rdlns=4 # .fast'q'
  elif [ "${exts: -1}" == "a" ]; then rdlns=2 # .fast'a'
  else num_reads=0; # not using fast[qa] filenames?!
  fi
  lines=$($CAT $seqfile | wc -l)
  count=$(($lines / $rdlns))
  log.info "$seqfile has '$count' reads."
  num_reads=$count
}

#=====================================================
# Func: trims file name's extensions, swaps with arg 2
# Args: [name] [new extension]
# Returns: filename with swapped extension
swap_file_ext( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  local name="${1:-"file.ext"}";
  local swap="${2:-".db"}";
  echo "${name%.*}${swap}";
}

#=====================================================
# Func: inserts text at top of the file (on new line)
# Args: <"line string"> <filepath>
# Returns: $?
insert_header( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  line="$1"
  file="$2"
  ex -sc "1i|${line}" -cx $2
  return $?
}

#=====================================================
# Func: filters tsv:
#   - sort into unique "target" sections
#   - calc percent of reads per target to entire sample
# Args: tsv-file sample-name sample-num-reads
# Globals: sets 'pct_reads' =
#     "sample_id target_id pct_matched% num_matched_reads num_sample_reads"
#     for each unique target 
filter_tsv( ) {
  print '' # blank
  log.info "Running '$FUNCNAME' with args: '$*'"
  local tsv="${1}";
  local sample_name="${2}";
  local sample_reads="${3}";
  local preads;
  declare -A targets 
  log.info "Parsing '$tsv' to calc % target matches to sample reads.."
  #while read -r query target pident alnlen qlen tlen qcov tcov; do
  #while read -r ${format_output//,/ }; do # s/,/ /g
  while read -r ${format_fields}; do
    if [[ -n $target ]]; then
      ((targets[$target]+=1))
    fi
  done < $tsv
  for t in ${!targets[@]}; do
    log.debug "target: $t = ${targets[$t]}";
    pct=$( echo "${targets[$t]} / $sample_reads *100" | bc -l )
    preads+="${sample_name}\t${t}\t${pct}%\t${targets[$t]}\t${sample_reads}\n"
  done
  pct_reads="$preads"
}


#=====================================================
# MM "object" => MMseqs2 operations, .init setup env
MM.init( ) {
  log.info "Initialising 'MMseqs2' environment."
  log.info "Running '$FUNCNAME'"

  methods=(
    'db_name'
    'aln_db_name'
    'createdb'
    'createindex'
    'search'
    'convert_tsv'
  ) # summary list
  requires=(
    'swap_file_ext'
    'log() funcs (log.info log.error ...)'
  ) # external funcs/vars  used in MM.* 'methods'

  # init processing vars
  MM_comp="--compressed $Compressed "
  export MM_comp
  MM_opts="$MM_comp --threads $Threads "
  export MM_opts
  MM_createdb_opts="$MM_comp --dont-split-seq-by-len 1 "
  export MM_createdb_opts

  #TODO: add cli opt for min identity
  min_ident="0.97"  # default 97%
  #TODO: add cli opt for fraction alignment/query&target coverage
  cov_fract="0.9"   # default 90% of align-length to query length
  cov_mode=0        # default 0: coverage of query and target
  MM_search_opts="  --min-seq-id $min_ident"
  MM_search_opts+=" --cov-mode $cov_mode"
  MM_search_opts+=" -c $cov_fract"
  MM_search_opts+=" --comp-bias-corr 0"
  MM_search_opts+=" --search-type 3"
  export MM_search_opts

  ## init results vars
  export        \
    sample_name \
    sample_db   \
    target_db   \
    align_db    \
    align_tsv   \
    num_reads   \
    pct_reads

  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}

#=====================================================
# Func: swaps filename extension for '.db'
MM.db_name( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  swap_file_ext "$1"; 
  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}
#=====================================================
# Func: swaps filename extension for '.align.db'
MM.aln_db_name( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  swap_file_ext "$1" ".align.db"; 
  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}

#=====================================================
# Func: creates mmseqs DB from passed 'seqfile'
#         using 'MM_createdb_opts'
# Args: <seq file path> [seq db path]
# Returns: $?
MM.createdb( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  local seqfile="${1}"
  [[ -r "$seqfile" ]] || log.die '"createdb" needs existing sequence file to make the db from.'
  local seq_db=${2:-$(MM.db_name $seqfile)}
  log.info "Creating db for '$seqfile' as '$seq_db'"
  mmseqs createdb $seqfile $seq_db $MM_createdb_opts
  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}

#=====================================================
# Func: creates index on DB for faster uses.
#         using 'MM_opts'
# Args: <seq db path>
# Returns: $?
MM.createindex( ) {
  log.info "Running '$FUNCNAME' with args: '$*'"
  local seq_db=${1:-$(MM.db_name $seqfile)}
  [[ -r "$seq_db" ]] || log.die '"createdb" needs existing sequence db to make the index.'
  log.info "Creating index on '$seq_db'"
  mmseqs createindex $seq_db $TmpDir --search-type 3 $MM_opts
  local state=$?
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}

#=====================================================
# Func: use 'mmseqs' search module on passed dbs
# requires 2 DB args, possible tmpdir, then 4th+ args are for search options
# Globals: sets 'align_db' = query_db + .align.db
MM.search( ) {
  print '' # blank
  log.info "Running '$FUNCNAME' with args: '$*'"
  local qdb="${1}"; shift;  #FIXME: check called with required args
  local tdb="${1}"; shift;
  local tmp="${1:-$TmpDir}"; shift;
  local opts="${*:-$MM_search_opts}"
  log.debug "Search opts: '$opts'"
  local aln="$( MM.aln_db_name $qdb )"
  log.info "Searching '$qdb' for matches in '$tdb'."
  mmseqs search $qdb $tdb $aln $tmp $MM_opts $opts
  local state=$?
  [ $state -eq 0 ] \
    && log.info "Created '$aln'" \
    || log.info "Problems creating '$aln'"
  align_db="$aln"
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}

#=====================================================
# Func: uses 'mmseqs' convertalis module on passed dbs
# Args: <sample_db> <target_db> <align_db>
# Globals: sets 'align_tsv' = sample_name.matches.tsv
MM.convert_tsv( ) {
  print '' # blank
  log.info "Running '$FUNCNAME' with args: '$*'"
  local qdb="${1}"; shift;  #FIXME: check called with required args
  local tdb="${1}"; shift;
  local aln="${1}"; shift;
  local tsv="${sample_name}.matches.tsv"
  local format_out="--format-output ${format_output}"
  log.info "Converting '$aln' to tsv."
  #log.debug "Running: mmseqs convertalis $qdb $tdb $aln $tsv $MM_opts $format_out"
  mmseqs convertalis $qdb $tdb $aln $tsv $MM_opts $format_out
  local state=$?
  #insert_header "$format_fields" "$tsv" # not added to file bc extra row
  log.info "TSV file: $tsv"
  align_tsv="$tsv"
  log.debug "'$FUNCNAME' exit status: '$state'"
  return $state
}


#################### Take action! ####################
declare -i continue=1 # set to 0 when error
log.info "Preparing to match specific sequences to sample file.";
MM.init # base settings

########## sample query setup ##########
log.debug "Next: Sample query setup"
#sample_path=$(dirname $sample_file)
sample_db_name="$( basename $(MM.db_name ${sample_file}) )"
sample_db="${TmpDir}/${sample_db_name}" ## Create sample db in TmpDir
log.debug "using sample db: $sample_db"
[ -r "$sample_db" ] || (MM.createdb "$sample_file" "$sample_db") || log.error "Problems with sample db '$sample_db'." && continue=0
#sample_db="$(MM.db_name ${sample_file})" ## Create sample db in sample's dir
#{ [ -r $sample_db ] || (MM.createdb "$sample_file") } || log.die "Problems with sample db '$sample_db'."

########## target db setup ##########
log.debug "Next: Target db setup"
#target_path=$(dirname target_file)
if [[ "${target_file}" =~ .*db ]] ;then
  log.info "target file is a pre-existing database"
  target_db="${target_file}"
else
  target_db_name="$(basename $(MM.db_name ${target_file}))"
  target_db="${TmpDir}/${target_db_name}" ## Create target db in TmpDir
fi
## Create target db if not pre-existing.
log.debug "using target_db: $target_db"
[ -r "$target_db" ] || (MM.createdb "$target_file" "$target_db") || log.error "Problems with target db '$target_db'." && continue=0
target_idx="${target_db}.idx"
[ -r "$target_idx" ] || (MM.createindex "$target_db") || log.error "Problems with target db index '$target_idx'." && continue=1 # still continue; idx optional

########## search me ##########
[ $continue -eq 1 ] && MM.search "$sample_db" "$target_db" "$TmpDir"        || log.error "Problems with mmseqs2.search db!" && continue=0 # sets 'align_db'
[ $continue -eq 1 ] && MM.convert_tsv "$sample_db" "$target_db" "$align_db" || log.error "Problems with mmseqs2.aligning!"  && continue=0 # sets 'align_tsv'
[ $continue -eq 1 ] && count_seq_reads "$sample_file"                       || log.error "Problems counting sample reads!"  && continue=0 # sets 'num_reads'
[ $continue -eq 1 ] && filter_tsv "$align_tsv" "$sample_name" "$num_reads"  || log.error "Problems percent aligned reads!"  && continue=0 # sets 'pct_reads'
[ -f "${align_tsv}" ] && printf "%b" "$pct_reads" >${align_tsv%.*}.pcts.tsv || log.error "Problems writing pct reads to tsv outfile!"   # writes pct tsv

########## write me ##########
print '' # blank
log.info "Sample: '${sample_name}' has the following percent spike reads:"
[ "$pct_reads" ] && print "$pct_reads" || log.error "Problems getting pct reads!"

########## remove me? ##########
[ "$KeepDbs" -eq 0 ] && \
  log.info "Removing alignment db's" && \
  (rm -v ${align_db}.db* || log.error "Problems removing alignment dbs!")

########## finish me ##########
log.info "'MMseqs2' sequence matching completed for sample '$sample_name'."
