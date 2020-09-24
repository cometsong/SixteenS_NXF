#!/usr/bin/env nextflow

import Helpers
import Logos

logo = new Logo()
println logo.show(logo.logo_mbcore1)

////////// PARAM Defaults \\\\\\\\\\
//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help = false

    params.read_dir = null //NOTE: must be passed by user!
    params.outdir = "${workflow.manifest.name}_results"

    // Email check:
    if (!params.email) { exit 1, "Must supply params.email address to send pipeline report."}

    // Configurable variable parameters
    params.flashMinOverlap = 30
    params.flashMaxOverlap = 150
    params.flashMaxMismatch = 0.1
    params.flashNumThreads = 1
    params.flashCompress = false

    params.trimMinLength = 100
    params.trimLeading = 20
    params.trimTrailing = 17
    params.trimSingleEnd = false

    params.classifyCutoff = 0.8

    //TODO: set default GOLD db to local or on container??
    params.uchimeDb = params.teamDbs + "/16S/rRNA16S.gold.foruchime.fasta"
    if ( file(params.uchimeDb).exists() ) {
        uchimeDb = params.uchimeDb
    } else {
        exit 1, "Uchime db ($params.uchimeDb) does not exist"
    }
    params.uchimeQuiet = false

    // Dehosting Params:
    params.bmtaggerConf = "${workflow.projectDir}/lib/bmtagger_se.conf"
    /* No other host choices configured. */
    params.host = 'mouse'
    if(params.host == 'human') {
        params.bmtaggerBitmask = "/projects/mbiomecore/dbs/hs38/hs38.fa.bitmask"
        params.bmtaggerSrprism = "/projects/mbiomecore/dbs/hs38/hs38.fa.srprism"
    } else if(params.host == 'mouse') {
        params.bmtaggerBitmask = "/projects/mbiomecore/dbs/mm9/mm9.fa.bitmask"
        params.bmtaggerSrprism = "/projects/mbiomecore/dbs/mm9/mm9.fa.srprism"
    } else {
        params.bmtaggerBitmask = "missing bmtagger database"
        params.bmtaggerSrprism = "missing bmtagger database"
        exit 1, "Dehosting 'host' ($params.host) is not configured."
    }
    // genome: homo sapiens: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
    // genome: mus musculus: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_000001635.8_GRCm38.p6/GCA_000001635.8_GRCm38.p6_genomic.fna.gz

    params.spikedRun = true
    params.spikesThreads = 4
    params.spikesCompressDbs = 1
    // TODO: Depends on `ln -s main-mbiome-dbs dir to scripts basedir`
    //params.spikesTargetDb = (process.executor == 'local') \
                          //? "$projectDir/data/dbs/zymo_spikes.db" \
                          //: "/projects/mbiomecore/dbs/spikes_db/zymo_spikes.db"
    params.spikesTargetDb = "/projects/mbiomecore/dbs/spikes_db/zymo_spikes.db"
    if ( file(params.spikesTargetDb).exists() ) {
        spikesTargetDb = params.spikesTargetDb
    } else {
        exit 1, "Spike Check 'target db' ($params.spikesTargetDb) does not exist"
    }

    params.keep_intermediate = false
}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
      The typical command for running the pipeline is as follows:
        nextflow run ${workflow.manifest.name} -profile sumner,singular --read_dir "raw_fastq_dir"
      OR if using a config file with these parameters:
        nextflow -c your_params.config run ${workflow.manifest.name} -profile sumner
    Mandatory:
      --read_dir    Path to paired fastq data.
    Optional:
      --host        Choice host genome. [ choices: 'human', 'mouse' ]
      --spikedRun   Check for spike-in sequences (e.g. Zymo spikes) [default: true]
      --outdir      The output directory where the results will be saved
                        [default: "${manifest.name}_results"]
      --email       The email address to send the pipeline report.
      -name         Name for the pipeline run. If not specified Nextflow will
                        automatically generate a random mnemonic.
      -profile      Environment config to use.
                        [ choices: standard (local), sumner, helix, winter ]
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
if ( ! file(params.outdir).exists() ) {
    outDir = file(params.outdir)
    outDir.mkdirs()
    if ( ! outDir.exists() ) {
        exit 1, "Parameter ERROR: Missing output directory ($params.outdir) is not found: check if path is correct."
    }
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
if(workflow.containerEngine) {
    summary['Container Engine'] = "$workflow.containerEngine"
}
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus"

// Pipeline Params:
summary['Parameters......'] = ''
if(params.email) {
    summary['. E-mail']     = params.email
    //summary['Notification Email'] = workflow.notification.to  //TODO: HOW TO Show this??
}
summary['. Output dir']     = params.outdir
summary['. Read Dir']       = params.read_dir
summary['. Host']           = params.host

//TODO: summary['. Other Params']           = params.host
summary['. Keep Tmp Files']  = params.keep_intermediate

//~~~~~ Project Code == First 4 Chars of Sample File Name ~~~~~\\
sample_one = file("${params.read_dir}/*_R1*.fastq*", checkIfExists: true)[0]
//if (! sample_one.exists) {
//    exit 1, "Samples in read_dir not found!?"
//}
fqName = sample_one.name
projCode = fqName.split("_")[0]
summary['Project Code'] = projCode //NOTE: presumes standard file nomenclature of raw fastq files!

//~~~~~ Flowcell == Basename of Read_Dir ~~~~~\\
flowcell = file(params.read_dir).name //NOTE: presumes standard folder nomenclature of raw fastq files!
summary['Flowcell'] = flowcell

//~~~~~ Count Number of Samples ~~~~~\\
sample_count = file("${params.read_dir}/*_R{1,2}*.fastq*").size() / 2
summary['Sample Count'] = sample_count

//~~~~~ Get Abspath of 'outdir' ~~~~~\\
abs_outdir = file(params.outdir).toAbsolutePath().toString()
//summary['Outdir abspath'] = abs_outdir

summary['Run Start Time']   = workflow.start
//~~~~~~~~~~~~~~~~~~~~~~ Print Summary Info ~~~~~~
println Summary.show(summary)


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Opening Reads Channel ~~~~~
Channel
    .fromFilePairs("${params.read_dir}/*_R{1,2}*fastq*")
    .set { read_pairs_ch }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Processes  ~~~~~
/** TODO Parse software version numbers for each container **/
// basic workflow tasks {
//        trimmomatic
//        flash
//        uchime4
//        dehost w/ bmtagger
//        RDP classifier
//        Spike seqs check
// } then all sample results to {
//        readcounts(*.QC.log)
//        spike summary
//        //upload readcounts to run_qc web app
//        gzip resulting fast[qa], log, tsv files
// }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Trimmomatic ~~~~~
// Check parameters are numbers
def number_vars = [ 'trimLeading', 'trimTrailing', 'trimMinLength', ]
CheckParams.is_numeric(params, number_vars)
process trimmomatic {
    tag "$sample_id"
    label 'trimms'
    publishDir "${params.outdir}/intermediate/${flowcell}/", enabled: params.keep_intermediate, \
        pattern: "*_trim.fastq"
    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.QC.log"
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move", \
        saveAs: { filename -> filename.endsWith("QC.log") ? null : filename }

    input:
    tuple sample_id, file(reads) from read_pairs_ch

    output:
    tuple sample_id, file("*_trim.fastq") into trimmed_reads_fq
    file("*raw.QC.log") into QC_log_raw
    file("*trim.QC.log") into QC_log_trim
    file("*.log") optional true

    script:
    log.info "Running trimmomatic PE,SE on sample_id: ${sample_id}"
    """
    trimmomatic \
    PE -phred33 \
    ${reads[0]} ${reads[1]} \
    ${sample_id}_R1_trim_pe.fastq ${sample_id}_R1_trim_unpaired.fastq \
    ${sample_id}_R2_trim_pe.fastq ${sample_id}_R2_trim_unpaired.fastq \
    MINLEN:${params.trimMinLength}

    trimmomatic \
    SE -phred33 \
    -trimlog ${sample_id}.trim.log \
    ${sample_id}_R1_trim_pe.fastq \
    ${sample_id}_R1_trim.fastq \
    HEADCROP:${params.trimLeading}

    trimmomatic \
    SE -phred33 \
    -trimlog ${sample_id}.trim.log \
    ${sample_id}_R2_trim_pe.fastq \
    ${sample_id}_R2_trim.fastq \
    HEADCROP:${params.trimTrailing}

    QC="QC_raw ${sample_id} \$(echo \$(zcat ${reads[0]}|wc -l)/4|bc)"
    echo "\${QC}" |tee -a ${sample_id}.1.raw.QC.log;
    QC="QC_trim ${sample_id} \$(echo \$(cat ${sample_id}_R1_trim.fastq|wc -l)/4|bc)"
    echo "\${QC}" |tee -a ${sample_id}.2.trim.QC.log;
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Flash Assembler ~~~~~
// Check parameters are numbers
number_vars = [ 'flashMinOverlap', 'flashMaxMismatch', 'flashMaxOverlap', 'flashNumThreads', ]
CheckParams.is_numeric(params, number_vars)
compress = params.flashCompress ? "--compress" : ""
process flash_assembler {
    tag "$sample_id"
    label 'flash'
    publishDir "${params.outdir}/intermediate/${flowcell}/", enabled: params.keep_intermediate, \
        pattern: "*.extendedFrags.fastq"
    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.QC.log"
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move", \
        saveAs: { filename -> filename.endsWith("QC.log") ? null : filename }

    cpus params.flashNumThreads

    input:
    tuple sample_id, file(fastq_pair) from trimmed_reads_fq

    //FIXME: !!!! outputs into ( flash_fq, asmbld_fq ) are somehow getting switched amond sample_id's ??
    output:
    tuple sample_id, file("*.extendedFrags.fastq") into ( flash_fq, asmbld_fq )
    file("*flash.QC.log") into QC_log_flash
    file("*.log") optional true

    script:
    log.info "Running flash assembler on sample_id: ${sample_id}"
    flashfq="${sample_id}.extendedFrags.fastq"
    """
    flash --min-overlap=${params.flashMinOverlap} \
          --max-overlap=${params.flashMaxOverlap} \
          --max-mismatch-density=${params.flashMaxMismatch} \
          --threads=${params.flashNumThreads} \
          --output-prefix=${sample_id} \
          ${compress} \
          ${fastq_pair[0]} ${fastq_pair[1]} \
          > ${sample_id}.flash.log 2>&1

    QC="QC_combined ${sample_id} \$(echo \$(cat ${flashfq}|wc -l)/4|bc)"
    echo "\${QC}"|tee -a ${sample_id}.3.flash.QC.log;
    """
}
process flash_fq2fa {
    tag "$sample_id"
    label 'fq2fa'

    input:
    tuple sample_id, file(in_fastq) from flash_fq

    output:
    tuple sample_id, file("*.extendedFrags.fasta") into flash_fasta_ch

    script:
    log.info "Running flash assembled fastq to fasta on sample_id: ${sample_id}"
    fasta = in_fastq.baseName + ".fasta"

    """
    fastq_to_fasta -Q33 -n \
        -i ${in_fastq} \
        -o ${fasta};
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Dechimera ~~~~~
if ( !file(params.uchimeDb).exists() ){ //FIXME: check exists uchimedb only once... here or up top
    log.error "uchimeDb file was not found. Provided value: ${params.uchimeDb}"
    exit 1
}
quiet = params.uchimeQuiet ? '--quiet' : ''
process dechimera {
    tag "$sample_id"
    label 'uchime'
    publishDir "${params.outdir}/intermediate/${flowcell}/", enabled: params.keep_intermediate, \
        pattern: "*.chimera"
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move"

    scratch true

    input:
    tuple sample_id, file(fasta_in) from flash_fasta_ch
    val db from Channel.value(params.uchimeDb)
    val quiet from quiet

    output:
    tuple sample_id, file("*.chimera") into chimera_out
    file("*.log") optional true

    //TODO: ¿implement plethora of other uchime options?
    script:
    log.info "Running dechimera on sample_id: ${sample_id}"
    """
    uchime \
      --input ${fasta_in} \
      --db ${db} \
      ${quiet} \
      --uchimealns ${sample_id}.alns \
      --uchimeout ${sample_id}.chimera \
      --log ${sample_id}.dechimera.log
    """
}
extractreads_byname_fq = "${workflow.projectDir}/bin/extractreads_byname_fq.pl"
process dechimed_2fqfa {
    tag "$sample_id"
    label 'fq2fa'
    publishDir "${params.outdir}/intermediate/${flowcell}/", enabled: params.keep_intermediate, \
        pattern: "*.clean.fast*"

    input:
    tuple sample_id, file(fastq_in) from asmbld_fq
    tuple sample_id, file(chimes)   from chimera_out

    output:
    tuple sample_id, file("${fastq_out}") into ( dechime_fq, dechime_count_qc )
    tuple sample_id, file("${fasta_out}") into dechime_fa

    shell:
    fastq_out = "${sample_id}.clean.fastq"
    fasta_out = "${sample_id}.clean.fasta"
    chime_ids = "${chimes}.seqids.list"

    log.info "Running dechimed_to_fastq_to_fasta on sample_id: ${sample_id}"
    '''
    awk -F'\t' '{if($17=="Y") print $2}' !{chimes} > !{chime_ids};
    { cat !{fastq_in} | extractreads_byname_fq.pl -v !{chime_ids} > !{fastq_out}; } \
    || echo "ERROR: problems script to extract listed sequences from fastq file!!";

    fastq_to_fasta -Q33 -n -i !{fastq_out} -o !{fasta_out};
    '''
}
process dechimed_qc {
    tag "$sample_id"
    label 'utils'

    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.QC.log"

    input:
    tuple sample_id, file(in_fq) from dechime_count_qc
    output:
    file("*uchime.QC.log") into QC_log_uchime

    script:
    """
    QC="QC_nonchimera ${sample_id} \$(echo \$(cat ${in_fq}|wc -l)/4|bc)"
    echo "\${QC}"|tee -a ${sample_id}.4.uchime.QC.log;
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Dehosting ~~~~~
process dehost {
    tag "$sample_id"
    label 'bmtagger'
    publishDir "${params.outdir}/dehost/${flowcell}/", pattern: "*.host.txt"
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move"

    scratch true

    input:
    tuple sample_id, file(fasta_in) from dechime_fa
    val bitmask from Channel.value(params.bmtaggerBitmask)
    val srprism from Channel.value(params.bmtaggerSrprism)

    output:
    tuple sample_id, file("${sample_id}.host.txt") into dehost_seqs
    file("*.log") optional true

    script:
    log.info "Running dehost on sample_id: ${sample_id}"
    """
    bmtagger_se.sh \
        -b ${bitmask} \
        -x ${srprism} \
        -q 0 \
        -1 ${fasta_in} \
        -C ${params.bmtaggerConf} \
        -o ${sample_id}.host.txt
    """
}
process dehost_2fq2fa {
    tag "$sample_id"
    label 'fq2fa'
    publishDir "${params.outdir}/clean/${flowcell}/", pattern: "*.dehost.fast[qa]", mode: "copy", overwrite: true

    input:
    tuple sample_id, file(fastq_in) from dechime_fq
    tuple sample_id, file(hosts)    from dehost_seqs

    output:
    tuple sample_id, file(fastq_out) into ( cleanfq_out, cleanfq_count_qc )
    tuple sample_id, file(fasta_out) into cleanfa_out

    script:
    fastq_out = fastq_in.baseName + ".dehost.fastq"
    fasta_out = fastq_in.baseName + ".dehost.fasta"

    log.info "Running dehost extract seqs into fastq/fasta on sample_id: ${sample_id}"
    """
    { cat ${fastq_in} | extractreads_byname_fq.pl -v ${hosts} > ${fastq_out}; } \
    || echo "ERROR: problems script to extract listed sequences from fastq file!!";

    fastq_to_fasta -Q33 -n -i ${fastq_out} -o ${fasta_out};
    """
}
process dehost_qc {
    tag "$sample_id"
    label 'utils'

    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.QC.log"
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move", \
        saveAs: { filename -> filename.endsWith("QC.log") ? null : filename }

    input:
    tuple sample_id, file(in_fq) from cleanfq_count_qc
    output:
    file("*dehost.QC.log") into QC_log_dehost

    script:
    """
    QC=\$( echo "QC_nonhost" ${sample_id} \$(echo \$(cat ${in_fq}|wc -l)/4|bc) )
    echo "\${QC}"|tee -a ${sample_id}.5.dehost.QC.log;
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Spike Seq Check ~~~~~
// Check parameters are numbers
number_vars = [ 'spikesThreads', 'spikesCompressDbs', ]
CheckParams.is_numeric(params, number_vars)
process spike_checks {
    tag "$sample_id"
    label 'spikes'
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move"

    errorStrategy 'ignore'
    //scratch true
    cpus params.spikesThreads

    when:
    params.spikedRun == true

    input:
    tuple sample_id, file(clean_fq) from cleanfq_out

    output:
    file("${sample_id}.matches.tsv") optional true into spikes_tsv
    file("*.log") optional true

    script:
    log.info "Running spike_checks on sample_id: ${sample_id}"
    threads  = params.spikesThreads
    compress = params.spikesCompressDbs
    targetdb = params.spikesTargetDb
    """
    match_target_seqs.sh \
        -t \$PWD \
        -p ${threads} \
        -c ${compress} \
        ${clean_fq} \
        ${targetdb}
    """
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Classify RDP ~~~~~
// Check parameters are numbers
number_vars = [ 'classifyCutoff', ]
CheckParams.is_numeric(params, number_vars)
process classify_rdp {
    tag "$sample_id"
    label 'rdpclass'
    publishDir "${params.outdir}/classify/${flowcell}/", pattern: "*.{classified,taxa}.csv", mode: "move", \
        saveAs: { filename -> filename.startsWith("cnadjusted") ? null : filename }
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move"

    input:
    tuple sample_id, file(fasta_in) from cleanfa_out

    output:
    val sample_id into sample_classified
    tuple sample_id, file("*.{classified,taxa}.csv") into classified_ch
    file("*.log") optional true

    script:
    class_out = sample_id + ".classified.tsv"
    class_csv = sample_id + ".classified.csv"
    hier_out  = sample_id + ".taxa.tsv"
    hier_csv  = sample_id + ".taxa.csv"
    cutoff    = params.classifyCutoff
    log.info "Running RDP classify on sample_id: ${sample_id} using cutoff ${cutoff}"
    """
    rdp_classifier      \
        -c ${cutoff}    \
        -f filterbyconf \
        -o ${class_out} \
        -h ${hier_out}  \
        ${fasta_in}
    # convert TSV to CSV
    tr "\t" "," < ${class_out} > ${class_csv} && rm ${class_out}
    tr "\t" "," < ${hier_out}  > ${hier_csv}  && rm ${hier_out}
    """
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Summation Channels ~~~~~
// Contents of samples' spikes match reports
spikes_tsv
    .collectFile(name: "${projCode}.${flowcell}.samples_completed.spike_list",
                 newLine: true,
                 tempDir: params.outdir)
    .set { samples_spike_final }

// List of samples ready for final gzip
sample_classified
    .collectFile(name: "${projCode}.${flowcell}.samples_completed.gzip_list",)
    .set { samples_gzip_final }

// samples' QC read counts
QC_log_raw.concat( QC_log_trim, QC_log_flash, QC_log_uchime, QC_log_dehost )
    .collectFile(name: "${projCode}.${flowcell}.samples_completed.qc_list",
                 tempDir: params.outdir)
    .set { samples_qc_final }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Summary Processes  ~~~~~
process spike_summary {
    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.tsv", mode: "copy"

    when:
    params.spikedRun == true

    input:
    file(spike_pcts) from samples_spike_final
    val(flow) from flowcell
    val(proj) from projCode

    output:
    file(spikes) into spike_summary_tsv

    script:
    spikes = "pipe_16S_spike_pcts-${flow}-${proj}.tsv"
    """
    [ -s ${spike_pcts} ] && \
    mv ${spike_pcts} ${spikes}
    """
}

process read_counts {
    label 'perl5'
    publishDir "${params.outdir}/qc/${flowcell}/", pattern: "*.csv", mode: "copy"

    input:
    val(flow) from flowcell
    val(proj) from projCode
    val(outdir) from abs_outdir
    file(qc_reads) from samples_qc_final

    output:
    file(qc_csv) into qc_log_csv

    // # QC read counts:
    // echo "QC_raw"        ${sample_id} \$(echo \$(zcat ${reads[0]}|wc -l)/4|bc)
    // echo "QC_trim"       ${sample_id} \$(echo \$(cat ${sample_id}_R1_trim.fastq|wc -l)/4|bc)
    // echo "QC_combined"   ${sample_id} \$(echo \$(cat ${flash_fq} |wc -l)/4|bc)
    // echo "QC_nonchimera" ${sample_id} \$(echo \$(cat ${fastq_out}|wc -l)/4|bc)
    // echo "QC_nonhost"    ${sample_id} \$(echo \$(cat ${fastq_out}|wc -l)/4|bc)

    //CIVET originals:
    // echo "QC_raw" $myprefix $(echo $(zless {in_1} |wc -l )/4|bc) | tee -a $myprefix.QC.log;
    // echo "QC_trim" $myprefix $(echo $(zless {out_1} |wc -l )/4|bc) | tee -a $myprefix.QC.log
    // echo "QC_combined" $myprefix $(echo $(wc -l $myY |cut -d" " -f1 )/4|bc) | tee -a $myprefix.QC.log
    // echo "QC_nonchimera" $myprefix $(echo $(wc -l {out_5} |cut -d" " -f1 )/4|bc) | tee -a $myprefix.QC.log
    // echo "QC_nonhost" $myprefix $(echo $(wc -l {out_1} |cut -d" " -f1 )/4|bc) | tee -a $myprefix.QC.log

    script:
    qcdir = "${outdir}/qc/${flowcell}"
    qc_csv = "pipe_16S_QC-${flow}-${proj}.csv";
    """
    mv ${qc_reads} QC.log;
    perl ${workflow.projectDir}/bin/qc_16S_reads.pl;
    perl ${workflow.projectDir}/bin/qc_16S_alert.pl;
    rm -rf ${qcdir}/*.QC.log; # */
    mv QC.log.csv ${qc_csv};
    """
}

process final_gzips {
    label 'utils'
    publishDir "${params.outdir}/logs/${flowcell}/", pattern: "*.log", mode: "move"

    errorStrategy 'ignore'
    echo true

    input:
    val(samples)  from samples_gzip_final
    val(flowcell) from flowcell
    val(outdir)   from abs_outdir

    output:
    file("*.log") optional true

    shell:
    '''
    for F in \
      !{outdir}/logs/!{flowcell}/*.log   \
      !{outdir}/*/!{flowcell}/*.fast[aq] \
      !{outdir}/spikes/!{flowcell}/*.tsv ; #*/
    do
      gzip -v9 $F || echo "ERR $F not found" \
      >> final_gzips.log 2>&1;
    done;
    '''
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Info ~ ~ ~ ~ ~ ~
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
    println "If workflow run successful, remember to delete the Working dir: '${workflow.workDir}'"
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// vim: set ft=groovy.nextflow ts=4 sw=0 tw=100 et fdm=syntax:
