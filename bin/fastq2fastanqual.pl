#!/usr/bin/perl
## copied from http://edwards.sdsu.edu/labsite/index.php/bas/309-converting-fastq-to-fasta-and-qual

use strict;

my $usage = "Use: $0 input.fastq\n";

if (scalar @ARGV != 1) { die $usage; }
if (! (open (IN, "<$ARGV[0]"))) { die "Can't open $ARGV[0]: $!\n"; }
if (-e "$ARGV[0].fasta") { die "$ARGV[0].fasta extists - please remove or rename\n"; }
if (-e "$ARGV[0].fasta.qual") { die "$ARGV[0].fasta.qual extists - please remove or rename\n"; }
if (! (open (FASTA, ">$ARGV[0].fasta"))) { die "Can't write to $ARGV[0].fasta: $!\n"; }
if (! (open (QUAL, ">$ARGV[0].fasta.qual"))) { die "Can't write to $ARGV[0].fasta.qual: $!\n"; }

my $id = "";
my $seq = "";
my $qual = "";
my $line = <IN>;
if ($line !~ /\@/) { die "$ARGV[0] does not look like a FASTQ file (it should start with \@)\n"; }
while ($line =~ s/^\@//o) {
  chomp ($line);
  my $id = $line;
  $line = <IN>;
  while ($line !~ s/^\+//o) {
    chomp ($line);
    $seq .= $line;
    $line =<IN> ; }
  chomp ($line);
  if (($line =~ /\S/o) && ($line ne $id)) { die "ID of $id not followed by same identifier for quality scores\n"; }
  $line = <IN>;
  while (length ($qual) < length ($seq)) {
    chomp ($line);
    $qual .= $line;
    $line = <IN>; }
  print FASTA ">$id\n$seq\n";
  print QUAL ">$id\n";
  my @a = split //o, $qual;
  foreach my $i (@a) {
#    my $q = 10 * log (1 + 10 ** ((ord ($i) - 64) / 10)) / log (10);
		my $q = ord ($i) - 33;	
    printf QUAL "%1.0f ", $q; }
  print QUAL "\n";
  $id = "";
  $seq = "";
  $qual = ""; }
