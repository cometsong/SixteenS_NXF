#!/usr/bin/perl

#lchen@jax.org
#read the QC.log.csv in current dir and output QC.log.alert.csv
# QC.log.alert.csv is based on QC.log.csv, with a column "alert",
# for 16S, it's based on number of non-chimera reads

#example QC.log.csv below
#  Sample_name,QC_raw,QC_trim,QC_combined,QC_nonchimera,QC_nonhost
#  HMP2_J08498_1_ST_T0_B0_0000_69_001_6011_AF28W_L001,571046,571046,106359,81167



my $very_few=5000;
my $too_few=2000;
open(F,"QC.log.csv");
my $header=<F>;
chomp $header;
open(OUT,">QC.log.alert.csv");
print OUT "#very few means less than $very_few non-host, too few means less than $too_few non-host\n";
print OUT $header.",Alert\n";
while(<F>){
    chomp;
    my @cols=split /,/;
    if($cols[5]<$too_few){
	print OUT $_.",TOO_FEW\n";
    }elsif($cols[5]<$very_few){
	print OUT $_.",VERY_FEW\n";
    }
}
close(F);
close(OUT);
