#!/usr/bin/perl

#lchen@jax.org
#script works with our 16S civet pipeline.
#Read QC.log in current dir and output QC.log.csv.
#extract QC information from individual samples'
#	QC output file, and put them into one csv. 
#The individual QC files are first concatenated together
#	into QC.log (outside this scrip), in the format of 
#sample_name   item_name value
#One of these lines becomes an entry with value "value""in
#	the final csv table under Column "item_name", for Row "sample_name".
#Has the flexibility to work with other pipeline.... 


open(F,"QC.log");
my @aitem;
my %hitem;
my %h;
while(<F>){
    chomp;
    next if(/^$/);
    my ($item ,$sample, $x)=split / /;
    if( !defined $hitem{$item}){
	$hitem{$item}=1;
	push @aitem, $item;
    }
}

seek(F,0,0);
open(OUT,">QC.log.csv");
print OUT "Sample_name";
foreach my $item (@aitem){
    print OUT ",";
    print OUT $item;
}
print OUT "\n";

while(<F>){
    chomp;
    next if(/^$/);
    my ($item ,$sample, $x)=split / /;
    if( !defined $h{$sample}){
	$h{$sample}={};
    }
    $h{$sample}->{$item}=$x;
#    print "debug1 $sample $item $x ".$h{$sample}->{$item}."\n";
}

foreach my $sample (keys %h){
    my $s=$sample;
    foreach my $item (@aitem){
    #	print "debug2 $sample $item\n";
        $s.=",";
        if(defined $h{$sample}->{$item}){
            $s.=$h{$sample}->{$item};
        }else{
            $s.="Missing";
        }
    }
    print OUT $s;
    print OUT "\n";
}
