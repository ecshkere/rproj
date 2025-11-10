#!/usr/bin/perl
($sample, $bam_file) = @ARGV;

if ($sample =~ /RNASEQ_(s\w+)/) { $sample = $1; }
elsif ($sample =~ /CHIPSEQ_(s\w+)/) { $sample = $1; }
elsif ($sample eq "Damarov") { $sample = "sDm"; }
elsif ($sample eq "E.Viacheslavovna") { $sample = "sAn"; }
elsif ($sample eq "N.Petrovna") { $sample = "sBn"; }

open(OUT,">$sample.stat");
print OUT join("\t", "id_sample", "chr", "pos", "DP", "ref", "QS_ref", "alt1", "QS_alt1", "alt2", "QS_alt2"), "\n";
$genome = "/media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa";

open(IN, "bcftools mpileup -f $genome -d 10000 -q 30 -Q 20 $bam_file --ff QCFAIL,SECONDARY |");

while ($IN=<IN>) {
	$IN=~/DP=(\d+)/;
	$dp=$1;
	# if ($dp>50) {
	if ($dp > 19) {
		@IN=split(/\t/,$IN);
		if($IN[4]=~/,/) {
			$IN[7]=~/QS=([\d\.,]+)/;
			@QS=split(/,/,$1);
			if($QS[1]>0.1 && $QS[1]<0.9) {
				@l=split(/,/,$IN[4]);
				print OUT	"$sample	$IN[0]	$IN[1]	$dp	$IN[3]	$QS[0]	$l[0]	$QS[1]	$l[1]	$QS[2]\n";
			}
		} 
	}
}
