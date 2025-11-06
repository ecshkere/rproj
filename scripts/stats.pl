#!/usr/bin/perl
# ($sample, $bam_file, $genome) = @ARGV;
($sample, $bam_file) = @ARGV;
# $sample =~ s/^[^_]*_[^_]*_//;
open(OUT,">$sample.stat");
print OUT join("\t", "id_sample", "chr", "pos", "DP", "ref", "QS_ref", "alt1", "QS_alt1", "alt2", "QS_alt2"), "\n";
# open(IN,"samtools view -bq 30 $samp |bcftools mpileup -f /media/leon/DATA/Genomes/grch38_snp/Human_hg38.fa -d 10000 -|");
$genome = "/media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa";

open(IN, "bcftools mpileup -f $genome -d 10000 -q 30 -Q 20 $bam_file |"); # DEFAULT

# open(IN, "bcftools mpileup -f $genome -d 10000 -q 30 -Q 20 $bam_file -R /media/leon/DISK2/icig/done/85_snps.tsv |");

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
