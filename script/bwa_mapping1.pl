#####################################################################################################
#	File Name: hisat_mapping_ZhangLiInHouse.pl
#	This program calculates gene expression abundance by tophat2, htseq and cufflinks software
#	> Author: Hua Yu
#	> Mail: huayu@genetics.ac.cn 
#	Created Time: Sun 08 April 2018 02:40:40 PM CST
#####################################################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$indexdir,$bwadir,$picarddir,$threads,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"indexdir=s" => \$indexdir,
	"bwadir=s" => \$bwadir,
	"picarddir=s" => \$picarddir,
	"threads=s" => \$threads,
	"help!" => \$help,
);

my @samples = `find "$inputdir" -type d -name "GSE*" -exec sh -c 'find {} -path "*/hisat2*/*" -prune -o -type f -name "SRR*" -print' \\;`;
print join("\n",@samples)."\n";
foreach my $sample_p1 (@samples){
	chomp $sample_p1;
	$sample_p1 =~ /.*\/(GSE.*?)\/(SRR.*?)\//;
    my $sample_id = $2;
    my $dataset_id = $1;

	if(!-e "$outputdir/$bwadir/$sample_id"){
		mkpath("$outputdir/$dataset_id/$bwadir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$dataset_id/$bwadir/$sample_id failed:\n$@";
			exit(1);
		}
	}
	
	open(SH,">$outputdir/$dataset_id/$bwadir/$sample_id/${sample_id}_mapping.sh") or die "$!\n";
	#    if(!-e "$outputdir/$dataset_id/$bwadir/$sample_id/accepted_hits.sam" || -z "$outputdir/$bwadir/$sample_id/accepted_hits.sam"){
	    #    print SH "fastq-dump --split-files $sample_p1 -O $inputdir/$dataset_id/$sample_id\n";
	#    }
	if(!-e "$outputdir/$dataset_id/$bwadir/$sample_id/accepted_hits.sam" || -z "$outputdir/$dataset_id/$bwadir/$sample_id/accepted_hits.sam"){
       print SH "bwa mem -t $threads $indexdir $inputdir/$dataset_id/$sample_id/${sample_id}_1.fastq $inputdir/$dataset_id/$sample_id/${sample_id}_2.fastq 1> $outputdir/$dataset_id/$bwadir/$sample_id/sample.sam 2> $outputdir/$dataset_id/$bwadir/$sample_id/sample.log\n";
	}
	close SH;
    
	   my $taskNum =`ps -aux | grep perl | grep bwa | wc -l`; 
 	   while($taskNum > 5){
 	       print "The num of task remaining $taskNum\n";
 	       sleep 30;
 	       print `date`;
 	       $taskNum = `ps -aux | grep perl | grep bwa | wc -l`;
 	   }
	my $out = system("sh $outputdir/$dataset_id/$bwadir/$sample_id/${sample_id}_mapping.sh 1>>$outputdir/$dataset_id/$bwadir/$sample_id/std.log 2>>$outputdir/$dataset_id/$bwadir/$sample_id/error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}							

# nohup perl bwa_mapping1.pl --inputdir /home/zmzhang/test --outputdir /home/zmzhang/test --indexdir /home/zmzhang/bwahumanindex/humanbwaIndex/GRCh38.p12.genome.fa --bwadir bwafile --threads 2 &
