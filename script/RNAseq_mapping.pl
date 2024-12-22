#####################################################################################################
#	File Name: hisat_mapping.pl
#	> Author:  Hua Yu
#	> Mail: 
#	Created Time: Sun 08 April 2024 02:40:40 PM CST
#####################################################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$indexdir,$hisat2dir,$picarddir,$htseqdir,$cufflinksdir,$stringtiedir,$gfffile,$threads,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"indexdir=s" => \$indexdir,
	"hisat2dir=s" => \$hisat2dir,
	"picarddir=s" => \$picarddir,
	"htseqdir=s" => \$htseqdir,
	"cufflinksdir=s" => \$cufflinksdir,
	"stringtiedir=s" => \$stringtiedir,
	"gfffile|gff=s" => \$gfffile,
	"threads=s" => \$threads,
	"help!" => \$help,
);

my @samples = `find "$inputdir" -type d -name "GSE*" -exec sh -c 'find {} -path "*/hisat2*/*" -prune -o -type f -name "SRR*" -print' \\;`;
print join("\n",@samples)."\n";
foreach my $sample_p1 (@samples){
	chomp $sample_p1;
	$sample_p1 =~ /.*\/(GSE.*?)\/(SRR.*?)\//;
	my $dataset_id = $1;
    my $sample_id = $2;
	
                	if(!-e "$outputdir/$stringtiedir/$sample_id"){
               		mkpath("$outputdir/$dataset_id/$stringtiedir/$sample_id",0644);
               		if($@){
               			print "Make path $outputdir/$dataset_id/$stringtiedir/$sample_id failed:\n";
               			exit(1);
               		}
               	}
               	if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id"){
               		mkpath("$outputdir/$dataset_id/$hisat2dir/$sample_id",0644);
               		if($@){
               			print "Make path $outputdir/$dataset_id/$hisat2dir/$sample_id failed:\n$@";
               			exit(1);
               		}
               	}
	
	open(SH,">$outputdir/$dataset_id/$hisat2dir/$sample_id/${sample_id}_expcal.sh") or die "$!\n";
	#			 if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
	#			 print SH "fastq-dump --split-files $sample_p1 -O $inputdir/$dataset_id/$sample_id\n";
	#			 }	
   
			 	         if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
			 	           		 print SH "hisat2 -x $indexdir -p $threads --dta --rg-id $sample_id --rg SM:$sample_id -1 $inputdir/$dataset_id/$sample_id/$sample_id\_1.fastq -2 $inputdir/$dataset_id/$sample_id/$sample_id\_2.fastq --summary-file $outputdir/$dataset_id/$hisat2dir/$sample_id/mapping_summary.txt -S $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits.sam\n";
			 	           	 }
			 	           	 if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sam"){
			 	           		 print SH "grep -v -E -w 'NH:i:2|NH:i:3|NH:i:4|NH:i:5|NH:i:6|NH:i:7|NH:i:8|NH:i:9|NH:i:10|NH:i:11|NH:i:12|NH:i:13|NH:i:14|NH:i:15|NH:i:16|NH:i:17|NH:i:18|NH:i:19|NH:i:20' $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits.sam > $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sam\n";
			 	           	 }
			 	           	 if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
			 	           	        print SH "samtools view -bS $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sam -o $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
			 	           	        print SH "samtools sort -@ $threads -o $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
			 	           	 }
			 	           	 if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
			 	           	        print SH "java -Xmx15g -jar $picarddir/picard.jar MarkDuplicates I=$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam O=$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam METRICS_FILE=$outputdir/$dataset_id/$hisat2dir/$sample_id/${sample_id}.metricsFile VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate\n";
			 	           	 }
			 	           	 if(!-e "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam" || -z "$outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam"){
			 	           		 print SH "samtools index $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam\n";
							  	 print SH "rm $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.bam\n";
							 	 	 print SH "rm $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.bam\n";
							  print SH "rm $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sam\n";
							  print SH "rm $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits.sam\n";
							  #	 print SH "rm $inputdir/$dataset_id/$sample_id/$sample_id\_1.fastq\n";
			    		 		  #	 print SH "rm $inputdir/$dataset_id/$sample_id/$sample_id\_2.fastq\n";
					  }
				   	 if(!-e "$outputdir/$dataset_id/$stringtiedir/$sample_id/transcripts.gtf" || -z "$outputdir/$dataset_id/$stringtiedir/$sample_id/transcripts.gtf"){
				        	 print SH "stringtie -p $threads -e -B -G $gfffile -A $outputdir/$dataset_id/$stringtiedir/$sample_id/gene_abund.tab -o $outputdir/$dataset_id/$stringtiedir/$sample_id/transcripts.gtf $outputdir/$dataset_id/$hisat2dir/$sample_id/accepted_hits_NHi1.sorted.unique.bam\n";
				    	 }
           	close SH;
    
		#	   my $taskNum =`ps -aux | grep fastq-dump | wc -l`; 
	 	#	   my $   while($taskNum > 30){
	 	#	   my $       print "The num of task remaining $taskNum\n";
	 	#	   my $       sleep 30;
	 	#	   my $       print `date`;
	 	#	   my $       $taskNum = `ps -aux | grep perl | grep hisat | wc -l`;
	 	#	   my $   }
	
		  	 my $taskNum =`ps -aux | grep perl | grep hisat | wc -l`; 
		        while($taskNum > 60){
		            print "The num of task remaining $taskNum\n";
		            sleep 30;
		            print `date`;
		            $taskNum = `ps -aux | grep perl | grep hisat | wc -l`;
		        }
	
	my $out = system("sh $outputdir/$dataset_id/$hisat2dir/$sample_id/${sample_id}_expcal.sh 1>>$outputdir/$dataset_id/$hisat2dir/$sample_id/std.log 2>>$outputdir/$dataset_id/$hisat2dir/$sample_id/error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}							

# perl hisat_mapping.pl --inputdir /home/huayu/eRNAHunter/CancerData/CMLData/human --outputdir /home/huayu/eRNAHunter/CancerData/CMLData/human --indexdir /home/huayu/genomeanno/humanhisat2Index/GRCh38 --picarddir /home/Yuhua/software/picard-2.18.2 --stringtiedir stringtiefile --hisat2dir hisat2file --gfffile /home/Yuhua/enhProj/genomeanno/gencode.v28.annotation.gff3 --threads 2
# perl RNAseq_mapping.pl --inputdir /home/zmzhang/testspace --outputdir /home/zmzhang/testspace --indexdir /home/storage/genomeanno/humanhisat2Index/GRCh38 --picarddir /home/storage/software/picard-2.18.2 --stringtiedir stringtieflie --hisat2dir hisat2file --gfffile  /home/storage/genomeanno/gencode.v44.mRNA.annotation.gff3 --threads 8
