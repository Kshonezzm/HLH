#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir,$outputdir,$ciriquantdir,$configfile,$threads,$bedfile,$help);

GetOptions(
	"inputdir|i=s" => \$inputdir,
	"outputdir|o=s" => \$outputdir,
	"ciriquantdir=s" => \$ciriquantdir,
	"configfile=s" => \$configfile,
	"threads=s" => \$threads,
	"bedfile=s" => \$bedfile,
	"help!" => \$help,
);


my @samples = `find $inputdir -type d -name "GSE*" -exec sh -c 'find {} -path "*/bwa*/*" -prune -o -type f -name "SRR*" -print' \\;`;
print join("\n",@samples)."\n";
foreach my $sample_p1 (@samples){
	chomp $sample_p1;
	$sample_p1 =~ /.*\/(GSE.*?)\/(SRR.*?)\//;
    my $sample_id = $2;
    my $dataset_id = $1;

	if(!-e "$outputdir/$ciriquantdir/$sample_id"){
		mkpath("$outputdir/$dataset_id/$ciriquantdir/$sample_id",0644);
		if($@){
			print "Make path $outputdir/$dataset_id/$ciriquantdir/$sample_id failed:\n$@";
			exit(1);
		}
	}
	
	open(SH,">$outputdir/$dataset_id/$ciriquantdir/$sample_id/${sample_id}_circRNAquantify.sh") or die "$!\n";
    print SH "CIRIquant -t 40 -1 $outputdir/$dataset_id/$sample_id/$sample_id\_1.fastq -2 $outputdir/$dataset_id/$sample_id/$sample_id\_2.fastq --config $configfile -o $outputdir/$dataset_id/$ciriquantdir/$sample_id -p $sample_id --bed $bedfile\n";
	close SH;
    
    my $taskNum =`ps -aux | grep CIRIquant | wc -l`; 
    while($taskNum > 10){
        print "The num of task remaining $taskNum\n";
        sleep 30;
        print `date`;
        $taskNum = `ps -aux | grep CIRIquant | wc -l`;
    }
	
	my $out = system("sh $outputdir/$dataset_id/$ciriquantdir/$sample_id/${sample_id}_circRNAquantify.sh 1>>$outputdir/$dataset_id/$ciriquantdir/$sample_id/std.log 2>>$outputdir/$dataset_id/$ciriquantdir/$sample_id/error.log &");
	if($out==0){
		print "The task of $sample_id is successfully submitted\n";
	}
}							

# perl quantifycircRNAbyCIRIquant.pl --inputdir /home/zmzhang/test --outputdir /home/zmzhang/test --ciriquantdir ciriquantfile --configfile /home/zmzhang/config_circAtlasv3.0.yml --bedfile /home/storage/ncRNAanno/human/circAtlasv3.0_circRNA.bed --threads 60
