#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;

my ($inputdir, $outputdir, $threads, $help);

GetOptions(
    "inputdir|i=s"  => \$inputdir,
    "outputdir|o=s" => \$outputdir,
    "threads=s"     => \$threads,
    "help!"         => \$help,
);

if ($help) {
    print "Usage: $0 --inputdir|-i <input_directory> --outputdir|-o <output_directory> --threads <number_of_threads>\n";
    exit;
}

sub monitor_tasks {
    my $limit = shift;
    my $task_name = shift;
    my $task_num = `ps -aux | grep $task_name | grep -v grep | wc -l`;
    chomp($task_num);
    while ($task_num >= $limit) {
        print "The number of $task_name tasks running: $task_num\n";
        sleep 30;
        print `date`;
        $task_num = `ps -aux | grep $task_name | grep -v grep | wc -l`;
        chomp($task_num);
    }
}

my @samfiles = `find $inputdir -name "accepted_hits.sam"`;
foreach my $samfile (@samfiles) {
    chomp $samfile;
    $samfile =~ /.*\/(.*)\/accepted_hits.sam/;
    my $sample_id = $1;
    my $samdir = "$inputdir/$sample_id";
    my $bamfile = "$outputdir/$sample_id/accepted_hits_NHi1.sorted.unique.bam";
    my $sorted_bamfile = "$outputdir/$sample_id/accepted_hits_NHi1.sorted.unique.sorted.bam";

    # Create output directory if it doesn't exist
    mkpath("$outputdir/$sample_id");

    # SAM to BAM conversion
    system("samtools view -bS $samfile -o $bamfile");

    # Sort BAM file
    system("samtools sort $bamfile -o $sorted_bamfile");

    # Index the sorted BAM file
    system("samtools index $sorted_bamfile");

    # Prepare the BAM to BigWig conversion script
    open(SH, ">$outputdir/$sample_id/${sample_id}_bamtobw.sh") or die "$!\n";
    print SH "bamCoverage -b $sorted_bamfile -of bigwig --binSize 10 --ignoreDuplicates --normalizeUsing BPM --numberOfProcessors $threads -o $outputdir/$sample_id/$sample_id.bw\n";
    close SH;

    monitor_tasks(4, 'bamCoverage');

    my $out = system("sh $outputdir/$sample_id/${sample_id}_bamtobw.sh 1>>$outputdir/$sample_id/std.log 2>>$outputdir/$sample_id/error.log &");
    if ($out == 0) {
        print "The task of $sample_id is successfully submitted\n";
    } else {
        print "There was an error in submitting the task for $sample_id\n";
    }
}

# Find all BAM files and create BigWig conversion scripts
my @bamfiles = `find $inputdir -name "accepted_hits_NHi1.sorted.unique.sorted.bam"`;
foreach my $bamfile (@bamfiles) {
    chomp $bamfile;
    $bamfile =~ /.*\/(.*)\/accepted_hits_NHi1.sorted.unique.sorted.bam/;
    my $sample_id = $1;
    open(SH,">$inputdir/$sample_id/${sample_id}_bamtobw.sh") or die "$!\n";
    print SH "samtools index $bamfile\n";
    print SH "bamCoverage -b $bamfile -of bigwig --binSize 10 --ignoreDuplicates --normalizeUsing BPM --numberOfProcessors $threads -o $inputdir/$sample_id/$sample_id.bw\n";
    close SH;

    monitor_tasks(4, 'bamCoverage');

    my $out = system("sh $inputdir/$sample_id/${sample_id}_bamtobw.sh 1>>$inputdir/$sample_id/std.log 2>>$inputdir/$sample_id/error.log &");
    if ($out == 0) {
        print "The task of $sample_id is successfully submitted\n";
    } else {
        print "There was an error in submitting the task for $sample_id\n";
    }
}

# perl bamtobw3.pl --inputdir /home/zlli/human/rawdata/Post-mortem/GSE144136/bwafile --outputdir /home/zlli/human/rawdata/Post-mortem/GSE144136/bwfile --threads 4