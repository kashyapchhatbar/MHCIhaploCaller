#!/usr/bin/env perl

# authors	= "Deepali Vasoya"
# copyright	= "Copyright 2016"
# version	= "0.1"
# maintainer	= "Deepali Vasoya"
# email	= "Deepali.Vasoya@staffmail.ed.ac.uk"
# status	= "Alpha"

# #################################################################

# mhcI_haplotype.pl

# Copyright (c) 2015 Deepali Vasoya

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This file is part of MHCI haplotype finder
#
# MHCI haplotype finder is the pipeline of perl scripts
# to analyse Illumina high throughput data of MHCI 
# alleles in cattle. 
#
# This script uses only pre-defined haplotype information 
# based on MHCI genes deposited in IPD database. 
# We are making it more generic tool in near future. 
# Please refer README.md for more detail on usage.

use warnings;
use Getopt::Long;
use Cwd;

my $command = shift or die "Please enter command:\npreAlign - to analyse the paired-end raw/quality trimmed fastq files.\nBlast - to align processed unique aligned reads on database\npostAlign - to analyse the blast results\nanalyseUnknowns - to find out which sample has unknown sequences\n\n";

my $inDir;
my $total_samples;
my $dir;
my $flashPath;
my $fastxPath;

my $blastPath;
my $db;
my $onlySelected;
my $primer;
my $workDir;
my $read1;
my $read2;

my $main_dir = getcwd;

if ($command eq "preAlign"){
	
	GetOptions(

		'fastqDir=s' => \$inDir, 
		'total_samples=i' => \$total_samples, 
		'prefix=s' => \$prefix,
		'dir=s' => \$dir, 
		'FlashPath=s' => \$flashPath,
		'fastxToolkitPath=s' => \$fastxPath,

	) or die "Usage:\n--fastqDir\t/Directory/of/Fastq/Files\n--total_samples\tTotal number of samples\n--prefix\tprefix of the fastq files, Boran, Cameroonian or HolsteinFriesian\n--dir\t/Full/path/to/project/directory\n--FlashPath\t/Full/Path/of/FLASH/software\n--fastxToolkitPath\t/Full/Path/of/fastx-toolkit/software\n\n";
	
	if (! defined $inDir or ! defined $total_samples or ! defined $dir or ! defined $flashPath){
		die "Usage:\n--fastqDir\t/Directory/of/Fastq/Files\n--total_samples\tTotal number of samples\n--prefix\tprefix of the fastq files, Boran, Cameroonian or HolsteinFriesian\n--dir\t/Full/path/to/project/directory\n--FlashPath\t/Full/Path/of/FLASH/software\n--fastxToolkitPath\t/Full/Path/of/fastx-toolkit/software\n\n";
	}

	for (my $n = 1; $n <= $total_samples; $n++){

		print "\nRunning Sample $n\n";

		$workDir = "$dir/Samples/$n";
		$read1 = "$inDir/$prefix"."-$n"."-R1.fastq.gz";
		$read2 = "$inDir/$prefix"."-$n"."-R2.fastq.gz";
		
		system ("mkdir -p $workDir");
		system ("mkdir -p $workDir/Flash");

		print "**** Running FLASH for Sample $n using $read1 and $read2 ****\n";
		system ("$flashPath/flash $read1 $read2 -m 20 -M 300 -o $n -d $workDir/Flash > $workDir/Flash/$n.flash.log");
		
		print "**** Converting fastq into fasta: Sample $n...\n";
		system ("$fastxPath/bin/fastq_to_fasta -r -v -Q33 -i $workDir/Flash/$n.extendedFrags.fastq -o $workDir/$n.extendedFrags.fasta > $workDir/$n.fastqToFasta.log");

		print "**** Sorting primers in Sample $n...\n";
		system ("perl $main_dir/scripts/sortPrimers.pl --input $workDir/$n.extendedFrags.fasta --outDir $workDir --prefix $n");

		print "**** Grouping unique sequences in Sample $n For1Rev2...\n";
		system ("perl $main_dir/scripts/findUniqueSeq.pl --input $workDir/$n.For1Rev2.fasta --outDir $workDir --cutoff 0.2 --prefix $n --primer For1Rev2");
		print "**** Grouping unique sequences in Sample $n For3Rev1...\n";
		system ("perl $main_dir/scripts/indUniqueSeq.pl --input $workDir/$n.For3Rev1.fasta --outDir $workDir --cutoff 0.2 --prefix $n --primer For3Rev1");
	
	}
	system("mkdir -p $dir/Reports");
	print "**** Generating reports for primer sorting and unique sequences...\n";
	system ("cat $dir/Samples/*/*.sortPrimer.log | sort -nk1 > $dir/Reports/sortPrimer.log");
	system ("cat $dir/Samples/*/*.For1Rev2.uniqueSeq.log | sort -nk1 > $dir/Reports/For1Rev2.uniqueSeq.log");
	system ("cat $dir/Samples/*/*.For3Rev1.uniqueSeq.log | sort -nk1 > $dir/Reports/For3Rev1.uniqueSeq.log");

}
elsif ($command eq "Blast"){
	
	GetOptions(
		
		'blastPath=s' => \$blastPath,
		'dir=s' => \$dir,
		'total_samples=i' => \$total_samples, 
		'onlySelected=s' => \$onlySelected,

	) or die "Usage:\n--blastPath\t/Full/Path/of/Blast/software\n--dir\t/Full/path/to/project/directory\n--total_samples\tTotal number of samples\n--onlySelected\talign only selected sequences which are above threshold, TRUE or FALSE\n\n";
	
	if (! defined $blastPath or ! defined $dir or ! defined $total_samples or ! defined $onlySelected){
		die "Usage:\n--blastPath\t/Full/Path/of/Blast/software\n--dir\t/Full/path/to/project/directory\n--total_samples\tTotal number of samples\n--onlySelected\talign only selected sequences which are above threshold, TRUE or FALSE\n\n";
	}
	if ($onlySelected eq "TRUE" or $onlySelected eq "true" or $onlySelected eq "True"){

		for ($n = 1; $n <= $total_samples; $n++){

			$workDir = "$dir/Samples/$n";
			print "Running BLAST for sample $n for $primer only selected sequences...\n";
			system("$blastPath/bin/blastall -p blastn -d $main_dir/db/For1Rev2.fasta -i $workDir/$n.For1Rev2.selected.fasta -m 8 -o $workDir/$n.For1Rev2.blast");
			system("$blastPath/bin/blastall -p blastn -d $main_dir/db/For3Rev1.fasta -i $workDir/$n.For3Rev1.selected.fasta -m 8 -o $workDir/$n.For3Rev1.blast");
		
		}
	}
	elsif ($onlySelected eq "FALSE" or $onlySelected eq "false" or $onlySelected eq "False"){

		for ($n = 1; $n <= $total_samples; $n++){

			$workDir = "$dir/Samples/$n";
			print "Running BLAST for sample $n for $primer...\n";
			system("$blastPath/bin/blastall -p blastn -d $main_dir/db/For1Rev2.fasta -i $workDir/$n.For1Rev2.unique.fasta -m 8 -o $workDir/$n.For1Rev2.blast");
			system("$blastPath/bin/blastall -p blastn -d $main_dir/db/For3Rev1.fasta -i $workDir/$n.For3Rev1.unique.fasta -m 8 -o $workDir/$n.For3Rev1.blast");
		
		}
	}

	else{

		die "Usage:\n--onlySelected\talign only selected sequences which are above threshold, TRUE or FALSE\n\n";
	
	}

}
elsif ($command eq "postAlign"){

	GetOptions(

		'dir=s' => \$dir,
		'total_samples=i' => \$total_samples, 

	) or die "Usage:\n--dir\t/Full/path/to/project/directory\n--total_samples\tTotal number of samples\n\n";
	
	if (! defined $dir or ! defined $total_samples){
		die "Usage:\n--dir\t/Full/path/to/project/directory\n--total_samples\tTotal number of samples\n\n";
	}

	for ($n = 1; $n <= $total_samples; $n++){	

		print "**** Analysing Blast of sample $n...\n";
		$workDir = "$dir/Samples/$n";

		system ("perl $main_dir/scripts/analyseBlast.pl --dir $workDir --prefix $n --primer For1Rev2");
		system ("perl $main_dir/scripts/analyseBlast.pl --dir $workDir --prefix $n --primer For3Rev1");

		print "**** Analysing known alleles and assigning haplotypes..\n";
		system ("perl $main_dir/scripts/analyseKnownGenes.pl --dir $workDir --prefix $n --primer For1Rev2");
		system ("perl $main_dir/scripts/analyseKnownGenes.pl --dir $workDir --prefix $n --primer For3Rev1");

		print "**** Generating combining reports...\n";
		system ("cat $dir/Samples/*/*.For1Rev2.unknown.fasta > $dir/Reports/Unknown.For1Rev2.fasta");
		system ("cat $dir/Samples/*/*.For3Rev1.unknown.fasta > $dir/Reports/Unknown.For3Rev1.fasta");
		system ("cat $dir/Samples/*/*.For1Rev2.analyseBlast.log | sort -nk1,1 > $dir/Reports/For1Rev2.analyseBlast.log.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.analyseBlast.log | sort -nk1,1 > $dir/Reports/For3Rev1.analyseBlast.log.txt");
		system ("cat $dir/Samples/*/*.For1Rev2.haplotypes.genes.txt | sort -nk1,1 > $dir/Reports/For1Rev2.haplotypes.genes.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.haplotypes.genes.txt | sort -nk1,1 > $dir/Reports/For3Rev1.haplotypes.genes.txt");
		system ("cat $dir/Samples/*/*.For1Rev2.haplotypes.txt | sort -nk1,1 > $dir/Reports/For1Rev2.haplotypes.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.haplotypes.txt | sort -nk1,1 > $dir/Reports/For3Rev1.haplotypes.txt");
		system ("cat $dir/Samples/*/*.For1Rev2.analyseKnownGenes.log | sort -nk1,1 > $dir/Reports/For1Rev2.analyseKnownGenes.log.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.analyseKnownGenes.log | sort -nk1,1 > $dir/Reports/For3Rev1.analyseKnownGenes.log.txt");
		system ("cat $dir/Samples/*/*.For1Rev2.extra.genes.txt | sort -nk1,1 > $dir/Reports/For1Rev2.extra.genes.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.extra.genes.txt | sort -nk1,1 > $dir/Reports/For3Rev1.extra.genes.txt");
		system ("cat $dir/Samples/*/*.For1Rev2.nonclassic.genes.txt | sort -nk1,1 > $dir/Reports/For1Rev2.nonclassic.genes.txt");
		system ("cat $dir/Samples/*/*.For3Rev1.nonclassic.genes.txt | sort -nk1,1 > $dir/Reports/For3Rev1.nonclassic.genes.txt");

	}
}

elsif ($command eq "analyseUnknowns"){

	GetOptions(

		'dir=s' => \$dir,
		'total_samples=i' => \$total_samples, 

	) or die "Usage:\n--dir\t/Full/path/to/main/directory\n--total_samples\tTotal number of samples\n\n";
	
	if (! defined $dir or ! defined $total_samples){
		die "Usage:\n--dir\t/Full/path/to/main/directory\n--total_samples\tTotal number of samples\n\n";
	}

	print "**** Grouping unique sequences from all unknown sequences for For1Rev2...\n";
	system ("perl $main_dir/scripts/findUniqueUnknown.pl --input $dir/Reports/Unknown.For1Rev2.fasta --output $dir/Reports/Unknown.For1Rev2.unique.fasta --prefix Foo");
	print "**** Grouping unique sequences from all unknown sequences for For3Rev1...\n";
	system ("perl $main_dir/scripts/findUniqueUnknown.pl --input $dir/Reports/Unknown.For3Rev1.fasta --output $dir/Reports/Unknown.For3Rev1.unique.fasta --prefix Bar");

	print "**** Allocating new unknown sequences ids to samples....\n";
	system ("perl $main_dir/scripts/allocateUnknownFastaIds.pl --dir $dir --input $dir/Reports/Unknown.For1Rev2.unique.fasta --output $dir/Reports/Unknown.For1Rev2.unique.occurance.report.txt --primer For1Rev2 --total_samples $total_samples");
	system ("perl $main_dir/scripts/allocateUnknownFastaIds.pl --dir $dir --input $dir/Reports/Unknown.For3Rev1.unique.fasta --output $dir/Reports/Unknown.For3Rev1.unique.occurance.report.txt --primer For3Rev1 --total_samples $total_samples");

	print "**** Generating reports of unknown sequences...\n";
	system ("cat ${dir}/Samples/*/*.For1Rev2.unknown.genes.txt | sort -nk1 > $dir/Reports/For1Rev2.unknown.genes.txt");
	system ("cat ${dir}/Samples/*/*.For3Rev1.unknown.genes.txt | sort -nk1 > $dir/Reports/For3Rev1.unknown.genes.txt");

}
else{
	
	die "Please enter command:\npreAlign - to analyse the paired-end raw fastq files.\nBlast - to align processed unique aligned reads on database\npostAlign - to analyse the blast results\nanalyseUnknowns - to find out which sample has unknown sequences\n\n";

}

















