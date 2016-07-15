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
# We are making it more generical in near future. 
# Please refer README.txt for more detail on usage.

use warnings;
use Getopt::Long;

my $in;
my $out;
my $dir;
my $primer;
my $total_samples;
GetOptions(
	'dir=s' => \$dir, 
	'input=s' => \$in, 
	'output=s' => \$out,
	'primer=s' => \$primer,
	'total_samples=i' => \$total_samples, 
) or die "Usage:\n--dir\tpath to main directory\n--input\tFasta file of unique unknown sequences\n--output\tfile to write list of unique unknown ids and their occurances in samples\n--primer\nprimer label\n--total_samples\tTotal number of samples\n";

if(! defined $dir or ! defined $in or ! defined $out or ! defined $primer or ! defined $total_samples){
	die "Usage:\n--dir\tpath to main directory\n--input\tFasta file of unique unknown sequences\n--output\tfile to write list of unique unknown ids and their occurances in samples\n--primer\nprimer label\n--total_samples\tTotal number of samples\n";
}

open(IN, "$in") or die "Cannot open $in in allocateUnknownFastaIds.pl\n";
open(OUT, ">$out") or die "Cannot write $out in allocateUnknownFastaIds.pl\n";

my $count = 0;
my $id;
my %all_genes;
my %seq;
my @list;
my $sample_unknown_fasta;
my $out_fasta;
my $out_list;
my @words;
my $per;
my $gene;
my %all_ids;
my %unknown_id;
my $sample = 1;

while(<IN>){
	chomp $_;
    if (/^>(\S+)/){
    	$id = $1;
		$count++;
		next;
	}
    else{
	$all_genes{$id} = $count;
	$seq{$_} = $id;
	#print "$seq{$_}\t$_\n";
    }
}

while($sample <= $total_samples){
	#print "Running $sample\n";
	$sample_unknown_fasta = "$dir/Samples/$sample/$sample.$primer.unknown.fasta";
	$out_fasta = "$dir/Samples/$sample/$sample.$primer.unknown.processed.fasta";
	$out_list = "$dir/Samples/$sample/$sample.$primer.unknown.genes.txt";

	open (FASTA, ">$out_fasta");
	open (LIST_OUT, ">$out_list");
	
	print LIST_OUT "$sample";
	
	open (UN, "$sample_unknown_fasta");
	while(<UN>){
		chomp $_;
		if (/^>(.+)/){
			@words = split("-", $1);
			$id = $1;
			$per = $words[3];
			next;
		}
		else{	
			$gene = $seq{$_};
			#print "$gene\n";
			if (exists $unknown_id{$gene}){
				$all_ids{$gene} = $all_ids{$gene}.",".$id;
				$unknown_id{$gene} = $unknown_id{$gene}.",".$sample;
				#print "$gene\n";
			}
			else{	
				$all_ids{$gene} = $id;
				$unknown_id{$gene} = $sample;
				#print "$gene\n";
				
			} 
			print FASTA ">$sample $gene:$per\n$_\n";
			push @list, "$gene:$per";
		#	print LIST "$seq{$_}:$per, ";
		}
	}
	foreach my $gene (@list){
		if ($gene ne ""){
			print LIST_OUT "\t$gene";
		}
	}
	print LIST_OUT "\n";
	close LIST_OUT;
	close UN;
	close FASTA;

	@list = undef;
	$per = undef;
	$sample++;
}

my $len;
foreach $gene (sort {$all_genes{$a} <=> $all_genes{$b}} keys %all_genes){
	@words = split(",", $unknown_id{$gene});
	$len = scalar(@words);
	#print "$gene: $len\n";
	print OUT "$gene\t$len\t$all_ids{$gene}\t$unknown_id{$gene}\n";
	#print OUT "$gene\t$unknown_id{$gene}\n";
}

