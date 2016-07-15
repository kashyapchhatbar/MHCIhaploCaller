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
my $dir;
my $cutoff;
my $prefix;
my $primer;
GetOptions(
	'input=s' => \$in, 
	'outDir=s' => \$dir, 
	'cutoff=s' => \$cutoff, 
	'prefix=s' => \$prefix,
	'primer=s' => \$primer,
) or die "Usage:\n--input\tfasta file to find unique sequences\n--outDir\tOutput direcotry\n--prefix\tprefix for output files\n--cutoff\tthreshold to select high abundant sequences\n--primer\tPrimer name\n\n";

if(! defined $in or ! defined $dir or ! defined $prefix or ! defined $cutoff or ! defined $primer){

	die "Usage:\n--input\tfasta file to find unique sequences\n--outDir\tOutput direcotry\n--prefix\tprefix for output files\n--cutoff\tthreshold to select high abundant sequences\n--primer\tPrimer name\n\n";

}

open (IN, "$in") or die "Cannot open $in in findUniqueSeq.pl\n";
open (OUT_Unique, ">$dir/$prefix.$primer.unique.fasta") or die "Cannot write $dir/$prefix.$primer.For1Rev2.fasta\n";
open (OUT_Selected_fasta, ">$dir/$prefix.$primer.selected.fasta") or die "Cannot write $dir/$prefix.$primer.For3Rev1.fasta\n";
open (OUT_selected_list, ">$dir/$prefix.$primer.selected.list.txt") or die "Cannot write $dir/$prefix.$primer.unknownPrimers.fasta\n";
open (LOG, ">$dir/$prefix.$primer.uniqueSeq.log") or die "Cannot write $dir/$prefix.$primer.uniqueSeq.log\n";

my $total = 0;
my %seq;
my %len;
my $s;
while (<IN>){
	chomp $_;
	if (/^>/){
		$total++;
		next;
	}
	else{
		$s = $_;
		if (exists $seq{$s}){
			$seq{$s} = $seq{$s} + 1;
		}
		else{
			$seq{$s} = 1;
			$len{$s} = length($s);
		}
	}
}

my $total_unique = keys(%seq);
my $threshold = ($total * $cutoff)/100;
my $count = 0;
my $per = 0;
my $count_selected = 0;
my $total_selected = 0;
foreach my $gene (sort {$seq{$b} <=> $seq{$a}} keys %seq){
	$count++;
	$per = sprintf("%.2f", ($seq{$gene} * 100)/$total);
	print OUT_Unique ">$count-$len{$gene}-$seq{$gene}-$per\n$gene\n";
	if ($seq{$gene} >= $threshold){
		print OUT_Selected_fasta ">$count-$len{$gene}-$seq{$gene}-$per\n$gene\n";
		print OUT_selected_list "$count-$len{$gene}-$seq{$gene}-$per\n";
		$count_selected++;
		$total_selected = $total_selected + $seq{$gene};
	}
}
$per = sprintf("%.2f", ($total_selected * 100)/$total);
print LOG "$prefix\t$total\t$total_unique\t$threshold\t$count_selected\t$total_selected\t$per\n";



