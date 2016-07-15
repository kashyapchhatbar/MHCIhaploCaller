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
my $prefix;
GetOptions(
	'input=s' => \$in, 
	'output=s' => \$out,
	'prefix=s' => \$prefix,
) or die "Usage:\n--input\tunknown fasta file of all unknown sequences\n--output\tunknown fasta file of unique sequences\n--prefix\tprefix for unique unknown sequences\n";

if(! defined $in or ! defined $out or ! defined $prefix){
	die "Usage:\n--input\tunknown fasta file of all unknown sequences\n--output\tunknown fasta file of unique sequences\n--prefix\tprefix for unique unknown sequences\n";
}

open(IN, "$in") or die "Cannot open $in in findUniqueUnknown.pl\n";
open(OUT, ">$out") or die "Cannot write $out in findUniqueUnknown.pl\n";

my $id;
my $sequence;
my %seq;
my %id;
my %length;

while (<IN>){
	chomp $_;
	if (/^>(.+)/){
		$id = $1;
		next;
	}
	else{
		$sequence = $_;
		if (exists $seq{$sequence}){
			$seq{$sequence} = $seq{$sequence} + 1;
			$id{$sequence} = $id{$sequence}.",".$id;
		}
		else{
			$length{$sequence} = length($_);
			$seq{$sequence} = 1;
			$id{$sequence} = $id;
		}
		$sequence = "";
		next;
	}
}

my $i = 0;
foreach my $sequence  (sort { $seq{$b} <=> $seq{$a} } keys %seq){
	$i++;
	print OUT ">$prefix$i\n$sequence\n";
}





