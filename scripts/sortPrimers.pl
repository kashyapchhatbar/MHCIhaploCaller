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
my $prefix;
GetOptions(
	'input=s' => \$in, 
	'outDir=s' => \$dir,
	'prefix=s' => \$prefix,
) or die "Usage:\n--input\tfastq file to sort primers\n--outDir\tOutput direcotry\n--prefix\tprefix for output files\n\n";

if(! defined $in or ! defined $dir or ! defined $prefix){
	die "Usage:\n--input\tfastq file to sort primers\n--outDir\tOutput direcotry\n--prefix\tprefix for output files\n\n";
}

open (IN, "$in") or die "Cannot open $in in sortPrimers.pl\n";
open (OUT_For1Rev2, ">$dir/$prefix.For1Rev2.fasta") or die "Cannot write $dir/$prefix.For1Rev2.fasta\n";
open (OUT_For3Rev1, ">$dir/$prefix.For3Rev1.fasta") or die "Cannot write $dir/$prefix.For3Rev1.fasta\n";
open (OUT_noPrimer, ">$dir/$prefix.unknownPrimers.fasta") or die "Cannot write $dir/$prefix.unknownPrimers.fasta\n";
open (LOG, ">$dir/$prefix.sortPrimer.log") or die "Cannot write $dir/$prefix.sortPrimer.log\n";

my $for1_1 = "GTCGGCTACGTGGACGAC";
my $for1_2 = "GTCGGCTATGTGGACGAC";
my $for1_3 = "GTTGGCTACGTGGACGAC";
my $for1_4 = "GTTGGCTATGTGGACGAC";

my $rev1_1 = "GCTCCGCAGACACCTGGAG";
my $rev1_2 = "GCTCCCCAGACACCTGGAG";
my $rev1_3 = "GCTCCGCAGATACCTGGAG";
my $rev1_4 = "GCTCCCCAGATACCTGGAG";

my $for3_1 = "GGGCCAGAGTATTGGGA";
my $for3_2 = "GGGCCCGAGTATTGGGA";
my $for3_3 = "GGGCCGGAGTATTGGGA";
my $for3_4 = "GGGCTAGAGTATTGGGA";
my $for3_5 = "GGGCTCGAGTATTGGGA";
my $for3_6 = "GGGCTGGAGTATTGGGA";

my $rev2_1 = "AGGAACTACGTGGAGGGCC";
my $rev2_2 = "AGGAACTACCTGGAGGGCC";
my $rev2_3 = "AGGAACTACGTCGAGGGCC";
my $rev2_4 = "AGGAACTACCTCGAGGGCC";

my $total = 0;
my $unknown = 0;
my $for1_rev2 = 0;
my $for3_rev1 = 0;
my $header;

while (<IN>){
	chomp $_;
	if (/^>/){
		$total++;
		$header = $_;
	}
	else{
		if ($_ =~ /^$for1_1(\w+)/ or $_ =~ /^$for1_2(\w+)/ or $_ =~ /^$for1_3(\w+)/ or $_ =~ /^$for1_4(\w+)/){
			if ($1 =~ /(\w+)$rev2_1$/ or $1 =~ /(\w+)$rev2_2$/ or $1 =~ /(\w+)$rev2_3$/ or $1 =~ /(\w+)$rev2_4$/){
				print OUT_For1Rev2 "$header\n$1\n";
				$for1_rev2++;
			}
			else{
				$unknown++;
                print OUT_noPrimer "$header\n$_\n";
			}
		}
		elsif ($_ =~ /^$for3_1(\w+)/ or $_ =~ /^$for3_2(\w+)/ or $_ =~ /^$for3_3(\w+)/ or $_ =~ /^$for3_4(\w+)/ or $_ =~ /^$for3_5(\w+)/ or $_ =~ /^$for3_6(\w+)/){
			if ($1 =~ /(\w+)$rev1_1$/ or $1 =~ /(\w+)$rev1_2$/ or $1 =~ /(\w+)$rev1_3$/ or $1 =~ /(\w+)$rev1_4$/){
                print OUT_For3Rev1 "$header\n$1\n";
                $for3_rev1++;
            }
			else{
				$unknown++;
                print OUT_noPrimer "$header\n$_\n";
			}
		}
		else{
			$unknown++;
			print OUT_noPrimer "$header\n$_\n";
		}
	}
}
print LOG "$prefix\t";
print LOG "$total\t";
my $known = $for1_rev2 + $for3_rev1;
my $per_for1 = sprintf("%.2f", (100 * $for1_rev2) / $known);
my $per_for3 = sprintf("%.2f", (100 * $for3_rev1) / $known);
print LOG "$unknown\t";
print LOG "$for1_rev2\t$per_for1%\t";
print LOG "$for3_rev1\t$per_for3%\n";



