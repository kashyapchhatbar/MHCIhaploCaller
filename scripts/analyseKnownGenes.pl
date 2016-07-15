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

my $dir;
my $prefix;
my $primer;

GetOptions(
	'dir=s' => \$dir, 
	'prefix=s' => \$prefix,
	'primer=s' => \$primer,
) or die "Usage:\n--dir\tdirectory of analysed blast output\n--prefix\tprefix for output files\n--primer\tprimer label\n\n";

if(! defined $dir or ! defined $prefix or ! defined $primer){
	die "Usage:\n--dir\tdirectory of analysed blast output\n--prefix\tprefix for output files\n--primer\tprimer label\n\n";
}

open(IN,"$dir/$prefix.$primer.known.genes.txt") or die "Cannot open $dir/$prefix.$primer.known.genes.txt in analyseKnownGenes.pl\n";
open(OUT_GENES,">$dir/$prefix.$primer.haplotypes.genes.txt") or die "Cannot open $dir/$prefix.$primer.haplotypes.genes.txt in analyseKnownGenes.pl\n";
open(OUT_HAPLO, ">$dir/$prefix.$primer.haplotypes.txt") or die "Cannot open $dir/$prefix.$primer.haplotypes.txt in analyseKnownGenes.pl\n";
open(OUT_NONCLASSIC, ">$dir/$prefix.$primer.nonclassic.genes.txt") or die "Cannot open $dir/$prefix.$primer.nonclassic.genes.txt in analyseKnownGenes.pl\n";
open(OUT_NONHAPLO, ">$dir/$prefix.$primer.extra.genes.txt") or die "Cannot open $dir/$prefix.$primer.extra.genes.txt in analyseKnownGenes.pl\n";
open(LOG, ">$dir/$prefix.$primer.analyseKnownGenes.log") or die "Cannot open $dir/$prefix.$primer.analyseKnownGenes.log in analyseKnownGenes.pl\n";

my %genes;
my %genes1;
my @haplo;
my $gene;
my $flag;
my $flag1;
my $name;
my $count_classic = 0;
my $count_nonclassic = 0;
my $header;
my %nonHaplo;

sub Haplotypes{	
		if  (exists $genes{"3*00201"} and exists $genes{"2*01201"} and exists $genes{"Roslin1.1"} and exists $genes{"RoslinSZ.1.1"}){
		push @haplo, "A10";
		$genes1{"3*00201"} = "A10";
		$genes1{"2*01201"} = "A10";
		$genes1{"Roslin1.1"} = "A10";
		$genes1{"RoslinSZ.1.1"} = "A10";
		print OUT_GENES "3*00201:".$genes{"3*00201"}.",2*01201:".$genes{"2*01201"}.",Roslin1.1:".$genes{"Roslin1.1"}.",RoslinSZ.1.1:".$genes{"RoslinSZ.1.1"}."\t";
	}
	elsif (exists $genes{"3*00201"} and exists $genes{"2*01201"} and exists $genes{"Roslin1.1"}){
		push @haplo, "A10*";
		$genes1{"3*00201"} = "A10*";
		$genes1{"2*01201"} = "A10*";
		$genes1{"Roslin1.1"} = "A10*";
		print OUT_GENES "3*00201:".$genes{"3*00201"}.",2*01201:".$genes{"2*01201"}.",Roslin1.1:".$genes{"Roslin1.1"}."\t";
	}
	elsif (exists $genes{"3*00201"} and exists $genes{"2*01201"} and exists $genes{"RoslinSZ.1.1"}){
		push @haplo, "A10**";
		$genes1{"3*00201"} = "A10**";
		$genes1{"2*01201"} = "A10**";
		$genes1{"RoslinSZ.1.1"} = "A10**";
		print OUT_GENES "3*00201:".$genes{"3*00201"}.",2*01201:".$genes{"2*01201"}.",RoslinSZ.1.1:".$genes{"RoslinSZ.1.1"}."\t";
	}
	elsif (exists $genes{"3*00201"} and exists $genes{"2*01201"}){
		push @haplo, "A10***";
		$genes1{"3*00201"} = "A10***";
		$genes1{"2*01201"} = "A10***";
		print OUT_GENES "3*00201:".$genes{"3*00201"}.",2*01201:".$genes{"2*01201"}."\t";
	}
	if  (exists $genes{"3*01701"} and exists $genes{"3*03301N"} and exists $genes{"2*01801"}){
		push @haplo, "A11";
		$genes1{"3*01701"} = "A11";
		$genes1{"3*03301N"} = "A11";
		$genes1{"2*01801"} = "A11";
		print OUT_GENES "3*01701:".$genes{"3*01701"}.",3*03301N:".$genes{"3*03301N"}.",2*01801:".$genes{"2*01801"}."\t";
	}
	elsif (exists $genes{"3*01701"} and exists $genes{"2*01801"}){
		push @haplo, "A11*";
		$genes1{"3*01701"} = "A11*";
		$genes1{"2*01801"} = "A11*";
		print OUT_GENES "3*01701:".$genes{"3*01701"}.",2*01801:".$genes{"2*01801"}."\t";
	}
	if  (exists $genes{"1*01901"} and exists $genes{"2*00801"}){
		push @haplo, "A12(w12B)";
		$genes1{"1*01901"} = "A12(w12B)";
		$genes1{"2*00801"} = "A12(w12B)";
		print OUT_GENES "1*01901:".$genes{"1*01901"}.",2*00801:".$genes{"2*00801"}."\t";
	}
	if  (exists $genes{"1*03101"} and exists $genes{"2*03201N"} and exists $genes{"Roslin1.2"}){
		push @haplo, "A13";
		$genes1{"1*03101"} = "A13";
		$genes1{"2*03201N"} = "A13";
		$genes1{"Roslin1.2"} = "A13";
		print OUT_GENES "1*03101:".$genes{"1*03101"}.",2*03201N:".$genes{"2*03201N"}.",Roslin1.2:".$genes{"Roslin1.2"}."\t";
	}
	elsif (exists $genes{"1*03101"} and exists $genes{"2*03201N"}){
		push @haplo, "A13*";
		$genes1{"1*03101"} = "A13*";
		$genes1{"2*03201N"} = "A13*";
		print OUT_GENES "1*03101:".$genes{"1*03101"}.",2*03201N:".$genes{"2*03201N"}."\t";
	}
	if  (exists $genes{"1*02301"} and exists $genes{"4*02401"} and exists $genes{"2*02501"} and exists $genes{"6*04001"}){
		push @haplo, "A14";
		$genes1{"1*02301"} = "A14";
		$genes1{"4*02401"} = "A14";
		$genes1{"2*02501"} = "A14";
		$genes1{"6*04001"} = "A14";
		print OUT_GENES "1*02301:".$genes{"1*02301"}.",4*02401:".$genes{"4*02401"}.",2*02501:".$genes{"2*02501"}.",6*04001:".$genes{"6*04001"}."\t";
	}
	if  (exists $genes{"1*00901"} and exists $genes{"4*02401"} and exists $genes{"2*02501"}){
		push @haplo, "A15";
		$genes1{"1*00901"} = "A15";
		$genes1{"4*02401"} = "A15";
		$genes1{"2*02501"} = "A15";
		print OUT_GENES "1*00901:".$genes{"1*00901"}.",4*02401:".$genes{"4*02401"}.",2*02501:".$genes{"2*02501"}."\t";
	}
	if  (exists $genes{"1*00902"} and exists $genes{"4*02401"} and exists $genes{"2*02501"} and exists $genes{"6*04001"}){
		push @haplo, "A15v";
		$genes1{"1*00902"} = "A15v";
		$genes1{"4*02401"} = "A15v";
		$genes1{"2*02501"} = "A15v";
		$genes1{"6*04001"} = "A15v";
		print OUT_GENES "1*00902:".$genes{"1*00902"}.",4*02401:".$genes{"4*02401"}.",2*02501:".$genes{"2*02501"}.",6*04001:".$genes{"6*04001"}."\t";
	}
	if  (exists $genes{"6*01402"} and exists $genes{"2*01601"}){
		push @haplo, "A19";
		$genes1{"6*01402"} = "A19";
		$genes1{"2*01601"} = "A19";
		print OUT_GENES "6*01402:".$genes{"6*01402"}.",2*01601:".$genes{"2*01601"}."\t";
	}
	elsif (!exists $genes{"Roslin1.9"} and !exists $genes{"Roslin1.10"} and !exists $genes{"Roslin3.12"} and !exists $genes{"Roslin3.36"} and !exists $genes{"Roslin3.37"} and exists $genes{"2*01601"}){
		push @haplo, "BF5";
		$genes1{"BF5"} = "BF5";
		$genes1{"2*01601"} = "BF5";
		print OUT_GENES "2*01601:".$genes{"2*01601"}."\t";
	}
	if  (exists $genes{"3*02702"} and exists $genes{"2*02601"}){
		push @haplo, "A20(v2)";
		$genes1{"3*02702"} = "A20(v2)";
		$genes1{"2*02601"} = "A20(v2)";
		print OUT_GENES "3*02702:".$genes{"3*02702"}.",2*02601:".$genes{"2*02601"}."\t";
	}
	if  (exists $genes{"1*02101"} and exists $genes{"2*02201"}){
		push @haplo, "A31";
		$genes1{"1*02101"} = "A31";
		$genes1{"2*02201"} = "A31";
		print OUT_GENES "1*02101:".$genes{"1*02101"}.",2*02201:".$genes{"2*02201"}."\t";
	}
	if  (exists $genes{"2*05401"} and exists $genes{"Roslin1.3"} and exists $genes{"Roslin1.4"}){
		push @haplo, "BF1";
		$genes1{"2*05401"} = "BF1";
		$genes1{"Roslin1.3"} = "BF1";
		$genes1{"Roslin1.4"} = "BF1";
		print OUT_GENES "2*05401:".$genes{"2*05401"}.",Roslin1.3:".$genes{"Roslin1.3"}.",Roslin1.4:".$genes{"Roslin1.4"}."\t";
	}
	if  (exists $genes{"3*03601"} and exists $genes{"3*03701"} and exists $genes{"Roslin1.5"}){
		push @haplo, "H5(New5)";
		$genes1{"3*03601"} = "H5(New5)";
		$genes1{"3*03701"} = "H5(New5)";
		$genes1{"Roslin1.5"} = "H5(New5)";
		print OUT_GENES "3*03601:".$genes{"3*03601"}.",3*03701:".$genes{"3*03701"}.",Roslin1.5:".$genes{"Roslin1.5"}."\t";
	}
	if  (exists $genes{"Roslin1.6"} and exists $genes{"4*02401"}){
		push @haplo, "HP1.1";
		$genes1{"Roslin1.6"} = "HP1.1";
		$genes1{"4*02401"} = "HP1.1";
		print OUT_GENES "Roslin1.6:".$genes{"Roslin1.6"}.",4*02401:".$genes{"4*02401"}."\t";
	}
	if  (exists $genes{"Roslin1.7"} and exists $genes{"Roslin1.8"} and exists $genes{"3*02702"} and exists $genes{"6*04001"}){
		push @haplo, "HP1.2";
		$genes1{"Roslin1.7"} = "HP1.2";
		$genes1{"Roslin1.8"} = "HP1.2";
		$genes1{"3*02702"} = "HP1.2";
		$genes1{"6*04001"} = "HP1.2";
		print OUT_GENES "Roslin1.7:".$genes{"Roslin1.7"}.",Roslin1.8:".$genes{"Roslin1.8"}.",3*02702:".$genes{"3*02702"}.",6*04001:".$genes{"6*04001"}."\t";
	}
	elsif (exists $genes{"Roslin1.7"} and exists $genes{"Roslin1.8"} and exists $genes{"3*02702"}){
		push @haplo, "HP1.2*";
		$genes1{"Roslin1.7"} = "HP1.2*";
		$genes1{"Roslin1.8"} = "HP1.2*";
		$genes1{"3*02702"} = "HP1.2*";
		print OUT_GENES "Roslin1.7:".$genes{"Roslin1.7"}.",Roslin1.8:".$genes{"Roslin1.8"}.",3*02702:".$genes{"3*02702"}."\t";
	}
	if  (exists $genes{"Roslin1.9"} and exists $genes{"Roslin1.10"} and exists $genes{"2*01602"} and exists $genes{"Roslin(For1Rev2).1.1"}){
		push @haplo, "HP1.3#";
		$genes1{"Roslin1.9"} = "HP1.3#";
		$genes1{"Roslin1.10"} = "HP1.3#";
		$genes1{"2*01602"} = "HP1.3#";
		$genes1{"Roslin(For1Rev2).1.1"} = "HP1.3#";
		print OUT_GENES "Roslin1.9:".$genes{"Roslin1.9"}.",Roslin1.10:".$genes{"Roslin1.10"}.",2*01602:".$genes{"2*01602"}.",Roslin(For1Rev2).1.1:".$genes{"Roslin(For1Rev2).1.1"}."\t";
	}
	elsif (exists $genes{"Roslin1.9"} and exists $genes{"Roslin1.10"} and exists $genes{"2*01602"}){
		push @haplo, "HP1.3";
		$genes1{"Roslin1.9"} = "HP1.3";
		$genes1{"Roslin1.10"} = "HP1.3";
		$genes1{"2*01602"} = "HP1.3";
		print OUT_GENES "Roslin1.9:".$genes{"Roslin1.9"}.",Roslin1.10:".$genes{"Roslin1.10"}.",2*01602:".$genes{"2*01602"}."\t";
	}
	if  (exists $genes{"Roslin2.1"} and exists $genes{"Roslin2.2"} and exists $genes{"3*00402"} and exists $genes{"Roslin(For3Rev1).2.1"}){
		push @haplo, "HP2.1#";
		$genes1{"Roslin2.1"} = "HP2.1#";
		$genes1{"Roslin2.2"} = "HP2.1#";
		$genes1{"3*00402"} = "HP2.1#";
		$genes1{"Roslin(For3Rev1).2.1"} = "HP2.1#";
		print OUT_GENES "Roslin2.1:".$genes{"Roslin2.1"}.",Roslin2.2:".$genes{"Roslin2.2"}.",3*00402:".$genes{"3*00402"}.",Roslin(For3Rev1).2.1:".$genes{"Roslin(For3Rev1).2.1"}."\t";
	}
	elsif (exists $genes{"Roslin2.1"} and exists $genes{"Roslin2.2"} and exists $genes{"3*00402"}){
		push @haplo, "HP2.1";
		$genes1{"Roslin2.1"} = "HP2.1";
		$genes1{"Roslin2.2"} = "HP2.1";
		$genes1{"3*00402"} = "HP2.1";
		print OUT_GENES "Roslin2.1:".$genes{"Roslin2.1"}.",Roslin2.2:".$genes{"Roslin2.2"}.",3*00402:".$genes{"3*00402"}."\t";
	}
	if  (exists $genes{"Roslin2.21"} and exists $genes{"Roslin2.22"}){
		push @haplo, "HP2.10";
		$genes1{"Roslin2.21"} = "HP2.10";
		$genes1{"Roslin2.22"} = "HP2.10";
		print OUT_GENES "Roslin2.21:".$genes{"Roslin2.21"}.",Roslin2.22:".$genes{"Roslin2.22"}."\t";
	}
	if  (exists $genes{"Roslin2.23"} and exists $genes{"Roslin2.19"} and exists $genes{"Roslin2.24"}){
		push @haplo, "HP2.11";
		$genes1{"Roslin2.23"} = "HP2.11";
		$genes1{"Roslin2.19"} = "HP2.11";
		$genes1{"Roslin2.24"} = "HP2.11";
		print OUT_GENES "Roslin2.23:".$genes{"Roslin2.23"}.",Roslin2.19:".$genes{"Roslin2.19"}.",Roslin2.24:".$genes{"Roslin2.24"}."\t";
	}
	elsif (exists $genes{"Roslin2.23"} and exists $genes{"Roslin2.24"}){
		push @haplo, "HP2.11*";
		$genes1{"Roslin2.23"} = "HP2.11*";
		$genes1{"Roslin2.24"} = "HP2.11*";
		print OUT_GENES "Roslin2.23:".$genes{"Roslin2.23"}.",Roslin2.24:".$genes{"Roslin2.24"}."\t";
	}
	if  (exists $genes{"Roslin2.25"} and exists $genes{"Roslin2.27"} and exists $genes{"Roslin2.26"}){
		push @haplo, "HP2.12";
		$genes1{"Roslin2.25"} = "HP2.12";
		$genes1{"Roslin2.27"} = "HP2.12";
		$genes1{"Roslin2.26"} = "HP2.12";
		print OUT_GENES "Roslin2.25:".$genes{"Roslin2.25"}.",Roslin2.27:".$genes{"Roslin2.27"}.",Roslin2.26:".$genes{"Roslin2.26"}."\t";
	}
	elsif (exists $genes{"Roslin2.25"} and exists $genes{"Roslin2.26"}){
		push @haplo, "HP2.12*";
		$genes1{"Roslin2.25"} = "HP2.12*";
		$genes1{"Roslin2.26"} = "HP2.12*";
		print OUT_GENES "Roslin2.25:".$genes{"Roslin2.25"}.",Roslin2.26:".$genes{"Roslin2.26"}."\t";
	}
	if  (exists $genes{"Roslin2.28"} and exists $genes{"4*02401"} and exists $genes{"2*02501"} and exists $genes{"6*04001"}){
		push @haplo, "HP2.13";
		$genes1{"Roslin2.28"} = "HP2.13";
		$genes1{"4*02401"} = "HP2.13";
		$genes1{"2*02501"} = "HP2.13";
		$genes1{"6*04001"} = "HP2.13";
		print OUT_GENES "Roslin2.28:".$genes{"Roslin2.28"}.",4*02401:".$genes{"4*02401"}.",2*02501:".$genes{"2*02501"}.",6*04001:".$genes{"6*04001"}."\t";
	}
	elsif (exists $genes{"Roslin2.28"} and exists $genes{"4*02401"} and exists $genes{"2*02501"}){
		push @haplo, "HP2.13*";
		$genes1{"Roslin2.28"} = "HP2.13*";
		$genes1{"4*02401"} = "HP2.13*";
		$genes1{"2*02501"} = "HP2.13*";
		print OUT_GENES "Roslin2.28:".$genes{"Roslin2.28"}.",4*02401:".$genes{"4*02401"}.",2*02501:".$genes{"2*02501"}."\t";
	}
	if  (exists $genes{"Roslin2.29"} and exists $genes{"Roslin2.30"} and exists $genes{"2*06001"}){
		push @haplo, "HP2.14";
		$genes1{"Roslin2.29"} = "HP2.14";
		$genes1{"Roslin2.30"} = "HP2.14";
		$genes1{"2*06001"} = "HP2.14";
		print OUT_GENES "Roslin2.29:".$genes{"Roslin2.29"}.",Roslin2.30:".$genes{"Roslin2.30"}.",2*06001:".$genes{"2*06001"}."\t";
	}
	if  (exists $genes{"Roslin2.31"} and exists $genes{"3*00102"} and exists $genes{"RoslinSZ.2.1"}){
		push @haplo, "HP2.15";
		$genes1{"Roslin2.31"} = "HP2.15";
		$genes1{"3*00102"} = "HP2.15";
		$genes1{"RoslinSZ.2.1"} = "HP2.15";
		print OUT_GENES "Roslin2.31:".$genes{"Roslin2.31"}.",3*00102:".$genes{"3*00102"}.",RoslinSZ.2.1:".$genes{"RoslinSZ.2.1"}."\t";
	}
	elsif (exists $genes{"Roslin2.31"} and exists $genes{"3*00102"}){
		push @haplo, "HP2.15*";
		$genes1{"Roslin2.31"} = "HP2.15*";
		$genes1{"3*00102"} = "HP2.15*";
		print OUT_GENES "Roslin2.31:".$genes{"Roslin2.31"}.",3*00102:".$genes{"3*00102"}."\t";
	}
	if  (exists $genes{"Roslin2.32"} and exists $genes{"3*05901"}){
		push @haplo, "HP2.16";
		$genes1{"Roslin2.32"} = "HP2.16";
		$genes1{"3*05901"} = "HP2.16";
		print OUT_GENES "Roslin2.32:".$genes{"Roslin2.32"}.",3*05901:".$genes{"3*05901"}."\t";
	}
	if  (exists $genes{"Roslin2.33"} and exists $genes{"Roslin2.12"} and exists $genes{"3*00402"}){
		push @haplo, "HP2.17";
		$genes1{"Roslin2.33"} = "HP2.17";
		$genes1{"Roslin2.12"} = "HP2.17";
		$genes1{"3*00402"} = "HP2.17";
		print OUT_GENES "Roslin2.33:".$genes{"Roslin2.33"}.",Roslin2.12:".$genes{"Roslin2.12"}.",3*00402:".$genes{"3*00402"}."\t";
	}
	if  (exists $genes{"Roslin2.34"} and exists $genes{"Roslin2.36"} and exists $genes{"Roslin2.35"}){
		push @haplo, "HP2.18";
		$genes1{"Roslin2.34"} = "HP2.18";
		$genes1{"Roslin2.36"} = "HP2.18";
		$genes1{"Roslin2.35"} = "HP2.18";
		print OUT_GENES "Roslin2.34:".$genes{"Roslin2.34"}.",Roslin2.36:".$genes{"Roslin2.36"}.",Roslin2.35:".$genes{"Roslin2.35"}."\t";
	}
	if  (exists $genes{"Roslin2.47"} and exists $genes{"Roslin2.48"} and exists $genes{"Roslin2.49"} and exists $genes{"Roslin(For3Rev1).2.4"}){
		push @haplo, "HP2.19#";
		$genes1{"Roslin2.47"} = "HP2.19#";
		$genes1{"Roslin2.48"} = "HP2.19#";
		$genes1{"Roslin2.49"} = "HP2.19#";
		$genes1{"Roslin(For3Rev1).2.4"} = "HP2.19#";
		print OUT_GENES "Roslin2.47:".$genes{"Roslin2.47"}.",Roslin2.48:".$genes{"Roslin2.48"}.",Roslin2.49:".$genes{"Roslin2.49"}.",Roslin(For3Rev1).2.4:".$genes{"Roslin(For3Rev1).2.4"}."\t";
	}
	elsif (exists $genes{"Roslin2.47"} and exists $genes{"Roslin2.48"} and exists $genes{"Roslin2.49"}){
		push @haplo, "HP2.19";
		$genes1{"Roslin2.47"} = "HP2.19";
		$genes1{"Roslin2.48"} = "HP2.19";
		$genes1{"Roslin2.49"} = "HP2.19";
		print OUT_GENES "Roslin2.47:".$genes{"Roslin2.47"}.",Roslin2.48:".$genes{"Roslin2.48"}.",Roslin2.49:".$genes{"Roslin2.49"}."\t";
	}
	if  (exists $genes{"Roslin2.3"} and exists $genes{"Roslin2.4"}){
		push @haplo, "HP2.2";
		$genes1{"Roslin2.3"} = "HP2.2";
		$genes1{"Roslin2.4"} = "HP2.2";
		print OUT_GENES "Roslin2.3:".$genes{"Roslin2.3"}.",Roslin2.4:".$genes{"Roslin2.4"}."\t";
	}
	if  (exists $genes{"Roslin2.37"} and exists $genes{"Roslin2.38"}){
		push @haplo, "HP2.20";
		$genes1{"Roslin2.37"} = "HP2.20";
		$genes1{"Roslin2.38"} = "HP2.20";
		print OUT_GENES "Roslin2.37:".$genes{"Roslin2.37"}.",Roslin2.38:".$genes{"Roslin2.38"}."\t";
	}
	if  (exists $genes{"Roslin2.39"} and exists $genes{"Roslin2.40"} and exists $genes{"Roslin2.41"}){
		push @haplo, "HP2.21";
		$genes1{"Roslin2.39"} = "HP2.21";
		$genes1{"Roslin2.40"} = "HP2.21";
		$genes1{"Roslin2.41"} = "HP2.21";
		print OUT_GENES "Roslin2.39:".$genes{"Roslin2.39"}.",Roslin2.40:".$genes{"Roslin2.40"}.",Roslin2.41:".$genes{"Roslin2.41"}."\t";
	}
	if (exists $genes{"Roslin2.42"} and exists $genes{"RoslinSZ.3.2"}){
		push @haplo, "HP2.22";
		$genes1{"Roslin2.42"} = "HP2.22";
		$genes1{"RoslinSZ.3.2"} = "HP2.22";
		print OUT_GENES "Roslin2.42:".$genes{"Roslin2.42"}.",RoslinSZ.3.2:".$genes{"RoslinSZ.3.2"}."\t";
	}
	elsif (exists $genes{"Roslin2.42"}){
		push @haplo, "HP2.22*";
		$genes1{"Roslin2.42"} = "HP2.22*";
		print OUT_GENES "Roslin2.42:".$genes{"Roslin2.42"}."\t";
	}
	if  (exists $genes{"Roslin2.39"} and exists $genes{"Roslin2.40"} and exists $genes{"Roslin2.43"}){
		push @haplo, "HP2.24";
		$genes1{"Roslin2.39"} = "HP2.24";
		$genes1{"Roslin2.40"} = "HP2.24";
		$genes1{"Roslin2.43"} = "HP2.24";
		print OUT_GENES "Roslin2.39:".$genes{"Roslin2.39"}.",Roslin2.40:".$genes{"Roslin2.40"}.",Roslin2.43:".$genes{"Roslin2.43"}."\t";
	}
	if  (exists $genes{"Roslin2.39"} and exists $genes{"Roslin2.44"} and exists $genes{"Roslin2.43"}){
		push @haplo, "HP2.25";
		$genes1{"Roslin2.39"} = "HP2.25";
		$genes1{"Roslin2.44"} = "HP2.25";
		$genes1{"Roslin2.43"} = "HP2.25";
		print OUT_GENES "Roslin2.39:".$genes{"Roslin2.39"}.",Roslin2.44:".$genes{"Roslin2.44"}.",Roslin2.43:".$genes{"Roslin2.43"}."\t";
	}
	if  (exists $genes{"Roslin2.45"} and exists $genes{"Roslin(For3Rev1).2.2"}){
		push @haplo, "HP2.26#";
		$genes1{"Roslin2.45"} = "HP2.26#";
		$genes1{"Roslin(For3Rev1).2.2"} = "HP2.26#";
		print OUT_GENES "Roslin2.45:".$genes{"Roslin2.45"}.",Roslin(For3Rev1).2.2:".$genes{"Roslin(For3Rev1).2.2"}."\t";
	}
	elsif (exists $genes{"Roslin2.45"}){
		push @haplo, "HP2.26";
		$genes1{"Roslin2.45"} = "HP2.26";
		print OUT_GENES "Roslin2.45:".$genes{"Roslin2.45"}."\t";
	}
	if  (exists $genes{"Roslin2.46"} and exists $genes{"Roslin(For3Rev1).2.3"}){
		push @haplo, "HP2.27#";
		$genes1{"Roslin2.46"} = "HP2.27#";
		$genes1{"Roslin(For3Rev1).2.3"} = "HP2.27#";
		print OUT_GENES "Roslin2.46:".$genes{"Roslin2.46"}.",Roslin(For3Rev1).2.3:".$genes{"Roslin(For3Rev1).2.3"}."\t";
	}
	elsif (exists $genes{"Roslin2.46"}){
		push @haplo, "HP2.27";
		$genes1{"Roslin2.46"} = "HP2.27";
		print OUT_GENES "Roslin2.46:".$genes{"Roslin2.46"}."\t";
	}
	if  (exists $genes{"1*02801"} and exists $genes{"4*02401"}){
		push @haplo, "HP2.23";
		$genes1{"1*02801"} = "HP2.23";
		$genes1{"4*02401"} = "HP2.23";
		print OUT_GENES "1*02801:".$genes{"1*02801"}.",4*02401:".$genes{"4*02401"}."\t";
	}
	if  (exists $genes{"Roslin2.5"} and exists $genes{"Roslin2.6"} and exists $genes{"Roslin2.7"}){
		push @haplo, "HP2.3";
		$genes1{"Roslin2.5"} = "HP2.3";
		$genes1{"Roslin2.6"} = "HP2.3";
		$genes1{"Roslin2.7"} = "HP2.3";
		print OUT_GENES "Roslin2.5:".$genes{"Roslin2.5"}.",Roslin2.6:".$genes{"Roslin2.6"}.",Roslin2.7:".$genes{"Roslin2.7"}."\t";
	}
	if  (exists $genes{"Roslin2.8"} and exists $genes{"Roslin2.9"} and exists $genes{"Roslin2.10"}){
		push @haplo, "HP2.4";
		$genes1{"Roslin2.8"} = "HP2.4";
		$genes1{"Roslin2.9"} = "HP2.4";
		$genes1{"Roslin2.10"} = "HP2.4";
		print OUT_GENES "Roslin2.8:".$genes{"Roslin2.8"}.",Roslin2.9:".$genes{"Roslin2.9"}.",Roslin2.10:".$genes{"Roslin2.10"}."\t";
	}
	if  (exists $genes{"Roslin2.11"} and exists $genes{"Roslin2.12"} and exists $genes{"3*00402"}){
		push @haplo, "HP2.5";
		$genes1{"Roslin2.11"} = "HP2.5";
		$genes1{"Roslin2.12"} = "HP2.5";
		$genes1{"3*00402"} = "HP2.5";
		print OUT_GENES "Roslin2.11:".$genes{"Roslin2.11"}.",Roslin2.12:".$genes{"Roslin2.12"}.",3*00402:".$genes{"3*00402"}."\t";
	}
	if  (exists $genes{"Roslin2.13"} and exists $genes{"Roslin(For3Rev1).2.2"}){
		push @haplo, "HP2.6#";
		$genes1{"Roslin2.13"} = "HP2.6#";
		$genes1{"Roslin(For3Rev1).2.2"} = "HP2.6#";
		print OUT_GENES "Roslin2.13:".$genes{"Roslin2.13"}.",Roslin(For3Rev1).2.2:".$genes{"Roslin(For3Rev1).2.2"}."\t";
	}
	elsif (exists $genes{"Roslin2.13"}){
		push @haplo, "HP2.6";
		$genes1{"Roslin2.13"} = "HP2.6";
		print OUT_GENES "Roslin2.13:".$genes{"Roslin2.13"}."\t";
	}
	if  (exists $genes{"Roslin2.14"} and exists $genes{"Roslin2.15"} and exists $genes{"Roslin2.16"}){
		push @haplo, "HP2.7";
		$genes1{"Roslin2.14"} = "HP2.7";
		$genes1{"Roslin2.15"} = "HP2.7";
		$genes1{"Roslin2.16"} = "HP2.7";
		print OUT_GENES "Roslin2.14:".$genes{"Roslin2.14"}.",Roslin2.15:".$genes{"Roslin2.15"}.",Roslin2.16:".$genes{"Roslin2.16"}."\t";
	}
	if  (exists $genes{"Roslin2.17"} and exists $genes{"Roslin2.19"} and exists $genes{"Roslin2.18"}){
		push @haplo, "HP2.8";
		$genes1{"Roslin2.17"} = "HP2.8";
		$genes1{"Roslin2.19"} = "HP2.8";
		$genes1{"Roslin2.18"} = "HP2.8";
		print OUT_GENES "Roslin2.17:".$genes{"Roslin2.17"}.",Roslin2.19:".$genes{"Roslin2.19"}.",Roslin2.18:".$genes{"Roslin2.18"}."\t";
	}
	elsif (exists $genes{"Roslin2.17"} and exists $genes{"Roslin2.18"}){
		push @haplo, "HP2.8*";
		$genes1{"Roslin2.17"} = "HP2.8*";
		$genes1{"Roslin2.18"} = "HP2.8*";
		print OUT_GENES "Roslin2.17:".$genes{"Roslin2.17"}.",Roslin2.18:".$genes{"Roslin2.18"}."\t";
	}
	if  (exists $genes{"Roslin2.20"} and exists $genes{"3*00102"} and exists $genes{"RoslinSZ.2.1"}){
		push @haplo, "HP2.9";
		$genes1{"Roslin2.20"} = "HP2.9";
		$genes1{"3*00102"} = "HP2.9";
		$genes1{"RoslinSZ.2.1"} = "HP2.9";
		print OUT_GENES "Roslin2.20:".$genes{"Roslin2.20"}.",3*00102:".$genes{"3*00102"}.",RoslinSZ.2.1:".$genes{"RoslinSZ.2.1"}."\t";
	}
	elsif (exists $genes{"Roslin2.20"} and exists $genes{"3*00102"}){
		push @haplo, "HP2.9*";
		$genes1{"Roslin2.20"} = "HP2.9*";
		$genes1{"3*00102"} = "HP2.9*";
		print OUT_GENES "Roslin2.20:".$genes{"Roslin2.20"}.",3*00102:".$genes{"3*00102"}."\t";
	}
	if  (exists $genes{"Roslin3.1"} and exists $genes{"Roslin3.2"} and exists $genes{"Roslin3.3"}){
		push @haplo, "HP3.1";
		$genes1{"Roslin3.1"} = "HP3.1";
		$genes1{"Roslin3.2"} = "HP3.1";
		$genes1{"Roslin3.3"} = "HP3.1";
		print OUT_GENES "Roslin3.1:".$genes{"Roslin3.1"}.",Roslin3.2:".$genes{"Roslin3.2"}.",Roslin3.3:".$genes{"Roslin3.3"}."\t";
	}
	if  (exists $genes{"Roslin3.15"} and exists $genes{"Roslin3.16"}){
		push @haplo, "HP3.10";
		$genes1{"Roslin3.15"} = "HP3.10";
		$genes1{"Roslin3.16"} = "HP3.10";
		print OUT_GENES "Roslin3.15:".$genes{"Roslin3.15"}.",Roslin3.16:".$genes{"Roslin3.16"}."\t";
	}
	if  (exists $genes{"Roslin3.17"}){
		push @haplo, "HP3.11";
		$genes1{"Roslin3.17"} = "HP3.11";
		print OUT_GENES "Roslin3.17:".$genes{"Roslin3.17"}."\t";
	}
	if  (exists $genes{"Roslin1.4"} and exists $genes{"Roslin3.18"} and exists $genes{"4*06301"}){
		push @haplo, "HP3.12";
		$genes1{"Roslin1.4"} = "HP3.12";
		$genes1{"Roslin3.18"} = "HP3.12";
		$genes1{"4*06301"} = "HP3.12";
		print OUT_GENES "Roslin1.4:".$genes{"Roslin1.4"}.",Roslin3.18:".$genes{"Roslin3.18"}.",4*06301:".$genes{"4*06301"}."\t";
	}
	if  (exists $genes{"Roslin3.19"}){
		push @haplo, "HP3.13";
		$genes1{"Roslin3.19"} = "HP3.13";
		print OUT_GENES "Roslin3.19:".$genes{"Roslin3.19"}."\t";
	}
	if  (exists $genes{"Roslin3.20"} and exists $genes{"Roslin(For3Rev1).3.2"} and exists $genes{"Roslin(For3Rev1).2.4"}){
		push @haplo, "HP3.14##";
		$genes1{"Roslin3.20"} = "HP3.14##";
		$genes1{"Roslin(For3Rev1).3.2"} = "HP3.14##";
		$genes1{"Roslin(For3Rev1).2.4"} = "HP3.14##";
		print OUT_GENES "Roslin3.20:".$genes{"Roslin3.20"}.",Roslin(For3Rev1).3.2:".$genes{"Roslin(For3Rev1).3.2"}.",Roslin(For3Rev1).2.4:".$genes{"Roslin(For3Rev1).2.4"}."\t";
	}
	elsif (exists $genes{"Roslin3.20"} and exists $genes{"RoslinSZ(For1Rev2).3.1"}){
		push @haplo, "HP3.14#";
		$genes1{"Roslin3.20"} = "HP3.14#";
		$genes1{"RoslinSZ(For1Rev2).3.1"} = "HP3.14#";
		print OUT_GENES "Roslin3.20:".$genes{"Roslin3.20"}.",RoslinSZ(For1Rev2).3.1:".$genes{"RoslinSZ(For1Rev2).3.1"}."\t";
	}
	if  (exists $genes{"Roslin3.5"} and exists $genes{"Roslin3.21"} and exists $genes{"3*03701"}){
		push @haplo, "HP3.15";
		$genes1{"Roslin3.5"} = "HP3.15";
		$genes1{"Roslin3.21"} = "HP3.15";
		$genes1{"3*03701"} = "HP3.15";
		print OUT_GENES "Roslin3.5:".$genes{"Roslin3.5"}.",Roslin3.21:".$genes{"Roslin3.21"}.",3*03701:".$genes{"3*03701"}."\t";
	}
	if  (exists $genes{"Roslin2.34"} and exists $genes{"Roslin3.22"} and exists $genes{"Roslin2.35"}){
		push @haplo, "HP3.16";
		$genes1{"Roslin2.34"} = "HP3.16";
		$genes1{"Roslin3.22"} = "HP3.16";
		$genes1{"Roslin2.35"} = "HP3.16";
		print OUT_GENES "Roslin2.34:".$genes{"Roslin2.34"}.",Roslin3.22:".$genes{"Roslin3.22"}.",Roslin2.35:".$genes{"Roslin2.35"}."\t";
	}
	if  (exists $genes{"Roslin3.24"} and exists $genes{"Roslin2.39"} and exists $genes{"Roslin2.40"}){
		push @haplo, "HP3.17";
		$genes1{"Roslin3.24"} = "HP3.17";
		$genes1{"Roslin2.39"} = "HP3.17";
		$genes1{"Roslin2.40"} = "HP3.17";
		print OUT_GENES "Roslin3.24:".$genes{"Roslin3.24"}.",Roslin2.39:".$genes{"Roslin2.39"}.",Roslin2.40:".$genes{"Roslin2.40"}."\t";
	}
	if  (exists $genes{"Roslin3.25"} and exists $genes{"Roslin3.26"} and exists $genes{"Roslin(For3Rev1).2.3"}){
		push @haplo, "HP3.18#";
		$genes1{"Roslin3.25"} = "HP3.18#";
		$genes1{"Roslin3.26"} = "HP3.18#";
		$genes1{"Roslin(For3Rev1).2.3"} = "HP3.18#";
		print OUT_GENES "Roslin3.25:".$genes{"Roslin3.25"}.",Roslin3.26:".$genes{"Roslin3.26"}.",Roslin(For3Rev1).2.3:".$genes{"Roslin(For3Rev1).2.3"}."\t";
		}
		elsif (exists $genes{"Roslin3.25"} and exists $genes{"Roslin3.26"}){
		push @haplo, "HP3.18";
		$genes1{"Roslin3.25"} = "HP3.18";
		$genes1{"Roslin3.26"} = "HP3.18";
		print OUT_GENES "Roslin3.25:".$genes{"Roslin3.25"}.",Roslin3.26:".$genes{"Roslin3.26"}."\t";
	}
	if  (exists $genes{"Roslin3.29"} and exists $genes{"Roslin3.27"} and exists $genes{"Roslin3.28"}){
		push @haplo, "HP3.19";
		$genes1{"Roslin3.29"} = "HP3.19";
		$genes1{"Roslin3.27"} = "HP3.19";
		$genes1{"Roslin3.28"} = "HP3.19";
		print OUT_GENES "Roslin3.29:".$genes{"Roslin3.29"}.",Roslin3.27:".$genes{"Roslin3.27"}.",Roslin3.28:".$genes{"Roslin3.28"}."\t";
	}
	if  (exists $genes{"Roslin3.5"} and exists $genes{"Roslin3.4"} and exists $genes{"3*03701"}){
		push @haplo, "HP3.2";
		$genes1{"Roslin3.5"} = "HP3.2";
		$genes1{"Roslin3.4"} = "HP3.2";
		$genes1{"3*03701"} = "HP3.2";
		print OUT_GENES "Roslin3.5:".$genes{"Roslin3.5"}.",Roslin3.4:".$genes{"Roslin3.4"}.",3*03701:".$genes{"3*03701"}."\t";
	}
	if  (exists $genes{"Roslin2.35"} and exists $genes{"Roslin3.30"} and exists $genes{"Roslin3.32"} and exists $genes{"Roslin3.31"} and exists $genes{"Roslin1.1"}){
		push @haplo, "HP3.20";
		$genes1{"Roslin2.35"} = "HP3.20";
		$genes1{"Roslin3.30"} = "HP3.20";
		$genes1{"Roslin3.32"} = "HP3.20";
		$genes1{"Roslin3.31"} = "HP3.20";
		$genes1{"Roslin1.1"} = "HP3.20";
		print OUT_GENES "Roslin2.35:".$genes{"Roslin2.35"}.",Roslin3.30:".$genes{"Roslin3.30"}.",Roslin3.32:".$genes{"Roslin3.32"}.",Roslin3.31:".$genes{"Roslin3.31"}.",Roslin1.1:".$genes{"Roslin1.1"}."\t";
	}
	elsif (exists $genes{"Roslin2.35"} and exists $genes{"Roslin3.30"} and exists $genes{"Roslin3.32"} and exists $genes{"Roslin3.31"}){
		push @haplo, "HP3.20*";
		$genes1{"Roslin2.35"} = "HP3.20*";
		$genes1{"Roslin3.30"} = "HP3.20*";
		$genes1{"Roslin3.32"} = "HP3.20*";
		$genes1{"Roslin3.31"} = "HP3.20*";
		print OUT_GENES "Roslin2.35:".$genes{"Roslin2.35"}.",Roslin3.30:".$genes{"Roslin3.30"}.",Roslin3.32:".$genes{"Roslin3.32"}.",Roslin3.31:".$genes{"Roslin3.31"}."\t";
	}
	if  (exists $genes{"Roslin3.33"} and exists $genes{"1*02101"}){
		push @haplo, "HP3.21";
		$genes1{"Roslin3.33"} = "HP3.21";
		$genes1{"1*02101"} = "HP3.21";
		print OUT_GENES "Roslin3.33:".$genes{"Roslin3.33"}.",1*02101:".$genes{"1*02101"}."\t";
	}
	if  (exists $genes{"Roslin3.34"} and exists $genes{"2*04701"} and exists $genes{"RoslinSZ.3.1"} and exists $genes{"Roslin(For3Rev1).3.1"}){
		push @haplo, "HP3.22#";
		$genes1{"Roslin3.34"} = "HP3.22#";
		$genes1{"2*04701"} = "HP3.22#";
		$genes1{"RoslinSZ.3.1"} = "HP3.22#";
		$genes1{"Roslin(For3Rev1).3.1"} = "HP3.22#";
		print OUT_GENES "Roslin3.34:".$genes{"Roslin3.34"}.",2*04701:".$genes{"2*04701"}.",RoslinSZ.3.1:".$genes{"RoslinSZ.3.1"}.",Roslin(For3Rev1).3.1:".$genes{"Roslin(For3Rev1).3.1"}."\t";
	}
	elsif (exists $genes{"Roslin3.34"} and exists $genes{"2*04701"} and exists $genes{"RoslinSZ.3.1"}){
		push @haplo, "HP3.22";
		$genes1{"Roslin3.34"} = "HP3.22";
		$genes1{"2*04701"} = "HP3.22";
		$genes1{"RoslinSZ.3.1"} = "HP3.22";
		print OUT_GENES "Roslin3.34:".$genes{"Roslin3.34"}.",2*04701:".$genes{"2*04701"}.",RoslinSZ.3.1:".$genes{"RoslinSZ.3.1"}."\t";
	}
	if  (exists $genes{"Roslin3.35"} and exists $genes{"Roslin(For1Rev2).3.2"}){
		push @haplo, "HP3.23#";
		$genes1{"Roslin3.35"} = "HP3.23#";
		$genes1{"Roslin(For1Rev2).3.2"} = "HP3.23#";
		print OUT_GENES "Roslin3.35:".$genes{"Roslin3.35"}.",Roslin(For1Rev2).3.2:".$genes{"Roslin(For1Rev2).3.2"}."\t";
	}
	elsif (exists $genes{"Roslin3.35"}){
		push @haplo, "HP3.23";
		$genes1{"Roslin3.35"} = "HP3.23";
		print OUT_GENES "Roslin3.35:".$genes{"Roslin3.35"}."\t";
	}
	if (exists $genes{"Roslin3.36"} and exists $genes{"Roslin3.37"} and exists $genes{"2*01601"} and exists $genes{"Roslin2.27"}){
		push @haplo, "HP3.24#";
		$genes1{"Roslin3.36"} = "HP3.24#";
		$genes1{"Roslin3.37"} = "HP3.24#";
		$genes1{"2*01601"} = "HP3.24#";
		$genes1{"Roslin2.27"} = "HP3.24#";
		print OUT_GENES "Roslin3.36:".$genes{"Roslin3.36"}.",Roslin3.37:".$genes{"Roslin3.37"}.",2*01601:".$genes{"2*01601"}.",Roslin2.27:".$genes{"Roslin2.27"}."\t";
	}
	elsif (exists $genes{"Roslin3.36"} and exists $genes{"Roslin3.37"} and exists $genes{"2*01601"}){
		push @haplo, "HP3.24";
		$genes1{"Roslin3.36"} = "HP3.24";
		$genes1{"Roslin3.37"} = "HP3.24";
		$genes1{"2*01601"} = "HP3.24";
		print OUT_GENES "Roslin3.36:".$genes{"Roslin3.36"}.",Roslin3.37:".$genes{"Roslin3.37"}.",2*01601:".$genes{"2*01601"}."\t";
	}
	if  (exists $genes{"Roslin3.38"} and exists $genes{"Roslin2.29"} and exists $genes{"2*06001"}){
		push @haplo, "HP3.25";
		$genes1{"Roslin3.38"} = "HP3.25";
		$genes1{"Roslin2.29"} = "HP3.25";
		$genes1{"2*06001"} = "HP3.25";
		print OUT_GENES "Roslin3.38:".$genes{"Roslin3.38"}.",Roslin2.29:".$genes{"Roslin2.29"}.",2*06001:".$genes{"2*06001"}."\t";
	}
	if  (exists $genes{"Roslin3.39"} and exists $genes{"Roslin3.23"} and exists $genes{"Roslin3.40"} and exists $genes{"Roslin3.41"}){
		push @haplo, "HP3.26";
		$genes1{"Roslin3.39"} = "HP3.26";
		$genes1{"Roslin3.23"} = "HP3.26";
		$genes1{"Roslin3.40"} = "HP3.26";
		$genes1{"Roslin3.41"} = "HP3.26";
		print OUT_GENES "Roslin3.39:".$genes{"Roslin3.39"}.",Roslin3.23:".$genes{"Roslin3.23"}.",Roslin3.40:".$genes{"Roslin3.40"}.",Roslin3.41:".$genes{"Roslin3.41"}."\t";
	}
	elsif (exists $genes{"Roslin3.39"} and exists $genes{"Roslin3.23"} and exists $genes{"Roslin3.40"}){
		push @haplo, "HP3.26*";
		$genes1{"Roslin3.39"} = "HP3.26*";
		$genes1{"Roslin3.23"} = "HP3.26*";
		$genes1{"Roslin3.40"} = "HP3.26*";
		print OUT_GENES "Roslin3.39:".$genes{"Roslin3.39"}.",Roslin3.23:".$genes{"Roslin3.23"}.",Roslin3.40:".$genes{"Roslin3.40"}."\t";
	}
	if  (exists $genes{"Roslin3.42"} and exists $genes{"Roslin3.43"} and exists $genes{"Roslin(For1Rev2).3.1"}){
		push @haplo, "HP3.27#";
		$genes1{"Roslin3.42"} = "HP3.27#";
		$genes1{"Roslin3.43"} = "HP3.27#";
		$genes1{"Roslin(For1Rev2).3.1"} = "HP3.27#";
		print OUT_GENES "Roslin3.42:".$genes{"Roslin3.42"}.",Roslin3.43:".$genes{"Roslin3.43"}.",Roslin(For1Rev2).3.1:".$genes{"Roslin(For1Rev2).3.1"}."\t";
	}
	elsif (exists $genes{"Roslin3.42"} and exists $genes{"Roslin3.43"}){
		push @haplo, "HP3.27";
		$genes1{"Roslin3.42"} = "HP3.27";
		$genes1{"Roslin3.43"} = "HP3.27";
		print OUT_GENES "Roslin3.42:".$genes{"Roslin3.42"}.",Roslin3.43:".$genes{"Roslin3.43"}."\t";
	}
	elsif (exists $genes{"Roslin3.42"} and exists $genes{"Roslin(For1Rev2).3.1"}){
		push @haplo, "HP3.27*#";
		$genes1{"Roslin3.42"} = "HP3.27*#";
		$genes1{"Roslin(For1Rev2).3.1"} = "HP3.27*#";
		print OUT_GENES "Roslin3.42:".$genes{"Roslin3.42"}.",Roslin(For1Rev2).3.1:".$genes{"Roslin(For1Rev2).3.1"}."\t";
	}
	elsif (exists $genes{"Roslin3.42"}){
		push @haplo, "HP3.27*";
		$genes1{"Roslin3.42"} = "HP3.27*";
		print OUT_GENES "Roslin3.42:".$genes{"Roslin3.42"}."\t";
	}
	if  (exists $genes{"Roslin3.44"} and exists $genes{"Roslin(For3Rev1).3.3"}){
		push @haplo, "HP3.28#";
		$genes1{"Roslin3.44"} = "HP3.28#";
		$genes1{"Roslin(For3Rev1).3.3"} = "HP3.28#";
		print OUT_GENES "Roslin3.44:".$genes{"Roslin3.44"}.",Roslin(For3Rev1).3.3:".$genes{"Roslin(For3Rev1).3.3"}."\t";
	}
	elsif (exists $genes{"Roslin3.44"}){
		push @haplo, "HP3.28";
		$genes1{"Roslin3.44"} = "HP3.28";
		print OUT_GENES "Roslin3.44:".$genes{"Roslin3.44"}."\t";
	}
	if  (exists $genes{"Roslin3.45"} and exists $genes{"3*03301N"} and exists $genes{"Roslin3.46"}){
		push @haplo, "HP3.29";
		$genes1{"Roslin3.45"} = "HP3.29";
		$genes1{"3*03301N"} = "HP3.29";
		$genes1{"Roslin3.46"} = "HP3.29";
		print OUT_GENES "Roslin3.45:".$genes{"Roslin3.45"}.",3*03301N:".$genes{"3*03301N"}.",Roslin3.46:".$genes{"Roslin3.46"}."\t";
	}
	if  (exists $genes{"Roslin3.6"} and exists $genes{"Roslin2.39"} and exists $genes{"Roslin2.40"}){
		push @haplo, "HP3.3";
		$genes1{"Roslin3.6"} = "HP3.3";
		$genes1{"Roslin2.39"} = "HP3.3";
		$genes1{"Roslin2.40"} = "HP3.3";
		print OUT_GENES "Roslin3.6:".$genes{"Roslin3.6"}.",Roslin2.39:".$genes{"Roslin2.39"}.",Roslin2.40:".$genes{"Roslin2.40"}."\t";
	}
	if  (exists $genes{"Roslin3.5"} and exists $genes{"Roslin3.48"} and exists $genes{"3*03701"}){
		push @haplo, "HP3.30";
		$genes1{"Roslin3.5"} = "HP3.30";
		$genes1{"Roslin3.48"} = "HP3.30";
		$genes1{"3*03701"} = "HP3.30";
		print OUT_GENES "Roslin3.5:".$genes{"Roslin3.5"}.",Roslin3.48:".$genes{"Roslin3.48"}.",3*03701:".$genes{"3*03701"}."\t";
	}
	if  (exists $genes{"Roslin3.7"} and exists $genes{"2*04701"} and exists $genes{"2*03201N"}){
		push @haplo, "HP3.4";
		$genes1{"Roslin3.7"} = "HP3.4";
		$genes1{"2*04701"} = "HP3.4";
		$genes1{"2*03201N"} = "HP3.4";
		print OUT_GENES "Roslin3.7:".$genes{"Roslin3.7"}.",2*04701:".$genes{"2*04701"}.",2*03201N:".$genes{"2*03201N"}."\t";
	}
	if  (exists $genes{"Roslin3.8"} and exists $genes{"Roslin3.9"}){
		push @haplo, "HP3.5";
		$genes1{"Roslin3.8"} = "HP3.5";
		$genes1{"Roslin3.9"} = "HP3.5";
		print OUT_GENES "Roslin3.8:".$genes{"Roslin3.8"}.",Roslin3.9:".$genes{"Roslin3.9"}."\t";
	}
	if  (exists $genes{"Roslin3.10"} and exists $genes{"Roslin3.11"} and exists $genes{"Roslin1.1"}){
		push @haplo, "HP3.6";
		$genes1{"Roslin3.10"} = "HP3.6";
		$genes1{"Roslin3.11"} = "HP3.6";
		$genes1{"Roslin1.1"} = "HP3.6";
		print OUT_GENES "Roslin3.10:".$genes{"Roslin3.10"}.",Roslin3.11:".$genes{"Roslin3.11"}.",Roslin1.1:".$genes{"Roslin1.1"}."\t";
	}
	elsif (exists $genes{"Roslin3.10"} and exists $genes{"Roslin3.11"}){
		push @haplo, "HP3.6*";
		$genes1{"Roslin3.10"} = "HP3.6*";
		$genes1{"Roslin3.11"} = "HP3.6*";
		print OUT_GENES "Roslin3.10:".$genes{"Roslin3.10"}.",Roslin3.11:".$genes{"Roslin3.11"}."\t";
	}
	if  (exists $genes{"Roslin3.12"} and exists $genes{"2*01601"}){
		push @haplo, "HP3.7";
		$genes1{"Roslin3.12"} = "HP3.7";
		$genes1{"2*01601"} = "HP3.7";
		print OUT_GENES "Roslin3.12:".$genes{"Roslin3.12"}.",2*01601:".$genes{"2*01601"}."\t";
	}
	if  (exists $genes{"Roslin3.13"} and exists $genes{"Roslin3.2"}){
		push @haplo, "HP3.8";
		$genes1{"Roslin3.13"} = "HP3.8";
		$genes1{"Roslin3.2"} = "HP3.8";
		print OUT_GENES "Roslin3.13:".$genes{"Roslin3.13"}.",Roslin3.2:".$genes{"Roslin3.2"}."\t";
	}
	if  (exists $genes{"Roslin3.14"} and exists $genes{"Roslin3.1"}){
		push @haplo, "HP3.9";
		$genes1{"Roslin3.14"} = "HP3.9";
		$genes1{"Roslin3.1"} = "HP3.9";
		print OUT_GENES "Roslin3.14:".$genes{"Roslin3.14"}.",Roslin3.1:".$genes{"Roslin3.1"}."\t";
	}
	if  (exists $genes{"Roslin2.1"} and exists $genes{"Roslin3.47"} and exists $genes{"3*00402"} and exists $genes{"Roslin(For3Rev1).2.1"}){
		push @haplo, "HP3.31#";
		$genes1{"Roslin2.1"} = "HP3.31#";
		$genes1{"Roslin3.47"} = "HP3.31#";
		$genes1{"3*00402"} = "HP3.31#";
		$genes1{"Roslin(For3Rev1).2.1"} = "HP3.31#";
		print OUT_GENES "Roslin2.1:".$genes{"Roslin2.1"}.",Roslin3.47:".$genes{"Roslin3.47"}.",3*00402:".$genes{"3*00402"}.",Roslin(For3Rev1).2.1:".$genes{"Roslin(For3Rev1).2.1"}."\t";
	}
	elsif (exists $genes{"Roslin2.1"} and exists $genes{"Roslin3.47"} and exists $genes{"3*00402"}){
		push @haplo, "HP3.31";
		$genes1{"Roslin2.1"} = "HP3.31";
		$genes1{"Roslin3.47"} = "HP3.31";
		$genes1{"3*00402"} = "HP3.31";
		print OUT_GENES "Roslin2.1:".$genes{"Roslin2.1"}.",Roslin3.47:".$genes{"Roslin3.47"}.",3*00402:".$genes{"3*00402"}."\t";
	}
	if (exists $genes{"Roslin2.5"} and exists $genes{"Roslin2.38"} and exists $genes{"Roslin2.7"}){
		push @haplo, "HP3.32";
		$genes1{"Roslin2.5"} = "HP3.32";
		$genes1{"Roslin2.38"} = "HP3.32";
		$genes1{"Roslin2.7"} = "HP3.32";
		print OUT_GENES "Roslin2.5:".$genes{"Roslin2.5"}.",Roslin2.38:".$genes{"Roslin2.38"}.",Roslin2.7:".$genes{"Roslin2.7"}."\t";
	}
}

sub ambiguous {
 	$gene = $_[0];
	if ($gene eq "3*00102"){
		$flag = 1;
		if (!exists $genes1{"3*00103"}){
			$flag1 = 1;
			$name = "3*00102|3*00103";
			$genes1{"3*00102"} = 1;
			$genes1{"3*00103"} = 1;
	}	}
	if ($gene eq "3*00103"){
		$flag = 1;
		if (!exists $genes1{"3*00102"}){
			$flag1 = 1;
			$name = "3*00102|3*00103";
			$genes1{"3*00102"} = 1;
			$genes1{"3*00103"} = 1;
	}	}
	if ($gene eq "2*00601"){
		$flag = 1;
		if (!exists $genes1{"2*00602"}){
			$flag1 = 1;
			$name = "2*00601|2*00602";
			$genes1{"2*00601"} = 1;
			$genes1{"2*00602"} = 1;
	}	}
	if ($gene eq "2*00602"){
		$flag = 1;
		if (!exists $genes1{"2*00601"}){
			$flag1 = 1;
			$name = "2*00601|2*00602";
			$genes1{"2*00601"} = 1;
			$genes1{"2*00602"} = 1;
	}	}
	if ($gene eq "2*02601"){
		$flag = 1;
		if (!exists $genes1{"2*02603"}){
			$flag1 = 1;
			$name = "2*02601|2*02603";
			$genes1{"2*02601"} = 1;
			$genes1{"2*02603"} = 1;
	}	}
	if ($gene eq "2*02603"){
		$flag = 1;
		if (!exists $genes1{"2*02601"}){
			$flag1 = 1;
			$name = "2*02601|2*02603";
			$genes1{"2*02601"} = 1;
			$genes1{"2*02603"} = 1;
	}	}
	if ($primer eq "For1Rev2"){
		if ($gene eq "3*00402"){
			$flag = 1;
			if (!exists $genes1{"3*05301"}){
				$flag1 = 1;
				$name = "3*00402|3*05301";
				$genes1{"3*00402"} = 1;
				$genes1{"3*05301"} = 1;
		}	}
		if ($gene eq "3*05301"){
			$flag = 1;
			if (!exists $genes1{"3*00402"}){
				$flag1 = 1;
				$name = "3*00402|3*05301";
				$genes1{"3*00402"} = 1;
				$genes1{"3*05301"} = 1;
		}	}
		if ($gene eq "Roslin3.46"){
			$flag = 1;
			if (!exists $genes1{"2*01801"}){
				$flag1 = 1;
				$name = "2*01801|2*01802|Roslin3.46";
				$genes1{"2*01801"} = 1;
				$genes1{"2*01802"} = 1;
				$genes1{"Roslin3.46"} = 1;
				return;
		}	}
		if ($gene eq "2*01801"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.46"} and !exists $genes1{"2*01802"}){
				$flag1 = 1;
				$name = "2*01801|2*01802|Roslin3.46";
				$genes1{"2*01801"} = 1;
				$genes1{"2*01802"} = 1;
				$genes1{"Roslin3.46"} = 1;
				return;
		}	}
		if ($gene eq "2*01802"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.46"} and !exists $genes1{"2*01801"}){
				$flag1 = 1;
				$name = "2*01801|2*01802|Roslin3.46";
				$genes1{"2*01801"} = 1;
				$genes1{"2*01802"} = 1;
				$genes1{"Roslin3.46"} = 1;
				return;
		}	}
		if ($gene eq "2*00801"){
			$flag = 1;
			if (!exists $genes1{"2*00802"}){
				$flag1 = 1;
				$name = "2*00801|2*00802";
				$genes1{"2*00801"} = 1;
				$genes1{"2*00802"} = 1;
		}	}
		if ($gene eq "2*00802"){
			$flag = 1;
			if (!exists $genes1{"2*00801"}){
				$flag1 = 1;
				$name = "2*00801|2*00802";
				$genes1{"2*00801"} = 1;
				$genes1{"2*00802"} = 1;
		}	}
		if ($gene eq "2*01601"){
			$flag = 1;
			if (!exists $genes1{"2*01602"}){
				$flag1 = 1;
				$name = "2*01601|2*01602";
				$genes1{"2*01601"} = 1;
				$genes1{"2*01602"} = 1;
		}	}
		if ($gene eq "2*01602"){
			$flag = 1;
			if (!exists $genes1{"2*01601"}){
				$flag1 = 1;
				$name = "2*01601|2*01602";
				$genes1{"2*01601"} = 1;
				$genes1{"2*01602"} = 1;
		}	}
		if ($gene eq "6*04101"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.30"}){
				$flag1 = 1;
				$name = "6*04101|Roslin2.30";
				$genes1{"6*04101"} = 1;
				$genes1{"Roslin2.30"} = 1;
		}	}
		if ($gene eq "Roslin2.30"){
			$flag = 1;
			if (!exists $genes1{"6*04101"}){
				$flag1 = 1;
				$name = "6*04101|Roslin2.30";
				$genes1{"6*04101"} = 1;
				$genes1{"Roslin2.30"} = 1;
		}	}
		if ($gene eq "1*06101"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.18"}){
				$flag1 = 1;
				$name = "1*06101|Roslin2.18";
				$genes1{"1*06101"} = 1;
				$genes1{"Roslin2.18"} = 1;
		}	}
		if ($gene eq "Roslin2.18"){
			$flag = 1;
			if (!exists $genes1{"1*06101"}){
				$flag1 = 1;
				$name = "1*06101|Roslin2.18";
				$genes1{"1*06101"} = 1;
				$genes1{"Roslin2.18"} = 1;
		}	}
		if ($gene eq "2*05401"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.8"}){
				$flag1 = 1;
				$name = "2*05401|Roslin3.8";
				$genes1{"2*05401"} = 1;
				$genes1{"Roslin3.8"} = 1;
		}	}
		if ($gene eq "Roslin3.8"){
			$flag = 1;
			if (!exists $genes1{"2*05401"}){
				$flag1 = 1;
				$name = "2*05401|Roslin3.8";
				$genes1{"2*05401"} = 1;
				$genes1{"Roslin3.8"} = 1;
		}	}
		if ($gene eq "Roslin2.20"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.18"}){
				$flag1 = 1;
				$name = "Roslin2.20|Roslin3.18";
				$genes1{"Roslin2.20"} = 1;
				$genes1{"Roslin3.18"} = 1;
		}	}
		if ($gene eq "Roslin3.18"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.20"}){
				$flag1 = 1;
				$name = "Roslin2.20|Roslin3.18";
				$genes1{"Roslin2.20"} = 1;
				$genes1{"Roslin3.18"} = 1;
		}	}
		if ($gene eq "Roslin2.40"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.1"}){
				$flag1 = 1;
				$name = "Roslin2.40|Roslin3.1";
				$genes1{"Roslin2.40"} = 1;
				$genes1{"Roslin3.1"} = 1;
		}	}
		if ($gene eq "Roslin3.1"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.40"}){
				$flag1 = 1;
				$name = "Roslin2.40|Roslin3.1";
				$genes1{"Roslin2.40"} = 1;
				$genes1{"Roslin3.1"} = 1;
		}	}
		if ($gene eq "Roslin2.5"){
			$flag = 1;
			if (!exists $genes1{"Unassigned2.6"}){
				$flag1 = 1;
				$name = "Roslin2.5|Unassigned2.6";
				$genes1{"Roslin2.5"} = 1;
				$genes1{"Unassigned2.6"} = 1;
		}	}
		if ($gene eq "Unassigned2.6"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.5"}){
				$flag1 = 1;
				$name = "Roslin2.5|Unassigned2.6";
				$genes1{"Roslin2.5"} = 1;
				$genes1{"Unassigned2.6"} = 1;
		}	}
		if ($gene eq "6*01402"){
			$flag = 1;
			if (!exists $genes1{"Roslin3.37"}){
				$flag1 = 1;
				$name = "6*01402|Roslin3.37";
				$genes1{"6*01402"} = 1;
				$genes1{"Roslin3.37"} = 1;
		}	}
		if ($gene eq "Roslin3.37"){
			$flag = 1;
			if (!exists $genes1{"6*01402"}){
				$flag1 = 1;
				$name = "6*01402|Roslin3.37";
				$genes1{"6*01402"} = 1;
				$genes1{"Roslin3.37"} = 1;
		}	}
	}
	if($primer eq "For3Rev1"){
		if ($gene eq "2*01801"){
			$flag = 1;
			if (!exists $genes1{"2*01802"}){
				$flag1 = 1;
				$name = "2*01801|2*01802";
				$genes1{"2*01801"} = 1;
				$genes1{"2*01802"} = 1;
		}	}
		if ($gene eq "2*01802"){
			$flag = 1;
			if (!exists $genes1{"2*01801"}){
				$flag1 = 1;
				$name = "2*01801|2*01802";
				$genes1{"2*01801"} = 1;
				$genes1{"2*01802"} = 1;
		}	}
		if ($gene eq "Roslin2.16"){
			$flag = 1;
			if (!exists $genes1{"3*00402"}){
				$flag1 = 1;
				$name = "3*00402|3*05301|Roslin2.16";
				$genes1{"3*00402"} = 1;
				$genes1{"3*05301"} = 1;
				$genes1{"Roslin2.16"} = 1;
				return;
		}	}
		if ($gene eq "3*00402"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.16"} and !exists $genes1{"3*05301"}){
				$flag1 = 1;
				$name = "3*00402|3*05301|Roslin2.16";
				$genes1{"3*00402"} = 1;
				$genes1{"3*05301"} = 1;
				$genes1{"Roslin2.16"} = 1;
				return;
		}	}
		if ($gene eq "3*05301"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.16"} and !exists $genes1{"3*00402"}){
				$flag1 = 1;
				$name = "3*00402|3*05301|Roslin2.16";
				$genes1{"3*00402"} = 1;
				$genes1{"3*05301"} = 1;
				$genes1{"Roslin2.16"} = 1;
				return;
		}	}
		if ($gene eq "2*03201N"){
			$flag = 1;
			if (!exists $genes1{"2*03202"}){
				$flag1 = 1;
				$name = "2*03201N|2*03202";
				$genes1{"2*03201N"} = 1;
				$genes1{"2*03202"} = 1;
		}	}
		if ($gene eq "2*03202"){
			$flag = 1;
			if (!exists $genes1{"2*03201N"}){
				$flag1 = 1;
				$name = "2*03201N|2*03202";
				$genes1{"2*03201N"} = 1;
				$genes1{"2*03202"} = 1;
		}	}
		if ($gene eq "1*02901"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.43"}){
				$flag1 = 1;
				$name = "1*02901|Roslin2.43";
				$genes1{"1*02901"} = 1;
				$genes1{"Roslin2.43"} = 1;
		}	}
		if ($gene eq "Roslin2.43"){
			$flag = 1;
			if (!exists $genes1{"1*02901"}){
				$flag1 = 1;
				$name = "1*02901|Roslin2.43";
				$genes1{"1*02901"} = 1;
				$genes1{"Roslin2.43"} = 1;
		}	}
		if ($gene eq "Roslin2.19"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.22"}){
				$flag1 = 1;
				$name = "Roslin2.19|Roslin2.22";
				$genes1{"Roslin2.19"} = 1;
				$genes1{"Roslin2.22"} = 1;
		}	}
		if ($gene eq "Roslin2.22"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.19"}){
				$flag1 = 1;
				$name = "Roslin2.19|Roslin2.22";
				$genes1{"Roslin2.19"} = 1;
				$genes1{"Roslin2.22"} = 1;
		}	}
		if ($gene eq "Roslin2.38"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.25"}){
				$flag1 = 1;
				$name = "Roslin2.38|Roslin2.25";
				$genes1{"Roslin2.38"} = 1;
				$genes1{"Roslin2.25"} = 1;
		}	}
		if ($gene eq "Roslin2.25"){
			$flag = 1;
			if (!exists $genes1{"Roslin2.38"}){
				$flag1 = 1;
				$name = "Roslin2.38|Roslin2.25";
				$genes1{"Roslin2.38"} = 1;
				$genes1{"Roslin2.25"} = 1;
		}	}
		if ($gene eq "Roslin1.10"){
			$flag = 1;
			if (!exists $genes1{"Roslin1.2"}){
				$flag1 = 1;
				$name = "Roslin1.10|Roslin1.2";
				$genes1{"Roslin1.10"} = 1;
				$genes1{"Roslin1.2"} = 1;
		}	}
		if ($gene eq "Roslin1.2"){
			$flag = 1;
			if (!exists $genes1{"Roslin1.10"}){
				$flag1 = 1;
				$name = "Roslin1.10|Roslin1.2";
				$genes1{"Roslin1.10"} = 1;
				$genes1{"Roslin1.2"} = 1;
		}	}
		if ($gene eq "RoslinSZ.1.1"){
			$flag = 1;
			if (!exists $genes1{"RoslinSZ.3.1"}){
				$flag1 = 1;
				$name = "RoslinSZ.1.1|RoslinSZ.3.1";
				$genes1{"RoslinSZ.1.1"} = 1;
				$genes1{"RoslinSZ.3.1"} = 1;
		}	}
		if ($gene eq "RoslinSZ.3.1"){
			$flag = 1;
			if (!exists $genes1{"RoslinSZ.1.1"}){
				$flag1 = 1;
				$name = "RoslinSZ.1.1|RoslinSZ.3.1";
				$genes1{"RoslinSZ.1.1"} = 1;
				$genes1{"RoslinSZ.3.1"} = 1;
		}	}
	}
}    
print OUT_GENES "$prefix\t";
print OUT_HAPLO "$prefix";
print OUT_NONCLASSIC "$prefix\t";
print OUT_NONHAPLO "$prefix\t";

my %nonClassic;

while(<IN>){
	chomp $_;
	@words = split("\t",$_);
	@header = split("-",$words[0]);
	if ($words[1] =~ /^NC/){
		$count_nonclassic++;
		$nonClassic{$words[1]} = $header[3];
	}
	else{
        $count_classic++;
		$genes{$words[1]} = $header[3];
	} 	
}

Haplotypes; #Call function to assign haplotypes based on genes found

foreach $gene (keys %genes){
 	if (! exists $genes1{$gene}){
 		$flag = 0;	
 		$flag1 = 0;
 		ambiguous($gene);
 		if ($flag == 0){
 			$nonHaplo{$gene} = $genes{$gene};
	 	}
	 	else{
	 		$count_classic--;
	 		if ($flag1 == 1){
	 			$nonHaplo{$name} = $genes{$gene};
	 		}
	 	}
}	}		
$count_haplo = scalar @haplo;
foreach $gene (keys %nonClassic){
	if ($gene eq "NC1-00201" and $nonClassic{$gene} eq $nonClassic{"NC1-00301"}){
		delete $nonClassic{"NC1-00301"};
		$count_nonclassic--;
	}
	elsif ($gene eq "NC1-00101" and $nonClassic{$gene} eq $nonClassic{"NC1-00101SV"}){
		delete $nonClassic{"NC1-00101SV"};
		$count_nonclassic--;
	}
	elsif ($gene eq "NC2-00101" and $nonClassic{$gene} eq $nonClassic{"NC2-00102"}){
		delete $nonClassic{"NC2-00102"};
		$count_nonclassic--;
	}
}
foreach $gene (sort {$nonClassic{$b} <=> $nonClassic{$a}} keys %nonClassic){
		print OUT_NONCLASSIC "$gene:$nonClassic{$gene}\t";
}
foreach $gene (sort {$nonHaplo{$b} <=> $nonHaplo{$a}} keys %nonHaplo){
	print OUT_NONHAPLO "$gene:$nonHaplo{$gene}\t";
}

$total_genes = $count_classic + $count_nonclassic;
print LOG "$prefix\t$total_genes\t$count_classic\t$count_nonclassic\n";
print OUT_GENES "\n";
print OUT_HAPLO "\t$count_haplo\t@haplo\n";
print OUT_NONCLASSIC "\n";
print OUT_NONHAPLO"\n";	




