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

use Getopt::Long;
use 5.006;
use Error;
use warnings;

my $dir;
my $prefix;
my $primer;
GetOptions(
	'dir=s' => \$dir, 
	'prefix=s' => \$prefix,
	'primer=s' => \$primer,
) or die "Usage:\n--dir\tdirectory of blast output\n--prefix\tprefix for output files\n--primer\tprimer label\n\n";

if(! defined $dir or ! defined $prefix or ! defined $primer){
	die "Usage:\n--dir\tdirectory of blast output\n--prefix\tprefix for output files\n--primer\tprimer label\n\n";
}
my $len;
my $len1;
my $len2;

if ($primer eq "For1Rev2"){
	$len = 378;
	$len1 = 376;
	$len2 = 375;
}
elsif ($primer eq "For3Rev1"){
	$len = 318;
	$len1 = 316;
	$len2 = 315;
}
else{
	die "Please provide proper primer name for analyseBlast.pl\n\n";
}

open(IN, "$dir/$prefix.$primer.blast") or die "Cannot open $dir/$prefix.$primer.blast in analyseBlast.pl\n";
open(FASTA, "$dir/$prefix.$primer.selected.fasta") or die "Cannot open $dir/$prefix.$primer.selected.fasta in analyseBlast.pl\n";
open(KNOWN, ">$dir/$prefix.$primer.known.genes.txt") or die "Cannot open $dir/$prefix.$primer.known.genes.txt in analyseBlast.pl\n";
open(UNKNOWN, ">$dir/$prefix.$primer.unknown.fasta") or die "Cannot open $dir/$prefix.$primer.unknown.fasta in analyseBlast.pl\n";
open(LEN9, ">$dir/$prefix.$primer.Len-9nc.diff.fasta") or die "Cannot open $dir/$prefix.$primer.Len-9nc.diff.fasta in analyseBlast.pl\n";
open(DETAIL, ">$dir/$prefix.$primer.selected.details.txt") or die "Cannot open $dir/$prefix.$primer.selected.details.txt in analyseBlast.pl\n";
open(LOG, ">$dir/$prefix.$primer.analyseBlast.log") or die "Cannot open $dir/$prefix.$primer.analyseBlast.log in analyseBlast.pl\n";

my $countSingleton = 0;
my $countKnown = 0;
my $countUnknown = 0;
my $countChimaera = 0;
my $countLen9 = 0;
my $countLen = 0;
my $countVar = 0;
my %genes = ();
my $id;
my %kn;
my @words;
my @info;
my %copies = ();
my %chimaera = ();
my $diff;

while (<FASTA>){
	chomp $_;
	if (/>(.+)/){
		$id = $1;
	}
	else{
		@words = split ("-", $id);
		$copies{$id} = $words[2];
		$genes{$id} = $_;	
		$lengths{$id} = $words[1];
}	}

while(<IN>){
	chomp $_;
	@words = split("\t", $_);
	if (exists $genes{$words[0]}){
		if ($words[0] ne $id){
			$id = $words[0];
			@info = split("-", $words[0]);
			if ($info[1] eq $words[3] and $words[2] eq "100.00" and $words[3] eq $len and $words[6] eq 1 and $words[7] eq $len and $words[8] eq 1 and $words[9] eq $len){
				$kn{$words[0]} = $words[1];	
			}
			elsif ($info[1] eq $words[3] and $words[1] eq "3*03301N" and $words[2] eq "100.00" and $words[3] eq $len1 and $words[6] eq 1 and $words[7] eq $len1 and $words[8] eq 1 and $words[9] eq $len1){
				$kn{$words[0]} = $words[1];
			}
			elsif ($words[1] eq "2*05701" or $words[1] eq "NC2-00101" or $words[1] eq "NC2-00102" or $words[1] eq "NC2-00103"){
				if ($info[1] eq $words[3] and $words[2] eq 100.00 and $words[3] eq $len2 and $words[6] eq 1 and $words[7] eq $len2 and $words[8] eq 1 and $words[9] eq $len2){
					$kn{$words[0]} = $words[1];
			}	}
			elsif ($info[1] eq $words[3] and $words[1] =~ /RoslinSZ/ and $words[2] eq "100.00" and $words[3] eq $info[1] and $words[6] eq 1 and $words[7] eq $info[1] and $words[8] eq 1 and $words[9] eq $info[1]){
				$kn{$words[0]} = $words[1];
			}	
			else{
				next;
		}	}
		else{
			if ($words[2] eq "100.00" and $words[3] eq $len and $words[6] eq 1 and $words[7] eq $len and $words[8] eq 1 and $words[9] eq $len){
				$kn{$words[0]} = $kn{$words[0]}."\t".$words[1];
			}
			elsif ($words[1] eq "3*03301N" and $words[2] eq "100.00" and $words[3] eq $len1 and $words[6] eq 1 and $words[7] eq $len1 and $words[8] eq 1 and $words[9] eq $len1){
				$kn{$words[0]} = $kn{$words[0]}."\t".$words[1];
			}
			elsif ($words[1] eq "2*05701" or $words[1] eq "NC2-00101" or $words[1] eq "NC2-00102" or $words[1] eq "NC2-00103"){
				if ($words[2] eq "100.00" and $words[3] eq $len2 and $words[6] eq 1 and $words[7] eq $len2 and $words[8] eq 1 and $words[9] eq $len2){
					$kn{$words[0]} = $kn{$words[0]}."\t".$words[1];
			}	}
			elsif ($info[1] eq $words[3] and $words[1] =~ /RoslinSZ/ and $words[2] eq "100.00" and $words[3] eq $info[1] and $words[6] eq 1 and $words[7] eq $info[1] and $words[8] eq 1 and $words[9] eq $info[1]){
				$kn{$words[0]} = $kn{$words[0]}."\t".$words[1];
			}	
	}	}
	else{
		last;
}	}

sub new {
	# my $invocant = shift;
	# my $class = ref($invocant) || $invocant;
	my $self = {@_};
	bless( $self );
	foreach my $index ( 1, 2 ) {
		my $haplo = "haplotype${index}";
		throw Error::Simple( "No haplotype${haplo}" )
		  unless defined( $self->{$haplo} );
	}
	my $len = length($self->{'haplotype1'});
	throw Error::Simple( "Haplotypes differ in length" )
		unless ( $len == length($self->{'haplotype2'}));

	throw Error::Simple( "Haplotypes are empty" )
		unless ( $len > 0);
		
	$self->{'haplotype1'} = uc($self->{'haplotype1'});
	$self->{'haplotype2'} = uc($self->{'haplotype2'});

	my $matcher = {};
	for (my $i = 0; $i < $len; $i++) {
		my $base = {};
		$base->{substr($self->{'haplotype1'}, $i, 1 )} += 1;
		$base->{substr($self->{'haplotype2'}, $i, 1 )} -= 1;
		
		$matcher->{$i} = $base;
	}
	$self->{'matcher'} = $matcher;
	$self->{'length'} = $len;
	
	return $self;
}
sub possible_chimaera {
	my $self = shift;
	my $candidate = shift;
	
	$candidate = uc($candidate);
	
	my %to_string = ( 0 => "", -1 => "A", 1 => "B");
	
	if ( length($candidate) != $self->{'length'} ) {
		return 0;
	}
	
	my $coded_string = "";
	for (my $i = 0; $i < $self->{'length'}; $i++) {
		my $base = substr($candidate, $i, 1);
		if (defined($self->{'matcher'}{$i}{$base})) {
			$coded_string .= $to_string{$self->{'matcher'}{$i}{$base}};
		}
		else {
			return 0;
		}
	}
	# print $coded_string, " - CODED\n";
	if ($coded_string =~ m/^(A+B+|B+A+)$/) {
		return 1;
	}
	return 0;	
}

my $match = 0;
my $sub1 = undef;
my $sub2 = undef;
foreach my $gene (sort{$copies{$a}<=>$copies{$b}} keys %copies){
	my $flag = 0;
	if ($copies{$gene} eq 1){
		print DETAIL "$gene\tSINGLETON\n";
		$countSingleton++;
		$flag = 1;
	}
	elsif (exists $kn{$gene}){
		@words = split("\t", $kn{$gene});
		foreach $known (@words){
			print KNOWN "$gene\t$known\n";
			print DETAIL "$gene\tKNOWN\t$known\n";
		}
		$countKnown++;
		$flag = 1;
	}
	elsif ($lengths{$gene} != $len){
		$diff = $len - $lengths{$gene};
		if ($lengths{$gene} <= ($len + 9) and $lengths{$gene} >= ($len - 9)){
			print DETAIL "$gene\tLENDIFF-9NC\t$diff\n";
			print LEN9 ">$gene\n$genes{$gene}\n";
			$countLen9++;
			$flag = 1;
		}
		else{
			print DETAIL "$gene\tLENDIFF\t$diff\n";
			$countLen++;
			$flag = 1;
		}
	}
	elsif ($lengths{$gene} eq $len){
		LOOP1: foreach my $gene1 (sort {$copies{$b} <=> $copies{$a}} keys %copies){
			if ($gene ne $gene1 and $lengths{$gene1} eq $len){
				foreach my $gene2 (sort {$copies{$b} <=> $copies{$a}} keys %copies){
					if ($gene ne $gene2 and $gene1 ne $gene2 and $lengths{$gene2} eq $len){
						# print "gene:$gene\tgene1:$gene1\tgene2:$gene2\n";
						my $check = new('haplotype1' => $genes{$gene1}, 'haplotype2' => $genes{$gene2});
						if ($check->possible_chimaera($genes{$gene})){
							if ($copies{$gene} < $copies{$gene1} and $copies{$gene} < $copies{$gene2}){
								$chimaera{$gene} = "$gene1\t$gene2";
								$flag = 1;
								$countChimaera++;
								print DETAIL "$gene\tCHIMAERA\t$chimaera{$gene}\n";
								# print "$gene $chimaera{$gene}\n";
								last LOOP1;
							}
						}
					}
				}
			}
		}
		if (! exists $chimaera{$gene} and $flag == 0){
			LOOP4: foreach $gene1 (sort{$copies{$b}<=>$copies{$a}} keys %copies){
				if (($copies{$gene1} / $copies{$gene}) >= 100){
					my $i = 0;
					my $match = 0;
					my $var_base = undef;
					my $variant_gene = undef;
					while ($i < $len){
						my $sub1 = substr ($genes{$gene}, $i, 1);
						my $sub2 = substr ($genes{$gene1}, $i, 1);
						if ($sub1 eq $sub2){
							$match++;
						}
						$i++;
					}
					if ((($copies{$gene1} / $copies{$gene}) >= 30 and $match eq $len-1) or (($copies{$gene1} / $copies{$gene}) >= 100 and $match eq $len-2)){
						$countVar++;
						$flag = 1;	
						print DETAIL "$gene\tVARIANT\t$gene1\t$kn{$gene1}\n";
						last LOOP4;
					}
				}
			}
			if ($flag == 0){
				print UNKNOWN ">$prefix\t$gene\n$genes{$gene}\n";
				print DETAIL "$gene\tUNKNOWN\n";
				$countUnknown++;
			}
		}
	}
}

$total = $countSingleton + $countKnown + $countUnknown + $countChimaera + $countLen9 + $countLen + $countVar;
print LOG "$prefix\t$total\t$countKnown\t$countUnknown\t$countChimaera\t$countLen9\t$countLen\t$countVar\t$countSingleton\n";





	
