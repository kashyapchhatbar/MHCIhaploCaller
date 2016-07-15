# MHCIhaploCaller

# Overview

MHCI haplotype caller is written to analyse high-throughput data of MHCI alleles in cattle. 
The reference paper is: Deepali Vasoya, Andy Law, Paolo Motta, Mingyan Yu, Adrian Muwonge, Elizabeth Cook, Xiaoying Li, Karen Bryson, Amanda MacCallam, Tatjana Sitt, PhilipToye, Barend Bronsvoort, Mick Watson, W. Ivan Morrison and Timothy Connelley. **_"Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations Using Next-Generation Sequencing."_**

MHCI haplotype caller is the pipeline of Perl scripts to analyse Illumina MiSeq data of MHCI alleles in cattle. This software uses only pre-defined haplotype information based on MHCI genes deposited in IPD database and novel haplotypes and MHCI allele found in our study. We continue working on it to make this tool more generic in terms of species, database and haplotypes. 

# Specificity

This pipeline is written to analyse 288 cattle - cohort of Holstein-Friesian, Boran and Cameroonian. Bovine MHCI primer sets amplifying the highly polymorphic exon 2 and exon 3 were designed which allow amplification and nearly unambiguous identification of all known MHCI alleles. There primers were used to generate amplicons (no greater than 500bp) from cohorts of 96 animals each time that were multiplexed and overlapping paired end sequencing data generated using Illumina MiSeq. There are 2 sets of primers designed: For1Rev2 and For3Rev1, spanning most of the hyper variable exons 2 and 3 that encode the alpha1 and alpha2 domains. The length of the amplicons generated by these primer sets were anticipated to be 378 and 318 by respectively and together cover 410bp. The database has the MHCI alleles sequences of the same size as the amplicons. The IDP database sequences are trimmed using the primer sequences. The newly found MHCI alleles are named with the prefix of Roslin, NCRoslin and Unassigned. More details can be found in paper.

# Installation/Download

* Please download latest version of MHCIhaploCaller.
* Unzip the folder.
* Enter into the folder.
* Run the mhcI_haplotype.pl script using the usage manual below. 

# Usage

mhcI_haplotype.pl is the pipeline with 4 different modules. They can only be run in fixed order as later steps are dependent on output delivered by previous step. 
The MiSeq data for which this software is written is available to download from here: PRJEB14552

```
$ perl mhcI_haplotype.pl
Please enter command:
preAlign - to analyse the paired-end raw/quality trimmed fastq files.
Blast - to align processed unique aligned reads on database
postAlign - to analyse the blast results
analyseUnknowns - to find out which sample has unknown sequences
```

NOTE: Make one project directory and use it (use the full path to directory) in all the steps to save the outputs. This directory will have sub-directories called Samples and Reports. Samples will have detailed analysis report for individual samples. Reports will have overall merged reports of all samples together.  

1. preAlign:
This is the first step of the pipeline. It takes the mandatory options shown below. Please provide the path to [FLASH](https://ccb.jhu.edu/software/FLASH/) and [fastx-toolkit]( http://hannonlab.cshl.edu/fastx_toolkit/). It won’t run without it. 

```
$ perl mhcI_haplotype.pl preAlign
Usage:
--fastqDir	/Directory/of/Fastq/Files
--total_samples	Total number of samples
--prefix	prefix of the fastq files, Boran, Cameroonian or HolsteinFriesian
--dir	/Full/path/to/project/directory
--FlashPath	/Full/Path/of/FLASH/software
--fastxToolkitPath	/Full/Path/of/fastx-toolkit/software
```

It takes the raw fastq files. The fastq file names should be as per the ENA data (PRJEB14552). The overlapping sequencing reads are converted into extended fragments using FLASH. Fastx-toolkit removes bad quality reads and converts into fasta format. The fasta files are then demultiplexed into two primer sets. Now the further steps process these two primer sets separately. The multiple copies of unique sequences are merged and ordered from high to low abundance. The header of these fasta files has the information on each sequence – ‘>index-length-frequency-%frequency’. It also calculates the threshold frequency of sequence at 0.2%. This unique sequences will go further for alignment on database. 

2. Blast
It uses the unique sequences received from preAlign and aligns the sequences on database. The database is in built. The full path to installed BLAST must be provided. It takes the following options:

```
$ perl mhcI_haplotype.pl Blast
Usage:
--blastPath	/Full/Path/of/Blast/software
--dir	/Full/path/to/project/directory
--total_samples	Total number of samples
--onlySelected	align only selected sequences which are above threshold, TRUE or FALSE
```

--onlySelected option is important to save the time. TRUE will only aligns the selected sequences which are above the frequency cutoff while FALSE will aligns all unique sequences. FALSE alignment will take extreme long time to run BLAST, so TRUE option is highly recommended.  

3. postAlign:
It uses the blast output. Based on alignment information on in-built database, it identifies the known alleles. It also removes the unwanted sequences such as chimaeras, 1 or 2 base variants of known alleles, +/- 9bp different from expteced length. The filtered sequences are then assigned as unknowns.

```
$ perl mhcI_haplotype.pl postAlign
Usage:
--dir	/Full/path/to/project/directory
--total_samples	Total number of samples
```

4. analyseUnknowns
It uses the unknown sequences and gives unique identifiers. It analyses the occurrences of each unknown in all samples and makes an informative list. The followings are the required options:

```
$ perl mhcI_haplotype.pl analyseUnknowns
Usage:
--dir	/Full/path/to/main/directory
--total_samples	Total number of samples
```

# License
Copyright (c) 2015 Deepali Vasoya

Permission is hereby granted, free of charge, to any person obtaining a copyof this software and associated documentation files (the "Software"), to dealin the Software without restriction, including without limitation the rightsto use, copy, modify, merge, publish, distribute, sublicense, and/or sellcopies of the Software, and to permit persons to whom the Software isfurnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in allcopies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS ORIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THEAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHERLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Cite
Deepali Vasoya, Andy Law, Paolo Motta, Mingyan Yu, Adrian Muwonge, Elizabeth Cook, Xiaoying Li, Karen Bryson, Amanda MacCallam, Tatjana Sitt, PhilipToye, Barend Bronsvoort, Mick Watson, W. Ivan Morrison and Timothy Connelley. **_"Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations Using Next-Generation Sequencing."_**
