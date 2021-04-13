------------------------------------------------------
--- drFAST:di-base read Fast Alignment Search Tool ---
---               VERSION 1.0.0.0                  ---
------------------------------------------------------

Table of Contents
=================

- Introduction
- Installation
- Usage
- Returned Files
- Parameters
- Quick Start
- Upcoming features
- Changes to previous version
- Additional Information


Introduction
============

drFAST is designed to map short reads generated with the AB Solid platform to 
reference genome assemblies; in a fast and memory-efficient manner.

Installation
============

On Unix systems, we recommend using GNU gcc 4.1.2 as your compiler and type 'make' to build.

Example:
        linux> make

if the compile is successful, output on the screen is as below:

gcc4 -c -o2 baseFAST.c -o baseFAST.o
gcc4 -c -o2 CommandLineParser.c -o CommandLineParser.o
gcc4 -c -o2 Common.c -o Common.o
gcc4 -c -o2 HashTable.c -o HashTable.o
gcc4 -c -o2 DrFAST.c -o DrFAST.o
gcc4 -c -o2 Output.c -o Output.o
gcc4 -c -o2 Reads.c -o Reads.o
gcc4 -c -o2 RefGenome.c -o RefGenome.o
gcc4 baseFAST.o CommandLineParser.o Common.o HashTable.o DrFAST.o Output.o Reads.o RefGenome.o  -o DrFAST -static -lz -lm
rm -rf *.o

Usage
=====


I.Indexing
	To index a reference genome like "refgen.fasta" run the following command:

		$>./drFAST --index refgen.fasta

	Upon the completion of the indexing phase, you can find "refgen.fasta.index" in the same directory as "refgen.fasta". 
	drFAST uses a window size of 12 (default) to make the index of the genome, this windows size can be modified with "--ws". 
	There is a restriction on the maximum of the window size as the	window size directly affects the memory usage.

		$>./drFAST --index refgen.fasta --ws 13	

NOTE: You can index more than one chromosom(segment, )

II.Mapping
	A. Single-end Reads

	To map single reads to a reference genome, run the following command. Use "--seq" to specify the input file. 
	refgen.fa and refgen.fa.index should be in the same folder.

		$>./drFAST --search refgen.fa --seq reads.fastq

	The reported locations will be saved into "output" by default. If you want to save it somewhere else, use "-o" 
	to specify another file. mrFAST can report the unmapped reads in fasta/fastq format.

		$>./drFAST --search refgen.fasta --seq reads.fastq -o myoutput 

	B. Paired-end Reads

	To map paired-end reads, use "--pe" option. The distance allowed between the paired-end reads should be specified with "--min" and "--max". 
	"--min" and "--max" specify the minmum and maximum of the inferred size (the distance between outer edges of the mapping mates).

		$>./drFAST --search refgen.fasta --pe --seq reads.fastq --min 150 --max 250 

	In order to get all the discordant mapping.
		$>./drFAST --search refgen.fasta --pe --seq reads.fastq --min 150 --max 250 --discordant-vh


Returned Files 
===============

A. Single-end Reads

		$>./drFAST --search refgen.fasta --seq reads.fastq -o myoutput.sam
	It will generate file called myoutput.sam in SAM format which contains all the read mapping and their locations.

B. Paired-end Reads

		$>./drFAST --search refgen.fasta --mp --seq reads.fastq --min 150 --max 250 --discordant-vh -o myoutput.sam --discordant-vh

	It will generate 5 files

		myoutput.sam: Contains the best concordant(distance between two read is between min and max value specified in command) 
			  and best discordant(distance between two read is not between min and max specified in command) pair end mappings (SAM format).
			
		myoutput__BEST.CONCORDANT: Contains the best concordant, concordant reads which has the minimum edit distance to 
					   reference and the distance is closer to mean of min and max value specified in the command (SAM format).

		myoutput__BEST.DISCONCORDANT: Contains the best discordant, discordant reads which has the minimum edit distance to 
					      reference and the distance is closer to mean of min and max value specified in the command (SAM format).

		myoutput__OEA1:Contains the reads in fastq format which specifies which reads their first pair (/1) didnt map but their 
				mate pair (/2) mapped to	the reference.

		myoutput__OEA2:Contains the reads in fastq format which specifies which reads their first pair (/1) did map but their 
				mate pair (/2) didnt mapped to the reference.

		myoutput__DIVET.vh: Contains all the discordant mapping reads which can be used later for Varation Hunter (other Structural Varation
		softwares), if is in .vh format.

Parameters
===========
General Options:
 -v|--version		Current Version.
 -h			Shows the help file.


Indexing Options:
 --index [file]		Generate an index from the specified fasta file. 
 -b			Indicates the indexing will be done in batch mode.
			The file specified in --index should contain the 
			list of fasta files.
 -ws [int]		Set window size for indexing (default:12-max:14).


Searching Options:
 --search [file]	Search the specified genome. Index file should be 
			in same directory as the fasta file.
 -b			Indicates the mapping will be done in batch mode. 
			The file specified in --search should contain the 
			list of fasta files.
 --mp 			Search will be done in Pairedend mode.
 --seq [file]		Input sequences in fasta/fastq format [file]. If 
			pairend reads are interleaved, use this option.
 --seq1 [file]		Input sequences in fasta/fastq format [file] (First 
			file). Use this option to indicate the first file of 
			pair-end reads. 
 --seq2 [file]		Input sequences in fasta/fastq format [file] (Second 
			file). Use this option to indicate the second file of 
			pair-end reads.  
 -o [file]		Output of the mapped sequences. The default is output.
 --seqcomp 		Indicates that the input sequences are compressed(gz).
 --outcomp 		Indicates that output file should be compressed(gz).
 -e [int]		edit distance (default 2).
 --min [int]		Min inferred distance allowed between two pairend sequences.
 --max [int]		Max inferred distance allowed between two pairend sequences.
 --discordant-vh	To generate all the discordant mapping for Variation Hunter 
			program
 --best			Returns the best mapping only
		
Quick Start
===========


I.Indexing
	
	$>./drFAST --index chr1_random 

	Output on the screen:
		Generating Index from chr1_random
		 - chr1_random 
		DONE in 1.55s!

	$>./drFAST --index hg18
		
		It will generate the index for hg18 (contains all the chromosomes in one file).	

II.Mapping
	
	A. Single-end Reads

	$>./drFAST --search example/chr1_random --seq example/36bp_20.txt -o example/output -e 2
	
	Output on the screen:

	20 sequences are read in 0.01. (0 discarded) [Mem:0.01 M]
	-----------------------------------------------------------------------------------------------------------
	|     Genome Name |    Loading Time |    Mapping Time | Memory Usage(M) |  Total Mappings    Mapped reads |
	-----------------------------------------------------------------------------------------------------------
	|     chr1_random |            0.45 |            0.04 |          147.39 |            2530              13 |
	-----------------------------------------------------------------------------------------------------------
	     Total:            0.46              0.04

	Total Time:                         0.51
	Total No. of Reads:                   20
	Total No. of Mappings:              2530
	Avg No. of locations verified:         0




	B. Paired-end Reads

	$>./drFAST --search example/chr1_random --seq example/36bp_20.txt -e 2 --mp --min 50 --max 200 --discordant-vh -o example/output

	Output on the screen:

	20 sequences are read in 0.00. (0 discarded) [Mem:0.01 M]
	-----------------------------------------------------------------------------------------------------------
	|     Genome Name |    Loading Time |    Mapping Time | Memory Usage(M) |  Total Mappings    Mapped reads |
	-----------------------------------------------------------------------------------------------------------
	|     chr1_random |            0.41 |            0.24 |          439.84 |               0               0 |
	10
	-----------------------------------------------------------------------------------------------------------
		     Total:            0.41              0.24

	Post Processing Time:               0.02 
	Total Time:                         0.66
	Total No. of Reads:                   20
	Total No. of Mappings:                 0
	Avg No. of locations verified:         0


Upcoming features
=================

	1 - OEA files contains the mapping location/ in SAM format
	2 - Auto-detect FASTQ offset (33 vs 64) and scale to 33 if 64-based fastq was the input


Changes to previous version
===========================

No need to run the Driver to transfer the genome from base (letter space) to  color space, this step is incorporated in indexing step.

Additional Information
======================

  This version of drFAST is based on the mrsFAST code developed by Faraz Hach (fhach AT cs DOT sfu DOT ca),
  and its predecessor mrFAST developed by (Fereydoun Hormozdiari and Can Alkan)
 
  Copyright (c) <2009 - 2010>, University of Washington, Simon Fraser University
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:
 
  Redistributions of source code must retain the above copyright notice, this list
  of conditions and the following disclaimer.
  - Redistributions in binary form must reproduce the above copyright notice, this
    list of conditions and the following disclaimer in the documentation and/or other
    materials provided with the distribution.
  - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
    used to endorse or promote products derived from this software without specific
    prior written permission.
 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
        Farhad Hormozdiari
        Faraz Hach
	Can Alkan
  Emails: 
        farhadh AT u DOT washington DOT edu / farhad DOT hormozdiari AT gmail DOT com
    	fhach AT cs DOT sfu DOT ca
 	calkan AT gs DOT washington DOT edu


