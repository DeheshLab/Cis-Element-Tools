#!/usr/bin/perl
#v4.cisfinder.pl by gkbenn
#the purpose of this program is to scan a list of input
#sequences (ie promoters) for entered motifs (ie a cis-element)
#the program produces 3 output files. The first provides a list of all genes with the element, 
#their sequences, and locations. The second provides a list of all genes with the searched element. 
#The third provides basic statistics about the search.

#To use this program, call it in the terminal, followed by the file with your sequence, 
#the length of each promoter you would like to search, and the element(s) to be searched.
# Example: perl v4cisfinder test.sequence.file.txt 2000 CCCGAA TTTTAA
#This searches the first 2000bp (from TSS) of each promoter in "test.sequence.file.txt" for CCCGAA and TTTTAA

#the main change with version 4 is the added ability to scan multiple motifs
#at once. This allows for searching for degenerate elements.

#IMPORTANT note: in order to screen the last sequence in the file, you
#must add "ENDLINE" as the final line in the sequence file

use strict; use warnings;

#my ($seqfile, $upstream) = @ARGV;
my $entry = @ARGV;
my $seqfile = $ARGV[0];
my $upstream = $ARGV[1];

open(OUT1, ">v4cis.$upstream.detail.output.txt");
open(OUT2, ">v4cis.$upstream.positive.list.txt");
open(OUT3, ">v4cis.$upstream.report.txt");

open(IN, "<$seqfile");

my $x = 0;						#variable used to add $line to $prom in 1st run through loop
my $i = 0;						#variable used to set $ATGnum = $line in 1st run through loop
my $y = 0;						#variable used to report each ATGnum only once
my $ATGnum = ();				#variable used to hold the ATG identifier
my $prom = ();					#string used to hold the promoter sequence
my $A = 0;						#variable to hold count of As
my $T = 0;						#variable to hold count of Ts
my $C = 0;						#variable to hold count of Cs
my $G = 0;						#variable to hold count of Gs
my $length = 0;					#variable to hold total length
my $numproms = 0;				#variable to hold number of promoters looked at
my $posproms = 0;				#variable to hold number of positives
my $numcis = 0;					#variable to hold number of elements found

while (<IN>) {
	chomp;
	my $line = $_;
	$line = uc($line);
	
### LOOP1:
	if (($line =~ m/^>/) or ($line =~ m/ENDLINE/)) {		#have we reached an info line/end of the file?
		if ($i == 0) {										#1st line of file enters this loop
			($ATGnum) = $line =~ m/^>(AT[\dMC]G\d{5})/;			#$ATGnum gets initial assignment here
			$i++;											#as this is 1st line, we now leave to build up $prom
		}
#look for motif in this loop		
		else {
#this section cuts the sequence down to the user-specified length
			$prom = substr($prom,-$upstream);
#this section of code collects information for later statistical analyis
			my $a = $prom =~ tr/A/A/;
			$A = $A + $a;
			my $t = $prom =~ tr/T/T/;
			$T = $T + $t;
			my $c = $prom =~ tr/C/C/;
			$C = $C + $c;
			my $g = $prom =~ tr/G/G/;
			$G = $G + $g;
			$length = $length + length($prom);
			$numproms++;			
#this section of code uses a sliding window to look for the specified motif
			my $n = 2;					#this is the count variable to pull the motifs out of the command line entry
			while ($n < $entry) {
				my $motif = $ARGV[$n];
				$motif = uc($motif);
				for (my $i = 0; $i < length($prom) - length($motif) + 1; $i++) {
					my $subseq = substr($prom, $i, length($motif));
					if ($subseq =~ m/$motif/) {
						my $dist = $i - length($prom);		#distance from translation start site
						if ($y == 0) {						#this prevents repetitive printing of $ATGnum
							print OUT1 ">$ATGnum\n";
							print OUT2 "$ATGnum\n";
							$y++;
							$posproms++;
						}
						print OUT1 "$motif\t$dist\n";
						$numcis++;					
					}
				}		
				$n++;		
			}
			($ATGnum) = $line =~ m/^>(AT[\dMC]G\d{5})/;	#set new value for $ATGnum as current line						
			$prom = ();								#empty $prom so it can take in the next sequence
			$x = 0;									#set $x to 0 so that $prom can initiate properly in LOOP2
			$y = 0;
		}
	}
	
### LOOP2:	
	else {
		$line =~ s/\s//;
		if ($x == 0) {
			$prom = $line;			#set initial value of $prom on first time through
			$x++;
		}
		else {
			$prom = join("","$prom","$line");	#adds $line to end of $prom
		}
	}
}

close IN;
close OUT1;
close OUT2;

#this section prints out the report file

my $pA = $A/$length;
my $pT = $T/$length;
my $pC = $C/$length;
my $pG = $G/$length;

print OUT3 "$seqfile cisfinder report file\n\n";
print OUT3	"Entered motifs found $numcis times in $posproms out of $numproms sequences searched\n\n";
print OUT3	"Nucleotide distribution for searched sequences:\n";
print OUT3	"A = $pA\n";
print OUT3	"T = $pT\n";
print OUT3	"C = $pC\n";
print OUT3	"G = $pG\n";

close OUT3;

	
	
