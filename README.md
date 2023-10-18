# RNA-to-DNA-Reverse-Transcription-Analysis
The primary objective of this is to develop a program that reverse transcribes RNA to DNA, computes nucleotide frequencies, and identifies differences in dinucleotide and trinucleotide composition between two RNA sequences. 
Key Features:

RNA to DNA Transcription: Implement a program to reverse transcribe RNA sequences into DNA for two input RNA files (RNA-sequence1.fna and RNA-sequence2.fna).

Nucleotide Frequency Analysis: Calculate di-nucleotide and tri-nucleotide frequencies for the original RNA sequences, presenting results in a specified three-column format, including absolute frequency and percentage.

Comparison of Nucleotide Composition: Utilize the program to identify differences in the percentage abundance of dinucleotide and trinucleotide compositions between the two RNA sequences. Detect pairs that differ by a significant margin (threefold or 3X times).

Data Files: The necessary input RNA sequence files (RNA-sequence1.fna and RNA-sequence2.fna) are available for download within the repository.

Required packages/Libraries: 
1, Pandas
- To install Pandas package, type the following command in command prompt
	py -m pip install "pandas"
- To import the package to start using it, use the following command
	import pandas as pd

2, Biopython
- To install Biopython package, type the following command in command prompt
	py -m pip install "biopython"
- To import the package you need, use the following commands
	from Bio import SeqIO
	from Bio.Seq import Seq

Description: 
- The first defined function reads the sequences from the fasta file
- The sequences are then reverse transcribed and the forward strand is printed. 
- The code laters defines the di and tri-nucleotides and calculates the absolute frequency and percentage in the given sequences.
- The code later compares the percentage fold difference between the di and tri nucleotides between the two sequences.
