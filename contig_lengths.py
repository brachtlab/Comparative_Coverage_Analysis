#!/usr/bin/python
from Bio import SeqIO
import sys
from sys import argv
outHandle=open("%s_lengths.csv"%argv[1],'w')
genome_handle=open("%s"%argv[1],"r")
outHandle.write("chr,length\n")
for seq_record in SeqIO.parse(genome_handle,"fasta"):
	seqlength=len(seq_record.seq)
	ident=seq_record.id
	outHandle.write("%s,%s\n"%(ident,seqlength))		
genome_handle.close()
outHandle.close()
