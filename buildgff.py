#!/usr/bin/python
from Bio import SeqIO
from Bio import Seq
import sys
from sys import argv
out=open("%s.gff"%argv[1],"w")
out.write("##gff-version 3\n")
out.write('##date Wed Aug 21 15:51:21 2013\n')
out.write('##source gbrowse gbgff gff3 dumper\n')
out.write('##sequence-region Contig7580.0:1..66022\n')
coordinates=open("%s"%argv[1],"r")#blast output file
j=0
for line in coordinates:
	MICcontig=line.split()[1]
	j+=1
	MACcontig=line.split()[0]
	t_start=int(line.split()[6])
	t_stop=int(line.split()[7])
	q_start=int(line.split()[8])
	q_stop=int(line.split()[9])
	evalue=line.split()[10]
	if t_stop < t_start:
		orientation='-'
	else:
		orientation='+'
	out.write("%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=test\n"%(MACcontig,MICcontig,MICcontig,t_start,t_stop,orientation))
out.close()
coordinates.close()
