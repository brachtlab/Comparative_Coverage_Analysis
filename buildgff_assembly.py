#!/usr/bin/python
from Bio import SeqIO
from Bio import Seq
import sys
from sys import argv
out=open("%s.gff"%argv[1],"w")
out.write("##gff-version 3\n")
out.write("##date Mon May 08 18:30:00 2020\n")
out.write("##source gbrowse gbgff gff3 dumper\n")
out.write("##sequence-region Contig7580.0:1..66022\n")
contigs=open("%s"%argv[1],"r")#assembly file
for seq_record in SeqIO.parse(contigs,"fasta"):
	MACcontig=seq_record.id
	MICcontig=seq_record.id
	t_start=1
	t_stop=len(seq_record.seq)
	out.write("%s\tmaker\ttranscript\t%s\t%s\t.\t.\t.\tID=test\n"%(MACcontig,t_start,t_stop))
out.close()
contigs.close()
