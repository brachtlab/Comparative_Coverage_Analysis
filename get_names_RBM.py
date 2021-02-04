#!/usr/bin/env python
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
import sys
from sys import argv
#get_names_RBM.py <contig_list.txt><RBM.txt>
contigHandle=open("%s"%argv[1],"r")
RBMHandle=open("%s"%argv[2],"r")
outHandle=open("%s_RBM.txt"%argv[1],"w")
contigDict={}
for line in contigHandle:
	contigDict[line.split()[0]]=1
for line in RBMHandle:
	q1=line.split()[0]
	s1=line.split()[1]
	if contigDict.has_key(q1):
		outHandle.write('%s\t%s\n'%(q1, s1))
		print "found it!"
		del contigDict[q1]
print "not found"
print contigDict
contigHandle.close()
RBMHandle.close()
outHandle.close()
