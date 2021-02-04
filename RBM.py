#!/usr/bin/python
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
import sys
from sys import argv
#RBM.py <blast_output_1> <blast_output_2> NOTE: does not require -max_target_seqs 1, instead takes top hit for each from blast output. For RBMs prints both query and subject. For one-way matches just puts the query and 'one-way'; if query is not in file a '0'.
b1Handle=open("%s"%argv[1],"r")
b2Handle=open("%s"%argv[2],"r")
outHandle=open("%s_%s_RBM.txt"%(argv[1],argv[2]),'w')
b1Dict={}
b2Dict={}
for line in b1Handle:
	q1=line.split()[0]
	s1=line.split()[1]
	if b1Dict.has_key(q1):
		pass
	else:
		b1Dict[q1]=s1
for line in b2Handle:
        q2=line.split()[0]
        s2=line.split()[1]
        if b2Dict.has_key(q2):
                pass
        else:
                b2Dict[q2]=s2
keys1=b1Dict.keys()
x=0
for k1 in keys1:
	if b2Dict.has_key(b1Dict[k1]):
		if b2Dict[b1Dict[k1]]==k1:
			x+=1
			outHandle.write('%s\t%s\n'%(k1,b1Dict[k1]))
			del b2Dict[b1Dict[k1]]#remove from dicts
			del b1Dict[k1]#remove
print '%s RBMs printed to file.'%x			
				
b1Handle.close()
b2Handle.close()
outHandle.close()
