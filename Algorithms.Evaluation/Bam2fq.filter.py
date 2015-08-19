#!/usr/bin/env python

#!python
#command='Bam2fq.filter.py --bam /mnt/EXT/Mills-scratch2/Xuefang/NA12878.Pacbio/NA12878_S1.bam --sv temp.sv.txt'
#sys.argv=command.split()
import os
import sys
import getopt
import re
import pickle
import time
import datetime
import random
import numpy
import glob
import numpy as np
opts,args=getopt.getopt(sys.argv[1:],'i:',['bam=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['--sv'])
pin=fin.readline().strip().split()
fin.close()

ref_sv=pin[0]
alt_sv=pin[1]
chrom=pin[2]
bps=pin[2:]
read_hash={}
fin=os.popen('samtools view %s %s:%d-%d'%(dict_opts['--bam'],bps[0],int(bps[1])-500,int(bps[-1])+500))
for line in fin:
	pin=line.strip().split()
	if not pin[0] in read_hash.keys():
		read_hash[pin[0]]=[]
	read_hash[pin[0]].append([pin[9],pin[10]])

fin.close()
fout1=dict_opts['--bam'].replace('.bam','.'+'_'.join([bps[0],bps[1],bps[-1]])+'.1.fq')
fout2=dict_opts['--bam'].replace('.bam','.'+'_'.join([bps[0],bps[1],bps[-1]])+'.2.fq')
fo1=open(fout1,'w')
fo2=open(fout2,'w')
rec=0
samp_name=dict_opts['--bam'].split('/')[-1].replace('.bam','')
for k1 in read_hash.keys():
	if len(read_hash[k1])==2:
		rec+=1
		print >>fo1, ' '.join(['@'+samp_name+'.'+str(rec),k1,'length=101'])
		print >>fo2, ' '.join(['@'+samp_name+'.'+str(rec),k1,'length=101'])
		print >>fo1, read_hash[k1][0][0]
		print >>fo2, read_hash[k1][1][0]
		print >>fo1, ' '.join(['+'+samp_name+'.'+str(rec),k1,'length=101'])
		print >>fo2, ' '.join(['+'+samp_name+'.'+str(rec),k1,'length=101'])
		print >>fo1, read_hash[k1][0][1]
		print >>fo2, read_hash[k1][1][1]

fo1.close()
fo2.close()


