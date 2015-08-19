#!/usr/bin/env python

#!python
#command='bed2rec.py -i /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Bed_Files_SV/SAMN02744161.overlap.DEL.bed --sv homo_del'
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
opts,args=getopt.getopt(sys.argv[1:],'i:',['sv=','bam=','ppre=','sv='])
dict_opts=dict(opts)

def sv_geno_readin():
	sv_in=dict_opts['--sv']
	if '_' in sv_in:
		geno=sv_in.split('_')[0]
		sv_type=sv_in.split('_')[1]
	else:
		sv_type=sv_in
		geno=0
	return [sv_type,geno]

def ref_alt_produce():
	sv_geno_info=sv_geno_readin()
	ref='a'
	if sv_geno_info[0]=='del':
		alt=''
	elif sv_geno_info[0]=='inv':
		alt='a^'
	elif sv_geno_info[0]=='dup':
		alt='aa'
	return alt

ref='a'
alt=ref_alt_produce()
sv_geno_info=sv_geno_readin()
genotype=sv_geno_info[1]

fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.bed','.rec'),'w')
for line in fin:
	pin=line.strip().split()
	if genotype==0:
		new_geno=pin[-1]
	else:
		new_geno=genotype
	if new_geno=='homo':
		alt_geno='/'.join([alt,alt])
	elif new_geno=='het':
		alt_geno='/'.join([ref,alt])
	ref_geno='/'.join([ref,ref])
	print >>fo,'\t'.join([ref_geno,alt_geno]+pin[:3])

fin.close()
fo.close()
 

