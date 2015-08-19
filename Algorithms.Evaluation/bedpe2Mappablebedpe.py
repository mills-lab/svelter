#!/usr/bin/python

#!Python
#command='bedpe2Mappablebedpe.py -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.het/Lumpy/bedpe/Simulate_het_RD30_pesr.DEL.bedpe --ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.Mappable'
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
from scipy.stats import scoreatpercentile
opts,args=getopt.getopt(sys.argv[1:],'i:r:',['ref=','path_in=','path_ref=','minIL=','maxIL='])
dict_opts=dict(opts)
ref_path='/'.join(dict_opts['--ref'].split('/')[:-1])
ref_prefix=dict_opts['--ref'].split('/')[-1]
ref_files=[]
for k1 in os.listdir(ref_path):
	if ref_prefix in k1 and k1.split('.')[-1]=='bed':
		ref_files.append(ref_path+'/'+k1)

ref_hash={}
for k1 in ref_files:
	chrom=k1.split('/')[-1].replace('.bed','').replace(ref_prefix+'.','')
	ref_hash[chrom]=[]
	fin=open(k1)
	for line in fin:
		pin=line.strip().split()
		ref_hash[chrom].append([int(pin[1]),int(pin[2])])
	fin.close()

file_in=dict_opts['-i']
file_out=file_in.replace('.bedpe','.Mappable.bedpe')
fin=open(file_in)
in_hash={}
chroms=[]
for line in fin:
	pin=line.strip().split()
	if not pin[0] in in_hash.keys():
		in_hash[pin[0]]=[]
		chroms.append(pin[0])
	in_hash[pin[0]].append([int(pin[1]),int(pin[2]),int(pin[4]),int(pin[5])])

fin.close()

fo=open(file_out,'w')
for k1 in chroms:
	if k1 in in_hash.keys() and k1 in ref_hash.keys():
		for k2 in in_hash[k1]:
			flag=0
			for k3 in ref_hash[k1]:
				if k3[0]-1<k2[0] and k3[1]+1>k2[3]:
					flag+=1
			if not flag==0:
				print >>fo, ' '.join([str(i) for i in [k1]+k2])

fo.close()




