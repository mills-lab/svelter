#!/usr/bin/python

#!Python
#command='Lumpy.bedpe2bedpe.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/Lumpy/bedpe/Simulate_het_RD30_pesr.bedpe'
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
opts,args=getopt.getopt(sys.argv[1:],'i:r:',['ref=','start=','end=','minIL=','maxIL='])
dict_opts=dict(opts)
lumpy_in=dict_opts['-i']
lum_hash={}
fin=open(lumpy_in)
for line in fin:
	pin=line.strip().split()
	if not pin[10].split(':')[1] in lum_hash.keys():
		lum_hash[pin[10].split(':')[1]]={}
	if not pin[0] in lum_hash[pin[10].split(':')[1]].keys():
		lum_hash[pin[10].split(':')[1]][pin[0]]=[]
	lum_hash[pin[10].split(':')[1]][pin[0]].append(pin[:6])

fin.close()
for k1 in lum_hash.keys():
	k2=k1[:3]
	fo=open(lumpy_in.replace('.bedpe','.'+k2+'.bedpe'),'w')
	for k3 in lum_hash[k1].keys():
		for k4 in lum_hash[k1][k3]:
			print >>fo, ' '.join(k4)
	fo.close()








