#!/usr/bin/env python

#!python
#command='Compare.Stat.Integrate.py -p ./'
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
import scipy
from scipy import stats

opts,args=getopt.getopt(sys.argv[1:],'p:',['ref=','bam=','ppre=','sv=','Path_Lumpy=','Path_Delly='])
dict_opts=dict(opts)
pathin=dict_opts['-p']
if not pathin[-1]=='/':
	pathin+='/'

data_in=[]
for k1 in os.listdir(pathin):
	if k1.split('.')[-1]=='stat' and k1.split('.')[-2]=='compare':
		fin=open(pathin+k1)
		fin.readline().strip().split()
		for line in fin:
			pin=line.strip().split()
			if pin[2]=='-1':
				pin[2]='1'
			if pin[3]=='-1':
				pin[3]='1'
			data_in.append(pin)
		fin.close()

fo=open(pathin+'All.Integrated.compare.stat','w')
for x in data_in:
	print >>fo, '\t'.join(x)

fo.close()

