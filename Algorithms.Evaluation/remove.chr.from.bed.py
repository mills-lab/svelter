#!/usr/bin/env python

#!python
#command='remove.chr.from.bed.py -i input.bed -o output.bed'
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
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=dict_opts['-o']
info_list=[]
fin=open(filein)
for line in fin:
	pin=line.strip().split()
	info_list.append([pin[0].replace('chr','')]+pin[1:])

fin.close()
fo=open(fileout,'w')
for x in info_list:
	print >>fo, '\t'.join(x)

fo.close()



