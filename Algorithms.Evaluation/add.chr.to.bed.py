#!/usr/bin/env python

#!python
#command='add.chr.to.bed.py -i input.bed'
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
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
filein=dict_opts['-i']
info_list=[]
fin=open(filein)
for line in fin:
	pin=line.strip().split()
	info_list.append(['chr'+pin[0]]+pin[1:])

fin.close()
fo=open(filein,'w')
for x in info_list:
	print >>fo, '\t'.join(x)

fo.close()

