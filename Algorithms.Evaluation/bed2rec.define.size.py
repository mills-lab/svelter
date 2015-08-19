#!/usr/bin/env python

#!python
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

import os
opts,args=getopt.getopt(sys.argv[1:],'i:',['size=','bam=','ppre=','sv='])
dict_opts=dict(opts)
k1=dict_opts['-i']
num_rec=int(dict_opts['--size'])
test=[]
fin=open(k1)
for line in fin:
	pin=line.strip().split()
	if not pin==[]:
			if pin[-1]=='homo':
		       		ref='a/a'
		       		alt='/'
			else:
		   		ref='a/a'
		   		alt='/a'
			test.append([ref,alt]+pin[:3])
fin.close()
test2=[[]]
for x in test:
		if len(test2[-1])<num_rec:
			test2[-1].append(x)
		else:
			test2.append([x])
for rec1 in range(len(test2)):
		fo=open('.'.join(k1.split('.')[:-1]+[str(rec1),'rec']),'w')
		for rec2 in test2[rec1]:
			print >>fo, '\t'.join(rec2)
		fo.close()

