#!/usr/bin/env python

#!python
#command='Order.NonValidated.cases.py -i NonValidated.Cases'
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
data_hash={}
extras=[]
fin=open(filein)
for line in fin:
    pin=line.strip().split()
    if len(pin)>1:
        if not float(pin[-1]) in data_hash.keys():
            data_hash[float(pin[-1])]=[]
        data_hash[float(pin[-1])].append(pin[0])
    else:
        extras.append(pin[0])

fin.close()
fo=open(filein,'w')
for x1 in sorted(data_hash.keys())[::-1]:
    for x2 in data_hash[x1]:
        print >>fo, '\t'.join([x2,str(x1)])

for x1 in extras:
            print >>fo, x1

fo.close()

