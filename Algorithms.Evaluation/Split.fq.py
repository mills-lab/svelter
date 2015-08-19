#!/usr/bin/python

#!python
#command='python Split.fq.py -i het.ref1.RD10.1.fq -n 200291000'
#sys.argv=command.split()[1:]
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
opts,args=getopt.getopt(sys.argv[1:],'i:n:',['ref='])
dict_opts=dict(opts)
max_num=int(dict_opts['-n'])
file_num=1
rec_num=0

fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.fq','.'+str(file_num)+'.fq'),'w')
for line in fin:
        pin=line.strip()
        rec_num+=1
        if rec_num<max_num:
                print >>fo, pin
                rec_num+=1
        else:
                fo.close()
                file_num+=1
                rec_num=0
                fo=open(dict_opts['-i'].replace('.fq','.'+str(file_num)+'.fq'),'w')
                print >>fo, pin
                rec_num+=1

fo.close()
fin.close()


