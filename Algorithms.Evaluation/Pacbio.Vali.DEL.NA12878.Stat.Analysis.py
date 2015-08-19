#!/usr/bin/env python

#!python
#command='Pacbio.Vali.DEL.Stat.Analysis.version2.py -p ./'
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

opts,args=getopt.getopt(sys.argv[1:],'p:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
out_path=dict_opts['-p']
if not out_path[-1]=='/':
    out_path+='/'

stats_mean=[]
stats_std=[]
stats_name=[]
total_name=[]
for x in os.listdir(out_path):
        if x.split('.')[-1]=='rsquare':
                total_name.append(x)
                stats1=0
                stats2=0
                fin=open(out_path+x)
                for line in fin:
                    pin=line.strip().split()
                    if max([float(pin[1]),float(pin[2])])>float(pin[0]):
                        stats1+=1
                    else:
                        stats2+=1
                fin.close()
                if stats1+stats2>0:
                    stats_mean.append(float(stats1)/float(stats1+stats2))
                    stats_name.append(x)

fo1=open(dict_opts['-p']+'Validated.cases','w')
fo2=open(dict_opts['-p']+'NonValidated.cases','w')
for x in range(len(stats_mean)):
    if stats_mean[x]>0.2:
        print >>fo1,'\t'.join([str(i) for i in [stats_name[x],stats_mean[x]]])
    else:
        print >>fo2,'\t'.join([str(i) for i in [stats_name[x],stats_mean[x]]])

for x in total_name:
    if not x in stats_name:
        print >>fo2,x

fo1.close()
fo2.close()



