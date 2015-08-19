#!/usr/bin/env python

#!python
#command='Pacbio.Vali.Rsquare.Stat.Analysis.py -p ./'
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
for x in os.listdir(out_path):
    if x.split('.')[-1]=='rsquare':
        stats=[[],[],[]]
        fin=open(out_path+x)
        for line in fin:
            pin=line.strip().split()
            stats[0].append(float(pin[0]))
            stats[1].append(float(pin[1]))
            stats[2].append(float(pin[2]))
        fin.close()
        if not stats[0]==[]:
            stats_mean.append([numpy.mean(stats[0]),numpy.mean(stats[1]),numpy.mean(stats[2])])
            stats_std.append([numpy.std(stats[0]),numpy.std(stats[1]),numpy.std(stats[2])])
            stats_name.append(x)

vali_path=out_path+'Validated/'
if not os.path.isdir(vali_path):
    os.system(r'''mkdir %s'''%(vali_path))

fo1=open(dict_opts['-p']+'Validated.cases','w')
fo2=open(dict_opts['-p']+'NonValidated.cases','w')
for x in range(len(stats_mean)):
    if numpy.max(stats_mean[x][1:2])-stats_mean[x][0]>0.1:
        print >>fo1,stats_name[x]
    else:
        print >>fo2,stats_name[x]

fo1.close()
fo2.close()








