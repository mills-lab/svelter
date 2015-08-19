#!/usr/bin/env python

#!python
#command='Extract.csv.from.rec.py --qc 0 -i /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/SVelter/NA12878.all.SV.rec'
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

opts,args=getopt.getopt(sys.argv[1:],'i:',['qc=','bam=','ppre=','sv='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fout=filein.replace('.rec','.csv.rec')
in_hash={}
fin=open(filein)
qc=float(dict_opts['--qc'])
for line in fin:
    pin=line.strip().split()
    if float(pin[-1])>qc:
        if not pin[0] in in_hash.keys():
            in_hash[pin[0]]={}
        if not pin[1] in in_hash[pin[0]].keys():
            in_hash[pin[0]][pin[1]]=[]
        in_hash[pin[0]][pin[1]].append(pin[2:-1])

fin.close()

def let2list(k1):
    out=[]
    for x in k1.split('/'):
        out.append([])
        for y in x:
            if not y=='^':
                out[-1].append(y)
            else:
                out[-1][-1]+=y
    return out

def combine_list(k1):
    ka=let2list(k1)
    out=[]
    for x in ka:
        out.append([])
        for y in x:
            if out[-1]==[]:
                out[-1].append([y])
            else:
                if not '^' in y and not '^' in out[-1][-1][-1]:
                    if ord(y[0])-ord(out[-1][-1][-1][0])==1:
                        out[-1][-1].append(y)
                    else:
                        out[-1].append([y])
                elif '^' in y and'^' in out[-1][-1][-1]:
                    if ord(y[0])-ord(out[-1][-1][-1][0])==-1:
                        out[-1][-1].append(y)
                    else:
                        out[-1].append([y])
                else:
                        out[-1].append([y])
    return out

def merge_list(k1):
    ka=combine_list(k1)
    out=[]
    for x in ka:
        out.append([])
        for y in x:
            if not '^' in y[0]:
                out[-1].append(''.join(y))
            else:
                out[-1].append(''.join([z[0] for z in y[::-1]]))
                out[-1][-1]+='^'
    return out

def simple_del_decide(k1,k2):
    flag=0
    for x in k2.split('/'):
        for y in range(len(x)-1):
            if not ord(x[y+1])-ord(x[y])>0:
                flag+=1
    if flag==0:return True
    else: return False

def simple_inv_decide(k1,k2):
    ka=merge_list(k2)
    out=[]
    for x in ka:
        out.append([])
        for y in x:
            out[-1].append(y.replace('^',''))
    out2='/'.join([''.join(x) for x in out])
    if k1==out2: return True
    else: return False

def simple_dup_decide(k1,k2):
    out=[]
    for x in k2.split('/'):
        out.append([])
        for y in x:
            if out[-1]==[]:
                out[-1].append(y)
            else:
                if not y==out[-1][-1]:
                    out[-1].append(y)
    k3='/'.join([''.join(x) for x in out])
    if k1==k3: return True
    else: return False

def main(fout):
    fo=open(fout,'w')
    for k1 in in_hash.keys():
        if not k1=='a/a':
            for k2 in in_hash[k1].keys():
                if not k1==k2:
                    if not simple_del_decide(k1,k2)==True:
                        if not simple_inv_decide(k1,k2)==True:
                            if not simple_dup_decide(k1,k2)==True:
                                for k3 in in_hash[k1][k2]:
                                    print >>fo, ' '.join([str(i) for i in [k1,k2]+k3])
    fo.close()

main(fout)









