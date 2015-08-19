#!/usr/bin/env python

#!python
#command='Bed.Intersect.py -f 0.5 -a SAMN02744161.Lumpy.DEL.bed -b SAMN02744161.DELLY.DEL.bed -o SAMN02744161.Lumpy.DELLY.DEL.bed'
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

opts,args=getopt.getopt(sys.argv[1:],'a:b:f:o:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
filein1=dict_opts['-a']
filein2=dict_opts['-b']
def read_in_hash(filein):
        out={}
        fin=open(filein)
        for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                        out[pin[0]]=[]
                out[pin[0]].append([int(i) for i in pin[1:3]]+pin[3:])
        fin.close()
        return out

def hash_compare(filein1,filein2,reciprocal):
        out={}
        out2={}
        for k1 in filein1.keys():
                if k1 in filein2.keys():
                        out[k1]=[]
			out2[k1]=[]
                        for k2 in filein1[k1]:
                                flag=0
                                for k3 in filein2[k1]:
                                        if k3[1]<k2[0]:continue
                                        elif k3[0]>k2[1]:continue
                                        else:
                                                if not float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max([k2[1]-k2[0],k3[1]-k2[0]])) <reciprocal:
                                                        flag+=1
                                if flag>0:
                                        out[k1].append(k2)
                                        out2[k1].append(k3)
        return [out,out2]

def hash_write(hash_intersect,out_file):
        fo=open(out_file,'w')
        for k1 in hash_intersect.keys():
                for k2 in hash_intersect[k1]:
                        print >>fo, '\t'.join([str(i) for i in [k1]+k2])
        fo.close()

hash1=read_in_hash(filein1)
hash2=read_in_hash(filein2)
reciprocal=float(dict_opts['-f'])
hash_intersect=hash_compare(hash1,hash2,reciprocal)
hash_write(hash_intersect[0],dict_opts['-a']+dict_opts['-b'])
hash_write(hash_intersect[1],dict_opts['-b']+dict_opts['-a'])


