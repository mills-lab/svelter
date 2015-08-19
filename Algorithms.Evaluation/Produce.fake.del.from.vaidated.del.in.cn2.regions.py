#!/usr/bin/env python

#!python
#Code to test specificity of PacVali code for simple del
#Stretegy: random put fake deletions into cn2 regions, then validate
#command='Produce.random.fake.del.file.from.validated.del.bed.py -i Validated.del.bed -o faked.del.bed'
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
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['server=','ref=','sv='])
dict_opts=dict(opts)

del_info_in=dict_opts['-i']
del_simu_out=dict_opts['-o']
cn2_region_in='/mnt/EXT/Mills-data/xuefzhao/SVelter/tools/CN2.hg19.bed'
cn2_hash={}
fin=open(cn2_region_in)
for line in fin:
        pin=line.strip().split()
        if not pin[0] in cn2_hash.keys():
                cn2_hash[pin[0]]=[]
        if int(pin[2])-int(pin[1])>2000:
                cn2_hash[pin[0]].append(pin[1:])

fin.close()

del_info=[]
fin=open(del_info_in)
for line in fin:
        pin=line.strip().split()
        del_info.append(int(pin[2])-int(pin[1]))

fin.close()
out_del=[]
for y in range(len(del_info)):
        print y
        x=del_info[y]
        chrom=random.choice(cn2_hash.keys())
        region=random.choice(cn2_hash[chrom])
        reg2=[int(i) for i in region]
        while True:
                if reg2[1]-reg2[0]-1100<x:
                        chrom=random.choice(cn2_hash.keys())
                        region=random.choice(cn2_hash[chrom])
                        reg2=[int(i) for i in region]
                else:
                        start=random.choice(range(reg2[0]+500,reg2[1]-500-x))
                        out_del.append([chrom,start,start+x])
                        break

fo=open(del_simu_out,'w')
for x in out_del:
        print >>fo, '\t'.join([str(i) for i in x])

fo.close()


