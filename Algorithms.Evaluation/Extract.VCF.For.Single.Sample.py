#!/usr/bin/python

#!python
#command='Extract.VCF.For.Single.Sample.py -i /nfs/remills-data/xuefzhao/projects.axiom/ALL.wgs.mergedSV.v3.20130502.svs.genotypes.vcf -s NA12878 '
#sys.argv=command.split()
import os
import sys
import getopt
import numpy
import re
import random
import pickle
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'i:s:',['ref=','ref_genome='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=dict_opts['-i'].replace('.vcf','.'+dict_opts['-s']+'.vcf')
fin=open(filein)
header=[]
def genotype_decide(pin):
  for x in pin[8].split(':'):
    if x=='GT':
      x_index=pin[8].index(x)
  gt_out='Error'
  gt=pin[9].split(':')[x_index]
  if '|' in gt:
    if sorted(gt.split('|'))==['0','0']: gt_out=gt_out
    elif sorted(gt.split('|'))==['0','1']: gt_out='het'
    elif sorted(gt.split('|'))==['1','1']: gt_out='homo'
  elif '/' in gt:
    if sorted(gt.split('/'))==['0','0']: gt_out=gt_out
    elif sorted(gt.split('/'))==['0','1']: gt_out='het'
    elif sorted(gt.split('/'))==['1','1']: gt_out='homo'
  return gt_out

fo_flag=0
for line in fin:
  pin=line.strip().split()
  if pin[0][:2]=='##':
    header.append(pin)
  elif pin[0][0]=='#':
    if not dict_opts['-s'] in pin:
      break
    else:
      samp_index=pin.index(dict_opts['-s'])
      fo=open(fileout,'w')
      fo_flag=1
      for k1 in header:
        print>>fo, '\t'.join(k1)
      print >>fo, '\t'.join(pin[:9]+[pin[samp_index]])
  else:
    if not genotype_decide(pin[:9]+[pin[samp_index]])=='Error':
      print >>fo, '\t'.join(pin[:9]+[pin[samp_index]])

fin.close()
if fo_flag==1:
  fo.close()




