#!/usr/bin/python

#!python
#command='TRA.comp.stat.integrate.py -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/TRA.Compare.Report/ --appdix .compare.Stat'
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
opts,args=getopt.getopt(sys.argv[1:],'p:',['appdix=','ref_genome='])
dict_opts=dict(opts)
if not dict_opts['-p'][-1]=='/':
  dict_opts['-p']+='/'

#fo=open(dict_opts['-p']+'Stat.Integrated'+dict_opts['--appdix'],'w')
for k1 in os.listdir(dict_opts['-p']):
  if dict_opts['--appdix'] in k1:
    fin=open(dict_opts['-p']+k1)
    fo=open(dict_opts['-p']+k1+'.Integrated','w')
    test_hash={}
    for line in fin:
      pin=line.strip().split()
      if not '='.join(sorted(pin[4:])) in test_hash.keys():
        test_hash['='.join(sorted(pin[4:]))]=[0, pin[:3]]
      test_hash['='.join(sorted(pin[4:]))][0]+=1
      if not pin[:3]==test_hash['='.join(sorted(pin[4:]))][-1]:
        test_hash['='.join(sorted(pin[4:]))].append(pin[:3])
    fin.close()
    test_order={}
    for x in test_hash.keys():
      if not test_hash[x][0] in test_order.keys():
        test_order[test_hash[x][0]]=[]
      test_order[test_hash[x][0]].append(x)
    for x in sorted(test_order.keys())[::-1]:
      for y in test_order[x]:
        print >>fo, ' '.join([str(i) for i in [test_hash[y][0]]+test_hash[y][1]+y.split('=')])
    fo.close()





