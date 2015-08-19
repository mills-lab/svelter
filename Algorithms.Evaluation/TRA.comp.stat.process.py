#!/usr/bin/python

#!python
#command='TRA.comp.stat.process.py -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/TRA.Compare.Report/ --appdix .compare.Stat'
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

fo=open(dict_opts['-p']+'Stat.Integrated'+dict_opts['--appdix'],'w')
for k1 in os.listdir(dict_opts['-p']):
  if dict_opts['--appdix'] in k1:
    fin=open(dict_opts['-p']+k1)
    keyName=k1.replace(dict_opts['--appdix'],'')
    Sens=[]
    Specs=[]
    for line in fin:
      pin=line.strip().split()
      Sens.append(float(pin[0])/float(pin[1]))
      Specs.append(float(pin[0])/float(pin[2]))
    fin.close()
    print >>fo, ' '.join([str(i) for i in [keyName,numpy.mean(Sens),numpy.mean(Specs),len([i for i in Sens if i==1.0])]])

fo.close()
