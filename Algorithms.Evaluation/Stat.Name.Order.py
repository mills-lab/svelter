#!/usr/bin/python

#!python
#command='Stat.Name.Order.py -i Stat.Integrated.compare.Stat'
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
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','ref_genome='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
hash={}
for line in fin:
  pin=line.strip().split()
  hash[pin[0]]=pin

fin.close()
fo=open(dict_opts['-i'],'w')
for k1 in sorted(hash.keys()):
  print >>fo, '\t'.join(hash[k1])

fo.close()
