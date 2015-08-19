#!/usr/bin/env python

#!python
#command='Integrate.compare.stat.py -p ./ -o ./Integrated.compare.stat'
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
opts,args=getopt.getopt(sys.argv[1:],'p:o:',['server=','ref=','sv='])
dict_opts=dict(opts)
import os
pathin=dict_opts['-p']
fileout=dict_opts['-o']
if not pathin[-1]=='/':
	pathin+='/'

data_hash={}
for k1 in os.listdir(pathin):
	if k1.split('.')[-1]=='stat' and k1.split('.')[-2]=='compare':
		fin=open(pathin+k1)
		pin=fin.readline().strip().split()
		for line in fin:
			pin=line.strip().split()
			if not pin[0].split('_')[-2] in data_hash.keys():
				data_hash[pin[0].split('_')[-2]]={}
			if not pin[0].split('_')[-1] in data_hash[pin[0].split('_')[-2]].keys():
				data_hash[pin[0].split('_')[-2]][pin[0].split('_')[-1]]=[]
			data_hash[pin[0].split('_')[-2]][pin[0].split('_')[-1]].append(pin)
		fin.close()

fo=open(fileout,'w')
for k1 in data_hash.keys():
	for k2 in data_hash[k1].keys():
		for k3 in data_hash[k1][k2]:
			if not float(k3[1])==float(k3[2])==float(k3[3])==float(k3[4])==float(1):
				print >>fo, '\t'.join(k3)

fo.close()


