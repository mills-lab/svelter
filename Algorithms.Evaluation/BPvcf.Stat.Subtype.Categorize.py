#!/usr/bin/python

#!Python
#command='BPvcf.Stat.Subtype.Categorize.py -p /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.comp/SVelter'
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
opts,args=getopt.getopt(sys.argv[1:],'p:x:',[])
dict_opts=dict(opts)
if not dict_opts['-p'][-1]=='/':
	dict_opts['-p']+='/'

all_stat_hash={}
for k1 in os.listdir(dict_opts['-p']):
	if k1.split('.')[-1]=='stat' and not 'Integrated' in k1:
		fin=open(dict_opts['-p']+k1)
		for line in fin:
			pin=line.strip().split()
			if not '_'.join(pin[-1].split('_')[-2:]) in all_stat_hash.keys():
				all_stat_hash['_'.join(pin[-1].split('_')[-2:])]={}
			if not k1 in all_stat_hash['_'.join(pin[-1].split('_')[-2:])].keys():
				all_stat_hash['_'.join(pin[-1].split('_')[-2:])][k1]=[[],[]]
			all_stat_hash['_'.join(pin[-1].split('_')[-2:])][k1][0].append(float(pin[0])/float(pin[1]))
			if not pin[2]=='0':
				all_stat_hash['_'.join(pin[-1].split('_')[-2:])][k1][1].append(float(pin[0])/float(pin[2]))
			else:
				all_stat_hash['_'.join(pin[-1].split('_')[-2:])][k1][1].append(1)
		fin.close()

fo=open(dict_opts['-p']+'BPLink.Integrated.Categorized'+'.stat','w')
for k1 in sorted(all_stat_hash.keys()):
	for k2 in sorted(all_stat_hash[k1].keys()):
		print >>fo, ' '.join([str(i) for i in [k1,numpy.mean(all_stat_hash[k1][k2][0]),numpy.mean(all_stat_hash[k1][k2][1]),k2]])

fo.close()









