#!/usr/bin/env python

#!Python
#usage: Stat2Integrated.Categorized.stat.py -p input/path/of/*.stat --ref simu.SV.rec
#eg:BPvcf.Stat.Subtype.Categorize.py -p /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp.het/Stat.compare/stats/ --ref /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.comp/comp_het.rec/comp_het.SV.rec
#sys.argv=command.split()
#This command integrates the sens and spec of each simulated CSV, and report them in 'BPLink.Integrated.Categorized.stat' under input path
#Input path is the 
import os
import sys
import getopt
import numpy
import re
import random
import pickle
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'p:x:',['ref='])
dict_opts=dict(opts)
if not dict_opts['-p'][-1]=='/':
	dict_opts['-p']+='/'

def read_in_ref_hash():
	ref_hash={}
	fin=open(dict_opts['--ref'])
	for line in fin:
		pin=line.strip().split()
		key1='_'.join([pin[0],pin[1]])
		if not key1 in ref_hash.keys():
			ref_hash[key1]=[]
		ref_hash[key1].append(pin[-2:])
	fin.close()
	return ref_hash

def stat_info_orignal_structure_modify(pin):
	key1='_'.join([pin[3].split('_')[0],pin[3].split('_')[ord(pin[4][0])-97+1]])
	key2=ref_hash[key1]
	pin.append('_'.join(key2[0]))
	return pin


ref_hash=read_in_ref_hash()


all_stat_hash={}
for k1 in os.listdir(dict_opts['-p']):
	if k1.split('.')[-1]=='stat' and not 'Integrated' in k1:
		fin=open(dict_opts['-p']+k1)
		for line in fin:
			pin=line.strip().split()
			pin=stat_info_orignal_structure_modify(pin)
			if not pin[-1] in all_stat_hash.keys():
				all_stat_hash[pin[-1]]={}
			if not k1 in all_stat_hash[pin[-1]].keys():
				all_stat_hash[pin[-1]][k1]=[[],[]]
			all_stat_hash[pin[-1]][k1][0].append(float(pin[0])/float(pin[1]))
			if not pin[2]=='0':
				all_stat_hash[pin[-1]][k1][1].append(float(pin[0])/float(pin[2]))
			else:
				all_stat_hash[pin[-1]][k1][1].append(1)
		fin.close()

fo=open(dict_opts['-p']+'BPLink.Integrated.Categorized'+'.stat','w')
for k1 in sorted(all_stat_hash.keys()):
	for k2 in sorted(all_stat_hash[k1].keys()):
		print >>fo, ' '.join([str(i) for i in [k1,numpy.mean(all_stat_hash[k1][k2][0]),numpy.mean(all_stat_hash[k1][k2][1]),k2]])

fo.close()









