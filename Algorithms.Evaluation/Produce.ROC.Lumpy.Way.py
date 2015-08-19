#!/usr/bin/python

#!Python
#command='Produce.ROC.Lumpy.Way.py --path_in /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.homo/Delly/bed_files/ --path_ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.homo/sv_new.rec/bed_files --appdix .Mappable.TRAFree.min100.max1000000000.bed'
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
from scipy.stats import scoreatpercentile
opts,args=getopt.getopt(sys.argv[1:],'i:r:',['appdix=','path_in=','path_ref=','minIL=','maxIL='])
dict_opts=dict(opts)
ref_path=dict_opts['--path_ref']
if not ref_path[-1]=='/':
	ref_path+='/'

in_path=dict_opts['--path_in']
if not in_path[-1]=='/':
	in_path+='/'

in_hash=[]
for k1 in os.listdir(in_path):
	if dict_opts['--appdix'] in k1:
		in_hash.append(in_path+k1)

out_stat={}
for k1 in in_hash:
	if os.path.isfile(ref_path+k1.split('/')[-1]):
		file1=k1
		file2=ref_path+k1.split('/')[-1]
		header=file1.split('/')[-1].split('.')[0]
		SV1=file1.split('/')[-1].split('.')[1]
		if not header in out_stat.keys():
			out_stat[header]={}
		if not SV1 in out_stat[header].keys():
			out_stat[header][SV1]={}
		fin1=open(file1)
		data_hash={}
		for line in fin1:
			pin1=line.strip().split()
			if not pin1[0] in data_hash.keys():
				data_hash[pin1[0]]=[]
			data_hash[pin1[0]].append([int(i) for i in pin1[1:]])
		fin1.close()
		fin2=open(file2)
		ref_hash={}
		for line in fin2:
			pin1=line.strip().split()
			if not pin1[0] in ref_hash.keys():
				ref_hash[pin1[0]]=[]
			ref_hash[pin1[0]].append([int(i) for i in pin1[1:]])
		fin2.close()
		for ka in data_hash.keys():
			if ka in ref_hash.keys():
				out_stat[header][SV1][ka]=[]
				flag2=0
				for kb in ref_hash[ka]:
					flag1=0
					for kc in data_hash[ka]:
						if abs(kc[0]-kb[0])<50 and abs(kc[1]-kb[1])<50:
							flag1+=1
							break
					flag2+=flag1
				out_stat[header][SV1][ka]=[flag2,len(ref_hash[ka]),len(data_hash[ka])]

fo=open(in_path+'ROC.Lumpy.Way','w')
for k1 in sorted(out_stat.keys()):
	for k2 in sorted(out_stat[k1].keys()):
		for k3 in out_stat[k1][k2].keys():
			print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+out_stat[k1][k2][k3]])

fo.close()


