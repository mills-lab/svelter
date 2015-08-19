#!/usr/bin/python

#!python
#command='TRA.pick.out.py -i /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/Delly/Delly_het_RD_10_DEL.DEL.Mappable.bed --ref /nfs/remills-scratch/datasets/Simulation.Xuefang/Simulate.het/sv.rec/het.het.TRA.rec'
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
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref='])
dict_opts=dict(opts)
ref_hash={}
fin=open(dict_opts['--ref'])
for line in fin:
	pin=line.strip().split()
	if not pin[0] in ref_hash.keys():
		ref_hash[pin[0]]={}
	if not int(pin[1]) in ref_hash[pin[0]].keys():
		ref_hash[pin[0]][int(pin[1])]=[]
	if not [pin[0]]+[int(i) for i in pin[1:4]] in ref_hash[pin[0]][int(pin[1])]:
		ref_hash[pin[0]][int(pin[1])].append([pin[0]]+[int(i) for i in pin[1:4]])

fin.close()

in_hash={}
chromos=[]
fin=open(dict_opts['-i'])
for line in fin:
	pin=line.strip().split()
	if not pin[0] in in_hash.keys():
		in_hash[pin[0]]={}
		chromos.append(pin[0])
	if not int(pin[1]) in in_hash[pin[0]].keys():
		in_hash[pin[0]][int(pin[1])]=[]
	in_hash[pin[0]][int(pin[1])].append([pin[0]]+[int(i) for i in pin[1:]])

fin.close()

def hash_to_list(in_hash):
	ref_list={}
	for k1 in in_hash.keys():
		ref_list[k1]=[]
		for k2 in sorted(in_hash[k1].keys()):
			for k3 in in_hash[k1][k2]:
				ref_list[k1].append(k3)
	return ref_list

def bp_check(list1,list2):
	#eg: list1=['1', 3615657, 3638623]; list2=['1', 3615658, 3638624, 3751340];
	flag=0
	rec=[]
	for x in list1[1:]:
		for y in list2[1:]:
			if abs(x-y)<51:
				flag+=1
				rec.append(list2.index(y))
	return [flag,rec]

def compare_list(in_list,ref_list):
	out={}
	out2=[]
	for k1 in in_list.keys():
		if k1 in ref_list.keys():
			out[k1]={}
			for x1 in ref_list[k1]:
				for x2 in in_list[k1]:
					if x2[-1]<x1[1]: continue
					elif x2[1]>x1[-1]:continue
					else:
						bp_result=bp_check(x2,x1)
						if bp_result[0]==2:
							if bp_result[1]==[1,2]:
								block='a'
							elif bp_result[1]==[2,3]:
								block='b'
							elif bp_result[1]==[1,3]:
								block='ab'
							if not '_'.join([str(i) for i in x1]) in out[k1].keys():
								out[k1]['_'.join([str(i) for i in x1])]=[]
							out[k1]['_'.join([str(i) for i in x1])].append(x2+[block])
							out2.append(x2)
	return [out,out2]		

in_list=hash_to_list(in_hash)
ref_list=hash_to_list(ref_hash)
comp_all=compare_list(in_list,ref_list)
comp_hash=comp_all[0]
comp_list=comp_all[1]
fo=open(dict_opts['-i'].replace('.bed','.TRAFree.bed'),'w')
for k1 in chromos:
	for k2 in in_list[k1]:
		if not k2 in comp_list:
			print >>fo, ' '.join([str(i) for i in k2])

fo.close()
fo2=open(dict_opts['-i'].replace('.bed','.TRARelated.bed'),'w')
for k1 in chromos:
	for k2 in comp_hash[k1].keys():
		for k3 in comp_hash[k1][k2]:
			print >>fo2, ' '.join([str(i) for i in [k2]+k3])

fo2.close()


