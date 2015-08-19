#!/usr/bin/env python

#!python
#command='Bed.Add.Score.py -a SVelter.Validated.cases.bed --pa /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/alt_bed/extra_SVs/SVelter.Stats/'
#sys.argv=command.split()
#script for validate large dels >3000, where not a lot of pacbio reads goes through entirely
import os
import sys
import getopt
import re
import numpy
opts,args=getopt.getopt(sys.argv[1:],'a:b:o:',['pa=','pb=','ppre=','sv='])
dict_opts=dict(opts)
patha=dict_opts['--pa']
if not patha[-1]=='/':
	patha+='/'

def file_name_readin(patha):
	file_hash_a={}
	for k1 in os.listdir(patha):
		if k1.split('.')[-1]=='rsquare':
			if not k1.split('.')[-2] in file_hash_a.keys():
				file_hash_a[k1.split('.')[-2]]=[]
			file_hash_a[k1.split('.')[-2]].append(patha+k1)
	return file_hash_a

def hash_readin(filea):
	out1={}
	fin=open(filea)
	for line in fin:
		pin=line.strip().split()
		if not pin[0] in out1.keys():
			out1[pin[0]]={}
		if not int(pin[1]) in out1[pin[0]].keys():
			out1[pin[0]][int(pin[1])]=[]
		if not int(pin[2]) in out1[pin[0]][int(pin[1])]:
			out1[pin[0]][int(pin[1])].append(int(pin[2]))
	fin.close()
	out2={}
	for k1 in out1.keys():
		out2[k1]=[]
		for k2 in sorted(out1[k1].keys()):
			for k3 in sorted(out1[k1][k2]):
				out2[k1].append([k2,k3])
	return out2

def rsquare_readin(file):
	fin=open(file)
	data=[[],[]]
	for line in fin:
		pin=line.strip().split()
		data_a=float(pin[0])
		data_b=max([float(pin[1]),float(pin[2])])
		if data_a<data_b:
			data[0].append(data_a)
			data[1].append(data_b)
	fin.close()
	return numpy.mean(data[1])-numpy.mean(data[0])

out=[]
file_hash_a=file_name_readin(patha)
hash_a=hash_readin(dict_opts['-a'])
for k1 in hash_a.keys():
	for k2 in range(len(hash_a[k1])):
			filea=file_hash_a['_'.join([str(i) for i in [k1]+hash_a[k1][k2]])]
			scorea=rsquare_readin(filea[0])
			hash_a[k1][k2]+=[scorea]

fo=open(dict_opts['-a']+'Score','w')
for k1 in hash_a.keys():
	for k2 in hash_a[k1]:
              #	print >>fo, '\t'.join([str(i) for i in [k1]+k2)
		print >>fo, '\t'.join([str(i) for i in [k1]+k2])

fo.close()

