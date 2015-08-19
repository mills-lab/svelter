#!/usr/bin/python

#!Python
#command='BPvcf.Stat.Analysis.py -i SVelter_comp_RD_20.BPLink.2.stat'
#sys.argv=command.split()
import os
import sys
import getopt
import numpy
import re
import random
opts,args=getopt.getopt(sys.argv[1:],'i:x:',[])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
out_hash={}
for line in fin:
	pin=line.strip().split()
	if not '_'.join(pin[-1].split('_')[-2:]) in out_hash.keys():
		out_hash['_'.join(pin[-1].split('_')[-2:])]=[]
	out_hash['_'.join(pin[-1].split('_')[-2:])].append(pin)

fin.close()
rec=0
for x in out_hash.keys():
	rec+=1
	fo=open(dict_opts['-i'].replace('.stat','.'+str(rec)+'.stat'),'w')
	for y in out_hash[x]:
		print >>fo, ' '.join(y)
	fo.close()

fo=open(dict_opts['-i'].replace('.stat','.integrated.stat'),'w')
for x in out_hash.keys():
	Sens=[]
	Spec=[]
	for y in out_hash[x]:
		if not y[2]=='0':
			Sens.append(float(y[0])/float(y[1]))
			Spec.append(float(y[0])/float(y[2]))
		else:
			print y
			Sens.append(0)
			Spec.append(1)
	print >>fo, ' '.join([str(i) for i in [numpy.mean(Sens),numpy.mean(Spec),x]])

fo.close()



