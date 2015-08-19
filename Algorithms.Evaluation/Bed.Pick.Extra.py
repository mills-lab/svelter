#!/usr/bin/env python

#!python
#command='Bed.Pick.Extra.Part.py -a SVelter.Validated.cases.bed -b SVelter.Pindel.cases.bed -o SVelter.Extra.cases.bed'
#sys.argv=command.split()
import os
import sys
import getopt
import re
import pickle
import time
import datetime

opts,args=getopt.getopt(sys.argv[1:],'a:b:o:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-a'])
all_data=[]
for line in fin:
	pin=line.strip().split()
	all_data.append('_'.join(pin))

fin.close()
fin=open(dict_opts['-b'])
all_data2=[]
for line in fin:
	pin=line.strip().split()
	all_data2.append('_'.join(pin))

fin.close()
fo=open(dict_opts['-o'],'w')
for x1 in all_data:
	if not x1 in all_data2:
		print >>fo, '\t'.join(x1.split('_'))

fo.close()

