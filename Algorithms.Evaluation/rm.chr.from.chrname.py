#!/usr/bin/env python

#!python
#command='rm.chr.from.chromname.py -i deletion.3col.bed'
#sys.argv=command.split()
import os
import sys
import getopt
import re
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
list=[]
for line in fin:
	pin=line.strip().split()
	pin[0]=pin[0].replace('chr','')
	list.append(pin)

fin.close()
fo=open(dict_opts['-i'],'w')
for k1 in list:
	print >>fo, '\t'.join(k1)

fo.close()


