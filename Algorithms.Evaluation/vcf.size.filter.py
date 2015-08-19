#!/usr/bin/env python

#!python
#command='python vcf.size.filter.py -i Pindel.het.RD.20.vcf  --size 10'
#sys.argv=command.split()[1:]
#vcf size filter
import os
import sys
import getopt
import numpy
import re
import random
import pickle
opts,args=getopt.getopt(sys.argv[1:],'i:',['size='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.vcf','.LargerThan'+dict_opts['--size']+'.vcf'),'w')
for line in fin:
	pin=line.strip().split()
	if not pin[0][0]=='#':
		if pin[6]=='PASS':
			start=int(pin[1])
			info=pin[7]
			flag=0
			for x in info.split(';'):
				if 'END' in x:
					end=int(x.split('=')[1])
					flag=1
			if flag==1:
				length=end-start
				if length>int(dict_opts['--size']):
					print >>fo, '\t'.join(pin)
	else:
					print >>fo, '\t'.join(pin)

fin.close()
fo.close()


