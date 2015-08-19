#!/usr/bin/python

#!python
#command='Sam.Header.RM.py -i SRR1283824.10.sam'
#sys.argv=command.split()
import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.sam','.head.rm'),'w')
for line in fin:
	pin=line.strip().split('\t')
	if not pin[0][0]=='@':
		print >>fo, '\t'.join(pin)

fin.close()
fo.close()



