#!/usr/bin/env python

#!Python
#Usage: bed2rec.homo.del.py -i input.bed
#sys.argv=command.split()
import os
import sys
import getopt
script_name=sys.argv[0]
opts,args=getopt.getopt(sys.argv[1:],'i:',['num='])
dict_opts=dict(opts)
filein=dict_opts['-i']
fileout=filein.replace('.bed','.homo.rec')
fin=open(filein)
fo=open(fileout,'w')
for line in fin:
	pin=line.strip().split()
	chrname=pin[0].replace('chr','')
	po=['a/a','/']+pin
	print >>fo, '\t'.join(po)

fin.close()
fo.close()


