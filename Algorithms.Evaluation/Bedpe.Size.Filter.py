#!/usr/bin/python

#!python
#command='Bedpe.Size.Filter.py -i filein --min 100 --max 1000000'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',['min=','max='])
dict_opts=dict(opts)
min=int(dict_opts['--min'])
max=int(dict_opts['--max'])
filein=dict_opts['-i']
fin=open(filein)
fo=open(filein.replace('.bed','.min'+str(min)+'.max'+str(max)+'.bed'),'w')
for line in fin:
	pin=line.strip().split()
	if int(pin[3])-int(pin[2])>min and int(pin[4])-int(pin[1])<max:
		print >>fo, ' '.join(pin)

fin.close()
fo.close()


