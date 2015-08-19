#!/usr/bin/python

#!python
#command='python Rec.TRA.Size.Filter.py -i filein --min 100 --max 1000000'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import numpy
opts,args=getopt.getopt(sys.argv[1:],'i:',['min=','max='])
dict_opts=dict(opts)
min=int(dict_opts['--min'])
max=int(dict_opts['--max'])
filein=dict_opts['-i']
fin=open(filein)
fo=open(filein.replace('.rec','.min'+str(min)+'.max'+str(max)+'.rec'),'w')
for line in fin:
  pin=line.strip().split()
  size=[int(pin[2])-int(pin[1]),int(pin[3])-int(pin[2])]
  if numpy.min(size)>min and numpy.min(size)<max:
    print >>fo, ' '.join(pin)

fin.close()
fo.close()

