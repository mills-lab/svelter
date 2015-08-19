#!/usr/bin/env python

#!python
#command='Split.rec.py -i file --size 30'
#sys.argv=command.split()
import os
import getopt
import re
import sys
opts,args=getopt.getopt(sys.argv[1:],'i:',['size=','bam=','ppre=','sv=','Path_Lumpy=','Path_Delly='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
test=[[]]
for line in fin:
	pin=line.strip().split()
	if len(test[-1])<int(dict_opts['--size']):
		test[-1].append(pin)
	else:
		test.append([pin])

fin.close()
for x in range(len(test)):
	fo=open(dict_opts['-i'].replace('.rec','.'+str(x)+'.rec'),'w')
	for y in test[x]:
		print >>fo, '\t'.join(y)
	fo.close()

