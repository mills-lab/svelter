#!/usr/bin/env python

#!python
#command='Vali.Cases2Bed.py -i Validated.cases'
#sys.argv=command.split()
#script for validate large dels >3000, where not a lot of pacbio reads goes through entirely
import os
import sys
import getopt
import re
import pickle
import time
import datetime
import random
import numpy
import glob
opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i']+'.bed','w')
for line in fin:
	pin=line.strip().split()
	print >>fo, '\t'.join(pin[0].split('.')[-2].split('_'))

fin.close()
fo.close()







