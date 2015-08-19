#!/usr/bin/env python

#!python
#command='Vali.Cases.size.filter.py -i Validated.cases  --min 100 --max 5000'
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
import numpy as np
import scipy
from scipy import stats

opts,args=getopt.getopt(sys.argv[1:],'i:',['min=','max=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i']+'.min'+dict_opts['--min']+'.max'+dict_opts['--max'],'w')
for line in fin:
	pin=line.strip().split()
	if int(pin[0].split('.')[-2].split('_')[2])-int(pin[0].split('.')[-2].split('_')[1])> int(dict_opts['--min']) and not int(pin[0].split('.')[-2].split('_')[2])-int(pin[0].split('.')[-2].split('_')[1])> int(dict_opts['--max']):
		print >>fo, ' '.join(pin)

fin.close()
fo.close()


