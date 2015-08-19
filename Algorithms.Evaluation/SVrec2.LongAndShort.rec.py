#!/usr/bin/env python

#!python
#command='SVrec2.LongAndShort.rec.py -i --size 5000 SVelter_SAMN02744161_flag1.SV.rec'
#sys.argv=command.split()
import os
import sys
import getopt
import re
import pickle
opts,args=getopt.getopt(sys.argv[1:],'i:',['size=','bam=','ppre=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo1=open(dict_opts['-i'].replace('.rec','.small.SV.rec'),'w')
fo2=open(dict_opts['-i'].replace('.rec','.large.SV.rec'),'w')
for line in fin:
        pin=line.strip().split()
        size=int(pin[-1])-int(pin[3])
        if size<int(dict_opts['--size']):
                print >>fo1, '\t'.join(pin)
        else:
                print >>fo2, '\t'.join(pin)

fo1.close()
fo2.close()


