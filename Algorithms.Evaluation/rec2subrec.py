#!/usr/bin/env python

#!Python
#sys.argv=command.split()
import os
import sys
import getopt
script_name=sys.argv[0]
opts,args=getopt.getopt(sys.argv[1:],'i:',['num=','ppre=','PathBam=','chr=','SamplePbs=', 'ref=','NullSplitLength=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ex=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
input_file=dict_opts['-i']
number=int(dict_opts['--num'])
fin=open(input_file)
rec=0
rec2=0
pre=input_file.replace('.rec','')
fo=open(pre+'.'+str(rec)+'.rec','w')
for line in fin:
	pin=line.strip().split()
	rec2+=1
	if rec2<number:
		print >>fo, '\t'.join(pin)
	else:
		rec2=0
		rec+=1
		fo.close()
		fo=open(pre+'.'+str(rec)+'.rec','w')
		print >>fo, '\t'.join(pin)

fo.close()


