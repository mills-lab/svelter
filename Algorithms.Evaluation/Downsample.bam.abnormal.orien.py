#!/usr/bin/env python

#!Python
#command='Downsample.bamfile.py -i NA12878_S1.chr11_95366462_95367193.all.bam --minIL 150 --maxIL 500 -o NA12878_S1.chr11_95366462_95367193.abn.sam'
#sys.argv=command.split()
#Suggestion: do not specify '-o' here. output files would be under: --ppre +'Breakpoints/'
import os
import sys
import getopt
import numpy
import re
import random
import time
import datetime
ILCutoff=0.95
RDCutoff=0.95
TBCutoff=0.9999
SplitCutoff=0.999
ABILCutoff=0.99
DRCutoff=0.99
SPLCutoff=0.85
#ToolMappingQ='/mnt/EXT/Mills-data/xuefzhao/software/bigWigSummary'
#FileMappingQ='/mnt/EXT/Mills-data/xuefzhao/projects/wgEncodeCrgMapabilityAlign36mer.bigWig'
opts,args=getopt.getopt(sys.argv[1:],'i:o:',['minIL=','maxIL=','ppre=','PathBam=','sample=','chr=','SamplePbs=', 'ref=','NullGenomeName=','NullSplitLength=','NullILCI=','NullRDCI=','NullTBCI=','NullILCff=','NullSPCff=','NullDRCff=','NullBed=','NullBedLength=','NullSampleLength=','NullSampleNumber=','ExcludeBed=','IncludeBed=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)
fbam=os.popen(r'''samtools view %s'''%(dict_opts['-i']))
def Reads_Direction_Detect(flag):
    #flag is the number on the second position of each read in bam file
	#	if int(pi[1])&2>0: #two both mapped
	#	elif int(pi[1])&4>0: #the read itself mapped, mate not
	#	elif int(pi[1])&8>0: #mate mapped, the read itself not
    flag2=int(flag)
    if int(flag2)&16==0: 
    		direct_1='+'
    elif int(flag2)&16>0:
    		direct_1='-'
    if int(flag2)&32==0:
    		direct_2='+'
    elif int(flag2)&32>0: 
    		direct_2='-'
    if flag2&8==0:
        return([direct_1,direct_2])
    else:
        return([direct_1,'0'])

QCAlign=20
os.system(r'''samtools view -H %s > %s'''%(dict_opts['-i'],dict_opts['-o']))
fo=open(dict_opts['-o'],'a')
while True:
    pbam=fbam.readline().strip().split()
    if not pbam: break
    if int(pbam[4])<int(QCAlign): continue             #fail quality control, skip
    if int(pbam[1])&4>0: continue           #the read was not mapped, skip
    DRtemp=Reads_Direction_Detect(pbam[1])
    absIL=abs(int(pbam[8]))
    signIL=0
    if DRtemp==['+','+']: signIL+=1
    if DRtemp==['-','-']: signIL+=1
    if not signIL==0:
    	print >>fo, '\t'.join(pbam)

fo.close()


