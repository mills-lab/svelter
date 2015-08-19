#!/usr/bin/env python

#!python
#command='vcf2simplevcf.py -i vcf_file'
#sys.argv=command.split()
import os
import sys
import getopt

opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv=','Path_Lumpy=','Path_Pindel=','Path_Delly='])
dict_opts=dict(opts)
fin=open(dict_opts['-i'])
fo=open(dict_opts['-i'].replace('.vcf','_simple.vcf'),'w')
fo2=open(dict_opts['-i'].replace('.vcf','_complex.vcf'),'w')
for line in fin:
	pin=line.strip().split()
	if pin[0][0]=='#':
		print >>fo, '\t'.join(pin)
		print >>fo2, '\t'.join(pin)
	else:
		if pin[2].split('_')[-1]=='S':
			print >>fo, '\t'.join(pin)
		elif pin[2].split('_')[-1]=='C':
			print >>fo2, '\t'.join(pin)

fin.close
fo.close()
fo2.close()

