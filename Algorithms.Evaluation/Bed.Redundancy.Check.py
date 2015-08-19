#!/usr/bin/python

#!python
#script usd to compare vcf files
#command='python SVelter.py -i /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.homo/SVelter/Simulate.homo.RD.10.RL.101.vcf'
#command='python Bed.Redundancy.Check -i Delly_NA12878_DEL.DEL.bed'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
opts,args=getopt.getopt(sys.argv[1:],'i:',)
dict_opts=dict(opts)
file_in_1=dict_opts['-i']
test=[]
fin=open(file_in_1)
for line in fin:
        pin=line.strip().split()
        if not pin in test:
                test.append(pin)

fin.close()
fo=open(file_in_1,'w')
for x in test:
        print >>fo, '\t'.join(x)

fo.close()

