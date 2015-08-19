#!/usr/bin/env python

#!python
#command='Produce.Pac.Vali.pbs.py -p ./ --server flux'
#sys.argv=command.split()
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
opts,args=getopt.getopt(sys.argv[1:],'p:',['server=','ref=','sv='])
dict_opts=dict(opts)
import os
bam_NA12878_axiom='/mnt/EXT/Mills-scratch2/Xuefang/CHM1/Pacbio/Merged.3.Bam/SAMN02744161.pacbio.sorted.bam'
#bam_NA12878_flux='/scratch/remills_flux/xuefzhao/NA12878.Pacbio/sorted_final_merged.bam'
ref_axiom='/mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta'
#ref_flux='/nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa'
ppre=dict_opts['-p']
if not ppre[-1]=='/':
        ppre+='/'

def write_axiom_pbs_header(fout,JobToDo):
        fo=open(fout,'w')
        print >>fo, '#!/bin/bash'
        print >>fo, ' '
        print >>fo, '#PBS -N '+JobToDo
        print >>fo, '#PBS -l mem=20gb,walltime=100:0:0,nodes=compute-4-3'
        print >>fo, '#PBS -m a'
        print >>fo, '#PBS -M xuefzhao@umich.edu'
        print >>fo, '#PBS -o '+JobToDo+'.log'
        print >>fo, '#PBS -e '+JobToDo+'.err'
        print >>fo, '#PBS -V'
        print >>fo, '#PBS -d .'
        fo.close()

def write_flux_pbs_header(fout,JobToDo):
        fo=open(fout,'w')#!/bin/bash
        print >>fo, ' '
        print >>fo, '#PBS -N '+JobToDo
        print >>fo, '#PBS -l mem=20gb,walltime=100:0:0,nodes=flux'
        print >>fo, '#PBS -q flux'
        print >>fo, '#PBS -m a'
        print >>fo, '#PBS -M xuefzhao@umich.edu'
        print >>fo, '#PBS -o '+JobToDo+'.log'
        print >>fo, '#PBS -e '+JobToDo+'.err'
        print >>fo, '#PBS -A remills_flux'
        print >>fo, '#PBS -V'
        print >>fo, '#PBS -d .'
        fo.close()

def chromos_readin(ref):
        fin=open(ref+'.fai')
        chromos=[]
        for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
        fin.close()
        return chromos

def path_mkdir(path):
        if not os.path.isdir(path):
                os.system(r'''mkdir %s'''%(path))

if dict_opts['--server']=='axiom':
        bam_NA12878=bam_NA12878_axiom
        ref=ref_axiom
        for k1 in os.listdir(ppre):
                if k1.split('.')[-1]=='rec':
                        fpbs=ppre+k1.replace('.rec','.pbs')
                        JobToDo=k1.replace('.rec','')
                        write_axiom_pbs_header(fpbs,JobToDo)
                        fo=open(fpbs,'a')
                        out_path=ppre+JobToDo
                        print >>fo, ' '.join(['mkdir',out_path])
                        print >>fo, ' '.join(['Pacbio.Vali.CSV.py --ppre',out_path,'--bam',bam_NA12878,'--ref',ref,'-i',ppre+k1,'--Path_Pindel  /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Pindel --Path_Delly /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Delly --Path_Lumpy /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Lumpy']) 
                        fo.close()



