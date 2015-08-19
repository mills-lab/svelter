#!/usr/bin/env python

#!python
#command='Pacbio.produce.ref.alt.py --fl 500 --ref /mnt/EXT/Mills-scratch/reference/hg19/genome.fa --sv  NA12878_S1.chr12_71532786_71533753.txt'
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
opts,args=getopt.getopt(sys.argv[1:],'i:',['fl=','ref=','sv='])
dict_opts=dict(opts)
fin=open(dict_opts['--sv'])
pin=fin.readline().strip().split()
fin.close()
flank_length=int(dict_opts['--fl'])
ref_sv=pin[0]
ref_sv='/'.join(['1'+ref_sv.split('/')[0]+'2','1'+ref_sv.split('/')[0]+'2'])
alt_sv=pin[1]
alt_sv='/'.join(['1'+alt_sv.split('/')[0]+'2','1'+alt_sv.split('/')[1]+'2'])
chrom=pin[2]
bps=pin[2:]
bps=[bps[0]]+[str(int(bps[1])-flank_length)]+bps[1:]+[str(int(bps[-1])+flank_length)]
ref_hash={}
rec=0

for x in ref_sv.split('/')[0]:
        rec+=1
        fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(dict_opts['--ref'],chrom,int(bps[rec]), int(bps[rec+1])))
        fref.readline().strip().split()
        seq=''
        while True:
                pref=fref.readline().strip().split()
                if not pref: break
                seq+=pref[0]
        fref.close()
        ref_hash[x]=seq

def alt_sv_to_list(alt_sv):
        out=[[],[]]
        for x in alt_sv.split('/')[0]:
                if not x=='^':
                        out[0].append(x)
                else:
                        out[0][-1]+=x
        for x in alt_sv.split('/')[1]:
                if not x=='^':
                        out[1].append(x)
                else:
                        out[1][-1]+=x
        return out

def complementary(seq):
        seq2=[]
        for i in seq:
                if i in 'ATGCN':
                        seq2.append('ATGCN'['TACGN'.index(i)])
                elif i in 'atgcn':
                        seq2.append('atgcn'['tacgn'.index(i)])
        return ''.join(seq2)

def reverse(seq):
        return seq[::-1]

def ref_hash_modi(ref_hash):
        out={}
        for x in ref_hash.keys():
                out[x]=ref_hash[x]
                out[x+'^']=reverse(complementary(ref_hash[x]))
        return out

def chop_x_for_fasta(x):
        if len(x)>60:
                lines=len(x)/60
                out=[]
                for k in range(lines):
                        out.append(x[(k*60):((k+1)*60)])
                out.append(x[((k+1)*60):])
                return out
        else:
                return [x]

ref_hash=ref_hash_modi(ref_hash)
alt_sv_list=alt_sv_to_list(alt_sv)
fo1=open(dict_opts['--sv'].replace('.txt','.ref.fa'),'w')
fo2=open(dict_opts['--sv'].replace('.txt','.alt1.fa'),'w')
fo3=open(dict_opts['--sv'].replace('.txt','.alt2.fa'),'w')
fo4=open(dict_opts['--sv'].replace('.txt','.alt.fa'),'w')
ref_seq=''.join([ref_hash[x] for x in ref_sv.split('/')[0]])
alt_1_seq=''.join([ref_hash[x] for x in alt_sv_list[0]])
alt_2_seq=''.join([ref_hash[x] for x in alt_sv_list[1]])
alt_seq=alt_1_seq+alt_2_seq
print >>fo1, '>'+'_'.join(bps)+'_ref'
for x in chop_x_for_fasta(ref_seq):
        print >>fo1, x

print >>fo2, '>'+'_'.join(bps)+'_alt1'
for x in chop_x_for_fasta(alt_1_seq):
        print >>fo2, x

print >>fo3, '>'+'_'.join(bps)+'_alt2'
for x in chop_x_for_fasta(alt_2_seq):
        print >>fo3, x

print >>fo4, '>'+'_'.join(bps)+'_alt'
for x in chop_x_for_fasta(alt_seq):
        print >>fo4, x


fo1.close()
fo2.close()
fo3.close()


