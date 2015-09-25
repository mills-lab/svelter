#!/usr/bin/env python

#!python
#sys.argv=command.split()
#command='SVelter.python --path_ref path1 --path_in path2 --appdix appdix'
#command='Produce.Pseudo.ROC.stats.py --path_ref sv_rec --path_in alt_sv --appdix .Mappable.min100.max100000.bed'
import os
import sys
import getopt
import numpy
import re
import random
import pickle
import time
import datetime
opts,args=getopt.getopt(sys.argv[1:],'o:',['path_ref=','path_in=','appdix='])
dict_opts=dict(opts)
path_in=dict_opts['--path_in']
path_ref=dict_opts['--path_ref']
if not path_ref[-1]=='/':
        path_ref+='/'

if not path_in[-1]=='/':
        path_in+='/'

def bed_read_in(bedfile):
        out={}
        fin=open(bedfile)
        for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                        out[pin[0]]=[]
                out[pin[0]].append([int(pin[1]),int(pin[2])])
        fin.close()
        return out

def hash_compare(ref_bed,samp_bed):
    for k1 in ref_bed.keys():
        if k1 in ref_bed.keys() and k1 in samp_bed.keys():
            stat_hash[k1]=[]
            missed_bed[k1]=[]
            flag1=0
            for k2 in ref_bed[k1]:
                flag2=0
                for k3 in samp_bed[k1]:
                    if k3[1]<k2[0]:
                        continue
                    elif k3[0]>k2[1]:
                        break
                    else:
                        if not float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/max([float(k2[1]-k2[0]),float(k3[1]-k3[0])])<0.5:
                            flag2+=1
                if flag2>0:
                    flag1+=1
                else:
                    missed_bed[k1].append(k2)
            stat_hash[k1]=[flag1,len(ref_bed[k1]),len(samp_bed[k1])]
        else:
            stat_hash[k1]=[0,len(ref_bed[k1]),0]


def write_stat(total_stat,fout):
        fo=open(fout,'w')
        print >>fo, ' '.join(['sample','sv','chrom','overlap','#Ref.SVs','#Samp.SVs'])
        for k1 in total_stat.keys():
                for k2 in total_stat[k1].keys():
                        for k3 in sorted(total_stat[k1][k2].keys()):
                                k4=total_stat[k1][k2][k3]
                                print >>fo, ' '.join([k1,k2,k3]+[str(i) for i in k4])
        fo.close()

def write_beds(total_missed,fout):
        fo=open(fout,'w')
        for k1 in total_missed.keys():
                for k2 in total_missed[k1].keys():
                        for k3 in sorted(total_missed[k1][k2].keys()):
                                for k4 in total_missed[k1][k2][k3]:
                                        print >>fo,' '.join([k3]+[str(i) for i in k4]+[k1,k2])
        fo.close()

#path_ref='/mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/Del.Bed/'
#path_in='/mnt/EXT/Mills-data/xuefzhao/projects.axiom/1000Genome.Validation/vcf_files/'
appdix=dict_opts['--appdix']
ref_hash={}
for k1 in os.listdir(path_ref):
        if appdix in k1:
                keya1=k1.split('.')[0]
                keya2=k1.replace(appdix,'').split('.')[-1].upper()
                if not keya1 in ref_hash.keys():
                        ref_hash[keya1]={}
                if not keya2 in ref_hash[keya1].keys():
                        ref_hash[keya1][keya2]=[]
                ref_hash[keya1][keya2].append(path_ref+k1)


in_hash={}
for k1 in os.listdir(path_in):
        if appdix in k1:
                keya1=k1.split('.')[0]
                keya2=k1.replace(appdix,'').split('.')[-1].upper()
                if not keya1 in in_hash.keys():
                        in_hash[keya1]={}
                if not keya2 in in_hash[keya1].keys():
                        in_hash[keya1][keya2]=[]
                in_hash[keya1][keya2].append(path_in+k1)

total_stat={}
total_missed={}
total_more={}
for k1 in in_hash.keys():
        total_stat[k1]={}
        total_missed[k1]={}
        total_more[k1]={}
        if k1 in ref_hash.keys():
                for k2 in in_hash[k1].keys():
                        if k2 in in_hash[k1].keys():
                                total_stat[k1][k2]={}
                                total_missed[k1][k2]={}
                                total_more[k1][k2]={}
                                ref_bed=bed_read_in(ref_hash[k1][k2][0])
                                samp_bed=bed_read_in(in_hash[k1][k2][0])
                                stat_hash={}
                                missed_bed={}
                                hash_compare(samp_bed,ref_bed)
                                total_more[k1][k2]=missed_bed
                                stat_hash={}
                                missed_bed={}
                                hash_compare(ref_bed,samp_bed)
                                total_stat[k1][k2]=stat_hash
                                total_missed[k1][k2]=missed_bed

fout=path_in+'Pseudo.ROC'+appdix.replace('.bed','')+'.Stats'
write_stat(total_stat,fout)
fout=path_in+'Missed.SVs'+appdix.replace('.bed','')+'.bed'
write_beds(total_missed,fout)
fout=path_in+'Extra.SVs'+appdix.replace('.bed','')+'.bed'
write_beds(total_more,fout)


