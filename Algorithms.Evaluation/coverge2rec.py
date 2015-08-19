#!/usr/bin/env python

#!python
#command='python coverge2rec.py --ref /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/reference/human_g1k_v37.fasta -p /mnt/EXT/Mills-scratch/datasets/Simulation.Xuefang/Simulate.homo/SVelter/bp_files.Simulate.homo.RD.10.RL.101.sorted.bam/Simulate.homo.RD.10.RL.101.sorted/SPCff3.CluCff3.AlignCff0.0'
#sys.argv=command.split()[1:]
import os
import sys
import getopt
import numpy
opts,args=getopt.getopt(sys.argv[1:],'p:',['ref=','max='])
dict_opts=dict(opts)
def bp_to_let(del_info_unit):
  flag=0
  for i in del_info_unit[0]:
      if i in chromos or not i.isdigit():
          flag+=1
  if not flag==0:
    letter=''.join([chr(i+97) for i in range(len(del_info_unit[0])-2*flag)])
    letters='/'.join([letter,letter])
    return letters
  else:
    return 0

def chromos_name_readin(ref_file):
  ref_index=ref_file+'.fai'
  ref=ref_file
  chromos=[]
  fin=open(ref_index)
  for line in fin:
    pin=line.strip().split()
    chromos.append(pin[0])
  fin.close()
  return chromos

if not dict_opts['-p'][-1]=='/':
  dict_opts['-p']+='/'

chromos=chromos_name_readin(dict_opts['--ref'])
fo=open(dict_opts['-p']+'All.SV.rec','w')
out_hash={}
for k1 in os.listdir(dict_opts['-p']):
  if k1.split('.')[-1]=='coverge':
    fin=open(dict_opts['-p']+k1)
    while True:
      pin1=fin.readline().strip().split()
      if not pin1: break
      pin2=fin.readline().strip().split()
      if not pin1: break
      pin3=fin.readline().strip().split()
      pin4=fin.readline().strip().split()
      pin5=fin.readline().strip().split()
      score=float(pin4[-1])-float(pin3[-1])
      let1=bp_to_let([pin1])
      let2=pin2[0]
      if not let1 in out_hash.keys():
        out_hash[let1]={}
      if not let2 in out_hash[let1].keys():
        out_hash[let1][let2]=[]
      out_hash[let1][let2].append(pin1+[score])
    fin.close()

for k1 in out_hash.keys():
  for k2 in out_hash[k1].keys():
    for k3 in out_hash[k1][k2]:
      print >>fo, ' '.join([str(i) for i in [k1,k2]+k3])

fo.close() 



