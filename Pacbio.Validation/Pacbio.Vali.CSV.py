#!/usr/bin/env python

#!python
#command='Pacbio.Vali.CSV.py --ppre /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/Pacbio.Vali.rec/ --bam /scratch/remills_flux/xuefzhao/NA12878.Pacbio/sorted_final_merged.bam --ref /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa -i /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/SVelter/CSV.PacVali/SVelter_NA12878.all.CSV.shorterthan5kb.PacVali.rec --Path_Delly /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/Delly --Path_Pindel pindel/ --Path_Lumpy /nfs/remills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/Lumpy'
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
import numpy as np
import scipy
from scipy import stats

def path_name_check(Path_Name):
  if not Path_Name[-1]=='/':
    Path_Name+='/'
  return Path_Name

def csv_info_readin(SVelter_in):
  hash1={}
  fin=open(SVelter_in)
  for line in fin:
    pin=line.strip().split()
    if not pin==[]:
      if not pin[2] in hash1.keys():
        hash1[pin[2]]={}
      if not int(pin[3]) in hash1[pin[2]].keys():
        hash1[pin[2]][int(pin[3])]={}
      if not int(pin[-1]) in hash1[pin[2]][int(pin[3])].keys():
        hash1[pin[2]][int(pin[3])][int(pin[-1])]=[]
      if not pin in hash1[pin[2]][int(pin[3])][int(pin[-1])]:
        hash1[pin[2]][int(pin[3])][int(pin[-1])].append(pin)
  fin.close()
  hash2={}
  for k1 in hash1.keys():
    hash2[k1]=[]
    for k2 in sorted(hash1[k1].keys()):
      for k3 in sorted(hash1[k1][k2].keys()):
        for k4 in hash1[k1][k2][k3]:
          hash2[k1].append([k2,k3]+k4)
  return hash2

def alt_SV_produce(k1):
  if 'DEL' in k1.split('.') or 'del' in k1.split('.') or 'Del' in k1.split('.'):
    return ''
  elif 'DUP' in k1.split('.') or 'dup' in k1.split('.') or 'Dup' in k1.split('.'):
    return 'aa'
  elif 'INV' in k1.split('.') or 'inv' in k1.split('.') or 'Inv' in k1.split('.'):
    return 'a^'
  else:
    return 'a'

def Other_Algorithm_result_readin(Path_Delly):
  out_hash={}
  ref='a/a'
  for k1 in os.listdir(Path_Delly):
    if k1.split('.')[-1]=='bed':
      fin=open(Path_Delly+k1)
      alt_allele=alt_SV_produce(k1)
      if not alt_allele=='a':
        for line in fin:
          pin=line.strip().split()
          if pin[-1]=='homo':
            alt='/'.join([alt_allele,alt_allele])
          else:
            alt='/'.join([alt_allele,'a'])
          if not pin[0] in out_hash.keys():
            out_hash[pin[0]]={}
          if not int(pin[1]) in out_hash[pin[0]].keys():
            out_hash[pin[0]][int(pin[1])]={}
          if not int(pin[2]) in out_hash[pin[0]][int(pin[1])].keys():
            out_hash[pin[0]][int(pin[1])][int(pin[2])]=[]
          out_hash[pin[0]][int(pin[1])][int(pin[2])].append([ref,alt])
      fin.close()
  out_hash2={}
  for k1 in out_hash.keys():
    out_hash2[k1]=[]
    for k2 in sorted(out_hash[k1].keys()):
      for k3 in sorted(out_hash[k1][k2].keys()):
        for k4 in out_hash[k1][k2][k3]:
          out_hash2[k1].append([k2,k3]+k4)
  return out_hash2

def compare_hash():
  out={}
  for k1 in SVelter_hash.keys():
    #if k1 in Delly_hash.keys() and k1 in Lumpy_hash.keys():
      out[k1]=[]
      for k2 in SVelter_hash[k1]:
        temp=[[],[],[]]
        for k3 in Delly_hash[k1]:
          if k3[1]<k2[0]: continue
          elif k3[0]>k2[1]: break
          else:
            if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
              if not k3 in temp[0]:
                temp[0].append(k3)
        for k3 in Lumpy_hash[k1]:
          if k3[1]<k2[0]: continue
          elif k3[0]>k2[1]: break
          else:
            if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
              if not k3 in temp[1]:
                temp[1].append(k3)
        for k3 in Pindel_hash[k1]:
          if k3[1]<k2[0]: continue
          elif k3[0]>k2[1]: break
          else:
            if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
              if not k3 in temp[1]:
                temp[2].append(k3)
        #if not temp[0]==[] and not temp[1]==[]:
        out[k1].append([[k2]]+temp)
  return out

def bps_check(bps):
  flag=0
  for x in range(len(bps)-2):
    if int(bps[x+2])-int(bps[x+1])>10**6:
      flag+=1
  return flag

def flank_length_decide(bps):
  if int(bps[-1])-int(bps[1])<100:
    flank_length=2*(int(bps[-1])-int(bps[1]))
  else:
    if int(bps[-1])-int(bps[1])<500:
      flank_length=int(bps[-1])-int(bps[1])
    else:
      flank_length=500
  return flank_length

def cigar_modify(cigar):
  temp=[[]]
  for x in cigar:
    if ord(x)>64 and ord(x)<91:
      temp[-1].append(x)
      temp.append([])
    else:
      if temp[-1]==[]:
        temp[-1].append(x)
      else:
        temp[-1][-1]+=x
  temp.remove([])
  for x in temp:
    if x[1]=='S' and not int(x[0])>min_length:
      x[1]='I'
  return ''.join(''.join(y) for y in temp)

def cigar_integrate(pin):
  pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
  cigar=cigar_modify(pin[5])
  pos_start=int(pin[3])
  cigars=[]
  for m in pcigar.finditer(cigar):
    cigars.append((m.groups()[0],m.groups()[1]))
  read_start=0
  read_end=0
  Seq_pos=[]
  Read_Pos=[]
  Map_Pos=[]
  rec_seq=0
  rec_map=int(pin[3])
  for x in cigars:
    if x[1]=='I':
      rec_seq1=rec_seq
      rec_seq=rec_seq1+int(x[0])
      Seq_pos.append(pin[9][rec_seq1:rec_seq])
      Read_Pos.append([rec_seq1,rec_seq])
      Map_Pos.append([rec_map,rec_map])
    if x[1]=='D' or x[1]=='N':
      Seq_pos.append('')
      Read_Pos.append([rec_seq,rec_seq])
      rec_map1=rec_map
      rec_map=rec_map1+int(x[0])
      Map_Pos.append([rec_map1,rec_map])
    if x[1]=='M':
      rec_seq1=rec_seq
      rec_seq=rec_seq1+int(x[0])
      rec_map1=rec_map
      rec_map=rec_map1+int(x[0])
      Seq_pos.append(pin[9][rec_seq1:rec_seq])
      Read_Pos.append([rec_seq1,rec_seq])
      Map_Pos.append([rec_map1,rec_map])
    if x[1]=='S':
      rec_seq1=rec_seq
      rec_seq=rec_seq1+int(x[0])
      rec_map1=rec_map
      rec_map=rec_map1
      Seq_pos.append(pin[9][rec_seq1:rec_seq])
      Read_Pos.append([rec_seq1,rec_seq])
      Map_Pos.append([rec_map1,rec_map])
  Cigar_group=[]
  for x in cigars:
    if int(x[0])>min_length and not x[1]=='M':
      Cigar_group.append([x])
      Cigar_group.append([])
    else:
      if Cigar_group==[]:
        Cigar_group.append([x])
      else:
        Cigar_group[-1]+=[x]
  Seq2_pos=[]
  Read2_Pos=[]
  Map2_Pos=[]
  rec=-1
  for x in Cigar_group:
    Seq2_pos.append([])
    Read2_Pos.append([])
    Map2_Pos.append([])
    for y in x:
      rec+=1
      Seq2_pos[-1].append(Seq_pos[rec])
      Read2_Pos[-1].append(Read_Pos[rec])
      Map2_Pos[-1].append(Map_Pos[rec])
  chopped_reads=[''.join(x) for x in Seq2_pos if not x==[]]
  chopped_Read_Pos=[[x[0][0],x[-1][1]] for x in Read2_Pos if not x==[]]
  chopped_Map_Pos=[[x[0][0],x[-1][1]] for x in Map2_Pos if not x==[]]
  return [chopped_reads,chopped_Read_Pos,chopped_Map_Pos,[x for x in Cigar_group if not x==[]]]

def modify_read(pin):
  start=int(pin[3])
  real_start=100
  real_end=int(bps[-1])-int(bps[1])+real_start
  pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
  cigar=pin[5]
  map_start=int(pin[3])
  pos_start=0
  cigars=[]
  for m in pcigar.finditer(cigar):
    cigars.append((m.groups()[0],m.groups()[1]))
  pos_rec=[]
  for x in cigars:
    if x[1]=='M':
      map_start+=int(x[0])
      pos_start+=int(x[0])
    if x[1]=='I':
      pos_start+=int(x[0])
    if x[1]=='D':
      map_start+=int(x[0])
    if x[1]=='S':
      pos_start+=int(x[0])
    if map_start>real_start:
      if pos_rec==[]:
        pos_rec.append(pos_start)
    if map_start>real_end:
      if len(pos_rec)==1:
        pos_rec.append(pos_start)
  return pin[9][pos_rec[0]:]

def Draw_dotplot(sam_in):
  fin=open(sam_in)
  test=0
  tread=''
  start=300
  pin_rec=''
  for line in fin:
    pin=line.strip().split()
    if not pin[0][0]=='@':
      if len(pin[9])>test:
        if int(pin[3])<start:
          tread=pin[9]
          test=len(pin[9])
          pin_rec=pin
  fin.close()
  fref2=open(sam_in.replace('.sam','.dotplot.2.fa'),'w')
  print >>fref2, '>'+sam_in.split('/')[-1].replace('.sam','')
  print >>fref2,modify_read(pin_rec)
  fref2.close()
  #os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,bps[0],int(bps[1]),int(bps[-1]),sam_in.replace('.sam','.dotplot.1.fa')))
  dotmatcher = "/home/remills/data/local/bin/dotmatcher"
  os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),sam_in.replace('.sam','.fa'),threshold))
  os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
  os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.svg')))
  os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.pdf')))
  os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),txt_in.replace('.txt','.ref.fa'),threshold))
  os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
  os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.2.svg')))
  os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.2.pdf')))

def cigar2reaadlength(cigar):
  import re
  pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
  cigars=[]
  for m in pcigar.finditer(cigar):
    cigars.append((m.groups()[0],m.groups()[1]))
  MapLen=0
  for n in cigars:
    if n[1]=='M' or n[1]=='D' or n[1]=='N':
      MapLen+=int(n[0])
  return MapLen

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

def rsquare_calcu(fo_ref,rsquare_ref):
  fo1=open(fo_ref+'.txt')
  temp_data1=[[],[]]
  for line in fo1:
    po1=line.strip().split()
    temp_data1[0].append(int(po1[0]))
    temp_data1[1].append(int(po1[1]))
  fo1.close()
  slope, intercept, r_value, p_value, std_err = stats.linregress(temp_data1[0],temp_data1[1])
  rsquare_ref.append(r_value**2)

def eu_dis_calcu_1(fo_ref,rsquare_ref,align_off,delta):
  fo1=open(fo_ref+'.txt')
  temp_data1=[[],[]]
  for line in fo1:
    po1=line.strip().split()
    temp_data1[0].append(int(po1[0])-align_off)
    temp_data1[1].append(int(po1[1]))
  fo1.close()
  #temp_data2=[[],[]]
  #for x in range(len(temp_data1[0])):
  #  if numpy.abs(temp_data1[0][x]-temp_data1[1][x])<delta:
  #    temp_data2[0].append(temp_data1[0][x])
  #    temp_data2[1].append(temp_data1[1][x])      
  #slope, intercept, r_value, p_value, std_err = stats.linregress(temp_data2[0],temp_data2[1])
  slope=1
  intercept=0
  if not temp_data1[0]==[]:
    out=sum([temp_data1[1][x]-temp_data1[0][x]*slope-intercept for x in range(len(temp_data1[0]))])
    rsquare_ref.append(out)
    #out=[1 for x in range(len(temp_data1[0])) if abs(temp_data1[1][x]-temp_data1[0][x]*slope-intercept)<delta ]
    #rsquare_ref.append(float(len(out))/float(len(temp_data1[0])))

def eu_dis_calcu_2(fo_ref,rsquare_ref,align_off,delta):
  fo1=open(fo_ref+'.txt')
  temp_data1=[[],[]]
  for line in fo1:
    po1=line.strip().split()
    temp_data1[0].append(int(po1[0])-align_off)
    temp_data1[1].append(int(po1[1]))
  fo1.close()
  rec1=1
  for x in range(len(temp_data1[0])):
    if abs(temp_data1[1][x]-temp_data1[0][x])<float(temp_data1[0][x])/10.0:
      rec1+=1
  rsquare_ref.append(float(len(temp_data1[0]))/float(rec1))

def dup_decide(structure):
  flag=0
  for x in structure:
    if not x=='^':
      if structure.count(x)>1:
        flag+=1
  return flag

def eu_dis_calcu(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
    structure1=info[2].split('/')[0]
    structure2=info[3].split('/')[0]
    structure3=info[3].split('/')[1]
    if sum([int(dup_decide(structure1)),int(dup_decide(structure2)),int(dup_decide(structure3))])>0:
        eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
        eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
        eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)
    else:
        eu_dis_calcu_2(fo_ref,rsquare_ref,y2,delta)
        eu_dis_calcu_2(fo1_alt,rsquare_alt1,y2,delta)
        eu_dis_calcu_2(fo2_alt,rsquare_alt2,y2,delta)

def remove_files(txt_file):
  os.system(r'''rm %s'''%(txt_file.replace('.txt','*.fa')))      

def bps_check(bps):
  flag=0
  for x in range(len(bps)-2):
    if int(bps[x+2])-int(bps[x+1])>10**6:
      flag+=1
  return flag

def cigar2alignstart(cigar,start,bps,flank_length):
  #eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
  import re
  pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
  cigars=[]
  for m in pcigar.finditer(cigar):
    cigars.append((m.groups()[0],m.groups()[1]))
  read_rec=0
  align_rec=start
  for x in cigars:
    if x[1]=='S':
      read_rec+=int(x[0])
    if x[1]=='M':
      read_rec+=int(x[0])
      align_rec+=int(x[0])
    if x[1]=='D':
      align_rec+=int(x[0])
    if x[1]=='I':
      read_rec+=int(x[0])
    if align_rec>int(bps[1])-flank_length: break
  return [read_rec,int(align_rec)-int(bps[1])+flank_length]

def cigar2alignstart_2(cigar,start,bps):
  #eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
  import re
  pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
  cigars=[]
  for m in pcigar.finditer(cigar):
    cigars.append((m.groups()[0],m.groups()[1]))
  read_rec=0
  align_rec=start
  for x in cigars:
    if x[1]=='S':
      read_rec+=int(x[0])
    if x[1]=='M':
      read_rec+=int(x[0])
      align_rec+=int(x[0])
    if x[1]=='D':
      align_rec+=int(x[0])
    if x[1]=='I':
      read_rec+=int(x[0])
    if align_rec>int(bps[1]): break
  return [read_rec,int(align_rec)-int(bps[1])]

def chop_pacbio_read(bps,flank_length,info):
    block_length={}
    for x in range(len(info[4:])-2):
        block_length[chr(97+x)]=int(info[x+6])-int(info[x+5])
    alA_len=numpy.sum([block_length[x] for x in info[3].split('/')[0] if not x=='^'])
    alB_len=numpy.sum([block_length[x] for x in info[3].split('/')[1] if not x=='^'])
    alRef_len=int(info[-1])-int(info[5])
    len_cff=max([alA_len,alB_len,alRef_len])
    if len_cff>10**4:
        len_cff=min([alA_len,alB_len,alRef_len])
    fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
    out=[]
    out2=[]
    #test=[]
    for line in fbam:
        pbam=line.strip().split()
        #test.append(int(pbam[3])-int(bps[1])+flank_length)
        if not pbam[0]=='@': 
          if int(pbam[3])<int(bps[1])-flank_length+1:
            align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
            align_start=align_info[0]
            miss_bp=align_info[1]
            align_pos=int(pbam[3])
            target_read=pbam[9][align_start:]
            if len(target_read)>flank_length+len_cff:
              out.append(target_read[:max([alA_len,alB_len,alRef_len])+2*flank_length])
              out2.append(miss_bp)
              #test.append(pbam[:9])
          elif int(pbam[3])<int(bps[1])+1:
            align_info=cigar2alignstart_2(pbam[5],int(pbam[3]),bps)
            align_start=align_info[0]
            miss_bp=align_info[1]+flank_length
            target_read=pbam[9][align_start:]
            if len(target_read)>len_cff:
              out.append(target_read[:max([alA_len,alB_len,alRef_len])+flank_length])
              out2.append(miss_bp)
              #test.append(pbam[:9])
    fbam.close()
    return [out,out2]

def Capitalize_ref(fa2):
  fin=open(fa2)
  data=[]
  for line in fin:
    pin=line.strip().split()
    data+=pin
  fin.close()
  for x in range(1,len(data)):
    data[x]=data[x].upper()
  fin=open(fa2,'w')
  for x in data:
    print >>fin,x
  fin.close()

def txt_rec_file_writing(out_path,SVelter_Name,info):
  txt_file=out_path+SVelter_Name+'.txt'
  fo=open(txt_file,'w')
  print >>fo, '\t'.join(info)
  fo.close()

def calcu_eu_dis_algorithm(SVelter_Name,info):
  #eg: info=k2[0][0]
    txt_rec_file_writing(out_path,SVelter_Name,info[2:])
    txt_file=out_path+SVelter_Name+'.txt'
    ref_sv=info[2]
    alt_sv=info[3]
    chrom=info[4]
    bps=info[4:]
    if (int(bps[-1])-int(bps[1]))>1000:
        delta=(int(bps[-1])-int(bps[1]))/10
    else:
        delta=50
    rsquare_ref=[]
    rsquare_alt1=[]
    rsquare_alt2=[]
    if bps_check(bps)==0:
        bl_len_hash={}
        for x in ref_sv.split('/')[0]:
          bl_len_hash[x]=int(bps[ord(x)-97+2])-int(bps[ord(x)-97+1])
        end_point=0
        for x in alt_sv.split('/')[0]:
          if not x=='^':
            end_point+=bl_len_hash[x]
        flank_length=flank_length_decide(bps)
        all_reads=chop_pacbio_read(bps,flank_length,info)
        read_hash=all_reads[0]
        miss_hash=all_reads[1]
        rec_len=0
        os.system(r'''Pacbio.produce.ref.alt.ref.py --fl %d --ref %s --sv %s'''%(flank_length,ref,txt_file))
        fa1=out_path+'temp.fa'
        Capitalize_ref(txt_file.replace('.txt','.ref.fa'))
        Capitalize_ref(txt_file.replace('.txt','.alt1.fa'))
        Capitalize_ref(txt_file.replace('.txt','.alt2.fa'))
        miss_rec=-1
        if not read_hash==[]:
            if len(read_hash)>20:
                new_read_hash=random.sample(range(len(read_hash)),20)
                read_hash2=[read_hash[i] for i in new_read_hash]
                miss_hash2=[miss_hash[i] for i in new_read_hash]
                read_hash=read_hash2
                miss_hash=miss_hash2
            for x in read_hash:
                    miss_rec+=1
                    y=x
                    y2=miss_hash[miss_rec]
                    fo=open(out_path+'temp.fa','w')
                    print >>fo, '>temp'
                    for z in chop_x_for_fasta(y):
                        print >>fo, z
                    fo.close()
                    fo_ref=txt_file.replace('.txt','.dotplot.ref')
                    fo1_alt=txt_file.replace('.txt','.dotplot.alt1')
                    fo2_alt=txt_file.replace('.txt','.dotplot.alt2')
                    os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,txt_file.replace('.txt','.ref.fa'),fo_ref))
                    os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,txt_file.replace('.txt','.alt1.fa'),fo1_alt))
                    os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,txt_file.replace('.txt','.alt2.fa'),fo2_alt))
                    eu_dis_calcu(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
                    if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
                          min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
                          rsquare_ref=rsquare_ref[:min_len]
                          rsquare_alt1=rsquare_alt1[:min_len]
                          rsquare_alt2=rsquare_alt2[:min_len]
                    if max([abs(rsquare_alt1[-1]),abs(rsquare_alt2[-1])])-rsquare_ref[-1] <rec_len:
                            rec_len=max([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
                            rec_start=y2
                            os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest'))
                            os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest'))
                            os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest'))
                            os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))                
            if not os.path.isfile(fo_ref+'longest'):
                  os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest'))
                  os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest'))
                  os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest'))
                  os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))                
            os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.ref.fa'),fo_ref+'longest.png'))
            os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.alt1.fa'),fo1_alt+'longest.png'))
            os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.alt2.fa'),fo2_alt+'longest.png'))
        remove_files(txt_file)
    out1=[[1 for i in rsquare_ref],[abs(float(rsquare_alt1[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0],[abs(float(rsquare_alt2[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0]]
    out_al1=[[],[],[]]
    out_al2=[[],[],[]]
    for x in range(len(out1[2])):
        if out1[1][x]>out1[2][x]:
          out_al1[0].append(out1[0][x])
          out_al1[1].append(out1[1][x])
          out_al1[2].append(out1[2][x])
        else:
          out_al2[0].append(out1[0][x])
          out_al2[1].append(out1[1][x])
          out_al2[2].append(out1[2][x])
    return [out1[0],out_al1[2]+out_al2[1]]

def multi_calcu_eu(Delly_Name,multi_info):
  #eg: multi_info=k2[1]
  Delly_Info=[[],[]]
  for k3 in multi_info:
    temp_info=k3+[k1]+[str(i) for i in k3[:2]]
    temp_result=calcu_eu_dis_algorithm(Delly_Name,temp_info)
    for x in range(len(Delly_Info)):
      Delly_Info[x]+=temp_result[x]
  return Delly_Info

def heter_info_Filter(SVelter_Info):
        SVelter_Info2=[[],[]]
        for x in range(len(SVelter_Info[0])):
            if SVelter_Info[0][x]>SVelter_Info[1][x]:
                SVelter_Info2[0].append(SVelter_Info[0][x])
                SVelter_Info2[1].append(SVelter_Info[1][x])
        if not SVelter_Info2==[[],[]]:
            return SVelter_Info2
        else:
            return SVelter_Info

opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv=','Path_Lumpy=','Path_Pindel=','Path_Delly='])
dict_opts=dict(opts)
start=0
window_size=8
bam_in=dict_opts['--bam']
ref=dict_opts['--ref']
min_length=50
flank_length=100
min_read_compare=20

Path_Delly=dict_opts['--Path_Delly']
Path_Delly=path_name_check(Path_Delly)
Path_Lumpy=dict_opts['--Path_Lumpy']
Path_Lumpy=path_name_check(Path_Lumpy)
Path_Pindel=dict_opts['--Path_Pindel']
Path_Pindel=path_name_check(Path_Pindel)
SVelter_in=dict_opts['-i']
out_path=dict_opts['--ppre']
out_path=path_name_check(out_path)
SVelter_hash=csv_info_readin(SVelter_in)
Delly_hash=Other_Algorithm_result_readin(Path_Delly)
Lumpy_hash=Other_Algorithm_result_readin(Path_Lumpy)
Pindel_hash=Other_Algorithm_result_readin(Path_Pindel)
comparison=compare_hash()
case_rec=0
chromos=[]
fref=open(dict_opts['--ref']+'.fai')
for line in fref:
    pref=line.strip().split()
    chromos.append(pref[0])

fref.close()
fo=open(SVelter_in+'.compare.stat','w')
print >>fo, ' '.join(['info','SVelter','Delly','Lumpy','Ref'])
for k1 in comparison.keys():
  if not k1 in ['chrX','chrY']:
    for k2 in comparison[k1]:
        print k2
        flag_run=0
        for k3 in k2[0][0][5:]:
            if k3 in chromos:
                flag_run+=1
        if flag_run==0:
            print k2
            case_rec+=1
            SVelter_Name='_'.join(['SVelter',str(case_rec)])
            Delly_Name='_'.join(['Delly',str(case_rec)])
            Lumpy_Name='_'.join(['Lumpy',str(case_rec)])
            Pindel_Name='_'.join(['Pindel',str(case_rec)])
            SVelter_Info=calcu_eu_dis_algorithm(SVelter_Name,k2[0][0])
            Delly_Info=multi_calcu_eu(Delly_Name,k2[1])
            Lumpy_Info=multi_calcu_eu(Lumpy_Name,k2[2])
            Pindel_Info=multi_calcu_eu(Pindel_Name,k2[3])
            #SVelter_Stat=numpy.mean([numpy.mean(SVelter_Info[1]),numpy.mean(SVelter_Info[2])])
            SVelter_Info2=heter_info_Filter(SVelter_Info)
            SVelter_Stat=numpy.mean(SVelter_Info2[1])
            if not Delly_Info[0]==[]:
                Delly_Info2=heter_info_Filter(Delly_Info)
                #Delly_Stat=numpy.mean([numpy.mean(Delly_Info[1]),numpy.mean(Delly_Info[2])])
                Delly_Stat=numpy.mean(Delly_Info2[1])
            else:
                Delly_Stat=1
            if not Lumpy_Info[0]==[]:
                Lumpy_Info2=heter_info_Filter(Lumpy_Info)
                Lumpy_Stat=numpy.mean(Lumpy_Info2[1])
            else:
              Lumpy_Stat=1
            if not Pindel_Info[0]==[]:
                Pindel_Info2=heter_info_Filter(Pindel_Info)
                Pindel_Stat=numpy.mean(Pindel_Info2[1])
            else:
              Pindel_Stat=1
            print >>fo, ' '.join([str(i) for i in ['_'.join([str(j) for j in [k1]+k2[0][0][:4]])]+[SVelter_Stat,Delly_Stat,Lumpy_Stat,Pindel_Stat]])

fo.close()


