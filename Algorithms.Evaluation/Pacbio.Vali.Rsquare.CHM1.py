#!/usr/bin/env python

#!python
#command='Pacbio.Vali.Rsquare.CHM1.py --bam /mnt/EXT/Mills-scratch2/Xuefang/CHM1/Pacbio/Merged.3.Bam/SAMN02744161.pacbio.sorted.bam --ref /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/reference/human_g1k_v37.fasta -i /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Ref_SV/inversions.3col.rec --ppre /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Pacbio.Vali/Ref.INV/'
#command='Pacbio.Vali.Rsquare.CHM1.py --bam /mnt/EXT/Mills-scratch2/Xuefang/CHM1/Pacbio/Merged.3.Bam/SAMN02744161.pacbio.sorted.bam --ref /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/reference/human_g1k_v37.fasta -i /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Ref_SV/deletion.3col.3.rec --ppre /mnt/EXT/Mills-scratch2/Xuefang/CHM1/IL500b/Pacbio.Vali/Ref.DEL/'
#command='Pacbio.Vali.Rsquare.CHM1.py --bam /mnt/EXT/Mills-scratch2/Xuefang/NA12878.Pacbio/sorted_final_merged.bam --ref /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa -i /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/NA12878/Pacbio.Vali.rec/CSV.1.rec --ppre /mnt/EXT/Mills-scratch2/Xuefang/NA12878.Pacbio/Pac.Vali.Ref.DEL/'
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

def eu_dis_calcu(fo_ref,rsquare_ref,align_off,delta):
    fo1=open(fo_ref+'.txt')
    temp_data1=[[],[]]
    for line in fo1:
        po1=line.strip().split()
        temp_data1[0].append(int(po1[0])-align_off)
        temp_data1[1].append(int(po1[1]))
    fo1.close()
    temp_data2=[[],[]]
    for x in range(len(temp_data1[0])):
        if numpy.abs(temp_data1[0][x]-temp_data1[1][x])<delta:
            temp_data2[0].append(temp_data1[0][x])
            temp_data2[1].append(temp_data1[1][x])            
    slope, intercept, r_value, p_value, std_err = stats.linregress(temp_data2[0],temp_data2[1])
    if not temp_data1[0]==[]:
        out=[1 for x in range(len(temp_data1[0])) if abs(temp_data1[1][x]-temp_data1[0][x]*slope-intercept)<delta ]
        rsquare_ref.append(float(len(out))/float(len(temp_data1[0])))
        print float(len(out))/float(len(temp_data1[0]))

def remove_files():
    os.system(r'''rm %s'''%(fo_ref+'.txt'))
    os.system(r'''rm %s'''%(fo1_alt+'.txt'))
    os.system(r'''rm %s'''%(fo2_alt+'.txt'))
    os.system(r'''rm %s'''%(fa1))                
    os.system(r'''rm %s'''%(fa2))                
    os.system(r'''rm %s'''%(fa3))                
    os.system(r'''rm %s'''%(fa4)) 
    os.system(r'''rm %s'''%(txt_file.replace('.txt','*.fa')))            
    os.system(r'''rm %s'''%(txt_file.replace('.txt','*.sam')))            

def bps_check(bps):
    flag=0
    for x in range(len(bps)-2):
        if int(bps[x+2])-int(bps[x+1])>10**6:
            flag+=1
    return flag

opts,args=getopt.getopt(sys.argv[1:],'i:',['ref=','bam=','ppre=','sv='])
dict_opts=dict(opts)
out_path=dict_opts['--ppre']
if not out_path[-1]=='/':
    out_path+='/'

if not os.path.isdir(out_path):
    os.system(r'''mkdir %s'''%(out_path))

filein=dict_opts['-i']
fin=open(filein)
sample_name=filein.split('/')[-1].split('.')[0]
case_hash={}
for line in fin:
    pin=line.strip().split()
    if not pin==[]:
        if not pin[0] in case_hash.keys():
            case_hash[pin[0]]={}
        if not pin[1] in case_hash[pin[0]].keys():
            case_hash[pin[0]][pin[1]]=[]
        case_hash[pin[0]][pin[1]].append(pin[2:])

fin.close()
start=0
window_size=10
delta=50
bam_in=dict_opts['--bam']
ref=dict_opts['--ref']
min_length=50
flank_length=100
min_read_compare=20
for k1 in case_hash.keys():
    for k2 in case_hash[k1].keys():
        for k3 in case_hash[k1][k2]:
            #if float(k3[-1])>start:
                case_name='_'.join(k3)
                fo=open(out_path+sample_name+'.'+case_name+'.txt','w')
                print >>fo, ' '.join([k1,k2]+k3)
                fo.close()
                txt_file=out_path+sample_name+'.'+case_name+'.txt'
                temp_bam=txt_file.replace('.txt','.bam')
                temp_sam=txt_file.replace('.txt','.sam')
                fin=open(txt_file)
                pin=fin.readline().strip().split()
                fin.close()
                ref_sv=pin[0]
                alt_sv=pin[1]
                chrom=pin[2]
                bps=pin[2:]
                if bps_check(bps)==0:
                    bl_len_hash={}
                    for x in ref_sv.split('/')[0]:
                        bl_len_hash[x]=int(bps[ord(x)-97+2])-int(bps[ord(x)-97+1])
                    end_point=0
                    for x in alt_sv.split('/')[0]:
                        if not x=='^':
                            end_point+=bl_len_hash[x]
                    if int(bps[-1])-int(bps[1])<100:
                        flank_length=2*(int(bps[-1])-int(bps[1]))
                    else:
                        if int(bps[-1])-int(bps[1])<500:
                            flank_length=int(bps[-1])-int(bps[1])
                        else:
                            flank_length=500
                    fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
                    fout=open(txt_file.replace('.txt','.clip.fa'),'w')
                    fout2=open(txt_file.replace('.txt','.nonclip.fa'),'w')
                    for line in fbam:
                        pbam=line.strip().split()
                        print >>fout2, '>'+pbam[0]
                        print >>fout2, pbam[9]
                        info_bam=cigar_integrate(pbam)
                        for x in range(len(info_bam[0])):
                            print >>fout,'>'+pbam[0]+'.'+str(x)
                            print >>fout,info_bam[0][x]
                    fbam.close()
                    fout.close()
                    fout2.close()
                    os.system(r'''Pacbio.produce.ref.alt.ref.py --fl %d --ref %s --sv %s'''%(flank_length,ref,txt_file))
                    os.system(r'''blasr %s %s -sam -out %s'''%(txt_file.replace('.txt','.nonclip.fa'),txt_file.replace('.txt','.alt1.fa'),txt_file.replace('.txt','.alt1.sam')))
                    sam_in=txt_file.replace('.txt','.alt1.sam')
                    fsam=open(sam_in)
                    rsquare_ref=[]
                    rsquare_alt1=[]
                    rsquare_alt2=[]
                    rec_len=0
                    read_hash={}
                    for line in fsam:
                        psam=line.strip().split()
                        if not psam[0][0]=='@':
                            if int(psam[3])<flank_length/2 and int(psam[3])+cigar2reaadlength(psam[5])>flank_length:
                                if not len(psam[9]) in read_hash.keys():
                                    read_hash[len(psam[9])]=[]
                                read_hash[len(psam[9])].append([psam[9],int(psam[3])])
                            if int(psam[3])<end_point-float(flank_length)*1.5 and int(psam[3])+cigar2reaadlength(psam[5])>end_point:
                                if not len(psam[9]) in read_hash.keys():
                                    read_hash[len(psam[9])]=[]
                                read_hash[len(psam[9])].append([psam[9],int(psam[3])])
                    fsam.close()
                    rec_num=0
                    if not read_hash=={}:
                        for x in sorted(read_hash.keys())[::-1]:
                            if x>1000:
                                delta=x/10
                            for y1 in read_hash[x]:
                                y=y1[0]
                                y2=y1[1]
                                if rec_num>min_read_compare: break
                                rec_num+=1
                                fo=open(out_path+'temp.fa','w')
                                print >>fo, '>temp'
                                for x in chop_x_for_fasta(y):
                                    print >>fo, x
                                fo.close()
                                fa1=out_path+'temp.fa'
                                fa2=txt_file.replace('.txt','.ref.fa')
                                fa3=txt_file.replace('.txt','.alt1.fa')
                                fa4=txt_file.replace('.txt','.alt2.fa')
                                fa5=txt_file.replace('.txt','.alt.fa')
                                fo_ref=txt_file.replace('.txt','.dotplot.ref')
                                fo1_alt=txt_file.replace('.txt','.dotplot.alt1')
                                fo2_alt=txt_file.replace('.txt','.dotplot.alt2')
                                os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,fa2,fo_ref))
                                os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,fa3,fo1_alt))
                                os.system(r'''dotdata.py %d %s %s %s'''%(window_size,fa1,fa4,fo2_alt))
                                if len(y)>rec_len:
                                    rec_len=len(y)
                                    rec_start=y1
                                    os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest'))
                                    os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest'))
                                    os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest'))
                                    os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))                                
                                eu_dis_calcu(fo_ref,rsquare_ref,y2,delta)
                                eu_dis_calcu(fo1_alt,rsquare_alt1,y2,delta)
                                eu_dis_calcu(fo2_alt,rsquare_alt2,y2,delta)
                        os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa2,fo_ref+'longest.png'))
                        os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa3,fo1_alt+'longest.png'))
                        os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa4,fo2_alt+'longest.png'))
                        os.system(r'''mv %s %s'''%(fo_ref+'longest',fo_ref+'longest.start.'+str(rec_start)))
                        os.system(r'''mv %s %s'''%(fo1_alt+'longest',fo1_alt+'longest.start.'+str(rec_start)))
                        os.system(r'''mv %s %s'''%(fo2_alt+'longest',fo2_alt+'longest.start.'+str(rec_start)))
                        fo=open(txt_file.replace('.txt','.rsquare'),'w')
                        for x in range(min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])):
                            print >>fo, ' '.join([str(y) for y in [rsquare_ref[x],rsquare_alt1[x],rsquare_alt2[x]]])
                        fo.close()
                        remove_files()



