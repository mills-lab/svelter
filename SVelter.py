#!/usr/bin/env python

#!Python
#Usage:
#SVelter.py [option] [Parametres]
#option:
#For debug use only
#command='SVelter.py --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.bam'
#sys.argv=command.split()
import os
import re
import sys
import glob
import getopt
import numpy
import scipy
import math
from math import sqrt,pi,exp
from scipy.stats import norm
import random
import pickle
import time
import datetime
import itertools
script_name=sys.argv[0]
if len(sys.argv)<2:
    print 'SVelter-0.1          Last Update:2015-08-20'
    print ''
    print 'SVelter.py Index should be run first:'
    print 'SVelter.py Index [parameters]'
    print 'Required Parameters:'
    print '--workdir, writable working directory.'
    print '--reference, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
    print '--exclude, absolute path of bed file indicating regions to be excluded from analysis. If not provided, no mappable regions will be excluded'
    print '--copyneutral,absolute, path of bed file indicating copy neutural regions based on which null statistical models would be built. If not provided, genome would be randomly sampled for null model'
    print '--svelter-path, folder which contains all SVelter scripts'
    print ''
    print 'To run SVelter.py:'
    print 'SVelter.py [options] [parameters]'
    print ' '
    print 'options:'
    print 'NullModel' 
    print 'BPSearch'  
    print 'BPIntegrate' 
    print 'SVPredict'
    print 'SVIntegrate'
    print ' '
    print 'Required Parameters:'
    print '--workdir, writable working directory.'
    print '--sample: input alignment file in bam format'
    print ' '
    print 'Optional Parameters:'
    print '--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
    print '--null-copyneutral-length, minimum length requirement for --copyneutral regions used to build null model (default: 2000)'
    print '--null-copyneutral-perc, percentage of regions from --copyneutral to utilize (default: 0.1)'
    print '--null-random-length, specify the length of random regions if --copyneutral parameter not used (default: 5000)'
    print '--null-random-num, specify the number of random regions if --copyneutral parameter not used (default: 10000)'
    print '--num-iteration, maximum number of iterations per structure will run in SV predicting step'
    print '--qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary'
    print '--qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/'
    print '--qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported (default: 0.0)'
    print '--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)'
    print '--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)'
    print '--qc-structure, minimum quality score of a resolved structure to be considered as PASS and included in the output vcf file'
    print '--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)'
    print '--prefix, output prefix for vcf and svelter files (default: input.vcf, input.svelter)'
    print '--ploidy, limit algorithm to specific zygosity (0:heterozygous only; 1:homozygous only; 2:both; default:2)'
    print ' '
    print 'for more details, please see: README'
else:
    function_name=sys.argv[1]
    if function_name=='Index':
        def File_Path_Modify(dict_opts):
            if not dict_opts.has_key('-p'):
                dict_opts['-p']='./'
            elif not dict_opts['-p'][-1]=='/':
                dict_opts['-p']+='/'
            return dict_opts
        def Read_Block_From_Position(bps,letters,bounce_length,position,end):
            bps2=[int(i) for i in bps]
            relative_bps=[i-min(bps2) for i in bps2]
            if end=='right':
                    if position-bounce_length<relative_bps[0]:
                        return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
                    elif position-bounce_length>relative_bps[-1]:
                        return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
                    else:
                        for j in range(len(letters)):
                            if position-bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
                                return [letters[j]]+[relative_bps[j]]+[relative_bps[j+1]]
            if end=='left':
                    if position+bounce_length<relative_bps[0]:
                        return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
                    elif position+bounce_length>relative_bps[-1]:
                        return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
                    else:
                        for j in range(len(letters)):
                            if position+bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
                                return([letters[j]]+[relative_bps[j]]+[relative_bps[j+1]])
        def Pick_CN2_Region(CN2_File,bamF,Length_Limit):
            fi=os.popen('samtools view -h %s'%(bamF))
            Chromosome=[]
            while True:
                pi=fi.readline().strip().split()
                if not pi[0][0]=='@': break
                if pi[0]=='@SQ':
                    Chromosome.append(pi[1].split(':')[1])
            CN2_Region={}
            for chrom in Chromosome:
                CN2_Region[chrom]=[]
            fcn2=open(CN2_File)
            while True:
                pcn2=fcn2.readline().strip().split()
                if not pcn2: break
                if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
                    CN2_Region[pcn2[0]].append(pcn2)
                elif int(pcn2[2])-int(pcn2[1])<Length_Limit: continue 
            return [Chromosome,CN2_Region]
        def Fasta_To_Sequence(fasta_file):
            ffa=open(fasta_file)
            seq=ffa.readline().strip().split()
            sequence=[]
            while True:
                seq=ffa.readline().strip().split()
                if not seq: break
                sequence.append(seq)
            seq2=sequence[0][0]
            for seq3 in sequence[1:]:
                seq2=''.join([seq2,seq3[0]])
            if 'N' in seq2:
                return 'ERROR!'
            else:
                return seq2
        def GC_Content_Calculate(seq2):
            NumAT=0
            NumGC=0
            for n in seq2:
                if n=='A' or n=='T':
                    NumAT+=1
                elif n=='G' or n=='C':
                    NumGC+=1
            return [NumGC,NumAT]            
        def Region_Coverage_Calculate(sam_file,Number_Of_Windows,Region_Info):
            fsam=open(sam_file)
            coverage={}
            num_of_wind=(int(Region_Info[2])-int(Region_Info[1]))/100
            for i in range(num_of_wind):
                coverage[i]=[int(Region_Info[1])+i*100,int(Region_Info[1])+i*100+99,0]
            coverage[-1]=[0,0,0]
            while True:
                psam=fsam.readline().strip().split()
                if not psam: break
                Read_Region=[int(psam[3]),int(psam[3])+len(psam[9])]
                left_block=(Read_Region[0]-int(Region_Info[1]))/100
                if left_block in coverage.keys():
                    left_length=coverage[left_block][1]-Read_Region[0]+1
                else:
                    left_length=0
                right_block=(Read_Region[1]-int(Region_Info[1]))/100
                if right_block in coverage.keys():
                    right_length=Read_Region[1]-coverage[right_block][0]
                else:
                    right_length=0
                if left_block<0:
                    left_block=-1
                if right_block<0:
                    right_block=-1
                if left_block>max(coverage.keys()):
                    left_block=-1
                if right_block>max(coverage.keys()):
                    right_block=-1
                if left_block==right_block:
                    coverage[left_block][-1]+=Read_Region[1]-Read_Region[0]
                else:
                    coverage[left_block][-1]+=left_length
                    coverage[right_block][-1]+=right_length
                    if right_block-left_block>1:
                        for k in range(left_block+1,right_block):
                            coverage[k][-1]+=100
            del coverage[-1]
            for k1 in coverage.keys():
                coverage[k1][-1]=float(coverage[k1][-1])/100
            return coverage 
        def cigar2reaadlength(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            MapLen=0
            for n in cigars:
                if n[1]=='M' or n[1]=='D' or n[1]=='N':
                    MapLen+=int(n[0])
            return MapLen
        def Reads_Direction_Detect(flag):
            flag2=int(flag)
            if int(flag2)&16==0: 
                    direct_1='+'
            elif int(flag2)&16>0:
                    direct_1='-'
            if int(flag2)&32==0:
                    direct_2='+'
            elif int(flag2)&32>0: 
                    direct_2='-'
            return([direct_1,direct_2])
        def cigar2split(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            MapLen=[]
            maplen=0
            if cigars[0][1]=='S':
                MapLen.append(0)
            for n in cigars[1:]:
                if n[1]=='M' or n[1]=='D' or n[1]=='N':
                    maplen+=int(n[0])
                if n[1]=='S': 
                    MapLen.append(maplen-1)
            return MapLen
        def cigar2splitlength(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            splitlen=[]
            for i in cigars:
                if i[1]=='S':
                    splitlen.append(int(i[0]))
            return splitlen
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def ex_file_read_in():
            if '--exclude' in dict_opts.keys():
                ex_file=dict_opts['--exclude']
            else:
                ex_file=workdir+'Support/Exclude.bed'
            return ex_file
        def Gap_Hash_Ref1_read_in(Gap_Refs):
            Gap_Hash_Ref1={}
            for Gap_Ref1 in Gap_Refs:
                fgap=open(Gap_Ref1)
                for line in fgap:
                    pgap=line.strip().split()
                    if not pgap[0] in Gap_Hash_Ref1.keys():
                        Gap_Hash_Ref1[pgap[0]]=[]
                    Gap_Hash_Ref1[pgap[0]].append(pgap[1:4])
                fgap.close()
            return Gap_Hash_Ref1
        def write_ExcludeBed(ExcludeBed):
            if not os.path.isfile(ExcludeBed):
                fo=open(ExcludeBed,'w')
                for chr_ex in chromos:
                    print >>fo, ' '.join([chr_ex,'0','0'])
                fo.close()
        def final_regions_decide(Gap_Hash_Ref1,hash_Cor,chrom):
            GapHash=Gap_Hash_Ref1[chrom]
            GapHash1={}
            for key1 in GapHash:
                GapHash1[int(key1[0])]=[int(key1[0]),int(key1[1])]
            GapHash=[]
            for key1 in sorted(GapHash1.keys()):
                GapHash.append(GapHash1[key1])
            for key1 in range(len(GapHash)-1):
                if GapHash[key1+1][0]<GapHash[key1][1]:
                    GapHash[key1]=[min(GapHash[key1]+GapHash[key1+1]), max(GapHash[key1]+GapHash[key1+1])]
                    GapHash[key1+1]=GapHash[key1]
            for key1 in GapHash:
                if GapHash.count(key1)>1:
                    for key2 in range(GapHash.count(key1)-1):
                        del GapHash[GapHash.index(key1)]
            flag1=0
            flag2=0
            Cor_Gap={}
            while True:
                if flag1==len(GapHash) or flag2==len(hash_Cor): break
                elif flag2==len(hash_Cor)-1 and not hash_Cor[flag2][1]>GapHash[flag1][0]: break
                else:
                    if not  GapHash[flag1][0]>hash_Cor[flag2][1] and not GapHash[flag1][1]<hash_Cor[flag2][0]:
                        if not flag2 in Cor_Gap.keys():
                            Cor_Gap[flag2]=[]
                        Cor_Gap[flag2].append(flag1)
                        flag1+=1
                    elif GapHash[flag1][0]>hash_Cor[flag2][1]:
                        flag2+=1
                    elif GapHash[flag1][1]<hash_Cor[flag2][0]:
                        flag1+=1
            Cor_Gap2=[]
            for key1 in range(len(hash_Cor)):
                if not key1 in Cor_Gap.keys():
                    Cor_Gap2.append(hash_Cor[key1])
                else:
                    CGapTemp=[hash_Cor[key1][0]]
                    for key2 in Cor_Gap[key1]:
                        CGapTemp+=GapHash[key2]
                    CGapTemp.append(hash_Cor[key1][1])
                    if CGapTemp[-2]>CGapTemp[-1]:
                        del CGapTemp[-1]
                        del CGapTemp[-1]
                    for key3 in range(len(CGapTemp)/2):
                        Cor_Gap2.append([CGapTemp[2*key3],CGapTemp[2*key3+1]])
            return Cor_Gap2
        def chromos_info_readin(ref_index):
            global chromos
            global Chromo_Length
            chromos=[]
            Chromo_Length={}
            fin=open(ref_index)
            for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
                Chromo_Length[pin[0]]=int(pin[1])
            fin.close()
        def Gap_Hash_Initiate(chromos):
            global Gap_Hash
            Gap_Hash={}
            for chr_ex in chromos:
                Gap_Hash[chr_ex]=[]
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','batch=','prefix=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        Code_path='/'.join(sys.argv[0].split('/')[:-1])+'/'
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print ''
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--reference, absolute path of reference genome. eg: .../SVelter/reference/genome.fa'
            print '--exclude, absolute path of bed file indicating regions to be excluded from analysis. If not provided, no mappable regions will be excluded.'
            print '--copyneutral,absolute, path of bed file indicating copy neutural regions based on which null statistical models would be built. If not provided, genome would be randomly sampled for null model.'
            print '--svelter-path, folder which contains all SVelter scripts.'
        else:
            if not '--workdir' in dict_opts.keys():
                print 'working directory not specified'
                print 'all temporal files would be writen under current directory'
                workdir='./'
                #print 'Error: please specify working directory using --workdir'
            else:
                workdir = path_modify(dict_opts['--workdir'])
            ref_file=0
            if not '--reference' in dict_opts.keys():
                print 'Error: please specify refrence genome using --reference'
            else:
                ref_file=dict_opts['--reference']
                ref_path='/'.join(ref_file.split('/')[:-1])+'/'
                ref_index=ref_file+'.fai'
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                    print 'Please index reference genome using samtools'
                else:
                    if not '--svelter-path' in dict_opts.keys():
                        print 'Error: please specify path of SVelter scripts using --svelter-path'                       
                    else:
                        time1=time.time()
                        ref_path=workdir+'reference/'
                        if not ref_path=='/'.join(ref_file.split('/')[:-1])+'/':
                            ref_path=workdir+'reference/'
                            if not os.path.isdir(ref_path):
                                os.system(r'''mkdir %s'''%(ref_path))
                            os.system(r'''ln -s %s %s'''%(dict_opts['--svelter-path']+'/SVelter*.r',ref_path))  
                            os.system(r'''ln -s %s %s'''%(ref_file,ref_path+'genome.fa'))
                            os.system(r'''ln -s %s %s'''%(ref_index,ref_path+'genome.fa.fai'))
                            if '--copyneutral' in dict_opts.keys():
                                os.system(r'''ln -s %s %s'''%(dict_opts['--copyneutral'],ref_path+'CN2.bed'))
                            if '--exclude' in dict_opts.keys():
                                os.system(r'''ln -s %s %s'''%(dict_opts['--exclude'],ref_path+'Exclude.bed'))
                                print 'symbolic link of reference genome built under '+ref_path
                                print 'symbolic link of CN file built under '+ref_path
                                print 'symbolic link of Exclude file built under '+ref_path                    
                                ref_file=ref_path+'genome.fa'
                                ref_index=ref_file+'.fai'
                                ExcludeBed=ref_path+'Exclude.bed'
                                chromos_info_readin(ref_index)
                        write_ExcludeBed(ExcludeBed)
                        Gap_Refs=[ExcludeBed]
                        Gap_Hash_Ref1=Gap_Hash_Ref1_read_in(Gap_Refs)
                        Gap_Hash_Initiate(chromos)
                        for chrom in chromos:
                            fout_Name='.'.join(ref_file.split('.')[:-1])+'.Mappable.'+chrom+'.bed'
                            fout_N2='.'.join(ref_file.split('.')[:-1])+'.GC_Content.'+chrom
                            if not os.path.isfile(fout_Name):
                                fref=os.popen(r'''samtools faidx %s %s:'''%(ref_file,chrom))
                                pref=fref.readline().strip().split()
                                while True:
                                    pref=fref.readline().strip().split()
                                    if not pref:break
                                    Gap_Hash[chrom].append(pref[0])
                                fref.close()
                                fout=open(fout_Name,'w')
                                fout2=open(fout_N2,'w')
                                hash_key=chrom
                                if not Gap_Hash[hash_key]==[]:
                                    hash_to_Seq=''.join(Gap_Hash[hash_key])
                                    hash_Cor=[]
                                    hash_cal=0
                                    if not hash_to_Seq[0] in ['N','n']:
                                        hash_Cor.append([0])
                                    for hts in hash_to_Seq[1:]:
                                        hash_cal+=1
                                        if hts=='N' or hts=='n':
                                            if len(hash_Cor)==0: continue
                                            else:
                                                if len(hash_Cor[-1])==1:
                                                    hash_Cor[-1].append(hash_cal)
                                        else:
                                            if len(hash_Cor)==0:
                                                hash_Cor.append([hash_cal])
                                            elif len(hash_Cor[-1])==2:
                                                hash_Cor.append([hash_cal])
                                    if len(hash_Cor[-1])==1:
                                        hash_Cor[-1].append(hash_cal)
                                    hakey2=hash_key
                                    if chrom in Gap_Hash_Ref1.keys():
                                        Cor_Gap2=final_regions_decide(Gap_Hash_Ref1,hash_Cor,chrom)
                                        for key1 in Cor_Gap2:
                                            if key1[1]-key1[0]>1000:
                                                print >>fout, ' '.join([hakey2, str(key1[0]),str(key1[1])])
                                                print >>fout2, ' '.join([hakey2, str(key1[0]),str(key1[1])])
                                                GC_out=[]
                                                for key2 in range((key1[1]-key1[0])/100):
                                                    GC_region=hash_to_Seq[(key1[0]+key2*100):(key1[0]+(key2+1)*100)]
                                                    GC_out.append(str(float(GC_region.count('g')+GC_region.count('G')+GC_region.count('c')+GC_region.count('C'))/100.0))
                                                GC_region=hash_to_Seq[(key1[0]+(key2+1)*100):key1[1]]
                                                if len(GC_region)>0:
                                                    GC_out.append(str(float(GC_region.count('g')+GC_region.count('G')+GC_region.count('c')+GC_region.count('C'))/float(len(GC_region))))
                                                print >>fout2, ' '.join(GC_out)
                                fout.close()
                                fout2.close()
                            fout=open(fout_Name)
                            total_out={}
                            for line in fout:
                                pout=line.strip().split()
                                if not pout[0] in total_out:
                                    total_out[pout[0]]=[]
                                total_out[pout[0]].append(pout[1:3])
                            fout.close()
                        time2=time.time()
                        print 'Reference Genome Indexed !'
                        print 'Time Consuming:'+str(time2-time1)
    if function_name=='NullModel':
        def Read_Block_From_Position(bps,letters,bounce_length,position,end):
            bps2=[int(i) for i in bps]
            relative_bps=[i-min(bps2) for i in bps2]
            if end=='right':
                    if position-bounce_length<relative_bps[0]:
                        return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
                    elif position-bounce_length>relative_bps[-1]:
                        return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
                    else:
                        for j in range(len(letters)):
                            if position-bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
                                return [letters[j]]+[relative_bps[j]]+[relative_bps[j+1]]
            if end=='left':
                    if position+bounce_length<relative_bps[0]:
                        return(['left']+[int(relative_bps[0])-flank,int(relative_bps[0])])
                    elif position+bounce_length>relative_bps[-1]:
                        return(['right']+[int(relative_bps[-1]),int(relative_bps[-1])+flank])
                    else:
                        for j in range(len(letters)):
                            if position+bounce_length in range(relative_bps[j],relative_bps[j+1]+1):
                                return([letters[j]]+[relative_bps[j]]+[relative_bps[j+1]])
        def Pick_CN2_Region(cn2_file,bamF,Length_Limit):
            fi=os.popen('samtools view -h %s'%(bamF))
            Chromosome=[]
            while True:
                pi=fi.readline().strip().split()
                if not pi[0][0]=='@': break
                if pi[0]=='@SQ':
                    Chromosome.append(pi[1].split(':')[1])
            fi.close()
            CN2_Region={}
            for chrom in Chromosome:
                CN2_Region[chrom]=[]
            fcn2=open(cn2_file)
            while True:
                pcn2=fcn2.readline().strip().split()
                if not pcn2: break
                if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
                    CN2_Region[pcn2[0]].append(pcn2)
                elif int(pcn2[2])-int(pcn2[1])<Length_Limit: continue 
            fcn2.close()
            return [Chromosome,CN2_Region]
        def Fasta_To_Sequence(fasta_file):
            ffa=open(fasta_file)
            ffa.readline().strip().split()
            sequence=[]
            while True:
                seq=ffa.readline().strip().split()
                if not seq: break
                sequence.append(seq)
            ffa.close()
            seq2=sequence[0][0]
            for seq3 in sequence[1:]:
                seq2+=seq3[0]
                if ''.join(['N' for x in range(100)]) in seq2: break
            if ''.join(['N' for x in range(100)]) in seq2:
                return 'ERROR!'
            else:
                return seq2
        def GC_Content_Calculate(seq2):
            NumAT=0
            NumGC=0
            for n in seq2:
                if n in ['A','T','a','t']:
                    NumAT+=1
                elif n in ['G','C','g','c']:
                    NumGC+=1
            total_num=NumAT+NumGC
            NumAT_new=int(float(NumAT)/float(total_num)*100)
            NumGC_new=int(float(NumGC)/float(total_num)*100)
            return [NumGC_new,NumAT_new]            
        def Region_Coverage_Calculate(sam_file,Number_Of_Windows,Region_Info):
            fsam=open(sam_file)
            coverage={}
            num_of_wind=(int(Region_Info[2])-int(Region_Info[1]))/Window_Size
            for i in range(num_of_wind):
                coverage[i]=[int(Region_Info[1])+i*Window_Size,int(Region_Info[1])+i*Window_Size+Window_Size-1,0]
            coverage[-1]=[0,0,0]
            while True:
                psam=fsam.readline().strip().split()
                if not psam: break
                Read_Region=[int(psam[3]),int(psam[3])+len(psam[9])]
                left_block=(Read_Region[0]-int(Region_Info[1]))/Window_Size
                if left_block in coverage.keys():
                    left_length=coverage[left_block][1]-Read_Region[0]+1
                else:
                    left_length=0
                right_block=(Read_Region[1]-int(Region_Info[1]))/Window_Size
                if right_block in coverage.keys():
                    right_length=Read_Region[1]-coverage[right_block][0]
                else:
                    right_length=0
                if left_block<0:
                    left_block=-1
                if right_block<0:
                    right_block=-1
                if left_block>max(coverage.keys()):
                    left_block=-1
                if right_block>max(coverage.keys()):
                    right_block=-1
                if left_block==right_block:
                    coverage[left_block][-1]+=Read_Region[1]-Read_Region[0]
                else:
                    coverage[left_block][-1]+=left_length
                    coverage[right_block][-1]+=right_length
                    if right_block-left_block>1:
                        for k in range(left_block+1,right_block):
                            coverage[k][-1]+=Window_Size
            fsam.close()
            del coverage[-1]
            for k1 in coverage.keys():
                coverage[k1][-1]=float(coverage[k1][-1])/Window_Size
            return coverage 
        def cigar2reaadlength(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            MapLen=0
            for n in cigars:
                if n[1]=='M' or n[1]=='D' or n[1]=='N':
                    MapLen+=int(n[0])
            return MapLen
        def Reads_Direction_Detect(preadpair):
            flag=preadpair[1]
            ILLength=int(preadpair[7])-int(preadpair[3])
            flag2=int(flag)
            if int(flag2)&16==0: 
                    direct_1='+'
            elif int(flag2)&16>0:
                    direct_1='-'
            if int(flag2)&32==0:
                    direct_2='+'
            elif int(flag2)&32>0: 
                    direct_2='-'
            if ILLength>0:
                return([direct_1,direct_2])
            else:
                return([direct_2,direct_1])
        def cigar2split(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            MapLen=[]
            maplen=0
            if cigars[0][1]=='S':
                MapLen.append(0)
            for n in cigars[1:]:
                if n[1]=='M' or n[1]=='D' or n[1]=='N':
                    maplen+=int(n[0])
                if n[1]=='S': 
                    MapLen.append(maplen-1)
            return MapLen
        def cigar2splitlength(cigar):
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            splitlen=[]
            for i in cigars:
                if i[1]=='S':
                    splitlen.append(int(i[0]))
            return splitlen
        def dict_opts_modify(dict_opts):
            global model_comp
            if not '--null-model' in dict_opts.keys():
                model_comp='S'
            else:
                if dict_opts['--null-model'] in ['S','Simple']:
                    model_comp='S'
                else:
                    model_comp='C'
            global ReadLength
            global ReadLength_Flag
            if '--read-length' in dict_opts.keys():
                ReadLength=int(dict_opts['--read-length'])
                ReadLength_Flag=1
            else:
                ReadLength=0
                ReadLength_Flag=0
            global QCAlign
            if '--qc-align' in dict_opts.keys():
                QCAlign=int(dict_opts['--qc-align'])
            else:
                QCAlign=20
            global QCSplit
            if '--qc-split' in dict_opts.keys():
                QCSplit=int(dict_opts['--qc-split'])
            else:
                QCSplit=20
            global NullSplitLen_perc
            if '--split-min-len' in dict_opts.keys():
                NullSplitLen_perc=int(dict_opts['--split-min-len'])
            else:
                NullSplitLen_perc=0.9
            global NullILCIs
            if '--NullILCI' in dict_opts.keys():
                NullILCIs=dict_opts['--NullILCI']
            else:
                NullILCIs=[0.025,0.05,0.95,0.975]
            global NullRDCIs
            if '--NullRDCI' in dict_opts.keys():
                NullRDCIs=dict_opts['--NullRDCI']
            else:
                NullRDCIs=[0.025,0.05,0.95,0.975]
            global NullTBCIs
            if '--NullTBCI' in dict_opts.keys():
                NullTBCIs=dict_opts['--NullTBCI']
            else:
                NullTBCIs=[0.0001,0.0005,0.9999,0.9995]
            global NullILCff
            if '--NullILCff'in dict_opts.keys():
                NullILCff=dict_opts['--NullILCff']
            else:
                NullILCff=0.999
            global NullSPCff
            if '--NullSPCff' in dict_opts.keys():
                NullSPCff=dict_opts['--NullSPCff']
            else:
                NullSPCff=0.999
            global NullDRCff
            if '--NullDRCff' in dict_opts.keys():
                NullDRCff=dict_opts['--NullDRCff']
            else:
                NullDRCff=0.999
            global KeepFile
            if '--keep-temp-files' in dict_opts.keys():
                KeepFile=dict_opts['--keep-temp-files']
            else:
                KeepFile='Yes'
            global KeepFigure
            if '--keep-temp-figs' in dict_opts.keys():
                KeepFigure=dict_opts['--keep-temp-figs']
            else:
                KeepFigure='No'
        def path_modify(path):
            if path=='':
                path='/'
            if not path[-1]=='/':
                path+='/'
            return path
        def IL_Stat_Calculate(InsertLength):
            TotalILNum=0
            for key in InsertLength.keys():
                    TotalILNum+=InsertLength[key]
            return TotalILNum
        def NullILCI_Calculate(InsertLength,TotalILNum):
            NullILCILeft=[]
            NullILCIRight=[]
            SubILNumleft=0
            SubILNumright=0
            NciLeft=0
            NciRight=0
            for keyn in range(len(sorted(InsertLength.keys()))):
                SubILNumleft+=InsertLength[sorted(InsertLength.keys())[keyn]]
                SubILNumright+=InsertLength[sorted(InsertLength.keys())[-(keyn+1)]]
                if NciLeft<len(NullILCIs)/2:
                        if SubILNumleft<NullILCIs[NciLeft]*float(TotalILNum): continue
                        if not SubILNumleft<NullILCIs[NciLeft]*float(TotalILNum): 
                                if len(NullILCILeft)==NciLeft:
                                        NullILCILeft.append(sorted(InsertLength.keys())[keyn])
                                        NciLeft+=1
                if NciRight<(len(NullILCIs)/2):
                        if SubILNumright<NullILCIs[NciRight]*float(TotalILNum): continue
                        if not SubILNumright<NullILCIs[NciRight]*float(TotalILNum):
                                if len(NullILCIRight)==NciRight:
                                        NullILCIRight.append(sorted(InsertLength.keys())[-(keyn+1)])
                                        NciRight+=1
                if NciLeft==len(NullILCIs)/2 and NciRight==len(NullILCIs)/2: break
            NullILCI=NullILCILeft+sorted(NullILCIRight)
            return NullILCI
        def SplitLenPNum_Calculate(SplitLength):
            SplitLengthPath=NullPath
            if not os.path.isdir(SplitLengthPath):
                    os.system(r'''mkdir %s'''%(SplitLengthPath))
            SplitLengthOutput=SplitLengthPath+'/'+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'.SplitLength'
            fslo=open(SplitLengthOutput,'w')
            print >> fslo, ' '.join(['Length_of_Split_Region', 'Time_of_Observation'])
            Total_Split_Reads=0
            for key in sorted(SplitLength.keys()):
                    print >> fslo, ' '.join([str(key),str(SplitLength[key])])
                    Total_Split_Reads+=SplitLength[key]
            fslo.close()
            Sub_Split_Reads=0
            SplitLenPNum=1
            if not NullSplitLen_perc<1:
                SplitLenPNum=int(NullSplitLen_perc)
            elif NullSplitLen_perc<1:
                for key in sorted(SplitLength.keys()):
                    Sub_Split_Reads+=SplitLength[key]
                    if float(Sub_Split_Reads)/float(Total_Split_Reads) > NullSplitLen_perc:
                        SplitLenPNum=key
                        break
            return SplitLenPNum
        def clean_files():
            os.system('''rm  %s'''%(InsertLenNullTemp))
            os.system('''rm  %s'''%(DRNullTemp))
            os.system('''rm  %s'''%(SplitNullTemp))
            os.system('''rm  %s'''%(ILNullTemp))
            os.system('''rm  %s'''%(TBNullTemp))
            os.system('''rm  %s'''%(RDNullTemp))
        def RDNullfigure2_Modify(RDNullfigure2):
            fin=open(RDNullfigure2)
            pin1=fin.readline().strip().split()
            pin2=fin.readline().strip().split()
            fin.close()
            fo=open(RDNullfigure2,'w')
            print >>fo, ' '.join(pin1)
            print >>fo, ' '.join([str(i) for i in [float(j)/Window_Size for j in pin2]])
            fo.close()
        def TBNullDensity_CleanUP(TBNullDensity):
            if 0 in TBNullDensity.keys():
                del TBNullDensity[0]
            TotalILNum=0
            for key in TBNullDensity.keys():
                TotalILNum+=TBNullDensity[key]
            Min_keys=0
            Max_keys=0
            for x in sorted(TBNullDensity.keys()):
                Min_keys+=TBNullDensity[x]
                if float(Min_keys)/float(TotalILNum)<0.1:
                    del TBNullDensity[x]
                else:
                    break
            for x in sorted(TBNullDensity.keys())[::-1]:
                Max_keys+=TBNullDensity[x]
                if float(Max_keys)/float(TotalILNum)<0.1:
                    del TBNullDensity[x]
                else:
                    break
            TotalILNum=0
            for key in TBNullDensity.keys():
                    TotalILNum+=TBNullDensity[key]
            return TotalILNum
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def SamplingPercentage_readin():
            if '--null-copyneutral-perc' in dict_opts.keys():
                SamplingPercentage=float(dict_opts['--null-copyneutral-perc'])
            else:
                SamplingPercentage=0.1
            return SamplingPercentage
        def chromos_read_in(ref_index):
            chromos=[]
            fin=open(ref_index)
            for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
            fin.close()
            return chromos
        def cn2_length_readin():
            if '--null-copyneutral-length' in dict_opts.keys():
                cn2_length=int(dict_opts['--null-copyneutral-length'])
            else:
                cn2_length=2000
            return cn2_length
        def cn2_region_write(cn2_file):
            if not '--null-random-length' in dict_opts.keys():
                dict_opts['--null-random-length']=5000
            else:
                dict_opts['--null-random-length']=int(dict_opts['--null-random-length'])
            if not '--null-random-num' in dict_opts.keys():
                dict_opts['--null-random-num']=10000
            else:
                dict_opts['--null-random-num']=int(dict_opts['--null-random-num'])
            cn2_length=dict_opts['--null-random-length']-100
            whole_genome={}
            fref=open(ref_index)
            for line in fref:
                pref=line.strip().split()
                whole_genome[pref[0]]=[int(pref[1])]
            fref.close()
            len_genome=0
            for i in whole_genome.keys():
                len_genome+=whole_genome[i][0]
            fo=open(cn2_file,'w')
            for i in whole_genome.keys():
                num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--null-random-num'])
                reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--null-random-length']) for j in range(num_i)]
                for j in sorted(reg_i):
                    print >>fo, ' '.join([i,str(j),str(j+dict_opts['--null-random-length']-1)])
            fo.close()
            SamplingPercentage=1
            return [cn2_length,SamplingPercentage]
        def genome_name_readin():
            if not '--NullGenomeName' in dict_opts.keys():
                genome_name='genome'
            else:
                genome_name=dict_opts['--NullGenomeName']
            return genome_name
        def NullPath_SetUp(out_path):
            NullPath=out_path+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/'
            if not os.path.isdir(NullPath):
                os.system(r'''mkdir %s'''%(NullPath))
            return NullPath
        def PathBP_SetUp(out_path):
            path_BP=out_path+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
            if not os.path.isdir(path_BP):
                os.system(r'''mkdir %s'''%(path_BP))
            return path_BP
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','batch=','prefix=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        dict_opts_modify(dict_opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print ''
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--sample, input alignment file in bam format'
            print ''
            print 'Optional Parameters:'
            print '--chromosome, name of chromosome to run. should match chromosome name in bam file'
            print '--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
            print '--null-copyneutral-length, minimum length requirement for --copyneutral regions used to build null model (default: 2000)'
            print '--null-copyneutral-perc, percentage of regions from --copyneutral to utilize (default: 0.1)'
            print '--null-random-length, specify the length of random regions if --copyneutral parameter not used (default: 5000)'
            print '--null-random-num, specify the number of random regions if --copyneutral parameter not used (default: 10000)'
            print '--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)'
            print '--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)'
            print '--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)'
        else:
            if not '--workdir' in dict_opts.keys():
                print 'working directory not specified'
                print 'all temporal files would be writen under current directory'
                workdir='./'
            else:
                workdir=path_modify(dict_opts['--workdir'])
            if not '--sample' in dict_opts.keys():
                print 'Error: please specify either input file using --sample'
            else:
                bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                bam_files=[dict_opts['--sample']]
                ref_path=workdir+'reference/'
                ref_file=ref_path+'genome.fa'
                ref_index=ref_file+'.fai'
                cn2_file=ref_path+'CN2.bed'
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                else:
                    SamplingPercentage=SamplingPercentage_readin()
                    chromos=chromos_read_in(ref_index)
                    if os.path.isfile(cn2_file):
                        cn2_length=cn2_length_readin()
                    else:
                        cn2_more_info=cn2_region_write(cn2_file)
                        cn2_length=cn2_more_info[0]
                        SamplingPercentage=cn2_more_info[1]
                    if not '--chromosome' in dict_opts.keys():
                        chr_flag=0
                    elif '--chromosome' in dict_opts.keys():
                        chr_flag=1
                        chrom_single=dict_opts['--chromosome']
                        if not chrom_single in chromos:
                            print 'Error: please make sure the chromosome defined by --chromosome is correct based on the reference genome'
                            chromos=[]
                        else:
                            chromos=[chrom_single]
                    if not chromos==[]:
                        genome_name=genome_name_readin()
                        out_path=workdir
                        NullPath=NullPath_SetUp(out_path)
                        path_BP=PathBP_SetUp(out_path)
                        print 'temp files produced under: '+workdir
                        Script_Path=workdir+'reference/'
                        for bamF in bam_files:
                            time1=time.time()
                            if ReadLength_Flag==0:
                                ReadLengthHash={}
                            outputfile=NullPath+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'.null'
                            fo=open(outputfile,'w')
                            print >>fo, ' '.join(['position','GCContent','ReadDepth','SplitReads','AbnormalDirection','ThroughBP'])
                            fo.close()
                            SplitLength={}
                            InsertLength={}
                            fcn2=open(cn2_file)
                            chr_cn2=[]
                            while True:
                                pcn2=fcn2.readline().strip().split()
                                if not pcn2: break
                                if not len(pcn2)==3: break
                                if pcn2[0] in chromos:
                                    if not pcn2[0] in chr_cn2:
                                        chr_cn2.append(pcn2[0])
                                    if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<10**6:
                                        if not random.choice(range(100))>SamplingPercentage*100:
                                            freadpairs=os.popen('''samtools view %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1])+100,int(pcn2[2])-100))
                                            while True:
                                                preadpair=freadpairs.readline().strip().split()
                                                if not preadpair: break
                                                if not int(preadpair[4])>QCAlign: continue
                                                if preadpair[8]=='0': continue
                                                if not abs(int(preadpair[8])) in InsertLength.keys():
                                                    InsertLength[abs(int(preadpair[8]))]=1
                                                else: 
                                                    InsertLength[abs(int(preadpair[8]))]+=1
                                                if 'S' in preadpair[5]:
                                                    SplitLen=cigar2splitlength(preadpair[5])
                                                    for s in SplitLen:
                                                        if not s in SplitLength.keys():
                                                            SplitLength[s]=1
                                                        else:
                                                            SplitLength[s]+=1
                                            freadpairs.close()
                            fcn2.close()
                            if not chr_cn2==[]:
                                SplitLenPNum=SplitLenPNum_Calculate(SplitLength)
                                TotalILNum=IL_Stat_Calculate(InsertLength)
                                NullILCI=NullILCI_Calculate(InsertLength,TotalILNum)
                                Window_Size=int(float(NullILCI[0])/3)
                                cn2_length=max([cn2_length,NullILCI[2]])
                                cn2_max_len=max(cn2_length*100,10**6)
                                ILNullDensity={}
                                RDNullDensity={}
                                DRNullDensity={}
                                TBNullDensity={}
                                SplitNullDensity={}
                                GC_Content={}
                                fcn2=open(cn2_file)
                                while True:
                                    pcn2=fcn2.readline().strip().split()
                                    if not pcn2: break
                                    if not len(pcn2)==3: break
                                    if pcn2[0] in chromos:
                                        if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<cn2_max_len:
                                            if random.choice(range(100))<SamplingPercentage*100:
                                                pos=[int(pcn2[1])+i*Window_Size for i in range(int(float(int(pcn2[2])-int(pcn2[1]))/Window_Size))]
                                                pos0=[pcn2[0]+'_'+str(i*Window_Size) for i in pos]
                                                RDNull=[0 for i in pos]
                                                SplitNull=[0 for i in pos]
                                                DRNull=[0 for i in pos]
                                                TBNull=[0 for i in pos]
                                                ILNull=[0 for i in pos]
                                                readInf=[]
                                                freadpairs=os.popen('''samtools view %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1]),int(pcn2[2])))
                                                while True:
                                                        preadpair=freadpairs.readline().strip().split()
                                                        if not preadpair: break
                                                        if not int(preadpair[4])>QCAlign: continue
                                                        if not int(preadpair[3])>pos[0]: continue
                                                        block_num=max(int(preadpair[3])-pos[0],0)/Window_Size
                                                        if not block_num<len(pos): continue
                                                        if ReadLength_Flag==0:
                                                            if not preadpair[9]=='*':
                                                                if not len(preadpair[9]) in ReadLengthHash.keys():
                                                                    ReadLengthHash[len(preadpair[9])]=1
                                                                else:
                                                                    ReadLengthHash[len(preadpair[9])]+=1
                                                        RDNull[block_num]+=cigar2reaadlength(preadpair[5])
                                                        if not int(preadpair[8])<NullILCI[0] and not int(preadpair[8])>NullILCI[-1]:
                                                                ILNull[block_num]+=1
                                                        if not preadpair[5].find('S')==-1:
                                                                splitpos=[i+int(preadpair[3]) for i in cigar2split(preadpair[5])]
                                                                splitlent=cigar2splitlength(preadpair[5])
                                                                for j in range(len(splitpos)):
                                                                        if not splitlent[j]<SplitLenPNum and splitpos[j] in pos:
                                                                                SplitNull[block_num]+=1
                                                        if preadpair[6]=='=':
                                                                if not Reads_Direction_Detect(preadpair)==['+', '-']:
                                                                        abdrpos=min([int(preadpair[3]),int(preadpair[7])])-int(pcn2[1])-100
                                                                        if abdrpos>-1 and abdrpos<(int(pcn2[2])-int(pcn2[1])-200):
                                                                                DRNull[abdrpos/Window_Size]+=1
                                                                if not preadpair[0] in readInf:
                                                                        for j in range(max((int(preadpair[3])-pos[0])/Window_Size,0), min((int(preadpair[3])+int(preadpair[8])-pos[0])/Window_Size,len(pos))):
                                                                                TBNull[j]+=1
                                                                        readInf.append(preadpair[0])
                                                freadpairs.close()
                                                if not sum(RDNull)==0:
                                                    fref=os.popen('''samtools faidx %s %s:%d-%d'''%(ref_file,pcn2[0],int(pcn2[1]),int(pcn2[1])+(int(pcn2[2])-int(pcn2[1]))/Window_Size*Window_Size))
                                                    tref=fref.readline().strip().split()
                                                    REFSEQUENCE=fref.readline().strip().split()
                                                    while True:
                                                            pref=fref.readline().strip().split()
                                                            if not pref: break
                                                            REFSEQUENCE=[''.join(REFSEQUENCE+pref)]
                                                    fref.close()
                                                    GCNull=[int(float(REFSEQUENCE[0][(Window_Size*i):(Window_Size*i+Window_Size)].count('G')+REFSEQUENCE[0][(Window_Size*i):(Window_Size*i+Window_Size)].count('C')+REFSEQUENCE[0][(Window_Size*i):(Window_Size*i+Window_Size)].count('g')+REFSEQUENCE[0][(Window_Size*i):(Window_Size*i+Window_Size)].count('c'))/float(Window_Size)*100) for i in range(len(REFSEQUENCE[0])/Window_Size)]
                                                    fo=open(outputfile,'a')
                                                    for k in range(len(pos)):
                                                            if not RDNull[k]==0:
                                                                    print >>fo, ' '.join([str(pos0[k]),str(GCNull[k]),str(RDNull[k]),str(SplitNull[k]),str(DRNull[k]),str(TBNull[k])])
                                                    for i in range(len(RDNull)):
                                                            if not GCNull[i] in GC_Content.keys():
                                                                    GC_Content[GCNull[i]]=[RDNull[i]]
                                                            if GCNull[i] in GC_Content.keys():
                                                                    GC_Content[GCNull[i]].append(RDNull[i])
                                                    fo.close()
                                                    for k in range(len(pos)):
                                                            if not ILNull[k] in ILNullDensity.keys():
                                                                    ILNullDensity[ILNull[k]]=1
                                                            elif ILNull[k] in ILNullDensity.keys():
                                                                    ILNullDensity[ILNull[k]]+=1
                                                            if not RDNull[k] in RDNullDensity.keys():
                                                                    RDNullDensity[RDNull[k]]=1
                                                            elif RDNull[k] in RDNullDensity.keys():
                                                                    RDNullDensity[RDNull[k]]+=1
                                                            if not SplitNull[k] in SplitNullDensity.keys():
                                                                    SplitNullDensity[SplitNull[k]]=1
                                                            elif SplitNull[k] in SplitNullDensity.keys():
                                                                    SplitNullDensity[SplitNull[k]]+=1
                                                            if not DRNull[k] in DRNullDensity.keys():
                                                                    DRNullDensity[DRNull[k]]=1
                                                            elif DRNull[k] in DRNullDensity.keys():
                                                                    DRNullDensity[DRNull[k]]+=1
                                                    for k in range(len(pos))[1:]:
                                                            if not TBNull[k] in TBNullDensity.keys():
                                                                    TBNullDensity[TBNull[k]]=1
                                                            elif TBNull[k] in TBNullDensity.keys():
                                                                    TBNullDensity[TBNull[k]]+=1
                                fcn2.close()
                                if 0 in RDNullDensity.keys():
                                        del RDNullDensity[0]
                                if 0 in TBNullDensity.keys():
                                        del TBNullDensity[0]
                                OverallRDDenominator=0
                                OverallRDNumeritor=0
                                if not RDNullDensity=={}:
                                    for key in RDNullDensity.keys():
                                        if not key==0:
                                            OverallRDNumeritor+=key*RDNullDensity[key]
                                            OverallRDDenominator+=RDNullDensity[key]
                                    OverallRDNullDensity=float(OverallRDNumeritor)/float(OverallRDDenominator)
                                    fbRD=open(outputfile)
                                    pbRD=fbRD.readline().strip().split()
                                    RD_Af_Adj={}
                                    for key in GC_Content.keys():
                                            GC_Content[key]=numpy.mean(GC_Content[key])
                                    while True:
                                            pbRD=fbRD.readline().strip().split()
                                            if not pbRD: break
                                            if not len(pbRD)==6 : break
                                            RDAfAdj=int(pbRD[2])*OverallRDNullDensity/GC_Content[int(pbRD[1])]
                                            if not int(RDAfAdj) in RD_Af_Adj.keys():
                                                    RD_Af_Adj[int(RDAfAdj)]=1
                                            elif int(RDAfAdj) in RD_Af_Adj.keys():
                                                    RD_Af_Adj[int(RDAfAdj)]+=1
                                    fbRD.close()
                                    RDMedian=numpy.median(RD_Af_Adj.keys())
                                    for key in RD_Af_Adj.keys():
                                            if key > RDMedian*10 or key==0:
                                                    del RD_Af_Adj[key]
                                    TotalRDNum=0
                                    for key in RD_Af_Adj.keys():
                                            TotalRDNum+=RD_Af_Adj[key]
                                    NullRDCILeft=[]
                                    NullRDCIRight=[]
                                    SubRDNumleft=0
                                    SubRDNumright=0
                                    NciLeft=0
                                    NciRight=0
                                    for keyn in range(len(sorted(RD_Af_Adj.keys()))):
                                            SubRDNumleft+=RD_Af_Adj[sorted(RD_Af_Adj.keys())[keyn]]
                                            SubRDNumright+=RD_Af_Adj[sorted(RD_Af_Adj.keys())[-(keyn+1)]]
                                            if NciLeft<len(NullRDCIs)/2:
                                                    if SubRDNumleft<NullRDCIs[NciLeft]*float(TotalRDNum): continue
                                                    if not SubRDNumleft<NullRDCIs[NciLeft]*float(TotalRDNum): 
                                                            if len(NullRDCILeft)==NciLeft:
                                                                    NullRDCILeft.append(sorted(RD_Af_Adj.keys())[keyn])
                                                                    NciLeft+=1
                                            if NciRight<(len(NullRDCIs)/2):
                                                    if SubRDNumright<NullRDCIs[NciRight]*float(TotalRDNum): continue
                                                    if not SubRDNumright<NullRDCIs[NciRight]*float(TotalRDNum):
                                                            if len(NullRDCIRight)==NciRight:
                                                                    NullRDCIRight.append(sorted(RD_Af_Adj.keys())[-(keyn+1)])
                                                                    NciRight+=1
                                            if NciLeft==len(NullRDCIs)/2 and NciRight==len(NullRDCIs)/2: break
                                    NullRDCI=NullRDCILeft+sorted(NullRDCIRight)
                                    TotalTBNum=TBNullDensity_CleanUP(TBNullDensity)
                                    NullTBCILeft=[]
                                    NullTBCIRight=[]
                                    SubTBNumleft=0
                                    SubTBNumright=0
                                    NciLeft=0
                                    NciRight=0
                                    for keyn in range(len(sorted(TBNullDensity.keys()))):
                                            SubTBNumleft+=TBNullDensity[sorted(TBNullDensity.keys())[keyn]]
                                            SubTBNumright+=TBNullDensity[sorted(TBNullDensity.keys())[-(keyn+1)]]
                                            if NciLeft<len(NullTBCIs)/2:
                                                    if SubTBNumleft<NullTBCIs[NciLeft]*float(TotalTBNum): continue
                                                    if not SubTBNumleft<NullTBCIs[NciLeft]*float(TotalTBNum): 
                                                            if len(NullTBCILeft)==NciLeft:
                                                                    NullTBCILeft.append(sorted(TBNullDensity.keys())[keyn])
                                                                    NciLeft+=1
                                            if NciRight<(len(NullTBCIs)/2):
                                                    if SubTBNumright<NullTBCIs[NciRight]*float(TotalTBNum): continue
                                                    if not SubTBNumright<NullTBCIs[NciRight]*float(TotalTBNum):
                                                            if len(NullTBCIRight)==NciRight:
                                                                    NullTBCIRight.append(sorted(TBNullDensity.keys())[-(keyn+1)])
                                                                    NciRight+=1
                                            if NciLeft==len(NullTBCIs)/2 and NciRight==len(NullTBCIs)/2: break
                                    NullTBCI=NullTBCILeft+sorted(NullTBCIRight)
                                    TotalILNum=0
                                    for key in ILNullDensity.keys():
                                            TotalILNum+=ILNullDensity[key]
                                    ILITotal=0
                                    for key in sorted(ILNullDensity.keys()):
                                            ILITotal+=ILNullDensity[key]
                                            if float(ILITotal)/float(TotalILNum)>NullILCff: break
                                    if sorted(ILNullDensity.keys()).index(key)>4 or len(ILNullDensity.keys())<4:
                                        ILIPoint=sorted(ILNullDensity.keys())[:sorted(ILNullDensity.keys()).index(key)+1]
                                        ILIPoint2=[float(ILNullDensity[i])/float(TotalILNum) for i in ILIPoint]
                                        ILIPoint3=[]
                                        for k in range(len(ILIPoint)):
                                            ILIPoint3.append(sum(ILIPoint2[:(k+1)]))
                                    else:
                                        ILIPoint=sorted(ILNullDensity.keys())[:4]
                                        ILIPoint2=[float(ILNullDensity[i])/float(TotalILNum) for i in ILIPoint]
                                        ILIPoint3=[]
                                        for k in range(len(ILIPoint)):
                                            ILIPoint3.append(sum(ILIPoint2[:(k+1)]))
                                    TotalSplitNum=0
                                    for key in SplitNullDensity.keys():
                                            TotalSplitNum+=SplitNullDensity[key]
                                    SplitITotal=0
                                    for key in sorted(SplitNullDensity.keys()):
                                            SplitITotal+=SplitNullDensity[key]
                                            if float(SplitITotal)/float(TotalSplitNum)>NullSPCff: break
                                    if sorted(SplitNullDensity.keys()).index(key)>4 or len(SplitNullDensity.keys())<4:
                                        SplitIPoint=sorted(SplitNullDensity.keys())[:sorted(SplitNullDensity.keys()).index(key)+1]
                                        SplitIPoint2=[float(SplitNullDensity[i])/float(TotalSplitNum) for i in SplitIPoint]
                                        SplitIPoint3=[]
                                        for k in range(len(SplitIPoint)):
                                            SplitIPoint3.append(sum(SplitIPoint2[:(k+1)]))
                                    else:
                                        SplitIPoint=sorted(SplitNullDensity.keys())[:4]
                                        SplitIPoint2=[float(SplitNullDensity[i])/float(TotalSplitNum) for i in SplitIPoint]
                                        SplitIPoint3=[]
                                        for k in range(len(SplitIPoint)):
                                            SplitIPoint3.append(sum(SplitIPoint2[:(k+1)]))
                                    TotalDRNum=0
                                    for key in DRNullDensity.keys():
                                            TotalDRNum+=DRNullDensity[key]
                                    DRITotal=0
                                    for key in sorted(DRNullDensity.keys()):
                                            DRITotal+=DRNullDensity[key]
                                            if float(DRITotal)/float(TotalDRNum)>NullDRCff: break
                                    if sorted(DRNullDensity.keys()).index(key)>4 or len(DRNullDensity.keys())<4:
                                        DRIPoint=sorted(DRNullDensity.keys())[:(sorted(DRNullDensity.keys()).index(key)+1)]
                                        DRIPoint2=[float(DRNullDensity[i])/float(TotalDRNum) for i in DRIPoint]
                                        DRIPoint3=[]
                                        for k in DRIPoint:
                                            DRIPoint3.append(sum(DRIPoint2[:(k+1)]))
                                    else:
                                        DRIPoint=sorted(DRNullDensity.keys())[:4]
                                        DRIPoint2=[float(DRNullDensity[i])/float(TotalDRNum) for i in DRIPoint]
                                        DRIPoint3=[]
                                        for k in DRIPoint:
                                            DRIPoint3.append(sum(DRIPoint2[:(k+1)]))
                                    outputfileStat=NullPath+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.Stats'
                                    fos=open(outputfileStat,'w')
                                    print >>fos, 'Insert Lenth CIs'
                                    print >>fos, ' '.join([str(i) for i in NullILCIs])
                                    print >>fos, ' '.join([str(i) for i in NullILCI])
                                    print >>fos, 'Read Depth CIs'
                                    print >>fos, ' '.join([str(i) for i in NullRDCIs])
                                    print >>fos, ' '.join([str(i) for i in NullRDCI])
                                    print >>fos, 'Number of Reads Going Through a Break Points CIs'
                                    print >>fos, ' '.join([str(i) for i in NullTBCIs])
                                    print >>fos, ' '.join([str(i) for i in NullTBCI])
                                    print >>fos, 'Number of Read Pairs with Aberrant Insert Length'
                                    print >>fos, ' '.join([str(i) for i in ILIPoint])
                                    print >>fos, ' '.join([str(i) for i in ILIPoint3])
                                    print >>fos, 'Number of Split Reads'
                                    print >>fos, ' '.join([str(i) for i in SplitIPoint])
                                    print >>fos, ' '.join([str(i) for i in SplitIPoint3])
                                    print >>fos, 'Number of Read Pairs with Aberrant Direction'
                                    print >>fos, ' '.join([str(i) for i in DRIPoint])
                                    print >>fos, ' '.join([str(i) for i in DRIPoint3])
                                    if ReadLength_Flag==0:
                                        ReadLengthTag=0
                                        ReadLengthOut=0
                                        for RLK1 in ReadLengthHash.keys():
                                            if ReadLengthHash[RLK1]>ReadLengthTag:
                                                ReadLengthOut=RLK1
                                                ReadLengthTag=ReadLengthHash[RLK1]
                                        print >>fos, 'Read Length Of Reads'+':'+str(ReadLengthOut)
                                    fos.close()
                                    outputfileIL=NullPath+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.density.null'
                                    foIL=open(outputfileIL,'w')
                                    print >>foIL, ' '.join(['InsertLength','Frequency'])
                                    for l in sorted(InsertLength.keys()):
                                            print >>foIL, ' '.join([str(l),str(InsertLength[l])])
                                    print >>foIL, ' '.join(['ReadDepth','Frequency'])
                                    for r in RD_Af_Adj.keys():
                                            print >>foIL, ' '.join([str(r),str(RD_Af_Adj[r])])
                                    print >>foIL, ' '.join(['ThroughBreakPoint','Frequency'])
                                    for t in TBNullDensity.keys():
                                            print >>foIL, ' '.join([str(t),str(TBNullDensity[t])])
                                    print >>foIL,' '.join(['BinPosition','GC_Content'])
                                    for b in GC_Content.keys():
                                        if not b==0:
                                            print >>foIL, ' '.join([str(b), str(GC_Content[b])])
                                    foIL.close()
                                    RFigureDRSplit=Script_Path+'SVelter1.NullModel.Figure.a.r'	
                                    InsertLenNullTemp=NullPath+'InsertLenNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fIL=open(InsertLenNullTemp,'w')
                                    for dr in ILNullDensity.keys():
                                            print >> fIL, ' '.join([str(dr),str(ILNullDensity[dr])])
                                    fIL.close()
                                    InsertLenNullfigure1='.'.join(InsertLenNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    InsertLenNullfigure2='.'.join(InsertLenNullTemp.split('.')[:-1])+'.2.jpg'
                                    if KeepFigure in ['no','N','No','n']:
                                        InsertLenNullfigure1=InsertLenNullfigure1.replace('.jpg','.na')
                                        InsertLenNullfigure2=InsertLenNullfigure2.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,InsertLenNullTemp,InsertLenNullfigure1,BoxPlotColor,InsertLenNullfigure2))
                                    DRNullTemp=NullPath+'DirectionNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fDR=open(DRNullTemp,'w')
                                    for dr in DRNullDensity.keys():
                                            print >> fDR, ' '.join([str(dr),str(DRNullDensity[dr])])
                                    fDR.close()
                                    DRNullfigure1='.'.join(DRNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    DRNullfigure2='.'.join(DRNullTemp.split('.')[:-1])+'.2.jpg'
                                    if KeepFigure in ['no','N','No','n']:
                                        DRNullfigure1=DRNullfigure1.replace('.jpg','.na')
                                        DRNullfigure2=DRNullfigure2.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,DRNullTemp,DRNullfigure1,BoxPlotColor,DRNullfigure2))
                                    SplitNullTemp=NullPath+'SplitNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fSP=open(SplitNullTemp,'w')
                                    for sp in SplitNullDensity.keys():
                                            print >> fSP, ' '.join([str(sp),str(SplitNullDensity[sp])])
                                    fSP.close()
                                    SplitNullfigure1='.'.join(SplitNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    SplitNullfigure2='.'.join(SplitNullTemp.split('.')[:-1])+'.2.jpg'
                                    if KeepFigure in ['no','N','No','n']:
                                        SplitNullfigure1=SplitNullfigure1.replace('.jpg','.na')
                                        SplitNullfigure2=SplitNullfigure2.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,SplitNullTemp,SplitNullfigure1,BoxPlotColor,SplitNullfigure2))
                                    if model_comp=='C':
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.d.r'	
                                    else:
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.d2.r'    
                                    RDNullTemp=NullPath+'RDNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fRD=open(RDNullTemp,'w')
                                    for rd in RD_Af_Adj.keys():
                                            print >> fRD, ' '.join([str(rd),str(RD_Af_Adj[rd])])
                                    fRD.close()
                                    RDNullfigure1='.'.join(RDNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    RDNullfigure2='.'.join(RDNullTemp.split('.')[:-1])+'.NegativeBinomial'
                                    if KeepFigure in ['no','N','No','n']:
                                        RDNullfigure1=RDNullfigure1.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,RDNullTemp,RDNullfigure1,BoxPlotColor,lineColor,RDNullfigure2))
                                    RDNullfigure2_Modify(RDNullfigure2)
                                    ILNullTemp=NullPath+'ILNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fIL=open(ILNullTemp,'w')
                                    for il in InsertLength.keys():
                                            print >> fIL, ' '.join([str(il),str(InsertLength[il])])
                                    fIL.close()
                                    ILNullfigure1='.'.join(ILNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    ILNullfigure2='.'.join(ILNullTemp.split('.')[:-1])+'.Bimodal'
                                    if KeepFigure in ['no','N','No','n']:
                                        ILNullfigure1=ILNullfigure1.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,ILNullTemp,ILNullfigure1,BoxPlotColor,lineColor,ILNullfigure2))
                                    TBNullTemp=NullPath+'TBNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fTB=open(TBNullTemp,'w')
                                    for tb in TBNullDensity.keys():
                                        print >> fTB, ' '.join([str(tb),str(TBNullDensity[tb])])
                                    fTB.close()
                                    TBNullfigure1='.'.join(TBNullTemp.split('.')[:-1]+['jpg'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    TBNullfigure2='.'.join(TBNullTemp.split('.')[:-1])+'.Bimodal'
                                    if KeepFigure in ['no','N','No','n']:
                                        TBNullfigure1=TBNullfigure1.replace('.jpg','.na')
                                    os.system('''Rscript %s %s %s %s %s %s'''%(RFigureDRSplit2,TBNullTemp,TBNullfigure1,BoxPlotColor,lineColor,TBNullfigure2))
                                    print 'Rscript'
                                    clean_files()
                                    Ref_Seq_File=ref_file
                                    Mini_CN2_Region=int(cn2_length)
                                    Length_Limit=int(cn2_length)
                                    CN2_Region={} #key of hash CN2_Region is the name of each chromosome
                                    for chrom in chromos:
                                            CN2_Region[chrom]={} #key of CN2_Region[chrom] is GC_content
                                            for con in range(101):
                                                    CN2_Region[chrom][con]=[]
                                    fcn2=open(cn2_file)
                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]
                                    while True:
                                            pcn2=fcn2.readline().strip().split()
                                            if not len(pcn2)==3: break
                                            Chromosome=pcn2[0]
                                            if Chromosome in CN2_Region.keys():
                                                if int(pcn2[2])-int(pcn2[1])<Length_Limit: continue
                                                if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(Ref_Seq_File,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
                                                    Seq1=Fasta_To_Sequence(fasta_file)
                                                    if Seq1=='ERROR!':continue
                                                    if not Seq1=='ERROR!':
                                                        sam_file=NullPath+temp_Name+'.sam'
                                                        os.system(r'''samtools view %s %s:%d-%d > %s'''%(bamF,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),sam_file))
                                                        Number_Of_Windows=len(Seq1)/Window_Size
                                                        GC_Content={}
                                                        for i in range(len(Seq1)/Window_Size+1)[1:]:				
                                                                Seq2=Seq1[(i-1)*Window_Size:i*Window_Size]
                                                                GC_Content[i]=GC_Content_Calculate(Seq2)
                                                        coverage=Region_Coverage_Calculate(sam_file,Number_Of_Windows,pcn2)
                                                        for j in GC_Content.keys():
                                                            if j in coverage.keys():
                                                                CN2_Region[Chromosome][GC_Content[j][0]].append(coverage[j][-1])
                                    fcn2.close()
                                    if os.path.isfile(NullPath+temp_Name+'.fa'):
                                        os.system(r'''rm %s'''%(NullPath+temp_Name+'.fa'))
                                    if os.path.isfile(NullPath+temp_Name+'.sam'):
                                        os.system(r'''rm %s'''%(NullPath+temp_Name+'.sam'))
                                    Output_File=NullPath+'RD_Stat/'+bamF.split('/')[-1].replace('.bam','')+'.'+genome_name+'_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
                                    Output_Path=NullPath+'RD_Stat/'
                                    if not os.path.isdir(Output_Path):
                                        os.system(r'''mkdir %s'''%(Output_Path))
                                    fo=open(Output_File,'w')
                                    print >>fo, ' '.join(chromos)
                                    print >>fo, ' '.join([str(i) for i in range(101)])
                                    for key_1 in chromos:
                                            for key_2 in range(101):
                                                    print >>fo, ':'.join(['@',','.join(str(j) for j in CN2_Region[key_1][key_2])])
                                    fo.close()
                            time2=time.time()
                            print 'Null Model Built for '+bamF
                            print bamF+str(time2-time1)
    if function_name=='BPSearch':
        def IL_Prob_Bimodal(IL_Length,ILStats):
            log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std1']))-(float(IL_Length)-float(ILStats['bimodal']['Mean1']))**2/2/float(ILStats['bimodal']['std1'])**2
            log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std2']))-(float(IL_Length)-float(ILStats['bimodal']['Mean2']))**2/2/float(ILStats['bimodal']['std2'])**2
            return log_y1*float(ILStats['bimodal']['Bimodal1'])+log_y2*float(ILStats['bimodal']['Bimodal2'])
        def IL_Prob_Normal(IL_Length,ILStats):
            log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['normal']['std']))-(float(IL_Length)-float(ILStats['normal']['Mean']))**2/2/float(ILStats['normal']['std'])**2
            return log_y
        def log_factorial(K):
            k=int(K)
            k_list=[0]
            for i in range(k+1)[2:]:
                k_list.append(k_list[i-2]+numpy.log(i))
            return k_list[-1]
        def RD_Prob_NegativeBinomial(ReadDepth,RDStats):
            P=1-float(RDStats['Mean'])/(float(RDStats['STD'])**2)
            R=float(RDStats['Mean'])*(1-P)/P
            log_y=log_factorial(ReadDepth+R-1)-log_factorial(ReadDepth)-log_factorial(R-1)+R*numpy.log(1-P)+ReadDepth*numpy.log(P)
            return log_y
        def RD_Prob_Normal(ReadDepth,RDStats):
            log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(RDStats['STD']))-(float(ReadDepth)-float(RDStats['Mean']))**2/2/float(RDStats['STD'])**2
            return log_y
        def TB_Prob_Bimodal(TB_Length,TBStats):
            log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std1']))-(float(TB_Length)-float(TBStats['bimodal']['Mean1']))**2/2/float(TBStats['bimodal']['std1'])**2
            log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std2']))-(float(TB_Length)-float(TBStats['bimodal']['Mean2']))**2/2/float(TBStats['bimodal']['std2'])**2
            return log_y1*float(TBStats['bimodal']['Bimodal1'])+log_y2*float(TBStats['bimodal']['Bimodal2'])
        def TB_Prob_Normal(TB_Length,TBStats):
            log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['normal']['std']))-(float(TB_Length)-float(TBStats['normal']['Mean']))**2/2/float(TBStats['normal']['std'])**2
            return log_y
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
        def Reads_Direction_Detect(flag):
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
        def cigar2split(cigar):
            MapLen=[]
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            if cigars[0][1]=='S' or cigars[0][1]=='H':
                MapLen.append(0)
                cigars=cigars[1:]
            maplen=0
            for n in cigars:
                if n[1]=='M' or n[1]=='D' or n[1]=='N':
                    maplen+=int(n[0])
                if n[1]=='S' or n[1]=='H': 
                    MapLen.append(maplen)
            return MapLen
        def cigar2splitlen(cigar):
            MapLen=[]
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            for i in cigars:
                if i[1]=='S' or i[1]=='H':
                    MapLen.append(int(i[0]))
            return MapLen
        def cigar2splitqual(cigar,qual):
            import re
            pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
            cigars=[]
            qual_out=[]
            for m in pcigar.finditer(cigar):
                cigars.append((m.groups()[0],m.groups()[1]))
            pos1=[0]
            for n in cigars:
                if n[1]=='S' or n[1]=='H':
                    pos1.append(pos1[-1]+int(n[0]))
                    pos1.append(pos1[-1])
                elif n[1]=='M' or n[1]=='D' or n[1]=='N':
                    pos1[-1]+=int(n[0])
            pos1.remove(pos1[-1])
            qual2=[]
            for i in range(len(pos1)/2):
                qual3=qual[pos1[2*i]:pos1[2*i+1]]
                for j in qual3:
                    qual2.append(ord(j))
            if not qual2==[]:
                qual_out.append(numpy.mean(qual2))
            else:
                qual_out.append(30)
            return numpy.max(qual_out)
        def clusterNums(data, ClusterLen, direction):
            if data==[]:
                return [[],[]]
            else:
                data.sort()
                out2=[]
                out3=[]
                out=[]
                for i in data:
                    if len(out)==0:
                        out.append(i)
                    else:
                        if i-out[-1]< ClusterLen:
                            out.append(i)
                        else:
                            if out[-1]-out[0]<ClusterLen:
                                if direction=='f':
                                    out2.append(out[-1])
                                    out3.append(len(out))
                                elif direction=='r':
                                    out2.append(out[0])
                                    out3.append(len(out))
                            else:
                                temp=[0]
                                lenTime=int((out[-1]-out[0])/ClusterLen)
                                lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                                lenInter2=sorted(lenInter)[::-1][:lenTime]
                                for j in range(len(lenInter)):
                                    if lenInter[j] in  lenInter2:
                                        temp.append(j+1)
                                temp.append(len(out))
                                if direction=='f':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][-1])
                                        out3.append(temp[k+1]-temp[k])
                                elif direction=='r':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][0])
                                        out3.append(temp[k+1]-temp[k])
                            out=[i]
                if out[-1]-out[0]<ClusterLen:
                    if direction=='f':
                        out2.append(out[-1])
                        out3.append(len(out))
                    elif direction=='r':
                        out2.append(out[0])
                        out3.append(len(out))
                else:
                    temp=[0]
                    lenTime=int((out[-1]-out[0])/ClusterLen)
                    lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                    lenInter2=sorted(lenInter)[::-1][:lenTime]
                    for j in range(len(lenInter)):
                        if lenInter[j] in  lenInter2:
                            temp.append(j+1)
                    temp.append(len(out))
                    if direction=='f':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][-1])
                            out3.append(temp[k+1]-temp[k])
                    elif direction=='r':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][0])
                            out3.append(temp[k+1]-temp[k])
                out2a=[]
                out3a=[]
                for x in range(len(out2)):
                    if out3[x] <2: continue
                    else:
                        out2a.append(out2[x])
                        out3a.append(out3[x])
                return [out2a,out3a]
        def clusterNum1(data, ClusterLen, direction):
            if data==[]:
                return [[],[]]
            else:
                data.sort()
                out2=[]
                out3=[]
                out=[]
                for i in data:
                    if len(out)==0:
                        out.append(i)
                    else:
                        if i-out[-1]< ClusterLen:
                            out.append(i)
                        else:
                            if out[-1]-out[0]<ClusterLen:
                                if direction=='f':
                                    out2.append(out[-1])
                                    out3.append(len(out))
                                elif direction=='r':
                                    out2.append(out[0])
                                    out3.append(len(out))
                            else:
                                temp=[0]
                                lenTime=int((out[-1]-out[0])/ClusterLen)
                                lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                                lenInter2=sorted(lenInter)[::-1][:lenTime]
                                for j in range(len(lenInter)):
                                    if lenInter[j] in  lenInter2:
                                        temp.append(j+1)
                                temp.append(len(out))
                                if direction=='f':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][-1])
                                        out3.append(temp[k+1]-temp[k])
                                elif direction=='r':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][0])
                                        out3.append(temp[k+1]-temp[k])
                            out=[i]
                if out[-1]-out[0]<ClusterLen:
                    if direction=='f':
                        out2.append(out[-1])
                        out3.append(len(out))
                    elif direction=='r':
                        out2.append(out[0])
                        out3.append(len(out))
                else:
                    temp=[0]
                    lenTime=int((out[-1]-out[0])/ClusterLen)
                    lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                    lenInter2=sorted(lenInter)[::-1][:lenTime]
                    for j in range(len(lenInter)):
                        if lenInter[j] in  lenInter2:
                            temp.append(j+1)
                    temp.append(len(out))
                    if direction=='f':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][-1])
                            out3.append(temp[k+1]-temp[k])
                    elif direction=='r':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][0])
                            out3.append(temp[k+1]-temp[k])
                return [out2,out3]
        def clusterNum2(data1,data2,direction,numCff):
            data1.sort()
            out={}
            datatemp=[]
            if direction=='f':
                if data2[1][0]>numCff:
                    datatemp.append([0,data2[0][0]])
                for i in range(len(data2[1]))[1:]:
                    if data2[1][i]>numCff:
                        datatemp.append([data2[0][i-1],data2[0][i]])
            elif direction=='r':
                for i in range(len(data2[1])-1):
                    if data2[1][i]>numCff:
                        datatemp.append([data2[0][i],data2[0][i+1]])
                if data2[1][-1]>numCff:
                    datatemp.append([data2[0][-1],data2[0][-1]+1000])
            if len(datatemp)>0:
                flag1=0
                flag2=0
                if direction=='f':
                    while True:
                        if not data1[flag2] <datatemp[flag1][0] and not data1[flag2] >datatemp[flag1][1]:
                            if not datatemp[flag1][1] in out.keys():
                                out[datatemp[flag1][1]]=[data1[flag2]]
                            else:
                                out[datatemp[flag1][1]]+=[data1[flag2]]
                            flag2+=1
                        elif data1[flag2] >datatemp[flag1][1]:
                            flag1+=1
                        elif data1[flag2]<datatemp[flag1][0]:
                            flag2+=1
                        if flag1==len(datatemp) or flag2==len(data1): break
                if direction=='r':
                    while True:
                        if not data1[flag2] <datatemp[flag1][0] and  not data1[flag2] >datatemp[flag1][1]:
                            if not datatemp[flag1][0] in out.keys():
                                out[datatemp[flag1][0]]=[data1[flag2]]
                            else:
                                out[datatemp[flag1][0]]+=[data1[flag2]]
                            flag2+=1
                        elif data1[flag2] >datatemp[flag1][1]:
                            flag1+=1
                        elif data1[flag2]<datatemp[flag1][0]:
                            flag2+=1
                        if flag1==len(datatemp) or flag2==len(data1): break
                return out
            else:
                return {}
        def clusterNum3(hash1,hash2):
            hash3={}
            for key in hash2.keys():
                hash3[key]=[0,0,0]
                for key2 in hash2[key]:
                    hash3[key][0]+=hash1[key2][0]
                    hash3[key][1]+=hash1[key2][1]
                    hash3[key][2]+=hash1[key2][2]
            return hash3
        def clusterNums4(hash, ClusterLen,direction):
            data1=clusterNum1(hash.keys(),ClusterLen,direction)[0]
            data2=[]
            if not data1==[]:
                if direction=='r':
                    data1=data1[::-1]
                temp=[[] for i in data1]
                if direction=='f':
                    rec_temp=0
                    for k1 in sorted(hash.keys()):
                        if k1<data1[rec_temp]:
                            temp[rec_temp].append(k1)
                        elif k1==data1[rec_temp]:
                            temp[rec_temp].append(k1)
                            rec_temp+=1
                elif direction=='r':
                    rec_temp=0
                    for k1 in sorted(hash.keys())[::-1]:
                        if k1>data1[rec_temp]:
                            temp[rec_temp].append(k1)
                        elif k1==data1[rec_temp]:
                            temp[rec_temp].append(k1)
                            rec_temp+=1
                data2=[0 for i in data1]
                rec_temp=-1
                for k1 in temp:
                    rec_temp+=1
                    for k2 in k1:
                        data2[rec_temp]+=hash[k2]
            return [data1,data2]
        def clusterNumSP(hash,ClusterLen,numCff):
            data=hash.keys()
            data.sort()
            out1=clusterNums(data, ClusterLen, 'f')
            out2=clusterNums(data, ClusterLen, 'r')
            out={}
            flag1=0
            flag2=0
            while True:
                if data[flag2]< out2[0][flag1]:
                    flag2+=1
                elif not data[flag2]< out2[0][flag1] and not data[flag2]> out1[0][flag1]:
                    if not out1[0][flag1] in out.keys():
                        out[out1[0][flag1]]=[out2[0][flag1],hash[data[flag2]]]
                    else:
                        out[out1[0][flag1]][1]+=hash[data[flag2]]
                    flag2+=1
                elif data[flag2]> out1[0][flag1]:
                    flag1+=1
                if flag1==len(out1[0]) or flag2==len(hash): break
            for i in out.keys():
                if not out[i][1]>numCff:
                    del out[i]
            return out 
        def clusterNumLN(hash,ClusterLen,numCff):
            out={}
            for key in hash.keys():
                out[key]={}
                key2={}
                for key3 in hash[key]:
                    if not key3[1] in key2.keys():
                        key2[key3[1]]=key3[0]
                    else:
                        key2[key3[1]]+=key3[0]
                for key4 in key2.keys():
                    if len(key2[key4])>numCff:
                        if key4=='ff' or key4=='rf':
                            temp1=clusterNums(key2[key4],ClusterLen,'f')
                            for temp2 in range(len(temp1[1])):
                                if temp1[1][temp2]>numCff:
                                    if not key4 in out[key].keys():
                                        out[key][key4]=temp1[0][temp2]
                        elif key4=='rr' or key4=='fr':
                            temp1=clusterNums(key2[key4],ClusterLen,'r')
                            for temp2 in range(len(temp1[1])):
                                if temp1[1][temp2]>numCff:
                                    if not key4 in out[key].keys():
                                        out[key][key4]=temp1[0][temp2]
                        elif key4==0:
                            if len(key2[key4])>numCff:
                                out[key][key4]=0
                        else:
                            if key4[-1]=='+':
                                temp=clusterNums(key2[key4],ClusterLen,'f')
                            else:
                                temp=clusterNums(key2[key4],ClusterLen,'r')
                            for temp2 in range(len(temp[1])):
                                if temp[1][temp2]>numCff:
                                    if not key4 in out[key].keys():
                                        if key4[-1]=='+':
                                            out[key][key4]=[temp[0],'f']
                                        elif key4[-1]=='-':
                                            out[key][key4]=[temp[0],'r']
            return out
        def clusterQC(hash, CluNum):
            out=[]
            for i in range(len(hash[1])):
                if not hash[1][i]<CluNum:
                    out.append(hash[0][i])
            return out
        def clusterSupVis2(dataRound, dataLeft, dataRight,keypick):
            flag1=0
            flag2=0
            out={}
            if not dataLeft==[] and not dataRight==[]:
                if keypick=='left':
                    while True:
                        if flag1==len(dataRound): break
                        else:
                            if flag2==len(dataLeft)-1 and dataRound[flag1]>dataRight[flag2]: break
                            else:
                                if dataRound[flag1]<dataLeft[flag2]:
                                    flag1+=1
                                elif dataRound[flag1]>dataRight[flag2]:
                                    flag2+=1
                                else:
                                    if not dataLeft[flag2] in out.keys():
                                        out[dataLeft[flag2]]=[dataRound[flag1]]
                                    else:
                                        out[dataLeft[flag2]]+=[dataRound[flag1]]
                                    flag1+=1
                if keypick=='right':
                    while True:
                        if flag1==len(dataRound): break
                        else:
                            if flag2==len(dataLeft)-1 and dataRound[flag1]>dataRight[flag2]: break
                            else:
                                if dataRound[flag1]<dataLeft[flag2]:
                                    flag1+=1
                                elif dataRound[flag1]>dataRight[flag2]:
                                    flag2+=1
                                else:
                                    if not dataRight[flag2] in out.keys():
                                        out[dataRight[flag2]]=[dataRound[flag1]]
                                    else:
                                        out[dataRight[flag2]]+=[dataRound[flag1]]
                                    flag1+=1
            return out
        def clusterSupVis(dataRou,dataCen,ClusterLen):
            flag1=0
            flag2=0
            dataRou.sort()
            dataCen.sort()
            out={}
            while True:
                if flag1==len(dataRou): break
                else:
                    if flag2==len(dataCen)-1:
                        if dataRou[flag1]>dataCen[flag2]+ClusterLen: break
                        elif dataRou[flag1]<dataCen[flag2]-ClusterLen: 
                            flag1+=1
                        else:
                            if not dataCen[flag2] in out.keys():
                                out[dataCen[flag2]]=[dataRou[flag1]]
                            else:
                                out[dataCen[flag2]]+=[dataRou[flag1]]
                            flag1+=1
                    else:
                        if dataRou[flag1]<dataCen[flag2]-ClusterLen:
                            flag1+=1
                        elif dataRou[flag1]>dataCen[flag2]+ClusterLen:
                            flag2+=1
                        else:
                            if not dataCen[flag2] in out.keys():
                                out[dataCen[flag2]]=[dataRou[flag1]]
                            else:
                                out[dataCen[flag2]]+=[dataRou[flag1]]
                            flag1+=1
            return out
        def clusterLN(data, LnHash, numCff):
            out={}
            for i in data:
                if not i in LnHash.keys():
                    continue
                else:
                    j=LnHash[i]
                    key=j[0]+j[2]
                    if not key in out.keys():
                        out[key]=[j[1]]
                    else:
                        out[key]+=[j[1]]
            for key in out.keys():
                if len(out[key])<numCff:
                    del out[key]
            return out
        def clusterSupVis3(dataRound,dataCen):
            out={}
            for key1 in dataCen:
                min_num=abs(dataRound[0]-key1)
                min_rec=dataRound[0]
                for key2 in dataRound[1:]:
                    alt_num=abs(key2-key1)
                    alt_rec=key2
                    if alt_num<min_num:
                        min_num=alt_num
                        min_rec=alt_rec
                out[key1]=min_rec
            return out
        def BPFilter_ReadNumThroughBP(BamInput,chrom,bpCoordinate,flank):
            flank=int(flank)
            bps=[bpCoordinate-flank,bpCoordinate+flank]
            if Bam_Flag==1:
                fin=os.popen(r'''samtools view %s chr%s:%d-%d'''%(BamInput,chrom,int(bps[0]),int(bps[1])))
            else:
                fin=os.popen(r'''samtools view %s %s:%d-%d'''%(BamInput,chrom,int(bps[0]),int(bps[1])))
            RDs=[0*i for i in range(bps[1]-bps[0])]
            while True:
                pbam=fin.readline().strip().split()
                if not pbam:break
                if int(pbam[4])<int(QCAlign): continue 
                if int(pbam[1])&4>0: continue
                ReadLen=cigar2reaadlength(pbam[5])
                if ReadLen<90:
                    pos_pbam=[int(pbam[3]),int(pbam[3])+ReadLen]
                else:
                    diff_pbam=(ReadLen-90)/2
                    pos_pbam=[int(pbam[3])+diff_pbam,int(pbam[3])+ReadLen-diff_pbam]
                for j in range(pos_pbam[0]-bps[0],pos_pbam[1]-bps[0]+1):
                    if j in range(len(RDs)):
                        RDs[j]+=1
            fin.close()
            min_through=float('inf')
            for i in range(len(RDs)):
                if RDs[i]<min_through:
                    min_through=RDs[i]
                    min_record=i
            if min_through<RD_Mean/2:
                return min_record+int(bps[0])
            else:
                return -1
        def clusterSupVis4(dataRou,left,right):
            out=[]
            for i in sorted(dataRou):
                if i<left:
                    continue
                if i>right:
                    break
                else:
                    out.append(i)
            return out
        def overlap_calcu(list1,list2):
            list1.sort()
            list2.sort()
            list_len=min(len(list1),len(list2))
            rec=0
            reca=0
            recb=0
            while True:
                if reca==len(list1) or recb==len(list2): break
                if list1[reca]==list2[recb]:
                    rec+=1
                    reca+=1
                    recb+=1
                elif list1[reca]<list2[recb]:
                    reca+=1
                elif list1[reca]>list2[recb]:
                    recb+=1
            return float(rec)/float(list_len)
        def lower_ILCff_calcu(IL_Null_Stat,CICff):
            fin=open(IL_Null_Stat)
            flag=0
            outhash={}
            totalnum=0
            for line in fin:
                pin=line.strip().split()
                if pin[0]=='ReadDepth' and  pin[1]=='Frequency':
                    flag=0
                if flag==1:
                    outhash[int(pin[0])]=int(pin[1])
                    totalnum+=int(pin[1])
                if pin[0]=='InsertLength' and pin[1]=='Frequency':
                    flag=1
            fin.close()
            subnum=0
            for k1 in sorted(outhash.keys()):
                subnum+=outhash[k1]
                if float(subnum)/float(totalnum)>CICff/2:
                    break
            return k1
        def bp_subgroup(list, Min_Distinguish_Len):
            if list==[]:
                return []
            else:
                out=[[list[0]]]
                list.sort()
                for x in range(len(list)-1):
                    if list[x+1]-list[x]>Min_Distinguish_Len:
                        out.append([list[x+1]])
                    else:
                        out[-1].append(list[x+1])
                return out
        def Fasta_To_Sequence(fasta_file):
            ffa=open(fasta_file)
            ffa.readline().strip().split()
            sequence=[]
            while True:
                seq=ffa.readline().strip().split()
                if not seq: break
                sequence.append(seq)
            seq2=sequence[0][0]
            for seq3 in sequence[1:]:
                seq2+=seq3[0]
                if ''.join(['N' for x in range(100)]) in seq2: break
            if ''.join(['N' for x in range(100)]) in seq2:
                return 'ERROR!'
            else:
                return seq2
        def file_start(filename):
            fRDind=open(filename,'w')
            fRDind.close()
        def GC_Content_Calculate(seq2):
            NumAT=0
            NumGC=0
            for n in seq2:
                if n=='A' or n=='T':
                    NumAT+=1
                elif n=='G' or n=='C':
                    NumGC+=1
            return [NumGC,NumAT]            
        def path_mkdir(path):
            if not os.path.isdir(path):
                    os.system(r'''mkdir %s'''%(path))
        def stat_file_name(bamF_Name,genome_name):
            global ILStat
            ILStat=NullPath+'ILNull.'+bamF_Name+'.'+genome_name+'.Bimodal'
            global IL_Null_Stat
            IL_Null_Stat=NullPath+bamF_Name+'.'+genome_name+'.density.null'
            global RDStat
            RDStat=NullPath+'RDNull.'+bamF_Name+'.'+genome_name+'.NegativeBinomial'
            global TBStat
            TBStat=NullPath+'TBNull.'+bamF_Name+'.'+genome_name+'.Bimodal'
            global SPLenStat
            SPLenStat=NullPath+bamF_Name+'.'+genome_name+'.SplitLength'
            global AllStat
            AllStat=NullPath+bamF_Name+'.'+genome_name+'.Stats'
        def SPLCff_Calculate(NullSplitLen_perc,SPLenStat):
            if NullSplitLen_perc>1:
                SPLCff=NullSplitLen_perc
            else:        
                fSPLStat=open(SPLenStat)
                fSPLStat.readline()
                SPLength={}
                totalSPlit=0
                subSPlit=0 
                while True:
                    pSPLStat=fSPLStat.readline().strip().split()
                    if not pSPLStat: break
                    SPLength[int(pSPLStat[0])]=int(pSPLStat[1])
                    totalSPlit+=int(pSPLStat[1])
                if not SPLength.keys()==[]:
                    for key in SPLength.keys():
                            subSPlit+=SPLength[key]
                            if float(subSPlit)/float(totalSPlit)>NullSplitLen_perc: 
                                    break
                    SPLCff=key-1
                    if SPLCff<ReadLength/10:
                        SPLCff=ReadLength/10
                else:
                    SPLCff=ReadLength/10
                if SPLCff<10:
                    SPLCff=10
                fSPLStat.close()
                return SPLCff
        def Define_Default_BPSearching():
            global ILCutoff
            ILCutoff=0.95
            global RDCutoff
            RDCutoff=0.95
            global TBCutoff
            TBCutoff=0.9999
            global SplitCutoff
            SplitCutoff=0.999
            global ABILCutoff
            ABILCutoff=0.99
            global DRCutoff
            DRCutoff=0.99
            global SPLCutoff
            SPLCutoff=0.85
            global Length_Limit
            Length_Limit=2000
            global model_comp
            if not '--null-model' in dict_opts.keys():
                model_comp='S'
            else:
                if dict_opts['--null-model'] in ['S','Simple']:
                    model_comp='S'
                else:
                    model_comp='C'
            global ToolMappingQ
            global FileMappingQ
            global align_QCflag
            if '--qc-map-tool' in dict_opts.keys() and '--qc-map-file' in dict_opts.keys():
                ToolMappingQ=dict_opts['--qc-map-tool']
                FileMappingQ=dict_opts['--qc-map-file']
                align_QCflag=1
            else:
                align_QCflag=0
            global BPAlignQC
            if '--BPAlignQC' in dict_opts.keys():
                BPAlignQC=float(dict_opts['--BPAlignQC'])
            else:
                BPAlignQC=0.0
            global QCAlign   
            if '--qc-align' in dict_opts.keys():
                QCAlign=int(dict_opts['--qc-align'])
            else:
                QCAlign=20
            global QCSplit   
            if '--qc-split' in dict_opts.keys():
                QCSplit=int(dict_opts['--qc-split'])
            else:
                QCSplit=20
            global NullSplitLen_perc
            if '--split-min-len' in dict_opts.keys():
                NullSplitLen_perc=float(dict_opts['--split-min-len'])
            else:
                NullSplitLen_perc=0.9
            global BPAlignQCFlank
            if '--BPAlignQCFlank' in dict_opts.keys():
                BPAlignQCFlank=int(dict_opts['--BPAlignQCFlank'])
            else:
                BPAlignQCFlank=500
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def chromos_read_in(ref_index):
            chromos=[]
            fin=open(ref_index)
            for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
            fin.close()
            return chromos
        def genome_name_readin():
            if '--NullGenomeName' in dict_opts.keys():
                genome_name=dict_opts['--NullGenomeName']
            else:
                genome_name='genome'
            return genome_name
        def ILStats_readin(ILStat):
            ILStats={}
            fILS=open(ILStat)
            if model_comp=='C':
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                ILStats['stat']={}
                for i in range(len(pILS1)):
                        ILStats['stat'][pILS1[i]]=pILS2[i]
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                ILStats['bimodal']={}
                for i in range(len(pILS1)):
                        ILStats['bimodal'][pILS1[i]]=pILS2[i]
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                for i in range(len(pILS1)):
                        ILStats['bimodal'][pILS1[i]]=pILS2[i]
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                ILStats['normal']={}
                for i in range(len(pILS1)):
                        ILStats['normal'][pILS1[i]]=pILS2[i]
            elif model_comp=='S':
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                ILStats['stat']={}
                for i in range(len(pILS1)):
                        ILStats['stat'][pILS1[i]]=pILS2[i]
                pILS1=fILS.readline().strip().split()
                pILS2=fILS.readline().strip().split()
                ILStats['normal']={}
                for i in range(len(pILS1)):
                        ILStats['normal'][pILS1[i]]=pILS2[i]
                ILStats['bimodal']={}
                for i in range(len(pILS1)):
                        ILStats['bimodal'][pILS1[i]]=pILS2[i]
                for i in range(len(pILS1)):
                        ILStats['bimodal'][pILS1[i]]=pILS2[i]
            fILS.close() 
            return ILStats
        def RDStats_readin(RDStat):    
            RDStats={}
            fRDS=open(RDStat)
            pRDS1=fRDS.readline().strip().split()
            pRDS2=fRDS.readline().strip().split()
            for i in range(len(pRDS1)):
                    RDStats[pRDS1[i]]=pRDS2[i]
            fRDS.close()
            return RDStats
        def TBStats_readin(TBStat):
            TBStats={}
            fTBS=open(TBStat)
            if model_comp=='C':
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                TBStats['stat']={}
                for i in range(len(pTBS1)):
                        TBStats['stat'][pTBS1[i]]=pTBS2[i]
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                TBStats['bimodal']={}
                for i in range(len(pTBS1)):
                        TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                for i in range(len(pTBS1)):
                        TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                TBStats['normal']={}
                for i in range(len(pTBS1)):
                        TBStats['normal'][pTBS1[i]]=pTBS2[i]
            elif model_comp=='S':
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                TBStats['stat']={}
                for i in range(len(pTBS1)):
                        TBStats['stat'][pTBS1[i]]=pTBS2[i]
                pTBS1=fTBS.readline().strip().split()
                pTBS2=fTBS.readline().strip().split()
                TBStats['normal']={}
                for i in range(len(pTBS1)):
                        TBStats['normal'][pTBS1[i]]=pTBS2[i]
                TBStats['bimodal']={}
                for i in range(len(pTBS1)):
                        TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
                for i in range(len(pTBS1)):
                        TBStats['bimodal'][pTBS1[i]]=pTBS2[i]
            fTBS.close()
            return TBStats
        def Null_Stats_Readin_One(NullPath,bamF,NullSplitLen_perc):
            fin=open(NullPath+bamF.split('/')[-1].replace('.bam','.')+genome_name+'.Stats')
            pin=fin.readline().strip().split()
            pin=fin.readline().strip().split()
            pin=fin.readline().strip().split()
            global Window_Size
            global ReadLength
            global sub_loc_size
            Window_Size=int(float(pin[0])/3)
            sub_loc_size=max(int(float(pin[3]))*100,10**6)
            for line in fin:
                pin=line.strip().split()
            ReadLength=int(pin[-1].split(':')[1])
            if NullSplitLen_perc>1:
                NullSplitLen_perc=float(NullSplitLen_perc)/float(ReadLength)
            fin.close()
        def split_loc_to_subloc(loc,sub_loc_size):
            if not loc[1]-loc[0]>2*sub_loc_size:
                loc2=[loc]
            else: 
                sublocNum=(loc[1]-loc[0])/sub_loc_size
                loc2=[[loc[0],loc[0]+sub_loc_size+int(ClusterLen2)]]
                for slNum in range(sublocNum)[1:]:
                    loc2.append([loc[0]+slNum*sub_loc_size-int(ClusterLen2),loc[0]+(slNum+1)*sub_loc_size+int(ClusterLen2)])
                if sublocNum>1:
                    loc2.append([loc[0]+(slNum+1)*sub_loc_size-int(ClusterLen2),loc[1]])
            return loc2
        def signIL_Decide(absIL,ILCffs):
            signIL=0
            if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                signIL+=1
            return signIL
        def pos_define_and_redefine(pbam,ReadLen,real_region,Window_Size):
            pos1=int(pbam[3])
            pos2=int(pbam[3])+ReadLen
            if pos2>real_region[0] and pos1<real_region[1]:
                if pos1<real_region[0] and pos2>real_region[0]:
                    pos1=real_region[0]
                if pos1<real_region[1] and pos2>real_region[1]:
                    pos2=real_region[1]
                block1=(pos1-real_region[0])/Window_Size
                RD_RealRegion[block1]+=ReadLen
            return [pos1,pos2]
        def QCFlag_Define(absIL,ILCffs,DRtemp,pbam,pbam1,fmini):
            QCFlag=0
            if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                QCFlag+=1
            if DRtemp==['+','-'] and int(pbam[8])<0:
                QCFlag+=1
            if DRtemp==['+','+'] or DRtemp==['-','-']:
                QCFlag+=1
            if int(pbam[8])==0:
               QCFlag+=1
            if not pbam[5].find('S')==-1:
                QCFlag+=1
            if not QCFlag==0:
                print>>fmini, pbam1
        def tempIL_Info_Add(clu_c_F,LinkCluMin):
            if not clu_c_F=={}:
                    for key2 in clu_c_F.keys():
                            key2b=key2-10
                            record=0
                            record2=0
                            for key3 in clu_c_F[key2]:
                                    if not key3 >key2b:
                                            record+=sum(abInfF[key3][:3])
                                    if abs(key2b-key3)<SPCluLen:
                                            record2+=abInfF[key3][3]
                            if not record+record2<LinkCluMin:
                                    tempIL[key2b]=[[record,record2,'f'],clu_c_F[key2]]
                            del clu_c_F[key2]
        def tempIL_Filter_One():
            if not tempIL=={}:
                for aqb in tempIL.keys():
                    if aqb<BPAlignQCFlank:
                        del tempIL[aqb]
                        continue
                    if align_QCflag==1:
                        LNSPFR_aqb=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,aqb-BPAlignQCFlank,aqb+BPAlignQCFlank))
                        tPairF_b=float(LNSPFR_aqb.read().strip()) 
                        if tPairF_b<float(BPAlignQC):
                            del tempIL[aqb]
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        CN2_Region={}
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--sample, input alignment file in bam format'
            print 'Optional Parameters:'
            print '--chromosome, name of chromosome to run. should match chromosome name in bam file'
            print '--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
            print '--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)'
            print '--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)'
            print '--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)'
            print '--qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary'
            print '--qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/'
            print '--qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported (default: 0.0)'
        else:
            Define_Default_BPSearching()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                if not '--sample' in dict_opts.keys():
                    print 'Error: please specify either input file using --sample'
                else:
                    if '--sample' in dict_opts.keys():
                        bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                        bam_files=[dict_opts['--sample']]
                    else:
                        bam_path=path_modify(dict_opts['--PathBam'])
                        bam_files=[]
                        for file1 in os.listdir(bam_path):
                            if file1.split('.')[-1]=='bam':
                                bam_files.append(bam_path+file1)
                    ref_path=workdir+'reference/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        chromos=chromos_read_in(ref_index)
                        if '--chromosome' in dict_opts.keys():
                            chrom_single=dict_opts['--chromosome']
                            if not chrom_single in chromos:
                                print 'Error: please make sure the chromosome defined by --chr is correct based on the reference genome'
                                chromos=[]
                            else:
                                chromos=[chrom_single]
                        if not chromos==[]:
                            genome_name=genome_name_readin()
                            out_path=workdir
                            print 'temp files produced under: '+workdir
                            BPPath=out_path+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
                            if not os.path.isdir(BPPath):
                                os.system(r'''mkdir %s'''%(BPPath))
                            NullPath=out_path+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/'
                            for bamF in bam_files:
                                time1=time.time()
                                Null_Stats_Readin_One(NullPath,bamF,NullSplitLen_perc)
                                for chrF in chromos:
                                    bamF_Name=bamF.split('/')[-1].replace('.bam','')
                                    floc_Name=BPPath+bamF_Name+'.'+chrF
                                    Refloc_name='.'.join(ref_file.split('.')[:-1])+'.Mappable.'+chrF+'.bed'
                                    stat_file_name(bamF_Name,genome_name)
                                    if os.path.isfile(Refloc_name):
                                        BamInput=bamF
                                        ILStats=ILStats_readin(ILStat)
                                        RDStats=RDStats_readin(RDStat)
                                        TBStats=TBStats_readin(TBStat)
                                        SPLCff=SPLCff_Calculate(NullSplitLen_perc,SPLenStat)
                                        fAllS=open(AllStat)
                                        pAllS=fAllS.readline().rstrip()
                                        ILCIs={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                                ILCIs[pAllS1[i]]=pAllS2[i]
                                        pAllS=fAllS.readline().rstrip()
                                        RDCIs={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                            RDCIs[pAllS1[i]]=pAllS2[i]
                                        pAllS=fAllS.readline().rstrip()
                                        TBCIs={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                            TBCIs[pAllS1[i]]=pAllS2[i]
                                        pAllS=fAllS.readline().rstrip()
                                        InsertLen={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                                InsertLen[pAllS1[i]]=pAllS2[i]
                                        for i in range(len(pAllS1)):
                                                if float(pAllS2[i])>ABILCutoff:
                                                        InsertLenMin=int(pAllS1[i])-1
                                                        break
                                        if InsertLenMin<5:
                                                InsertLenMin=5
                                        pAllS=fAllS.readline().rstrip()
                                        SplitReads={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                                SplitReads[pAllS1[i]]=pAllS2[i]
                                        for i in range(len(pAllS1)):
                                                if float(pAllS2[i])>SplitCutoff:
                                                        SplitMin=int(pAllS1[i])-1
                                                        break
                                        if SplitMin<5:
                                                SplitMin=5
                                        pAllS=fAllS.readline().rstrip()
                                        AbDirection={}
                                        pAllS1=fAllS.readline().strip().split()
                                        pAllS2=fAllS.readline().strip().split()
                                        for i in range(len(pAllS1)):
                                                AbDirection[pAllS1[i]]=pAllS2[i]
                                        for i in range(len(pAllS1)):
                                            if float(pAllS2[i])>DRCutoff:
                                                DRMin=int(pAllS1[i])-1
                                                break
                                        if DRMin<5:
                                                DRMin=5
                                        fbam=os.popen('''samtools view -H %s'''%(BamInput))
                                        if '--BPSPCff' in dict_opts.keys():
                                            BPSPCff=int(float(dict_opts['--BPSPCff']))
                                        else:
                                            BPSPCff=int(round(1.2*float(TBStats['stat']['Median'])/float(10)))
                                        if BPSPCff<3:
                                            BPSPCff=3
                                        if '--BPLNCff' in dict_opts.keys():
                                            BPLNCff=int(float(dict_opts['--BPLNCff']))
                                        else:
                                            BPLNCff=int(round(1.8*float(TBStats['stat']['Median'])/float(10)))
                                        if BPLNCff<3:
                                            BPLNCff=3
                                        SPCluMin=BPSPCff
                                        LnCluMin=BPLNCff
                                        if '--SPCluLen' in dict_opts.keys():
                                            SPCluLen=int(dict_opts['--SPCluLen'])
                                        else:
                                            SPCluLen=5
                                        SubLnCluMin=max([LnCluMin,SPCluMin])
                                        LinkCluMin=min([LnCluMin,SPCluMin])
                                        ClusterLen=float(ILStats['normal']['Mean'])+2*float(ILStats['normal']['STD'])-ReadLength
                                        ClusterLen2=int(ClusterLen/Window_Size+1)*Window_Size
                                        Min_Distinguish_Len=Window_Size
                                        subLnClusterLen=ClusterLen/2
                                        if not '-S' in dict_opts.keys():
                                            CICff=0.003
                                            ILCffs=[lower_ILCff_calcu(IL_Null_Stat,CICff),int(float(ILStats['normal']['Mean'])+3*float(ILStats['normal']['STD']))]
                                        else:
                                            CICff=0.01/float(dict_opts['-S'])
                                            ILCffs=[lower_ILCff_calcu(IL_Null_Stat,CICff),int(float(ILStats['normal']['Mean'])+float(dict_opts['-S'])*float(ILStats['normal']['STD']))]
                                        BPOutputa=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.SPs'
                                        fout=open(BPOutputa,'w')
                                        fout.close()
                                        BPOutputb=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.LNs'
                                        fout=open(BPOutputb,'w')
                                        fout.close()
                                        BPOutputd=floc_Name+'.'+'.'.join(['SPCff'+str(SPCluMin),'CluCff'+str(LnCluMin),'AlignCff'+str(BPAlignQC)])+'.chromLNs'
                                        fout=open(BPOutputd,'w')
                                        fout.close()
                                        BPOutpute='/'.join(BPOutputd.split('/')[:-1])+'/'+'.'.join(BPOutputd.split('/')[-1].split('.')[:-6]+BPOutputd.split('/')[-1].split('.')[-5:])
                                        if not os.path.isfile(BPOutpute):
                                            fout=open(BPOutpute,'w')
                                            fout.close()
                                        abtLink={}
                                        floc=open(Refloc_name)
                                        loc_rec={}
                                        for line in floc:
                                            ploc=line.strip().split()
                                            if not ploc[0] in loc_rec.keys():
                                                loc_rec[ploc[0]]=[]
                                                loc_rec[ploc[0]].append([int(ploc2) for ploc2 in ploc[1:]])
                                            else:
                                                loc_rec[ploc[0]].append([int(ploc2) for ploc2 in ploc[1:]])
                                        floc.close()
                                        if not loc_rec=={}:
                                            test_mul_RP=[]
                                            test_mul_SP=[]
                                            chrom=chrF
                                            mini_fout_Name=BPPath+bamF_Name+'.mini.'+loc_rec.keys()[0]+'.sam'
                                            mini_fout_N2=BPPath+bamF_Name+'.mini.'+chrom+'.bam'
                                            mini_fout_N3=BPPath+bamF_Name+'.mini.'+chrom+'.sorted'
                                            mini_fout_N4=BPPath+bamF_Name+'.mini.'+chrom+'.sorted.bam'
                                            if not os.path.isfile(mini_fout_N4):
                                                os.system(r''' samtools view -H %s -o %s'''%(BamInput,mini_fout_Name)) 
                                                RD_index_Path=NullPath+'RD_Stat/'
                                                if not os.path.isdir(RD_index_Path):
                                                    os.system(r'''mkdir %s'''%(RD_index_Path))
                                                RD_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.RD.index'
                                                fRDind=open(RD_index_File,'w')
                                                fRDind.close()
                                                for loc in loc_rec[chrom]:
                                                    loc2=split_loc_to_subloc(loc,sub_loc_size)
                                                    for real_region in loc2:
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/Window_Size+1)]
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            pos1=int(pbam[3])
                                                            pos2=int(pbam[3])+ReadLen
                                                            if pos2>real_region[0] and pos1<real_region[1]:
                                                                if pos1<real_region[0] and pos2>real_region[0]:
                                                                    pos1=real_region[0]
                                                                if pos1<real_region[1] and pos2>real_region[1]:
                                                                    pos2=real_region[1]
                                                                block1=(pos1-real_region[0])/Window_Size
                                                                RD_RealRegion[block1]+=ReadLen
                                                            if int(pbam[4])>int(QCAlign):             #fail quality control, skip
                                                                absIL=abs(int(pbam[8]))
                                                                QCFlag=0
                                                                link_flag=0
                                                                if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    QCFlag+=1
                                                                if DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    QCFlag+=1
                                                                if DRtemp==['+','+'] or DRtemp==['-','-']:
                                                                    QCFlag+=1
                                                                if int(pbam[8])==0:
                                                                   QCFlag+=1
                                                                if not pbam[5].find('S')==-1:
                                                                    QCFlag+=1
                                                                if not QCFlag==0:
                                                                    print>>fmini, pbam1
                                                        fbam.close()
                                                        fmini.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/Window_Size)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/Window_Size*Window_Size ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/Window_Size*Window_Size))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()
                                                os.system(r'''samtools view -h -Sb %s -o %s'''%(mini_fout_Name,mini_fout_N2))
                                                os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
                                                os.system(r'''samtools index %s '''%(mini_fout_N4))
                                                os.system(r'''rm %s'''%(mini_fout_N2))
                                                os.system(r'''rm %s'''%(mini_fout_Name))
                                            temp_IL_Rec={}
                                            Link_IL_Rec={}
                                            for loc in loc_rec[chrom]: 
                                                loc2=split_loc_to_subloc(loc,sub_loc_size)
                                                for real_region in loc2:
                                                    if real_region[1]-real_region[0]<Window_Size:continue
                                                    else:
                                                        fbam=os.popen('''samtools view %s %s:%d-%d'''%(mini_fout_N4,chrom,real_region[0],real_region[1]))
                                                        abInfF={}
                                                        abInfR={}
                                                        abLink={}
                                                        abInf3={}
                                                        LinkFR={}
                                                        LinkFF={}
                                                        LinkRR={}
                                                        LinkRF={}
                                                        LinkNM={}
                                                        LinkNM['F']=[]
                                                        LinkNM['R']=[]
                                                        LinkSP={}
                                                        abLinkSP={}
                                                        test=[]
                                                        LinkSPTemp={}
                                                        while True:
                                                            pbam=fbam.readline().strip().split()
                                                            if not pbam: break
                                                            DRtemp=Reads_Direction_Detect(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            absIL=abs(int(pbam[8]))
                                                            signIL=signIL_Decide(absIL,ILCffs)
                                                            posF=[]
                                                            if not pbam[5].find('S')==-1 or not pbam[5].find('H')==-1:
                                                                ClippedLen=cigar2splitlen(pbam[5])
                                                                ClippedPos=cigar2split(pbam[5])
                                                                ClippedQual=cigar2splitqual(pbam[5],pbam[10])
                                                                if ClippedQual>QCSplit:
                                                                        ClipAbsPos=[]
                                                                        for c in range(len(ClippedLen)):
                                                                            if ClippedLen[c]>SPLCff:
                                                                                    pos1=int(pbam[3])+ClippedPos[c]
                                                                                    posF.append(pos1)
                                                                                    pos2=int(pbam[7])
                                                                                    if DRtemp[0]=='+':
                                                                                        if not pos1 in abInfF.keys():
                                                                                            abInfF[pos1]=[0,0,0,1]
                                                                                        else:
                                                                                            abInfF[pos1][3]+=1
                                                                                    elif DRtemp[0]=='-':
                                                                                        if not pos1 in abInfR.keys():
                                                                                            abInfR[pos1]=[0,0,0,1]
                                                                                        else:
                                                                                            abInfR[pos1][3]+=1
                                                                if not pbam[0] in LinkSPTemp.keys():
                                                                    LinkSPTemp[pbam[0]]=[[]]
                                                                else:
                                                                    LinkSPTemp[pbam[0]]+=[[]]
                                                                for c in ClippedPos:
                                                                    pos1=int(pbam[3])+c
                                                                    LinkSPTemp[pbam[0]][-1]+=[pos1,DRtemp[0]]
                                                                LinkSPTemp[pbam[0]].append('S')
                                                            if posF==[]:
                                                                if DRtemp[0]=='+':
                                                                    pos1=int(pbam[3])+len(pbam[9])
                                                                elif DRtemp[0]=='-':
                                                                    pos1=int(pbam[3])
                                                            else:
                                                                pos1=posF[0]
                                                            pos2=int(pbam[7])
                                                            if pbam[0] in LinkSPTemp.keys():
                                                                LinkSPTemp[pbam[0]]+=[[pos1,DRtemp[0]]]
                                                            else:
                                                                LinkSPTemp[pbam[0]]=[[pos1,DRtemp[0]]]
                                                            if DRtemp==['+','-'] and int(pbam[8])>0:
                                                                    if signIL>0:
                                                                            if not pos1 in LinkFR.keys():
                                                                                    LinkFR[pos1]=[pos2]
                                                                            else:
                                                                                    LinkFR[pos1]+=[pos2]
                                                                            if not pos1 in abInfF.keys():
                                                                                    abInfF[pos1]=[1,0,0,0]
                                                                            else:
                                                                                    abInfF[pos1][0]+=1
                                                                            if not pos2 in abInfR.keys():
                                                                               abInfR[pos2]=[1,0,0,0]
                                                                            else:
                                                                                    abInfR[pos2][0]+=1
                                                            elif DRtemp==['+','-'] and int(pbam[8])<0:
                                                                    if not pos1 in LinkFR.keys():
                                                                            LinkFR[pos1]=[pos2]
                                                                    else:
                                                                            LinkFR[pos1]+=[pos2]
                                                                    if not pos1 in abInfF.keys():
                                                                            abInfF[pos1]=[0,1,0,0]
                                                                    else:
                                                                            abInfF[pos1][1]+=1
                                                                    if not pos2 in abInfR.keys():
                                                                       abInfR[pos2]=[0,1,0,0]
                                                                    else:
                                                                            abInfR[pos2][1]+=1
                                                                    if signIL>0:
                                                                            abInfF[pos1][0]+=1
                                                                            abInfR[pos2][0]+=1
                                                            elif DRtemp==['+','+'] and not int(pbam[8])==0:
                                                                    if not pos1 in abInfF.keys():
                                                                            abInfF[pos1]=[0,1,0,0]
                                                                    else:
                                                                            abInfF[pos1][1]+=1
                                                                    if signIL>0:
                                                                            abInfF[pos1][0]+=1
                                                                    if not pbam[0] in abtLink.keys():
                                                                            abtLink[pbam[0]]=[pos1]
                                                                    elif pbam[0] in abtLink.keys():
                                                                            abtLink[pbam[0]].append(pos1)
                                                                            if not min(abtLink[pbam[0]]) in LinkFF.keys():
                                                                                    LinkFF[min(abtLink[pbam[0]])]=[max(abtLink[pbam[0]])]
                                                                            else:
                                                                                    LinkFF[min(abtLink[pbam[0]])]+=[max(abtLink[pbam[0]])]
                                                                            del abtLink[pbam[0]]
                                                            elif DRtemp==['-','-'] and int(pbam[8])>0:
                                                                    if not min(pos1,pos2) in LinkRR.keys():
                                                                            LinkRR[min(pos1,pos2)]=[max(pos1,pos2)]
                                                                    else:
                                                                            LinkRR[min(pos1,pos2)]+=[max(pos1,pos2)]
                                                                    if not pos1 in abInfR.keys():
                                                                            abInfR[pos1]=[0,1,0,0]
                                                                    else:
                                                                            abInfR[pos1][1]+=1                                    
                                                                    if not pos2 in abInfR.keys():
                                                                            abInfR[pos2]=[0,1,0,0]
                                                                    else:
                                                                            abInfR[pos2][1]+=1
                                                                    if signIL>0:
                                                                            abInfR[pos1][0]+=1
                                                                            abInfR[pos2][0]+=1
                                                            elif int(pbam[8])==0:
                                                                    if int(pbam[1])&8>0:
                                                                            if DRtemp[0]=='+':
                                                                                    LinkNM['F'].append(pos1)
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abInfF.keys():
                                                                                                    abInfF[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfF[pos1][2]+=1
                                                                            elif DRtemp[0]=='-':
                                                                                    LinkNM['R'].append(pos1)
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abInfR.keys():
                                                                                                    abInfR[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfR[pos1][2]+=1
                                                                    if not pbam[6]=='=':
                                                                            if DRtemp[0]=='+':
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abLink.keys():
                                                                                                    abLink[pos1]=['f',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                            else:
                                                                                                    abLink[pos1]+=['f',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                            if not pos1 in abInfF.keys():
                                                                                                    abInfF[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfF[pos1][2]+=1
                                                                            elif DRtemp[0]=='-':
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abLink.keys():
                                                                                                    abLink[pos1]=['r',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                            else:
                                                                                                    abLink[pos1]+=['r',int(pbam[7]),pbam[6]+'_'+DRtemp[1]]
                                                                                            if not pos1 in abInfR.keys():
                                                                                                    abInfR[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfR[pos1][2]+=1
                                                        for k1 in LinkSPTemp.keys():
                                                            if not 'S' in LinkSPTemp[k1]:
                                                                del LinkSPTemp[k1]
                                                            else:
                                                                if len(LinkSPTemp[k1])==2:
                                                                    del LinkSPTemp[k1]
                                                        for k1 in LinkSPTemp.keys():
                                                            tempk1=[]
                                                            for k2 in LinkSPTemp[k1]:
                                                                if not k2=='S':
                                                                    tempk1+=[k2]
                                                            for k2 in range(len(tempk1[0])/2):
                                                                for k3 in range(len(tempk1[1])/2):
                                                                    if tempk1[0][2*k2]<tempk1[1][2*k3]:
                                                                        tempk2=[tempk1[0][2*k2],tempk1[0][2*k2+1],tempk1[1][2*k3],tempk1[1][2*k3+1]]
                                                                        if [tempk2[1],tempk2[3]]==['+','-']:
                                                                            if not tempk2[0] in LinkFR.keys():
                                                                                LinkFR[tempk2[0]]=[]
                                                                            LinkFR[tempk2[0]].append(tempk2[2])
                                                                        if [tempk2[1],tempk2[3]]==['+','+']:
                                                                            if not tempk2[0] in LinkFF.keys():
                                                                                LinkFF[tempk2[0]]=[]
                                                                            LinkFF[tempk2[0]].append(tempk2[2])
                                                                        if [tempk2[1],tempk2[3]]==['-','-']:
                                                                            if not tempk2[0] in LinkRR.keys():
                                                                                LinkRR[tempk2[0]]=[]
                                                                            LinkRR[tempk2[0]].append(tempk2[2])
                                                                        if [tempk2[1],tempk2[3]]==['-','+']:
                                                                            if not tempk2[0] in LinkRF.keys():
                                                                                LinkRF[tempk2[0]]=[]
                                                                            LinkRF[tempk2[0]].append(tempk2[2])
                                                        for k1 in abInfF.keys():
                                                                if not abInfF[k1][-1]==0:
                                                                        if not k1 in LinkSP.keys():
                                                                                LinkSP[k1]=0
                                                                        LinkSP[k1]+=abInfF[k1][-1]
                                                        for k1 in abInfR.keys():
                                                                if not abInfR[k1][-1]==0:
                                                                        if not k1 in LinkSP.keys():
                                                                                LinkSP[k1]=0
                                                                        LinkSP[k1]+=abInfR[k1][-1]                            
                                                        for k1 in LinkFR.keys():
                                                                for k2 in LinkFR[k1]:
                                                                        if not k2 in LinkRF.keys():
                                                                                LinkRF[k2]=[]
                                                                        LinkRF[k2].append(k1)
                                                        out_pair_bp=[]
                                                        out_single_bp=[]
                                                        SP4S=[]
                                                        if not LinkSP=={}:
                                                            SP2S=[]
                                                            linkSP_First=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'f')
                                                            for k1 in range(len(linkSP_First[0])):
                                                                    if linkSP_First[1][k1]==1:
                                                                            if linkSP_First[0][k1] in abInfF.keys():
                                                                                    if abInfF[linkSP_First[0][k1]]==[0,0,0,1]:
                                                                                            del abInfF[linkSP_First[0][k1]]
                                                                            if linkSP_First[0][k1] in abInfR.keys():
                                                                                    if abInfR[linkSP_First[0][k1]]==[0,0,0,1]:
                                                                                            del abInfR[linkSP_First[0][k1]]
                                                            linkSPF=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'f')[0]
                                                            linkSPR=clusterNums(sorted(LinkSP.keys()), SPCluLen, 'r')[0]
                                                            linkSPFR=clusterSupVis2(sorted(LinkSP.keys()),sorted(linkSPR),sorted(linkSPF),'left')
                                                            for k1 in linkSPFR.keys():
                                                                    qual_num=0
                                                                    qual_rec=[]
                                                                    for k2 in linkSPFR[k1]:
                                                                            qual_num+=LinkSP[k2]
                                                                            qual_rec.append(LinkSP[k2])
                                                                    if qual_num<SPCluMin:
                                                                            del linkSPFR[k1]
                                                                    else:
                                                                            SP2S.append(linkSPFR[k1][qual_rec.index(max(qual_rec))])
                                                            if not SP2S==[]:
                                                                    SP3S=[]
                                                                    key=[]
                                                                    if align_QCflag==1:
                                                                        for key1 in SP2S:
                                                                            QCLNSP_flag=0
                                                                            for aqb in [key1]:
                                                                                if aqb>BPAlignQCFlank:
                                                                                    LNSPFR_aqb=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,aqb-BPAlignQCFlank,aqb+BPAlignQCFlank))
                                                                                    LNSPFR_score=float(LNSPFR_aqb.read().strip())  
                                                                                    if LNSPFR_score<float(BPAlignQC):
                                                                                            QCLNSP_flag+=1
                                                                            if not QCLNSP_flag==0:
                                                                                    SP3S.append(key1)
                                                                    SP4S=[]
                                                                    for k1 in SP2S:
                                                                            if not k1 in SP3S:
                                                                                    SP4S.append(k1)
                                                                    SP4S.sort()
                                                                    for i in range(len(SP4S)-1):
                                                                            if SP4S[i+1]-SP4S[i]<10:
                                                                                    SP4S[i]=SP4S[i+1]
                                                                    SP5S=[]
                                                                    for i in SP4S:
                                                                            if not i in SP5S:
                                                                                    SP5S.append(i)
                                                                    SP4S=SP5S
                                                                    if not SP4S==[]:
                                                                            LinkSPF2=clusterSupVis2(sorted(abInfF.keys()), [i-ClusterLen for i in sorted(SP4S)], [i+10 for i in sorted(SP4S)],'right')
                                                                            for k1x in LinkSPF2.keys():
                                                                                key1=k1x-10
                                                                                LinkSPF2[key1]=[i for i in LinkSPF2[k1x]]
                                                                                del LinkSPF2[k1x]
                                                                            for k1x in LinkSPF2.keys():
                                                                                test1=bp_subgroup(LinkSPF2[k1x],Min_Distinguish_Len)
                                                                                if len(test1)>1:
                                                                                    for k2x in test1[:-1]:
                                                                                        new_core_ele=0
                                                                                        for k3x in k2x:
                                                                                            if k3x in abInfF.keys():
                                                                                                new_core_ele+=numpy.sum(abInfF[k3x])
                                                                                        if new_core_ele>BPLNCff:
                                                                                            if not max(k2x) in LinkSPF2.keys():
                                                                                                LinkSPF2[max(k2x)]=k2x
                                                                                            else:
                                                                                                LinkSPF2[max(k2x)]+=k2x
                                                                                            if not max(k2x) in SP4S:
                                                                                                SP4S.append(max(k2x))
                                                                            LinkSPR2=clusterSupVis2(sorted(abInfR.keys()), [i-10 for i in sorted(SP4S)], [i+ClusterLen for i in sorted(SP4S)],'left')
                                                                            for k1x in LinkSPR2.keys():
                                                                                key1=k1x+10
                                                                                LinkSPR2[key1]=[i for i in LinkSPR2[k1x]]
                                                                                del LinkSPR2[k1x]
                                                                            for k1x in LinkSPR2.keys():
                                                                                test1=bp_subgroup(LinkSPR2[k1x],Min_Distinguish_Len)
                                                                                if len(test1)>1:
                                                                                    for k2x in test1[1:]:
                                                                                        new_core_ele=0
                                                                                        for k3x in k2x:
                                                                                            if k3x in abInfR.keys():
                                                                                                new_core_ele+=numpy.sum(abInfR[k3x])
                                                                                        if new_core_ele>BPLNCff:
                                                                                            if not min(k2x) in LinkSPR2.keys():
                                                                                                LinkSPR2[min(k2x)]=k2x
                                                                                            else:
                                                                                                LinkSPR2[min(k2x)]+=k2x
                                                                                            if not min(k2x) in SP4S:
                                                                                                SP4S.append(min(k2x))
                                                                            for k1x in LinkSPF2.keys():
                                                                                temp_rec=LinkSPF2[k1x]
                                                                                LinkSPF2[k1x]=[LinkSPF2[k1x],[]]
                                                                                for k2 in temp_rec:
                                                                                    if k2 in LinkFR.keys():
                                                                                        LinkSPF2[k1x][1]+=LinkFR[k2]
                                                                                    if k2 in LinkFF.keys():
                                                                                        LinkSPF2[k1x][1]+=LinkFF[k2]
                                                                            for k1x in LinkSPR2.keys():
                                                                                temp_rec=LinkSPR2[k1x]
                                                                                LinkSPR2[k1x]=[LinkSPR2[k1x],[]]
                                                                                for k2 in temp_rec:
                                                                                    if k2 in LinkRR.keys():
                                                                                        LinkSPR2[k1x][1]+=LinkRR[k2]
                                                                                    if k2 in LinkRF.keys():
                                                                                        LinkSPR2[k1x][1]+=LinkRF[k2]
                                                                            LinkSP_To_Link={}
                                                                            for k1x in SP4S:
                                                                                if k1x in LinkSPF2.keys():
                                                                                    LinkSP_To_Link[k1x]=[[],[]]
                                                                                    LinkSP_To_Link[k1x][0]+=LinkSPF2[k1x][0]
                                                                                    LinkSP_To_Link[k1x][1]+=LinkSPF2[k1x][1]
                                                                                if k1x in LinkSPR2.keys():
                                                                                    if not k1x in LinkSP_To_Link.keys():
                                                                                        LinkSP_To_Link[k1x]=[[],[]]
                                                                                    LinkSP_To_Link[k1x][0]+=LinkSPR2[k1x][0]
                                                                                    LinkSP_To_Link[k1x][1]+=LinkSPR2[k1x][1]
                                                                                if k1x in LinkSP_To_Link.keys():
                                                                                    LinkSP_To_Link[k1x][0].sort()
                                                                                    LinkSP_To_Link[k1x][1].sort()
                                                                            for k1x in sorted(LinkSP_To_Link.keys()):
                                                                                for k2 in LinkSP_To_Link[k1x][1]:
                                                                                    if k2 in LinkSP_To_Link.keys():
                                                                                        if not sorted([k1x,k2]) in out_pair_bp and not k1x==k2:
                                                                                            out_pair_bp.append(sorted([k1x,k2]))
                                                                                        elif k1x==k2:
                                                                                            if not k1x in out_single_bp:
                                                                                                out_single_bp.append(k1x)
                                                                            for k1x in sorted(LinkSP_To_Link.keys()):
                                                                                if not LinkSP_To_Link[k1x][1]==[]:
                                                                                    for k2 in sorted(LinkSP_To_Link.keys())[sorted(LinkSP_To_Link.keys()).index(k1x):]:
                                                                                        if not LinkSP_To_Link[k2][1]==[]:
                                                                                            if k2==k1x: 
                                                                                                if not k1x in out_single_bp:
                                                                                                    out_single_bp.append(k1x)
                                                                                            else:
                                                                                                overlap_rec=[overlap_calcu(LinkSP_To_Link[k1x][1],LinkSP_To_Link[k2][0]),overlap_calcu(LinkSP_To_Link[k1x][0],LinkSP_To_Link[k2][1])]
                                                                                                if overlap_rec[0]+overlap_rec[1]>0:
                                                                                                    if not sorted([k1x,k2]) in out_pair_bp and not k1x==k2:
                                                                                                        out_pair_bp.append(sorted([k1x,k2]))
                                                                                                    elif k1x==k2:
                                                                                                        out_single_bp.append(k1x)
                                                                            for k1x in sorted(LinkSP_To_Link.keys()):
                                                                                if LinkSP_To_Link[k1x][1]==[]:
                                                                                    out_single_bp.append(k1x)
                                                                                    del LinkSP_To_Link[k1x]
                                                                            for k1x in out_pair_bp:
                                                                                for k2 in k1x:
                                                                                    if k2 in out_single_bp:
                                                                                        del out_single_bp[out_single_bp.index(k2)]
                                                                                    if k2 in LinkSP_To_Link.keys():
                                                                                        del LinkSP_To_Link[k2]
                                                                            for k1x in out_single_bp:
                                                                                if k1x in LinkSP_To_Link.keys():
                                                                                    del LinkSP_To_Link[k1x]
                                                                    SP4S.sort()
                                                        LinkSP_To_Link={}
                                                        tempIL={}
                                                        LinkIL={}
                                                        if not abInfF=={}:
                                                            clu_a_F=clusterNums(abInfF.keys(), ClusterLen, 'f')[0]
                                                            if not clu_a_F==[]:
                                                                clu_b_F=clusterNums(abInfF.keys(), ClusterLen, 'r')[0]
                                                                clu_c_F=clusterSupVis2(sorted(abInfF.keys()), clu_b_F, [caf+10 for caf in clu_a_F],'right')
                                                                tempIL_Info_Add(clu_c_F,LinkCluMin)
                                                        if not abInfR=={}:
                                                            clu_a_R=clusterNums(abInfR.keys(), ClusterLen, 'r')[0]
                                                            if not clu_a_R==[]:
                                                                clu_b_R=clusterNums(abInfR.keys(), ClusterLen, 'f')[0]
                                                                clu_c_R=clusterSupVis2(sorted(abInfR.keys()), [car-10 for car in clu_a_R], clu_b_R,'left')
                                                                tempIL_Info_Add(clu_c_R,LinkCluMin)
                                                        tempIL_Filter_One()
                                                        for k1 in tempIL.keys():
                                                            temp_mate_F={}
                                                            temp_mate_R={}
                                                            temp_mate_NewChr={}
                                                            info_mate=0
                                                            if tempIL[k1][0][2]=='f':
                                                                for k2 in tempIL[k1][1]:
                                                                    if k2 in LinkFF.keys():
                                                                        for k3 in LinkFF[k2]:
                                                                            if k3 in abInfF.keys():
                                                                                temp_mate_F[k3]=sum(abInfF[k3])
                                                                    if k2 in LinkFR.keys():
                                                                        for k3 in LinkFR[k2]:
                                                                            if k3 in abInfR.keys():
                                                                                temp_mate_R[k3]=sum(abInfR[k3])
                                                                    if k2 in abLink.keys():
                                                                        if not abLink[k2][2] in temp_mate_NewChr.keys():
                                                                            temp_mate_NewChr[abLink[k2][2]]={}
                                                                        if not abLink[k2][0] in temp_mate_NewChr[abLink[k2][2]].keys():
                                                                            temp_mate_NewChr[abLink[k2][2]][abLink[k2][0]]=[]
                                                                        temp_mate_NewChr[abLink[k2][2]][abLink[k2][0]].append(abLink[k2][1])
                                                            elif tempIL[k1][0][2]=='r':
                                                                for k2 in tempIL[k1][1]:
                                                                    if k2 in LinkRF.keys():
                                                                        for k3 in LinkRF[k2]:
                                                                            if k3 in abInfF.keys():
                                                                                temp_mate_F[k3]=sum(abInfF[k3])
                                                                    if k2 in LinkRR.keys():
                                                                        for k3 in LinkRR[k2]:
                                                                            if k3 in abInfR.keys():
                                                                                temp_mate_R[k3]=sum(abInfR[k3])
                                                                    if k2 in abLink.keys():
                                                                        if not abLink[k2][2] in temp_mate_NewChr.keys():
                                                                            temp_mate_NewChr[abLink[k2][2]]={}
                                                                        if not abLink[k2][0] in temp_mate_NewChr[abLink[k2][2]].keys():
                                                                            temp_mate_NewChr[abLink[k2][2]][abLink[k2][0]]=[]
                                                                        temp_mate_NewChr[abLink[k2][2]][abLink[k2][0]].append(abLink[k2][1])
                                                            if clusterQC(clusterNums4(temp_mate_F, ClusterLen, 'f'),LinkCluMin)==[] and clusterQC(clusterNums4(temp_mate_R, ClusterLen, 'r'),LinkCluMin)==[] and len(temp_mate_NewChr.keys())>5:
                                                                del tempIL[k1]
                                                            else:
                                                                for k1x in temp_mate_F.keys():
                                                                    info_mate+=temp_mate_F[k1x]
                                                                for k1x in temp_mate_R.keys():
                                                                    info_mate+=temp_mate_R[k1x]
                                                                for k1x in temp_mate_NewChr.keys():
                                                                    for k2x in temp_mate_NewChr[k1x].keys():
                                                                        info_mate+=len(temp_mate_NewChr[k1x][k2x])
                                                                if not info_mate<LinkCluMin:
                                                                    LinkIL[k1]=[[],[]]
                                                                    if not temp_mate_F=={}:
                                                                        LinkIL[k1][0]=clusterQC(clusterNums4(temp_mate_F, ClusterLen, 'f'),LinkCluMin)
                                                                    if not temp_mate_R=={}:
                                                                        LinkIL[k1][1]=clusterQC(clusterNums4(temp_mate_R, ClusterLen, 'r'),LinkCluMin)
                                                                    for xa in temp_mate_NewChr.keys():
                                                                        for xb in temp_mate_NewChr[xa].keys():
                                                                            LinkILTemp=clusterQC(clusterNums(temp_mate_NewChr[xa][xb],ClusterLen,'f'),LinkCluMin)
                                                                            if not LinkILTemp==[]:
                                                                                LinkIL[k1].append([xa,xb]+LinkILTemp)
                                                        for k1 in LinkIL.keys():
                                                                for k2 in LinkIL[k1]:
                                                                        for k3 in k2:
                                                                            if k3>BPAlignQCFlank:
                                                                                if align_QCflag==1:
                                                                                    tPairF_QC=0
                                                                                    tPairF_a=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,k3-BPAlignQCFlank,k3+BPAlignQCFlank))
                                                                                    tPairF_b=float(tPairF_a.read().strip())  
                                                                                    if not tPairF_b<float(BPAlignQC):
                                                                                            if not [min([k3,k1]),max([k3,k1])] in out_pair_bp and not k1==k3:
                                                                                                out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                                                                            elif k1==k3:
                                                                                                out_single_bp.append(k1)
                                                                                else:
                                                                                            if not [min([k3,k1]),max([k3,k1])] in out_pair_bp and not k1==k3:
                                                                                                    out_pair_bp.append([min([k3,k1]),max([k3,k1])])        
                                                                                            elif  k1==k3:     
                                                                                                out_single_bp.append(k1)                                  
                                                        if not out_pair_bp==[]:
                                                            out_pair_bp_2=[]
                                                            for x in out_pair_bp:
                                                                out_pair_bp_flag=0
                                                                for y in x:
                                                                    if not type(y)==type(1):
                                                                        out_pair_bp_flag+=1
                                                                if out_pair_bp_flag==0:
                                                                    out_pair_bp_2.append(x)
                                                            out_pair_bp=out_pair_bp_2
                                                            temp_out_pair_bp=[]
                                                            out_BPmodify={}
                                                            for k1 in out_pair_bp:
                                                                for k2 in k1:
                                                                    if type(k2)==type(1):
                                                                        out_BPmodify[k2]=[]
                                                            if not SP4S==[]:
                                                                LBSP_tempIL=clusterSupVis3(sorted(SP4S),sorted(out_BPmodify.keys()))
                                                                for k1 in out_pair_bp:
                                                                    temp_k1=[]
                                                                    for k2 in k1:
                                                                        if type(k2)==type(1):
                                                                            k3=LBSP_tempIL[k2]
                                                                            if abs(k2-k3)<ClusterLen:
                                                                                temp_k1.append(k3)
                                                                            else:
                                                                                temp_k1.append(k2)
                                                                    if len(temp_k1)>1:
                                                                        temp_out_pair_bp.append(temp_k1)
                                                            if not temp_out_pair_bp==[]:
                                                                out_pair_bp=temp_out_pair_bp
                                                            for k1 in out_pair_bp:
                                                                for k2 in k1:
                                                                    if k2 in SP4S:
                                                                        del SP4S[SP4S.index(k2)]
                                                            out_pair_modify={}
                                                            for i in out_pair_bp:
                                                                if not i[0] in out_pair_modify.keys():
                                                                    out_pair_modify[i[0]]=[]
                                                                if not i[1] in out_pair_modify[i[0]]:
                                                                    out_pair_modify[i[0]].append(i[1])
                                                            if len(out_pair_modify)>1:
                                                                while True: 
                                                                    if len(out_pair_modify)==1: break
                                                                    out_pair_qc=[]
                                                                    for i in range(len(sorted(out_pair_modify.keys()))-1):
                                                                        out_pair_qc.append(sorted(out_pair_modify.keys())[i+1]-sorted(out_pair_modify.keys())[i])
                                                                    if min(out_pair_qc)>50:break
                                                                    else:
                                                                        out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))+1]]+=out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                                        del out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                            for k1 in out_pair_modify.keys():
                                                                while True:
                                                                    if len(out_pair_modify[k1])==1: break
                                                                    out_pair_modify[k1].sort()
                                                                    out_pair_qc=[]
                                                                    for i in range(len(out_pair_modify[k1])-1):
                                                                        out_pair_qc.append(out_pair_modify[k1][i+1]-out_pair_modify[k1][i])
                                                                    if min(out_pair_qc)>50:break
                                                                    else:
                                                                        out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))+1]=out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))]
                                                                    for i in out_pair_modify[k1]:
                                                                        if out_pair_modify[k1].count(i)>1:
                                                                            for j in range(out_pair_modify[k1].count(i)-1):
                                                                                del out_pair_modify[k1][out_pair_modify[k1].index(i)]
                                                            out_pair_numrec={}
                                                            for k1 in out_pair_modify.keys():
                                                                out_pair_numrec[k1]=[]
                                                                for k2 in [k1]+out_pair_modify[k1]:
                                                                    if k2 in tempIL.keys() and not k2 in LinkSP.keys():
                                                                        out_pair_numrec[k1].append(len(tempIL[k2][1]))
                                                                    elif k2 in LinkSP.keys() and not k2 in tempIL.keys():
                                                                        out_pair_numrec[k1].append(LinkSP[k2])
                                                                    elif k2 in LinkSP.keys() and k2 in tempIL.keys():
                                                                        out_pair_numrec[k1].append(LinkSP[k2]+len(tempIL[k2][1]))
                                                                    else:
                                                                        out_pair_numrec[k1].append(0)
                                                            fout=open(BPOutputb,'a')
                                                            for i in sorted(out_pair_modify.keys()):
                                                                for j in out_pair_modify[i]:                            
                                                                    print >>fout, ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                                    if i in tempIL.keys():
                                                                        del tempIL[i]
                                                                    if j in tempIL.keys():
                                                                        del tempIL[j]
                                                            fout.close()
                                                            fout=open(BPOutputa,'a')
                                                            for i in sorted(out_single_bp+LinkSP_To_Link.keys()):
                                                                num=0
                                                                if i in abInfF.keys():
                                                                    num+=sum(abInfF[i])
                                                                if i in abInfR.keys():
                                                                    num+=sum(abInfR[i])
                                                                print >>fout, ' '.join([str(j) for j in [i,num]])
                                                            fout.close()
                                                        for i in tempIL.keys():
                                                            temp_IL_Rec[i]=tempIL[i]
                                                            temp_mate_F=[]
                                                            temp_mate_R=[]
                                                            for k2 in tempIL[i][1]:
                                                                if k2 in LinkFF.keys():
                                                                        temp_mate_F+=LinkFF[k2]
                                                                if k2 in LinkFR.keys():
                                                                        temp_mate_R+=LinkFR[k2]
                                                            temp_IL_Rec[i].append(temp_mate_F) 
                                                            temp_IL_Rec[i].append(temp_mate_R) 
                                                        tempLNF=[]
                                                        tempLNR=[]
                                                        for key in abLink.keys():
                                                                for key2 in range(len(abLink[key])/3):
                                                                        if abLink[key][3*key2]=='f':
                                                                                tempLNF.append(key)
                                                                        else:
                                                                                tempLNR.append(key)
                                                        for key in abLinkSP.keys():
                                                                for key2 in range(len(abLinkSP[key])/3):
                                                                        if abLinkSP[key][3*key2]=='f':
                                                                                tempLNF.append(key)
                                                                        else:
                                                                                tempLNR.append(key)
                                                        FtLNFR=clusterNums(tempLNF+tempLNR,ClusterLen,'f')[0]
                                                        RtLNFR=clusterNums(tempLNF+tempLNR,ClusterLen,'r')[0]
                                                        FRtLNFR=clusterSupVis2(sorted(tempLNF+tempLNR),RtLNFR,FtLNFR,'left')
                                                        for key1 in FRtLNFR.keys():
                                                                t_LNFR=[]
                                                                for key2 in FRtLNFR[key1]:
                                                                        if key2 in abLink.keys():
                                                                                t_LNFR+=abLink[key2]
                                                                        else:
                                                                                if abs(key2-key1)<SPCluLen or abs(key2-max(FRtLNFR[key1]))<SPCluLen:
                                                                                        t_LNFR+=abLinkSP[key2]
                                                                if len(t_LNFR)/3<LinkCluMin:
                                                                        del FRtLNFR[key1]
                                                                else:
                                                                    if not t_LNFR in FRtLNFR[key1]:
                                                                        FRtLNFR[key1].append(t_LNFR)
                                                        FRtLNFRb=[]
                                                        for key1 in FRtLNFR.keys():
                                                            t1_LN=[]
                                                            t2_LN=[] 
                                                            t3_LN={}
                                                            t4out=[]
                                                            for key2 in FRtLNFR[key1][:-1]:
                                                                    if key2 in abLink.keys():
                                                                        if not key2 in t1_LN:
                                                                            for key3 in range(len(abLink[key2])/3):
                                                                                if not abLink[key2][3*key3:3*(key3+1)] in t2_LN:
                                                                                    t2_LN.append(abLink[key2][3*key3:3*(key3+1)])
                                                                                    t1_LN.append(key2)
                                                            for key2 in t2_LN:
                                                                    if not key2[-1].split('_')[0] in t3_LN.keys():
                                                                            t3_LN[key2[-1].split('_')[0]]={}
                                                                            t3_LN[key2[-1].split('_')[0]]['a']=[]
                                                                            t3_LN[key2[-1].split('_')[0]]['b']=[]
                                                                            t3_LN[key2[-1].split('_')[0]]['c']=[]
                                                                            t3_LN[key2[-1].split('_')[0]]['d']=[]
                                                                    t3_LN[key2[-1].split('_')[0]]['a'].append(key2[0])
                                                                    t3_LN[key2[-1].split('_')[0]]['b'].append(key2[1])
                                                                    t3_LN[key2[-1].split('_')[0]]['c'].append(key2[-1].split('_')[1])
                                                                    t3_LN[key2[-1].split('_')[0]]['d'].append(key2)                    
                                                            for key2 in t3_LN.keys():
                                                                if clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'f'), LinkCluMin)==[]:
                                                                        del t3_LN[key2]
                                                                else:
                                                                        t4LN=clusterSupVis2(t3_LN[key2]['b'],clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'r'), LinkCluMin),clusterQC(clusterNums(t3_LN[key2]['b'], ClusterLen, 'f'), LinkCluMin), 'left')                               
                                                                        for key5 in t4LN.keys():
                                                                                t4LNa=[]
                                                                                t4LNb=[]
                                                                                t4LNc=[]
                                                                                t4LNd=[]
                                                                                t4out=[]
                                                                                for key6 in t4LN[key5]:
                                                                                        t4LNa.append(t3_LN[key2]['a'][t3_LN[key2]['b'].index(key6)])
                                                                                        t4LNc.append(t3_LN[key2]['c'][t3_LN[key2]['b'].index(key6)])
                                                                                        t4LNd.append(key6)
                                                                                        t4LNb.append(t1_LN[t2_LN.index(t3_LN[key2]['d'][t3_LN[key2]['b'].index(key6)])])
                                                                                if not 'f' in t4LNa or float(t4LNa.count('r'))/float(t4LNa.count('f'))>5:
                                                                                        t4out+=[chrom,'r',min(t4LNb)]
                                                                                elif not 'r' in t4LNa or float(t4LNa.count('f'))/float(t4LNa.count('r'))>5:
                                                                                        t4out+=[chrom,'f',max(t4LNb)]
                                                                                else:
                                                                                        t4out+=[chrom,min(t4LNb),max(t4LNb)]
                                                                                if not '+' in t4LNc or float(t4LNc.count('-'))/float(t4LNc.count('+'))>5:
                                                                                        t4out+=[key2,'-',min(t4LNd)]
                                                                                elif not '-' in t4LNc or float(t4LNc.count('+'))/float(t4LNc.count('-'))>5:
                                                                                        t4out+=[key2,'+',max(t4LNd)]
                                                                                else:
                                                                                        t4out+=[key2,min(t4LNd),max(t4LNd)]
                                                            if not t4out==[] and not t4out in FRtLNFRb:
                                                                    FRtLNFRb.append(t4out)             
                                                        if not FRtLNFRb==[]:
                                                                fout=open(BPOutputd,'a')
                                                                for keyfrt in FRtLNFRb:
                                                                        print >>fout, ' '.join([str(keyfr2) for keyfr2 in keyfrt])
                                                                fout.close()
                                            Link_IL_Rec={}
                                            for k1 in temp_IL_Rec.keys():
                                                temp_mate_F=temp_IL_Rec[k1][2]
                                                temp_mate_R=temp_IL_Rec[k1][3]
                                                if not len(temp_IL_Rec[k1][1])+len(temp_IL_Rec[k1][2])+len(temp_IL_Rec[k1][3])<LnCluMin:
                                                    Link_IL_Rec[k1]=[clusterQC(clusterNums(temp_mate_F, ClusterLen, 'f'),LinkCluMin),
                                                                    clusterQC(clusterNums(temp_mate_R, ClusterLen, 'r'),LinkCluMin)]
                                                    del temp_IL_Rec[k1]
                                                else:continue
                                            for k1 in Link_IL_Rec.keys():
                                                for k2 in Link_IL_Rec[k1]:
                                                        for k3 in k2:
                                                            if k3>BPAlignQCFlank:
                                                                if align_QCflag==1:
                                                                    tPairF_QC=0
                                                                    tPairF_a=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom,k3-BPAlignQCFlank,k3+BPAlignQCFlank))
                                                                    tPairF_b=float(tPairF_a.read().strip())  
                                                                    if not tPairF_b<float(BPAlignQC):
                                                                            if not [min([k3,k1]),max([k3,k1])] in out_pair_bp:
                                                                                    out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                                                else:
                                                                            if not [min([k3,k1]),max([k3,k1])] in out_pair_bp:
                                                                                    out_pair_bp.append([min([k3,k1]),max([k3,k1])])
                                            if not out_pair_bp==[]:
                                                temp_out_pair_bp=[]
                                                out_BPmodify={}
                                                for k1 in out_pair_bp:
                                                    for k2 in k1:
                                                        out_BPmodify[k2]=[]
                                                if not SP4S==[]:
                                                    LBSP_tempIL=clusterSupVis3(sorted(SP4S),sorted(out_BPmodify.keys()))
                                                    for k1 in out_pair_bp:
                                                        temp_k1=[]
                                                        for k2 in k1:
                                                            k3=LBSP_tempIL[k2]
                                                            if abs(k2-k3)<ClusterLen:
                                                                temp_k1.append(k3)
                                                            else:
                                                                temp_k1.append(k2)
                                                        temp_out_pair_bp.append(temp_k1)
                                                out_pair_bp=temp_out_pair_bp
                                                for k1 in out_pair_bp:
                                                    for k2 in k1:
                                                        if k2 in SP4S:
                                                            del SP4S[SP4S.index(k2)]
                                                out_pair_modify={}
                                                for i in out_pair_bp:
                                                    if not i[0] in out_pair_modify.keys():
                                                        out_pair_modify[i[0]]=[]
                                                    if not i[1] in out_pair_modify[i[0]]:
                                                        out_pair_modify[i[0]].append(i[1])
                                                if len(out_pair_modify)>1:
                                                    while True: 
                                                        if len(out_pair_modify)==1: break
                                                        out_pair_qc=[]
                                                        for i in range(len(sorted(out_pair_modify.keys()))-1):
                                                            out_pair_qc.append(sorted(out_pair_modify.keys())[i+1]-sorted(out_pair_modify.keys())[i])
                                                        if min(out_pair_qc)>50:break
                                                        else:
                                                            out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))+1]]+=out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                            del out_pair_modify[sorted(out_pair_modify.keys())[out_pair_qc.index(min(out_pair_qc))]]
                                                for k1 in out_pair_modify.keys():
                                                    while True:
                                                        if len(out_pair_modify[k1])==1: break
                                                        out_pair_modify[k1].sort()
                                                        out_pair_qc=[]
                                                        for i in range(len(out_pair_modify[k1])-1):
                                                            out_pair_qc.append(out_pair_modify[k1][i+1]-out_pair_modify[k1][i])
                                                        if min(out_pair_qc)>50:break
                                                        else:
                                                            out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))+1]=out_pair_modify[k1][out_pair_qc.index(min(out_pair_qc))]
                                                        for i in out_pair_modify[k1]:
                                                            if out_pair_modify[k1].count(i)>1:
                                                                for j in range(out_pair_modify[k1].count(i)-1):
                                                                    del out_pair_modify[k1][out_pair_modify[k1].index(i)]
                                                out_pair_numrec={}
                                                for k1 in out_pair_modify.keys():
                                                    out_pair_numrec[k1]=[]
                                                    for k2 in [k1]+out_pair_modify[k1]:
                                                        if k2 in tempIL.keys() and not k2 in LinkSP.keys():
                                                            out_pair_numrec[k1].append(len(tempIL[k2][1]))
                                                        elif k2 in LinkSP.keys() and not k2 in tempIL.keys():
                                                            out_pair_numrec[k1].append(LinkSP[k2])
                                                        elif k2 in LinkSP.keys() and k2 in tempIL.keys():
                                                            out_pair_numrec[k1].append(LinkSP[k2]+len(tempIL[k2][1]))
                                                        else:
                                                            out_pair_numrec[k1].append(0)
                                                fout=open(BPOutputb,'a')
                                                for i in sorted(out_pair_modify.keys()):
                                                    for j in out_pair_modify[i]:                            
                                                        print >>fout, ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                        if i in tempIL.keys():
                                                            del tempIL[i]
                                                        if j in tempIL.keys():
                                                            del tempIL[j]
                                                fout.close()
                                            fout=open(BPOutputa,'a')
                                            for i in sorted(temp_IL_Rec.keys()):
                                                print >>fout, ' '.join([str(i),str(sum(temp_IL_Rec[i][0][:2]))])    
                                            fout.close()
                                            time2=time.time()
                                            os.system(r'''cat %s >> %s'''%(BPOutputd,BPOutpute))
                                            os.system(r'''rm %s'''%(BPOutputd))
                                time2=time.time()
                                print 'BPSearch Complete for '+bamF
                                print 'Time Consuming: '+str(time2-time1)
    if function_name=='BPIntegrate':
        def clusterSupVis2(dataRound, dataLeft, dataRight):
            flag1=0
            flag2=0
            out={}
            while True:
                if flag1==len(dataRound): break
                else:
                    if flag2==len(dataLeft)-1 and dataRound[flag1]>dataRight[flag2]: break
                    else:
                        if dataRound[flag1]<dataLeft[flag2]:
                            flag1+=1
                        elif dataRound[flag1]>dataRight[flag2]:
                            flag2+=1
                        else:
                            if not dataLeft[flag2] in out.keys():
                                out[dataLeft[flag2]]=[dataRound[flag1]]
                            else:
                                out[dataLeft[flag2]]+=[dataRound[flag1]]
                            flag1+=1
            return out
        def SearchClosest(dataRound, hashcentre):
            data1=sorted(dataRound)
            data2=sorted(hashcentre.keys())
            for data in data2:
                if data1[0]-data>0:
                    disa=data1[0]-data
                    hashcentre[data]=[data1[0],disa]
                else:
                    for flag2 in range(len(data1)-1):
                        disa=data1[flag2]-data
                        disb=data1[flag2+1]-data
                        if not disa*disb>0:
                            if numpy.abs(disa)<numpy.abs(disb):
                                hashcentre[data]=[data1[flag2],disa]
                            else:
                                hashcentre[data]=[data1[flag2+1],disb]
        def SearchClosest2(DataRou1, DataCen1,cutoff):
            dataRou=sorted(DataRou1)
            dataCen=sorted(DataCen1)
            out={}
            out2=[]
            for k1 in dataCen:
                for k2 in dataRou:
                    if k2>k1 and dataRou[dataRou.index(k2)-1]<k1:
                        k3=dataRou[dataRou.index(k2)-1]
                        if abs(k2-k1)>abs(k3-k1):
                            out[k1]=k3
                        else:
                            out[k1]=k2
                        if abs(out[k1]-k1)>cutoff:
                            out2.append(k1)
                            del out[k1]
                        break
            return [out,out2]
        def clusterNums(Data1, ClusterLen, direction):
            if not Data1==[]:
                data=sorted(Data1)
                out2=[]
                out3=[]
                out=[]
                for i in data:
                    if len(out)==0:
                        out.append(i)
                    else:
                        if i-out[-1]< ClusterLen:
                            out.append(i)
                        else:
                            if out[-1]-out[0]<ClusterLen:
                                if direction=='f':
                                    out2.append(out[-1])
                                    out3.append(len(out))
                                elif direction=='r':
                                    out2.append(out[0])
                                    out3.append(len(out))
                            else:
                                temp=[0]
                                lenTime=int((out[-1]-out[0])/ClusterLen)
                                lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                                lenInter2=sorted(lenInter)[::-1][:lenTime]
                                for j in range(len(lenInter)):
                                    if lenInter[j] in  lenInter2:
                                        temp.append(j+1)
                                temp.append(len(out))
                                if direction=='f':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][-1])
                                        out3.append(temp[k+1]-temp[k])
                                elif direction=='r':
                                    for k in range(len(temp)-1):
                                        out2.append(out[temp[k]:temp[k+1]][0])
                                        out3.append(temp[k+1]-temp[k])
                            out=[i]
                if out[-1]-out[0]<ClusterLen:
                    if direction=='f':
                        out2.append(out[-1])
                        out3.append(len(out))
                    elif direction=='r':
                        out2.append(out[0])
                        out3.append(len(out))
                else:
                    temp=[0]
                    lenTime=int((out[-1]-out[0])/ClusterLen)
                    lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                    lenInter2=sorted(lenInter)[::-1][:lenTime]
                    for j in range(len(lenInter)):
                        if lenInter[j] in  lenInter2:
                            temp.append(j+1)
                    temp.append(len(out))
                    if direction=='f':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][-1])
                            out3.append(temp[k+1]-temp[k])
                    elif direction=='r':
                        for k in range(len(temp)-1):
                            out2.append(out[temp[k]:temp[k+1]][0])
                            out3.append(temp[k+1]-temp[k])
                return [out2,out3]
            else:
                return [[],[]]
        def clusterQC(hash, CluNum):
            out=[]
            for i in range(len(hash[1])):
                if not hash[1][i]<CluNum:
                    out.append(hash[0][i])
            return out
        def clusterSupVis3(dataRound,dataCen):
            out={}
            if not dataCen==[] and not dataRound==[]:
                for key1 in dataCen:
                    min_num=abs(dataRound[0]-key1)
                    min_rec=dataRound[0]
                    for key2 in dataRound[1:]:
                        alt_num=abs(key2-key1)
                        alt_rec=key2
                        if alt_num<min_num:
                            min_num=alt_num
                            min_rec=alt_rec
                    out[key1]=min_rec
            return out
        def clusterSupVis(DataRou1,DataCen1,ClusterLen):
            dataRou=sorted(DataRou1)
            dataCen=sorted(DataCen1)
            flag1=0
            flag2=0
            out={}
            if not dataRou==[] and not dataCen==[]:
                while True:
                    if flag1==len(dataRou): break
                    else:
                        if flag2==len(dataCen)-1:
                            if dataRou[flag1]>dataCen[flag2]+ClusterLen: break
                            elif dataRou[flag1]<dataCen[flag2]-ClusterLen: 
                                flag1+=1
                            else:
                                if not dataCen[flag2] in out.keys():
                                    out[dataCen[flag2]]=[dataRou[flag1]]
                                else:
                                    out[dataCen[flag2]]+=[dataRou[flag1]]
                                flag1+=1
                        else:
                            if dataRou[flag1]<dataCen[flag2]-ClusterLen:
                                flag1+=1
                            elif dataRou[flag1]>dataCen[flag2]+ClusterLen:
                                flag2+=1
                            else:
                                if not dataCen[flag2] in out.keys():
                                    out[dataCen[flag2]]=[dataRou[flag1]]
                                else:
                                    out[dataCen[flag2]]+=[dataRou[flag1]]
                                flag1+=1
            return out
        def LNd_hash_readin(chrom_LN):
            LNd_hash={}
            fin=open(bps_in_path+chrom_LN)
            for line in fin:
                pin=line.strip().split()
                temp=sorted([pin[0],pin[3]])
                pin2=[temp[0],pin[pin.index(temp[0])+1],pin[pin.index(temp[0])+2],temp[1],pin[pin.index(temp[1])+1],pin[pin.index(temp[1])+2]]
                if not '_'.join([pin2[0],pin2[3]]) in LNd_hash.keys():
                    tempkey='_'.join([pin2[0],pin2[3]])
                    LNd_hash[tempkey]=[{},{}]
                if not pin2[1] in LNd_hash[tempkey][0].keys():
                    LNd_hash[tempkey][0][pin2[1]]=[]
                LNd_hash[tempkey][0][pin2[1]]+=pin2[4:6]
                if not pin2[2] in LNd_hash[tempkey][0].keys():
                    LNd_hash[tempkey][0][pin2[2]]=[]
                LNd_hash[tempkey][0][pin2[2]]+=pin2[4:6]
                if not pin2[4] in LNd_hash[tempkey][1].keys():
                    LNd_hash[tempkey][1][pin2[4]]=[]
                LNd_hash[tempkey][1][pin2[4]]+=pin2[1:3]
                if not pin2[5] in LNd_hash[tempkey][1].keys():
                    LNd_hash[tempkey][1][pin2[5]]=[]
                LNd_hash[tempkey][1][pin2[5]]+=pin2[1:3]
            fin.close()
            for k1 in LNd_hash.keys():
                for k2 in LNd_hash[k1][0].keys():
                    if k2 in ['f','r','+','-']:
                        del LNd_hash[k1][0][k2]
                    else:
                        temp=[]
                        for k3 in LNd_hash[k1][0][k2]:
                            if not k3 in ['f','r','+','-'] and not k3 in temp:
                                temp.append(k3)
                        LNd_hash[k1][0][k2]=temp
                for k2 in LNd_hash[k1][1].keys():
                    if k2 in ['f','r','+','-']:
                        del LNd_hash[k1][1][k2]
                    else:
                        temp=[]
                        for k3 in LNd_hash[k1][1][k2]:
                            if not k3 in ['f','r','+','-'] and not k3 in temp:
                                temp.append(k3)
                        LNd_hash[k1][1][k2]=temp
            tempall={}
            temp1=[]
            for k1 in LNd_hash.keys():
                temp1=[]
                for k2 in LNd_hash[k1][0].keys():
                    temp2=[k2]
                    for k3 in LNd_hash[k1][0][k2]:
                        for k4 in LNd_hash[k1][0].keys():
                            if not k4==k2 and k3 in LNd_hash[k1][0][k4]:
                                if not k4 in temp2:
                                    temp2.append(k4)
                    if not temp2==[k2]:
                        temp2.sort()
                        temp3=[]
                        for k3 in temp2:
                            for k4 in LNd_hash[k1][0][k3]:
                                if not k4 in temp3:
                                    temp3.append(k4)
                        temp3.sort()
                        if not [k1.split('_')[0]]+temp2+[k1.split('_')[1]]+temp3 in temp1:
                            temp1.append([k1.split('_')[0]]+temp2+[k1.split('_')[1]]+temp3)
                for k2 in LNd_hash[k1][1].keys():
                    temp2=[k2]
                    for k3 in LNd_hash[k1][1][k2]:
                        for k4 in LNd_hash[k1][1].keys():
                            if not k4==k2 and k3 in LNd_hash[k1][1][k4]:
                                if not k4 in temp2:
                                    temp2.append(k4)
                    if not temp2==[k2]:
                        temp2.sort()
                        temp3=[]
                        for k3 in temp2:
                            for k4 in LNd_hash[k1][1][k3]:
                                if not k4 in temp3:
                                    temp3.append(k4)
                        temp3.sort()
                        if not [k1.split('_')[1]]+temp2+[k1.split('_')[0]]+temp3 in temp1:
                            temp1.append([k1.split('_')[1]]+temp2+[k1.split('_')[0]]+temp3)
                tempall[k1]=temp1
                for k2 in tempall[k1]:
                    if k2[0]==k1.split('_')[0]:
                        for k3 in k2[1:]:
                            if k3 in k1.split('_'): break
                            if k3 in LNd_hash[k1][0].keys():
                                del LNd_hash[k1][0][k3]
                    if k2[0]==k1.split('_')[1]:
                        for k3 in k2[1:]:
                            if k3 in k1.split('_'): break
                            if k3 in LNd_hash[k1][1].keys():
                                del LNd_hash[k1][1][k3]
                temp2=[]
                clu1=clusterNums(sorted([int(i) for i in LNd_hash[k1][0].keys()]), 1000, 'f')
                cen1=clusterQC(clu1, 2)
                if not cen1==[]:
                    sup1=clusterSupVis(sorted([int(i) for i in LNd_hash[k1][0].keys()]),cen1,1000)
                    for k3 in sup1.keys():
                        temp3=[k1.split('_')[0]]+[str(i) for i in sorted(sup1[k3])]+[k1.split('_')[1]]
                        temp4=[]
                        for k4 in sup1[k3]:
                            temp4+=LNd_hash[k1][0][str(k4)]
                        temp4.sort()
                        temp3+=temp4
                        if not temp3 in temp2:
                            temp2.append(temp3)
                clu2=clusterNums(sorted([int(i) for i in LNd_hash[k1][1].keys()]), 1000, 'f')
                cen2=clusterQC(clu2, 2)
                if not cen2==[]:
                    sup2=clusterSupVis(sorted([int(i) for i in LNd_hash[k1][1].keys()]),cen2,1000)
                    for k3 in sup2.keys():
                        temp3=[k1.split('_')[1]]+[str(i) for i in sorted(sup2[k3])]+[k1.split('_')[0]]
                        temp4=[]
                        for k4 in sup2[k3]:
                            temp4+=LNd_hash[k1][1][str(k4)]
                        temp4.sort()
                        temp3+=temp4
                        if not temp3 in temp2:
                            temp2.append(temp3)
                tempall[k1]+=temp2
                for k2 in temp2:
                    if k2[0]==k1.split('_')[0]:
                        for k3 in k2[1:]:
                            if k3 in k1.split('_'): break
                            if k3 in LNd_hash[k1][0].keys():
                                del LNd_hash[k1][0][k3]
                    if k2[0]==k1.split('_')[1]:
                        for k3 in k2[1:]:
                            if k3 in k1.split('_'): break
                            if k3 in LNd_hash[k1][1].keys():
                                del LNd_hash[k1][1][k3]
            return [LNd_hash,tempall]
        def LNe_hash_Filter(LNe_hash):
            out={}
            for k1 in LNe_hash.keys():
                tempall=[]
                for k2 in LNe_hash[k1]:
                    khash={}
                    kchr=k2[0]
                    khash[kchr]=[]
                    for k3 in k2[1:]:
                        if k3 in allchromos or not type(k3) == type(1):
                            khash[k3]=[]
                            kchr=k3
                        elif type(k3) == type(1) and int(k3)<23:
                            khash[k3]=[]
                            kchr=k3                    
                        else:
                            if type(k3) == type(1) and int(k3)>22:
                                khash[kchr].append(int(k3))
                    temp=[]
                    if len(khash.keys())>1:
                        for k3 in sorted(khash.keys()):
                            temp+=[k3]+sorted(khash[k3])
                    if not temp in tempall and not temp==[]:
                        tempall.append(temp)
                if not tempall==[]:
                    out[k1]=tempall
            return out
        def LNd_hash_pair(LNd_hash):
            out={}
            for k1 in LNd_hash.keys():
                if k1.split('_')[0] in SP_Info.keys() and k1.split('_')[1] in SP_Info.keys():
                    out[k1]=[]
                    cen1=[int(i) for i in LNd_hash[k1][0].keys()]
                    if not cen1==[]:
                        clu1=clusterSupVis(sorted(SP_Info[k1.split('_')[0]].keys()),cen1,10000)
                        for k2 in cen1:
                            if k2 in clu1.keys():
                                temp=[]
                                for k3 in clu1[k2]:
                                    temp.append(abs(k3-k2))
                                clu1[k2]=clu1[k2][temp.index(min(temp))]
                            else:
                                clu1[k2]=k2+100
                    for k2 in LNd_hash[k1][0].keys():
                        k3=sorted([int(i) for i in LNd_hash[k1][0][k2]])
                        if len(k3)>1:
                            continue
                        else:
                            if not k3==[]:
                                clu2=clusterSupVis(sorted(SP_Info[k1.split('_')[0]].keys()),k3,10000)
                                for k4 in k3:
                                    if k4 in clu2.keys():
                                        temp=[]
                                        for k3 in clu2[k4]:
                                            temp.append(abs(k3-k4))
                                        clu2[k4]=clu2[k4][temp.index(min(temp))]
                                    else:
                                        clu2[k4]=k4+100
                                LNd_hash[k1][0][k2]=[clu2.keys()[0],clu2[clu2.keys()[0]]]
                    if not cen1==[]:
                        for k2 in clu1.keys():
                            temp=[k1.split('_')[0]]+sorted([k2,clu1[k2]])+[k1.split('_')[1]]+sorted(LNd_hash[k1][0][str(k2)])
                            if not temp in out[k1]:
                                out[k1].append(temp)
            return out
        def SP_Info_Cluster(SP_Info):
            for S_Chr in SP_Info.keys():
                bpCluster=clusterNums(sorted(SP_Info[S_Chr].keys()), 11, 'f')
                bpMilti=[]
                for k1 in range(len(bpCluster[1])):
                    if bpCluster[1][k1]>1:
                        bpMilti.append(bpCluster[0][k1])
                if not bpMilti==[]:
                    bpClu2=clusterSupVis(sorted(SP_Info[S_Chr].keys()),sorted(bpMilti),11)
                    for k1 in bpClu2.keys():
                        temp=[]
                        for k2 in bpClu2[k1]:
                            temp.append(SP_Info[S_Chr][k2])
                        newk=bpClu2[k1][temp.index(max(temp))]
                        for k2 in bpClu2[k1]:
                            if not k2==newk:
                                SP_Info[S_Chr][newk]+=SP_Info[S_Chr][k2]
                                del SP_Info[S_Chr][k2]
                else:
                    bpClu2=[]
            del bpCluster
            del bpClu2
        def SP_info_ReadIn(bps_hash,bps_in_path):
            single_chromos=bps_hash[S_Sample][S_Para].keys()
            SP_links={}
            SP_Link1={}
            SP_Info={}
            SP_Link2={}
            SP_Link3={}
            for S_Chr in single_chromos:
                SP_links[S_Chr]={}
                SP_Link1[S_Chr]=[]
                SP_Info[S_Chr]={}
                filein1=bps_hash[S_Sample][S_Para][S_Chr][0]
                fin=open(bps_in_path+filein1)
                for line in fin:
                    pin=line.strip().split()
                    if not int(pin[0]) in SP_Info[S_Chr].keys():
                        SP_Info[S_Chr][int(pin[0])]=int(pin[1])
                    else:
                        SP_Info[S_Chr][int(pin[0])]=max(int(pin[1]),SP_Info[S_Chr][int(pin[0])])
                    if not int(pin[2]) in SP_Info[S_Chr].keys():
                        SP_Info[S_Chr][int(pin[2])]=int(pin[3])
                    else:
                        SP_Info[S_Chr][int(pin[2])]=max(int(pin[3]),SP_Info[S_Chr][int(pin[2])])
                    if not pin[0]==pin[2]:
                        SP_Link1[S_Chr]+=[pin[0],pin[2]]
                fin.close()
                filein1=bps_hash[S_Sample][S_Para][S_Chr][1]
                fin=open(bps_in_path+filein1)
                for line in fin:
                        pin=line.strip().split()
                        if not int(pin[0]) in SP_Info[S_Chr].keys():
                            SP_Info[S_Chr][int(pin[0])]=int(pin[1])
                        else:
                            SP_Info[S_Chr][int(pin[0])]=max(int(pin[1]),SP_Info[S_Chr][int(pin[0])])
                fin.close()
            SP_Info_Cluster(SP_Info)
            for S_Chr in SP_Link1.keys():
                SP_Link2[S_Chr]=[int(i) for i in SP_Link1[S_Chr]]
            for S_Chr in single_chromos:
                bpCluster=clusterSupVis(SP_Info[S_Chr].keys(),SP_Link2[S_Chr],10)
                SP_Link3[S_Chr]=[]
                for k1 in bpCluster.keys():
                    if len(bpCluster[k1])==1 and bpCluster[k1][0]==k1:
                        del bpCluster[k1]
                for k1 in range(len(SP_Link2[S_Chr])):
                    k2=SP_Link2[S_Chr][k1]
                    if k2 in bpCluster.keys():
                        if len(bpCluster[k2])==1 and bpCluster[k2][0]==k2:
                            SP_Link3[S_Chr]+=[k2]
                        else:
                            if len(bpCluster[k2])==1:
                                SP_Link3[S_Chr]+=[bpCluster[k2][0]]
                            else:
                                temp=[]
                                for k3 in bpCluster[k2]:
                                    temp.append(SP_Info[S_Chr][k3])
                                k1_replace=bpCluster[k2][temp.index(max(temp))]
                                SP_Link3[S_Chr]+=[k1_replace]
                    else:
                        SP_Link3[S_Chr]+=[k2]
                SP_links[S_Chr]={}
                for k1 in range(len(SP_Link3[S_Chr])/2):
                    minNum=min([SP_Link3[S_Chr][2*k1],SP_Link3[S_Chr][2*k1+1]])
                    maxNum=max([SP_Link3[S_Chr][2*k1],SP_Link3[S_Chr][2*k1+1]])
                    if not minNum in SP_links[S_Chr].keys():
                        SP_links[S_Chr][minNum]=[]
                    SP_links[S_Chr][minNum]+=[maxNum]
            return [SP_links,SP_Info,SP_Link3]
        def SP_Info_Sort(SP_links):
            out={}
            for k1 in SP_links.keys():
                out[k1]=[]
                for k2 in sorted(SP_links[k1].keys()):
                    if not sorted([k2]+SP_links[k1][k2]) in out[k1]:
                        out[k1].append(sorted([k2]+SP_links[k1][k2]))
            return out
        def BP_Pair_Collaps(SP_list_Chr,SP_Info_Chr):
            out1={}
            out2=[]
            temp=[]
            out=[]
            Cff1=10**6
            for k1 in SP_list_Chr:
                k2=[]
                for k3 in k1:
                    if not k3 in k2:
                        k2.append(k3)
                k2.sort()
                if len(k2)>1:
                    if len(k2)>4:
                        out1[min(k2)]=sorted(k2)
                    else:
                        if k2[-1]-k2[0]<Cff1:
                            temp.append(k2)
                        else:
                            out2.append(k2)
            if not temp==[]:
                temp2=temp[0]
                flag1=0
                while True:
                    flag1+=1
                    if flag1==len(temp): break
                    if temp[flag1][0]-temp[flag1-1][-1]<500:
                        temp2+=temp[flag1]
                    else:
                        out.append(temp2)
                        temp2=temp[flag1]
            Cff2=5*(10**6)
            out3=[]
            if not out2==[]:
                for k1 in out2:
                    k1.sort()
                    k2=[[k1[0]]]
                    for k3 in k1[1:]:
                        if k3-k2[-1][-1]<Cff2:
                            k2[-1].append(k3)
                        else:
                            k2.append([k3])
                    k5=[]
                    for k3 in k2:
                        if len(k3)>1: 
                            k5+=k3
                        else: 
                            k4=SearchClosest2(SP_Info_Chr.keys(),k3,Cff2)[0]
                            if not k4=={}:
                                k5+=sorted([k4.keys()[0],k4[k4.keys()[0]]])
                            else:
                                k5+=k3
                    out3.append(k5)
            if not out+out3==[]:
                for k1 in out+out3:
                    if not min(k1) in out1.keys():
                        out1[min(k1)]=[]
                    out1[min(k1)]+=k1
            reorder=[]
            for k1 in sorted(out1.keys()):
                reorder.append(out1[k1])
            return reorder
        def Left_SP_Cluster(Left_SP_1,Cff):
            Left_SP=sorted(Left_SP_1)
            out=[[Left_SP[0]]]
            for k1 in Left_SP[1:]:
                if k1-out[-1][-1]<Cff:
                    out[-1].append(k1)
                else:
                    out.append([k1])
            out1=[]
            out2=[]
            for k1 in out:
                if len(k1)>1:
                    out1.append(k1)
                else:
                    out2+=k1
            return [out1,out2]
        def pair_link(SP_LI_in):
            singletons={}
            SP_links=SP_LI_in[0]
            SP_Info=SP_LI_in[1]
            SP_Link1=SP_LI_in[2]
            SP_Single={}
            SP_link2={}
            SP_link4={}
            for k1 in SP_Info.keys():
                SP_Single[k1]=[]
                for k2 in SP_Info[k1].keys():
                    if SP_Info[k1][k2]<2:
                        del SP_Info[k1][k2]
                    else:
                        if not k2 in SP_Link1[k1]:
                            SP_Single[k1].append(k2)
                SP_assign=SearchClosest2(SP_Link1[k1], SP_Single[k1],5000)
                Left_SP=SP_assign[1]
                Paired_SP=SP_assign[0]
                if not Left_SP==[]:
                    SingletonPair=Left_SP_Cluster(Left_SP,5000)
                    singletons[k1]=SingletonPair[1]
                    for k2 in Paired_SP.keys():
                        if SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2] in SP_links[k1].keys():
                            tempk=SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2]
                        elif SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2+1] in SP_links[k1].keys():
                            tempk=SP_Link1[k1][SP_Link1[k1].index(Paired_SP[k2])/2*2+1]
                        else: continue
                        SP_links[k1][tempk].append(k2)
                    for k2 in SingletonPair[0]:
                        if not min(k2) in SP_links[k1].keys():
                            SP_links[k1][min(k2)]=[]
                        SP_links[k1][min(k2)]+=k2
            SP_list=SP_Info_Sort(SP_links)
            collapsed_SP_list={}
            for k1 in SP_list.keys():
                collapsed_SP_list[k1]=BP_Pair_Collaps(SP_list[k1],SP_Info[k1])   
            return [collapsed_SP_list,singletons]
        def single_split(singleton_chr,out_single):
            if len(singleton_chr)<5:
                out_single.append(singleton_chr)
            else:
                single_inter=[]
                for kx in range(len(singleton_chr[1:])):
                    single_inter.append(singleton_chr[kx+1]-singleton_chr[kx])
                max_index=1
                max_rec=single_inter[1]
                for kx in range(2,len(single_inter)-1):
                    if single_inter[kx]>max_rec:
                        max_rec=single_inter[kx]
                        max_index=kx
                temp=[singleton_chr[:(max_index+1)],singleton_chr[(max_index+1):]]
                for temp_list in temp:
                    single_split(temp_list,out_single)
        def singleton_split(singleton_chr):
            out_single=[]
            SingletonChr=sorted(singleton_chr)
            single_split(SingletonChr,out_single)
            return out_single
        def write_bp_1(SP_link3,missed_pairs,rec):
            fout_rec=[]
            max_num_rec=[[]]
            max_num=int(dict_opts['--batch'])
            for k1 in SP_link3.keys():
                for k2 in SP_link3[k1]:
                    k2.sort()
                    k3=[k2[0]]
                    for k4 in k2:
                        if not k4 in k3 and k4-k3[-1]>10:
                            k3.append(k4)
                    k5=k3
                    if len(k5)>1 and not len(k5)>8:
                        if len(max_num_rec[-1])<max_num:
                            max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])
                        elif len(k5)==1:
                            max_num_rec.append([])
                            max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])                                                
                    else:
                        singletons[k1]+=k5
                for k5 in missed_pairs[k1]:
                    if len(max_num_rec[-1])<max_num:
                            max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])
                    else:
                        max_num_rec.append([])
                        max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k5])                                                
            for k1 in singletons.keys():      
                if len(singletons[k1])>1:
                    out_single=singleton_split(singletons[k1])
                    for k2 in out_single:
                        if len(max_num_rec[-1])<max_num:
                            max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k2])
                        else:
                            max_num_rec.append([])
                            max_num_rec[-1].append([str(k1)]+[str(k3) for k3 in k2])
            for k1 in max_num_rec:
                rec+=1
                fout=para_sample_bps_fold+S_Sample+'.'+str(rec)+'.'+'txt'
                if not fout in fout_rec:
                    fout_rec.append(fout)
                if not os.path.isfile(fout):
                    fo=open(fout,'w')
                else:
                    fo=open(fout,'a')
                for k2 in k1:
                    print >>fo, ' '.join(k2)
                    print >>fo, ' '
                fo.close()
            for fout in fout_rec:
                fin=open(fout)
                test=[]
                for line in fin:
                    pin=line.strip().split()
                    if not pin==[]:
                        if not pin in test:
                            test.append(pin)
                fin.close()
                fo=open(fout,'w')
                for k1 in test:
                    print >>fo, ' '.join(k1)
                    print >>fo, ' '
                fo.close()
            return rec
        def bp_lists_length_decide(bps):
            flag=0
            for x in bps:
                if len(x)>8:
                    flag+=1
            if flag==0:
                return 'TRUE'
            else:
                return 'FALSE'
        def bp_list_separate(bps):
            while True:
                if bp_lists_length_decide(bps)=='TRUE':
                    break
                elif bp_lists_length_decide(bps)=='FALSE':
                    for bps2 in bps:
                        if len(bps2)>8:
                            inter=[bps2[x+1]-bps2[x] for x in range(len(bps2)-1)]
                            inter1=max(inter[1:-1])
                            bps+=[bps2[:inter.index(inter1)+1],bps2[inter.index(inter1)+1:]]
                            del bps[bps.index(bps2)]
        def bps_list_chop(SP_link4):
            for x in SP_link4.keys():
                for y in SP_link4[x]:
                    if len(y)>8:
                        bps_y=[y]
                        bp_list_separate(bps_y)
                        SP_link4[x]+=bps_y
        def write_bp_2(SP_link3,missed_pairs):
            fout_rec=[]
            for k1 in SP_link3.keys():
                fout=para_sample_bps_fold+S_Sample+'.'+str(k1)+'.'+'txt'
                if not fout in fout_rec:
                    fout_rec.append(fout)
                if not os.path.isfile(fout):
                    fo=open(fout,'w')
                else:
                    fo=open(fout,'a')
                for k2 in SP_link3[k1]:
                    k2.sort()
                    k3=[k2[0]]
                    for k4 in k2:
                        if not k4 in k3 and k4-k3[-1]>10:
                            k3.append(k4)
                    k5=k3
                    if len(k5)>1 and not len(k5)>8:
                        print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k5])
                        print >>fo, ' '
                    elif len(k5)==1:
                        if k1 in singletons.keys():
                            singletons[k1]+=k5
                        else:
                            singletons[k1]=k5
                for k2 in missed_pairs[k1]:
                    k2.sort()
                    print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                    print >>fo, ' '
                fo.close()
            for k1 in singletons.keys():
                fout=para_sample_bps_fold+S_Sample+'.'+str(k1)+'.'+'txt'
                if not fout in fout_rec:
                    fout_rec.append(fout)
                if not os.path.isfile(fout):
                    fo=open(fout,'w')
                else:
                    fo=open(fout,'a')
                if len(singletons[k1])>1:
                    out_single=singleton_split(singletons[k1])
                    for k2 in out_single:
                        print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                        print >>fo, ' '
                fo.close()
            for fout in fout_rec:
                fin=open(fout)
                test=[]
                for line in fin:
                    pin=line.strip().split()
                    if not pin==[]:
                        if not pin in test:
                            test.append(pin)
                fin.close()
                fo=open(fout,'w')
                for k1 in test:
                    print >>fo, ' '.join(k1)
                    print >>fo, ' '
                fo.close()
        def write_bp_3(SP_link3,missed_pairs):
            fout_rec=[]
            fout=para_sample_bps_fold+S_Sample+'.txt'
            if not fout in fout_rec:
                fout_rec.append(fout)        
            if not os.path.isfile(fout):
                fo=open(fout,'w')
            else:
                fo=open(fout,'a')
            fo.close()
            for k1 in SP_link3.keys():
                fo=open(fout,'a')
                for k2 in SP_link3[k1]:
                    k2.sort()
                    k3=[k2[0]]
                    for k4 in k2:
                        if not k4 in k3 and k4-k3[-1]>10:
                            k3.append(k4)
                    k5=k3
                    if len(k5)>1 and not len(k5)>8:
                        print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k5])
                        print >>fo, ' '
                    elif len(k5)==1:
                        if k1 in singletons.keys():
                            singletons[k1]+=k5
                        else:
                            singletons[k1]=k5
                for k2 in missed_pairs[k1]:
                    k2.sort()
                    print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                    print >>fo, ' '
                fo.close()
            for k1 in singletons.keys():
                fo=open(fout,'a')
                if len(singletons[k1])>1:
                    out_single=singleton_split(singletons[k1])
                    for k2 in out_single:
                        print >>fo, ' '.join([str(k1)]+[str(k3) for k3 in k2])
                        print >>fo, ' '
                fo.close()
            for fout in fout_rec:
                fin=open(fout)
                test=[]
                for line in fin:
                    pin=line.strip().split()
                    if not pin==[]:
                        if not pin in test:
                            test.append(pin)
                fin.close()
                fo=open(fout,'w')
                for k1 in test:
                    print >>fo, ' '.join(k1)
                    print >>fo, ' '
                fo.close()
        def missed_pair_check(SP_Check,SP_link3):
            out={}
            out2={}
            for k1 in SP_Check.keys():
                rec1=0
                rec2=0
                temp={}
                temp1=SP_Check[k1]
                temp1.sort()
                temp2=SP_link3[k1]
                for x in temp2:
                    if not x[0] in temp.keys():
                        temp[x[0]]=[]
                    temp[x[0]].append(x)
                temp3=[]
                for x in sorted(temp.keys()):
                    temp3+=temp[x]
                temp4={}
                while True:
                    if rec1==len(temp1) or rec2==len(temp3): break
                    if temp3[rec2][-1]-temp3[rec2][0]<10**7:
                        if temp1[rec1]<temp3[rec2][0]-10**4:
                            rec1+=1
                        elif temp1[rec1]>temp3[rec2][1]+10**4:
                            rec2+=1
                        else:
                            if not '_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])]) in temp4.keys():
                                temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])]=[]
                            if not temp3[rec2] in temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])]:
                                temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])].append(temp3[rec2])
                            if not temp1[rec1] in temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])][0]:
                                temp4['_'.join([str(temp3[rec2][0]),str(temp3[rec2][-1])])].append(temp1[rec1])
                            rec1+=1
                    else:
                        rec2+=1
                out[k1]=[]
                temp5=[]
                for x in temp4.keys():
                    if len(temp4[x])==1:
                        out[k1].append(sorted(temp4[x][0]))
                        temp5+=sorted(temp4[x][0])
                    else:
                        out[k1].append(sorted(temp4[x][0]+temp4[x][1:]))
                        temp5+=sorted(temp4[x][0]+temp4[x][1:])
                temp6=[]
                for x in temp1:
                    if not x in temp5:
                        temp6.append(x)
                temp6.sort()
                temp6b=[]
                for x in temp3:
                    flagb=0
                    for y in x:
                        if not y in temp5:
                            flagb+=1
                    if not flagb==0:
                        temp6b.append(x)
                out[k1]+=temp6b
                data_cen=clusterNums(temp6, 10**6, 'r')[0]
                out2[k1]=[]
                if not data_cen==[]:
                    temp7=[[data_cen[0]]]
                    rec1=0
                    for x in temp6:
                        if x<temp7[-1][0]:
                            temp7[-1].append(x)
                        else:
                            rec1+=1
                            if rec1<len(data_cen):
                                temp7.append([data_cen[rec1]])  
                    for x in temp7:
                        if len(x)>1:
                            out2[k1].append(x)
            return [out,out2]
        min_length=100
        def link_SP_link3(SP_link3):
            out={}
            for k1 in SP_link3.keys():
                temp2=[]
                temp1={}
                for k2 in SP_link3[k1]:
                    if not k2[0] in temp1.keys():
                        temp1[k2[0]]=[]
                    temp1[k2[0]].append(k2)
                for k2 in sorted(temp1.keys()):
                    for k3 in temp1[k2]:
                        temp2.append(k3)
                out[k1]=[]
                for k2 in temp2:
                    if out[k1]==[]:
                        out[k1].append(k2)
                    else:
                        if k2[0]<out[k1][-1][-1]+10 and out[k1][-1][-1]-out[k1][-1][0]<2*10**6:
                            out[k1][-1]+=k2
                        else:
                            out[k1].append(k2)
            out2={}
            for k1 in out.keys():
                out2[k1]=[]
                for k2 in out[k1]:
                    k2.sort()
                    temp=[k2[0]]
                    for k3 in k2[1:]:
                        if not k3 in temp and k3-temp[-1]>100:
                            temp.append(k3)
                    if len(temp)>1:
                        out2[k1].append(temp)
            for k1 in out2.keys():
                temp=[]
                for k3 in out2[k1]:
                    temp2=[k3[0]]
                    for k4 in k3[1:]:
                        if k4-temp2[-1]>min_length:
                            temp2.append(k4)
                        else:
                            if SP_Info[k1][k4]>SP_Info[k1][temp2[-1]]:
                                temp2[-1]=k4
                            else:
                                continue
                    if len(temp2)>1:
                        temp.append(temp2)
                out2[k1]=temp
            return out2
        def Define_Default_BPIntegrate():
            global ReadLength
            if not '--read-length' in dict_opts.keys():
                ReadLength=101
            else:
                ReadLength=int(dict_opts['--read-length'])
            global ToolMappingQ
            global FileMappingQ
            global align_QCflag
            if '--qc-map-tool' in dict_opts.keys() and '--qc-map-file' in dict_opts.keys():
                ToolMappingQ=dict_opts['--qc-map-tool']
                FileMappingQ=dict_opts['--qc-map-file']
                align_QCflag=1
            else:
                align_QCflag=0
            global BPalignQC
            if '--BPalignQC' in dict_opts.keys():
                BPalignQC=float(dict_opts['--BPalignQC'])
            else:
                BPalignQC=0.2
            global QCAlign
            if '--qc-align' in dict_opts.keys():
                QCAlign=int(dict_opts['--qc-align'])
            else:
                QCAlign=20
            global QCSplit
            if '--qc-split' in dict_opts.keys():
                QCSplit=int(dict_opts['--qc-split'])
            else:
                QCSplit=20
            global Null_SplitLen_perc
            if '--split-min-len' in dict_opts.keys():
                Null_SplitLen_perc=float(dict_opts['--split-min-len'])
            else:
                Null_SplitLen_perc=0.1
            global BPalignQCFlank
            if '--BPalignQCFlank' in dict_opts.keys():
                BPalignQCFlank=int(dict_opts['--BPalignQCFlank'])
            else:
                BPalignQCFlank=500
            global para_filter
            para_filter=[]
            if '--BPSPCff' in dict_opts.keys() and '--BPLNCff' in dict_opts.keys() and '--BPalignQC' in dict_opts.keys():
                para_filter=['SPCff'+dict_opts['--BPSPCff']+'.CluCff'+dict_opts['--BPLNCff']+'.AlignCff'+dict_opts['--BPalignQC']]
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def chromos_read_in(ref_index):
            chromos=[]
            fin=open(ref_index)
            for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
            fin.close()
            return chromos
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print ' '
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--sample, input alignment file in bam format'
            print ' '
            print 'Optional Parameters:'
            print '--chromosome, name of chromosome to run. should match chromosome name in bam file'
            print '--batch, specify number of structures in each separate file (if 0, output files will be calssified by chromosomes; default, all BP clustered will be integrated in one txt file)'
        else:
            Define_Default_BPIntegrate()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                out_path=workdir
                print 'temp files produced under: '+workdir
                bps_in_path=workdir+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
                if not '--sample' in dict_opts.keys():
                    print 'Error: please specify either input file using --sample'
                else:
                    if '--sample' in dict_opts.keys():
                        bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                        bam_files=[dict_opts['--sample']]
                        bam_names=[dict_opts['--sample'].split('/')[-1].replace('.bam','')]
                    else:
                        bam_path=dict_opts['--PathBam']
                        if not bam_path[-1]=='/':
                            bam_path+='/'
                        bam_files=[]
                        bam_names=[]
                        for file in os.listdir(bam_path):
                            if file.split('.')[-1]=='bam':
                                bam_files.append(bam_path+file)
                                bam_names.append(file.replace('.bam',''))
                    ref_path=workdir+'reference/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        chromos=chromos_read_in(ref_index)
                        allchromos=chromos
                        if '--chromosome' in dict_opts.keys():
                            chrom_single=dict_opts['--chromosome']
                            if not chrom_single in chromos:
                                print 'Error: please make sure the chromosome defined by --chr is correct based on the reference genome'
                                chromos=[]
                            else:
                                chromos=[chrom_single]
                        if not chromos==[]:
                            time1=time.time()
                            bps_hash={}
                            for i in bam_names:
                                bps_hash[i]={}
                            bps_folder=out_path+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'
                            if not os.path.isdir(bps_folder):
                                os.system(r'''mkdir %s'''%(bps_folder))
                            for i in bam_names:
                                bps_hash[i]={}
                                for file1 in os.listdir(bps_in_path):
                                    if file1.split('.')[-1]=='LNs':
                                        if i in file1: 
                                            keyj='.'.join(file1.split('.')[-5:-1])
                                            #single_chr_name=file1.split('.')[-6]
                                            single_chr_name='.'.join(file1.replace('.'.join(dict_opts['--sample'].split('/')[-1].split('.')[:-1])+'.','').split('.')[:-5])
                                            if not keyj in bps_hash[i].keys():
                                                bps_hash[i][keyj]={}
                                            if single_chr_name in chromos:
                                                bps_hash[i][keyj][single_chr_name]=[]
                                                bps_hash[i][keyj][single_chr_name].append(file1)
                                                bps_hash[i][keyj][single_chr_name].append(file1.replace('LNs','SPs'))
                            for i in bps_hash.keys():
                                if not os.path.isdir(bps_folder+i):
                                    os.system(r'''mkdir %s'''%(bps_folder+i))
                                for j in bps_hash[i].keys():
                                    if not os.path.isdir(bps_folder+i+'/'+j):
                                        os.system(r'''mkdir %s'''%(bps_folder+i+'/'+j))
                            for S_Sample in bps_hash.keys():
                                if para_filter==[]:
                                    para_filter=sorted(bps_hash[S_Sample].keys())[::-1]
                                for S_Para in para_filter:
                                    rec=0
                                    para_sample_bps_fold=bps_folder+S_Sample+'/'+S_Para+'/'
                                    SP_LI_in=SP_info_ReadIn(bps_hash,bps_in_path)
                                    SP_links=SP_LI_in[0]
                                    SP_Info=SP_LI_in[1]
                                    SP_Check=SP_LI_in[2]
                                    SP_link3a=pair_link(SP_LI_in)
                                    singletons=SP_link3a[1]
                                    SP_link3=SP_link3a[0]
                                    missed_pairs=missed_pair_check(SP_Check,SP_link3)
                                    SP_link4=link_SP_link3(missed_pairs[0])
                                    bps_list_chop(SP_link4)
                                    missed_pairs=link_SP_link3(missed_pairs[1])
                                    devide_start=10**6
                                    if not '--batch' in dict_opts.keys():
                                        write_bp_3(SP_link4,missed_pairs)
                                    else:
                                        if dict_opts['--batch']=='0':
                                            write_bp_2(SP_link4,missed_pairs)
                                        else:
                                            rec=write_bp_1(SP_link4,missed_pairs,rec)
                                    S_Chr=bps_hash[S_Sample][S_Para].keys()[0] 
                                    chrom_LN=bps_hash[S_Sample][S_Para][S_Chr][1].replace('.'+S_Chr+'.SPCff','.SPCff').replace('SPs','chromLNs')
                                    LNall_hash=LNd_hash_readin(chrom_LN)
                                    LNd_hash=LNall_hash[0]
                                    LNe_hash=LNe_hash_Filter(LNall_hash[1])
                                    LNd_hash=LNd_hash_pair(LNd_hash)
                                    LN_out=[]
                                    for k1 in LNe_hash.keys():
                                        for k2 in LNe_hash[k1]:
                                            temp=[]
                                            for k3 in k2:
                                                if k3 in chromos:
                                                    temp.append([k3])
                                                else:
                                                    if not temp==[]:
                                                        temp[-1].append(str(k3))
                                            LN_out.append(temp)
                                    for k1 in LNd_hash.keys():
                                        for k2 in LNd_hash[k1]:
                                            temp=[]
                                            for k3 in k2:
                                                if k3 in chromos:
                                                    temp.append([k3])
                                                else:
                                                    if not temp==[]:
                                                        temp[-1].append(str(k3))
                                            LN_out.append(temp)
                                    if not '--batch' in dict_opts.keys():
                                        fout=para_sample_bps_fold+S_Sample+'.txt'
                                        if not os.path.isfile(fout):
                                            fo=open(fout,'w')
                                        else:
                                            fo=open(fout,'a')
                                        for k1 in LN_out:
                                            for k2 in k1:
                                                print >>fo, ' '.join(k2)
                                            print >>fo, ' '
                                        fo.close()
                                    else:
                                        if dict_opts['--batch']=='0':
                                            fout=para_sample_bps_fold+S_Sample+'.'+'LN'+'.'+'txt'
                                            if not os.path.isfile(fout):
                                                fo=open(fout,'w')
                                            else:
                                                fo=open(fout,'a')
                                            for k1 in LN_out:
                                                for k2 in k1:
                                                    print >>fo, ' '.join(k2)
                                                print >>fo, ' '
                                            fo.close()
                                        else:
                                            max_num=int(dict_opts['--batch'])
                                            LN_out2=[[]]
                                            for k1 in LN_out:
                                                if len(LN_out2[-1])<max_num:
                                                    LN_out2[-1].append(k1)
                                                else:
                                                    LN_out2.append([])
                                                    LN_out2[-1].append(k1)
                                            for k1 in LN_out2:
                                                rec+=1
                                                fout=para_sample_bps_fold+S_Sample+'.'+str(rec)+'.'+'txt'
                                                if not os.path.isfile(fout):
                                                    fo=open(fout,'w')
                                                else:
                                                    fo=open(fout,'a')
                                                for k2 in k1:
                                                    for k3 in k2:
                                                        print >>fo, ' '.join(k3)
                                                    print >>fo, ' '
                                                fo.close()
                            time2=time.time()
                            print 'BPIntegrate Complete !'
                            print 'Time Consuming: '+str(time2-time1)
    if function_name=='SVPredict':
        def pdf_calculate(x,alpha,mean1,mean2,std1,std2,upper_limit,lower_limit,Penalty_For_InsertLengthZero):
            Alpha=numpy.min([alpha,1-alpha])
            if mean1<mean2:
                Mean1=mean1
                Std1=std1
                Mean2=mean2
                Std2=std2
            elif mean1>mean2: 
                Mean1=mean2
                Std1=std2
                Mean2=mean1
                Std2=std1
            if not Alpha==0:
                if x<upper_limit and x>lower_limit:
                    return math.log(Alpha/math.sqrt(2*math.pi*math.pow(Std1,2))*math.exp(-math.pow((x-Mean1),2)/(2*math.pow(Std1,2)))+(1-Alpha)/math.sqrt(2*math.pi*math.pow(Std2,2))*math.exp(-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))))
                elif x>=upper_limit:
                    test1=math.pow((x-Mean1),2)/(2*math.pow(Std1,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
                    if test1<200:
                        return math.log(Alpha)-0.5*math.log(2*math.pi*math.pow(Std1,2))-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))+math.log(1+(1-Alpha)*Std1/(Alpha*Std2)*math.exp(test1))
                    elif test1>200:
                        return math.log(1-Alpha)-0.5*math.log(2*math.pi*math.pow(Std2,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
                elif x<=lower_limit and not x==0:
                    test2=-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))+math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
                    if test2<200:
                        return math.log(1-Alpha)-0.5*math.log(2*math.pi*math.pow(Std2,2))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))+math.log(1+Alpha*Std2/((1-Alpha)*Std1)*math.exp(test2))
                    if test2>200:
                        return  math.log(Alpha)-0.5*math.log(2*math.pi*math.pow(Std1,2))-math.pow((x-Mean1),2)/(2*math.pow(Std1,2))
                elif x==0:
                    return Penalty_For_InsertLengthZero
            else:
                if not x==0:
                    return math.log(1-Alpha)-math.log(math.sqrt(2*math.pi*math.pow(Std2,2)))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
                elif x==0:
                    return Penalty_For_InsertLengthZero
        def standard_pdf_IL_calculate(x,alpha,mean1,mean2,std1,std2,upper_limit,lower_limit,Penalty_For_InsertLengthZero):
            if alpha==1:
                return standard_norm_pdf_solver(x,mean1,std1)
            elif alpha==0:
                return standard_norm_pdf_solver(x,mean2,std2)
            else:
                return 'Error'
        def bimodal_cdf_solver(cdf,alpha,mean1,mean2,std1,std2):
            if cdf>0.5:
                x=numpy.min([mean1,mean2])
                while True:
                    x+=0.1
                    fc=alpha*(norm.cdf((x-mean1)/std1))+(1-alpha)*(norm.cdf((x-mean2)/std2))
                    if abs(fc-cdf)<0.001:
                        return x
            elif cdf<0.5:
                x=numpy.max([mean1,mean2])
                while True:
                    x-=0.1
                    fc=alpha*(norm.cdf((x-mean1)/std1))+(1-alpha)*(norm.cdf((x-mean2)/std2))
                    if abs(fc-cdf)<0.001:
                        return x
        def norm_cdf_solver(cdf,mean,std):
            x=mean
            if cdf>0.5:
                while True:
                    x+=0.1
                    fc=(norm.cdf((x-mean)/std))
                    if abs(fc-cdf)<0.001:
                        return x
            elif cdf<0.5:
                while True:
                    x-=0.1
                    fc=(norm.cdf((x-mean)/std))
                    if abs(fc-cdf)<0.001:
                        return x
        def norm_pdf_solver(x,mean,std):
            return -0.5*math.log(2*math.pi*math.pow(std,2))-math.pow((x-mean),2)/(2*math.pow(std,2))
        def standard_norm_pdf_solver(x,mean,std):
            x_new=(x-mean)/std
            return -0.5*math.log(2*math.pi)-math.pow((x_new),2)/(2)
        def Prob_Possion(Number, Mean):
            if not Mean==0:
                probs=[]
                from scipy.stats import poisson
                lamda=Mean
                K=Number
                K2=Float_2_Integer(K)
                pdf=K2*math.log(lamda)-lamda-math.log(math.factorial(K2))
                return pdf
            elif Mean==0:
                return -Number
        def Prob_NB(Number, Mean, Variance):
            if not Mean==0:
                P=1-Mean/Variance
                R=Float_2_Integer(Mean*(1-P)/P)
                K2=Float_2_Integer(Number)
                if P>0 and R>1 and K2+R-1>0 :
                    log_pdf=math.log(math.factorial(K2+R-1))-math.log(math.factorial(K2))-math.log(math.factorial(R-1))+R*math.log(1-P)+K2*math.log(P)
                else:
                    log_pdf=Prob_Possion(Number, Mean)
                return log_pdf
            elif Mean==0:
                return -Number
        def Prob_Norm(Number, Mean, Variance):
            return math.log(1/sqrt(2*numpy.pi*Variance)*exp(-(Number-Mean)**2/2/Variance))
        def cdf_solver_application(Insert_Len_Stat,cdf):
            fstat=open(Insert_Len_Stat)
            temp=fstat.readline()
            temp=fstat.readline()
            temp=fstat.readline()
            if model_comp=='S':
                data1=fstat.readline().strip().split()
                fstat.close()
                flank_out=int(round(norm_cdf_solver(float(cdf),float(data1[1]),float(data1[2]))))
            elif model_comp=='C':
                data1=fstat.readline().strip().split()
                fstat.readline()
                data2=fstat.readline().strip().split()
                flank_out=int(round(bimodal_cdf_solver(float(cdf),float(data1[0]),float(data1[1]),float(data2[1]),float(data1[2]),float(data2[2]))))
                fstat.close()
            return flank_out
        def GC_Stat_ReadIn(BamN,GC_Stat_Path,affix):
            GC_Stat_File=GC_Stat_Path+'/'+BamN+'.'+genome_name+affix
            f_GC_stat=open(GC_Stat_File)
            CN2_Region={}
            Chromos=f_GC_stat.readline().strip().split()
            GC_Content=f_GC_stat.readline().strip().split()
            for key_1 in Chromos:
                    CN2_Region[key_1]={}
                    for key_2 in GC_Content:
                            CN2_Region[key_1][int(key_2)]=f_GC_stat.readline().strip().split()
            f_GC_stat.close()
            return [CN2_Region,Chromos,GC_Content]
        def Reads_Direction_Detect(flag):
            flag2=int(flag)
            if int(flag2)&16==0: 
                    direct_1='+'
            elif int(flag2)&16>0:
                    direct_1='-'
            if int(flag2)&32==0:
                    direct_2='+'
            elif int(flag2)&32>0: 
                    direct_2='-'
            return([direct_1,direct_2])
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
        def Reads_block_assignment_2(bps,letters,block1,block2,flank):
            bps2=[int(i) for i in bps]
            relative_bps=[i-numpy.min(bps2) for i in bps2]
            if Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)==Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank):
                return Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)
            elif not Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)==Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank):
                length_left=Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)[-1]-numpy.min([block1,block2])
                length_right=-Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank)[-2]+numpy.max([block1,block2])
                if not length_left<length_right:
                    return Read_Block_From_Position(bps,letters,0,numpy.min([block1,block2]),'left',flank)
                elif length_left<length_right:
                    return Read_Block_From_Position(bps,letters,0,numpy.max([block1,block2]),'right',flank)
        def letters_bps_produce(letters,bps,flank):
            letters_bps={}
            letters_relative_bps={}
            letters_bps['left']=[bps[0]-flank,bps[0]]
            letters_relative_bps['left']=[-flank,0]
            for i in range(len(bps)-1):
                letters_relative_bps[letters[i]]=[bps[i]-bps[0],bps[i+1]-bps[0]]
                letters_bps[letters[i]]=[bps[i],bps[i+1]]
            letters_bps['right']=[bps[-1],bps[-1]+flank]
            letters_relative_bps['right']=[bps[-1]-bps[0],bps[-1]-bps[0]+flank]
            return [letters_bps,letters_relative_bps]
        def Reads_block_assignment_3(bps,letters,letters_relative_bps,block1,block2,flank):
            bps2=[int(i) for i in bps]
            relative_bps=[-flank]+[i-numpy.min(bps2) for i in bps2]+[bps2[-1]-bps2[0]+flank]
            relative_letters=['left']+letters+['right']
            meanbl=int(numpy.mean([block1,block2]))
            if meanbl<relative_bps[0]:
                return ['leftError',relative_bps[0],relative_bps[0]]
            if not meanbl<relative_bps[-1]:
                return ['rightError', relative_bps[-1], relative_bps[-1]]
            else:
                resultbl=[relative_letters[i] for i in range(len(relative_bps)-1) if meanbl<relative_bps[i+1] and not meanbl<relative_bps[i]]
                return resultbl+letters_relative_bps[resultbl[0]]
        def Reads_block_assignment_1(bps,letters,position):
            if position<bps[0]-2*flank or position>bps[-1]+2*flank:
                return '0'
            else:
                if position<bps[0]+1 and position>bps[0]-2*flank-1:
                    return letters[0]
                elif position>bps[-1]-1 and position<bps[-1]+2*flank+1:
                    return letters[-1]
                else:
                    for i in range(len(bps)-1):
                        if not position<bps[i] and not position>bps[i+1]:
                            return letters[i]
        def RD_Index_ReadIn(ppre_Path,BamN, chromo, region):
            if not ppre_Path[-1]=='/':
                ppre_Path+='/'
            path_in=ppre_Path+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/RD_Stat/'
            file_in=BamN+'.'+chromo+'.RD.index'
            fin=open(path_in+file_in)
            pos1=int(region[0])
            pos2=int(region[1])
            while True:
                pin1=fin.readline().strip().split()
                if not pin1: break
                pin2=fin.readline().strip().split()
                reg1=int(pin1[0].split(':')[1].split('-')[0])
                reg2=int(pin1[0].split(':')[1].split('-')[1])
                if not pos1<reg1 and not pos2>reg2:
                    break
        def Full_Info_of_Reads_Product_3(Initial_Bam,temp_bp,temp_let,bamChr,target_region,Chr_Link):
            Letter_Double={}
            Pair_ThroughBP=[]
            Double_Read_ThroughBP=[]
            Single_Read_ThroughBP=[]
            blackList=[]
            fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,target_region[0]-flank,target_region[-1]+flank))
            num_of_reads=0
            while True:
                pbam=fbam.readline().strip().split()
                if not pbam: break
                if int(pbam[1])&4>0: continue
                if int(pbam[1])&1024>0:continue
                if not int(pbam[4])>QCAlign or int(pbam[1])&512>0:
                    blackList.append(pbam[0])
                    continue
                if pbam[0] in blackList: continue
                num_of_reads+=1
                if int(pbam[1])&8>0 or not pbam[6]=='=':
                    pos1=int(pbam[3])+low_qual_edge
                    pos2=int(pbam[3])+cigar2reaadlength(pbam[5])-low_qual_edge
                    block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
                    block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
                    if block1==block2:
                        BlockCov[block1]+=cigar2reaadlength(pbam[5])
                    else:
                        reg1a=temp_bp[temp_let.index(block1)]
                        reg1b=temp_bp[temp_let.index(block1)+1]
                        reg2a=temp_bp[temp_let.index(block2)]
                        reg2b=temp_bp[temp_let.index(block2)+1]
                        rela_1=pos1-low_qual_edge-temp_bp[temp_let.index(block1)]
                        rela_2=pos2+low_qual_edge-temp_bp[temp_let.index(block2)]
                        Single_Read_ThroughBP.append([block1,rela_1,block2,rela_2,pbam[5]])
                    if not pbam[6]=='=':               
                        if not pbam[0] in Chr_Link:
                            Chr_Link[pbam[0]]=[pbam[1:9]]
                        else:
                            Chr_Link[pbam[0]]+=[pbam[1:9]]
                elif int(pbam[1])&8==0:
                    if pbam[6]=='=':
                        if not pbam[0] in Letter_Double.keys():
                            Letter_Double[pbam[0]]=[pbam[:9]]
                        else:
                            if not pbam[:9] in Letter_Double[pbam[0]]:
                                Letter_Double[pbam[0]]+=[pbam[:9]]
                                if int(Letter_Double[pbam[0]][0][3])<int(Letter_Double[pbam[0]][1][3]):
                                    pos1=int(Letter_Double[pbam[0]][0][3])+low_qual_edge
                                    pos2=int(Letter_Double[pbam[0]][1][3])+cigar2reaadlength(Letter_Double[pbam[0]][1][5])-low_qual_edge
                                else:
                                    pos1=int(Letter_Double[pbam[0]][1][3])+low_qual_edge
                                    pos2=int(Letter_Double[pbam[0]][0][3])+cigar2reaadlength(Letter_Double[pbam[0]][0][5])-low_qual_edge
                                block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
                                block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
                                if block1==block2:
                                    BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][0][5])
                                    del Letter_Double[pbam[0]]
                                    blackList.append(pbam[0])
            fbam.close()
            for key in Letter_Double.keys():
                if key in blackList:
                    del Letter_Double[key]
                    continue
                if len(Letter_Double[key])==2:
                    pos1=int(Letter_Double[key][0][3])
                    pos2=int(Letter_Double[key][1][3])
                    if not pos1>pos2:
                        pos1=int(Letter_Double[key][0][3])
                        pos1b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                        pos2=int(Letter_Double[key][1][3])
                        pos2b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double[key][0][1])
                    elif pos1>pos2:
                        pos1=int(Letter_Double[key][1][3])
                        pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                        pos2=int(Letter_Double[key][0][3])
                        pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double[key][1][1])
                    block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
                    block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
                    block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
                    block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
                    rela_1=pos1-temp_bp[temp_let.index(block1)]
                    rela_2=pos2-temp_bp[temp_let.index(block2)]
                    rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
                    rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
                    if block1==block1b and block2==block2b:
                        Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
                    else:
                        Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
                    del Letter_Double[key]
                elif len(Letter_Double[key])==1:
                    if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
                        BlockCov[Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
                        del Letter_Double[key]
            Initial_DR_Penal=0
            for j in Pair_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            for j in Double_Read_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            for j in Pair_ThroughBP:
                Initial_Cov[j[0]]+=j[2]-j[1]
                Initial_Cov[j[3]]+=j[5]-j[4]
            for j in Single_Read_ThroughBP:
                Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                Initial_Cov[j[2]]+=j[3]
            for j in Double_Read_ThroughBP:
                if j[0]==j[2]:
                    Initial_Cov[j[0]]+=j[3]-j[1]
                else:
                    Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                    Initial_Cov[j[2]]+=j[3]
                if j[4]==j[6]:
                    Initial_Cov[j[4]]+=j[7]-j[5]
                else:
                    Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
                    Initial_Cov[j[6]]+=j[7]
            for j in Pair_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
            for j in Double_Read_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
            return [Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,num_of_reads,Initial_DR_Penal]
        def Full_Info_of_Reads_Product(Initial_Bam,bps,total_bps,total_letters,bamChr,flank,QCAlign,ReadLength,chr_link):
            temp_bp=total_bps
            temp_let=total_letters
            BlockCov={}
            for j in temp_let:
                BlockCov[j]=0
            Letter_Double={}
            Pair_ThroughBP=[]
            Double_Read_ThroughBP=[]
            Single_Read_ThroughBP=[]
            blackList=[]
            fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,bps[0]-flank,bps[-1]+flank))
            while True:
                pbam=fbam.readline().strip().split()
                if not pbam: break
                if int(pbam[1])&4>0: continue
                if int(pbam[1])&1024>0:continue
                if int(pbam[1])&512>0:
                    blackList.append(pbam[0])
                    continue
                if not int(pbam[4])>QCAlign:
                    continue
                if pbam[0] in blackList: continue
                if int(pbam[1])&8>0 or not pbam[6]=='=':
                    pos1=int(pbam[3])+low_qual_edge
                    pos2=int(pbam[3])+cigar2reaadlength(pbam[5])-low_qual_edge
                    block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
                    block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
                    if block1==block2:
                        BlockCov[block1]+=cigar2reaadlength(pbam[5])
                    else:
                        rela_1=pos1-low_qual_edge-temp_bp[temp_let.index(block1)]
                        rela_2=pos2+low_qual_edge-temp_bp[temp_let.index(block2)]
                        Single_Read_ThroughBP.append([block1,rela_1,block2,rela_2,pbam[5]])
                    if not pbam[6]=='=':
                        if not pbam[0] in chr_link.keys():
                            chr_link[pbam[0]]=[pbam[1:9]]
                        else:
                            chr_link[pbam[0]]+=[pbam[1:9]]
                elif int(pbam[1])&8==0:
                    if pbam[6]=='=':
                        if not pbam[0] in Letter_Double.keys():
                            Letter_Double[pbam[0]]=[pbam[:9]]
                        else:
                            if not pbam[:9] in Letter_Double[pbam[0]]:
                                Letter_Double[pbam[0]]+=[pbam[:9]]
                                if int(Letter_Double[pbam[0]][0][3])<int(Letter_Double[pbam[0]][1][3]):
                                    pos1=int(Letter_Double[pbam[0]][0][3])+low_qual_edge
                                    pos2=int(Letter_Double[pbam[0]][1][3])+cigar2reaadlength(Letter_Double[pbam[0]][1][5])-low_qual_edge
                                else:
                                    pos1=int(Letter_Double[pbam[0]][1][3])+low_qual_edge
                                    pos2=int(Letter_Double[pbam[0]][0][3])+cigar2reaadlength(Letter_Double[pbam[0]][0][5])-low_qual_edge
                                block1=Reads_block_assignment_1(temp_bp,temp_let,pos1)
                                block2=Reads_block_assignment_1(temp_bp,temp_let,pos2)
                                if block1==block2:
                                    BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][0][5])
                                    BlockCov[block1]+=cigar2reaadlength(Letter_Double[pbam[0]][1][5])
                                    del Letter_Double[pbam[0]]
                                    blackList.append(pbam[0])
            fbam.close()
            for key in Letter_Double.keys():
                if key in blackList:
                    del Letter_Double[key]
                    continue
                if len(Letter_Double[key])==2:
                    pos1=int(Letter_Double[key][0][3])
                    pos2=int(Letter_Double[key][1][3])
                    if not pos1>pos2:
                        pos1=int(Letter_Double[key][0][3])
                        pos1b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                        pos2=int(Letter_Double[key][1][3])
                        pos2b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double[key][0][1])
                    elif pos1>pos2:
                        pos1=int(Letter_Double[key][1][3])
                        pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                        pos2=int(Letter_Double[key][0][3])
                        pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double[key][1][1])
                    block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
                    block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
                    block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
                    block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
                    rela_1=pos1-temp_bp[temp_let.index(block1)]
                    rela_2=pos2-temp_bp[temp_let.index(block2)]
                    rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
                    rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
                    if block1==block1b and block2==block2b:
                        Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
                    else:
                        Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
                    del Letter_Double[key]
                elif len(Letter_Double[key])==1:
                    if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][7]))==0:
                        if Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
                            BlockCov[Reads_block_assignment_1(temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
                            del Letter_Double[key]
            Initial_DR_Penal=0
            for j in Pair_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            for j in Double_Read_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            Initial_Cov={}
            for j in temp_let:
                Initial_Cov[j]=0
            for j in Pair_ThroughBP:
                Initial_Cov[j[0]]+=j[2]-j[1]
                Initial_Cov[j[3]]+=j[5]-j[4]
            for j in Single_Read_ThroughBP:
                Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                Initial_Cov[j[2]]+=j[3]
            for j in Double_Read_ThroughBP:
                if j[0]==j[2]:
                    Initial_Cov[j[0]]+=j[3]-j[1]
                else:
                    Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                    Initial_Cov[j[2]]+=j[3]
                if j[4]==j[6]:
                    Initial_Cov[j[4]]+=j[7]-j[5]
                else:
                    Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
                    Initial_Cov[j[6]]+=j[7]
            Initial_IL=[]
            for j in Pair_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
            for j in Double_Read_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
            Initial_ILPenal=[]
            for j in Initial_IL:
                Initial_ILPenal+=[standard_pdf_IL_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
            return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov,Letter_Double]
        def Single_Rec_Read_Locate(Letter_Double_rec,temp_bp, temp_let):
            Pair_ThroughBP=[]
            Double_Read_ThroughBP=[]
            Single_Read_ThroughBP=[]
            Initial_IL=[]
            BlockCov={}
            Initial_Cov={}
            Initial_DR_Penal=0
            for j in temp_let:
                BlockCov[j]=0
            for key in Letter_Double_rec.keys():
                if len(Letter_Double_rec[key])==1:
                    pos1=int(Letter_Double_rec[key][0][3]) 
                    pos2=int(Letter_Double_rec[key][0][7])
                    bamChr=Letter_Double_rec[key][0][2]
                    fbamtemp=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,bamChr,pos2,pos2+ReadLength))
                    while True:
                        pbam=fbamtemp.readline().strip().split()
                        if not pbam: break
                        flag=0
                        if pbam[0]==key:
                            Letter_Double_rec[key]+=[pbam[:9]]
                            flag+=1
                        if flag==1:
                            break
                    fbamtemp.close()         
            for key in Letter_Double_rec.keys():
                if len(Letter_Double_rec[key])==2:
                    pos1=int(Letter_Double_rec[key][0][3])
                    pos2=int(Letter_Double_rec[key][1][3])
                    if not pos1>pos2:
                        pos1=int(Letter_Double_rec[key][0][3])
                        pos1b=pos1+cigar2reaadlength(Letter_Double_rec[key][0][5])
                        pos2=int(Letter_Double_rec[key][1][3])
                        pos2b=pos2+cigar2reaadlength(Letter_Double_rec[key][1][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double_rec[key][0][1])
                    elif pos1>pos2:
                        pos1=int(Letter_Double_rec[key][1][3])
                        pos1b=pos2+cigar2reaadlength(Letter_Double_rec[key][1][5])
                        pos2=int(Letter_Double_rec[key][0][3])
                        pos2b=pos1+cigar2reaadlength(Letter_Double_rec[key][0][5])
                        direct_temp=Reads_Direction_Detect(Letter_Double_rec[key][1][1])
                    if not pos1<temp_bp[0]-flank+1 and not pos2b>temp_bp[-1]+flank-1:
                        block1=Reads_block_assignment_1(temp_bp,temp_let,pos1+low_qual_edge)
                        block2=Reads_block_assignment_1(temp_bp,temp_let,pos2+low_qual_edge)
                        block1b=Reads_block_assignment_1(temp_bp,temp_let,pos1b-low_qual_edge)
                        block2b=Reads_block_assignment_1(temp_bp,temp_let,pos2b-low_qual_edge)
                        rela_1=pos1-temp_bp[temp_let.index(block1)]
                        rela_2=pos2-temp_bp[temp_let.index(block2)]
                        rela_1b=pos1b-temp_bp[temp_let.index(block1b)]
                        rela_2b=pos2b-temp_bp[temp_let.index(block2b)]
                        if block1==block1b==block2==block2:
                            BlockCov[block1]+=cigar2reaadlength(Letter_Double_rec[key][0][5])
                        else:                       
                            if block1==block1b and block2==block2b:
                                Pair_ThroughBP.append([block1,rela_1,rela_1b, block2,rela_2,rela_2b]+direct_temp)
                            else:
                                Double_Read_ThroughBP.append([block1,rela_1,block1b,rela_1b, block2,rela_2,block2b,rela_2b]+direct_temp)
                        del Letter_Double_rec[key]
            for j in Pair_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            for j in Double_Read_ThroughBP:
                if not j[-2:]==['+', '-']:
                    Initial_DR_Penal+=1
            for j in temp_let:
                Initial_Cov[j]=0
            for j in Pair_ThroughBP:
                Initial_Cov[j[0]]+=j[2]-j[1]
                Initial_Cov[j[3]]+=j[5]-j[4]
            for j in Single_Read_ThroughBP:
                Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                Initial_Cov[j[2]]+=j[3]
            for j in Double_Read_ThroughBP:
                if j[0]==j[2]:
                    Initial_Cov[j[0]]+=j[3]-j[1]
                else:
                    Initial_Cov[j[0]]+=temp_bp[temp_let.index(j[0])+1]-temp_bp[temp_let.index(j[0])]-j[1]
                    Initial_Cov[j[2]]+=j[3]
                if j[4]==j[6]:
                    Initial_Cov[j[4]]+=j[7]-j[5]
                else:
                    Initial_Cov[j[4]]+=temp_bp[temp_let.index(j[4])+1]-temp_bp[temp_let.index(j[4])]-j[5]
                    Initial_Cov[j[6]]+=j[7]
            Initial_IL=[]
            for j in Pair_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[3])]-temp_bp[temp_let.index(j[0])]-j[1]+j[5])
            for j in Double_Read_ThroughBP:
                Initial_IL.append(temp_bp[temp_let.index(j[6])]-temp_bp[temp_let.index(j[0])]-j[1]+j[7])
            Initial_ILPenal=[]
            for j in Initial_IL:
                Initial_ILPenal+=[standard_pdf_IL_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
            return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov]
        def letter_rearrange(bps2):
            chr_letter_bp={}
            let_start=96
            for i in bps2:
                if not i[0] in chr_letter_bp.keys():
                    chr_letter_bp[i[0]]={}
                for j in range(len(i))[1:-1]:
                    chr_letter_bp[i[0]][chr(let_start+j)]=[]
                    if int(i[j+1])-int(i[j])<10*flank:
                        chr_letter_bp[i[0]][chr(let_start+j)]+=[int(i[j]),int(i[j+1])]
                    else:
                        chr_letter_bp[i[0]][chr(let_start+j)]+=[int(i[j]),int(i[j])+flank,int(i[j+1])-flank,int(i[j+1])]
                let_start+=len(i)-2
            return chr_letter_bp
        def letter_GC_ReadIn(chr_letter_bp):
            block_GC_temp={}
            filein=ref_prefix+'.GC_Content'
            block_range={}
            GC_hash={}
            test_flag=0
            for i in chr_letter_bp.keys():
                if not os.path.isfile(filein+'.'+str(i)):
                    test_flag+=1
            if test_flag==0:
                for i in chr_letter_bp.keys():
                    GC_hash[i]={}
                    block_range[i]=[]
                    for j in chr_letter_bp[i].keys():
                        block_range[i]+=chr_letter_bp[i][j]
                    block_range[i]=[min(block_range[i]),max(block_range[i])]
                    fin=open(filein+'.'+str(i))
                    while True:
                        pin=fin.readline().strip().split()
                        if not pin: break
                        pin2=fin.readline().strip().split()
                        if 'chr' in block_range.keys()[0]:
                            if not 'chr' in pin[0]:
                                pin[0]='chr'+pin[0]
                        elif not 'chr' in block_range.keys()[0]:
                            if 'chr' in pin[0]:
                                pin[0]=pin[0][3:]
                        if pin[0] in block_range.keys():
                            if not int(pin[2])<block_range[pin[0]][0] and not int(pin[1])>block_range[pin[0]][1]:
                                GC_hash[pin[0]][pin[1]+'-'+pin[2]]=pin2
                    fin.close()
                for k1 in chr_letter_bp.keys():
                    block_GC_temp[k1]={}
                    for k2 in GC_hash[k1].keys():
                        bl2=[int(k2.split('-')[0]),int(k2.split('-')[1])]
                        for k3 in chr_letter_bp[k1].keys():
                            if min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])<bl2[1]+1:
                                block_GC_temp[k1][k3]=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/Window_Size:(max(chr_letter_bp[k1][k3])-bl2[0])/Window_Size+1]
                            elif min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])>bl2[1]:
                                if not k3 in block_GC_temp[k1].keys():
                                    block_GC_temp[k1][k3]=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/Window_Size:]
                                else:
                                    block_GC_temp[k1][k3]+=GC_hash[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/Window_Size:]                        
                            elif min(chr_letter_bp[k1][k3])<bl2[0] and max(chr_letter_bp[k1][k3])>bl2[0]-1:
                                if not k3 in block_GC_temp[k1].keys():
                                    block_GC_temp[k1][k3]=GC_hash[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/Window_Size+1]
                                else:
                                    block_GC_temp[k1][k3]+=GC_hash[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/Window_Size+1]                      
                            elif min(chr_letter_bp[k1][k3])<bl2[0]+1 and max(chr_letter_bp[k1][k3])>bl2[1]-1:
                                if not k3 in block_GC_temp[k1].keys():
                                    block_GC_temp[k1][k3]=GC_hash[k1][k2]
                                else:
                                    block_GC_temp[k1][k3]+=GC_hash[k1][k2]                     
                for k1 in block_GC_temp.keys():
                    for k2 in block_GC_temp[k1].keys():
                        if not block_GC_temp[k1][k2]==[]:
                            block_GC_temp[k1][k2]=numpy.mean([float(k3) for k3 in block_GC_temp[k1][k2]])
                        else:
                            return 'error'
                return block_GC_temp
            else:
                return 'error'
        def letter_RD_ReadIn(chr_letter_bp):
            test_flag=0
            for k1 in chr_letter_bp.keys():
                filein=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+k1+'.RD.index'
                if not os.path.isfile(filein):
                    test_flag+=1
            if test_flag==0:
                out={}
                RD_hash={}
                block_range={}
                for i in chr_letter_bp.keys():
                    RD_hash[i]={}
                    out[i]={}
                    block_range[i]=[]
                    for j in chr_letter_bp[i].keys():
                        block_range[i]+=chr_letter_bp[i][j]
                    block_range[i]=[min(block_range[i]),max(block_range[i])]
                for k1 in chr_letter_bp.keys():
                    filein=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+k1+'.RD.index'
                    fin=open(filein)
                    while True:
                        pin=fin.readline().strip().split()
                        if not pin: break
                        pin2=fin.readline().strip().split()
                        bl2=[int(pin[0].split(':')[1].split('-')[0]),int(pin[0].split(':')[1].split('-')[1])]
                        if not bl2[1]<block_range[k1][0]+1 and not bl2[0]>block_range[k1][1]-1:
                            RD_hash[k1][str(bl2[0])+'-'+str(bl2[1])]=pin2
                    fin.close()
                for k1 in chr_letter_bp.keys():
                    for k2 in RD_hash[k1].keys():
                        bl2=[int(k2.split('-')[0]),int(k2.split('-')[1])]
                        for j in sorted(chr_letter_bp[k1].keys()):
                            if not j in out[k1].keys():
                                out[k1][j]=[]
                            if len(chr_letter_bp[k1][j])==4:
                                bl1=chr_letter_bp[k1][j][1:-1]
                                if bl1[0]>bl2[0]-1 and bl1[1]<bl2[1]+1:
                                    out[k1][j]+=RD_hash[k1][k2][(bl1[0]-bl2[0])/Window_Size:(bl1[1]-bl2[0])/Window_Size+1]
                                elif bl1[0]>bl2[0]-1 and bl1[1]>bl2[1]:
                                    out[k1][j]+=RD_hash[k1][k2][(bl1[0]-bl2[0])/Window_Size:]
                                elif bl1[0]<bl2[0] and bl1[1]<bl2[1]+1:
                                    out[k1][j]+=RD_hash[k1][k2][:(bl1[1]-bl2[0])/Window_Size+1]
                                elif bl1[0]<bl2[0] and bl1[1]>bl2[1]:
                                    out[k1][j]+=RD_hash[k1][k2]
                for k1 in out.keys():
                    for k2 in out[k1].keys():
                        if out[k1][k2]==[]:
                            out[k1][k2]=0
                        else:
                            out[k1][k2]=numpy.mean([float(k3) for k3 in out[k1][k2]])
                return out
            else:
                return 'error'
        def letter_range_report(chr_letter_bp):
            ks_range={}
            for k1 in chr_letter_bp.keys():
                ks_range[k1]=[]
                for k2 in chr_letter_bp[k1].keys():
                    ks_range[k1]+=chr_letter_bp[k1][k2]
                ks_range[k1]=[min(ks_range[k1]),max(ks_range[k1])]
            for k1 in chr_letter_bp.keys():
                chr_letter_bp[k1]['left']=[ks_range[k1][0]-flank,ks_range[k1][0]]
                chr_letter_bp[k1]['right']=[ks_range[k1][1],ks_range[k1][1]+flank]
        def block_Read_From_Bam(chr_letter_bp):
            blocks_out={}
            for k1 in chr_letter_bp.keys():
                k_block=[]
                blocks_out[k1]=[]
                for k2 in chr_letter_bp[k1].keys():
                    if not k2 in ['left','right']:
                        k_block.append(k2)
                blocks_out[k1].append(chr_letter_bp[k1]['left'][:]+['left'])
                k_block.sort()
                for k2 in k_block+['right']:
                    if len(chr_letter_bp[k1][k2])==2:
                        blocks_out[k1][-1]+=chr_letter_bp[k1][k2]+[k2]
                    elif len(chr_letter_bp[k1][k2])==4:
                        blocks_out[k1][-1]+=chr_letter_bp[k1][k2][:2]+[k2]
                        blocks_out[k1].append(chr_letter_bp[k1][k2][2:]+[k2])
            for k1 in blocks_out.keys():
                for k2 in blocks_out[k1]:
                    k2.sort()
            return blocks_out
        tolerance_bp=10
        def pos_block_assign(block_bps_chr,read_pos):
            read_new=[]
            for i in read_pos:
                if type(i)==type(1):
                    read_new.append(i)
            for j in range(len(read_new)/2):
                read_new[2*j]=read_new[2*j]+tolerance_bp
                read_new[2*j+1]=read_new[2*j+1]-tolerance_bp
            for i in read_new:
                if type(i)==type(1):
                    for j in block_bps_chr.keys():
                        if j=='left':
                            if i<block_bps_chr[j][1]:
                                read_pos.append(j)
                        elif j=='right':
                            if i>block_bps_chr[j][0]-1:
                                read_pos.append(j)
                        else:
                            if block_bps_chr[j][0]-1<i and block_bps_chr[j][1]>i:
                                read_pos.append(j)
        def block_Info_ReadIn(chr_letter_bp,blocks_read_in):
            block_bps={}
            block_rds={}
            for k1 in chr_letter_bp.keys():
                block_bps[k1]={}
                block_rds[k1]={}
                for k2 in chr_letter_bp[k1].keys():
                    block_bps[k1][k2]=[min(chr_letter_bp[k1][k2]),max(chr_letter_bp[k1][k2])]
                    block_rds[k1][k2]=0
            Pair_ThroughBP={}
            Double_Read_ThroughBP={}
            Single_Read_ThroughBP={}
            total_rec={}
            rd_low_qual={}
            for k1 in chr_letter_bp.keys():
                Pair_ThroughBP[k1]=[]
                Double_Read_ThroughBP[k1]=[]
                Single_Read_ThroughBP[k1]=[]
                rd_low_qual[k1]={}
                for k2 in blocks_read_in[k1]:
                    k2a=[]
                    k2b=[]
                    for k3 in k2:
                        if type(k3)==type(1):
                            k2a.append(k3)
                        else:
                            k2b.append(k3)
                    fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(Initial_Bam,k1,min(k2a)-flank,max(k2a)+flank))
                    blackList=[]
                    temp_rec={} 
                    temp_rec_LowQual={}     
                    while True:
                        pbam=fbam.readline().strip().split()
                        if not pbam: break
                        if int(pbam[1])&4>0: continue
                        if int(pbam[1])&1024>0:continue
                        if int(pbam[1])&512>0:
                            blackList.append(pbam[0])
                            continue
                        if pbam[0] in blackList: continue
                        if not int(pbam[4])>QCAlign:
                            if not pbam[0] in temp_rec_LowQual.keys():
                                temp_rec_LowQual[pbam[0]]=[]
                            if not pbam[1:9] in temp_rec_LowQual[pbam[0]]:
                                temp_rec_LowQual[pbam[0]]+=[pbam[1:9]]
                        else:
                            if not pbam[0] in temp_rec.keys():
                                temp_rec[pbam[0]]=[]
                            if not pbam[1:9] in temp_rec[pbam[0]]:
                                temp_rec[pbam[0]]+=[pbam[1:9]]
                    fbam.close()
                    flank_region=[]
                    for k3 in k2b:
                        flank_region+=block_bps[k1][k3]
                    flank_region=[min(flank_region),max(flank_region)]
                    for k3 in temp_rec_LowQual.keys():
                        for k4 in temp_rec_LowQual[k3]:
                            read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
                            pos_block_assign(block_bps[k1],read_pos)
                            if read_pos[-1]==read_pos[-2]:
                                if not read_pos[-1] in rd_low_qual[k1].keys():
                                    rd_low_qual[k1][read_pos[-1]]=0
                                rd_low_qual[k1][read_pos[-1]]+=(read_pos[1]-read_pos[0])
                            else:
                                if not read_pos[-2] in rd_low_qual[k1].keys():
                                    rd_low_qual[k1][read_pos[-2]]=0
                                if not read_pos[-1] in rd_low_qual[k1].keys():
                                    rd_low_qual[k1][read_pos[-1]]=0
                                rd_low_qual[k1][read_pos[-2]]+=block_bps[k1][read_pos[-2]][1]-read_pos[0]
                                rd_low_qual[k1][read_pos[-1]]+=-block_bps[k1][read_pos[-1]][0]+read_pos[1]
                    for k3 in temp_rec.keys():
                        if len(temp_rec[k3])>2:
                            test_rec=[int(temp_rec[k3][0][7])]
                            test_rec2=[temp_rec[k3][0]]
                            test_let=0
                            for k4 in temp_rec[k3][1:]:
                                delflag=0
                                for k5 in test_rec:
                                    if int(k4[7])+k5==0:
                                        test_let+=1
                                        k6=k3+chr(96+test_let)
                                        temp_rec[k6]=[test_rec2[test_rec.index(k5)],k4]
                                        del test_rec2[test_rec.index(k5)]
                                        del test_rec[test_rec.index(k5)]
                                        delflag+=1
                                if delflag==0:
                                    test_rec.append(int(k4[7]))
                                    test_rec2.append(k4)
                            temp_rec[k3]=test_rec2 
                    for k3 in temp_rec.keys():
                        if len(temp_rec[k3])==1:
                            del_flag=0
                            k4=temp_rec[k3][0]
                            read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
                            mate_pos=[int(k4[6]),int(k4[6])+ReadLength]
                            if 'left' in k2b and mate_pos[1]<flank_region[0]:
                                del_flag+=1
                            elif 'right' in k2b and mate_pos[0]>flank_region[0]:
                                del_flag+=1
                            if del_flag>0:
                                del temp_rec[k3]
                                pos_block_assign(block_bps[k1],read_pos)
                                if read_pos[-1]==read_pos[-2]:
                                    block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                                else:
                                    Single_Read_ThroughBP[k1].append(read_pos)
                            else:
                                if not k3 in total_rec.keys():
                                    total_rec[k3]=[k4]
                                else:
                                    total_rec[k3]+=[k4]
                        elif len(temp_rec[k3])==2: 
                            if int(temp_rec[k3][0][7])==0 or int(temp_rec[k3][1][7])==0:
                                continue
                            if int(temp_rec[k3][0][7])+int(temp_rec[k3][1][7])==0 and int(temp_rec[k3][0][7])<0:
                                temp_rec[k3]=[temp_rec[k3][1],temp_rec[k3][0]]
                            read_pos=[int(temp_rec[k3][0][2]),int(temp_rec[k3][0][2])+cigar2reaadlength(temp_rec[k3][0][4]),int(temp_rec[k3][1][2]),int(temp_rec[k3][1][2])+cigar2reaadlength(temp_rec[k3][1][4])]+Reads_Direction_Detect(temp_rec[k3][0][0])
                            if read_pos[0]>read_pos[2]:
                                read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
                            pos_block_assign(block_bps[k1],read_pos)
                            if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
                                block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                                block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
                            elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
                                Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
                            else:
                                Double_Read_ThroughBP[k1].append(read_pos)
                            del temp_rec[k3]
            for k3 in total_rec.keys():
                if len(total_rec[k3])==1: 
                    del_flag=0
                    k4=total_rec[k3][0]
                    read_pos=[int(k4[2]),int(k4[2])+cigar2reaadlength(k4[4])]
                    mate_pos=[int(k4[6]),int(k4[6])+ReadLength]
                    if 'left' in k2b and mate_pos[1]<flank_region[0]:
                        del_flag+=1
                    elif 'right' in k2b and mate_pos[0]>flank_region[0]:
                        del_flag+=1
                    elif not mate_pos[1]<flank_region[0] and not mate_pos[0]>flank_region[1]:
                        del_flag+=1
                    if del_flag>0:
                        del total_rec[k3]
                        pos_block_assign(block_bps[k1],read_pos)
                        if read_pos[-1]==read_pos[-2]:
                            block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                        else:
                            Single_Read_ThroughBP[k1].append(read_pos)
                elif len(total_rec[k3])==2:
                    read_pos=[int(total_rec[k3][0][2]),int(total_rec[k3][0][2])+cigar2reaadlength(total_rec[k3][0][4]),int(total_rec[k3][1][2]),int(total_rec[k3][1][2])+cigar2reaadlength(total_rec[k3][1][4])]+Reads_Direction_Detect(total_rec[k3][0][0])
                    if read_pos[0]>read_pos[2]:
                        read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
                    pos_block_assign(block_bps[k1],read_pos)
                    if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
                        block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                        block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
                    elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
                        Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
                    else:
                        Double_Read_ThroughBP[k1].append(read_pos)
                    del total_rec[k3]
            direction_penal=0
            block_rd2={}
            block_bp2=block_bps
            for k1 in block_rds.keys():
                block_rd2[k1]={}
                for k2 in block_rds[k1].keys():
                    block_rd2[k1][k2]=0
            for i2 in Pair_ThroughBP.keys():
                for i in Pair_ThroughBP[i2]:
                    if not i[4:6]==['+','-']:
                        direction_penal+=1
                    block_rd2[i2][i[6]]+=i[1]-i[0]
                    block_rd2[i2][i[7]]+=i[3]-i[2]
            for i2 in Double_Read_ThroughBP.keys():
                for i in Double_Read_ThroughBP[i2]:
                    if i[6]==i[7]:
                        block_rd2[i2][i[6]]+=i[1]-i[0]
                        block_rd2[i2][i[8]]+=-i[2]+block_bp2[i2][i[8]][1]
                        block_rd2[i2][i[9]]+=i[3]-block_bp2[i2][i[9]][0]
                    elif i[8]==i[9]:
                        block_rd2[i2][i[8]]+=i[3]-i[2]
                        block_rd2[i2][i[6]]+=-i[0]+block_bp2[i2][i[6]][1]
                        block_rd2[i2][i[7]]+=i[1]-block_bp2[i2][i[7]][0]
                    else:
                        block_rd2[i2][i[6]]+=-i[0]+block_bp2[i2][i[6]][1]
                        block_rd2[i2][i[7]]+=i[1]-block_bp2[i2][i[7]][0]
                        block_rd2[i2][i[8]]+=-i[2]+block_bp2[i2][i[8]][1]
                        block_rd2[i2][i[9]]+=i[3]-block_bp2[i2][i[9]][0]
            for i2 in Single_Read_ThroughBP.keys():
                for i in Single_Read_ThroughBP[i2]:
                    block_rd2[i2][i[2]]+=-i[0]+block_bp2[i2][i[2]][1]
                    block_rd2[i2][i[3]]+=i[1]-block_bp2[i2][i[3]][0]
            for k1 in rd_low_qual.keys():
                for k2 in rd_low_qual[k1].keys():
                    block_rds[k1][k2]+=rd_low_qual[k1][k2]
            return [block_rds,block_rd2,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP]
        def total_rd_calcu(letter_RD2,letter_GC,chr_letter_bp,block_rd2):
            out={}
            for k1 in block_rd2.keys():
                for k2 in block_rd2[k1].keys():
                    if not k2 in ['left','right']:
                        out[k2]=[]
            for k2 in letter_RD2.keys():
                out[k2].append(letter_RD2[k2])
            for k1 in block_rd2.keys():
                for k2 in block_rd2[k1].keys():
                    if not k2 in ['left','right']:
                        if not chr_letter_bp[k1][k2][-1]==chr_letter_bp[k1][k2][0]:
                            out[k2].append(float(block_rd2[k1][k2])/float(chr_letter_bp[k1][k2][-1]-chr_letter_bp[k1][k2][0]))
            for k2 in out.keys():
                out[k2]=numpy.sum(out[k2])
            GC2={}
            for k1 in letter_GC.keys():
                for k2 in letter_GC[k1].keys():
                    if not k2 in ['left','right']:
                        GC2[k2]=letter_GC[k1][k2]
            out2=GC_RD_Adj_hash(GC_Median_Num,GC_Overall_Median_Num,k1,GC2,out)
            return out2
        def GC_RD_Adj_hash(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
            Coverage_af_Adj={}
            Overall_Median_Coverage=float(GC_Overall_Median_Num)
            for key_1 in GC_Content.keys():
                    Be_Adj_RD=Coverage[key_1]
                    GC_Con=int(round(GC_Content[key_1]*100))
                    if GC_Con in GC_Median_Num.keys():
                            Median_Coverage=GC_Median_Num[GC_Con]
                            Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
                    elif not GC_Con in GC_Median_Num.keys():
                            Af_Adj_RD=Overall_Median_Coverage
                    Coverage_af_Adj[key_1]=Af_Adj_RD
            return  Coverage_af_Adj
        def DR_Penal_Calcu(read_info):
            DR_Penal=0
            for i2 in read_info[2].keys():
                for i in read_info[2][i2]:
                    if i[4:6]==['+', '+'] or i[4:6]==['-', '-'] or i[4:6]==['-', '+']:
                        DR_Penal+=1
            for i2 in read_info[3].keys():
                for i in read_info[3][i2]:
                    if i[4:6]==['+', '+'] or i[4:6]==['-', '-'] or i[4:6]==['-', '+']:
                        DR_Penal+=1
            return DR_Penal
        def IL_Penal_Calcu(read_info):
            Initial_IL=[]
            for j2 in read_info[2].keys():
                for j in read_info[2][j2]:
                    Initial_IL.append(j[3]-j[1])
            for j2 in read_info[3].keys():
                for j in read_info[3][j2]:
                    Initial_IL.append(j[3]-j[1])
            Initial_ILPenal=[]
            for j in Initial_IL:
                Initial_ILPenal+=[standard_pdf_IL_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
            return Initial_ILPenal
        def original_bp_let_produce(chr_letter_bp,bps2):
            all_bp=[int(i) for i in bps2[0][1:]]
            all_let=[chr(97+i) for i in range(len(all_bp)-1)]
            origin=int(bps2[0][1])
            for i in bps2[1:]:
                temp_ks=[]
                for k2 in chr_letter_bp[i[0]].keys():
                    if not k2 in ['left','right']:
                        temp_ks.append(k2)
                for k2 in sorted(temp_ks):
                    all_let.append(k2)
                    all_bp.append(all_bp[-1]+chr_letter_bp[i[0]][k2][-1]-chr_letter_bp[i[0]][k2][0])
            return [all_bp,all_let] 
        def rela_Pair_ThroughBP(chr_letter_bp,Pair_ThroughBP):
            out=[]
            for k1 in Pair_ThroughBP.keys():
                for k2 in Pair_ThroughBP[k1]:
                    rela=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                                k2[1]-chr_letter_bp[k1][k2[6]][0],
                          k2[7],k2[2]-chr_letter_bp[k1][k2[7]][0],
                                k2[3]-chr_letter_bp[k1][k2[7]][0],k2[4],k2[5]]
                    out.append(rela)
            return out 
        def rela_Pair_Double_Read_ThroughBP(chr_letter_bp,Double_Read_ThroughBP):
            out=[]
            for k1 in Double_Read_ThroughBP.keys():
                 for k2 in Double_Read_ThroughBP[k1]:
                    if k2[8]==k2[9]:
                        rela1=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[6],k2[1]-chr_letter_bp[k1][k2[6]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        rela2=[k2[7],k2[0]-chr_letter_bp[k1][k2[7]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        temp_out=[rela1,rela2]
                    elif k2[6]==k2[7]:
                        rela1=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[8],k2[3]-chr_letter_bp[k1][k2[8]][0],k2[4],k2[5]]
                        rela2=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[9],k2[2]-chr_letter_bp[k1][k2[9]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        temp_out=[rela1,rela2]
                    else:
                        rela1=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[6],k2[1]-chr_letter_bp[k1][k2[6]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[8],k2[3]-chr_letter_bp[k1][k2[8]][0],k2[4],k2[5]]
                        rela2=[k2[7],k2[0]-chr_letter_bp[k1][k2[7]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[8],k2[3]-chr_letter_bp[k1][k2[8]][0],k2[4],k2[5]]
                        rela3=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[6],k2[1]-chr_letter_bp[k1][k2[6]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        rela4=[k2[7],k2[0]-chr_letter_bp[k1][k2[7]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[9],k2[2]-chr_letter_bp[k1][k2[9]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        temp_out=[rela1,rela2,rela3,rela4]
                    out.append(temp_out)
            return out
        def read_Pair_Single_Read_ThroughBP(chr_letter_bp,Single_Read_ThroughBP):
            out=[]
            for k1 in Single_Read_ThroughBP.keys():
                for k2 in Single_Read_ThroughBP[k1]:
                    rela=[k2[2],k2[0]-chr_letter_bp[k1][k2[2]][0],
                          k2[3],k2[1]-chr_letter_bp[k1][k2[3]][0]]
                    out.append(rela)
            return out
        def Full_Info_of_Reads_Integrate(bps2):
            chr_letter_bp=letter_rearrange(bps2)
            letter_GC=letter_GC_ReadIn(chr_letter_bp)
            letter_RD=letter_RD_ReadIn(chr_letter_bp)
            letter_range_report(chr_letter_bp)
            blocks_read_in=block_Read_From_Bam(chr_letter_bp)
            read_info=block_Info_ReadIn(chr_letter_bp,blocks_read_in)
            block_rds=read_info[0]
            block_rd2=read_info[1]
            letter_RD2={}
            for k1 in letter_RD.keys():
                for k2 in letter_RD[k1].keys():
                        if len(chr_letter_bp[k1][k2])==4:
                                letter_RD2[k2]=letter_RD[k1][k2]*(chr_letter_bp[k1][k2][2]-chr_letter_bp[k1][k2][1])/(chr_letter_bp[k1][k2][3]-chr_letter_bp[k1][k2][0])
                        else:
                                letter_RD2[k2]=letter_RD[k1][k2]
            for k1 in block_rds.keys():
                for k2 in block_rds[k1].keys():
                        if not k2 in ['left','right']:
                            if not chr_letter_bp[k1][k2][-1]==chr_letter_bp[k1][k2][0]:
                                letter_RD2[k2]+=float(block_rds[k1][k2])/float(chr_letter_bp[k1][k2][-1]-chr_letter_bp[k1][k2][0])
            Pair_ThroughBP=rela_Pair_ThroughBP(chr_letter_bp,read_info[2])
            Double_Read_ThroughBP=rela_Pair_Double_Read_ThroughBP(chr_letter_bp,read_info[3])
            Single_Read_ThroughBP=read_Pair_Single_Read_ThroughBP(chr_letter_bp,read_info[4])
            Initial_RD=total_rd_calcu(letter_RD2,letter_GC,chr_letter_bp,block_rd2)
            DR_Penal=DR_Penal_Calcu(read_info)
            IL_Penal=IL_Penal_Calcu(read_info)
            letter_GC_out={}
            for k1 in letter_GC.keys():
                    for k2 in letter_GC[k1].keys():
                            letter_GC_out[k2]=letter_GC[k1][k2]
            return [letter_RD2,Initial_RD,DR_Penal,numpy.mean(IL_Penal),Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,letter_GC_out]+original_bp_let_produce(chr_letter_bp,bps2)
        def GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
            Coverage_af_Adj=[]
            GC_Content=[float(i) for i in GC_Content]
            Coverage=[float(i) for i in Coverage]
            Overall_Median_Coverage=float(GC_Overall_Median_Num)
            for key_1 in range(len(GC_Content)):
                    Be_Adj_RD=Coverage[key_1]
                    GC_Con=GC_Content[key_1]
                    GC_Con=int(round(GC_Con*Window_Size))
                    if GC_Con in GC_Median_Num.keys():
                            Median_Coverage=GC_Median_Num[GC_Con]
                            Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
                    elif not GC_Con in GC_Median_Num.keys():
                            Af_Adj_RD=Overall_Median_Coverage
                    Coverage_af_Adj+=[Af_Adj_RD]
            return Coverage_af_Adj
        def GC_RD_Correction(chrbam):
            ref_index=ref_prefix+'GC_Content'
            cov_index=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/RD_Stat/'+BamN+'.'+chrbam+'.RD.index'
            fref=open(ref_index)
            fcov=open(cov_index)
            GC=[]
            RD=[]
            while True:
                pref1=fref.readline().strip().split()
                if not pref1: break
                pref2=fref.readline().strip().split()
                if pref1[0]==chrbam: 
                    GC.append([int(pref1[1]),int(pref1[2])]+pref2)
            while True:
                pcov1=fcov.readline().strip().split()
                if not pcov1: break
                pcov2=fcov.readline().strip().split()
                reg1=int(pcov1[0].split(':')[1].split('-')[0])
                reg2=int(pcov1[0].split(':')[1].split('-')[1])
                RD.append([reg1,reg2]+pcov2)
            fref.close()
            fcov.close()
            rd_rec=0
            RD2={}
            for i in GC:
                RD2[i[0]]=[]
                region1=i[:2]
                while True:
                    if rd_rec==len(RD): break
                    region2=RD[rd_rec]
                    if region2[1]-1<region1[1]:
                        RD2[i[0]].append(region2[:2])
                        rd_rec+=1
                    else:
                        break
            RD3=[]
            rd_rec=-1
            for i in GC:
                rd_rec+=1
                RD3.append(i[:2])
                j=RD2[i[0]][0]
                j2=RD[rd_rec]
                RD3[-1]+=j2[2:]     
                for j in RD2[i[0]][1:]:
                    rd_rec+=1
                    j2=RD[rd_rec]
                    RD3[-1]+=j2[12:]
            RD4=[]
            for i in range(len(RD3)):
                RD4.append(RD3[i][:2])
                RD4[-1]+=GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrbam,GC[i][2:],RD3[i][2:])
            return [GC,RD3,RD4]
        def GC_index(chrbam):
            ref_index=ref_prefix+'GC_Content'
            fref=open(ref_index)
            Cov=[]
            RD=[]
            while True:
                pref1=fref.readline().strip().split()
                if not pref1: break
                pref2=fref.readline().strip().split()
                if pref1[0]==chrbam: 
                    Cov.append([int(pref1[1]),int(pref1[2])]+pref2)
            return Cov
        def Reads_Direction_Detect(flag):
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
        def Reads_block_assignment(bps,letters,block):
            if numpy.mean(block[:2])<bps[0]:
                return 'left'
            elif not numpy.mean(block[:2])<bps[-1]:
                return 'right'
            else:
                for i in letters:
                    if not numpy.mean(block[:2])<bps[letters.index(i)] and numpy.mean(block[:2])<bps[letters.index(i)+1]:
                        return i
        def clusterNums(data, ClusterLen, direction):
            data.sort()
            out2=[]
            out3=[]
            out=[]
            for i in data:
                if len(out)==0:
                    out.append(i)
                else:
                    if i-out[-1]< ClusterLen:
                        out.append(i)
                    else:
                        if out[-1]-out[0]<ClusterLen:
                            if direction=='f':
                                out2.append(out[-1])
                                out3.append(len(out))
                            elif direction=='r':
                                out2.append(out[0])
                                out3.append(len(out))
                        else:
                            temp=[0]
                            lenTime=int((out[-1]-out[0])/ClusterLen)
                            lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                            lenInter2=sorted(lenInter)[::-1][:lenTime]
                            for j in range(len(lenInter)):
                                if lenInter[j] in  lenInter2:
                                    temp.append(j+1)
                            temp.append(len(out))
                            if direction=='f':
                                for k in range(len(temp)-1):
                                    out2.append(out[temp[k]:temp[k+1]][-1])
                                    out3.append(temp[k+1]-temp[k])
                            elif direction=='r':
                                for k in range(len(temp)-1):
                                    out2.append(out[temp[k]:temp[k+1]][0])
                                    out3.append(temp[k+1]-temp[k])
                        out=[i]
            if out[-1]-out[0]<ClusterLen:
                if direction=='f':
                    out2.append(out[-1])
                    out3.append(len(out))
                elif direction=='r':
                    out2.append(out[0])
                    out3.append(len(out))
            else:
                temp=[0]
                lenTime=int((out[-1]-out[0])/ClusterLen)
                lenInter=[out[j+1]-out[j] for j in range(len(out)-1)]
                lenInter2=sorted(lenInter)[::-1][:lenTime]
                for j in range(len(lenInter)):
                    if lenInter[j] in  lenInter2:
                        temp.append(j+1)
                temp.append(len(out))
                if direction=='f':
                    for k in range(len(temp)-1):
                        out2.append(out[temp[k]:temp[k+1]][-1])
                        out3.append(temp[k+1]-temp[k])
                elif direction=='r':
                    for k in range(len(temp)-1):
                        out2.append(out[temp[k]:temp[k+1]][0])
                        out3.append(temp[k+1]-temp[k])
            return [out2,out3]
        def clusterQC(hash, CluNum):
            out=[]
            for i in range(len(hash[1])):
                if not hash[1][i]<CluNum:
                    out.append(hash[0][i])
            return out
        def Single_Read_Assort_For_insert(Full_Info,bp_list,flank):
            relative_bps=[i-bp_list[0] for i in bp_list]
            letter_list=[chr(97+i) for i in range(len(bp_list)-1)]
            Block_and_Reads={}
            Block_and_Reads['left']=[]
            Block_and_Reads['right']=[]
            SingleR_Through=Full_Info[6]
            Pair_Through=Full_Info[4]
            Read_Through=Full_Info[5]
            for block in letter_list:
                    Block_and_Reads[block]=[]
            for j in Pair_Through:
                Block_and_Reads[j[0]]=[j[1:3],j[3:]]
                Block_and_Reads[j[3]]=[j[4:6],j[:3]+j[6:8]]
            for j in Read_Through:
                Block_and_Reads[j[0]]=[]
            for key in Full_Info_of_Reads.keys():
                    read_left=[int(i) for i in Full_Info_of_Reads[key][:2]]+[Full_Info_of_Reads[key][-2]]
                    read_right=[int(i) for i in Full_Info_of_Reads[key][2:4]]+[Full_Info_of_Reads[key][-1]]
                    assign_left=Reads_block_assignment_2(relative_bps,letter_list,read_left[0],read_left[1],flank)
                    assign_right=Reads_block_assignment_2(relative_bps,letter_list,read_right[0],read_right[1],flank)
                    New_Info=['_'.join([assign_left[0],str(int(co)-assign_left[1])]) for co in Full_Info_of_Reads[key][:2]]+['_'.join([assign_right[0],str(int(co)-assign_right[1])]) for co in Full_Info_of_Reads[key][2:4]]+Full_Info_of_Reads[key][4:]
                    Block_and_Reads[assign_left[0]][key]=New_Info
                    Block_and_Reads[assign_right[0]][key]=New_Info
            return Block_and_Reads
        def complementary(seq):
            seq2=[]
            for i in seq:
                if i in 'ATGCN':
                    seq2.append('ATGCN'['TACGN'.index(i)])
                elif i in 'atgcn':
                    seq2.append('atgcn'['tacgn'.index(i)])
            return ''.join(seq2)        
        def Insert_Seq_Pool_Prod_2(original_bp_list,ori_1_Seq,flank):
            ini_letters=['left']+['I'+chr(97+i) for i in range(len(original_bp_list)-1)]+['right']+['I'+chr(97+i)+'^' for i in range(len(original_bp_list)-1)]
            relative_bps=[0]+[j-original_bp_list[0]+flank for j in original_bp_list]+[original_bp_list[-1]+flank-original_bp_list[0]+flank]
            Insert_Seq_Pool={}
            for k in range(len(original_bp_list)+1):
                Insert_Seq_Pool[ini_letters[k]]=ori_1_Seq[relative_bps[k]:relative_bps[k+1]]
            for k in range(len(original_bp_list)+1,len(ini_letters)):
                Insert_Seq_Pool[ini_letters[k]]=complementary(ori_1_Seq[relative_bps[k-len(original_bp_list)]:relative_bps[k+1-len(original_bp_list)]])
            return Insert_Seq_Pool
        def Median_Pick(number_list):
            if 2*(len(number_list)/2)+1==len(number_list):
                return float(sorted(number_list)[len(number_list)/2])
            elif 2*(len(number_list)/2)==len(number_list):
                return float(sorted(number_list)[len(number_list)/2-1]+sorted(number_list)[len(number_list)/2])/float(2)
        def Block_Assign_To_Letters(bp_list,letter_list,flank):
            number_of_blocks=(numpy.max(bp_list)-numpy.min(bp_list)+2*flank)/Window_Size+1
            blocks={}
            bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
            relative_bp_list=[i-numpy.min(bp_list_new) for i in bp_list_new]
            bp_length=[(bp_list_new[i+1]-bp_list_new[i]) for i in range(len(bp_list_new)-1)]
            letter_list_new=['left']+letter_list+['right']
            bp_blocks=[[letter_list_new[j]]+range(relative_bp_list[j]/Window_Size,relative_bp_list[j+1]/Window_Size+1) for j in range(len(relative_bp_list)-1)]
            blocks_bp={}
            for i in range(number_of_blocks):
                blocks_bp[i+1]=[bp_list_new[0]+i*Window_Size,bp_list_new[0]+i*Window_Size+99]
                for j in bp_blocks:
                    if i in j:
                        blocks_bp[i+1].append(j[0])
            blocks_bp[0]=[blocks_bp[1][0]-Window_Size,blocks_bp[1][0]-1,'0']
            blocks_bp[number_of_blocks+1]=[blocks_bp[number_of_blocks][1]+1,blocks_bp[number_of_blocks][1]+Window_Size,'0']
            return blocks_bp
        def GC_Content_Calculate(seq2):
            NumAT=0
            NumGC=0
            for n in seq2:
                if n=='A' or n=='T' or n=='a' or n=='t':
                    NumAT+=1
                elif n=='G' or n=='C' or n=='g' or n=='c':
                    NumGC+=1
            return [NumGC,NumAT]            
        def Float_2_Integer(Number):
            if Number-int(Number)<0.5:
                K2=int(Number)
            elif not Number-int(Number)<0.5:
                K2=int(Number)+1
            return K2
        def c_Coverage_Calculate_InfoList(Full_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
            bp_M=[i-original_bp_list[0] for i in bp_MP[0]]
            bp_P=[i-original_bp_list[0] for i in bp_MP[1]]
            M_New_bp=[bp_M[0]-flank]+bp_M+[bp_M[-1]+flank]
            P_New_bp=[bp_P[0]-flank]+bp_P+[bp_P[-1]+flank]
            M_coverage=Block_Assign_To_Letters(bp_MP[0],letter_MP[0],flank)
            P_coverage=Block_Assign_To_Letters(bp_MP[1],letter_MP[1],flank)
            for key in M_coverage.keys():
                M_coverage[key].append(0)
            for key in P_coverage.keys():
                P_coverage[key].append(0)       
            for key in Half_Info.keys():
                Half=Half_Info[key]
                if Half[0]<-flank-Window_Size: continue
                else:
                    if Half[-1]=='M':
                        M_coverage[(Half[0]-(M_New_bp[0]))/Window_Size+1][-1]+=1
                    elif Half[-1]=='P':
                        P_coverage[(Half[0]-(P_New_bp[0]))/Window_Size+1][-1]+=1
            return [M_coverage,P_coverage]
        def c_GCContent_Calculate_InfoList(Ori_1_Seq,original_bp_list,flank):
            region_length=original_bp_list[-1]-original_bp_list[0]+2*flank
            region_length_new=(region_length/Window_Size+1)*Window_Size-2*flank
            Number_Of_Blocks=len(Ori_1_Seq)/Window_Size
            GC_Content={}
            for i in range(Number_Of_Blocks):
                GC_Content[i+1]=GC_Content_Calculate(Ori_1_Seq[i*Window_Size:(i+1)*Window_Size])[0]
            return GC_Content
        def left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj,flank):
            left_blocks=range(flank/Window_Size+1)[1:]
            left_Changes=[Af_GCRD_Adj[j][3]-Through_GCRD_Adj[j][3] for j in left_blocks]
            return numpy.mean(left_Changes)
        def All_Block_RD(Initial_block_RD,Af_GCRD_Adj,Af_block_RD,Af_Letter,flank):
            All_Letters=['left']+[chr(97+i) for i in range(len(Initial_block_RD)-1)]
            CNm=[1]+[0 for j in range(len(Initial_block_RD)-1)]
            CNp=[1]+[0 for j in range(len(Initial_block_RD)-1)]
            k=Af_Letter[0]
            for m in k:
                CNm[ord(m[0])-96]+=1
            k=Af_Letter[1]
            for m in k:
                CNp[ord(m[0])-96]+=1
            RDm=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[0],flank))/2]+[0 for j in range(len(Initial_block_RD)-1)]
            RDp=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[1],flank))/2]+[0 for j in range(len(Initial_block_RD)-1)]
            RDs=[RDm,RDp]
            for p in range(len(Af_Letter)):
                for q in range(len(Af_Letter[p])):
                    RDs[p][ord(Af_Letter[p][q][0])-96]+=Af_block_RD[p][q]
            for r in range(len(Initial_block_RD))[1:]:
                if CNm[r]==CNp[r]:
                    RDs[0][r]+=Initial_block_RD[r]/2
                    RDs[1][r]+=Initial_block_RD[r]/2
                elif CNm[r]==0 and not CNp[r]==0:
                    RDs[1][r]+=Initial_block_RD[r]
                elif CNp[r]==0 and not CNm[r]==0:
                    RDs[0][r]+=Initial_block_RD[r]
                else:
                    RDs[0][r]+=Initial_block_RD[r]*CNm[r]/(CNp[r]+CNm[r])
                    RDs[1][r]+=Initial_block_RD[r]*CNp[r]/(CNp[r]+CNm[r])
            CNs=[CNm,CNp]
            return [CNs,RDs]
        def All_Block_RD_2(Initial_block_RD,Af_block_RD,Af_Letter,bps,flank):
            RDs=[[],[]]
            CNs=[[],[]]
            for let in [chr(97+i) for i in range(len(bps)-1)]:
                CNs[0].append(Af_Letter[0].count(let)+Af_Letter[0].count(let+'^'))
                CNs[1].append(Af_Letter[1].count(let)+Af_Letter[1].count(let+'^'))
                if not CNs[0][-1]+CNs[1][-1]==0:
                    RDs[0].append(Initial_block_RD[ord(let)-96]*CNs[0][-1]/(CNs[0][-1]+CNs[1][-1]))
                    RDs[1].append(Initial_block_RD[ord(let)-96]*CNs[1][-1]/(CNs[0][-1]+CNs[1][-1]))
                if CNs[0][-1]+CNs[1][-1]==0:
                    RDs[0].append(0)
                    RDs[1].append(0)
            for key in  Af_block_RD[0].keys():
                if not key=='left' and not key=='right':
                    RDs[0][ord(key.split('_')[0])-97]+=float(Af_block_RD[0][key])/float(bps[ord(key.split('_')[0])-96]-bps[ord(key.split('_')[0])-97])*Window_Size
            for key in  Af_block_RD[1].keys():
                if not key=='left' and not key=='right':
                    RDs[1][ord(key.split('_')[0])-97]+=float(Af_block_RD[1][key])/float(bps[ord(key.split('_')[0])-96]-bps[ord(key.split('_')[0])-97])*Window_Size
            CNs[0]=[1]+CNs[0]
            CNs[1]=[1]+CNs[1]
            RDs[0]=[Af_block_RD[0]['left']+Initial_block_RD[0]/2]+RDs[0]
            RDs[1]=[Af_block_RD[1]['left']+Initial_block_RD[0]/2]+RDs[1]
            return [CNs,RDs]
        def basic_block_RD(Initial_block_RD,original_letters,Af_Letter):
            blockRD={}
            for i in range(len(original_letters)):
                blockRD[original_letters[i]]=Initial_block_RD[i]
            after_letter=[]
            for j in Af_Letter:
                for k in j:
                    after_letter.append(k[0])
            afterletter=sorted(after_letter)
            blockCN={}
            for m in after_letter:
                if not m in blockCN.keys():
                    blockCN[m]=1
                elif m in blockCN.keys():
                    blockCN[m]+=1
            Af_rd=[]
            for n in Af_Letter:
                Af_rd.append([])
                for o in n:
                    Af_rd[-1].append(blockRD[o[0]]/blockCN[o[0]])
            return Af_rd
        def IL_Stat_Simp(StatFile):
            fstat=open(StatFile)
            temp=fstat.readline()
            temp=fstat.readline()
            temp=fstat.readline()
            p1=fstat.readline().strip().split()
            Mean1=float(p1[1])
            STD1=float(p1[2])
            fstat.close()
            return [[Mean1,0,STD1,1,1,0],[1,Mean1,STD1]]
        def IL_Stat(StatFile):  
            fstat=open(StatFile)
            temp=fstat.readline()
            temp=fstat.readline()
            temp=fstat.readline()
            p1=fstat.readline().strip().split()
            Prop1=float(p1[0])
            Mean1=float(p1[1])
            STD1=float(p1[2])
            fstat.readline()
            p1=fstat.readline().strip().split()
            Prop2=float(p1[0])
            Mean2=float(p1[1])
            STD2=float(p1[2])
            p1=fstat.readline().strip().split()
            p1=fstat.readline().strip().split()
            Prop3=float(p1[0])
            Mean3=float(p1[1])
            STD3=float(p1[2])
            fstat.close()
            return [[Mean1,Mean2,STD1,STD2,Prop1,Prop2],[Prop3,Mean3,STD3]]
        def BPs_Coverage(Af_Letter,original_bp_list,original_letters,Letter_Through,Af_Info,flank):
            blocklen={}
            for i in range(len(original_bp_list)-1):
                blocklen[original_letters[i]]=original_bp_list[i+1]-original_bp_list[i]
            blocklen['left']=flank
            blocklen['right']=flank
            tempM=[blocklen[j[0]] for j in Af_Letter[0]]
            tempP=[blocklen[j[0]] for j in Af_Letter[1]]
            Af_BPs=[[-flank,0]+[sum(tempM[:(k+1)]) for k in range(len(tempM))],[-flank,0,]+[sum(tempP[:(k+1)]) for k in range(len(tempP))]]
            Af_BPs=[Af_BPs[0]+[Af_BPs[0][-1]+flank],Af_BPs[1]+[Af_BPs[1][-1]+flank]]
            Af_BP_Through=[[0 for i in range(len(Af_BPs[0]))],[0 for i in range(len(Af_BPs[1]))]]
            for key in Af_Info.keys():
                if Af_Info[key][6]=='M':
                    tempbps=Af_BPs[0]
                    leftMost=numpy.min([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
                    rightMost=numpy.max([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
                    for m in range(len(tempbps)-1):
                        if tempbps[m+1]>leftMost and tempbps[m]<leftMost:
                            for n in range(m,len(tempbps)-1):
                                if tempbps[n+1]>rightMost and tempbps[n]<rightMost:
                                    for p in range(m+1,n+1):
                                        if len(Af_Info[key])==7:
                                            Af_BP_Through[0][p]+=1
                                        elif len(Af_Info[key])==8:
                                            Af_BP_Through[0][p]+=float(Af_Info[key][7])
                if Af_Info[key][6]=='P':
                    tempbps=Af_BPs[1]
                    leftMost=numpy.min([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
                    rightMost=numpy.max([numpy.mean(Af_Info[key][:2]),numpy.mean(Af_Info[key][2:4])])
                    for m in range(len(tempbps)-1):
                        if tempbps[m+1]>leftMost and tempbps[m]<leftMost:
                            for n in range(m,len(tempbps)-1):
                                if tempbps[n+1]>rightMost and tempbps[n]<rightMost:
                                    for p in range(m+1,n+1):
                                        if len(Af_Info[key])==7:
                                            Af_BP_Through[1][p]+=1
                                        elif len(Af_Info[key])==8:
                                            Af_BP_Through[1][p]+=float(Af_Info[key][7])
            return  [Af_BP_Through[0][1:-1],Af_BP_Through[1][1:-1]]
        def RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chromo,GC_Hash,RD_List,Let_List):
            Coverage_af_Adj=[]
            Letters=[['left']+Let_List[0]+['right'],['left']+Let_List[1]+['right']]
            Overall_Median_Coverage=float(GC_Overall_Median_Num)
            Coverage_af_Adj=RD_List
            Theo_RD=GC_Overall_Median_Coverage[str(Chromo)]
            Theo_Var=GC_Var_Coverage[str(Chromo)]
            Theo_Std=GC_Mean_Coverage[str(Chromo)]
            Prob_out=[]
            if Let_List==[[], []]:
                for i in Initial_GCRD_Adj.keys():
                    if not i in ['left','right']:
                        Prob_out.append(standard_norm_pdf_solver(Initial_GCRD_Adj[i]*2,0,Theo_Std))
                        #Prob_out.append(Prob_Norm(Initial_GCRD_Adj[i],0,Theo_Var))
            else:
                for i in Coverage_af_Adj:
                    for j in i:
                        #Prob_out.append(Prob_Norm(j*2,Theo_RD,Theo_Std))
                        Prob_out.append(standard_norm_pdf_solver(j*2,Theo_RD,Theo_Std))
            return numpy.mean(Prob_out)
        def delete_block_produce(Letter_List):
            delete_block=[x.upper() for x in Letter_List]
            return delete_block
        def delete_BPs_produce(Letter_List):
            delete_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
            return delete_BPs
        def insert_block_produce(Letter_List,Letter_List_origin):
            Insert_Pool=[]
            for i in Letter_List_origin:
                if not i in Letter_List and not i+'^' in Letter_List:
                    Insert_Pool.append(i)
            return(Insert_Pool)
        def insert_BPs_produce(Letter_List,Letter_List_origin):
            Insert_Pool=[]
            for i in Letter_List_origin:
                if not i in Letter_List and not i+'^' in Letter_List:
                    Insert_Pool.append(i)   
            Insert_Pool2=[[(ord(j)-96),(ord(j)-95)] for j in Insert_Pool]
            return(Insert_Pool2)
        def invert_block_produce(Letter_List):
            invert_block=[x.upper() for x in Letter_List]
            return invert_block
        def invert_BPs_produce(Letter_List):
            invert_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
            return invert_BPs
        def CopyPaste_block_produce(Letter_List):
            CopyPaste_block=[x.upper() for x in Letter_List]
            return CopyPaste_block
        def CopyPaste_BPs_produce(Letter_List):
            CopyPaste_BPs=[tuple(range(len(Letter_List)+2)[1:]) for j in range(len(Letter_List)+1)]
            return CopyPaste_BPs
        def CutPaste_block_produce(Letter_List):
            CutPaste_block=[x.upper() for x in Letter_List]
            return CutPaste_block
        def CutPaste_BPs_produce(Letter_List):
            BPs=[i for i in range(len(Letter_List)+2)[1:]]
            Letter_List=[chr(96+i) for i in BPs][:-1]
            CutPaste_BPs=[BPs[0:delete_BPs_produce(Letter_List)[ord(j[0])-97][0]-1]+BPs[delete_BPs_produce(Letter_List)[ord(j[0])-97][1]:len(BPs)] for j in Letter_List]
            return CutPaste_BPs
        def Move_Choose(Move_Sample_Pool,Ploidy,Move_probs):
            import random
            random_Pool=range(10**4)
            random_Prob=[0]+Move_probs
            random_Prob=[i*10**4 for i in random_Prob]
            for i in range(len(random_Prob[1:])):
                random_Prob[i+1]+=random_Prob[i]
            Move_numA=random.choice(random_Pool)
            for i in range(len(random_Prob)-1):
                if Move_numA<random_Prob[i+1] and not Move_numA<random_Prob[i]:
                    Move_M=Move_Sample_Pool[i]
            if Ploidy==0:
                    Move_P='x'
            elif Ploidy==1:
                    Move_P='x'
            elif Ploidy==2:
                Move_numA=random.choice(random_Pool)
                for i in range(len(random_Prob)-1):
                    if Move_numA<random_Prob[i+1] and not Move_numA<random_Prob[i]:
                        Move_P=Move_Sample_Pool[i]
            return [Move_M,Move_P]
        def Move_Choice_procedure_2(Move,Letter_List,Letter_List_origin,ChrAllele):
            from random import choice
            if len(Letter_List)==0 and not Move=='insert':
                return 'ERROR!'
            else:
                if Move=='delete':
                    Move_Pool2=[]
                    for k in range(len(delete_block_produce(Letter_List))):
                        maternal_Move=[ChrAllele,'1',str(delete_BPs_produce(Letter_List)[k][0]),str(delete_BPs_produce(Letter_List)[k][1]),'del']
                        Move_Pool2.append(maternal_Move)
                elif Move=='invert':
                    Move_Pool2=[]
                    for k in range(len(invert_block_produce(Letter_List))):
                        maternal_Move=[ChrAllele,'1',str(invert_BPs_produce(Letter_List)[k][0]),str(invert_BPs_produce(Letter_List)[k][1]),'inv']
                        Move_Pool2.append(maternal_Move)
                elif Move=='insert':
                    Move_Pool2=[]
                    insert_p=random.choice(range(len(Letter_List)+2)[1:])
                    insert_b=Letter_List_origin
                    insert_bl=[[ord(insert_b1)-96,ord(insert_b1)-95] for insert_b1 in insert_b]
                    for ip in insert_bl:
                        maternal_Move=[ChrAllele,insert_p,ip[0],ip[1],'ins']
                        Move_Pool2.append(maternal_Move)
                elif Move=='x':
                    Move_Pool2=[]
            return Move_Pool2
        def BPList_Invert(BP_List,Command):
            return BP_List
        def BPList_Invert_Letter(Letter_List,Command):
            letters_origin=Letter_List
            letters_before=letters_origin[0:(int(Command[2])-1)]
            letters_invert=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
            letters_after=letters_origin[(int(Command[3])-1):]
            letters_invert_after=[]
            for j in range(len(letters_invert)):
                letters_invert_after.append(letters_invert[len(letters_invert)-1-j])
            for k in range(len(letters_invert_after)):
                if not letters_invert_after[k][-1]=='^':
                    letters_invert_after[k]+='^'
                elif letters_invert_after[k][-1]=='^':
                    letters_invert_after[k]=letters_invert_after[k][:-1]
            return letters_before+letters_invert_after+letters_after        
        def BPList_Delete(BP_List,Command):
            bp1=BP_List[int(Command[2])-1]
            bp2=BP_List[int(Command[3])-1]
            delete=BP_List[(int(Command[2])-1):int(Command[3])]
            delete_length=delete[-1]-delete[0]
            delete_after=BP_List[int(Command[3]):]
            for i in range(len(delete_after)):
                delete_after[i]=delete_after[i]-delete_length
            return BP_List[0:(int(Command[2])-1)]+[delete[0]]+delete_after  
        def BPList_Delete_Letter(Letter_List,Command):
            letters_origin=Letter_List
            letters_delete=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
            return letters_origin[0:(int(Command[2])-1)]+letters_origin[(int(Command[3]))-1:]       
        def BPList_Insert_Pool_Sequence(BP_List_origin,Ref_Sequence_Path_File,flank):
            fi=open(Ref_Sequence_Path_File)
            temp=fi.readline().strip().split()
            seq1=[]
            while True:
                seq_temp=fi.readline().rstrip()
                if not seq_temp: break
                else:
                    seq1.append(seq_temp)
            seq_ori=''.join(seq1)
            relative_List=map(lambda a:a-BP_List_origin[0]+flank,BP_List_origin)
            Letter_List_origin=[chr(x+97) for x in range(len(BP_List_origin)-1)]
            seq_letter={}
            for i in range(len(BP_List_origin)-1):
                seq_letter[Letter_List_origin[i]]=seq_ori[relative_List[i]:relative_List[i+1]]
                seq_letter[Letter_List_origin[i]+'^']=complementary(seq_ori[relative_List[i]:relative_List[i+1]])
            fi.close()
            return seq_letter
        def BPList_Insert_Pool_Length(BP_List_origin):
            Letter_List_origin=[chr(x+97) for x in range(len(BP_List_origin)-1)]
            seq_letter={}
            for i in range(len(BP_List_origin)-1):
                seq_letter[Letter_List_origin[i]]=int(BP_List_origin[i+1])-int(BP_List_origin[i])
                seq_letter[Letter_List_origin[i]+'^']=int(BP_List_origin[i+1])-int(BP_List_origin[i])
            return seq_letter
        def BPList_Insert(BP_List,Command,BP_List_origin):
            Insert_block=chr(int(Command[2])+96)
            IL=BPList_Insert_Pool_Length(BP_List_origin)[Insert_block]
            Insert_Before=BP_List[0:int(Command[1])]
            insert_diff=[x+IL for x in BP_List[(int(Command[1])-1):]]
            return Insert_Before+insert_diff
        def BPList_Insert_Letter(Letter_List,Command):
            Insert_block=chr(int(Command[2])+96)
            Letters_Before=Letter_List[0:(int(Command[1])-1)]
            Letters_After=Letter_List[(int(Command[1])-1):]
            return Letters_Before+[Insert_block]+Letters_After
        def BPList_CopyPaste(BP_List,Command):
            if int(Command[1])<=int(Command[2]):
                Block_Before=BP_List[0:int(Command[1])]
                Block_ip_bp1=BP_List[(int(Command[1])-1):int(Command[2])]
                Block_bp1_bp2=BP_List[(int(Command[2])-1):int(Command[3])]
                Block_After=BP_List[(int(Command[3])-1):]
                Block_bp1_bp2_cp=[]
                Dist1=Block_bp1_bp2[0]-Block_Before[-1]
                for i in range(len(Block_bp1_bp2)):
                    Block_bp1_bp2_cp.append(Block_bp1_bp2[i]-Dist1)
                Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
                for j in range(len(Block_ip_bp1)):
                    Block_ip_bp1[j]=Block_ip_bp1[j]+Dist2
                for j in range(len(Block_bp1_bp2)):
                    Block_bp1_bp2[j]=Block_bp1_bp2[j]+Dist2
                for j in range(len(Block_After)):
                    Block_After[j]=Block_After[j]+Dist2
                return Block_Before[:-1]+Block_bp1_bp2_cp[:-1]+Block_ip_bp1[:-1]+Block_bp1_bp2[:-1]+Block_After
            if int(Command[1])>=int(Command[3]):
                Block_Before=BP_List[0:int(Command[2])]
                Block_bp1_bp2=BP_List[(int(Command[2])-1):int(Command[3])]
                Block_bp2_ip=BP_List[(int(Command[3])-1):int(Command[1])]
                Block_After=BP_List[(int(Command[1])-1):]
                Block_bp1_bp2_cp=[]
                Dist1=Block_After[0]-Block_bp1_bp2[0]
                for i in range(len(Block_bp1_bp2)):
                    Block_bp1_bp2_cp.append(Block_bp1_bp2[i]+Dist1)
                Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
                for j in range(len(Block_After)):
                    Block_After[j]=Block_After[j]+Dist2
                return Block_Before[:-1]+Block_bp1_bp2[:-1]+Block_bp2_ip[:-1]+Block_bp1_bp2_cp[:-1]+Block_After 
        def BPList_CopyPaste_Letter(Letter_List,Command):
            letters_origin=Letter_List
            letters_before=letters_origin[0:(int(Command[1])-1)]
            letters_copy=letters_origin[(int(Command[2])-1):(int(Command[3])-1)]
            letters_after=letters_origin[(int(Command[1])-1):]
            return letters_before+letters_copy+letters_after    
        def BPList_CutPaste(BP_List,Command):
            BP_Line=Command
            if int(BP_Line[1])<int(BP_Line[2]):
                Block_Before=BP_List[0:int(BP_Line[1])]
                Block_ip_bp1=BP_List[(int(BP_Line[1])-1):int(BP_Line[2])]
                Block_bp1_bp2=BP_List[(int(BP_Line[2])-1):int(BP_Line[3])]
                Block_After=BP_List[(int(BP_Line[3])-1):]
                Dist1=Block_bp1_bp2[0]-Block_Before[-1]
                for i in range(len(Block_bp1_bp2)):
                    Block_bp1_bp2[i]=Block_bp1_bp2[i]-Dist1
                Dist2=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
                for j in range(len(Block_ip_bp1)):
                    Block_ip_bp1[j]=Block_ip_bp1[j]+Dist2
                return Block_Before[:-1]+Block_bp1_bp2[:-1]+Block_ip_bp1[:-1]+Block_After
            if int(BP_Line[1])>int(BP_Line[3]):
                Block_Before=BP_List[0:int(BP_Line[2])]
                Block_bp1_bp2=BP_List[(int(BP_Line[2])-1):int(BP_Line[3])]
                Block_bp2_ip=BP_List[(int(BP_Line[3])-1):int(BP_Line[1])]
                Block_After=BP_List[(int(BP_Line[1])-1):]
                Dist1=Block_bp1_bp2[-1]-Block_bp1_bp2[0]
                for i in range(len(Block_bp2_ip)):
                    Block_bp2_ip[i]=Block_bp2_ip[i]-Dist1
                Dist2=Block_After[0]-Block_bp1_bp2[-1]
                for j in range(len(Block_bp1_bp2)):
                    Block_bp1_bp2[j]=Block_bp1_bp2[j]+Dist2
                return Block_Before[:-1]+Block_bp2_ip[:-1]+Block_bp1_bp2[:-1]+Block_After       
        def BPList_CutPaste_Letter(Letter_List,Command):
            Letter_Line=Command
            if int(Letter_Line[1])<int(Letter_Line[2]): 
                Block_Before=Letter_List[:(numpy.min(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1)]
                Block_After=Letter_List[(numpy.max(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1):]        
                Block_cut=Letter_List[(int(Letter_Line[2])-1):(int(Letter_Line[3])-1)]
                Block_ip_cut=Letter_List[(int(Letter_Line[1])-1):(int(Letter_Line[2])-1)]
                return Block_Before+Block_cut+Block_ip_cut+Block_After
            if  int(Letter_Line[1])>int(Letter_Line[3]):
                Block_Before=Letter_List[:(numpy.min(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1)]
                Block_After=Letter_List[(numpy.max(int(Letter_Line[1]),int(Letter_Line[2]),int(Letter_Line[3]))-1):]        
                Block_cut=Letter_List[(int(Letter_Line[2])-1):(int(Letter_Line[3])-1)]
                Block_ip_cut=Letter_List[(int(Letter_Line[3])-1):(int(Letter_Line[1])-1)]
                return Block_Before+Block_ip_cut+Block_cut+Block_After  
        def BPList_X(BP_List,Command):
            return BP_List
        def BPList_X_Letter(Letter_List,Command):
            return Letter_List
        def BPList_Rearrange(BP_List,Command,BP_List_origin):
            if Command[-1]=='del' or Command[-1]=='delete':
                return BPList_Delete(BP_List,Command)
            elif Command[-1]=='inv' or Command[-1]=='invert':
                return BPList_Invert(BP_List,Command)   
            elif Command[-1]=='ins' or Command[-1]=='insert':
                return BPList_Insert(BP_List,Command,BP_List_origin)
            elif Command[-1]=='copy+paste' or Command[-1]=='CopyPaste':
                return BPList_CopyPaste(BP_List,Command)
            elif Command[-1]=='cut+paste' or Command[-1]=='CutPaste':
                return BPList_CutPaste(BP_List,Command)
            elif Command[-1]=='x' or Command[-1]=='X':
                return BPList_X(BP_List,Command)
        def LetterList_Rearrange(Letter_List,Command,BP_List_origin):
            if Command[-1]=='del' or Command[-1]=='delete':
                return BPList_Delete_Letter(Letter_List,Command)
            elif Command[-1]=='inv' or Command[-1]=='invert':
                return BPList_Invert_Letter(Letter_List,Command)    
            elif Command[-1]=='ins' or Command[-1]=='insert':
                return BPList_Insert_Letter(Letter_List,Command)
            elif Command[-1]=='copy+paste' or Command[-1]=='CopyPaste':
                return BPList_CopyPaste_Letter(Letter_List,Command)
            elif Command[-1]=='cut+paste' or Command[-1]=='CutPaste':
                return BPList_CutPaste_Letter(Letter_List,Command)
            elif Command[-1]=='x' or Command[-1]=='X':
                return BPList_X_Letter(Letter_List,Command)
        def Rearrange_RefSeq(Ref_Sequence_Path_File,BP_List_origin,Letter_List_After,flank):
            Seq_pool=BPList_Insert_Pool_Sequence(BP_List_origin,Ref_Sequence_Path_File,flank)
            Seq_Out=[Seq_pool[x] for x in Letter_List_After]
            return(''.join(Seq_Out))
        def Sequence_Rearrange_3(Insert_Seq,Af_Letter,original_bp_list,flank):
            seql=[]
            seql.append(Insert_Seq['left'])
            for let in Af_Letter:
                seql.append(Insert_Seq['I'+let])
            seql.append(Insert_Seq['right'])
            return ''.join(seql)
        def MP_Letter_Through_Read_Rearrange_3(Full_Info_of_Read,Af_Letter,original_bp):
            Full_Info_of_Read=[int(i) for i in Full_Info_of_Read[:4]]+Full_Info_of_Read[4:]
            Left_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
            Right_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
            if numpy.max(Left_Read[:2])>numpy.max(Right_Read[:2]):
                Left_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
                Right_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
            Original_Letter=[chr(i+97) for i in range(len(original_bp)-1)]
            Letter_Len={}
            for i in range(len(Original_Letter)):
                Letter_Len[Original_Letter[i]]=original_bp[i+1]-original_bp[i]
            bp_MP=[[0],[0]]
            for j in range(len(Af_Letter)):
                for k in Af_Letter[j]:
                    if k=='le' or k=='ri':
                        bp_MP[j].append(bp_MP[j][-1]+Letter_Len[k])
                    else:
                        bp_MP[j].append(bp_MP[j][-1]+Letter_Len[k[0]])
            Af_M_let=['left']+Af_Letter[0]+['right']
            Af_M_bpp=[-flank]+bp_MP[0]+[bp_MP[0][-1]+flank]
            Af_P_let=['left']+Af_Letter[1]+['right']
            Af_P_bpp=[-flank]+bp_MP[1]+[bp_MP[1][-1]+flank]
            LeftM_New=[]
            RightM_New=[]
            for i in range(len(Af_M_let)):
                if Left_Read[3][0]==Af_M_let[i]:
                    LeftM_New.append([Af_M_bpp[i]+Left_Read[3][1],Af_M_bpp[i]+Left_Read[3][2],Left_Read[2]])
                elif Left_Read[3][0]+'^'==Af_M_let[i]:
                    LeftM_New.append([Af_M_bpp[i+1]-Left_Read[3][2],Af_M_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
                elif Left_Read[3][0]==Af_M_let[i]+'^':
                    LeftM_New.append([Af_M_bpp[i+1]-Left_Read[3][2],Af_M_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
                if Right_Read[3][0]==Af_M_let[i]:
                    RightM_New.append([Af_M_bpp[i]+Right_Read[3][1],Af_M_bpp[i]+Right_Read[3][2],Right_Read[2]])
                elif Right_Read[3][0]+'^'==Af_M_let[i]:
                    RightM_New.append([Af_M_bpp[i+1]-Right_Read[3][2],Af_M_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
                elif Right_Read[3][0]==Af_M_let[i]+'^':
                    RightM_New.append([Af_M_bpp[i+1]-Right_Read[3][2],Af_M_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
            LeftP_New=[]
            RightP_New=[]
            for i in range(len(Af_P_let)):
                if Left_Read[3][0]==Af_P_let[i]:
                    LeftP_New.append([Af_P_bpp[i]+Left_Read[3][1],Af_P_bpp[i]+Left_Read[3][2],Left_Read[2]])
                elif Left_Read[3][0]+'^'==Af_P_let[i]:
                    LeftP_New.append([Af_P_bpp[i+1]-Left_Read[3][2],Af_P_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
                elif Left_Read[3][0]==Af_P_let[i]+'^':
                    LeftP_New.append([Af_P_bpp[i+1]-Left_Read[3][2],Af_P_bpp[i+1]-Left_Read[3][1],oppo_direct(Left_Read[2])])
                if Right_Read[3][0]==Af_P_let[i]:
                    RightP_New.append([Af_P_bpp[i]+Right_Read[3][1],Af_P_bpp[i]+Right_Read[3][2],Right_Read[2]])
                elif Right_Read[3][0]+'^'==Af_P_let[i]:
                    RightP_New.append([Af_P_bpp[i+1]-Right_Read[3][2],Af_P_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
                elif Right_Read[3][0]==Af_P_let[i]+'^':
                    RightP_New.append([Af_P_bpp[i+1]-Right_Read[3][2],Af_P_bpp[i+1]-Right_Read[3][1],oppo_direct(Right_Read[2])])
            IN_Read_New=[]
            OUT_Read_New=[]
            for l in LeftM_New:
                for r in RightM_New:
                    if not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])<Cut_Lower and not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])>Cut_Upper and not l[2]==r[2]:
                        IN_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['M'])
                    else:
                        OUT_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['M'])
            for l in LeftP_New:
                for r in RightP_New:
                    if not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])<Cut_Lower and not numpy.max(l[:2]+r[:2])-numpy.min(l[:2]+r[:2])>Cut_Upper and not l[2]==r[2]:
                        IN_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['P'])
                    else:
                        OUT_Read_New.append(l[:2]+r[:2]+[l[2]]+[r[2]]+['P'])
            if not len(IN_Read_New)==0:
                return ['IN']+IN_Read_New
            elif len(IN_Read_New)==0 and not len(OUT_Read_New)==0:
                return ['OUT']+OUT_Read_New
            else:
                return ['NON',[-flank/2,-flank/2,-flank/2,-flank/2,'+','-','M'],[-flank/2,-flank/2,-flank/2,-flank/2,'+','-','P']]
        def MP_Letter_Through_Read_Rearrange_4(Full_Info_of_Read,Letter_Len,Af_BP_rela,Af_Letter_rela,original_bp):
            Full_Info_of_Read=[int(i) for i in Full_Info_of_Read[:4]]+Full_Info_of_Read[4:]
            Left_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
            Right_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
            if numpy.max(Left_Read[:2])>numpy.max(Right_Read[:2]):
                Left_Read=Full_Info_of_Read[2:4]+[Full_Info_of_Read[5],Full_Info_of_Read[9:12]]
                Right_Read=Full_Info_of_Read[:2]+[Full_Info_of_Read[4],Full_Info_of_Read[6:9]]
            Left_Reads=[]
            Right_Reads=[]
            Left_Info=Left_Read[3]
            Right_Info=Right_Read[3]
            if Left_Read[3][0]=='left':
                Left_Reads.append(Left_Read[:3]+['M',Af_Letter_rela[0][0]])
                Left_Reads.append(Left_Read[:3]+['P',Af_Letter_rela[1][0]])
            if Left_Read[3][0]=='right':
                Left_Reads.append([Af_BP_rela[0][-2]+j for j in Left_Read[3][1:]]+[Left_Read[2]]+['M',Af_Letter_rela[0][-1]])
                Left_Reads.append([Af_BP_rela[1][-2]+j for j in Left_Read[3][1:]]+[Left_Read[2]]+['P',Af_Letter_rela[1][-1]])
            if not Left_Read[3][0] in ['left','right']:
                for i in range(len(Af_Letter_rela)):
                    for j in range(len(Af_Letter_rela[i]))[1:-1]:
                        if Af_Letter_rela[i][j].split('_')[0][0]==Left_Read[3][0]:
                            if not Af_Letter_rela[i][j].split('_')[0][-1]=='^':
                                Left_Reads.append([Af_BP_rela[i][j]+k for k in Left_Read[3][1:]]+[Left_Read[2],['M','P'][i],Af_Letter_rela[i][j]])
                            if Af_Letter_rela[i][j].split('_')[0][-1]=='^':
                                Left_Reads.append(sorted([Af_BP_rela[i][j+1]-k for k in Left_Read[3][1:]])+[oppo_direct(Left_Read[2]),['M','P'][i],Af_Letter_rela[i][j]])
            if Right_Read[3][0]=='left':
                Right_Reads.append(Right_Read[:3]+['M',Af_Letter_rela[0][0]])
                Right_Reads.append(Right_Read[:3]+['P',Af_Letter_rela[1][0]])
            if Right_Read[3][0]=='right':
                Right_Reads.append([Af_BP_rela[0][-2]+j for j in Right_Read[3][1:]]+[Right_Read[2]]+['M',Af_Letter_rela[0][-1]])
                Right_Reads.append([Af_BP_rela[1][-2]+j for j in Right_Read[3][1:]]+[Right_Read[2]]+['P',Af_Letter_rela[1][-1]])
            if not Right_Read[3][0] in ['left','right']:
                for i in range(len(Af_Letter_rela)):
                    for j in range(len(Af_Letter_rela[i]))[1:-1]:
                        if Af_Letter_rela[i][j].split('_')[0][0]==Right_Read[3][0]:
                            if not Af_Letter_rela[i][j].split('_')[0][-1]=='^':
                                Right_Reads.append([Af_BP_rela[i][j]+k for k in Right_Read[3][1:]]+[Right_Read[2],['M','P'][i],Af_Letter_rela[i][j]])
                            if Af_Letter_rela[i][j].split('_')[0][-1]=='^':
                                Right_Reads.append(sorted([Af_BP_rela[i][j+1]-k for k in Right_Read[3][1:]])+[oppo_direct(Right_Read[2]),['M','P'][i],Af_Letter_rela[i][j]])
            In_Pair_Reads=['IN']
            Out_Pair_Reads=['OUT']
            for m in Left_Reads:
                for n in Right_Reads:
                    if m[3]==n[3]:
                        Pair_Read=m[:2]+n[:2]+[m[2],n[2],m[3],m[4],n[4]]
                        if (numpy.max(Pair_Read[:4])-numpy.min(Pair_Read[:4]))>Cut_Lower and (numpy.max(Pair_Read[:4])-numpy.min(Pair_Read[:4]))<Cut_Upper and not Pair_Read[4]==Pair_Read[5]:
                            In_Pair_Reads.append(Pair_Read) 
                        else:   
                            Out_Pair_Reads.append(Pair_Read)    
            output=[]
            if not len(In_Pair_Reads)==1:
                for i in In_Pair_Reads[1:]:
                    output.append(i+[Left_Info[0],Left_Info[1]/Window_Size,float(Left_Info[1]/Window_Size*Window_Size+Window_Size-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/Window_Size,float(Left_Info[2]-Left_Info[2]/Window_Size*Window_Size)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/Window_Size,float(Right_Info[1]/Window_Size*Window_Size+Window_Size-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/Window_Size,float(Right_Info[2]-Right_Info[2]/Window_Size*Window_Size)/float(Right_Info[2]-Right_Info[1])])
            else:
                if not len(Out_Pair_Reads)==1:
                    output1=[j for j in Out_Pair_Reads[1:] if not j[4]==j[5]]
                    output2=[j for j in Out_Pair_Reads[1:] if j[4]==j[5]]
                    if not len(output1)==0:
                        pdfs=[standard_pdf_IL_calculate(numpy.max(output1[i][:4])-numpy.min(output1[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(output1))]
                        for j in range(len(pdfs)):
                            if pdfs[j]==numpy.max(pdfs):
                                output.append(output1[j]+[Left_Info[0],Left_Info[1]/Window_Size,float(Left_Info[1]/Window_Size*Window_Size+Window_Size-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/Window_Size,float(Left_Info[2]-Left_Info[2]/Window_Size*Window_Size)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/Window_Size,float(Right_Info[1]/Window_Size*Window_Size+Window_Size-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/Window_Size,float(Right_Info[2]-Right_Info[2]/Window_Size*Window_Size)/float(Right_Info[2]-Right_Info[1])])
                    else:
                        pdfs=[standard_pdf_IL_calculate(numpy.max(output2[i][:4])-numpy.min(output2[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(output2))]
                        for j in range(len(pdfs)):
                            if pdfs[j]==numpy.max(pdfs):
                                output.append(output2[j]+[Left_Info[0],Left_Info[1]/Window_Size,float(Left_Info[1]/Window_Size*Window_Size+Window_Size-Left_Info[1])/float(Left_Info[2]-Left_Info[1]),Left_Info[2]/Window_Size,float(Left_Info[2]-Left_Info[2]/Window_Size*Window_Size)/float(Left_Info[2]-Left_Info[1])]+[Right_Info[0],Right_Info[1]/Window_Size,float(Right_Info[1]/Window_Size*Window_Size+Window_Size-Right_Info[1])/float(Right_Info[2]-Right_Info[1]),Right_Info[2]/Window_Size,float(Right_Info[2]-Right_Info[2]/Window_Size*Window_Size)/float(Right_Info[2]-Right_Info[1])])
                else:
                    output.append([-flank/2,-flank/2,-flank/2,-flank/2,'+','-','M','left_0','left_0','left',1,1,0,0,'left',1,1,0,0])
                    output.append([-flank/2,-flank/2,-flank/2,-flank/2,'+','-','P','left_0','left_0','left',1,1,0,0,'left',1,1,0,0])
            return  output
        def oppo_direct(x):
            if x=='+':
                return '-'
            elif x=='-':
                return '+'
        def Letter_Through_Rearrange_3(IL_Statistics,Letter_Through,Af_Letter,original_bp_list):
            New_Info_of_Reads={}
            for key in Letter_Through.keys():
                after_key=MP_Letter_Through_Read_Rearrange_3(Letter_Through[key],Af_Letter,original_bp_list)
                if len(after_key)==2:
                    New_Info_of_Reads[key]=after_key[1]
                elif len(after_key)>2:
                    if after_key[0]=='IN' or after_key[0]=='NON':
                        for i in range(len(after_key))[1:]:
                            New_Info_of_Reads[key+'_'+str(i)]=after_key[i]+[str(float(1)/float(len(after_key)-1))]
                    elif after_key[0]=='OUT':
                        pdfs=[standard_pdf_IL_calculate(numpy.max(after_key[i][:4])-numpy.min(after_key[i][:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for i in range(len(after_key))[1:]]
                        best_Out=[after_key[j+1] for j in range(len(pdfs)) if pdfs[j]==numpy.max(pdfs)]
                        if len(best_Out)==1:
                            New_Info_of_Reads[key]=best_Out[0]
                        elif len(best_Out)>1:
                            for i in range(len(best_Out)+1)[1:]:
                                New_Info_of_Reads[key+'_'+str(i)]=best_Out[i-1]+[str(float(1)/float(len(best_Out)))]
            return New_Info_of_Reads
        def complement(direction):
            if direction=='+':
                return '-'
            elif direction=='-':
                return '+'
            else:
                return 'error'
        min_resolution=70
        def candidate_QC_Control(Read_List):
            if Read_List==[]:
                return []
            else:
                Qual_Filter_1=[]
                for j in Read_List:
                    if not j[1]-j[0]>ReadLength+min_resolution and j[1]-j[0]>0 and not j[3]-j[2]>ReadLength+min_resolution and j[3]-j[2]>0:
                        Qual_Filter_1.append(j)
                if not Qual_Filter_1==[]:
                    if len(Qual_Filter_1)==1:
                        Qual_Filter_1[0]+=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_1]
                        return Qual_Filter_1
                    else:
                        Qual_Filter_2=[]
                        for j2 in Qual_Filter_1:
                            if j2[-2:]==['+','-']:
                                  Qual_Filter_2.append(j2)
                        if not Qual_Filter_2==[]:
                            if len(Qual_Filter_2)==1:
                                Qual_Filter_2[0]+=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                return Qual_Filter_2
                            else:
                                Qual_Filter_3=[]
                                Qual_IL=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                for jq in range(len(Qual_IL)):
                                    if Qual_IL[jq]==max(Qual_IL) and not Qual_Filter_1[jq] in Qual_Filter_3:
                                        Qual_Filter_3.append(Qual_Filter_1[jq]+[max(Qual_IL)])
                                return Qual_Filter_3                    
                        else:
                            Qual_Filter_2=Qual_Filter_1
                            if len(Qual_Filter_2)==1:
                                Qual_Filter_2[0]+=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                return Qual_Filter_2
                            else:
                                Qual_Filter_3=[]
                                Qual_IL=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                for jq in range(len(Qual_IL)):
                                    if Qual_IL[jq]==max(Qual_IL) and not Qual_Filter_1[jq] in Qual_Filter_3:
                                        Qual_Filter_3.append(Qual_Filter_1[jq]+[max(Qual_IL)])
                                return Qual_Filter_3 
                else:
                    return []
        def candidate_QC_Control2(M_Read_List,P_Read_List):
            Qual_Filter_1=[]
            for i in M_Read_List:
                Qual_Filter_1.append(i+['m'])
            for i in P_Read_List:
                Qual_Filter_1.append(i+['p'])
            Qual_Filter_2=[]
            for i in Qual_Filter_1:
                if i[-4:-2]==['+','-']:
                    Qual_Filter_2.append(i)
            if not Qual_Filter_2==[]:
                Qual_Filter_3=[]
                IL_Qual=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                for j in range(len(IL_Qual)):
                    if IL_Qual[j]==max(IL_Qual) and not Qual_Filter_2[j] in Qual_Filter_3:
                        Qual_Filter_3.append(Qual_Filter_2[j])
            else:
                Qual_Filter_2=Qual_Filter_1
                Qual_Filter_3=[]
                IL_Qual=[standard_pdf_IL_calculate(max(j3[:4])-min(j3[:4]),IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                for j in range(len(IL_Qual)):
                    if IL_Qual[j]==max(IL_Qual) and not Qual_Filter_2[j] in Qual_Filter_3:
                        Qual_Filter_3.append(Qual_Filter_2[j])
            return Qual_Filter_3
        def Cov_Cal_Block(pos,bp,cov,perc):
            for j in range(len(bp)-2):
                if not pos[0]<bp[j] and pos[0]<bp[j+1]:
                    if not pos[1]<bp[j] and pos[1]<bp[j+1]:
                        cov[j]+=(pos[1]-pos[0])*perc
                    elif not pos[1]<bp[j+1] and pos[1]<bp[j+2]:
                        cov[j]+=(bp[j+1]-pos[0])*perc
                        cov[j+1]+=(pos[1]-bp[j+1])*perc
                    elif not pos[1]<temp_bp[0][j+2] and pos[1]<temp_bp[0][j+3]:
                        cov[j]+=(bp[j+1]-pos[0])*perc
                        cov[j+1]+=(bp[j+2]-bp[j+1])*perc
                        cov[j+2]+=(pos[1]-bp[j+2])*perc
            j=len(bp)-2
            if not pos[0]<bp[j] and pos[0]<bp[j+1]:
                if not pos[1]<bp[j] and pos[1]<bp[j+1]:
                    cov[j]+=(pos[1]-pos[0])*perc
                else:
                    cov[j]+=(bp[j+1]-pos[0])*perc
        def penal_calculate(Map_All,temp_bp, Af_Letter,Af_BP,RD_within_B,letters_numbers,BlockGC,NoMapPenal):
            out_rd=[[0 for i in temp_bp[0][:-1]],[0 for i in temp_bp[1][:-1]]]
            IL_Rec={}
            DR_Penal=0
            out_tb=[[0 for i in temp_bp[0]],[0 for i in temp_bp[1]]]
            IL_All=[]
            for i in Map_All:
                if len(i)>4:
                    IL_All.append(i[6])
            IL_All_Mean=numpy.mean(IL_All)
            for i in Map_All:
                if len(i)>4:
                    if not i[6] in IL_Rec.keys() and not abs(i[6]/IL_All_Mean)>10:
                        IL_Rec[i[6]]=i[8]
                    elif i[6] in IL_Rec.keys() and not abs(i[6]/IL_All_Mean)>10:
                        IL_Rec[i[6]]+=i[8]
                    if not i[4:6]==['+','-']:
                        DR_Penal+=1
                    if i[7]=='m':
                        i_block=[]
                        for k in i[:4]:
                            if k<temp_bp[0][1]:
                                i_block.append(0)
                            elif k>temp_bp[0][-2]-1:
                                i_block.append(len(temp_bp[0])-2)
                            else:
                                for j in range(len(temp_bp[0])-1)[1:-1]:
                                    if temp_bp[0][j]-1<k and temp_bp[0][j+1]>k:
                                        i_block.append(j)
                        if i_block[0]==i_block[1] and i_block[2]==i_block[3]:
                            out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
                            out_rd[0][i_block[2]]+=(i[3]-i[2])*i[-1]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[0][k2]+=i[8]/(i_block[2]-i_block[1])
                        elif not i_block[0]==i_block[1] and i_block[2]==i_block[3]:
                            out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
                            out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
                            out_rd[0][i_block[2]]+=(i[3]-i[2])*i[-1]
                            out_tb[0][i_block[1]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[0][k2]+=i[8]
                        elif i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
                            out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
                            out_rd[0][i_block[2]]+=(temp_bp[0][i_block[2]+1]-i[2])*i[-1]
                            out_rd[0][i_block[3]]+=(i[3]-temp_bp[0][i_block[3]])*i[-1]
                            out_tb[0][i_block[3]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[0][k2]+=i[8]
                        elif not i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
                            out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
                            out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
                            out_rd[0][i_block[2]]+=(temp_bp[0][i_block[2]+1]-i[2])*i[-1]
                            out_rd[0][i_block[3]]+=(i[3]-temp_bp[0][i_block[3]])*i[-1]
                            out_tb[0][i_block[1]]+=i[8]
                            out_tb[0][i_block[3]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[0][k2]+=i[8]
                    if i[7]=='p':
                        i_block=[]
                        for k in i[:4]:
                            if k<temp_bp[1][1]:
                                i_block.append(0)
                            elif k>temp_bp[1][-2]-1:
                                i_block.append(len(temp_bp[1])-2)
                            else:
                                for j in range(len(temp_bp[1])-1)[1:-1]:
                                    if temp_bp[1][j]-1<k and temp_bp[1][j+1]>k:
                                        i_block.append(j)
                        if i_block[0]==i_block[1] and i_block[2]==i_block[3]:
                            out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
                            out_rd[1][i_block[2]]+=(i[3]-i[2])*i[-1]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[1][k2]+=i[8]
                        elif not i_block[0]==i_block[1] and i_block[2]==i_block[3]:
                            out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
                            out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
                            out_rd[1][i_block[2]]+=(i[3]-i[2])*i[-1]
                            out_tb[1][i_block[1]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[1][k2]+=i[8]
                        elif i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
                            out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
                            out_rd[1][i_block[2]]+=(temp_bp[1][i_block[2]+1]-i[2])*i[-1]
                            out_rd[1][i_block[3]]+=(i[3]-temp_bp[1][i_block[3]])*i[-1]
                            out_tb[1][i_block[3]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[1][k2]+=i[8]
                        elif not i_block[0]==i_block[1] and not i_block[2]==i_block[3]:
                            out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
                            out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
                            out_rd[1][i_block[2]]+=(temp_bp[1][i_block[2]+1]-i[2])*i[-1]
                            out_rd[1][i_block[3]]+=(i[3]-temp_bp[1][i_block[3]])*i[-1]
                            out_tb[1][i_block[1]]+=i[8]
                            out_tb[1][i_block[3]]+=i[8]
                            if i[4:6]==['+', '-'] and i[6]>Penalty_For_InsertLengthZero:
                                for k2 in range(i_block[1]+1,i_block[2]+1):
                                    out_tb[1][k2]+=i[8]
                else:
                    if i[2]=='m':
                        i_block=[]
                        for k in i[:2]:
                            if k<temp_bp[0][1]:
                                i_block.append(0)
                            elif k>temp_bp[0][-2]-1:
                                i_block.append(len(temp_bp[0])-2)
                            else:
                                for j in range(len(temp_bp[0])-1)[1:-1]:
                                    if temp_bp[0][j]-1<k and temp_bp[0][j+1]>k:
                                        i_block.append(j)
                        if i_block[0]==i_block[1]:
                            out_rd[0][i_block[0]]+=(i[1]-i[0])*i[-1]
                        elif not i_block[0]==i_block[1]:
                            out_rd[0][i_block[0]]+=(temp_bp[0][i_block[0]+1]-i[0])*i[-1]
                            out_rd[0][i_block[1]]+=(i[1]-temp_bp[0][i_block[1]])*i[-1]
                    if i[2]=='p':
                        i_block=[]
                        for k in i[:2]:
                            if k<temp_bp[1][1]:
                                i_block.append(0)
                            elif k>temp_bp[1][-2]-1:
                                i_block.append(len(temp_bp[1])-2)
                            else:
                                for j in range(len(temp_bp[1])-1)[1:-1]:
                                    if temp_bp[1][j]-1<k and temp_bp[1][j+1]>k:
                                        i_block.append(j)
                        if i_block[0]==i_block[1]:
                            out_rd[1][i_block[0]]+=(i[1]-i[0])*i[-1]
                        elif not i_block[0]==i_block[1]:
                            out_rd[1][i_block[0]]+=(temp_bp[1][i_block[0]+1]-i[0])*i[-1]
                            out_rd[1][i_block[1]]+=(i[1]-temp_bp[1][i_block[1]])*i[-1]
            block_bps_chr={}
            block_bps_chr['m']={}
            block_bps_chr['p']={}
            if not Penalty_For_InsertLengthZero in IL_Rec.keys():
                IL_Rec[Penalty_For_InsertLengthZero]=NoMapPenal
            else:
                IL_Rec[Penalty_For_InsertLengthZero]+=NoMapPenal
            IL_Penal=0
            IL_Weight=0
            for i in IL_Rec.keys():
                IL_Penal+=i*IL_Rec[i]
                IL_Weight+=IL_Rec[i]
            if not IL_Weight==0:
                IL_Output=IL_Penal/IL_Weight
            else:
                IL_Output=0
            Num_Read_TB=[out_tb[0][1:-1],out_tb[1][1:-1]]
            TB_Pena_2_out=0
            Num_total_TB=[]
            for x in Num_Read_TB:
                Num_total_TB+=x
            if numpy.sum(Num_total_TB)>0:
                pvalue=scipy.stats.chisquare(Num_total_TB)[1]
            else:
                pvalue=0.0
            if pvalue>0:
                TB_Pena_2_out=numpy.log(pvalue)
            else:
                TB_Pena_2_out=-100000000
            Af_Block_Len=[[flank]+[Af_BP[0][i+1]-Af_BP[0][i] for i in range(len(Af_BP[0])-1)]+[flank],[flank]+[Af_BP[1][i+1]-Af_BP[1][i] for i in range(len(Af_BP[1])-1)]+[flank]]
            out_rd=[[out_rd[0][i]/Af_Block_Len[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]/Af_Block_Len[1][i] for i in range(len(out_rd[1]))]]
            out_rd_new=[[(RD_within_B['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[0][0],
            (RD_within_B['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[0][-1]],
            [(RD_within_B['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[1][0],
            (RD_within_B['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[1][-1]]]
            out_rd=[[out_rd_new[0][0]]+out_rd[0][1:-1]+[out_rd_new[0][-1]],[out_rd_new[1][0]]+out_rd[1][1:-1]+[out_rd_new[1][-1]]]
            out_rd_within=[[RD_within_B[Af_Letter[0][i]]/letters_numbers[0][i] for i in range(len(Af_Letter[0]))],[RD_within_B[Af_Letter[1][i]]/letters_numbers[1][i] for i in range(len(Af_Letter[1]))]]
            out_rd_within[0]=[0]+out_rd_within[0]+[0]
            out_rd_within[1]=[0]+out_rd_within[1]+[0]   
            cov_bp2=[[out_rd[0][i]+out_rd_within[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]+out_rd_within[1][i] for i in range(len(out_rd[1]))]]
            Cov_GC=[[BlockGC[k] for k in Af_Letter[0]],[BlockGC[k] for k in Af_Letter[1]]]
            adj_cov_bp=[GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrom_N,Cov_GC[0],cov_bp2[0][1:-1]),GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,chrom_N,Cov_GC[1],cov_bp2[1][1:-1])]
            return [IL_Output,adj_cov_bp,DR_Penal,TB_Pena_2_out,Num_total_TB]
        def TB_Cal_BP(pos,bp,tblist,perc):
            for i in range(len(bp)-2)[1:]:
                if not pos[0]<bp[i] and pos[1]<bp[i+1]:
                    break
                elif not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+1] and pos[1]<bp[i+2]:
                    tblist[i+1]+=perc
                elif not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+2] and pos[1]<bp[i+3]:
                    tblist[i+1]+=perc
                    tblist[i+2]+=perc
            i=len(bp)-2
            if not pos[0]<bp[i] and pos[0]<bp[i+1] and not pos[1]<bp[i+1]:
                tblist[i+1]+=perc
        def Be_Info_1_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
            be_info_1=Be_Info[0]
            for j in be_info_1:
                jMapPenam=0
                j_m_new=[]
                if j[0] in temp_letter[0] and j[3] in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]]:
                        for kb in Let_BP_Info['m'][j[3]]:
                            j_m_temp=[j[1]+ka[0],j[2]+ka[0],j[4]+kb[0],j[5]+kb[0]]
                            if j_m_temp[0]>j_m_temp[2]:
                                j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[j[-1],j[-2]]
                            else:
                                j_m_temp+=[j[-2],j[-1]]
                            j_m_new.append(j_m_temp)
                if j[0]+'^' in temp_letter[0] and j[3] in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]+'^']:
                        for kb in Let_BP_Info['m'][j[3]]:
                            j_m_temp=[ka[1]-j[2],ka[1]-j[1],kb[0]+j[4],kb[0]+j[5]]
                            if j_m_temp[0]>j_m_temp[2]:
                                j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[j[-1],complement(j[-2])]
                            else:
                                j_m_temp+=[complement(j[-2]),j[-1]]
                            j_m_new.append(j_m_temp)
                if j[0] in temp_letter[0] and j[3]+'^' in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]]:
                        for kb in  Let_BP_Info['m'][j[3]+'^']:
                            j_m_temp=[j[1]+ka[0],j[2]+ka[0],kb[1]-j[5],kb[1]-j[4]]
                            if j_m_temp[0]>j_m_temp[2]:
                                j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[complement(j[-1]),j[-2]]
                            else:
                                j_m_temp+=[j[-2],complement(j[-1])]
                            j_m_new.append(j_m_temp)
                if j[0]+'^' in temp_letter[0] and j[3]+'^' in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]+'^']:
                        for kb in  Let_BP_Info['m'][j[3]+'^']:
                            j_m_temp=[ka[1]-j[2],ka[1]-j[1],kb[1]-j[5],kb[1]-j[4]]
                            if j_m_temp[0]>j_m_temp[2]:
                                j_m_temp=j_m_temp[2:4]+j_m_temp[:2]+[complement(j[-1]),complement(j[-2])]
                            else:
                                j_m_temp+=[complement(j[-2]),complement(j[-1])]
                            j_m_new.append(j_m_temp)
                j_m_3a=candidate_QC_Control(j_m_new)
                if j_m_3a==[]:
                    jMapPenam+=1
                j_p_new=[]
                if j[0] in temp_letter[1] and j[3] in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]]:
                        for kb in Let_BP_Info['p'][j[3]]:
                            j_p_temp=[j[1]+ka[0],j[2]+ka[0],j[4]+kb[0],j[5]+kb[0]]
                            if j_p_temp[0]>j_p_temp[2]:
                                j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[j[-1],j[-2]]
                            else:
                                j_p_temp+=[j[-2],j[-1]]
                            j_p_new.append(j_p_temp)
                if j[0]+'^' in temp_letter[1] and j[3] in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]+'^']:
                        for kb in Let_BP_Info['p'][j[3]]:
                            j_p_temp=[ka[1]-j[2],ka[1]-j[1],kb[0]+j[4],kb[0]+j[5]]
                            if j_p_temp[0]>j_p_temp[2]:
                                j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[j[-1],complement(j[-2])]
                            else:
                                j_p_temp+=[complement(j[-2]),j[-1]]
                            j_p_new.append(j_p_temp)
                if j[0] in temp_letter[1] and j[3]+'^' in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]]:
                        for kb in Let_BP_Info['p'][j[3]+'^']:
                            j_p_temp=[j[1]+ka[0],j[2]+ka[0],kb[1]-j[5],kb[1]-j[4]]
                            if j_p_temp[0]>j_p_temp[2]:
                                j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[complement(j[-1]),j[-2]]
                            else:
                                j_p_temp+=[j[-2],complement(j[-1])]
                            j_p_new.append(j_p_temp)
                if j[0]+'^' in temp_letter[1] and j[3]+'^' in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]+'^']:
                        for kb in  Let_BP_Info['p'][j[3]+'^']:
                            j_p_temp=[ka[1]-j[2],ka[1]-j[1],kb[1]-j[5],kb[1]-j[4]]
                            if j_p_temp[0]>j_p_temp[2]:
                                j_p_temp=j_p_temp[2:4]+j_p_temp[:2]+[complement(j[-1]),complement(j[-2])]
                            else:
                                j_p_temp+=[complement(j[-2]),complement(j[-1])]
                            j_p_new.append(j_p_temp)
                j_p_3a=candidate_QC_Control(j_p_new)
                if j_p_3a==[]:
                    jMapPenam+=1
                if jMapPenam==2:
                    Total_Cov_For_Pen[j[0]]+=j[2]-j[1]
                    Total_Cov_For_Pen[j[3]]+=j[5]-j[4]
                    NoMapPenal+=2
                elif jMapPenam==1:
                    if j_m_3a==[]:
                        Map_P+=[jp3+['p']+[float(1)/float(len(j_p_3a))] for jp3 in j_p_3a]
                    elif j_p_3a==[]:
                        Map_M+=[jp3+['m']+[float(1)/float(len(j_m_3a))] for jp3 in j_m_3a]
                else:
                    j_mp_4a=candidate_QC_Control2(j_m_3a,j_p_3a)
                    if not j_mp_4a==[]:
                        Map_Both+=[j4+[float(1)/float(len(j_mp_4a))] for j4 in j_mp_4a]
                    else:
                        Total_Cov_For_Pen[j[0]]+=j[2]-j[1]
                        Total_Cov_For_Pen[j[3]]+=j[5]-j[4]
                        NoMapPenal+=2
            return NoMapPenal
        def Be_Info_2_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
            be_info_2=Be_Info[1]
            for j2 in be_info_2:
                j2_m_new=[]
                j2_p_new=[]
                for j in j2:
                    jMapPenam=0
                    j_m_new=[]
                    if j[0] in temp_letter[0] and j[2] in temp_letter[0] and j[4] in temp_letter[0] and j[6] in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]]:
                            for kb in Let_BP_Info['m'][j[2]]:
                                for kc in Let_BP_Info['m'][j[4]]:
                                    for kd in Let_BP_Info['m'][j[6]]:
                                        j_info_new=[ka[0]+j[1],kb[0]+j[3],kc[0]+j[5],kd[0]+j[7]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],j[-2]])
                                        else:
                                            j_m_new.append(j_info_new+[j[-2],j[-1]])
                    if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0] and j[4] in temp_letter[0] and j[6] in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]+'^']:
                            for kb in Let_BP_Info['m'][j[2]+'^']:
                                for kc in Let_BP_Info['m'][j[4]]:
                                    for kd in Let_BP_Info['m'][j[6]]:
                                        j_info_new=[kb[1]-j[3],ka[1]-j[1],kc[0]+j[5],kd[0]+j[7]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],complement(j[-2])])
                                        else:
                                            j_m_new.append(j_info_new+[complement(j[-2]),j[-1]])
                    if j[0] in temp_letter[0] and j[2] in temp_letter[0] and j[4]+'^' in temp_letter[0] and j[6]+'^' in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]]:
                            for kb in Let_BP_Info['m'][j[2]]:
                                for kc in Let_BP_Info['m'][j[4]+'^']:
                                    for kd in Let_BP_Info['m'][j[6]+'^']:
                                        j_info_new=[ka[0]+j[1],kb[0]+j[3],kd[1]-j[7],kc[1]-j[5]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),j[-2]])
                                        else:
                                            j_m_new.append(j_info_new+[j[-2],complement(j[-1])])
                    if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0] and j[4]+'^' in temp_letter[0] and j[6]+'^' in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]+'^']:
                            for kb in Let_BP_Info['m'][j[2]+'^']:
                                for kc in Let_BP_Info['m'][j[4]+'^']:
                                    for kd in Let_BP_Info['m'][j[6]+'^']:
                                        j_info_new=[kb[1]-j[3],ka[1]-j[1],kd[1]-j[7],kc[1]-j[5]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_m_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),complement(j[-2])])
                                        else:
                                            j_m_new.append(j_info_new+[complement(j[-2]),complement(j[-1])])
                    j_p_new=[]
                    if j[0] in temp_letter[1] and j[2] in temp_letter[1] and j[4] in temp_letter[1] and j[6] in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]]:
                            for kb in Let_BP_Info['p'][j[2]]:
                                for kc in Let_BP_Info['p'][j[4]]:
                                    for kd in Let_BP_Info['p'][j[6]]:
                                        j_info_new=[ka[0]+j[1],kb[0]+j[3],kc[0]+j[5],kd[0]+j[7]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],j[-2]])
                                        else:
                                            j_p_new.append(j_info_new+[j[-2],j[-1]])
                    if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1] and j[4] in temp_letter[1] and j[6] in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]+'^']:
                            for kb in Let_BP_Info['p'][j[2]+'^']:
                                for kc in Let_BP_Info['p'][j[4]]:
                                    for kd in Let_BP_Info['p'][j[6]]:
                                        j_info_new=[kb[1]-j[3],ka[1]-j[1],kc[0]+j[5],kd[0]+j[7]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[j[-1],complement(j[-2])])
                                        else:
                                            j_p_new.append(j_info_new+[complement(j[-2]),j[-1]])
                    if j[0] in temp_letter[1] and j[2] in temp_letter[1] and j[4]+'^' in temp_letter[1] and j[6]+'^' in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]]:
                            for kb in Let_BP_Info['p'][j[2]]:
                                for kc in Let_BP_Info['p'][j[4]+'^']:
                                    for kd in Let_BP_Info['p'][j[6]+'^']:
                                        j_info_new=[ka[0]+j[1],kb[0]+j[3],kd[1]-j[7],kc[1]-j[5]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),j[-2]])
                                        else:
                                            j_p_new.append(j_info_new+[j[-2],complement(j[-1])])
                    if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1] and j[4]+'^' in temp_letter[1] and j[6]+'^' in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]+'^']:
                            for kb in Let_BP_Info['p'][j[2]+'^']:
                                for kc in Let_BP_Info['p'][j[4]+'^']:
                                    for kd in Let_BP_Info['p'][j[6]+'^']:
                                        j_info_new=[kb[1]-j[3],ka[1]-j[1],kd[1]-j[7],kc[1]-j[5]]
                                        if j_info_new[0]>j_info_new[2]:
                                            j_p_new.append(j_info_new[2:4]+j_info_new[:2]+[complement(j[-1]),complement(j[-2])])
                                        else:
                                            j_p_new.append(j_info_new+[complement(j[-2]),complement(j[-1])])
                    j2_m_new+=j_m_new
                    j2_p_new=j_p_new
                j_m_3a=candidate_QC_Control(j_m_new)
                j_p_3a=candidate_QC_Control(j_p_new)
                if j_m_3a==[]:
                    jMapPenam+=1
                if j_p_3a==[]:
                    jMapPenam+=1
                if jMapPenam==2:
                    if j[0]==j[2]:
                        Total_Cov_For_Pen[j[0]]+=j[3]-j[1]
                    else:
                        Total_Cov_For_Pen[j[0]]+=Be_BP_Letter[j[0]]-j[1]
                        Total_Cov_For_Pen[j[2]]+=j[3]
                    if j[4]==j[6]:
                        Total_Cov_For_Pen[j[4]]+=j[7]-j[5]
                    else:
                        Total_Cov_For_Pen[j[4]]+=Be_BP_Letter[j[4]]-j[5]
                        Total_Cov_For_Pen[j[6]]+=j[7]
                    NoMapPenal+=2
                elif jMapPenam==1:
                    if j_m_3a==[]:
                        Map_P+=[jp3+['p']+[float(1)/float(len(j_p_3a))] for jp3 in j_p_3a]
                    elif j_p_3a==[]:
                        Map_M+=[jp3+['m']+[float(1)/float(len(j_m_3a))] for jp3 in j_m_3a]
                else:
                    j_mp_4a=candidate_QC_Control2(j_m_3a,j_p_3a)
                    if not j_mp_4a==[]:
                        Map_Both+=[j4+[float(1)/float(len(j_mp_4a))] for j4 in j_mp_4a]
                    else:
                        if j[0]==j[2]:
                            Total_Cov_For_Pen[j[0]]+=j[3]-j[1]
                        else:
                            Total_Cov_For_Pen[j[0]]+=Be_BP_Letter[j[0]]-j[1]
                            Total_Cov_For_Pen[j[2]]+=j[3]
                        if j[4]==j[6]:
                            Total_Cov_For_Pen[j[4]]+=j[7]-j[5]
                        else:
                            Total_Cov_For_Pen[j[4]]+=Be_BP_Letter[j[4]]-j[5]
                            Total_Cov_For_Pen[j[6]]+=j[7]
                        NoMapPenal+=2
            return NoMapPenal
        def Be_Info_3_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
            be_info_3=Be_Info[2]
            for j in be_info_3:
                j_m_new=[]
                if j[0] in temp_letter[0] and j[2] in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]]:
                        for kb in Let_BP_Info['m'][j[2]]:
                            temp_single=[ka[0]+j[1],kb[0]+j[3]]
                            if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
                                j_m_new.append(temp_single)
                if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0]:
                    for ka in Let_BP_Info['m'][j[0]+'^']:
                        for kb in Let_BP_Info['m'][j[2]+'^']:
                            temp_single=[kb[1]-j[3],ka[1]-j[1]]
                            if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
                                j_m_new.append(temp_single)
                j_p_new=[]
                if j[0] in temp_letter[1] and j[2] in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]]:
                        for kb in Let_BP_Info['p'][j[2]]:
                            temp_single=[ka[0]+j[1],kb[0]+j[3]]
                            if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
                                j_p_new.append(temp_single)
                if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1]:
                    for ka in Let_BP_Info['p'][j[0]+'^']:
                        for kb in Let_BP_Info['p'][j[2]+'^']:
                            temp_single=[kb[1]-j[3],ka[1]-j[1]]
                            if not temp_single[1]-temp_single[0]>ReadLength*1.2 and temp_single[1]-temp_single[0]>0:
                                j_p_new.append(temp_single)
                if not j_m_new+j_p_new==[]:
                    for j2 in j_m_new:
                        Map_Both.append(j2+['m',float(1)/float(len(j_m_new+j_p_new))])
                    for j2 in j_p_new:
                        Map_Both.append(j2+['p',float(1)/float(len(j_m_new+j_p_new))])
                else:
                    Total_Cov_For_Pen[j[0]]=Be_BP_Letter[j[0]]-j[1]
                    Total_Cov_For_Pen[j[2]]=j[3]
                    NoMapPenal+=1
            return NoMapPenal
        def Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC,RD_within_B):
            Total_Cov_For_Pen={}
            for key in RD_within_B.keys():
                Total_Cov_For_Pen[key]=0
            Map_M=[]
            Map_P=[]
            Map_Both=[]
            Let_BP_Info={}
            Let_BP_Info['m']={}
            Let_BP_Info['p']={}
            temp_letter=[['left']+Af_Letter[0]+['right'],['left']+Af_Letter[1]+['right']]
            temp_bp=[[Af_BP[0][0]-flank]+Af_BP[0]+[Af_BP[0][-1]+flank],[Af_BP[1][0]-flank]+Af_BP[1]+[Af_BP[1][-1]+flank]]
            for j1 in range(len(temp_letter[0])):
                j=temp_letter[0][j1]
                if not j in Let_BP_Info['m'].keys():
                        Let_BP_Info['m'][j]=[[temp_bp[0][j1],temp_bp[0][j1+1]]]
                else:
                        Let_BP_Info['m'][j]+=[[temp_bp[0][j1],temp_bp[0][j1+1]]]
            for j1 in range(len(temp_letter[1])):
                j=temp_letter[1][j1]
                if not j in Let_BP_Info['p'].keys():
                    Let_BP_Info['p'][j]=[[temp_bp[1][j1],temp_bp[1][j1+1]]]
                else:
                    Let_BP_Info['p'][j]+=[[temp_bp[1][j1],temp_bp[1][j1+1]]]
            letters_numbers=[[Af_Letter[0].count(i[0])+Af_Letter[1].count(i[0])+Af_Letter[0].count(i[0]+'^')+Af_Letter[1].count(i[0]+'^') for i in Af_Letter[0]],[Af_Letter[0].count(i[0])+Af_Letter[1].count(i[0])+Af_Letter[0].count(i[0]+'^')+Af_Letter[1].count(i[0]+'^') for i in Af_Letter[1]]]
            NoMapPenal=0
            IL_Rec={}
            DR_Rec=0
            cov_bp=[[0 for i in range(len(temp_letter[0]))],[0 for i in range(len(temp_letter[1]))]]
            cov_bp2=[]
            NoMapPenal=Be_Info_1_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
            NoMapPenal=Be_Info_2_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
            NoMapPenal=Be_Info_3_rearrange(Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
            best_structure_sign_flag=0
            for key in Total_Cov_For_Pen.keys():
                if Total_Cov_For_Pen[key]==0:
                    del Total_Cov_For_Pen[key]
                else:
                    Total_Cov_For_Pen[key]/=float(Be_BP_Letter[key])
            for key in RD_within_B.keys():
                if not key[-1]=='^' and not key in ['left','right','left^', 'right^']:
                    if not key in Af_Letter[0]+Af_Letter[1] and not key+'^' in Af_Letter[0]+Af_Letter[1]:
                        if not key in Total_Cov_For_Pen.keys():
                            Total_Cov_For_Pen[key]=0
                        Total_Cov_For_Pen[key]+=RD_within_B[key]
            if NoMapPenal>0:
                best_structure_sign_flag+=1
            for key1 in Total_Cov_For_Pen.keys():
                if Total_Cov_For_Pen[key1]>2.58*GC_Std_Coverage[chrom_N]:
                    best_structure_sign_flag+=1
            if not Map_M+Map_P+Map_Both==[]:
                penals=penal_calculate(Map_M+Map_P+Map_Both,temp_bp,Af_Letter,Af_BP,RD_within_B,letters_numbers,BlockGC,NoMapPenal)
                if penals[2]>0:
                    best_structure_sign_flag+=1
                return penals[:-1]+[NoMapPenal,Total_Cov_For_Pen,best_structure_sign_flag]+[penals[-1]]
            else:
                return 0
        def Move_Decide_2(IL_List,RD_List,GC_Overall_Median_Coverage,GC_Var_Coverage):
            IL_Weight=1
            RD_Weight=1
            regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
            T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
            T2_Penal=[math.exp(k) for k in T_Penal]
            Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
            indicator=float(random.choice(range(1000)))/1000
            for j in range(len(IL_List)):
                cdf=numpy.sum(Normalized_Penal[:j+1])
                if cdf>indicator or cdf==indicator:
                    return [j,IL_List[j]*IL_Weight+RD_List[j]*RD_Weight]
        def Move_Decide_3(IL_List,RD_List,GC_Overall_Median_Coverage,GC_Var_Coverage):
            IL_Weight=1
            RD_Weight=1
            regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
            T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
            T2_Penal=[math.exp(k) for k in T_Penal]
            Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
            j=Normalized_Penal.index(max(Normalized_Penal))
            return [j,IL_List[j]*IL_Weight+RD_List[j]*RD_Weight]
        def Af_Letter_QC(Af_Letter,Copy_num_estimate):
            Copy_Num_Real={}
            for i in Copy_num_estimate.keys():
                Copy_Num_Real[i]=0
            for i in Af_Letter:
                for j in i:
                    Copy_Num_Real[j[0]]+=1
            FLAG=0
            for i in Copy_num_estimate.keys():
                if abs(Copy_num_estimate[i]-Copy_Num_Real[i])>1:
                    FLAG+=1
            return FLAG
        def getPwd(length):
            s=''.join([chr(97+i) for i in range(length)])
            sSet=set(s)
            if len(sSet)<length:return 0
            result=tmp=list(sSet)
            for i in range(length-1):
                result=createPwd(result,tmp)
            return result
        def createPwd(sList,kList):
            assert len(sList)*len(kList)>0
            upper,sList=len(sList),sList*len(kList)
            for i in range(len(kList)):
                for j in range(upper):
                    sList[i*upper+j]+=kList[i]
            return sList
        def changeLet(num):
            if num>1:
                let1=getPwd(num)
                let2=[]
                for k1 in let1:
                    let2.append([])
                    for k2 in k1:
                        let2[-1].append(chr(97+(ord(k2)-97)%2))
                let3=[]
                for k1 in let2:
                    if not k1 in let3:
                        let3.append(k1)
                let4=[]
                for k1 in let3:
                    temp=[]
                    for k2 in k1:
                        if k2=='a':
                            temp.append('a')
                        else:
                            temp.append('a^')
                    let4.append(temp)
                return let4
            elif num==1:
                return[['a'],['a^']]
            elif num==0:
                return [[]]
        def struc_propose_single_block(num):
            if num<1:
                return 'error'
            else:
                struc_rec=[]
                for i in range(num/2+1):
                    rec_a=changeLet(i)
                    rec_b=changeLet(num-i)
                    for k1 in rec_a:
                        for k2 in rec_b:
                            if not sorted([k1,k2]) in struc_rec:
                                struc_rec.append(sorted([k1,k2]))
                return struc_rec
        def bp_to_let(del_info_unit,chromos):
            flag=0
            for i in del_info_unit:
                if i  in chromos:
                    flag+=1
            letter=''.join([chr(i+97) for i in range(len(del_info_unit)-2*flag)])
            letters='/'.join([letter,letter])
            return letters
        def write_best_letter(bps_all,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters):
            fo=open(output_Score_File,'a')
            time2=time.time()
            Best_Letter_2=[]
            if not Score_rec_hash=={}:
                temp1=Best_Let_modify(Best_Letter_Rec,Best_Score_Rec,Score_rec_hash)
                Best_Letter_Rec=temp1[0]
                Best_Score_Rec=temp1[1]
            for bestletter in Best_Letter_Rec:
                if not sorted(bestletter) in Best_Letter_2:
                    Best_Letter_2.append(sorted(bestletter))
            bps3=[]
            for bps in bps_all:
                bps3+=bps
            for bestletter in Best_Letter_2:
                if not '/'.join([''.join(original_letters),''.join(original_letters)])=='/'.join([''.join(bestletter[0]),''.join(bestletter[1])]):
                    print >>fo, ' '.join([str(bp_ele) for bp_ele in bps3])
                    #print ' '.join([str(bp_ele) for bp_ele in bps3])
                    print >>fo, '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
                    #print '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
                    print >>fo, 'Theoretical Best Score: '+str(Best_IL_Score+Best_RD_Score)
                    print >>fo, 'Current Best Scure: '+str(Best_Score_Rec)
                    print >>fo, 'Time Consuming:'+str(datetime.timedelta(seconds=(time2-time1)))
            fo.close()
        def zero_RD_Process(run_flag):
            Best_Letter_Rec=[[[], []]]
            Af_Letter=[[],[]]
            Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
            Best_Score_Rec=Best_IL_Score+Best_RD_Score
            run_flag+=1
            return([Best_Letter_Rec,Best_Score_Rec,run_flag])
        def Af_TB_Penal_Less_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
            Af_TB_Rec=0
            Af_TB_Rec_list=Af_Info_all[-1]
            Af_TB_Penal_a=Af_Info_all[4]
            for x in Af_TB_Rec_list:
                if x <(IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5])/ReadLength*GC_Overall_Median_Num/float(10):
                    Af_TB_Rec+=1
            return -float(Af_TB_Penal_a)/float(num_of_reads)-float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
        def Af_TB_Penal_More_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
            Af_TB_Penal=Af_Info_all[3]
            Af_TB_Penal_a=Af_Info_all[4]
            return -float(Af_TB_Penal_a)/float(num_of_reads)-Af_TB_Penal
        def Af_Rearrange_Info_Collect(Letter_Candidates):
            P_IL=[]
            P_RD=[]
            P_DR=[]
            P_TB=[]
            Letter_Rec=[]
            BP_Rec=[]
            for Af_Letter in Letter_Candidates:
                Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
                for i in Af_Letter[0]:
                        Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
                for i in Af_Letter[1]:
                        Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
                Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
                if not Af_Info_all==0:
                    Letter_Rec.append(Af_Letter)
                    BP_Rec.append(Af_BP)
                    Af_IL_Penal=Af_Info_all[0]
                    Af_RD_Rec=Af_Info_all[1]
                    Af_DR_Penal=(Af_Info_all[2])**2
                    Af_TB_Penal_a=Af_Info_all[4]
                    Af_TB_Penal=Af_TB_Penal_Less_Important_Caldu(Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num)
                    Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
                    for key in Af_Info_all[5].keys():
                            Af_RD_Penal+=standard_norm_pdf_solver(Af_Info_all[5][key],0,GC_Mean_Coverage[chrom_N]/2)
                    P_IL.append(Af_IL_Penal)
                    P_RD.append(Af_RD_Penal)
                    P_DR.append(Af_DR_Penal/num_of_read_pairs)
                    P_TB.append(Af_TB_Penal)
            if len(P_IL)==0:
                return 'Error'
            else:
                ILTemp=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                RDTemp=[P_RD[i]*(1-P_TB[i]) for i in range(len(P_RD))] 
                return [ILTemp,RDTemp,Letter_Rec,BP_Rec]
        def Af_Rearrange_Info_Collect_M(M_Move_Choices):    
            P_IL=[]
            P_RD=[]
            P_DR=[]
            P_TB=[]
            Letter_Rec=[]
            BP_Rec=[]
            if not M_Move_Choices==[]:
                for m in [['2m','1','1','1','X']]+M_Move_Choices:
                        p=[str(Chr)+'p','1','1','1','X']
                        Move_MP=[m,p]
                        Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
                        Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
                        if Ploidy==1:
                            Af_Letter[1]=Af_Letter[0]
                            Af_BP[1]=Af_BP[0]
                        if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
                        if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
                        letter_num_flag=0
                        for key in Block_CN_Upper.keys():
                                if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
                                        letter_num_flag+=1
                        if not letter_num_flag==0: continue
                        Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
                        if Af_Info_all==0:continue
                        Letter_Rec.append(Af_Letter)
                        BP_Rec.append(Af_BP)
                        Af_IL_Penal=Af_Info_all[0]
                        Af_RD_Rec=Af_Info_all[1]
                        Af_DR_Penal=(Af_Info_all[2])**2
                        Af_TB_Penal_a=Af_Info_all[4]
                        Af_TB_Rec=Af_Info_all[3]
                        Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                        Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
                        for key in Af_Info_all[5].keys():
                            Af_RD_Penal+=standard_norm_pdf_solver(Af_Info_all[5][key],0,GC_Mean_Coverage[chrom_N]/2)
                            #Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
                        P_IL.append(Af_IL_Penal)
                        P_RD.append(Af_RD_Penal)
                        P_DR.append(Af_DR_Penal/num_of_read_pairs)
                        P_TB.append(Af_TB_Penal)
            if len(P_IL)==0: 
                return 'Error'
            else:
                ILTemp=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                RDTemp=[P_RD[i]*(1-P_TB[i]) for i in range(len(P_RD))] 
                return [ILTemp,RDTemp,Letter_Rec,BP_Rec]
        def Af_Rearrange_Info_Collect_P(P_Move_Choices):
            P_IL=[]
            P_RD=[]
            P_DR=[]
            P_TB=[]
            Letter_Rec=[]
            BP_Rec=[]
            if not P_Move_Choices==[]:
                for p in [['2p','1','1','1','X']]+P_Move_Choices:
                    m=[str(Chr)+'m','1','1','1','X']
                    Move_MP=[m,p]
                    Af_BP=[BPList_Rearrange(Be_BP[0],m,original_bp_list),BPList_Rearrange(Be_BP[1],p,original_bp_list)]
                    Af_Letter=[LetterList_Rearrange(Be_Letter[0],m,original_bp_list),LetterList_Rearrange(Be_Letter[1],p,original_bp_list)]
                    if Ploidy==1:
                        Af_Letter[1]=Af_Letter[0]
                        Af_BP[1]=Af_BP[0]
                    if not Af_Letter_QC(Af_Letter,Copy_num_estimate)==0:continue
                    if not Best_Score_Rec==0 and Af_Letter in Best_Letter_Rec: continue
                    letter_num_flag=0
                    for key in Block_CN_Upper.keys():
                            if (Af_Letter[0]+Af_Letter[1]).count(key)>Block_CN_Upper[key]:
                                    letter_num_flag+=1
                    if not letter_num_flag==0: continue
                    Af_Info_all=Letter_Through_Rearrange_4(IL_Statistics,Be_Info,Af_Letter,Af_BP,BlockGC2,RD_within_B)
                    if Af_Info_all==0: continue
                    Letter_Rec.append(Af_Letter)
                    BP_Rec.append(Af_BP)
                    Af_IL_Penal=Af_Info_all[0]
                    Af_RD_Rec=Af_Info_all[1]
                    Af_DR_Penal=(Af_Info_all[2])**2
                    Af_TB_Penal_a=Af_Info_all[4]
                    Af_TB_Rec=Af_Info_all[3]
                    Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                    Af_RD_Penal=RD_Adj_Penal(GC_Median_Coverage,GC_Overall_Median_Num,Chr,BlockGC2,Af_RD_Rec,Af_Letter)
                    for key in Af_Info_all[5].keys():
                        Af_RD_Penal+=standard_norm_pdf_solver(Af_Info_all[5][key],0,GC_Mean_Coverage[chrom_N]/2)
                        #Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
                    P_IL.append(Af_IL_Penal)
                    P_RD.append(Af_RD_Penal)
                    P_DR.append(Af_DR_Penal/num_of_read_pairs)
                    P_TB.append(Af_TB_Penal)
            if len(P_IL)==0:return 'Error'
            else:
                ILTemp=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                RDTemp=[P_RD[i]*(1-P_TB[i]) for i in range(len(P_RD))] 
                return [ILTemp,RDTemp,Letter_Rec,BP_Rec]
        def score_rec_hash_Modify_for_short_del(Score_rec_hash):
            Score_rec_hash_new={}
            for x in sorted(Score_rec_hash.keys())[::-1][:1]:
                Score_rec_hash_new[x]=Score_rec_hash[x]
            for x in sorted(Score_rec_hash.keys())[::-1][1:]:
                Score_rec_hash_new[x-1.1]=Score_rec_hash[x]
            return Score_rec_hash_new
        def one_RD_Process(run_flag,Score_rec_hash):
            if Ploidy==2:
                Letter_Candidates=[[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]]
            elif Ploidy==1:
                Letter_Candidates=[i for i in [[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]] if i[0]==i[1]]
            elif Ploidy==0:
                Letter_Candidates=[i for i in [[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a'], []],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]] if ['a'] in i]
            IL_RD_Temp_Info=Af_Rearrange_Info_Collect(Letter_Candidates)
            ILTemp=IL_RD_Temp_Info[0]
            RDTemp=IL_RD_Temp_Info[1]
            Letter_Rec=IL_RD_Temp_Info[2]
            BP_Rec=IL_RD_Temp_Info[3]
            Regu_IL=ILTemp
            Regu_RD=RDTemp
            if not ILTemp==[]:
                DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
                run_flag+=1
                for x in range(len(Letter_Rec)):
                    xy=ILTemp[x]+RDTemp[x]
                    if not xy in Score_rec_hash.keys():
                        Score_rec_hash[xy]=[]
                    Score_rec_hash[xy].append(Letter_Rec[x])               
            else:
                Best_Letter_Rec=[]
                Best_Score_Rec=100
                run_flag+=1
            Score_rec_hash2=score_rec_hash_Modify_for_short_del(Score_rec_hash)
            return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash2])
        def two_RD_Process(run_flag,Score_rec_hash):
            Letter_Candidates=struc_propose_single_block(2)+struc_propose_single_block(3)+struc_propose_single_block(4)+struc_propose_single_block(5)   
            if Ploidy==2:
                Letter_Candidates=Letter_Candidates
            elif Ploidy==1:
                Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
            elif Ploidy==0:
                Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
            IL_RD_Temp_Info=Af_Rearrange_Info_Collect(Letter_Candidates)
            ILTemp=IL_RD_Temp_Info[0]
            RDTemp=IL_RD_Temp_Info[1]
            Letter_Rec=IL_RD_Temp_Info[2]
            BP_Rec=IL_RD_Temp_Info[3]
            Regu_IL=ILTemp
            Regu_RD=RDTemp
            if not ILTemp==[]:
                DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
                run_flag+=1
                for x in range(len(Letter_Rec)):
                    xy=ILTemp[x]+RDTemp[x]
                    if not xy in Score_rec_hash.keys():
                        Score_rec_hash[xy]=[]
                    Score_rec_hash[xy].append(Letter_Rec[x])               
            else:
                Best_Letter_Rec=[]
                Best_Score_Rec=100
                run_flag+=1                     
            return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash])
        def few_RD_Process(run_flag,Score_rec_hash):
            Letter_Candidates=struc_propose_single_block(copy_num_a)+struc_propose_single_block(copy_num_b)
            if Ploidy==2:
                Letter_Candidates=Letter_Candidates
            elif Ploidy==1:
                Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
            elif Ploidy==0:
                Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
            IL_RD_Temp_Info=Af_Rearrange_Info_Collect(Letter_Candidates)
            ILTemp=IL_RD_Temp_Info[0]
            RDTemp=IL_RD_Temp_Info[1]
            Letter_Rec=IL_RD_Temp_Info[2]
            BP_Rec=IL_RD_Temp_Info[3]
            Regu_IL=ILTemp
            Regu_RD=RDTemp
            if not ILTemp==[]:
                DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
                run_flag+=1
                for x in range(len(Letter_Rec)):
                    xy=ILTemp[x]+RDTemp[x]
                    if not xy in Score_rec_hash.keys():
                        Score_rec_hash[xy]=[]
                    Score_rec_hash[xy].append(Letter_Rec[x])               
            else:
                Best_Letter_Rec=[]
                Best_Score_Rec=100
                run_flag+=1                     
            return([Best_Letter_Rec,Best_Score_Rec,run_flag,Score_rec_hash])
        def many_RD_Process(run_flag):
            Best_Letter_Rec=[[['a' for i in range(copy_num_a/2)],['a' for i in range(copy_num_a/2)]]]
            Best_Score_Rec=100
            run_flag+=1     
            return([Best_Letter_Rec,Best_Score_Rec,run_flag])
        def two_block_a(cpNum_a,cpNum_b):
            out=[]
            for k1 in cpNum_a:
                for k2 in cpNum_b:
                    out.append(['a' for k in range(int(k1))]+['b' for k in range(int(k2))])
            out2=[]
            for k1 in out:
                t1=[]
                ta=[]
                tb=[]
                for k2 in range(k1.count('a')+1):
                    ta.append([['a' for k in range(k2)],['a' for k in range(k1.count('a')-k2)]])
                for k2 in range(k1.count('b')+1):
                    tb.append([['b' for k in range(k2)],['b' for k in range(k1.count('b')-k2)]])
                for k2 in ta:
                    for k3 in tb:
                        if not sorted([k2[0]+k3[0],k2[1]+k3[1]]) in t1:
                            t1.append(sorted([k2[0]+k3[0],k2[1]+k3[1]]))
                out2+=t1
            return out2
        def two_block_b(all_Strucs):
            out=[]
            for k1 in all_Strucs:
                k1a=[]
                k1b=[]
                rect=0
                if k1[0]==[]:
                    k1a2=[[]]
                else:
                    for k2 in range(len(k1[0])+1):
                        all_lets=k1[0]+['^' for ki in range(k2)]
                        for k3 in itertools.permutations(all_lets,len(all_lets)):
                            if not k3[0]=='^' and not k3 in k1a:
                                k3flag=0
                                for k4 in range(len(k3[:-1])):
                                    if k3[k4+1]==k3[k4]=='^':
                                        k3flag+=1
                                if k3flag==0:
                                    k1a.append(k3)
                    k1a2=[]
                    for k2 in k1a:
                        tk2=[]
                        for k3 in k2:
                            if k3=='^':
                                tk2[-1]+=k3
                            else:
                                tk2.append(k3)
                        k1a2.append(tk2)
                if k1[1]==[]:
                    k1b2=[[]]
                else:
                    for k2 in range(len(k1[1])+1):
                        all_lets=k1[1]+['^' for ki in range(k2)]
                        for k3 in itertools.permutations(all_lets,len(all_lets)):
                            if not k3[0]=='^' and not k3 in k1b:
                                k3flag=0
                                for k4 in range(len(k3[:-1])):
                                    if k3[k4+1]==k3[k4]=='^':
                                        k3flag+=1
                                if k3flag==0:
                                    k1b.append(k3)
                    k1b2=[]
                    for k2 in k1b:
                        tk2=[]
                        for k3 in k2:
                            if k3=='^':
                                tk2[-1]+=k3
                            else:
                                tk2.append(k3)
                        k1b2.append(tk2)
                for ka in k1a2:
                    for kb in k1b2:
                        out.append([ka,kb])
            return out
        def struc_produce_two_block(Copy_num_estimate):
            cpNum_a=[j for j in [Copy_num_estimate['a']+i for i in [-1,0,1]] if j>-1]
            cpNum_b=[j for j in [Copy_num_estimate['b']+i for i in [-1,0,1]] if j>-1]
            all_Strucs=two_block_a(cpNum_a,cpNum_b)
            all_Str2=two_block_b(all_Strucs)
            return all_Str2
        def two_block_RD_Process(run_flag):
            Letter_Candidates=struc_produce_two_block(Copy_num_estimate)
            if Ploidy==2:
                Letter_Candidates=Letter_Candidates
            elif Ploidy==1:
                Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
            elif Ploidy==0:
                Letter_Candidates=[i for i in Letter_Candidates if ['a','b'] in i]
            IL_RD_Temp_Info=Af_Rearrange_Info_Collect(Letter_Candidates)
            ILTemp=IL_RD_Temp_Info[0]
            RDTemp=IL_RD_Temp_Info[1]
            Letter_Rec=IL_RD_Temp_Info[2]
            BP_Rec=IL_RD_Temp_Info[3]
            Regu_IL=ILTemp
            Regu_RD=RDTemp
            if not ILTemp==[]:
                DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                Best_Score_Rec=Regu_IL[DECISION_Score[0]]+Regu_RD[DECISION_Score[0]]
                run_flag+=1
            else:
                Best_Letter_Rec=[]
                Best_Score_Rec=100
                run_flag+=1
            return([Best_Letter_Rec,Best_Score_Rec,run_flag])
        def qual_check_bps2(bps2):
            flag=0
            if len(bps2)==1 and len(bps2[0])==3:
                for x1 in bps2:
                    for x2 in range(len(x1)-2):
                        if int(x1[x2+2])-int(x1[x2+1])<Window_Size: 
                            flag=1
            if flag==0:
                return 'right'
            else:
                return 'error'
        def let_reclust(vec_in):
            if vec_in==[]:
                return []
            else:
                k2e=[]
                k2e=[vec_in[0]]
                for k3 in range(len(vec_in)-1):
                    if '^' in vec_in[k3+1]:
                        if '^' in vec_in[k3] and ord(vec_in[k3][0])-ord(vec_in[k3+1][0])==1:
                            k2e[-1]+=vec_in[k3+1]
                        else:
                            k2e.append(vec_in[k3+1])
                    else:
                        if ord(vec_in[k3+1][0])-ord(vec_in[k3][0])==1 and not '^' in vec_in[k3]:
                            k2e[-1]+=vec_in[k3+1]
                        else:
                            k2e.append(vec_in[k3+1])
                k2f=[]
                for k3 in k2e:
                    if '^' in k3:
                        k5=''
                        for k4 in range(len(k3)/2):
                            k5+=k3[2*k4]
                        k6=k5[::-1]+'^'
                        if not k6 in k2f:
                            k2f.append(k6)
                    else:
                        k2f.append(k3)
                return k2f 
        def simple_flag_SA(k1,k2):
            if k2=='':
                return[[i for i in k1],[],[],0]
            else:
                temp=[]
                break_flag=0
                for i in k2:
                    if not i=='^':
                        temp.append(i)
                    else:
                        temp[-1]+=i
                temp2=[temp[0]]
                for i in range(len(temp[1:])):
                    if not '^' in temp[i] and not '^' in temp[i+1] and ord(temp[i+1])-ord(temp[i])==1:
                        temp2[-1]+=temp[i+1]
                    elif '^' in temp[i] and '^' in temp[i+1] and ord(temp[i+1][0])-ord(temp[i][0])==-1:
                        temp2[-1]=temp[i+1][0]+temp2[-1]
                    else:
                        temp2.append(temp[i+1]) 
                outdel=[]
                outinv=[]
                outdup=[]
                outtra=0
                for i in range(len(temp2)):
                    j=temp2[i]
                    if '^' in j:
                        if not j.replace('^','') in outinv:
                            outinv.append(j.replace('^',''))
                        temp2[i]=j.replace('^','')
                temp3=''.join(temp2)
                for i in range(len(temp3)-1):
                    if not temp3[i+1]=='/' and not temp3[i]=='/':
                        if ord(temp3[i+1])-ord(temp3[i])<0:
                            outtra=1
                if not temp3==k1:
                    temp4=[]
                    for i in temp3:
                        if temp3.count(i)>1:
                            if not i in outdup:
                                outdup.append(i)
                        if not i in temp4:
                            temp4.append(i)
                    if not ''.join(temp4)==k1:
                        for i in k1:
                            if not i in temp4:
                                outdel.append(i)
                if not outdup==[]:
                    dupuni=unit_produce(outdup)
                    outdup2=[]
                    k3=k2
                    for i in dupuni:
                        ia=i
                        ib=''.join([j+'^' for j in i[::-1]])
                        if len(i)>1:
                            if temp2.count(ia)+temp2.count(ib)>2:
                                outdup2.append([i,temp2.count(ia)+temp2.count(ib)])
                                k3=k3.replace(ia,'')
                                k3=k3.replace(ib,'')
                        elif len(i)==1:
                            if k3.count(ia)+k3.count(ib)>1:
                                outdup2.append([i,k3.count(ia)])
                                k3=k3.replace(ia,'')
                                k3=k3.replace(ib,'')
                else:
                    outdup2=[]
                outdup3=[]
                for i in outdup2:
                    flag=0
                    for j in outdup2:
                        if not j==i:
                            if len(j[0])>len(i[0]) and i[0] in j[0]:
                                flag+=1
                    if flag==0:
                        outdup3.append(i)       
                if len(outdup3)>1:
                    outdup3.sort()
                    outdup4=[outdup3[0]]
                    for i in range(len(outdup3)-1):
                        if outdup3[i+1][-1]==outdup3[i][-1] and ord(outdup3[i+1][0][0])-ord(outdup3[i][0][-1])==1:
                            outdup4[-1][0]+=outdup3[i+1][0]
                        else:
                            outdup4.append(outdup3[i+1])
                else:
                    outdup4=outdup3
                return [outdel,outinv,outdup4,outtra]
        def unit_produce(list):
            temp1=[sorted(list)[0]]
            for k1 in sorted(list)[1:]:
                if ord(k1)-ord(temp1[-1][-1])==1:
                    temp1[-1]+=k1
                else:
                    temp1.append(k1)
            temp2=[]
            for k1 in temp1:
                for k2 in range(len(k1)+1)[1:]:
                    for k3 in range(len(k1)-k2+1):
                        temp2.append(k1[k3:(k3+k2)])
            return temp2[::-1]
        def simpler_struc_pick(list_of_structures,original_letters):
            for x in list_of_structures:
                x.sort()
            out=[]
            for x in list_of_structures:
                out.append(0)
                for y in x:
                    temp=simple_flag_SA(''.join(original_letters),''.join(y))
                    temp2={}
                    for x in temp[0]:
                        if not x in temp2.keys():
                            temp2[x]=1
                        else:
                            temp2[x]+=1
                    for x in temp[1]:
                        if not x in temp2.keys():
                            temp2[x]=1
                        else:
                            temp2[x]+=1
                    for x in temp[2]:
                        if not x[0] in temp2.keys():
                            temp2[x[0]]=1
                        else:
                            temp2[x[0]]+=1
                    for x in temp2.keys():
                        out[-1]+=float(temp2[x]-1)/float(2)+1
            out2=[list_of_structures[i] for i in range(len(out)) if out[i]==min(out)]
            return out2
        def Best_Let_modify(Best_Letter_Rec,Best_Score_Rec,Score_rec_hash):
            if not Score_rec_hash=={}:
                temp1=[Best_Score_Rec]
                for x in sorted(Score_rec_hash.keys())[::-1][:2]:
                    if not x in temp1 and x-Best_Score_Rec> -1:
                        temp1.append(x)
                temp2=[]
                for x in temp1:
                    temp2+=Score_rec_hash[x]
                temp3=simpler_struc_pick(temp2,original_letters)
                if len(temp3)==1:
                    Best_Letter_Rec=temp3
                    for x in temp1:
                        if temp3[0] in Score_rec_hash[x]:
                            new_best_score=x
                    return [Best_Letter_Rec,new_best_score]
                else:
                    thash=[]
                    for x in temp3:
                        for y in temp1:
                            if x in Score_rec_hash[y]:
                                thash.append(y)
                    thash2=[temp3[i] for i in range(len(temp3)) if thash[i]==min(thash)]
                    Best_Letter_Rec=thash2
                    new_best_score=min(thash)
                    return [Best_Letter_Rec,new_best_score] 
            else:
                return [Best_Letter_Rec,Best_Score_Rec] 
        def ori_let_Modi(Be_Info, ori_let2):
            Ab_Dir_L1_a=[]
            Ab_Dir_L1_b=[]
            Ab_Dir_L2_a=[]
            Ab_Dir_L2_b=[]
            norm_dir={}
            for y in Be_Info[0]:
                if y[-2:]==['+','+']:
                    Ab_Dir_L1_a.append(y)
                elif y[-2:]==['-','-']:
                    Ab_Dir_L2_a.append(y)
            for y in Be_Info[1]:
                if y[-2:]==['+','+']:
                    Ab_Dir_L1_b.append(y)
                elif y[-2:]==['-','-']:
                    Ab_Dir_L2_b.append(y)
            vote_inv_hash={}
            for x in Ab_Dir_L1_a:
                if not x[3] in vote_inv_hash.keys():
                    vote_inv_hash[x[3]]=1
                else:
                    vote_inv_hash[x[3]]+=1
            for x in Ab_Dir_L1_b:
                if not x[4] in vote_inv_hash.keys():
                    vote_inv_hash[x[4]]=0.5
                else:
                    vote_inv_hash[x[4]]+=0.5
                if not x[6] in vote_inv_hash.keys():
                    vote_inv_hash[x[6]]=0.5
                else:
                    vote_inv_hash[x[6]]+=0.5
            for x in Ab_Dir_L2_a:
                if not x[0] in vote_inv_hash.keys():
                    vote_inv_hash[x[0]]=1
                else:
                    vote_inv_hash[x[0]]+=1
            for x in Ab_Dir_L2_b:
                if not x[0] in vote_inv_hash.keys():
                    vote_inv_hash[x[0]]=0.5
                else:
                    vote_inv_hash[x[0]]+=0.5
                if not x[2] in vote_inv_hash.keys():
                    vote_inv_hash[x[2]]=0.5
                else:
                    vote_inv_hash[x[2]]+=0.5
            new_let2=[]
            for x in ori_let2:
                if x in vote_inv_hash.keys() and vote_inv_hash[x]>2:
                    new_let2.append(x+'^')
                else:
                    new_let2.append(x)
            new_group=[]
            for x in new_let2:
                if new_group==[]:
                    new_group.append([x])
                else:
                    if x.count('^')==new_group[-1][-1].count('^'):
                        new_group[-1].append(x)
                    else:
                        new_group.append([x])
            out=[]
            for x in new_group:
                if '^' in x[0] and len(x)>1:
                    out+=x[::-1]
                else:
                    out+=x
            return out
        def Define_Default_SVPredict():
            global Penalty_For_InsertLengthZero
            Penalty_For_InsertLengthZero=-20 #Toy example,decides later 
            if not '/' in dict_opts['--bp-file']:
                dict_opts['--bp-file']='./'+dict_opts['--bp-file']
            global model_comp
            if not '--null-model' in dict_opts.keys():
                model_comp='S'
            else:
                if dict_opts['--null-model'] in ['S','Simple']:
                    model_comp='S'
                else:
                    model_comp='C'
            global Ploidy
            if '--ploidy' in dict_opts.keys():
                Ploidy=int(dict_opts['--ploidy'])
            else:
                Ploidy=2
            global QCAlign
            if '--qc-align' in dict_opts.keys():
                QCAlign=int(dict_opts['--qc-align'])
            else:
                QCAlign=20
            global genome_name
            if '--NullGenomeName' in dict_opts.keys():
                genome_name=dict_opts['--NullGenomeName']
            else:
                genome_name='genome'
            global Trail_Number
            if '--num-iteration' in dict_opts.keys():
                Trail_Number=int(dict_opts['--num-iteration'])
            else:
                Trail_Number=100000
            global IL_Weight
            global RD_Weight
            global DR_Weight
            global TB_Weight
            IL_Weight=1
            RD_Weight=1
            DR_Weight=5
            TB_Weight=5
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def GC_RD_Info_Complete(ref_file):
            global refFlag
            refFlag=0
            refgenome=ref_file
            reftest=open(refgenome)
            preftest=reftest.readline().strip().split()
            if len(preftest[0])>3:
                refFlag=1
            reftest.close()
            global chrom_N
            global chrom_X
            global chrom_Y  
            if refFlag==0:
                chrom_N='N'
                chrom_X='X'
                chrom_Y='Y'
            else:
                chrom_N='chrN'
                chrom_X='chrX'
                chrom_Y='chrY'      
            for x in [chrom_N,chrom_X,chrom_Y]:
                if not x in GC_Median_Coverage.keys():
                    GC_Median_Coverage[x]={}
                for i in ChrN_Median_Coverage.keys():
                    GC_Median_Coverage[x][i]=numpy.mean(ChrN_Median_Coverage[i])
            for x in [chrom_N,chrom_X,chrom_Y]:
                GC_Overall_Median_Coverage[x]=numpy.mean([GC_Overall_Median_Coverage[key] for key in GC_Overall_Median_Coverage.keys()])
                GC_Var_Coverage[x]=numpy.mean([GC_Var_Coverage[key] for key in GC_Var_Coverage.keys()])
                GC_Mean_Coverage[x]=numpy.mean([GC_Mean_Coverage[key] for key in GC_Mean_Coverage.keys()])
                GC_Std_Coverage[x]=numpy.mean([GC_Std_Coverage[key] for key in GC_Std_Coverage.keys()])
            for x in Chromosome:
                if not x in GC_Median_Coverage.keys():
                    GC_Median_Coverage[x]={}
                for i in ChrN_Median_Coverage.keys():
                    GC_Median_Coverage[x][i]=numpy.mean(ChrN_Median_Coverage[i])
                if not x in GC_Overall_Median_Coverage.keys():
                    GC_Overall_Median_Coverage[x]=GC_Overall_Median_Coverage[chrom_N]
                if not x in GC_Var_Coverage.keys():
                    GC_Var_Coverage[x]=GC_Var_Coverage[chrom_N]
                if not x in GC_Mean_Coverage.keys():
                    GC_Mean_Coverage[x]=GC_Mean_Coverage[chrom_N]
                if not x in GC_Std_Coverage.keys():
                    GC_Std_Coverage[x]=GC_Std_Coverage[chrom_N]
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1        Last Update:2015-08-20'
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--sample, input alignment file in bam format'
            print '--bp-file, input txt file containing clustered bps.'
            print ' '
            print 'Optional Parameters:'
            print '--num-iteration, maximum number of iterations per structure will run'
            print '--ploidy, limit algorithm to specific zygosity (0:heterozygous only; 1:homozygous only; 2:both; default:2)'
            print '--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
            print '--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)'
        else:
            Define_Default_SVPredict()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                if not '--bp-file' in dict_opts.keys():
                    print 'Error: please specify input txt file using : --bp-file'
                else:
                    if not '-o' in dict_opts.keys():
                        dict_opts['-o']='/'.join(dict_opts['--bp-file'].split('/')[:-1])
                    if not dict_opts['-o'][-1]=='/':
                        dict_opts['-o']+='/'
                    ref_path=workdir+'reference/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    ref_ppre=ref_path
                    ref_prefix='.'.join(ref_file.split('.')[:-1])
                    chromos_all=[]
                    refFile=ref_file
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        refIndFile=refFile+'.fai'
                        fin=open(refIndFile)
                        for line in fin:
                            pin=line.strip().split()
                            chromos_all.append(pin[0])
                        fin.close()
                        if not '--sample' in dict_opts.keys():
                            print 'Error: please specify either input file using --sample'
                        else:
                            time1=time.time()
                            BamN=dict_opts['--sample'].split('/')[-1].replace('.bam','')
                            Input_File=dict_opts['--bp-file']
                            Insert_Len_Stat=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/ILNull.'+BamN+'.'+genome_name+'.Bimodal'
                            if not os.path.isfile(Insert_Len_Stat):
                                print 'wrong workdir defined'
                            ReadLenFin=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/'+BamN+'.'+genome_name+'.Stats'
                            if not os.path.isfile(ReadLenFin):
                                print 'wrong workdir defined'
                            else:
                                fin=open(ReadLenFin)
                                pin=fin.readline().strip().split()
                                pin=fin.readline().strip().split()
                                pin=fin.readline().strip().split()
                                Window_Size=int(pin[0])/3                   
                                for line in fin:
                                    pin=line.strip().split()
                                fin.close()
                                ReadLength=int(pin[-1].split(':')[-1])
                            Initial_Bam_Name=BamN+'.bam'
                            Initial_Bam=dict_opts['--sample']
                            flank=cdf_solver_application(Insert_Len_Stat,0.95)
                            Cut_Lower=cdf_solver_application(Insert_Len_Stat,0.005)
                            Cut_Upper=cdf_solver_application(Insert_Len_Stat,0.995)
                            if model_comp=='C':
                                IL_Stat_all=IL_Stat(Insert_Len_Stat)
                                IL_Statistics=IL_Stat_all[0]
                                IL_Normal_Stat=IL_Stat_all[1]
                            elif model_comp=='S':
                                IL_Stat_all=IL_Stat_Simp(Insert_Len_Stat)
                                IL_Statistics=IL_Stat_all[0]
                                IL_Normal_Stat=IL_Stat_all[1]
                            IL_Estimate=IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5]
                            IL_SD=((IL_Statistics[2]*IL_Statistics[4])**2+(IL_Statistics[3]*IL_Statistics[5])**2)**(0.5)
                            IL_Penal_Two_End_Limit=min([standard_pdf_IL_calculate(IL_Estimate-3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
                            ,standard_pdf_IL_calculate(IL_Estimate+3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)])
                            low_qual_edge=5
                            fi=open(Input_File)
                            bps_hash={}
                            bps_temp=[]
                            break_flag=0
                            for line in fi:
                                pi=line.strip().split()
                                if pi==[] or len(pi)<3:
                                    if bps_temp==[]:
                                        continue
                                    else:
                                        bp_key=0
                                        for l1 in bps_temp:
                                            bp_key+=len(l1)
                                        if not bp_key in bps_hash.keys():
                                            bps_hash[bp_key]=[]
                                        bps_hash[bp_key].append(bps_temp)
                                        bps_temp=[]
                                        if bp_key<3: 
                                            pi=pi
                                else:
                                    bps_temp.append(pi)
                            fi.close()
                            bps_hash_inter={}
                            for k1 in bps_hash.keys():
                                bps_hash_inter[k1]=[]
                                for k2 in bps_hash[k1]:
                                    if not k2 in bps_hash_inter[k1]:
                                        bps_hash_inter[k1].append(k2)
                            bps_hash=bps_hash_inter
                            output_Score_File=dict_opts['-o']+'_'.join(dict_opts['--bp-file'].split('/')[-1].split('.')[:-1]+['flag'+str(Ploidy)])+'_'+'MP'+str(QCAlign)+'.coverge'
                            fo=open(output_Score_File,'w')
                            fo.close()
                            for bpsk1 in sorted(bps_hash.keys()):
                                for bps2 in bps_hash[bpsk1]:
                                    for i in bps2:
                                        if len(i)<3:
                                            i.append(str(int(i[-1])+Window_Size))
                            GC_Stat_Path=dict_opts['--workdir']+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/RD_Stat'
                            Affix_GC_Stat='_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
                            GC_Stat=GC_Stat_ReadIn(BamN,GC_Stat_Path,Affix_GC_Stat)
                            GC_Content_Coverage=GC_Stat[0]
                            Chromosome=GC_Stat[1]
                            Coverage=[int(k) for k in GC_Stat[2][1:]]
                            GC_Overall_Median_Coverage={}
                            GC_Overall_Median_Num=[]
                            GC_Median_Coverage={}
                            GC_Median_Num={}
                            GC_Mean_Coverage={}
                            GC_Std_Coverage={}
                            GC_Var_Coverage={}
                            for a in Chromosome:
                                GC_Overall_temp=[]
                                for b in Coverage:
                                    if not b in GC_Median_Num.keys():
                                        GC_Median_Num[b]=[]
                                    if len(GC_Content_Coverage[a][b][0])==2: continue
                                    elif len(GC_Content_Coverage[a][b][0])>2:
                                            num_list=[float(c) for c in GC_Content_Coverage[a][b][0][2:].split(',')]
                                            if not sum(num_list)==0:
                                                GC_Median_Num[b]+=num_list
                                                GC_Overall_Median_Num+=num_list
                                                GC_Overall_temp=GC_Overall_temp+num_list
                                                if not Median_Pick(num_list)==0.0:
                                                    if not a in GC_Median_Coverage.keys():
                                                        GC_Median_Coverage[a]={}
                                                    GC_Median_Coverage[a][b]=Median_Pick(num_list)
                                if len(GC_Overall_temp)==0: continue
                                if sum(GC_Overall_temp)==0.0: continue
                                elif len(GC_Overall_temp)>0: 
                                        GC_Overall_Median_Coverage[a]=Median_Pick(GC_Overall_temp)
                                        GC_Mean_Coverage[a]=numpy.mean(GC_Overall_temp)
                                        GC_Std_Coverage[a]=numpy.std(GC_Overall_temp)
                                        GC_Var_Coverage[a]=(GC_Std_Coverage[a])**2
                            GC_Overall_Median_Num=Median_Pick(GC_Overall_Median_Num)
                            for a in GC_Median_Num.keys():
                                if GC_Median_Num[a]==[]:
                                    GC_Median_Num[a]=GC_Overall_Median_Num
                                else:
                                    GC_Median_Num[a]=Median_Pick(GC_Median_Num[a])
                            ChrN_Median_Coverage={}
                            for i in GC_Median_Coverage.keys():
                                for j in GC_Median_Coverage[i].keys():
                                    if not j in ChrN_Median_Coverage.keys():
                                        ChrN_Median_Coverage[j]=[GC_Median_Coverage[i][j]]
                                    else:
                                        ChrN_Median_Coverage[j]+=[GC_Median_Coverage[i][j]]
                            GC_RD_Info_Complete(ref_file)
                            for bpsk1 in sorted(bps_hash.keys()):
                                for bps2 in bps_hash[bpsk1]:
                                    if qual_check_bps2(bps2)=='right':
                                        Chromo=bps2[0][0]
                                        K_RD=GC_Std_Coverage[str(Chromo)]/GC_Mean_Coverage[str(Chromo)]
                                        K_IL=IL_Normal_Stat[2]/IL_Normal_Stat[1]
                                        K_RD_new=1
                                        K_IL_new=(K_IL/K_RD)**2
                                        IL_GS=Prob_Norm(IL_Normal_Stat[1],IL_Normal_Stat[1],IL_Normal_Stat[2]**2)
                                        RD_GS=Prob_Norm(GC_Mean_Coverage[str(Chromo)],GC_Mean_Coverage[str(Chromo)],GC_Std_Coverage[str(Chromo)]**2)
                                        for i in bps2:
                                            temp2=[int(j) for j in i[1:]]
                                            k=[i[0]]+sorted(temp2)
                                            k2=k[:2]
                                            for k3 in temp2:
                                                if not k3 in k2 and k3-k2[-1]>10:
                                                    k2.append(k3)
                                            if len(k2)>2:
                                                bps2[bps2.index(i)]=k2
                                            else:
                                                del bps2[bps2.index(i)]
                                        original_bps_all=[]
                                        for obas in bps2:
                                            original_bps_all+=obas
                                        original_structure=bp_to_let(original_bps_all,chromos_all)
                                        chr_letter_tbp=letter_rearrange(bps2)
                                        letter_tGC=letter_GC_ReadIn(chr_letter_tbp)
                                        if letter_tGC=='error': continue
                                        letter_tRD=letter_RD_ReadIn(chr_letter_tbp)
                                        if letter_tRD=='error': continue
                                        chr_letter_bp={}
                                        letter_GC={}
                                        letter_RD={}
                                        for k1 in chr_letter_tbp.keys():
                                            chr_letter_bp[k1]={}
                                            letter_GC[k1]={}
                                            letter_RD[k1]={}
                                            for k2 in chr_letter_tbp[k1].keys():
                                                if k2 in letter_tGC[k1].keys() and k2 in letter_tRD[k1].keys() and not math.isnan(letter_tRD[k1][k2]) and not math.isnan(letter_tGC[k1][k2]):
                                                    chr_letter_bp[k1][k2]=chr_letter_tbp[k1][k2]
                                                    letter_GC[k1][k2]=letter_tGC[k1][k2]
                                                    letter_RD[k1][k2]=letter_tRD[k1][k2]
                                        left_keys=[]
                                        for k1 in chr_letter_bp.keys():
                                            for k2 in chr_letter_bp[k1].keys():
                                                left_keys.append(k2)
                                        if not left_keys==[]:
                                            bamFlag=0
                                            bamtest=os.popen(r'''samtools view -H %s'''%(dict_opts['--sample']))
                                            for line in bamtest:
                                                pbamtest=line.strip().split()
                                                if pbamtest[0]=='@SQ' and pbamtest[1].split(':')[0]=='SN':
                                                    if len(pbamtest[1].split(':')[1])>3:
                                                        bamFlag+=1
                                                    break
                                            bps3={}
                                            for k1 in chr_letter_bp.keys():
                                                bps3[k1]={}
                                                for k2 in chr_letter_bp[k1].keys():
                                                    bps3[k1][chr_letter_bp[k1][k2][0]]=[chr_letter_bp[k1][k2][0],chr_letter_bp[k1][k2][-1]]
                                            bps4={}
                                            for k1 in bps3.keys():
                                                if not bps3[k1]=={}:
                                                    bps4[k1]=[[k1]+bps3[k1][sorted(bps3[k1].keys())[0]]]
                                                    for k2 in range(len(bps3[k1].keys())-1):
                                                        if bps3[k1][sorted(bps3[k1].keys())[k2+1]][0]==bps3[k1][sorted(bps3[k1].keys())[k2]][-1]:
                                                            bps4[k1][-1]+=[bps3[k1][sorted(bps3[k1].keys())[k2+1]][-1]]
                                                        else:
                                                            bps4[k1].append([bps3[k1][sorted(bps3[k1].keys())[k2+1]]])
                                            bps2=[]
                                            for k1 in bps4.keys():
                                                for k2 in bps4[k1]:
                                                    bps2.append(k2)
                                            Chr=bps2[0][0]
                                            Through_BP_Minimum=GC_Overall_Median_Coverage[str(Chr)]/10
                                            if len(Chr)>3:
                                                refChr0=Chr[3:]
                                                refChr1=Chr
                                            else:
                                                refChr0=Chr
                                                refChr1='chr'+Chr
                                            if refFlag==0:
                                                refChr=refChr0
                                            else:
                                                refChr=refChr1
                                            if bamFlag==0:
                                                bamChr=refChr0
                                            else:
                                                bamChr=refChr1
                                            bamtest.close()
                                            local_Minimum=[]
                                            block2_bin=[]
                                            GC_blocks_index=[]
                                            Full_Info=Full_Info_of_Reads_Integrate(bps2)
                                            RD_within_B=Full_Info[0]
                                            RD_within_B['left']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
                                            RD_within_B['right']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
                                            for i in RD_within_B.keys():
                                                RD_within_B[i+'^']=RD_within_B[i]
                                            Initial_GCRD_Adj=Full_Info[1]
                                            Initial_GCRD_Adj['left']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
                                            Initial_GCRD_Adj['right']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
                                            Best_IL_Score=0
                                            for j in range(Cut_Lower,Cut_Upper+1):
                                                Single_ILScore=standard_pdf_IL_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
                                                Best_IL_Score+=Single_ILScore*exp(Single_ILScore)
                                            Best_RD_Score=0
                                            let_chr_rec={}
                                            for i in chr_letter_bp.keys():
                                                for j in chr_letter_bp[i].keys():
                                                    if j in left_keys:
                                                        let_chr_rec[j]=i
                                            for i in let_chr_rec.keys():
                                                Theo_RD=GC_Overall_Median_Coverage[str(let_chr_rec[i])]
                                                Theo_Var=GC_Var_Coverage[str(let_chr_rec[i])]
                                                for j in range(int(Theo_RD/2),int(Theo_RD/2*3+1)):
                                                    single_ProbNB=Prob_Norm(j,Theo_RD,Theo_Var)
                                                    Best_RD_Score+=single_ProbNB*exp(single_ProbNB)
                                            Block_CN_Upper={}
                                            Copy_num_estimate={}
                                            for i in Initial_GCRD_Adj.keys():
                                                if not i in ['left','right']:
                                                    Copy_num_estimate[i]=round(Initial_GCRD_Adj[i]*2/GC_Mean_Coverage[Chr])
                                                    if Initial_GCRD_Adj[i]<float(GC_Mean_Coverage[Chr])/10.0:
                                                        Copy_num_estimate[i]=-1
                                            Copy_num_Check=[]
                                            for CNE in Copy_num_estimate.keys():
                                                if Copy_num_estimate[CNE]>5:
                                                    Copy_num_Check.append(CNE)
                                            if Copy_num_Check==[]:
                                                median_CN=GC_Overall_Median_Coverage[chrom_N]/2
                                                for key in Initial_GCRD_Adj.keys():
                                                    if not key in ['left','right']:
                                                        Block_CN_Upper[key]=Initial_GCRD_Adj[key]/median_CN+2
                                                Initial_DR=Full_Info[2]
                                                Initial_IL=Full_Info[3]
                                                BlockGC=Full_Info[7]
                                                BlockGC['left']=0.476
                                                BlockGC['right']=0.476
                                                BlockGC2={}
                                                for key_B_GC in BlockGC.keys():
                                                    BlockGC2[key_B_GC]=BlockGC[key_B_GC]
                                                    BlockGC2[key_B_GC+'^']=BlockGC[key_B_GC]    
                                                original_letters=Full_Info[9]
                                                original_bp_list=Full_Info[8]
                                                Be_BP_Letter={}
                                                for let_key in original_letters:
                                                    Be_BP_Letter[let_key]=original_bp_list[original_letters.index(let_key)+1]-original_bp_list[original_letters.index(let_key)]
                                                ori_let2=[]
                                                for i in original_letters:
                                                    ori_let2.append(i)
                                                for i in original_letters:
                                                    if Copy_num_estimate[i]<1:
                                                        ori_let2.remove(i)
                                                    elif Copy_num_estimate[i]>3:
                                                        letter_copy=int(Copy_num_estimate[i]/2)
                                                        for j in range(letter_copy)[1:]:
                                                            ori_let2.append(i)
                                                ori_bp2=[original_bp_list[0]]
                                                for i in ori_let2:
                                                    ori_bp2.append(ori_bp2[-1]+Be_BP_Letter[i])
                                                Initial_TB=0
                                                Initial_Move_Prob=[1.0/3,1.0/3,1.0/3]
                                                Pair_Through=Full_Info[4]
                                                Read_Through=Full_Info[5]
                                                SingleR_Through=Full_Info[6]
                                                bp_MP=[original_bp_list,original_bp_list]
                                                letter_MP=[original_letters,original_letters]
                                                Be_BP=[ori_bp2,ori_bp2]
                                                Be_BP_Letter['left']=flank
                                                Be_BP_Letter['right']=flank
                                                for let_key in Be_BP_Letter.keys():
                                                    Be_BP_Letter[let_key+'^']=Be_BP_Letter[let_key]
                                                num_of_read_pairs=1
                                                for k1 in Be_BP_Letter.keys():
                                                        if not k1[-1]=='^' and not k1 in ['left','right']:
                                                                num_of_read_pairs+=Be_BP_Letter[k1]*RD_within_B[k1]/2/ReadLength
                                                num_of_read_pairs+=len(Full_Info[4])+len(Full_Info[5])+len(Full_Info[6])
                                                Be_Info=[Pair_Through,Read_Through,SingleR_Through]
                                                ori_let2=ori_let_Modi(Be_Info,ori_let2)
                                                Be_Letter=[ori_let2,ori_let2]
                                                Move_Step=0
                                                best_iterations=0
                                                Best_Score=float("-inf")
                                                Best_Letter=[]
                                                Best_BPs=[]
                                                score_record=[]
                                                total_read_Pairs_Num=(original_bp_list[-1]-original_bp_list[0])*GC_Mean_Coverage[Chr]/2/ReadLength
                                                num_of_reads=total_read_Pairs_Num
                                                Best_Score_Rec=0
                                                Score_rec_hash={}
                                                break_Iteration_Flag=0
                                                run_flag=0
                                                Best_Letter_Rec=[]
                                                if len(Full_Info[9])==1:
                                                    if Full_Info[1]['a']<GC_Mean_Coverage[Chr]/4 and Full_Info[2]<3:
                                                        Run_Result=zero_RD_Process(run_flag)
                                                        Best_Letter_Rec=Run_Result[0]
                                                        Best_Score_Rec=Run_Result[1]
                                                        run_flag=Run_Result[2]
                                                        Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                    else:
                                                        if Full_Info[1]['a']<GC_Mean_Coverage[Chr]:
                                                            Run_Result=one_RD_Process(run_flag,Score_rec_hash)
                                                            Best_Letter_Rec=Run_Result[0]
                                                            Best_Score_Rec=Run_Result[1]
                                                            run_flag=Run_Result[2]
                                                            Score_rec_hash=Run_Result[3]
                                                        else:
                                                            if Full_Info[1]['a']<2*GC_Mean_Coverage[Chr]:
                                                                Run_Result=two_RD_Process(run_flag,Score_rec_hash)
                                                                Best_Letter_Rec=Run_Result[0]
                                                                Best_Score_Rec=Run_Result[1]
                                                                run_flag=Run_Result[2]
                                                                Score_rec_hash=Run_Result[3]
                                                            else:
                                                                copy_num_a=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))
                                                                copy_num_b=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))+1
                                                                if copy_num_b<6:
                                                                    Run_Result=few_RD_Process(run_flag,Score_rec_hash)
                                                                    Best_Letter_Rec=Run_Result[0]
                                                                    Best_Score_Rec=Run_Result[1]
                                                                    run_flag=Run_Result[2]
                                                                    Score_rec_hash=Run_Result[3]
                                                                else:
                                                                    Run_Result=many_RD_Process(run_flag)
                                                                    Best_Letter_Rec=Run_Result[0]
                                                                    Best_Score_Rec=Run_Result[1]
                                                                    run_flag=Run_Result[2]
                                                                    Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                elif len(Full_Info[9])==2:
                                                    bl2_flag=0
                                                    for keyCNE in Copy_num_estimate.keys():
                                                        if not Copy_num_estimate[keyCNE]<2:
                                                            bl2_flag+=1
                                                    if bl2_flag==0:
                                                        Run_Result=two_block_RD_Process(run_flag)
                                                        Best_Letter_Rec=Run_Result[0]
                                                        Best_Score_Rec=Run_Result[1]
                                                        run_flag=Run_Result[2]
                                                        Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                if run_flag==0:
                                                    speed_test=10
                                                    t1_sptest=time.time()
                                                    while True:
                                                        if Move_Step>speed_test: break
                                                        Move_Step+=1
                                                        M_Key_0='_'.join(['Step',str(Move_Step),'M'])
                                                        P_Key_0='_'.join(['Step',str(Move_Step),'P'])
                                                        Move_Sample_Pool=['delete','invert','insert']
                                                        Move_M_P=Move_Choose(Move_Sample_Pool,Ploidy,Initial_Move_Prob)
                                                        M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
                                                        P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
                                                        if M_Move_Choices=='ERROR!' or P_Move_Choices=='ERROR!':
                                                                Move_Step-=1
                                                                continue
                                                        if not M_Move_Choices=='ERROR!' and not P_Move_Choices=='ERROR!':
                                                                IL_RD_Temp_Info=Af_Rearrange_Info_Collect_M(M_Move_Choices)
                                                                if IL_RD_Temp_Info=='Error': continue
                                                                ILTemp=IL_RD_Temp_Info[0]
                                                                RDTemp=IL_RD_Temp_Info[1]
                                                                Letter_Rec=IL_RD_Temp_Info[2]
                                                                BP_Rec=IL_RD_Temp_Info[3]
                                                                Regu_IL=ILTemp
                                                                Regu_RD=RDTemp
                                                                DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                                                                DECISION=DECISION_Score[0]
                                                                S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
                                                                Be_Letter=Letter_Rec[DECISION]
                                                                Be_BP=BP_Rec[DECISION]
                                                                if not S_DECISION in Score_rec_hash.keys():
                                                                    Score_rec_hash[S_DECISION]=[]
                                                                Score_rec_hash[S_DECISION].append(Be_Letter)
                                                                if S_DECISION>Best_Score:
                                                                        Best_Letter=[Be_Letter]
                                                                        Best_BPs=[Be_BP]
                                                                        Best_Score=S_DECISION
                                                                        best_iterations=0
                                                                elif S_DECISION==Best_Score:
                                                                        if not Be_Letter in Best_Letter:
                                                                                Best_Letter+=[Be_Letter]
                                                                                Best_BPs+=[Be_BP]
                                                                        best_iterations+=1
                                                                else:
                                                                    best_iterations+=1
                                                                score_record.append(S_DECISION)
                                                                IL_RD_Temp_Info=Af_Rearrange_Info_Collect_P(P_Move_Choices)
                                                                if IL_RD_Temp_Info=='Error': continue
                                                                ILTemp=IL_RD_Temp_Info[0]
                                                                RDTemp=IL_RD_Temp_Info[1]
                                                                Letter_Rec=IL_RD_Temp_Info[2]
                                                                BP_Rec=IL_RD_Temp_Info[3]
                                                                Regu_IL=ILTemp
                                                                Regu_RD=RDTemp
                                                                DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                                                                DECISION=DECISION_Score[0]
                                                                S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
                                                                Be_Letter=Letter_Rec[DECISION]
                                                                Be_BP=BP_Rec[DECISION]
                                                                if not S_DECISION in Score_rec_hash.keys():
                                                                    Score_rec_hash[S_DECISION]=[]
                                                                Score_rec_hash[S_DECISION].append(Be_Letter)
                                                                if S_DECISION>Best_Score:
                                                                        Best_Letter=[Be_Letter]
                                                                        Best_BPs=[Be_BP]
                                                                        Best_Score=S_DECISION
                                                                        best_iterations=0
                                                                elif S_DECISION==Best_Score:
                                                                        if not Be_Letter in Best_Letter:
                                                                                Best_Letter+=[Be_Letter]
                                                                                Best_BPs+=[Be_BP]
                                                                        best_iterations+=1
                                                                else:
                                                                        best_iterations+=1
                                                                score_record.append(S_DECISION)
                                                    t2_sptest=time.time()
                                                    if t2_sptest-t1_sptest<10 or bpsk1<4:
                                                        while True:
                                                                if Move_Step>Trail_Number: break
                                                                if best_iterations>100: 
                                                                        if Best_Score_Rec==0:
                                                                                best_iterations=0
                                                                                Best_Score_Rec=Best_Score
                                                                                Best_Letter_Rec=Best_Letter
                                                                                Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                                Best_BPs_Rec=Best_BPs
                                                                                Be_Letter=Best_Letter[0]
                                                                                Be_BP=Best_BPs[0]
                                                                                Best_Score-=100
                                                                        else:
                                                                                if Best_Score<Best_Score_Rec:
                                                                                    break_Iteration_Flag=1
                                                                                elif Best_Score==Best_Score_Rec:
                                                                                    break_Iteration_Flag=1
                                                                                    for i in Best_Letter:
                                                                                        if not i in Best_Letter_Rec:
                                                                                            Best_Letter_Rec.append(i)
                                                                                else:
                                                                                        best_iterations=0
                                                                                        Best_Score_Rec=Best_Score
                                                                                        Best_Letter_Rec=Best_Letter
                                                                                        Best_BPs_Rec=Best_BPs
                                                                                        Be_Letter=Best_Letter[0]
                                                                                        Be_BP=Best_BPs[0]
                                                                                        Best_Score-=100
                                                                if break_Iteration_Flag>0:
                                                                        break
                                                                Move_Step+=1
                                                                M_Key_0='_'.join(['Step',str(Move_Step),'M'])
                                                                P_Key_0='_'.join(['Step',str(Move_Step),'P'])
                                                                Move_Sample_Pool=['delete','invert','insert']
                                                                Move_M_P=Move_Choose(Move_Sample_Pool,Ploidy,Initial_Move_Prob)
                                                                M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
                                                                P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
                                                                if M_Move_Choices=='ERROR!' or P_Move_Choices=='ERROR!':
                                                                        Move_Step-=1
                                                                        continue
                                                                if not M_Move_Choices=='ERROR!' and not P_Move_Choices=='ERROR!':
                                                                        IL_RD_Temp_Info=Af_Rearrange_Info_Collect_M(M_Move_Choices)
                                                                        if IL_RD_Temp_Info=='Error': continue
                                                                        ILTemp=IL_RD_Temp_Info[0]
                                                                        RDTemp=IL_RD_Temp_Info[1]
                                                                        Letter_Rec=IL_RD_Temp_Info[2]
                                                                        BP_Rec=IL_RD_Temp_Info[3]
                                                                        Regu_IL=ILTemp
                                                                        Regu_RD=RDTemp
                                                                        DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                                                                        DECISION=DECISION_Score[0]
                                                                        S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
                                                                        Be_Letter=Letter_Rec[DECISION]
                                                                        Be_BP=BP_Rec[DECISION]
                                                                        if not S_DECISION in Score_rec_hash.keys():
                                                                            Score_rec_hash[S_DECISION]=[Be_Letter]
                                                                        else: 
                                                                            if not Be_Letter in Score_rec_hash[S_DECISION]:
                                                                                Score_rec_hash[S_DECISION].append(Be_Letter)
                                                                        if S_DECISION>Best_Score:
                                                                                Best_Letter=[Be_Letter]
                                                                                Best_BPs=[Be_BP]
                                                                                Best_Score=S_DECISION
                                                                                best_iterations=0
                                                                        elif S_DECISION==Best_Score:
                                                                                if not Be_Letter in Best_Letter:
                                                                                        Best_Letter+=[Be_Letter]
                                                                                        Best_BPs+=[Be_BP]
                                                                                best_iterations+=1
                                                                        else:
                                                                                best_iterations+=1
                                                                        score_record.append(S_DECISION)
                                                                        IL_RD_Temp_Info=Af_Rearrange_Info_Collect_P(P_Move_Choices)
                                                                        if IL_RD_Temp_Info=='Error': continue
                                                                        ILTemp=IL_RD_Temp_Info[0]
                                                                        RDTemp=IL_RD_Temp_Info[1]
                                                                        Letter_Rec=IL_RD_Temp_Info[2]
                                                                        BP_Rec=IL_RD_Temp_Info[3]
                                                                        Regu_IL=ILTemp
                                                                        Regu_RD=RDTemp
                                                                        DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Overall_Median_Num,GC_Var_Coverage)
                                                                        DECISION=DECISION_Score[0]
                                                                        S_DECISION=Regu_IL[DECISION]+Regu_RD[DECISION]
                                                                        Be_Letter=Letter_Rec[DECISION]
                                                                        Be_BP=BP_Rec[DECISION]
                                                                        if not S_DECISION in Score_rec_hash.keys():
                                                                            Score_rec_hash[S_DECISION]=[Be_Letter]
                                                                        else: 
                                                                            if not Be_Letter in Score_rec_hash[S_DECISION]:
                                                                                Score_rec_hash[S_DECISION].append(Be_Letter)
                                                                        if S_DECISION>Best_Score:
                                                                                Best_Letter=[Be_Letter]
                                                                                Best_BPs=[Be_BP]
                                                                                Best_Score=S_DECISION
                                                                                best_iterations=0
                                                                        elif S_DECISION==Best_Score:
                                                                                if not Be_Letter in Best_Letter:
                                                                                        Best_Letter+=[Be_Letter]
                                                                                        Best_BPs+=[Be_BP]
                                                                                best_iterations+=1
                                                                        else:
                                                                                best_iterations+=1
                                                                        score_record.append(S_DECISION)
                                                    else:
                                                        gaps=[]
                                                        bps2_new=[]
                                                        for k1 in bps2:
                                                            gaps.append([])
                                                            for k2 in range(len(k1)-2):
                                                                gaps[-1].append(int(k1[k2+2])-int(k1[k2+1]))
                                                        for k1 in range(len(gaps)):
                                                            bps2_new.append([])
                                                            chr_rec=bps2[k1][0]
                                                            rec1=1
                                                            for k2 in range(len(gaps[k1])):
                                                                if gaps[k1][k2]==max(gaps[k1]):
                                                                    bps2_new[-1].append([chr_rec]+bps2[k1][rec1:(k2+2)])
                                                                    bps2_new[-1].append([chr_rec]+bps2[k1][(k2+1):(k2+3)])
                                                                    rec1=k2+2
                                                            bps2_new[-1].append([chr_rec]+bps2[k1][rec1:])
                                                        for k1 in bps2_new:
                                                            for k2 in k1:
                                                                bps_hash[max(bps_hash.keys())].append([k2])
                                                        Best_Letter_Rec=[]
                                                        Best_Score_Rec=100
                                                struc_to_remove=[]
                                                for bestletter in Best_Letter_Rec:
                                                    if '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])==original_structure:
                                                        struc_to_remove.append(bestletter)
                                                Best_Letter_Rec=[i for i in Best_Letter_Rec if not i in struc_to_remove]
                                                if Best_Letter_Rec==[] and Best_Score_Rec==100:
                                                    continue
                                                else:
                                                    write_best_letter(bps2,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters)
                                            else:
                                                Score_rec_hash={}
                                                bps_new={}
                                                Initial_DR=Full_Info[2]
                                                Initial_IL=Full_Info[3]
                                                BlockGC=Full_Info[7]
                                                BlockGC['left']=0.476
                                                BlockGC['right']=0.476
                                                BlockGC2={}
                                                for key_B_GC in BlockGC.keys():
                                                    BlockGC2[key_B_GC]=BlockGC[key_B_GC]
                                                    BlockGC2[key_B_GC+'^']=BlockGC[key_B_GC]    
                                                original_letters=Full_Info[9]
                                                original_bp_list=Full_Info[8]
                                                for bl in Copy_num_Check:
                                                    for blk1 in chr_letter_bp.keys():
                                                        for blk2 in sorted(chr_letter_bp[blk1].keys()):
                                                            if blk2==bl:
                                                                bps2_temp=[blk1]+[chr_letter_bp[blk1][blk2][0],chr_letter_bp[blk1][blk2][-1]]
                                                                copy_num_a=int(Copy_num_estimate[bl])/2
                                                                copy_num_b=Copy_num_estimate[bl]-copy_num_a
                                                                Best_Letter_Rec=[[['a' for i in range(copy_num_a)],['a' for i in range(copy_num_a)]]]
                                                                Best_Score_Rec=100
                                                                write_best_letter([bps2_temp],Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters)
                                                for blk1 in chr_letter_bp.keys():
                                                    bps_new[blk1]=[]
                                                    for blk2 in sorted(chr_letter_bp[blk1].keys()):
                                                        if not blk2 in Copy_num_Check:
                                                            bps_new[blk1].append([chr_letter_bp[blk1][blk2][0],chr_letter_bp[blk1][blk2][-1]])
                                                bps_new_2=[]
                                                for k1 in bps_new.keys():
                                                    for k2 in bps_new[k1]:
                                                        if bps_new_2==[]:
                                                            bps_new_2.append([k1]+k2)
                                                        else:
                                                            if k1==bps_new_2[-1][0] and k2[0]==bps_new_2[-1][-1]:
                                                                bps_new_2[-1]+=k2[1:]
                                                            else:
                                                                bps_new_2.append([k1]+k2)
                                                for k1 in bps_new_2:
                                                    bps_hash[max(bps_hash.keys())].append([k1])
                            time2=time.time()
                            print 'SVPredict Complete !'
                            print 'Time Consuming: '+str(time2-time1)
    if function_name=='SVIntegrate':
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
        def let_to_bp_list(letter_list,bp_hash):
            out=[]
            if letter_list==[]:
                return []
            else:
                for k1 in letter_list:
                    out.append([])
                    for k2 in k1:
                        if not k2=='^':
                            if out[-1]==[]:
                                out[-1]+=[bp_hash[k2][0],bp_hash[k2][2],bp_hash[k2][4]]
                            else:
                                if bp_hash[k2][0]==out[-1][0]:
                                    out[-1]+=[bp_hash[k2][2],bp_hash[k2][4]]
                                else:
                                    out.append([bp_hash[k2][0],bp_hash[k2][2],bp_hash[k2][4]])
                out2=[]
                for k1 in out:
                    out2.append([k1[0],k1[1],k1[-1]])
                return out2
        def let_to_lets(letter_list):
            out=[]
            for k1 in letter_list:
                if k1==[]:
                    out.append(k1)
                else:
                    k2=[]
                    if not '^' in k1[0]:
                        for k3 in k1:
                            if k2==[]:
                                k2.append(k3)
                            else:
                                if ord(k3)-ord(k2[-1][-1])==1:
                                    k2[-1]+=k3
                                else:
                                    k2.append(k3)
                    else:
                        for k3 in k1:
                            if k2==[]:
                                k2.append(k3)
                            else:
                                if -ord(k3[0])+ord(k2[-1][0])==1:
                                    k2[-1]=k3[0]+k2[-1]
                                else:
                                    k2.append(k3)
                    out.append(k2)
            return out
        def read_in_structures(filein):
            fin=open(filein)
            while True:
                pin1=fin.readline().strip().split()
                if not pin1: break
                if pin1[0]=='Total': break
                pin2=fin.readline().strip().split()
                pin3=fin.readline().strip().split()
                pin4=fin.readline().strip().split()
                pin5=fin.readline().strip().split()
                if pin3[0]=='Theoretical' and pin4[0]=='Current' and pin5[0]=='Time':
                    let1=bp_to_let([pin1])
                    if not let1==0:
                        let2='/'.join(sorted(pin2[0].split('/')))
                        if not let1 in sv_info.keys():
                            sv_info[let1]={}
                        if not let2 in sv_info[let1].keys():
                            sv_info[let1][let2]=[]
                        if not pin1 in sv_info[let1][let2]:
                            sv_info[let1][let2].append(pin1+[float(pin4[-1])-float(pin3[-1])])
            fin.close()
        def let_recombine(let_list):
            out=[]
            if not let_list==[]:
                out=[let_list[0]]
                for i in let_list[1:]:
                    if ord(i)-ord(out[-1][-1])==1:
                        out[-1]+=i
                    else:
                        out.append(i)
            return out
        def inv_flag(k1,k2):
            if '^' in k2:
                return 1
            else:
                return 0
        def inv_flag_SA(k1,k2):
            out=0
            if '^' in k2:
                if k2.replace('^','')==k1:
                    out+=1
            return out
        def unit_produce(list):
            temp1=[sorted(list)[0]]
            for k1 in sorted(list)[1:]:
                if ord(k1)-ord(temp1[-1][-1])==1:
                    temp1[-1]+=k1
                else:
                    temp1.append(k1)
            temp2=[]
            for k1 in temp1:
                for k2 in range(len(k1)+1)[1:]:
                    for k3 in range(len(k1)-k2+1):
                        temp2.append(k1[k3:(k3+k2)])
            return temp2[::-1]
        def dup_flag(k1,k2):
            test=0
            for i in k2.split('/')[0]:
                if k2.split('/')[0].count(i)>1:
                    test=1
            for i in k2.split('/')[1]:
                if k2.split('/')[1].count(i)>1:
                    test=1
            return test
        def dup_flag_SA(k1,k2):
            test=0
            for i in k2:
                if k2.count(i)>1:
                    if k2.replace(''.join([i for j in range(k2.count(i))]),i)==k1:
                        test=1
            return test
        def tra_flag(k1,k2):
            test=0
            k2a=[i for i in k2.split('/')[0]]
            k2b=[i for i in k2.split('/')[1]]
            for i in range(len(k2a)-1):
                if ord(k2a[i+1])-ord(k2a[i])<0:
                    test=1
            for i in range(len(k2b)-1):
                if ord(k2b[i+1])-ord(k2b[i])<0:
                    test=1
            return test
        def tra_flag_SA(k1,k2):
            test=0
            k2a=[i for i in k2]
            for i in range(len(k2a)-1):
                if ord(k2a[i+1])-ord(k2a[i])<0:
                    test=1
            return test
        def del_flag(k1,k2):
            if inv_flag(k1,k2)+dup_flag(k1,k2)+tra_flag(k1,k2)==0:
                return 1
            else:
                return 0
        def del_flag_SA(k1,k2):
            out=0
            if not '^' in k2:
                flagdup=0
                for i in k2:
                    if k2.count(i)>2:
                        flagdup+=1
                if flagdup==0:
                    flagtra=0
                    for i in range(len(k2)-1):
                        if ord(k2[i+1])-ord(k2[i])<1:
                            flagtra+=1
                    if flagtra==0:
                        if not k1==k2:
                            out=1   
            return out
        def allele_DUP_decide(x1):
            flag_Dup=0
            for x2 in x1:
                if x1.count(x2)>1:
                    flag_Dup+=1
            if flag_Dup==0:
                return 'NO_DUP'
            else:
                return 'DUP'
        def allele_TRA_decide(x1):
            flag_TRA=0
            if len(x1)<2:
                return 'NO_TRA'
            else:
                for x2 in range(len(x1)-1):
                    if ord(x1[x2+1])-ord(x1[x2])<0:
                        flag_TRA+=1
                if flag_TRA==0:
                    return 'NO_TRA'
                else:
                    return 'TRA'
        def simple_DEL_decide(k1,k2):
            out='Not_Simple_DEL' 
            if not k2==k1:
                if not '^' in k2:
                    flag_Dup=0
                    for x1 in k2.split('/'):
                        if allele_DUP_decide(x1)=='DUP':
                            flag_Dup+=1
                    if flag_Dup==0:
                        flag_TRA=0
                        for x1 in k2.split('/'):
                            if allele_TRA_decide(x1)=='TRA':
                                flag_TRA+=1
                        if flag_TRA==0:
                            out='Simple_DEL'
            return out
        def simple_DEL_type(k1,k2):
            if simple_DEL_decide(k1,k2)=='Simple_DEL':
                del_part=[[],[]]
                for x1 in k1.split('/')[0]:
                    if not x1 in k2.split('/')[0]:
                        del_part[0].append(x1)
                for x1 in k1.split('/')[1]:
                    if not x1 in k2.split('/')[1]:
                        del_part[1].append(x1)
                return del_part
            else:
                return 'Error'
        def simple_INV_decide(k1,k2):
            out='Not_Simple_INV' 
            if not k2==k1:
                if '^' in k2:
                    k2_new=letter_recombine(k2)
                    if '/'.join([''.join([i.replace('^','') for i in k2_new[0]]),''.join([i.replace('^','') for i in k2_new[1]])])==k1:
                        out='Simple_INV'
            return out
        def letter_separate(k2):
            out=[[],[]]
            for x1 in k2.split('/')[0]:
                if not x1=='^':
                    out[0].append(x1)
                else:
                    out[0][-1]+=x1
            for x1 in k2.split('/')[1]:
                if not x1=='^':
                    out[1].append(x1)
                else:
                    out[1][-1]+=x1
            return out
        def letter_recombine(k2):
            k2_new=letter_separate(k2)
            out=[[],[]]
            if k2_new[0]==[]:
                out=[[],[k2_new[1][0]]]
            if k2_new[1]==[]:
                out=[[k2_new[0][0]],[]]
            if not k2_new[0]==[] and not k2_new[1]==[]:
                out=[[k2_new[0][0]],[k2_new[1][0]]]
            for x1 in k2_new[0]:
                if not '^' in x1 and not '^' in out[0][-1]:
                    if ord(x1)-ord(out[0][-1][-1])==1:
                        out[0][-1]+=x1
                    else:
                        out[0].append(x1)
                elif '^' in x1 and '^' in out[0][-1]:
                    if ord(x1[0])-ord(out[0][-1][0])==-1:
                        out[0][-1]=x1[0]+out[0][-1]
                    else:
                        out[0].append(x1)
                else:
                    out[0].append(x1)
            for x1 in k2_new[1]:
                if not '^' in x1 and not '^' in out[1][-1]:
                    if ord(x1)-ord(out[1][-1][-1])==1:
                        out[1][-1]+=x1
                    else:
                        out[1].append(x1)
                elif '^' in x1 and '^' in out[1][-1]:
                    if ord(x1[1])-ord(out[1][-1][1])==-1:
                        out[1][-1]=x1[1]+out[1][-1]
                    else:
                        out[1].append(x1)
                else:
                    out[1].append(x1)
            if not out[0]==[]:
                del out[0][0]
            if not out[1]==[]:
                del out[1][0]
            return out
        def simple_INV_type(k1,k2):
            k2_new=letter_separate(k2)
            if simple_INV_decide(k1,k2)=='Simple_INV':
                out=[[],[]]
                for x1 in k2_new[0]:
                    if '^' in x1:
                        out[0].append(x1)
                for x1 in k2_new[1]:
                    if '^' in x1:
                        out[1].append(x1)
                out2=[[x1.replace('^','') for x1 in out[0]],[x1.replace('^','') for x1 in out[1]]]
                return out2
            else:
                return 'Error'
        def simple_DUP_decide(k1,k2):
            out='Not_Simple_DUP' 
            if not k2==k1:
                if not '^' in k2:
                    k3=[[],[]]
                    for x1 in k2.split('/')[0]:
                        if k3[0]==[]:
                            k3[0].append(x1)
                        else:
                            if not x1 in k3[0][0]:
                                k3[0][0]+=x1
                            else:
                                if not x1==k3[0][0][-1]:
                                    index1=k3[0][0].index(x1)-k3[0][0].index(k3[0][0][-1])-(ord(x1)-ord(k3[0][0][-1]))
                                    index2=1-(ord(x1)-ord(k3[0][0][-1]))
                                    if abs(index1)>abs(index2): 
                                        k3[0][0]=k3[0][0].replace(x1,'')
                                        k3[0][0]+=x1
                    for x1 in k2.split('/')[1]:
                        if k3[1]==[]:
                            k3[1].append(x1)
                        else:
                            if not x1 in k3[1][0]:
                                k3[1][0]+=x1
                            else:
                                if not x1==k3[1][0][-1]:
                                    index1=k3[1][0].index(x1)-k3[1][0].index(k3[1][0][-1])-(ord(x1)-ord(k3[1][0][-1]))
                                    index2=1-(ord(x1)-ord(k3[1][0][-1]))
                                    if abs(index1)>abs(index2): 
                                        k3[1][0]=k3[1][0].replace(x1,'')
                                        k3[1][0]+=x1
                        if k3==letter_recombine(k1):
                            out='Simple_DUP'
            return out
        def index_element(k1,ele):
            out=[]
            for x in range(len(k1)-len(ele)+1):
                if k1[x:(x+len(ele))]==ele:
                    out.append(x)
            return out
        def simple_DUP_type(k1,k2):
            if simple_DUP_decide(k1,k2)=='Simple_DUP':
                k3=letter_recombine(k2)
                out=[[],[]]
                for x1 in k3[0]:
                    index_pos=index_element(k2.split('/')[0],x1)
                    if len(index_pos)>1:
                        tandem_flag=0
                        for y1 in range(len(index_pos)-1):
                            if index_pos[y1+1]-index_pos[y1]>len(x1):
                                tandem_flag+=1
                        if tandem_flag==0:
                            if not x1+'_Tandem' in out[0]:
                                out[0].append(x1+'_Tandem')
                        else:
                            if not x1+'_Disperse' in out[0]:
                                out[0].append(x1+'_Disperse')
                for x1 in k3[1]:
                    index_pos=index_element(k2.split('/')[1],x1)
                    if len(index_pos)>1:
                        tandem_flag=0
                        for y1 in range(len(index_pos)-1):
                            if index_pos[y1+1]-index_pos[y1]>len(x1):
                                tandem_flag+=1
                        if tandem_flag==0:
                            if not x1+'_Tandem' in out[1]:
                                out[1].append(x1+'_Tandem')
                        else:
                            if not x1+'_Disperse' in out[1]:
                                out[1].append(x1+'_Disperse')
                return out  
            else:
                return 'Error'
        def simple_TRA_decide(k1,k2):
            out='Not_Simple_TRA' 
            if '/'.join([''.join(sorted(k2.split('/')[0])),''.join(sorted(k2.split('/')[1]))])==k1:
                out='simple_TRA'
            return out
        def simple_flag_SA(k1,k2):
            if k2=='':
                return[[i for i in k1],[],[],0]
            else:
                temp=[]
                break_flag=0
                for i in k2:
                    if not i=='^':
                        temp.append(i)
                    else:
                        temp[-1]+=i
                temp2=[temp[0]]
                for i in range(len(temp[1:])):
                    if not '^' in temp[i] and not '^' in temp[i+1] and ord(temp[i+1])-ord(temp[i])==1:
                        temp2[-1]+=temp[i+1]
                    elif '^' in temp[i] and '^' in temp[i+1] and ord(temp[i+1][0])-ord(temp[i][0])==-1:
                        temp2[-1]=temp[i+1][0]+temp2[-1]
                    else:
                        temp2.append(temp[i+1]) 
                outdel=[]
                outinv=[]
                outdup=[]
                outtra=0
                for i in range(len(temp2)):
                    j=temp2[i]
                    if '^' in j:
                        if not j.replace('^','') in outinv:
                            outinv.append(j.replace('^',''))
                        temp2[i]=j.replace('^','')
                temp3=''.join(temp2)
                for i in range(len(temp3)-1):
                    if not temp3[i+1]=='/' and not temp3[i]=='/':
                        if ord(temp3[i+1])-ord(temp3[i])<0:
                            outtra=1
                if not temp3==k1:
                    temp4=[]
                    for i in temp3:
                        if temp3.count(i)>1:
                            if not i in outdup:
                                outdup.append(i)
                        if not i in temp4:
                            temp4.append(i)
                    if not ''.join(temp4)==k1:
                        for i in k1:
                            if not i in temp4:
                                outdel.append(i)
                if not outdup==[]:
                    dupuni=unit_produce(outdup)
                    outdup2=[]
                    k3=k2
                    for i in dupuni:
                        ia=i
                        ib=''.join([j+'^' for j in i[::-1]])
                        if len(i)>1:
                            if temp2.count(ia)+temp2.count(ib)>2:
                                outdup2.append([i,temp2.count(ia)+temp2.count(ib)])
                                k3=k3.replace(ia,'')
                                k3=k3.replace(ib,'')
                        elif len(i)==1:
                            if k3.count(ia)+k3.count(ib)>1:
                                outdup2.append([i,k3.count(ia)])
                                k3=k3.replace(ia,'')
                                k3=k3.replace(ib,'')
                else:
                    outdup2=[]
                outdup3=[]
                for i in outdup2:
                    flag=0
                    for j in outdup2:
                        if not j==i:
                            if len(j[0])>len(i[0]) and i[0] in j[0]:
                                flag+=1
                    if flag==0:
                        outdup3.append(i)       
                if len(outdup3)>1:
                    outdup3.sort()
                    outdup4=[outdup3[0]]
                    for i in range(len(outdup3)-1):
                        if outdup3[i+1][-1]==outdup3[i][-1] and ord(outdup3[i+1][0][0])-ord(outdup3[i][0][-1])==1:
                            outdup4[-1][0]+=outdup3[i+1][0]
                        else:
                            outdup4.append(outdup3[i+1])
                else:
                    outdup4=outdup3
                return [outdel,outinv,outdup4,outtra]
        def sv_rec_2(sv_info):
            for k1ab in sv_info.keys():
                for k2ab in sv_info[k1ab].keys():
                    if not k2ab==k1ab:
                        k1aba=k1ab.split('/')[0]
                        k2aba=k2ab.split('/')[0]
                        k2abb=k2ab.split('/')[1]
                        flaga=[]
                        flagb=[]
                        test=[[],[]]
                        if flaga==[] and not k1aba==k2aba:
                            if k2aba=='':
                                csv1=[[i for i in k1aba],[],[],0]
                            else:
                                csv1=simple_flag_SA(k1aba,k2aba)
                            add_csv_info(csv1,1,k1ab,k2ab)
                        if flagb==[] and not k1aba==k2abb:
                            if k2abb=='':
                                csv1=[[i for i in k2abb],[],[],0]
                            else:
                                csv1=simple_flag_SA(k1aba,k2abb)
                            add_csv_info(csv1,2,k1ab,k2ab)
        def sv_rec(sv_info):
            for k1ab in sv_info.keys():
                for k2ab in sv_info[k1ab].keys():
                    if not k2ab==k1ab:
                        if del_flag(k1ab,k2ab)==1:
                            delM=[]
                            delP=[]
                            for i in k1ab.split('/')[0]:
                                if not i in k2ab.split('/')[0]:
                                    delM.append(i)
                                if not i in k2ab.split('/')[1]:
                                    delP.append(i)
                            for k3 in sv_info[k1ab][k2ab]:
                                del_info_add(k3,[delM,delP])
                        else:
                            if inv_flag(k1ab,k2ab)+dup_flag(k1ab,k2ab)==0:
                                tra_info_add(k1ab,k2ab)
                                if del_flag_SA(k1ab.split('/')[0],k2ab.split('/')[0])==1:
                                    delM=[]
                                    delP=[]
                                    for i in k1ab.split('/')[0]:
                                        if not i in k2ab.split('/')[0]:
                                            delM.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        del_info_add(k3,[delM,delP])
                                if del_flag_SA(k1ab.split('/')[1],k2ab.split('/')[1])==1:
                                    delM=[]
                                    delP=[]
                                    for i in k1ab.split('/')[0]:
                                        if not i in k2ab.split('/')[1]:
                                            delP.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        del_info_add(k3,[delM,delP])
                            else:
                                k1aba=k1ab.split('/')[0]
                                k2aba=k2ab.split('/')[0]
                                k2abb=k2ab.split('/')[1]
                                flaga=[]
                                flagb=[]
                                if del_flag_SA(k1aba,k2aba)==1:#simple del on one allele
                                    delM=[]
                                    delP=[]
                                    for i in k1ab.split('/')[0]:
                                        if not i in k2ab.split('/')[0]:
                                            delM.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        del_info_add(k3,[delM,delP])
                                    flaga.append('del')
                                if del_flag_SA(k1aba,k2abb)==1:#simple del on one allele
                                    delM=[]
                                    delP=[]
                                    for i in k1ab.split('/')[0]:
                                        if not i in k2ab.split('/')[1]:
                                            delP.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        del_info_add(k3,[delM,delP])
                                    flagb.append('del')
                                if dup_flag_SA(k1aba,k2aba)==1:#simple dup on one allele
                                    dupM=[]
                                    dupP=[]
                                    for i in k1aba:
                                        if k2aba.count(i)>1:
                                            dupM.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        dup_info_add(k3,[dupM,dupP])
                                    flaga.append('dup')
                                if dup_flag_SA(k1aba,k2abb)==1:#simple dup on one allele
                                    dupM=[]
                                    dupP=[]
                                    for i in k1aba:
                                        if k2abb.count(i)>1:
                                            dupP.append(i)
                                    for k3 in sv_info[k1ab][k2ab]:
                                        dup_info_add(k3,[dupM,dupP])
                                    flagb.append('dup')
                                if inv_flag_SA(k1aba,k2aba)==1:#simple inv on one allele
                                    invM=[]
                                    invP=[]
                                    for i in range(len(k2aba)):
                                        if k2aba[i]=='^':
                                            invM.append(k2aba[i-1])
                                    for k3 in sv_info[k1ab][k2ab]:
                                        inv_info_add(k3,[invM,invP])
                                    flaga.append('inv')
                                if inv_flag_SA(k1aba,k2abb)==1:#simple inv on one allele
                                    invM=[]
                                    invP=[]
                                    for i in range(len(k2abb)):
                                        if k2abb[i]=='^':
                                            invP.append(k2abb[i-1])
                                    for k3 in sv_info[k1ab][k2ab]:
                                        inv_info_add(k3,[invM,invP])
                                    flagb.append('inv')
                                if flaga==[] and not k1aba==k2aba:
                                    csv1=simple_flag_SA(k1aba,k2aba)
                                    add_csv_info(csv1,1,k1ab,k2ab)
                                if flagb==[] and not k1aba==k2abb:
                                    csv1=simple_flag_SA(k1aba,k2abb)
                                    add_csv_info(csv1,2,k1ab,k2ab)
        def dup_type_decide(dup_let,flag_sex,k1,k2):
            if flag_sex==1:
                k1x=k1.split('/')[0]
                k2x=k2.split('/')[0]
            else:
                k1x=k1.split('/')[1]
                k2x=k2.split('/')[1]
            out=[]     
            for x in dup_let:
                out.append([])
                for y in x:
                    pos_index=[z for z in range(len(k2x)) if k2x[z]==y[0]]
                    inter_index=[pos_index[z+1]-pos_index[z] for z in range(len(pos_index)-1)]
                    if 1 in inter_index:
                        out[-1].append('Tandem')
                        inter_index.remove(1)
                    if not inter_index==[]:
                        out[-1].append('Disperse')
            return out
        def add_csv_info(csv1,flag_sex,k1,k2):
            if flag_sex==1:
                del_let=[csv1[0],[]]
                inv_let=[csv1[1],[]]
                dup_let=[csv1[2],[]]
            else:
                del_let=[[],csv1[0]]
                inv_let=[[],csv1[1]]
                dup_let=[[],csv1[2]]
            if simple_DEL_decide(k1,k2)=='Simple_DEL':
                for k3 in sv_info[k1][k2]:
                    del_info_add(k3,del_let)
            elif simple_DUP_decide(k1,k2)=='Simple_DUP':
                dup_subtype=simple_DUP_type(k1,k2)
                dis_dup_let=[[],[]]
                tan_dup_let=[[],[]]
                for x in range(2):
                    if not dup_subtype[x]==[]:
                        for y in dup_subtype[x]:
                            if y.split('_')[1]=='Disperse':
                                dis_dup_let[x].append([y.split('_')[0],k2.split('/')[x].count(y.split('_')[0])])
                            else:
                                tan_dup_let[x].append([y.split('_')[0],k2.split('/')[x].count(y.split('_')[0])])
                if not tan_dup_let==[[],[]]:
                    for k3 in sv_info[k1][k2]:
                        dup_info_2_add(k3,tan_dup_let)
                if not dis_dup_let==[[],[]]:
                    for k3 in sv_info[k1][k2]:
                        disperse_dup_info_2_add(k3,dis_dup_let)
            elif simple_INV_decide(k1,k2)=='Simple_INV':
                for k3 in sv_info[k1][k2]:
                    inv_info_add(k3,inv_let)
            elif simple_TRA_decide(k1,k2)=='simple_TRA':
                tra_info_add(k1,k2)
            else:
                dup_csv_subtype=dup_type_decide(dup_let,flag_sex,k1,k2)
                for k3 in sv_info[k1][k2]:
                    del_csv_info_add(k3,del_let)
                    inv_csv_info_add(k3,inv_let)
                    dup_csv_info_add(k3,dup_let,dup_csv_subtype)
                if csv1[3]==1:
                    tra_csv_info_add(k1,k2)
        def bp_to_hash(bp_list,sv_let):
            bp_hash={}
            block_rec=0
            block_hash=[]
            sv_let=[i[0] for i in sv_let]
            for a3 in bp_list:
                if a3 in chromos or not a3.isdigit():
                    block_hash.append([a3])
                else:
                    block_hash[-1].append(a3)
            for a3 in block_hash:
                for a4 in range(len(a3)-2):
                    bp_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                    block_rec+=1
            out=[]
            if not sv_let==[]:
                if len(sv_let)==1:
                    out=[bp_hash[sv_let[0]]]
                else:
                    out.append(bp_hash[sv_let[0]])
                    for ka in range(len(sv_let)-1):
                        if ord(sv_let[ka+1])-ord(sv_let[ka])==1 and bp_hash[sv_let[ka+1]][0]==bp_hash[sv_let[ka]][0]:
                            out[-1]+=bp_hash[sv_let[ka+1]][1:]
                        else:
                            out.append(bp_hash[sv_let[ka+1]])
            out2=[]
            for ka in out:
                out2.append([ka[0],int(ka[1]),int(ka[-1])])
            return out2
        def del_info_reorganize(k1,k2):
            del_let=[[],[]]
            for k3 in k1.split('/')[0]:
                if not k3 in k2.split('/')[0]:
                    del_let[0].append(k3)
            for k3 in k1.split('/')[1]:
                if not k3 in k2.split('/')[1]:
                    del_let[1].append(k3)
            for k3 in sv_info[k1][k2]:
                del_bp=[]
                if not del_let[0]==[]:
                    del_bp.append(bp_to_hash(k3,del_let[0]))
                else:
                    del_bp.append([])
                if not del_let[1]==[]:
                    del_bp.append(bp_to_hash(k3,del_let[1]))
                else:
                    del_bp.append([])            
                if del_bp[0]==del_bp[1]:
                    for k4 in del_bp[0]:
                        if not k4[0] in del1.keys():
                            del1[k4[0]]=[]
                        if not [int(k4[1]),int(k4[-1]),'hom'] in del1[k4[0]]:
                            del1[k4[0]].append([int(k4[1]),int(k4[-1]),'hom'])
                else:
                    for k5 in del_bp:
                        for k4 in k5:
                            if not k4[0] in del1.keys():
                                del1[k4[0]]=[]
                            if not [int(k4[1]),int(k4[-1]),'het'] in del1[k4[0]]:
                                del1[k4[0]].append([int(k4[1]),int(k4[-1]),'het'])
        def del_info_add(k3,del_let):
            tempa=bp_to_hash(k3[:-1],del_let[0])
            tempb=bp_to_hash(k3[:-1],del_let[1])
            for k1 in tempa:
                if k1 in tempb:
                    tempc='hom'
                    tempb.remove(k1)
                else:
                    tempc='heta'
                if not k1[0] in del1.keys():
                    del1[k1[0]]=[]
                del1[k1[0]].append(k1[1:]+[tempc,k3[-1],'_'.join(k3[:-1]+['S'])])
            for k1 in tempb:
                if not k1[0] in del1.keys():
                    del1[k1[0]]=[]
                del1[k1[0]].append(k1[1:]+['hetb',k3[-1],'_'.join(k3[:-1]+['S'])])
        def del_csv_info_add(k3,del_let):
            tempa=bp_to_hash(k3[:-1],del_let[0])
            tempb=bp_to_hash(k3[:-1],del_let[1])
            for k1 in tempa:
                if k1 in tempb:
                    tempc='hom'
                    tempb.remove(k1)
                else:
                    tempc='heta'
                if not k1[0] in del1.keys():
                    del1[k1[0]]=[]
                del1[k1[0]].append(k1[1:]+[tempc,k3[-1],'_'.join(k3[:-1]+['C'])])
            for k1 in tempb:
                if not k1[0] in del1.keys():
                    del1[k1[0]]=[]
                del1[k1[0]].append(k1[1:]+['hetb',k3[-1],'_'.join(k3[:-1]+['C'])])
        def dup_info_add(k3,dup_let):
            for k2x in dup_let:
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4])
                    for k5 in temp:
                        if not k5[0] in dup1.keys():
                            dup1[k5[0]]=[]
                        dup1[k5[0]].append(k5[1:]+[k3[-1],'_'.join(k3[:-1])])
        def dup_info_2_add(k3,dup_let):
            temprec=-1
            for k2x in dup_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
                    for k5 in temp:
                        if not k5[0] in dup1.keys():
                            dup1[k5[0]]=[]
                        if k4[1]>1:
                            dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])
        def dup_csv_info_add(k3,dup_let,dup_csv_subtype):
            temprec=-1
            dup_index_1=-1
            for k2x in dup_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                dup_index_1+=1
                dup_index_2=-1
                for k4 in k2x:
                    dup_index_2+=1
                    dup_subtype_current=dup_csv_subtype[dup_index_1][dup_index_2]
                    if dup_subtype_current=='Tandem':
                        temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
                        for k5 in temp:
                            if not k5[0] in disperse_dup.keys():
                                disperse_dup[k5[0]]=[]
                            if k4[1]>1:
                                disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])
                    elif dup_subtype_current=='Disperse':
                        temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
                        for k5 in temp:
                            if not k5[0] in disperse_dup.keys():
                                disperse_dup[k5[0]]=[]
                            if k4[1]>1:
                                disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])
        def disperse_dup_info_2_add(k3,dup_let):
            temprec=-1
            for k2x in dup_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
                    for k5 in temp:
                        if not k5[0] in disperse_dup.keys():
                            disperse_dup[k5[0]]=[]
                        if k4[1]>1:
                            disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S']),k4[1]])
        def dup_csv_info_2_add(k3,dup_let):
            temprec=-1
            for k2x in dup_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4[0]])
                    for k5 in temp:
                        if not k5[0] in dup1.keys():
                            dup1[k5[0]]=[]
                        if k4[1]>1:
                            dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['C']),k4[1]])
        def inv_info_add(k3,inv_let):
            temprec=-1
            for k2x in inv_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4])
                    for k5 in temp:
                        if not k5[0] in inv1.keys():
                            inv1[k5[0]]=[]
                        inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['S'])])
        def inv_csv_info_add(k3,inv_let):
            temprec=-1
            for k2x in inv_let:
                temprec+=1
                hetx=['heta','hetb'][temprec]
                for k4 in k2x:
                    temp=bp_to_hash(k3[:-1],[i for i in k4])
                    for k5 in temp:
                        if not k5[0] in inv1.keys():
                            inv1[k5[0]]=[]
                        inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],'_'.join(k3[:-1]+['C'])])
        def tra_csv_info_add(k1,k2):
            for k3 in sv_info[k1][k2]:
                SV_ID='_'.join([str(i) for i in k3]+['C'])
                if not SV_ID in tra1.keys():
                    tra1[SV_ID]={}
                k2a=k2.split('/')[0]
                k2b=k2.split('/')[1]
                bp_hash={}
                block_rec=0
                block_hash=[]
                for a3 in k3[:-1]:
                    if a3 in chromos or not a3.isdigit():
                        block_hash.append([a3])
                    else:
                        block_hash[-1].append(a3)
                for a3 in block_hash:
                    for a4 in range(len(a3)-2):
                        bp_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                        block_rec+=1
                for a3 in bp_hash.keys():
                    temp=[]
                    for a4 in bp_hash[a3][1:]:
                        temp.append(int(a4)-1)
                        temp.append(int(a4))
                    bp_hash[a3][1:]=temp
                bp_hash['left']=[bp_hash[k1[0]][0],bp_hash[k1[0]][1],bp_hash[k1[0]][2]]
                bp_hash['right']=[bp_hash[k1[-1]][0],bp_hash[k1[-1]][3],bp_hash[k1[-1]][4]]
                ref_allele={}
                for a3 in bp_hash.keys():
                    ref_allele[a3]=[bp_hash[a3][0]]
                    for a4 in bp_hash[a3][1:]:
                        ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
                if not k2a==k1.split('/')[0] and del_flag_SA(k1.split('/')[0],k2a)==0:
                    flag1=0#flag1==0:w/o inversion in the alt structure
                    if '^' in k2a:
                        flag1+=1
                    flag2=0#flag2==0:w/o duplication in the alt structure
                    for j in k2a:
                        if k2a.count(j)>1:
                            flag2+=1
                    flag3=0 #flag3==0: w/o translocation
                    if len(k2a)>1:
                        for i in range(len(k2a)-1):
                            if not ord(k2a[i+1])>ord(k2a[i]):
                                flag3+=1
                    if flag1+flag2+flag3==0:
                        heta_Del_block=[]
                        for a1 in k1.split('/')[0]:
                            if not a1 in k2a:
                                heta_Del_block.append(a1)   
                        if not 'a' in tra1[SV_ID].keys():
                            tra1[SV_ID]['a']=[]
                        block_hash=[]
                        del_hash={}
                        block_rec=0
                        for a3 in a2[0]:
                            if a3 in chromos:
                                block_hash.append([a3])
                            else:
                                block_hash[-1].append(a3)
                        for a3 in block_hash:
                            for a4 in range(len(a3)-2):
                                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                                block_rec+=1                                
                        if not heta_Del_block==[]:
                            a_heta=0
                            heta_Del_new=[heta_Del_block[0]]
                            while True:
                                a_heta+=1
                                if a_heta==len(heta_Del_block):break
                                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
                                    heta_Del_new[-1]+=heta_Del_block[a_heta]
                                else:
                                    heta_Del_new.append(heta_Del_block[a_heta])
                            for a3 in heta_Del_new:
                                a4=a3[0]
                                tra1[SV_ID]['a'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
                                a4=a3[-1]
                                tra1[SV_ID]['a'][-1].append(int(del_hash[a4][2])-1)
                    else:
                        if not 'a' in tra1[SV_ID].keys():
                            tra1[SV_ID]['a']=[]
                        t1=[]
                        for a3 in k2a:
                            if not a3=='^':
                                t1.append(a3)
                            else:
                                t1[-1]+=a3
                        t2=[t1[0]]
                        for a3 in t1[1:]:
                            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
                                t2[-1]+=a3
                            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
                                t2[-1]+=a3
                            else:
                                t2.append(a3)
                        a3='left'
                        a4=t2[0]
                        l_chr=bp_hash[a3][0]
                        r_chr=bp_hash[a4[0]][0]
                        if not '^' in a4:
                            if not a4[0]==k1[0]:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
                        elif '^' in a4:
                            tra1[SV_ID]['a'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
                            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                        for t3 in range(len(t2)-1):
                            a3=t2[t3]
                            a4=t2[t3+1]
                            l_chr=bp_hash[a3[0]][0]
                            r_chr=bp_hash[a4[0]][0]
                            if not '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
                            elif '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
                            elif not '^' in a3 and '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                            elif '^' in a3 and '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
                        if len(t2)>1:
                            a3=t2[t3+1]
                        else:
                            a3=t2[0]
                        a4='right'
                        l_chr=bp_hash[a3[0]][0]
                        r_chr=bp_hash[a4][0]
                        if not '^' in a3:
                            if not a3[-1]==k1[-1]:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
                        if '^' in a3:
                            tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
                            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])
                if not k2b==k1.split('/')[1] and del_flag_SA(k1.split('/')[1],k2b)==0:
                    flag1=0#flag1==0:w/o inversion in the alt structure
                    if '^' in k2b:
                        flag1+=1
                    flag2=0#flag2==0:w/o duplication in the alt structure
                    for j in k2b:
                        if k2b.count(j)>1:
                            flag2+=1
                    flag3=0 #flag3==0: w/o translocation
                    if len(k2b)>1:
                        for i in range(len(k2b)-1):
                            if not ord(k2b[i+1])>ord(k2b[i]):
                                flag3+=1
                    if flag1+flag2+flag3==0:
                        heta_Del_block=[]
                        for a1 in k1.split('/')[1]:
                            if not a1 in k2b:
                                heta_Del_block.append(a1)
                        if not 'b' in tra1[SV_ID].keys():       
                            tra1[SV_ID]['b']=[]
                        block_hash=[]
                        del_hash={}
                        block_rec=0
                        for a3 in a2[0]:
                            if a3 in chromos:
                                block_hash.append([a3])
                            else:
                                block_hash[-1].append(a3)
                        for a3 in block_hash:
                            for a4 in range(len(a3)-2):
                                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                                block_rec+=1                                
                        if not heta_Del_block==[]:
                            a_heta=0
                            heta_Del_new=[heta_Del_block[0]]
                            while True:
                                a_heta+=1
                                if a_heta==len(heta_Del_block):break
                                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
                                    heta_Del_new[-1]+=heta_Del_block[a_heta]
                                else:
                                    heta_Del_new.append(heta_Del_block[a_heta])
                            for a3 in heta_Del_new:
                                a4=a3[0]
                                tra1[SV_ID]['b'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
                                a4=a3[-1]
                                tra1[SV_ID]['b'][-1].append(int(del_hash[a4][2])-1)
                    else:
                        if not 'b' in tra1[SV_ID].keys():
                            tra1[SV_ID]['b']=[]
                        t1=[]
                        for a3 in k2b:
                            if not a3=='^':
                                t1.append(a3)
                            else:
                                t1[-1]+=a3
                        t2=[t1[0]]
                        for a3 in t1[1:]:
                            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
                                t2[-1]+=a3
                            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
                                t2[-1]+=a3
                            else:
                                t2.append(a3)
                        a3='left'
                        a4=t2[0]
                        l_chr=bp_hash[a3][0]
                        r_chr=bp_hash[a4[0]][0]
                        if not '^' in a4:
                            if not a4[0]==k1[0]:
                                 tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
                                 tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
                        elif '^' in a4:
                            tra1[SV_ID]['b'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
                            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                        for t3 in range(len(t2)-1):
                            a3=t2[t3]
                            a4=t2[t3+1]
                            l_chr=bp_hash[a3[0]][0]
                            r_chr=bp_hash[a4[0]][0]
                            if not '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
                            elif '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
                            elif not '^' in a3 and '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                            elif '^' in a3 and '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
                        if len(t2)>1:
                            a3=t2[t3+1]
                        else:
                            a3=t2[0]                                
                        a4='right'
                        l_chr=bp_hash[a3[0]][0]
                        r_chr=bp_hash[a4][0]
                        if not '^' in a3:
                            if not a3[-1]==k1[-1]:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
                        if '^' in a3:
                            tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
                            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])
        def tra_info_add(k1,k2):
            for k3 in sv_info[k1][k2]:
                SV_ID='_'.join([str(i) for i in k3]+['S'])
                if not SV_ID in tra1.keys():
                    tra1[SV_ID]={}
                k2a=k2.split('/')[0]
                k2b=k2.split('/')[1]
                bp_hash={}
                block_rec=0
                block_hash=[]
                for a3 in k3[:-1]:
                    if a3 in chromos or not a3.isdigit():
                        block_hash.append([a3])
                    else:
                        block_hash[-1].append(a3)
                for a3 in block_hash:
                    for a4 in range(len(a3)-2):
                        bp_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                        block_rec+=1
                for a3 in bp_hash.keys():
                    temp=[]
                    for a4 in bp_hash[a3][1:]:
                        temp.append(int(a4)-1)
                        temp.append(int(a4))
                    bp_hash[a3][1:]=temp
                bp_hash['left']=[bp_hash[k1[0]][0],bp_hash[k1[0]][1],bp_hash[k1[0]][2]]
                bp_hash['right']=[bp_hash[k1[-1]][0],bp_hash[k1[-1]][3],bp_hash[k1[-1]][4]]
                ref_allele={}
                for a3 in bp_hash.keys():
                    ref_allele[a3]=[bp_hash[a3][0]]
                    for a4 in bp_hash[a3][1:]:
                        ref_allele[a3].append(ref_base_readin(ref,bp_hash[a3][0],a4))
                if not k2a==k1.split('/')[0] and del_flag_SA(k1.split('/')[0],k2a)==0:
                    flag1=0#flag1==0:w/o inversion in the alt structure
                    if '^' in k2a:
                        flag1+=1
                    flag2=0#flag2==0:w/o duplication in the alt structure
                    for j in k2a:
                        if k2a.count(j)>1:
                            flag2+=1
                    flag3=0 #flag3==0: w/o translocation
                    if len(k2a)>1:
                        for i in range(len(k2a)-1):
                            if not ord(k2a[i+1])>ord(k2a[i]):
                                flag3+=1
                    if flag1+flag2+flag3==0:
                        heta_Del_block=[]
                        for a1 in k1.split('/')[0]:
                            if not a1 in k2a:
                                heta_Del_block.append(a1)  
                        if not 'a' in tra1[SV_ID].keys():      
                            tra1[SV_ID]['a']=[]
                        block_hash=[]
                        del_hash={}
                        block_rec=0
                        for a3 in a2[0]:
                            if a3 in chromos:
                                block_hash.append([a3])
                            else:
                                block_hash[-1].append(a3)
                        for a3 in block_hash:
                            for a4 in range(len(a3)-2):
                                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                                block_rec+=1                                
                        if not heta_Del_block==[]:
                            a_heta=0
                            heta_Del_new=[heta_Del_block[0]]
                            while True:
                                a_heta+=1
                                if a_heta==len(heta_Del_block):break
                                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
                                    heta_Del_new[-1]+=heta_Del_block[a_heta]
                                else:
                                    heta_Del_new.append(heta_Del_block[a_heta])
                            for a3 in heta_Del_new:
                                a4=a3[0]
                                tra1[SV_ID]['a'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
                                a4=a3[-1]
                                tra1[SV_ID]['a'][-1].append(int(del_hash[a4][2])-1)
                    else:
                        if not 'a' in tra1[SV_ID].keys():
                            tra1[SV_ID]['a']=[]
                        t1=[]
                        for a3 in k2a:
                            if not a3=='^':
                                t1.append(a3)
                            else:
                                t1[-1]+=a3
                        t2=[t1[0]]
                        for a3 in t1[1:]:
                            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
                                t2[-1]+=a3
                            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
                                t2[-1]+=a3
                            else:
                                t2.append(a3)
                        a3='left'
                        a4=t2[0]
                        l_chr=bp_hash[a3][0]
                        r_chr=bp_hash[a4[0]][0]
                        if not '^' in a4:
                            if not a4[0]==k1[0]:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
                        elif '^' in a4:
                            tra1[SV_ID]['a'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
                            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                        for t3 in range(len(t2)-1):
                            a3=t2[t3]
                            a4=t2[t3+1]
                            l_chr=bp_hash[a3[0]][0]
                            r_chr=bp_hash[a4[0]][0]
                            if not '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
                            elif '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
                            elif not '^' in a3 and '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                            elif '^' in a3 and '^' in a4:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
                        if len(t2)>1:
                            a3=t2[t3+1]
                        else:
                            a3=t2[0]
                        a4='right'
                        l_chr=bp_hash[a3[0]][0]
                        r_chr=bp_hash[a4][0]
                        if not '^' in a3:
                            if not a3[-1]==k1[-1]:
                                tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
                                tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
                        if '^' in a3:
                            tra1[SV_ID]['a'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
                            tra1[SV_ID]['a'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])
                if not k2b==k1.split('/')[1] and del_flag_SA(k1.split('/')[1],k2b)==0:
                    flag1=0#flag1==0:w/o inversion in the alt structure
                    if '^' in k2b:
                        flag1+=1
                    flag2=0#flag2==0:w/o duplication in the alt structure
                    for j in k2b:
                        if k2b.count(j)>1:
                            flag2+=1
                    flag3=0 #flag3==0: w/o translocation
                    if len(k2b)>1:
                        for i in range(len(k2b)-1):
                            if not ord(k2b[i+1])>ord(k2b[i]):
                                flag3+=1
                    if flag1+flag2+flag3==0:
                        heta_Del_block=[]
                        for a1 in k1.split('/')[1]:
                            if not a1 in k2b:
                                heta_Del_block.append(a1)   
                        if not 'b' in tra1[SV_ID].keys():     
                            tra1[SV_ID]['b']=[]
                        block_hash=[]
                        del_hash={}
                        block_rec=0
                        for a3 in a2[0]:
                            if a3 in chromos:
                                block_hash.append([a3])
                            else:
                                block_hash[-1].append(a3)
                        for a3 in block_hash:
                            for a4 in range(len(a3)-2):
                                del_hash[chr(97+block_rec)]=[a3[0],a3[a4+1],a3[a4+2]]
                                block_rec+=1                                
                        if not heta_Del_block==[]:
                            a_heta=0
                            heta_Del_new=[heta_Del_block[0]]
                            while True:
                                a_heta+=1
                                if a_heta==len(heta_Del_block):break
                                if ord(heta_Del_block[a_heta])-ord(heta_Del_block[a_heta-1])==1 and del_hash[heta_Del_block[a_heta]][0]==del_hash[heta_Del_block[a_heta-1]][0]:
                                    heta_Del_new[-1]+=heta_Del_block[a_heta]
                                else:
                                    heta_Del_new.append(heta_Del_block[a_heta])
                            for a3 in heta_Del_new:
                                a4=a3[0]
                                tra1[SV_ID]['b'].append(['DEL',del_hash[a4][0],int(del_hash[a4][1]),ref_allele[a4][2]])
                                a4=a3[-1]
                                tra1[SV_ID]['b'][-1].append(int(del_hash[a4][2])-1)
                    else:
                        if not 'b' in tra1[SV_ID].keys():
                            tra1[SV_ID]['b']=[]
                        t1=[]
                        for a3 in k2b:
                            if not a3=='^':
                                t1.append(a3)
                            else:
                                t1[-1]+=a3
                        t2=[t1[0]]
                        for a3 in t1[1:]:
                            if not '^' in a3 and not '^' in t2[-1] and ord(a3)-ord(t2[-1][-1])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-1]][0]:
                                t2[-1]+=a3
                            elif '^' in a3 and '^' in t2[-1] and ord(t2[-1][-2])-ord(a3[0])==1 and bp_hash[a3[0]][0]==bp_hash[t2[-1][-2]][0]:
                                t2[-1]+=a3
                            else:
                                t2.append(a3)
                        a3='left'
                        a4=t2[0]
                        l_chr=bp_hash[a3][0]
                        r_chr=bp_hash[a4[0]][0]
                        if not '^' in a4:
                            if not a4[0]==k1[0]:
                                 tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3][1])+']'+ref_allele[a4[0]][2]])
                                 tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+'['+r_chr+':'+str(bp_hash[a4[0]][2])+'['])
                        elif '^' in a4:
                            tra1[SV_ID]['b'].append([r_chr, bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3][1])+']'])
                            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3][1],ref_allele[a3][1],ref_allele[a3][1]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                        for t3 in range(len(t2)-1):
                            a3=t2[t3]
                            a4=t2[t3+1]
                            l_chr=bp_hash[a3[0]][0]
                            r_chr=bp_hash[a4[0]][0]
                            if not '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['])
                            elif '^' in a3 and not '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][2],ref_allele[a4[0]][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4[0]][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4[0]][0]+':'+str(bp_hash[a4[0]][2])+'['+ref_allele[a3[-2]][2]])
                            elif not '^' in a3 and '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'])
                            elif '^' in a3 and '^' in a4:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4[0]][3],ref_allele[a4[0]][3],ref_allele[a4[0]][3]+'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2], ']'+r_chr+':'+str(bp_hash[a4[0]][3])+']'+ref_allele[a3[-2]][2]])
                        if len(t2)>1:
                            a3=t2[t3+1]
                        else:
                            a3=t2[0]                                
                        a4='right'
                        l_chr=bp_hash[a3[0]][0]
                        r_chr=bp_hash[a4][0]
                        if not '^' in a3:
                            if not a3[-1]==k1[-1]:
                                tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],']'+l_chr+':'+str(bp_hash[a3[-1]][3])+']'+ref_allele[a4][2]])
                                tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-1]][3],ref_allele[a3[-1]][3],ref_allele[a3[-1]][3]+'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['])
                        if '^' in a3:
                            tra1[SV_ID]['b'].append([r_chr,bp_hash[a4][2],ref_allele[a4][2],'['+l_chr+':'+str(bp_hash[a3[-2]][2])+'['+ref_allele[a4][2]])
                            tra1[SV_ID]['b'].append([l_chr,bp_hash[a3[-2]][2],ref_allele[a3[-2]][2],'['+bp_hash[a4][0]+':'+str(bp_hash[a4][2])+'['+ref_allele[a3[-2]][2]])
        def let_reclust(vec_in):
            if vec_in==[]:
                return []
            else:
                k2e=[]
                k2e=[vec_in[0]]
                for k3 in range(len(vec_in)-1):
                    if '^' in vec_in[k3+1]:
                        if '^' in vec_in[k3] and ord(vec_in[k3][0])-ord(vec_in[k3+1][0])==1:
                            k2e[-1]+=vec_in[k3+1]
                        else:
                            k2e.append(vec_in[k3+1])
                    else:
                        if ord(vec_in[k3+1][0])-ord(vec_in[k3][0])==1 and not '^' in vec_in[k3]:
                            k2e[-1]+=vec_in[k3+1]
                        else:
                            k2e.append(vec_in[k3+1])
                k2f=[]
                for k3 in k2e:
                    if '^' in k3:
                        k5=''
                        for k4 in range(len(k3)/2):
                            k5+=k3[2*k4]
                        k6=k5[::-1]+'^'
                        if not k6 in k2f:
                            k2f.append(k6)
                    else:
                        k2f.append(k3)
                return k2f 
        def dup_let_recombind(vec_in):
            if vec_in==[]:
                return []
            else:
                vec2=sorted(vec_in)
                vec=[[vec2[0]]]
                for ka in vec2[1:]:
                    if ord(ka)-ord(vec[-1][-1])==1:
                        vec[-1].append(ka)
                    else:
                        vec.append([ka])
                vec3=[]
                for ka in vec:
                    if len(ka)==1:
                        vec3.append(ka)
                    else:
                        for kb in range(2,len(ka)+1):
                            for kc in ka[:(1-kb)]:
                                vec3.append([])
                                for kd in range(kb):
                                    vec3[-1].append(ka[ka.index(kc)+kd])                
                vec4=[''.join(i) for i in vec3]
                return vec4
        def comp_info_reorganize(k1,k2):
            del_let=[[],[]]
            dup_let=[[],[]]
            inv_let=[[],[]]
            tra_let=[[],[]]
            k2a=k2.split('/')[0]
            k2b=k2.split('/')[1]
            k2c=[]
            k2d=[]
            for k3 in k2a:
                if not k3=='^':
                    k2c.append(k3)
                else:
                    k2c[-1]+=k3
            for k3 in k2b:
                if not k3=='^':
                    k2d.append(k3)
                else:
                    k2d[-1]+=k3
            for k3 in k1.split('/')[0]:
                if k2a.count(k3)==0:
                    del_let[0].append(k3)
                if k2b.count(k3)==0:
                    del_let[1].append(k3)
                if k2a.count(k3)>1:
                    dup_let[0].append(k3)
                if k2b.count(k3)>1:
                    dup_let[1].append(k3)
            k2e=let_reclust(k2c)
            k2f=let_reclust(k2d)
            k2g=dup_let_recombind(dup_let[0])
            k2h=dup_let_recombind(dup_let[1])
            k2i=[]
            k2j=[]
            for k3 in k2g:
                flag1=0
                for k4 in k2e:
                    if k3 in k4:
                        flag1+=1
                if flag1>1:
                    k2i.append(k3)
            for k3 in dup_let[0]:
                if k2e.count(k3[0])+k2e.count(k3[0]+'^')>0:
                    if not k3[0] in k2i:
                        k2i.append(k3[0])
            for k3 in k2h:
                flag1=0
                for k4 in k2e:
                    if k3 in k4:
                        flag1+=1
                if flag1>1:
                    k2j.append(k3)
            for k3 in dup_let[1]:
                if k2e.count(k3[0])+k2e.count(k3[0]+'^')>0:
                    if not k3[0] in k2j:
                        k2j.append(k3[0])
            k2m=[]
            for k3 in k2e:
                if k3[-1]=='^':
                    k2m.append(k3)
            k2n=[]
            for k3 in k2f:
                if k3[-1]=='^':
                    k2n.append(k3)
            for k3 in sv_info[k1][k2]:
                del_info_add(k3,del_let)
                dup_info_add(k3,[k2i,k2j])
                inv_info_add(k3,[k2m,k2n])
        def hash_collaps():
            for k1 in sv_out.keys():
                for k2 in sv_out[k1].keys():
                    if len(sv_out[k1][k2])>1:
                        temp=[]
                        temp2=[]
                        for k3 in sv_out[k1][k2]:
                            if not k3[:-1] in temp:
                                temp.append(k3[:-1])
                                temp2.append([k3[-1]])
                            else:
                                temp2[temp.index(k3[:-1])].append(k3[-1])
                        for k3 in range(len(temp2)):
                            if len(temp2[k3])>1:
                                if sorted([temp2[k3][0].split(':')[0],temp2[k3][1].split(':')[0]])==['0|1', '1|0']:
                                    if not ':' in temp2[k3][0]:
                                        temp2[k3]=['1|1']
                                    else:   
                                        temp2[k3]=['1|1:'+str(int(temp2[k3][0].split(':')[1])+int(temp2[k3][1].split(':')[1]))]
                        temp3=[]
                        for k3 in range(len(temp2)):
                            temp3.append(temp[k3]+temp2[k3])
                        sv_out[k1][k2]=temp3
        def end_cordi_calcu(pin):
            chromo=pin[0]
            start=int(pin[1])
            end=0
            for x in pin[7].split(';'):
                if 'END' in x.split('='):
                    end=int(x.split('=')[1])
            if not end==0:
                return [chromo,start,end]
            else:
                return 'Error!'
        def hash_collaps2():
            temp={}
            for k1 in sv_out.keys():
                temp[k1]={}
                for k2 in sv_out[k1].keys():
                    for k3 in sv_out[k1][k2]:
                        pos=end_cordi_calcu(k3)
                        if not pos[1] in temp[k1].keys():
                            temp[k1][pos[1]]={}
                        if not pos[2] in temp[k1][pos[1]].keys():
                            temp[k1][pos[1]][pos[2]]=[]
                        temp[k1][pos[1]][pos[2]].append([k1,k2,k3])
            out={}
            for k1 in temp.keys():
                out[k1]={}
                for k2 in temp[k1].keys():
                    if len(temp[k1][k2])>1:
                        flag=1
                        for k3 in temp[k1][k2].keys():
                            for k4 in temp[k1][k2][k3]:
                                    if not k4[2][4]=='<DUP>':
                                        flag=0
                        if flag==1:
                            for k4 in temp[k1][k2].keys():
                                if not k4==max(temp[k1][k2].keys()):
                                    for k5 in temp[k1][k2][k4]:
                                        del sv_out[k1][k5[2][2]][sv_out[k1][k5[2][2]].index(k5[2])]
                                        if sv_out[k1][k5[2][2]]==[]:
                                            del sv_out[k1][k5[2][2]]
        def hash_collaps3():
            for k1 in sv_out.keys():
                for k2 in sv_out[k1].keys():
                    if len(sv_out[k1][k2])>1:
                        temp1=[]
                        temp2=[]
                        for k3 in range(len(sv_out[k1][k2])):
                            if not sv_out[k1][k2][k3][:5]+sv_out[k1][k2][k3][6:-1] in temp1:
                                temp1.append(sv_out[k1][k2][k3][:5]+sv_out[k1][k2][k3][6:-1])
                                temp2.append(sv_out[k1][k2][k3])
                            else:
                                continue
                        sv_out[k1][k2]=temp2
        def hash_reorder():
            for ka1 in del1.keys():
                if not ka1 in sv_out.keys():
                    sv_out[ka1]={}
                for ka2 in del1[ka1]:
                    REF_AL='N'
                    Pass_Sign='PASS'
                    if ka2[3]<score_Cff:
                        Pass_Sign='LowQual'
                    if ka2[2]=='heta':
                        GenoType='1|0'
                    elif ka2[2]=='hetb':
                        GenoType='0|1'
                    elif ka2[2]=='homo':
                        GenoType='1|1'
                    ka_new=[ka1,ka2[0],ka2[-1],REF_AL,'<DEL>',ka2[3],Pass_Sign,'SVTYPE=DEL;END='+str(ka2[1]),'GT',GenoType]
                    if not ka2[-1] in sv_out[ka1].keys():
                        sv_out[ka1][ka2[-1]]=[]
                    if not ka_new in sv_out[ka1][ka2[-1]]:
                        sv_out[ka1][ka2[-1]].append(ka_new)
            for ka1 in inv1.keys():
                if not ka1 in sv_out.keys():
                    sv_out[ka1]={}
                for ka2 in inv1[ka1]:
                    REF_AL='N'
                    Pass_Sign='PASS'
                    if ka2[3]<score_Cff:
                        Pass_Sign='LowQual'
                    if ka2[2]=='heta':
                        GenoType='1|0'
                    elif ka2[2]=='hetb':
                        GenoType='0|1'
                    elif ka2[2]=='homo':
                        GenoType='1|1'
                    ka_new=[ka1,ka2[0],ka2[-1],REF_AL,'<INV>',ka2[3],Pass_Sign,'SVTYPE=INV;END='+str(ka2[1]),'GT',GenoType]
                    if not ka2[-1] in sv_out[ka1].keys():
                        sv_out[ka1][ka2[-1]]=[]
                    if not ka_new in sv_out[ka1][ka2[-1]]:
                        sv_out[ka1][ka2[-1]].append(ka_new)
            for ka1 in dup1.keys():
                if not ka1 in sv_out.keys():
                    sv_out[ka1]={}
                for ka2 in dup1[ka1]:
                    REF_AL='N'
                    CopyNumber=str(ka2[-1])
                    Pass_Sign='PASS'
                    if ka2[3]<score_Cff:
                        Pass_Sign='LowQual'
                    if ka2[2]=='heta':
                        GenoType='1|0'
                    elif ka2[2]=='hetb':
                        GenoType='0|1'
                    elif ka2[2]=='homo':
                        GenoType='1|1'
                    ka_new=[ka1,ka2[0],ka2[-2],REF_AL,'<DUP:TANDEM>',ka2[3],Pass_Sign,'SVTYPE=DUP;END='+str(ka2[1]),'GT:CN',GenoType+':'+CopyNumber]
                    if not ka2[-2] in sv_out[ka1].keys():
                        sv_out[ka1][ka2[-2]]=[]
                    if not ka_new in sv_out[ka1][ka2[-2]]:
                        sv_out[ka1][ka2[-2]].append(ka_new)
            for ka1 in disperse_dup.keys():
                if not ka1 in sv_out.keys():
                    sv_out[ka1]={}
                for ka2 in disperse_dup[ka1]:
                    REF_AL='N'
                    CopyNumber=str(ka2[-1])
                    Pass_Sign='PASS'
                    if ka2[3]<score_Cff:
                        Pass_Sign='LowQual'
                    if ka2[2]=='heta':
                        GenoType='1|0'
                    elif ka2[2]=='hetb':
                        GenoType='0|1'
                    elif ka2[2]=='homo':
                        GenoType='1|1'
                    ka_new=[ka1,ka2[0],ka2[-2],REF_AL,'<DUP>',ka2[3],Pass_Sign,'SVTYPE=DUP;END='+str(ka2[1]),'GT:CN',GenoType+':'+CopyNumber]
                    if not ka2[-2] in sv_out[ka1].keys():
                        sv_out[ka1][ka2[-2]]=[]
                    if not ka_new in sv_out[ka1][ka2[-2]]:
                        sv_out[ka1][ka2[-2]].append(ka_new)
            for ka1 in tra1.keys():
                ks1=ka1.split('_')[0]
                ks2='_'.join(ka1.split('_')[:-2]+[ka1.split('_')[-1]])
                SV_Score=float(ka1.split('_')[-2])
                Pass_Sign='PASS'
                if SV_Score<score_Cff:
                    Pass_Sign='LowQual'
                if not ks1 in sv_out.keys():
                    sv_out[ks1]={}
                if not ks2 in sv_out[ks1].keys():
                    sv_out[ks1][ks2]=[]
                for ka2 in tra1[ka1].keys():
                    hetx='het'+ka2
                    if ka2=='a':
                        GenoType='1|0'
                    elif ka2=='b':
                        GenoType='0|1'
                    for ka3 in tra1[ka1][ka2]:
                        ka_new=ka3[:2]+[ks2,ka3[2]]+ka3[3:]+[SV_Score,Pass_Sign,'SVTYPE=TRA','GT',GenoType]
                        if not ka_new in sv_out[ks1][ks2]:
                            sv_out[ks1][ks2].append(ka_new)
        def write_VCF_header(output_file):
            fo=open(output_file,'w')
            print>>fo, '##fileformat=VCFv4.1'
            print>>fo,'##fileDate='+time.strftime("%Y%m%d")
            print>>fo,'##reference=hg19'
            print>>fo,'##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">'
            print>>fo,'##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">'
            print>>fo,'##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">'
            print>>fo,'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'
            print>>fo,'##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">'
            print>>fo,'##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">'
            print>>fo,'##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">'
            print>>fo,'##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">'
            print>>fo,'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
            print>>fo,'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
            print>>fo,'##FILTER=<ID=LowQual,Description="Score of final structural - Theoretical Score <-50">'
            print>>fo,'##ALT=<ID=DEL,Description="Deletion">'
            print>>fo,'##ALT=<ID=DUP,Description="Duplication">'
            print>>fo,'##ALT=<ID=INV,Description="Inversion">'
            print>>fo,'##ALT=<ID=TRA,Description="Translocation">'
            print>>fo,'##ALT=<ID=INS,Description="Insertion">'
            print>>fo,'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
            print>>fo,'##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">'
            print>>fo,'##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">'
            print>>fo,'##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">'
            print>>fo,'\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',output_file.split('/')[-1].replace('.vcf','')])
            fo.close()
        def ref_base_readin(ref,chromo,pos):
            fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,chromo,str(pos),str(pos)))
            tre=fref.readline().strip().split()
            REF_AL=fref.readline().strip().split()
            if not REF_AL==[]:
                return REF_AL[0]
            else:
                return 'N'
        def write_VCF_main(output_file):
            fo=open(output_file,'a')
            sv_reorganize={}
            for k1 in sv_out.keys():
                sv_reorganize[k1]={}
                for k2 in sv_out[k1].keys():
                    start=int(k2.split('_')[1])
                    if not start in sv_reorganize[k1].keys():
                        sv_reorganize[k1][start]={}
                    SVtemp_a=[]
                    SVtemp_b=[]
                    for k3 in sv_out[k1][k2]:
                        if not k3[:-1] in SVtemp_a:
                            SVtemp_a.append(k3[:-1])
                            SVtemp_b.append([k3[-1]])
                        else:
                            SVtemp_b[SVtemp_a.index(k3[:-1])].append(k3[-1])
                    SVtemp=[]
                    sv_reorganize[k1][start][k2]=[]
                    for k3 in range(len(SVtemp_a)):
                        if len(SVtemp_b[k3])==2 and SVtemp_b[k3] in [['0|1', '1|0'],['1|0', '0|1']]:
                            SVtemp_b[k3]=['1|1']
                    for k3 in range(len(SVtemp_a)):
                        for k4 in SVtemp_b[k3]:
                            sv_reorganize[k1][start][k2].append(SVtemp_a[k3]+[k4])
            for k1 in chromos:
                if k1 in sv_reorganize.keys():
                    for k2 in sorted(sv_reorganize[k1].keys()):
                        for k3 in sorted(sv_reorganize[k1][k2].keys()):
                            for k4 in sv_reorganize[k1][k2][k3]:
                                if k4[3]=='N':
                                    k4[3]=ref_base_readin(ref,k4[0],k4[1])
                                print >>fo, '\t'.join([str(i) for i in k4])
            fo.close()
        def ROC_produce_files(ref_file,samp_file):
            ref_hash={}
            samp_hash={}
            out={}
            for i in chromos:
                ref_hash[i]=[]
                samp_hash[i]=[]
            fin=open(ref_file)
            for line in fin:
                pin=line.strip().split()
                ref_hash[pin[0]].append([int(i) for i in pin[1:3]])
            fin.close()
            fin=open(samp_file)
            for line in fin:
                pin=line.strip().split()
                samp_hash[pin[0]].append([int(i) for i in pin[1:3]])
            fin.close()
            for k1 in chromos:
                flag1=0
                if not ref_hash[k1]==[]:
                    out[k1]=[]
                    for k2 in ref_hash[k1]:
                        flag2=0
                        for k3 in samp_hash[k1]:
                            if k3[1]<k2[0]: continue
                            elif k3[0]>k2[1]: continue
                            else:
                                if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max(k2[1]-k2[0],k3[1]-k3[0]))>0.5:
                                    flag2+=1
                        if flag2>0:
                            flag1+=1
                    out[k1]=[flag1,len(ref_hash[k1]),len(samp_hash[k1]),float(flag1)/float(len(ref_hash[k1]))]
            return out
        def MissedSV_Produce_files(ref_file,samp_file):
            ref_hash={}
            samp_hash={}
            out={}
            for i in chromos:
                ref_hash[i]=[]
                samp_hash[i]=[]
            fin=open(ref_file)
            for line in fin:
                pin=line.strip().split()
                ref_hash[pin[0]].append([int(i) for i in pin[1:3]])
            fin.close()
            fin=open(samp_file)
            for line in fin:
                pin=line.strip().split()
                samp_hash[pin[0]].append([int(i) for i in pin[1:3]])
            fin.close()
            for k1 in chromos:
                flag1=0
                if not ref_hash[k1]==[]:
                    out[k1]=[]
                    for k2 in ref_hash[k1]:
                        flag2=0
                        for k3 in samp_hash[k1]:
                            if k3[1]<k2[0]: continue
                            elif k3[0]>k2[1]: continue
                            else:
                                if float(sorted(k2+k3)[2]-sorted(k2+k3)[1])/float(max(k2[1]-k2[0],k3[1]-k3[0]))>0.5:
                                    flag2+=1
                        if flag2>0:
                            flag1+=1
                        else:
                            out[k1].append(k2)
            return out
        def bed_writing(filename,hash):
            fo=open(filename,'w')
            for k3 in chromos:
                if k3 in hash.keys():
                    for k4 in hash[k3]:
                        print >>fo, ' '.join([str(i) for i in k4])
            fo.close()
        def ROC_writing(filename,hash):
            fo=open(filename,'w')
            for k1 in hash.keys():
                for k2 in hash[k1].keys():
                    for k3 in chromos:
                        if k3 in hash[k1][k2].keys():
                            print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+hash[k1][k2][k3]])
            fo.close()
        def MissSV_writing(filename,hash):
            fo=open(filename,'w')
            for k1 in hash.keys():
                for k2 in hash[k1].keys():
                    for k3 in chromos:
                        if k3 in hash[k1][k2].keys():
                            for k4 in hash[k1][k2][k3]:
                                print >>fo, ' '.join([str(i) for i in [k3]+k4+[k1,k2]])
            fo.close()
        def MissSV_Compare(File1,File2):
            hash1={}
            hash2={}
            for k1 in chromos:
                hash1[k1]={}
                hash2[k1]={}
            fin=open(File1)
            for line in fin:
                pin=line.strip().split()
                if not pin[3] in hash1[pin[0]].keys():
                    hash1[pin[0]][pin[3]]={}
                if not pin[4].upper() in hash1[pin[0]][pin[3]].keys():
                    hash1[pin[0]][pin[3]][pin[4].upper()]=[]
                hash1[pin[0]][pin[3]][pin[4].upper()].append([pin[1],pin[2]])
            fin.close()
            fin=open(File2)
            for line in fin:
                pin=line.strip().split()
                if not pin[3] in hash2[pin[0]].keys():
                    hash2[pin[0]][pin[3]]={}
                if not pin[4].upper() in hash2[pin[0]][pin[3]].keys():
                    hash2[pin[0]][pin[3]][pin[4].upper()]=[]
                hash2[pin[0]][pin[3]][pin[4].upper()].append([pin[1],pin[2]])
            fin.close()
            hash3={}
            for k1 in hash1.keys():
                hash3[k1]={}
                for k2 in hash1[k1].keys():
                    hash3[k1][k2]={}
                    if k2 in hash2[k1].keys():
                        for k3 in hash1[k1][k2].keys():
                            hash3[k1][k2][k3]=[]
                            if k3 in hash2[k1][k2].keys():
                                for k4 in hash1[k1][k2][k3]:
                                    if not k4 in hash2[k1][k2][k3]:
                                        hash3[k1][k2][k3].append(k4)
                            else:
                                hash3[k1][k2][k3]=hash1[k1][k2][k3]
                    else:
                        hash3[k1][k2]=hash1[k1][k2]
            fo=open(File1+'.vs.'+File2.split('/')[-1],'w')
            for k1 in chromos:
                if k1 in hash3.keys():
                    for k2 in hash3[k1].keys():
                        for k3 in hash3[k1][k2].keys():
                            for k4 in hash3[k1][k2][k3]:
                                print >>fo, ' '.join([str(i) for i in [k1]+k4+[k2,k3]])
            fo.close()
        def dup_collaps(dup1):
            out={}
            for k1 in dup1.keys():
                out[k1]={}
                for k2 in dup1[k1]:
                    if not k2[0] in out[k1].keys():
                        out[k1][k2[0]]={}
                    if not k2[1] in out[k1][k2[0]].keys():
                        out[k1][k2[0]][k2[1]]=[]
                    out[k1][k2[0]][k2[1]].append(k2)
            out=dup_redun_refine(out)
            out=dup_homo_refine(out)
            out=dup_name_refine(out)
            out2=dup_collaps_refine2(out)
            out_final={}
            for k1 in out2.keys():
                out_final[k1]=[]
                for k2 in sorted(out2[k1].keys()):
                    for k3 in sorted(out2[k1][k2].keys()):
                        for k4 in out2[k1][k2][k3]:
                            out_final[k1].append(k4)
            return out_final
        def dup_redun_refine(out):
            for k1 in out.keys():
                for k2 in sorted(out[k1].keys()):
                    for k3 in sorted(out[k1][k2].keys()):
                        if len(out[k1][k2][k3])>1:
                            temp_dup=[]
                            temp_dup2=[]
                            for k4 in out[k1][k2][k3]:
                                if not k4[:3]+k4[4:-1] in temp_dup:
                                    temp_dup.append(k4[:3]+k4[4:-1])
                                    temp_dup2.append(k4)
                            out[k1][k2][k3]=temp_dup2
            return out
        def dup_homo_refine(out):
            for k1 in out.keys():
                for k2 in sorted(out[k1].keys()):
                    for k3 in sorted(out[k1][k2].keys()):
                        if len(out[k1][k2][k3])>1:
                            temp_dup=[]
                            temp_dup2=[]
                            for k4 in out[k1][k2][k3]:
                                if not k4[:2]+k4[4:] in temp_dup:
                                    temp_dup.append(k4[:2]+k4[4:])
                                    temp_dup2.append(k4)
                                else:
                                    temp_dup2[temp_dup.index(k4[:2]+k4[4:])][2]='homo'
                            out[k1][k2][k3]=temp_dup2
            return out
        def dup_name_refine(out):
            for k1 in out.keys():
                for k2 in sorted(out[k1].keys()):
                    for k3 in sorted(out[k1][k2].keys()):
                        if len(out[k1][k2][k3])>1:
                            out[k1][k2][k3]=[out[k1][k2][k3][0]]
            return out
        def dup_collaps_refine(out):
            outout={}
            for k1 in out.keys():
                outout[k1]={}
                out2=[]
                for k2 in sorted(out[k1].keys()):
                    for k3 in sorted(out[k1][k2].keys()):
                        out2.append([k2,k3])
                out3=[out2[0]]
                for k2 in out2[1:]:
                    add_flag=0
                    for k3 in range(len(out3)):
                        if k2[0]==out3[k3][-1]:
                            out3[k3]+=k2
                            add_flag=1
                    if add_flag==0:
                        out3.append(k2)
                for k2 in out3:
                    if len(k2)>2:
                        rec_temp=[]
                        for k3 in range(len(k2)/2):
                            rec_temp+=out[k1][k2[k3*2]][k2[k3*2+1]]
                        rec_t2=[rec_temp[0]]
                        for k3 in range(len(rec_temp)-1):
                            if abs(rec_temp[k3+1][-1]-rec_temp[k3][-1])<2:
                                rec_t2[-1][1]=rec_temp[k3+1][1]
                            else:
                                rec_t2.append(rec_temp[k3+1])
                    else:
                        rec_t2=out[k1][k2[0]][k2[1]]
                    for k3 in rec_t2:
                        if not k3[0] in outout[k1].keys():
                            outout[k1][k3[0]]={}
                        if not k3[1] in outout[k1][k3[0]].keys():
                            outout[k1][k3[0]][k3[1]]=[]
                        outout[k1][k3[0]][k3[1]].append(k3)
            return outout
        def dup_collaps_refine2(out):
            out2=dup_collaps_refine(out)
            for k1 in out2.keys():
                for k2 in out2[k1].keys():
                    for k3 in out2[k1][k2].keys():
                        if len(out2[k1][k2][k3])>1:
                            test=[]
                            for k4 in out2[k1][k2][k3]:
                                test.append(k4[2])
                            if 'homo' in test:
                                out2[k1][k2][k3]=[out2[k1][k2][k3][test.index('homo')]]
                            else:
                                out2[k1][k2][k3]=[out2[k1][k2][k3][0]]
            return out2
        def SV_Info_Write_svelter(sv_info):
            temp1={}
            for k1 in sv_info.keys():
                for k2 in sv_info[k1].keys():
                    for k3 in sv_info[k1][k2]:
                        if not k3[0] in temp1.keys():
                            temp1[k3[0]]={}
                        if not int(k3[1]) in temp1[k3[0]].keys():
                            temp1[k3[0]][int(k3[1])]={}
                        if not int(k3[-2]) in temp1[k3[0]][int(k3[1])].keys():
                            temp1[k3[0]][int(k3[1])][int(k3[-2])]=[]
                        temp1[k3[0]][int(k3[1])][int(k3[-2])].append(k3+[k1,k2])
            fo=open(output_file.replace('.vcf','.svelter'),'w')
            print >>fo, '\t'.join(['chr','start','end','bp_info','ref','alt','score'])
            for k1 in chromos:
                if k1 in temp1.keys():
                    for k2 in sorted(temp1[k1].keys()):
                        for k3 in sorted(temp1[k1][k2].keys()):
                            for k4 in temp1[k1][k2][k3]:
                                chrom_svelter=k1
                                bp_start_svelter=k2
                                bp_end_svelter=k3
                                bps_info_svelter=':'.join(k4[:-3])
                                struc_ref_svelter=k4[-2]
                                struc_alt_svelter=k4[-1]
                                score_svelter=k4[-3]
                                print >>fo, '\t'.join([str(i) for i in [chrom_svelter,bp_start_svelter,bp_end_svelter,bps_info_svelter,struc_ref_svelter,struc_alt_svelter,score_svelter]])
            fo.close()
        def Define_Default_SVIntegrate():
            global score_Cff
            if not '--qc-structure' in dict_opts:
                score_Cff=-20
            else:
                score_Cff=int(dict_opts['--qc-structure'])
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['help=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--input-path, path of .coverage files produced by SVelter.py SVPredict'
            print '--prefix, output prefix for vcf and svelter files (default: input.vcf, input.svelter)'
            print 'Optional Parameters:'
            print '--qc-structure, minimum quality score of a resolved structure to be considered as PASS and included in the output vcf file'
        else:
            Define_Default_SVIntegrate()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                if not '--input-path' in dict_opts.keys():
                    print 'Error: please specify path of input .coverge files using --input-path'
                else:
                    if '--input-path' in dict_opts.keys():
                        if not dict_opts['--input-path'][-1]=='/':
                            dict_opts['--input-path']+='/'
                        InputPath=[dict_opts['--input-path']]
                    else:
                        InputPath=[]
                        for fi1 in os.listdir(workdir+'bp_files'):
                            path1=workdir+'bp_files/'+fi1+'/'
                            if os.path.isdir(path1):
                                for fi2 in os.listdir(path1):
                                    path2=path1+fi2+'/'
                                    InputPath.append(path2)
                    ref_path=workdir+'reference/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed '
                    else:
                        if not '--prefix' in dict_opts.keys():
                            print 'Error: please specify output file using --prefix'
                        else:
                            time1=time.time()
                            output_file=dict_opts['--prefix']+'.vcf'
                            ref=ref_file
                            chromos=[]
                            fin=open(ref_index)
                            for line in fin:
                                pin=line.strip().split()
                                chromos.append(pin[0])
                            fin.close()
                            for path2 in InputPath:
                                sv_info={}
                                for k3 in os.listdir(path2):
                                    if k3.split('.')[-1]=='coverge':
                                        read_in_structures(path2+k3)
                                SV_Info_Write_svelter(sv_info)
                                dup1={}
                                disperse_dup={}
                                inv1={}
                                del1={}
                                tra1={}
                                sv_rec_2(sv_info)
                                dup1=dup_collaps(dup1)
                                sv_out={}
                                hash_reorder()
                                hash_collaps()
                                hash_collaps2()
                                hash_collaps3()
                                write_VCF_header(output_file)
                                write_VCF_main(output_file)
                            time2=time.time()
                            print 'SVIntegrate Complete !'
                            print 'Time Consuming: '+str(time2-time1)
    if not function_name in ['Index','NullModel','BPSearch','BPIntegrate','SVPredict','SVIntegrate']:
        def Define_Default_AllInOne():
            if '--core' in dict_opts.keys():
                global pool
                pool = Pool(processes=int(dict_opts['--core']))
            global model_comp
            if not '--null-model' in dict_opts.keys():
                model_comp='S'
            else:
                if dict_opts['--null-model'] in ['S','Simple']:
                    model_comp='S'
                else:
                    model_comp='C'
            global QCAlign
            if '--qc-align' in dict_opts.keys():
                QCAlign=int(dict_opts['--qc-align'])
            else:
                QCAlign=20
            global QCSplit
            if '--qc-split' in dict_opts.keys():
                QCSplit=int(dict_opts['--qc-split'])
            else:
                QCSplit=20
            global NullSplitLen_perc
            if '--split-min-len' in dict_opts.keys():
                NullSplitLen_perc=int(dict_opts['--split-min-len'])
            else:
                NullSplitLen_perc=0.9
            global KeepFile
            if '--keep-temp-files' in dict_opts.keys():
                KeepFile=dict_opts['--keep-temp-files']
            else:
                KeepFile='No'
            global KeepFigure
            if '--keep-temp-figs' in dict_opts.keys():
                KeepFigure=dict_opts['--keep-temp-figs']
            else:
                KeepFigure='No'
            global Trail_Number
            if '--num-iteration' in dict_opts.keys():
                Trail_Number=int(dict_opts['--num-iteration'])
            else:
                Trail_Number=10000
            global Ploidy
            if '--ploidy' in dict_opts.keys():
                Ploidy=int(dict_opts['--ploidy'])
            else:
                Ploidy=2
            global ILCff_STD_Time
            if '-S' in dict_opts.keys():
                ILCff_STD_Time=int(dict_opts['-S'])
            else:
                ILCff_STD_Time=3
        def Indicator_Readin(temp_inter):
            ft=open(temp_inter)
            F_test='ERROR'
            for line in ft:
                pin=line.strip().split()
                if pin[0]=='SVelter':
                    F_test=pin[1]
            ft.close()
            return F_test
        def path_modify(path):
            if not path[-1]=='/':
                path+='/'
            return path
        def check_scripts(Code_path):
            flag=0
            out=[]
            Code0_file=Code_path+'SVelter0.Ref.Index.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter1.NullModel.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter1.NullModel.Figure.a.r'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter1.NullModel.Figure.d.r'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter1.NullModel.Figure.d2.r'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter2.BP.Searching.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter3.BPIntegrate.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter4.StructureResolvation.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            Code0_file=Code_path+'SVelter5.result.integrate.py'
            if not os.path.isfile(Code0_file):
                flag+=1
                out.append(Code0_file)
            return out
        def SamplingPercentage_read_in():
            if '--null-copyneutral-perc' in dict_opts.keys():
                SamplingPercentage=float(dict_opts['--null-copyneutral-perc'])
            else:
                SamplingPercentage=0.001
            return SamplingPercentage
        def cn2_file_read_in():
            if '--copyneutral' in dict_opts.keys():
                cn2_file=dict_opts['--copyneutral']
            else:
                cn2_file=workdir+'reference/CN2.bed'
            return cn2_file
        def ex_file_read_in():
            if '--exclude' in dict_opts.keys():
                ex_file=dict_opts['--exclude']
            else:
                ex_file=workdir+'reference/Exclude.bed'
            return ex_file
        def cn2_length_read_in():
            if '--null-copyneutral-length' in dict_opts.keys():
                cn2_length=int(dict_opts['--null-copyneutral-length'])
            else:
                cn2_length=2000
            return cn2_length
        def chromos_read_in(ref_index):
            whole_genome={}
            fref=open(ref_index)
            for line in fref:
                pref=line.strip().split()
                whole_genome[pref[0]]=[int(pref[1])]
            fref.close()
            return whole_genome
        def run_SVelter0_chrom(chrom_name):
            os.system(r'''%s --workdir %s --ref %s --ex %s --sample %s --chr %s'''%(Code0_file,workdir,ref_file,ex_file,sin_bam_file,chrom_name))
        def run_SVelter1_chrom(sin_bam_file):
            os.system(r'''%s %s --keep-temp-files %s --keep-temp-figs %s --null-model %s --workdir %s --sample %s'''%(Code_File,Code1_Function,KeepFile,KeepFigure,model_comp,workdir,sin_bam_file)) 
        def run_SVelter1_Single_chrom(sin_bam_file,chromos_single):
            os.system(r'''%s %s --keep-temp-files %s --keep-temp-figs %s --null-model %s --workdir %s --sample %s --chromosome %s'''%(Code_File,Code1_Function,KeepFile,KeepFigure,model_comp,workdir,sin_bam_file,chromos_single)) 
        def run_SVelter2_chrom(chrom_name,sin_bam_file,ILCff_STD_Time):
            os.system(r'''%s %s --chromosome %s --workdir %s --sample %s --null-model %s -S %s'''%(Code_File,Code2_Function,chrom_name,workdir,sin_bam_file,model_comp,ILCff_STD_Time))
            print chrom_name+' done!'
        def run_SVelter3_chrom(sin_bam_file):
            os.system(r'''%s %s --batch %s --workdir %s --sample %s'''%(Code_File,Code3_Function,dict_opts['--batch'],workdir,sin_bam_file)) 
        def run_SVelter4_chrom(txt_name,sin_bam_file):
            os.system(r'''%s %s --workdir %s --bp-file %s --sample %s --num-iteration %s --ploidy %s --null-model %s'''%(Code_File,Code4_Function,workdir,txt_name,sin_bam_file,str(Trail_Number),str(Ploidy),model_comp))
            print txt_name+' done!'
        def run_SVelter5_chrom(path2,out_vcf):
            os.system(r'''%s %s --workdir %s --input-path %s --prefix %s'''%(Code_File,Code5_Function,workdir,path2,out_vcf))
        def Code_Files_Define():
            global Code_File
            global Code0_Function
            global Code1_Function
            global Code2_Function
            global Code3_Function
            global Code4_Function
            global Code5_Function
            global RCode_Path
            global Code1a_file
            global Code1d_file
            global Code1d2_file
            Code_File=script_name
            Code0_Function='Index'
            Code1_Function='NullModel'
            Code2_Function='BPSearch'
            Code3_Function='BPIntegrate'
            Code4_Function='SVPredict'
            Code5_Function='SVIntegrate'
            RCode_Path=workdir+'reference/'
            Code1a_file=RCode_Path+'SVelter1.NullModel.Figure.a.r'
            Code1d_file=RCode_Path+'SVelter1.NullModel.Figure.d.r'
            Code1d2_file=RCode_Path+'SVelter1.NullModel.Figure.d2.r'
        opts,args=getopt.getopt(sys.argv[1:],'o:h:S:',['help=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print 'Usage:'
            print 'SVelter.py [options] [parameters]'
            print 'options:'
            print 'NullModel' 
            print 'BPSearch'  
            print 'BPIntegrate' 
            print 'SVPredict'
            print 'SVIntegrate'
            print 'Required Parameters:'
            print '--workdir, writable working directory.'
            print '--sample, input alignment file in bam format'
            print ' '
            print 'Optional Parameters:'
            print '--prefix, output prefix for vcf and svelter files (default: input.vcf, input.svelter)'
            print '--num-iteration, maximum number of iterations per structure will run'
            print '--ploidy, limit algorithm to specific zygosity (0:heterozygous only; 1:homozygous only; 2:both; default:2)'
            print '--null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
            print '--null-copyneutral-length, minimum length requirement for --copyneutral regions used to build null model (default: 2000)'
            print '--null-copyneutral-perc, percentage of regions from --copyneutral to utilize (default: 0.1)'
            print '--null-random-length, specify the length of random regions if --copyneutral parameter not used (default: 5000)'
            print '--null-random-num, specify the number of random regions if --copyneutral parameter not used (default: 10000)'
            print '--qc-align, minimum alignment quality required for mapped reads in bam file (default: 20)'
            print '--qc-split, minimum alighment of clipped parts of reads considered as a soft clip (default: 20)'
            print '--split-min-len, the minumum length of clip read considered as split; (default:10% of read length)'
            print '--qc-structure, minimum quality score of a resolved structure to be considered as PASS and included in the output vcf file'
            print '--qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary'
            print '--qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/'
            print '--qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported (default: 0.0)'
            print '--batch, specify number of structures in each separate file (if 0, output files will be calssified by chromosomes; default, all BP clustered will be integrated in one txt file)'
        else:
            Define_Default_AllInOne()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                if not os.path.isdir(workdir):
                    print 'Error: working directory does not exit!'
                Code_Files_Define()
                if not '--sample' in dict_opts.keys() and not '--samplePath' in dict_opts.keys():
                    print 'Error: please specify input file using --sample'
                else:
                    if '--sample' in dict_opts.keys():
                        bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                        bam_files=[dict_opts['--sample']]
                    else:
                        bam_path=path_modify(dict_opts['--samplePath'])
                        bam_files=[]
                        for file in os.listdir(bam_path):
                            if file.split('.')[-1]=='bam':
                                bam_files.append(bam_path+file)
                    ref_path=workdir+'reference/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed '
                    else:
                        whole_genome=chromos_read_in(ref_index)
                        len_genome=0
                        for i in whole_genome.keys():
                            len_genome+=whole_genome[i][0]
                        chromos=whole_genome.keys()
                        chr_name_check=0
                        fin=open(ref_index)
                        chr_ref_check=[]
                        for line in fin:
                            pin=line.strip().split()
                            chr_ref_check.append(pin[0])
                        fin.close()
                        for filein_bam in bam_files:
                            chr_bam_check=[]
                            fin=os.popen(r'''samtools view -H %s'''%(filein_bam))
                            for line in fin:
                                pin=line.strip().split()
                                if pin[0]=='@SQ':
                                    chr_bam_check.append(pin[1].split(':')[1])
                            fin.close()
                        if not chr_ref_check==chr_bam_check:
                            print 'Warning: please make sure the reference file matches the bam'
                        chr_flag=0
                        if 'chr' in chr_ref_check[0]:
                            chr_flag=1
                        SamplingPercentage=float(SamplingPercentage_read_in())
                        cn2_file=cn2_file_read_in()
                        ex_file=ex_file_read_in()
                        cn2_length=int(cn2_length_read_in())
                        Gap_Refs=[ex_file]
                        if not os.path.isfile(cn2_file):
                            cn2_path='/'.join(cn2_file.split('/')[:-1])+'/'
                            if not os.path.isdir(cn2_path): 
                                os.system(r'''mkdir %s'''%(cn2_path))
                            if not '--null-random-length' in dict_opts.keys():
                                dict_opts['--null-random-length']=5000
                            else:
                                dict_opts['--null-random-length']=int(dict_opts['--null-random-length'])
                            if not '--null-random-num' in dict_opts.keys():
                                dict_opts['--null-random-num']=10000
                            else:
                                dict_opts['--null-random-num']=int(dict_opts['--null-random-num'])
                            cn2_length=dict_opts['--null-random-length']-100
                            fo=open(cn2_file,'w')
                            for i in sorted(whole_genome.keys()):
                                num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--null-random-num'])
                                reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--null-random-length']) for j in range(num_i)]
                                for j in sorted(reg_i):
                                    print >>fo, ' '.join([i,str(j),str(j+dict_opts['--null-random-length']-1)])
                            fo.close()
                            SamplingPercentage=1
                        if not os.path.isfile(ex_file):
                                fo=open(ex_file,'w')
                                for chr_ex in chromos:
                                    print >>fo, ' '.join([chr_ex,'0','0'])
                                fo.close()
                        if '--prefix' in dict_opts.keys():
                            out_vcf=dict_opts['--prefix']+'.vcf'
                            out_svelter=dict_opts['--prefix']+'.svelter'
                        else:
                            out_vcf=workdir+dict_opts['--sample'].split('/')[-1].replace('.bam','.vcf')
                            out_svelter=workdir+dict_opts['--sample'].split('/')[-1].replace('.bam','.svelter')
                            print 'Warning: output file is not specified'
                            print 'output file: '+out_vcf
                            print 'output file: '+out_svelter
                        temp_inter_replace=0
                        if '--chromosome' in dict_opts.keys():
                            chrom_single=dict_opts['--chromosome']
                            if not chrom_single in chromos:
                                print 'Error: please make sure the chromosome defined by --chr is correct based on the reference genome'
                                chromos=[]
                            else:
                                chromos=[chrom_single]
                        for sin_bam_file in bam_files:
                            running_time=[]
                            print ' '
                            print 'Step1: Running null parameters for '+sin_bam_file.split('/')[-1].replace('.bam','')+' ...'
                            time1=time.time()
                            if len(chromos)>1:
                                run_SVelter1_chrom(sin_bam_file)
                            elif len(chromos)==1:
                                run_SVelter1_Single_chrom(sin_bam_file,chromos[0])
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Null model built for '+sin_bam_file.split('/')[-1].replace('.bam','')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step2: Searching for BreakPoints of sample '+sin_bam_file.split('/')[-1].replace('.bam','')+' ...'
                            time1=time.time()
                            for x in chromos:
                                print x
                                run_SVelter2_chrom(x,sin_bam_file,ILCff_STD_Time)
                                print time.time()-time1
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Break points searching done for sample:'+sin_bam_file.split('/')[-1].replace('.bam','')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step3: Integrating breakpoints ... '
                            if not '--batch' in dict_opts.keys():
                                dict_opts['--batch']='0'
                            time1=time.time()
                            run_SVelter3_chrom(sin_bam_file)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Break points cluster done for sample:'+sin_bam_file.split('/')[-1].replace('.bam','')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step4: Resolving structure ... '
                            time1=time.time()
                            for k1 in os.listdir(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'):
                                if k1==dict_opts['--sample'].split('/')[-1].replace('.bam',''):
                                    path1=workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'+k1+'/'
                                    for k2 in os.listdir(path1):
                                        path2=path1+k2+'/'
                                        all_txt_files=[]
                                        for k3 in os.listdir(path2):
                                            if k3.split('.')[-1]=='txt':
                                                all_txt_files.append(path2+k3)
                                        for x in all_txt_files:
                                            run_SVelter4_chrom(x,sin_bam_file)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Structure resolved !'
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            for k1 in os.listdir(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'):
                                path1=workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'+k1+'/'
                                for k2 in os.listdir(path1):
                                    path2=path1+k2+'/'
                                    print 'Step5: Integrating results in VCF file: '+out_vcf+' ... '
                                    time1=time.time()
                                    run_SVelter5_chrom(path2,out_vcf)
                                    time2=time.time() 
                                    running_time.append(time2-time1)
                                    if temp_inter_replace==0:
                                        print out_vcf+' completed! '
                                        print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print 'Total Running Time:'+' '.join([str(i) for i in running_time])
                        if os.path.isfile(out_vcf):
                            NullPath=workdir+'NullModel.'+dict_opts['--sample'].split('/')[-1]
                            BPPath=workdir+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]
                            TXTPath=workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]
                            RefPath=workdir+'reference/'
                            os.system(r'''rm -r %s'''%(NullPath))
                            os.system(r'''rm -r %s'''%(BPPath))
                            os.system(r'''rm -r %s'''%(TXTPath))
                            os.system(r'''rm -r %s'''%(RefPath))

