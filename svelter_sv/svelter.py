#!/usr/bin/env python

#!Python
#Usage:
#SVelter.py [option] [Parametres]
#option:
#For debug use only
#command='SVelter.py SVPredict --deterministic-flag 1 --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.chr10.bam'
#command='SVelter.py SVPredict --deterministic-flag 1 --workdir /scratch/remills_flux/xuefzhao/NA12878.NGS/hg19 --sample /scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/alignment/NA12878_S1.chr10.bam --bp-file /scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/bp_files.NA12878_S1.chr10.bam/NA12878_S1.chr10.txt'
#command='SVelter_Add_cram.03312016.py PredefinedBP --input-bed /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/het_RD10_INV.INV.bed --workdir /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/predefinedBP_Test/ --sample /mnt/EXT/Mills-scratch2/Xuefang/NA12878.NGS/alignment/NA12878_S1.bam'
#sys.argv=command.split()

from __future__ import absolute_import

import os
import re
import sys
from svelter_sv import readme
script_name=sys.argv[0]
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

if len(sys.argv)<2:
    readme.print_default_parameters()
else:
    from svelter_sv.function import *
    function_name=sys.argv[1]    
    if function_name=='Clean':
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','ref-index=','help=','batch=','prefix=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','segdup=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_clean()
        else:
            workdir=path_modify(dict_opts['--workdir'])
            clean_svelter_set(workdir+'reference_SVelter/')
            print 'SVelter Cleared Up!'
    if function_name=='Setup':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['support=','deterministic-flag=','ref-index=','help=','batch=','prefix=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','segdup=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        Code_path='/'.join(sys.argv[0].split('/')[:-1])+'/'
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_setup()
        else:
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
            def final_regions_decide(Gap_Hash_Ref1,hash_Cor,chrom):
                GapHash=Gap_Hash_Ref1[chrom]
                GapHash2=calculate_interval_region(hash_Cor,Chromo_Length,chrom)
                GapHash+=GapHash2
                temp_hash={}
                for k1 in GapHash:
                    if not k1[0] in temp_hash.keys():
                        temp_hash[k1[0]]={}
                    if not k1[1] in temp_hash[k1[0]].keys():
                        temp_hash[k1[0]][k1[1]]=[]
                temp_list=[]
                for k1 in sorted(temp_hash.keys()):
                    for k2 in sorted(temp_hash[k1].keys()):
                        temp_list.append([k1,k2])
                temp2_list=[]
                temp2_list.append(temp_list[0])
                for k1 in temp_list[1:]:
                    if k1[0]-temp2_list[-1][1]<1000:
                        if k1[1]>temp2_list[-1][1]:
                            temp2_list[-1][1]=k1[1]
                    else:
                        temp2_list.append(k1)
                temp3_list=[]
                for k1 in temp2_list:
                    if not k1 in temp3_list:
                        temp3_list.append(k1)
                return calculate_interval_region(temp3_list,Chromo_Length,chrom)
            def Gap_Hash_Ref_filter(Gap_Hash_Ref1,Chromo_Length):
                for x in Gap_Hash_Ref1.keys():
                    if len(Gap_Hash_Ref1[x])==1:
                        if Gap_Hash_Ref1[x][0][0]==0:
                            if Gap_Hash_Ref1[x][0][1]==Chromo_Length[x]:
                                del Gap_Hash_Ref1[x]
            def Gap_Hash_Ref1_read_in(Gap_Refs):
                Gap_Hash_Ref1={}
                for Gap_Ref1 in Gap_Refs:
                    fgap=open(Gap_Ref1)
                    for line in fgap:
                        pgap=line.strip().split()
                        if pgap[0] in chromos:
                            if not pgap[0] in Gap_Hash_Ref1.keys():
                                Gap_Hash_Ref1[pgap[0]]=[]
                            Gap_Hash_Ref1[pgap[0]].append(pgap[1:4])
                    fgap.close()
                Gap_Hash_Ref2=bed_hash_short(Gap_Hash_Ref1,Chromo_Length)
                for x in chromos:
                    if not x in Gap_Hash_Ref2.keys():
                        Gap_Hash_Ref2[x]=[[0,0]]
                return Gap_Hash_Ref2
            def Global_para_declear_setup():
                global chromos
                global Chromo_Length
                global Gap_Hash
            def write_ExcludeBed(ExcludeBed):
                if not os.path.isfile(ExcludeBed):
                    fo=open(ExcludeBed,'w')
                    for chr_ex in chromos:
                        print >>fo, ' '.join([chr_ex,'0','0'])
                    fo.close()
            if not '--workdir' in dict_opts.keys():
                print 'working directory not specified'
                print 'all temporal files would be writen under current directory'
                workdir='./'
                #print 'Error: please specify working directory using --workdir'
            else:
                workdir = path_modify(dict_opts['--workdir'])
                path_mkdir(workdir)
            ref_file=0
            if not '--reference' in dict_opts.keys():
                print 'Error: please specify refrence genome using --reference'
            else:
                Global_para_declear_setup()
                ref_file=dict_opts['--reference']
                ref_path='/'.join(ref_file.split('/')[:-1])+'/'
                ref_index=ref_file+'.fai'
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                    print 'Please index reference genome using samtools'
                else:
                    #if not '--svelter-path' in dict_opts.keys():
                    #    print 'Error: please specify path of SVelter scripts using --svelter-path'                       
                    #else:
                        time1=time.time()
                        ref_path=workdir+'reference_SVelter/'
                        if not ref_path=='/'.join(ref_file.split('/')[:-1])+'/':
                            ref_path=workdir+'reference_SVelter/'
                            path_mkdir(ref_path)
                            if not ref_file[0]=='/':
                                print 'Error: refrence should be specified using absolute path ! '
                            os.symlink(ref_file,ref_path+'genome.fa')
                            os.symlink(ref_index,ref_path+'genome.fa.fai')
                            if '--ref-index' in dict_opts.keys():
                                if os.path.isdir(dict_opts['--ref-index']):
                                    ref_index_path=path_modify(dict_opts['--ref-index'])
                                    for ref_index_file in os.listdir(ref_index_path):
                                        if ref_index_file.split('.')[-1]=='GC_Content':
                                            if ref_index_path[0]=='/':
                                                os.symlink(ref_index_path+ref_index_file,ref_path+'genome.GC_Content')
                                            else:
                                                os.system(r'''cp %s %s'''%(ref_index_path+ref_index_file,ref_path+'genome.GC_Content'))
                                        if ref_index_file.split('.')[-1]=='bed' and ref_index_file.split('.')[-2]=='Mappable':
                                            if ref_index_path[0]=='/':
                                                os.symlink(ref_index_path+ref_index_file,ref_path+'genome.Mappable.bed')
                                            else:
                                                os.system(r'''cp %s %s'''%(ref_index_path+ref_index_file,ref_path+'genome.Mappable.bed'))
                            if '--support' in dict_opts.keys():
                                support_path=path_modify(dict_opts['--support'])
                                for k1 in os.listdir(support_path):
                                    if 'SVelter' in k1 and k1.split('.')[-1]=='r':
                                        if support_path[0]=='/':
                                            os.symlink(support_path+k1,ref_path+k1)
                                        else:
                                            os.system(r'''cp %s %s'''%(support_path+k1,ref_path))                                        
                                    if 'CN2' in k1:
                                        if not '--copyneutral' in dict_opts.keys():
                                            dict_opts['--copyneutral']=support_path+k1
                                    if 'Exclude' in k1:
                                        if not '--exclude' in dict_opts.keys():
                                            dict_opts['--exclude']=support_path+k1
                                    if 'Segdup' in k1:
                                        if not '--segdup' in dict_opts.keys():
                                            dict_opts['--segdup']=support_path+k1
                            if '--copyneutral' in dict_opts.keys():
                                if dict_opts['--copyneutral'][0]=='/':
                                    os.symlink(dict_opts['--copyneutral'],ref_path+'CN2.bed')
                                else:
                                    os.system(r'''cp %s %s'''%(dict_opts['--copyneutral'],ref_path+'CN2.bed'))
                            else:
                                [whole_genome,len_genome]=calculate_len_genome(ref_file)
                                random_produce_cn2_region(ref_path+'CN2.bed',whole_genome,len_genome,dict_opts)
                            if '--exclude' in dict_opts.keys():
                                if dict_opts['--exclude'][0]=='/':
                                    os.symlink(dict_opts['--exclude'],ref_path+'Exclude.bed')
                                else:
                                    os.system(r'''cp %s %s'''%(dict_opts['--exclude'],ref_path+'Exclude.bed'))
                            else:
                                chromos=chromos_readin_list(ref_file)
                                random_produce_exclude_region(ref_path+'Exclude.bed',chromos)
                            if '--segdup' in dict_opts.keys():
                                if dict_opts['--segdup'][0]=='/':
                                    os.symlink(dict_opts['--segdup'],ref_path+'Segdup.bed')
                                else:
                                    os.system(r'''cp %s %s'''%(dict_opts['--segdup'],ref_path+'Segdup.bed'))
                        ref_file=ref_path+'genome.fa'
                        ref_index=ref_file+'.fai'
                        ExcludeBed=ref_path+'Exclude.bed'
                        [chromos,Chromo_Length]=chromos_info_readin(ref_index)
                        write_ExcludeBed(ExcludeBed)
                        fout_Name='.'.join(ref_file.split('.')[:-1])+'.Mappable.bed'
                        fout_N2='.'.join(ref_file.split('.')[:-1])+'.GC_Content'
                        if not os.path.isfile(fout_Name):
                            Gap_Refs=[ExcludeBed]
                            Gap_Hash_Ref1=Gap_Hash_Ref1_read_in(Gap_Refs)
                            Gap_Hash_Ref_filter(Gap_Hash_Ref1,Chromo_Length)
                            Gap_Hash=Gap_Hash_Initiate(chromos)
                            file_initiate(fout_Name)
                            file_initiate(fout_N2)                        
                            for chrom in chromos:
                                fref=os.popen(r'''samtools faidx %s %s:'''%(ref_file,chrom))
                                pref=fref.readline().strip().split()
                                while True:
                                    pref=fref.readline().strip().split()
                                    if not pref:break
                                    Gap_Hash[chrom].append(pref[0])
                                fref.close()
                                fout=open(fout_Name,'a')
                                fout2=open(fout_N2,'a')
                                hash_key=chrom
                                if not Gap_Hash[hash_key]==[]:
                                    hash_Cor=[]
                                    hash_cal=0
                                    if not ''.join(set(Gap_Hash[hash_key][0])) in ['N','n','Nn','nN']:
                                        hash_Cor.append([0])
                                    for hts in Gap_Hash[hash_key][1:]:
                                        hash_cal+=len(hts)
                                        if ''.join(set(hts)) in ['N','n','Nn','nN']:
                                            if len(hash_Cor)==0: continue
                                            else:   
                                                if len(hash_Cor[-1])==1:
                                                    hash_Cor[-1].append(hash_cal)
                                        else:
                                            if len(hash_Cor)==0:
                                                hash_Cor.append([hash_cal])
                                            elif len(hash_Cor[-1])==2:
                                                hash_Cor.append([hash_cal])
                                    hash_Cor=hash_Cor_modify(hash_Cor,Chromo_Length,hash_key,chrom)
                                    hash_to_Seq=''.join(Gap_Hash[hash_key])
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
                        time2=time.time()
                        print 'Suppport files completely set !'
                        print 'Time Consuming:'+str(time2-time1)
    if function_name=='NullModel':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','out-path=','help=','long-insert=','batch=','prefix=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
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
                KeepFile='No'
            global KeepFigure
            if '--keep-temp-figs' in dict_opts.keys():
                KeepFigure=dict_opts['--keep-temp-figs']
            else:
                KeepFigure='Yes'
        dict_opts_modify(dict_opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_nullmodel()
        else:
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
            def clean_files():
                os.system('''rm  %s'''%(InsertLenNullTemp))
                os.system('''rm  %s'''%(DRNullTemp))
                os.system('''rm  %s'''%(SplitNullTemp))
                os.system('''rm  %s'''%(ILNullTemp))
                os.system('''rm  %s'''%(TBNullTemp))
                os.system('''rm  %s'''%(RDNullTemp))
            def global_para_declaration_nullmodel():
                global bam_path
                global bam_files
                global bam_file_appdix
                global cn2_file
                global len_genome
                global NullPath
                global ref_path
                global ref_file
                global ref_file
                global ref_index
                global SamplingPercentage
                global whole_genome
                bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                bam_files=[dict_opts['--sample']]
                bam_file_appdix=dict_opts['--sample'].split('.')[-1]
                ref_path=workdir+'reference_SVelter/'
                ref_file=ref_path+'genome.fa'
                ref_index=ref_file+'.fai'
                cn2_file=ref_path+'CN2.bed'
            if not '--workdir' in dict_opts.keys():
                print 'working directory not specified'
                print 'all temporal files would be writen under current directory'
                workdir='./'
            else:
                workdir=path_modify(dict_opts['--workdir'])
            if not '--sample' in dict_opts.keys():
                print 'Error: please specify either input file using --sample'
            else:
                global_para_declaration_nullmodel()
                if not os.path.isfile(ref_index):
                    print 'Error: reference genome not indexed'
                else:
                    SamplingPercentage=SamplingPercentage_readin(dict_opts)
                    [whole_genome,len_genome]=calculate_len_genome(ref_file)
                    chromos=chromos_readin_list(ref_file)
                    if os.path.isfile(cn2_file):
                        cn2_length=cn2_length_readin(dict_opts)
                    else:
                        [cn2_length,SamplingPercentage,whole_genome,len_genome]=cn2_region_write(cn2_file,ref)
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
                        genome_name=genome_name_readin(dict_opts)
                        NullPath=NullPath_SetUp(workdir,dict_opts)
                        path_BP=PathBP_SetUp(NullPath)
                        print 'temp files produced under: '+workdir
                        Script_Path=workdir+'reference_SVelter/'
                        for bamF in bam_files:
                            time1=time.time()
                            if ReadLength_Flag==0:
                                ReadLengthHash={}
                            outputfile=NullPath+bamF.split('/')[-1].replace('.'+bam_file_appdix,'')+'.'+genome_name+'.null'
                            fo=open(outputfile,'w')
                            print >>fo, ' '.join(['position','GCContent','ReadDepth','SplitReads','AbnormalDirection','ThroughBP'])
                            fo.close()
                            SplitLength={}
                            InsertLength={}
                            fcn2=open(cn2_file)
                            chr_cn2=[]
                            cn2_regions=[]
                            while True:
                                pcn2=fcn2.readline().strip().split()
                                if not pcn2: break
                                if not len(pcn2)==3: break
                                if pcn2[0] in chromos:
                                    if not pcn2[0] in chr_cn2:
                                        chr_cn2.append(pcn2[0])
                                    if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<10**6:
                                        if not random.choice(range(100))>SamplingPercentage*100:
                                            cn2_regions.append(pcn2)
                            fcn2.close()
                            if chr_cn2==[]:
                                whole_genome=chromos_read_in(ref_file)
                                cn2_regions=random_pick_cn2_region(cn2_file,whole_genome,chromos,len_genome,dict_opts)
                                chr_cn2=chromos
                            for pcn2 in cn2_regions:
                                freadpairs=os.popen('''samtools view -F 256 %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1])+100,int(pcn2[2])-100))
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
                            if not chr_cn2==[]:
                                SplitLenPNum=SplitLenPNum_Calculate(SplitLength,NullPath,bamF,bam_file_appdix,genome_name,NullSplitLen_perc)
                                TotalILNum=IL_Stat_Calculate(InsertLength)
                                NullILCI=NullILCI_Calculate(InsertLength,TotalILNum,NullILCIs)
                                Window_Size=int(float(NullILCI[0])/3)
                                cn2_length=max([cn2_length,NullILCI[2]])
                                cn2_max_len=max(cn2_length*100,10**6)
                                ILNullDensity={}
                                RDNullDensity={}
                                DRNullDensity={}
                                TBNullDensity={}
                                SplitNullDensity={}
                                GC_Content={}
                                #fcn2=open(cn2_file)
                                #while True:
                                for pcn2 in cn2_regions:
                                    #pcn2=fcn2.readline().strip().split()
                                    #if not pcn2: break
                                    if not len(pcn2)==3: break
                                    if pcn2[0] in chromos:
                                        #if (int(pcn2[2])-int(pcn2[1]))>cn2_length and (int(pcn2[2])-int(pcn2[1]))<cn2_max_len:
                                        #    if random.choice(range(100))<SamplingPercentage*100:
                                                pos=[int(pcn2[1])+i*Window_Size for i in range(int(float(int(pcn2[2])-int(pcn2[1]))/Window_Size))]
                                                pos0=[pcn2[0]+'_'+str(i*Window_Size) for i in pos]
                                                RDNull=[0 for i in pos]
                                                SplitNull=[0 for i in pos]
                                                DRNull=[0 for i in pos]
                                                TBNull=[0 for i in pos]
                                                ILNull=[0 for i in pos]
                                                readInf=[]
                                                freadpairs=os.popen('''samtools view -F 256 %s %s:%d-%d'''%(bamF,pcn2[0],int(pcn2[1]),int(pcn2[2])))
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
                                #fcn2.close()
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
                                            if int(pbRD[1]) in GC_Content.keys():
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
                                    [TotalTBNum,TBNullDensity]=TBNullDensity_CleanUP(TBNullDensity)
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
                                    InsertLenNullfigure1='.'.join(InsertLenNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    InsertLenNullfigure2='.'.join(InsertLenNullTemp.split('.')[:-1])+'.2.pdf'
                                    if KeepFigure in ['no','N','No','n']:
                                        InsertLenNullfigure1=InsertLenNullfigure1.replace('.pdf','.na')
                                        InsertLenNullfigure2=InsertLenNullfigure2.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,InsertLenNullTemp,InsertLenNullfigure1,BoxPlotColor,InsertLenNullfigure2))
                                    DRNullTemp=NullPath+'DirectionNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fDR=open(DRNullTemp,'w')
                                    for dr in DRNullDensity.keys():
                                            print >> fDR, ' '.join([str(dr),str(DRNullDensity[dr])])
                                    fDR.close()
                                    DRNullfigure1='.'.join(DRNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    DRNullfigure2='.'.join(DRNullTemp.split('.')[:-1])+'.2.pdf'
                                    if KeepFigure in ['no','N','No','n']:
                                        DRNullfigure1=DRNullfigure1.replace('.pdf','.na')
                                        DRNullfigure2=DRNullfigure2.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,DRNullTemp,DRNullfigure1,BoxPlotColor,DRNullfigure2))
                                    SplitNullTemp=NullPath+'SplitNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fSP=open(SplitNullTemp,'w')
                                    for sp in SplitNullDensity.keys():
                                            print >> fSP, ' '.join([str(sp),str(SplitNullDensity[sp])])
                                    fSP.close()
                                    SplitNullfigure1='.'.join(SplitNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    SplitNullfigure2='.'.join(SplitNullTemp.split('.')[:-1])+'.2.pdf'
                                    if KeepFigure in ['no','N','No','n']:
                                        SplitNullfigure1=SplitNullfigure1.replace('.pdf','.na')
                                        SplitNullfigure2=SplitNullfigure2.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s'''%(RFigureDRSplit,SplitNullTemp,SplitNullfigure1,BoxPlotColor,SplitNullfigure2))
                                    if model_comp=='C':
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.b.r' 
                                    else:
                                        RFigureDRSplit2=Script_Path+'SVelter1.NullModel.Figure.c.r'    
                                    RDNullTemp=NullPath+'RDNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fRD=open(RDNullTemp,'w')
                                    for rd in RD_Af_Adj.keys():
                                            print >> fRD, ' '.join([str(rd),str(RD_Af_Adj[rd])])
                                    fRD.close()
                                    RDNullfigure1='.'.join(RDNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    RDNullfigure2='.'.join(RDNullTemp.split('.')[:-1])+'.NegativeBinomial'
                                    if KeepFigure in ['no','N','No','n']:
                                        RDNullfigure1=RDNullfigure1.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s %s %d'''%(RFigureDRSplit2,RDNullTemp,RDNullfigure1,BoxPlotColor,lineColor,RDNullfigure2,Window_Size))
                                    RDNullfigure2_Modify(RDNullfigure2,Window_Size)
                                    ILNullTemp=NullPath+'ILNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fIL=open(ILNullTemp,'w')
                                    for il in InsertLength.keys():
                                            print >> fIL, ' '.join([str(il),str(InsertLength[il])])
                                    fIL.close()
                                    ILNullfigure1='.'.join(ILNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    ILNullfigure2='.'.join(ILNullTemp.split('.')[:-1])+'.Bimodal'
                                    if KeepFigure in ['no','N','No','n']:
                                        ILNullfigure1=ILNullfigure1.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s %s %d'''%(RFigureDRSplit2,ILNullTemp,ILNullfigure1,BoxPlotColor,lineColor,ILNullfigure2,Window_Size))
                                    TBNullTemp=NullPath+'TBNull.'+'.'.join(bamF.split('/')[-1].split('.')[:-1])+'.'+genome_name+'.temp'
                                    fTB=open(TBNullTemp,'w')
                                    for tb in TBNullDensity.keys():
                                        print >> fTB, ' '.join([str(tb),str(TBNullDensity[tb])])
                                    fTB.close()
                                    TBNullfigure1='.'.join(TBNullTemp.split('.')[:-1]+['pdf'])
                                    BoxPlotColor='blue'
                                    lineColor='red'
                                    TBNullfigure2='.'.join(TBNullTemp.split('.')[:-1])+'.Bimodal'
                                    if KeepFigure in ['no','N','No','n']:
                                        TBNullfigure1=TBNullfigure1.replace('.pdf','.na')
                                    os.system('''Rscript %s %s %s %s %s %s %d'''%(RFigureDRSplit2,TBNullTemp,TBNullfigure1,BoxPlotColor,lineColor,TBNullfigure2,Window_Size))
                                    clean_files()
                                    Ref_Seq_File=ref_file
                                    Mini_CN2_Region=int(cn2_length)
                                    Length_Limit=int(cn2_length)
                                    CN2_Region={} #key of hash CN2_Region is the name of each chromosome
                                    for chrom in chromos:
                                            CN2_Region[chrom]={} #key of CN2_Region[chrom] is GC_content
                                            for con in range(101):
                                                    CN2_Region[chrom][con]=[]
                                    #fcn2=open(cn2_file)
                                    temp_Name='temp.Null1.'+bamF.split('/')[-1]
                                    #while True:
                                    for pcn2 in cn2_regions:
                                            #pcn2=fcn2.readline().strip().split()
                                            if not len(pcn2)==3: break
                                            Chromosome=pcn2[0]
                                            if Chromosome in CN2_Region.keys():
                                                #if int(pcn2[2])-int(pcn2[1])<Length_Limit: continue
                                                #if not int(pcn2[2])-int(pcn2[1])<Length_Limit:
                                                    fasta_file=NullPath+temp_Name+'.fa'
                                                    os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(Ref_Seq_File,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),fasta_file))                    
                                                    Seq1=Fasta_To_Sequence_nullmodel(fasta_file)
                                                    if Seq1=='ERROR!':continue
                                                    if not Seq1=='ERROR!':
                                                        sam_file=NullPath+temp_Name+'.sam'
                                                        os.system(r'''samtools view -F 256 %s %s:%d-%d > %s'''%(bamF,str(pcn2[0]),int(pcn2[1]),int(pcn2[2]),sam_file))
                                                        Number_Of_Windows=len(Seq1)/Window_Size
                                                        GC_Content={}
                                                        for i in range(len(Seq1)/Window_Size+1)[1:]:                
                                                                Seq2=Seq1[(i-1)*Window_Size:i*Window_Size]
                                                                GC_Content[i]=GC_Content_Calculate(Seq2)
                                                        coverage=Region_Coverage_Calculate(sam_file,Number_Of_Windows,pcn2,Window_Size)
                                                        for j in GC_Content.keys():
                                                            if j in coverage.keys():
                                                                CN2_Region[Chromosome][GC_Content[j][0]].append(coverage[j][-1])
                                    #fcn2.close()
                                    if os.path.isfile(NullPath+temp_Name+'.fa'):
                                        os.system(r'''rm %s'''%(NullPath+temp_Name+'.fa'))
                                    if os.path.isfile(NullPath+temp_Name+'.sam'):
                                        os.system(r'''rm %s'''%(NullPath+temp_Name+'.sam'))
                                    Output_File=NullPath+'RD_Stat/'+bamF.split('/')[-1].replace('.'+bam_file_appdix,'')+'.'+genome_name+'_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
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
                            print 'Time Consuming: '+str(time2-time1)
    if function_name=='BPSearch':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','out-path=','help=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration=','BPSPCff=','BPLNCff='])
        dict_opts=dict(opts)
        CN2_Region={}
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_bpsearch()
        else:
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
                    BPAlignQC=0.2
                global QCAlign   
                if '--qc-align' in dict_opts.keys():
                    QCAlign=int(dict_opts['--qc-align'])
                else:
                    QCAlign=20
                global QC_RDCalculate_Cff
                QC_RDCalculate_Cff=10
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
            def global_para_declaration():
                global BPPath
                global NullPath
                global workdir
                global bam_path
                global bam_files
                global bam_files_appdix
                global ref_path
                global ref_file
                global ref_index
                global Window_Size
                global ReadLength
                global sub_loc_size
                workdir=path_modify(dict_opts['--workdir'])
                bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                bam_files=[dict_opts['--sample']]
                bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                ref_path=workdir+'reference_SVelter/'
                ref_file=ref_path+'genome.fa'
                ref_index=ref_file+'.fai'
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
            Define_Default_BPSearching()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                if not '--sample' in dict_opts.keys():
                    print 'Error: please specify either input file using --sample'
                else:
                    global_para_declaration()
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        chromos=chromos_read_in(ref_file)
                        if '--chromosome' in dict_opts.keys():
                            chrom_single=dict_opts['--chromosome']
                            if not chrom_single in chromos:
                                print 'Error: please make sure the chromosome defined by --chr is correct based on the reference genome'
                                chromos=[]
                            else:
                                chromos=[chrom_single]
                        if not chromos==[]:
                            genome_name=genome_name_readin(dict_opts)
                            print 'temp files produced under: '+workdir
                            if '--out-path' in dict_opts.keys():
                                BPPath=path_modify(dict_opts['--out-path'])
                            else:
                                BPPath=workdir+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
                            if not os.path.isdir(BPPath):
                                os.system(r'''mkdir %s'''%(BPPath))
                            NullPath='/'.join(BPPath.split('/')[:-2])+'/'+'.'.join(['NullModel']+BPPath.split('/')[-2].split('.')[1:])+'/'
                            for bamF in bam_files:
                                time1=time.time()                          
                                [Window_Size,ReadLength,sub_loc_size]=Null_Stats_Readin_One(NullPath,bamF,NullSplitLen_perc,genome_name,bam_files_appdix)
                                for chrF in chromos:
                                    bamF_Name=bamF.split('/')[-1].replace('.'+bam_files_appdix,'')
                                    floc_Name=BPPath+bamF_Name+'.'+chrF
                                    Refloc_name='.'.join(ref_file.split('.')[:-1])+'.Mappable.bed'
                                    stat_file_name(bamF_Name,genome_name)
                                    if os.path.isfile(Refloc_name):
                                        BamInput=bamF
                                        ILStats=ILStats_readin(ILStat)
                                        RDStats=RDStats_readin(RDStat)
                                        TBStats=TBStats_readin(TBStat)
                                        SPLCff=SPLCff_Calculate(NullSplitLen_perc,SPLenStat,ReadLength)
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
                                        [InsertLenMin,SplitMin,DRMin,BPSPCff,BPLNCff,SPCluLen]=[5,5,5,3,3,5]
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
                                            #BPSPCff=int(round(0.8*float(TBStats['stat']['Median'])/float(10)))
                                            BPSPCff=int(round(2*float(RDStats['Median'])/float(10)))
                                        if BPSPCff<3:
                                            BPSPCff=3
                                        if '--BPLNCff' in dict_opts.keys():
                                            BPLNCff=int(float(dict_opts['--BPLNCff']))
                                        else:
                                            #BPLNCff=int(round(1.2*float(TBStats['stat']['Median'])/float(10)))
                                            BPLNCff=int(round(3*float(RDStats['Median'])/float(10)))
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
                                        ClusterLen=ClusterLen_Calculation(ILStats,model_comp,ReadLength)
                                        ClusterLen2=int(ClusterLen/Window_Size+1)*Window_Size
                                        Min_Distinguish_Len=Window_Size
                                        subLnClusterLen=ClusterLen/2
                                        if not '-S' in dict_opts.keys():
                                            dict_opts['-S']=3
                                        ILCffs=IL_CI_Decide(ILStats,int(dict_opts['-S']),model_comp)
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
                                        if not loc_rec=={} and chrF in loc_rec.keys():
                                            test_mul_RP=[]
                                            test_mul_SP=[]
                                            chrom=chrF
                                            mini_fout_Name=BPPath+bamF_Name+'.mini.'+chrom+'.sam'
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
                                                    loc2=split_loc_to_subloc(loc,sub_loc_size,ClusterLen2)
                                                    for real_region in loc2:
                                                        fmini=open(mini_fout_Name,'a')
                                                        fRDind=open(RD_index_File,'a')
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/Window_Size+1)]
                                                        fbam=os.popen('''samtools view -F 256 %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if not int(pbam[4])>int(QC_RDCalculate_Cff): continue             #fail quality control, skip
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect_flag(pbam[1])
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
                                                os.system(r'''samtools view -h -Sb -F 256 %s -o %s'''%(mini_fout_Name,mini_fout_N2))
                                                os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
                                                os.system(r'''samtools index %s '''%(mini_fout_N4))
                                                os.system(r'''rm %s'''%(mini_fout_N2))
                                                os.system(r'''rm %s'''%(mini_fout_Name))
                                            temp_IL_Rec={}
                                            Link_IL_Rec={}
                                            for loc in loc_rec[chrom]: 
                                                loc2=split_loc_to_loc2(loc,ClusterLen)
                                                for real_region in loc2:
                                                        fbam=os.popen('''samtools view -F 256 %s %s:%d-%d'''%(mini_fout_N4,chrom,real_region[0],real_region[1]))
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
                                                            if int(pbam[4])<int(QCAlign): continue             #fail quality control, skip
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect_flag(pbam[1])
                                                            ReadLen=cigar2reaadlength(pbam[5])
                                                            absIL=abs(int(pbam[8]))
                                                            signIL=0
                                                            posF=[]
                                                            if absIL<int(ILCffs[0]) or absIL>int(ILCffs[1]):
                                                                    signIL+=1
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
                                                            #else:
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
                                                                            #pos1=int(pbam[3])+ReadLen
                                                                            #pos2=int(pbam[7])
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
                                                                    #pos1=int(pbam[3])+ReadLen
                                                                    #pos2=int(pbam[7])
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
                                                                    #pos1=int(pbam[3])+ReadLen
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
                                                                    #pos1=int(pbam[3])
                                                                    #pos2=int(pbam[7])
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
                                                                                    #pos1=int(pbam[3])+ReadLen
                                                                                    LinkNM['F'].append(pos1)
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abInfF.keys():
                                                                                                    abInfF[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfF[pos1][2]+=1
                                                                            elif DRtemp[0]=='-':
                                                                                    #pos1=int(pbam[3])
                                                                                    LinkNM['R'].append(pos1)
                                                                                    if pos1>real_region[0] and pos1<real_region[1]:
                                                                                            if not pos1 in abInfR.keys():
                                                                                                    abInfR[pos1]=[0,0,1,0]
                                                                                            else:
                                                                                                    abInfR[pos1][2]+=1
                                                                    if not pbam[6]=='=':
                                                                            if DRtemp[0]=='+':
                                                                                    #pos1=int(pbam[3])+ReadLen
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
                                                                                    #pos1=int(pbam[3])
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
                                                                                            #new_core.append(max(k2x))
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
                                                                                            #new_core.append(min(k2x))
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
                                                        if not abInfF=={}:
                                                                clu_a_F=clusterNums(abInfF.keys(), ClusterLen, 'f')[0]
                                                                if not clu_a_F==[]:
                                                                        clu_b_F=clusterNums(abInfF.keys(), ClusterLen, 'r')[0]
                                                                        clu_c_F=clusterSupVis2(sorted(abInfF.keys()), clu_b_F, [caf+10 for caf in clu_a_F],'right')
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
                                                        if not abInfR=={}:
                                                                clu_a_R=clusterNums(abInfR.keys(), ClusterLen, 'r')[0]
                                                                if not clu_a_R==[]:
                                                                        clu_b_R=clusterNums(abInfR.keys(), ClusterLen, 'f')[0]
                                                                        clu_c_R=clusterSupVis2(sorted(abInfR.keys()), [car-10 for car in clu_a_R], clu_b_R,'left')
                                                                        if not clu_c_R=={}:
                                                                                for key2 in clu_c_R.keys():  
                                                                                        key2b=key2+10
                                                                                        record=0
                                                                                        record2=0
                                                                                        for key3 in clu_c_R[key2]:
                                                                                                if not key3<key2b:
                                                                                                        record+= sum(abInfR[key3][:3])   
                                                                                                if abs(key3-key2b)<SPCluLen:
                                                                                                        record2+=abInfR[key3][3]
                                                                                        if not record+record2<LinkCluMin:
                                                                                                tempIL[key2b]=[[record,record2,'r'],clu_c_R[key2]]
                                                                                        del clu_c_R[key2]
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
                                                        LinkIL={}
                                                        for k1 in tempIL.keys():
                                                            temp_mate_F={}
                                                            temp_mate_R={}
                                                            info_mate=0
                                                            if tempIL[k1][0][2]=='f':
                                                                for k2 in tempIL[k1][1]:
                                                                    if k2 in LinkFF.keys():
                                                                        for k3 in LinkFF[k2]:
                                                                            if k3 in abInfF.keys():
                                                                                temp_mate_F[k3]=sum(abInfF[k3])
                                                                                #info_mate+=sum(abInfF[k3])
                                                                    if k2 in LinkFR.keys():
                                                                        for k3 in LinkFR[k2]:
                                                                            if k3 in abInfR.keys():
                                                                                temp_mate_R[k3]=sum(abInfR[k3])
                                                                                #info_mate+=sum(abInfR[k3])
                                                            elif tempIL[k1][0][2]=='r':
                                                                for k2 in tempIL[k1][1]:
                                                                    if k2 in LinkRF.keys():
                                                                        for k3 in LinkRF[k2]:
                                                                            if k3 in abInfF.keys():
                                                                                temp_mate_F[k3]=sum(abInfF[k3])
                                                                                #info_mate+=sum(abInfF[k3])
                                                                    if k2 in LinkRR.keys():
                                                                        for k3 in LinkRR[k2]:
                                                                            if k3 in abInfR.keys():
                                                                                temp_mate_R[k3]=sum(abInfR[k3])
                                                                                #info_mate+=sum(abInfR[k3])
                                                            for k1x in temp_mate_F.keys():
                                                                info_mate+=temp_mate_F[k1x]
                                                            for k1x in temp_mate_R.keys():
                                                                info_mate+=temp_mate_R[k1x]
                                                            if not info_mate<LinkCluMin:
                                                                LinkIL[k1]=[[],[]]
                                                                if not temp_mate_F=={}:
                                                                    LinkIL[k1][0]=clusterQC(clusterNums4(temp_mate_F, ClusterLen, 'f'),LinkCluMin)
                                                                if not temp_mate_R=={}:
                                                                    LinkIL[k1][1]=clusterQC(clusterNums4(temp_mate_R, ClusterLen, 'r'),LinkCluMin)
                                                                else:continue
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
                                                            out_pair_bp_temp=out_pair_bp_check(out_pair_bp,3)
                                                            out_pair_bp=out_pair_bp_temp
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
                                                            out_pair_modify=out_pair_modify_check(out_pair_modify,3)
                                                            for i in sorted(out_pair_modify.keys()):
                                                                for j in out_pair_modify[i]:                            
                                                                    print >>fout, ' '.join([chrom,str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                                    #print ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
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
                                                                print >>fout, ' '.join([str(j) for j in [chrom,i,num]])
                                                                #print ' '.join([str(j) for j in [i,num]])
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
                                                        print >>fout, ' '.join([chrom,str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                        #print ' '.join([str(i),str(out_pair_numrec[i][0]),str(j),str(out_pair_numrec[i][out_pair_modify[i].index(j)+1])])
                                                        if i in tempIL.keys():
                                                            del tempIL[i]
                                                        if j in tempIL.keys():
                                                            del tempIL[j]
                                                fout.close()
                                            fout=open(BPOutputa,'a')
                                            for i in sorted(temp_IL_Rec.keys()):
                                                print >>fout, ' '.join([chrom,str(i),str(sum(temp_IL_Rec[i][0][:2]))])    
                                            fout.close()
                                            time2=time.time()
                                            LN_Filter(BPOutputb,BPOutputa,workdir)
                                            os.system(r'''cat %s >> %s'''%(BPOutputd,BPOutpute))
                                            os.system(r'''rm %s'''%(BPOutputd))
                                            print 'BPSearch Complete for '+bamF+'.'+chrF
                                            print 'Time Consuming: '+str(time2-time1)
    if function_name=='BPSearch_Predefined':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','out-path=','input-bed=','help=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration=','BPSPCff=','BPLNCff='])
        dict_opts=dict(opts)
        CN2_Region={}
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            print 'SVelter-0.1          Last Update:2015-08-20'
            print 'Required Parameters:'
            print '    --sample, input alignment file in bam format'
            print '    --workdir, writable working directory.'
            print 'Optional Parameters:'
            print '    --chromosome, name of chromosome to run. should match chromosome name in bam file'
            print '    --null-model, specify which stat model to be fitted on each parameter. if --null-model==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
            print '    --qc-align, minimum alignment quality required for mapped reads in bam file [20)'
            print '    --qc-split, minimum alighment of clipped parts of reads considered as a soft clip [20)'
            print '    --split-min-len, the minumum length of clip read considered as split; (default:10% of read length)'
            print '    --qc-map-tool, the tool extracts mappability information from a bigWig file,avaliable from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary'
            print '    --qc-map-file, .bigWig file used to decide local genomic mappability, avaliable from: ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/'
            print '    --qc-map-cutoff, the minimum mapping quality required for a breakpoint to be reported [0.0)'
        else:
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
                    BPAlignQC=0.2
                global QCAlign   
                if '--qc-align' in dict_opts.keys():
                    QCAlign=int(dict_opts['--qc-align'])
                else:
                    QCAlign=20
                global QC_RDCalculate_Cff
                QC_RDCalculate_Cff=10
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
            def global_para_declaration():
                print 'temp files produced under: '+workdir
                global bam_path
                global bam_files
                global bam_files_appdix
                global ref_path
                global ref_file
                global ref_index
                global input_bed
                global Window_Size
                global ReadLength
                global sub_loc_size
                bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                bam_files=[dict_opts['--sample']]
                bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                ref_path=workdir+'reference_SVelter/'
                ref_file=ref_path+'genome.fa'
                ref_index=ref_file+'.fai'
                if '--input-bed' in dict_opts.keys():
                    input_bed=dict_opts['--input-bed']
                global BPPath
                if '--out-path' in dict_opts.keys():
                    BPPath=path_modify(dict_opts['--out-path'])
                else:
                    BPPath=workdir+'BreakPoints.'+'.'.join(dict_opts['--sample'].split('/')[-1].split('.')[:-1])+'.predefinedBP.'+'.'.join(input_bed.split('/')[-1].split('.')[:-1])+'/'
                path_mkdir(BPPath)
                global NullPath
                NullPath='/'.join(BPPath.split('/')[:-2])+'/'+'.'.join(['NullModel']+BPPath.split('/')[-2].split('.')[1:])+'/'
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
            Define_Default_BPSearching()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                workdir=path_modify(dict_opts['--workdir'])
                if not '--sample' in dict_opts.keys():
                    print 'Error: please specify either input file using --sample'
                else:
                    global_para_declaration()
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        chromos=chromos_read_in(ref_file)
                        if '--chromosome' in dict_opts.keys():
                            chrom_single=dict_opts['--chromosome']
                            if not chrom_single in chromos:
                                print 'Error: please make sure the chromosome defined by --chr is correct based on the reference genome'
                                chromos=[]
                            else:
                                chromos=[chrom_single]
                        if not chromos==[]:
                            genome_name=genome_name_readin(dict_opts)
                            for bamF in bam_files:
                                time1=time.time()
                                [Window_Size,ReadLength,sub_loc_size]=Null_Stats_Readin_One(NullPath,bamF,NullSplitLen_perc,genome_name,bam_files_appdix)
                                for chrF in chromos:
                                    bamF_Name=bamF.split('/')[-1].replace('.'+bam_files_appdix,'')
                                    floc_Name=BPPath+bamF_Name+'.'+chrF
                                    Refloc_name='.'.join(ref_file.split('.')[:-1])+'.Mappable.bed'
                                    stat_file_name(bamF_Name,genome_name)
                                    if os.path.isfile(Refloc_name):
                                        BamInput=bamF
                                        ILStats=ILStats_readin(ILStat)
                                        RDStats=RDStats_readin(RDStat)
                                        TBStats=TBStats_readin(TBStat)
                                        SPLCff=SPLCff_Calculate(NullSplitLen_perc,SPLenStat,ReadLength)
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
                                            #BPSPCff=int(round(0.8*float(TBStats['stat']['Median'])/float(10)))
                                            BPSPCff=int(round(2*float(RDStats['Median'])/float(10)))
                                        if BPSPCff<3:
                                            BPSPCff=3
                                        if '--BPLNCff' in dict_opts.keys():
                                            BPLNCff=int(float(dict_opts['--BPLNCff']))
                                        else:
                                            #BPLNCff=int(round(1.2*float(TBStats['stat']['Median'])/float(10)))
                                            BPLNCff=int(round(3*float(RDStats['Median'])/float(10)))
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
                                        ClusterLen=ClusterLen_Calculation(ILStats,model_comp,ReadLength)
                                        ClusterLen2=int(ClusterLen/Window_Size+1)*Window_Size
                                        Min_Distinguish_Len=Window_Size
                                        subLnClusterLen=ClusterLen/2
                                        if not '-S' in dict_opts.keys():
                                            dict_opts['-S']=3
                                        ILCffs=IL_CI_Decide(ILStats,int(dict_opts['-S']),model_comp)
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
                                        if not loc_rec=={} and chrF in loc_rec.keys():
                                            test_mul_RP=[]
                                            test_mul_SP=[]
                                            chrom=chrF
                                            RD_index_Path=NullPath+'RD_Stat/'
                                            path_mkdir(RD_index_Path)
                                            RD_index_File=RD_index_Path+bamF_Name+'.'+chrF+'.RD.index'
                                            if not os.path.isfile(RD_index_File):
                                                file_setup(RD_index_File)
                                                for loc in loc_rec[chrom]:
                                                    loc2=split_loc_to_subloc(loc,sub_loc_size,ClusterLen2)
                                                    for real_region in loc2:
                                                        fRDind=open(RD_index_File,'a')
                                                        print >>fRDind,chrom+':'+str(real_region[0])+'-'+str(real_region[1])                
                                                        RD_RealRegion=[0 for i in range((real_region[1]-real_region[0])/Window_Size+1)]
                                                        fbam=os.popen('''samtools view -F 256 %s %s:%d-%d'''%(BamInput,chrom,real_region[0],real_region[1]))
                                                        while True:
                                                            pbam1=fbam.readline().strip()
                                                            if not pbam1: break
                                                            pbam=pbam1.split()
                                                            if not int(pbam[4])>int(QC_RDCalculate_Cff): continue             #fail quality control, skip
                                                            if int(pbam[1])&4>0: continue           #the read was not mapped, skip
                                                            DRtemp=Reads_Direction_Detect_flag(pbam[1])
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
                                                        fbam.close()
                                                        for rfrr in range(len(RD_RealRegion[:-1])):
                                                            RD_RealRegion[rfrr]=str(float(RD_RealRegion[rfrr])/Window_Size)
                                                        if real_region[1]-real_region[0]-(real_region[1]-real_region[0])/Window_Size*Window_Size ==0:
                                                            del RD_RealRegion[-1]
                                                        else:
                                                            RD_RealRegion[-1]=str(float(RD_RealRegion[-1])/float(real_region[1]-real_region[0]-(real_region[1]-real_region[0])/Window_Size*Window_Size))
                                                        print >>fRDind, ' '.join(RD_RealRegion)
                                                        fRDind.close()
    if function_name=='BPIntegrate':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','help=','bp-path=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_bpintegrate()
        else:
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
            def global_para_declaration():  
                global workdir  
                workdir=path_modify(dict_opts['--workdir'])
                print 'temp files produced under: '+workdir
                global bps_in_path
                if not '--bp-path' in dict_opts.keys():
                    bps_in_path=workdir+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
                else:
                    bps_in_path=path_modify(dict_opts['--bp-path'])
                if not bps_in_path[0]=='/' and not bps_in_path[:2]=='./':
                    bps_in_path='./'+bps_in_path
                global S_Sample
                global chromo_name
                global LN_list
                global all_SPs
            min_length=100
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
            Define_Default_BPIntegrate()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                global_para_declaration()
                if not '--sample' in dict_opts.keys():
                    print 'Error: please specify either input file using --sample'
                else:
                    bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                    bam_files=[dict_opts['--sample']]
                    bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                    bam_names=[dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'')]
                    ref_path=workdir+'reference_SVelter/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        chromos=chromos_read_in(ref_file)
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
                            bps_folder='/'.join(bps_in_path.split('/')[:-2])+'/'+'.'.join(['bp_files']+bps_in_path.split('/')[-2].split('.')[1:])+'/'
                            path_mkdir(bps_folder)
                            for i in bam_names:
                                bps_hash[i]={}
                                for file1 in os.listdir(bps_in_path):
                                    if file1.split('.')[-1]=='LNs':
                                        key_bps_hash='.'.join(file1.split('.')[:-1])
                                        if not key_bps_hash in bps_hash[i].keys():
                                            bps_hash[i][key_bps_hash]=[]
                                        bps_hash[i][key_bps_hash].append(file1)
                                        bps_hash[i][key_bps_hash].append(key_bps_hash+'.SPs')
                            for S_Sample in bps_hash.keys():
                                [LN_list,all_SPs]=[{},{}]
                                for chromo_name in bps_hash[S_Sample].keys():
                                    [LN_list,all_SPs]=SP_LN_info_ReadIn(LN_list,all_SPs,bps_in_path+bps_hash[S_Sample][chromo_name][0],bps_in_path+bps_hash[S_Sample][chromo_name][1])
                                for chromo_name in LN_list.keys():
                                    if not LN_list[chromo_name]==[]:
                                        if chromo_name in all_SPs.keys():
                                            unique_SPs=SP_Info_Merge(all_SPs[chromo_name])
                                        else:
                                            all_SPs[chromo_name]={}
                                            unique_SPs=[]
                                        modified_LNs=LN_Info_Correct(LN_list[chromo_name],unique_SPs)
                                        multi_removed_LNs=multi_trans_detect(modified_LNs)
                                        LN_LN_Merge_0=merge_LNs_into_LNs(multi_removed_LNs)
                                        LN_LN_Merge=LN_Merge_Final_Check(LN_LN_Merge_0)
                                        if not '--batch' in dict_opts.keys():
                                            write_bp_1a(LN_LN_Merge,bps_folder,chromo_name,S_Sample)
                                        else:
                                            if dict_opts['--batch']=='0':
                                                write_bp_2a(LN_LN_Merge,bps_folder,chromo_name,S_Sample)
                                            else:
                                                file_length=int(dict_opts['--batch'])
                                                file_index=write_bp_3a(LN_LN_Merge,bps_folder,file_length,chromo_name,S_Sample)
                                LN_bps_write(bps_hash,bps_folder,S_Sample,dict_opts,chromos,allchromos,bps_in_path)
                            time2=time.time()
                            print 'BPIntegrate Complete !'
                            print 'Time Consuming: '+str(time2-time1)
    if function_name=='SVPredict':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','help=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','input-bed=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_svpredict()
        else:
            import os
            import sys
            import getopt
            import random
            import scipy
            import math
            import numpy
            import pickle
            from math import sqrt,pi,exp
            import scipy
            from scipy.stats import norm
            import time
            import datetime
            import itertools
            def Af_Rearrange_Info_Collect(GC_para_dict,BP_para_dict,Be_BP_Letter,Be_Info,Letter_Candidates):
                P_IL=[]
                P_RD=[]
                P_DR=[]
                P_TB=[]
                Letter_Rec=[]
                BP_Rec=[]
                for Af_Letter in Letter_Candidates:
                    Af_BP=[[BP_para_dict['original_bp_list'][0]],[BP_para_dict['original_bp_list'][0]]]
                    for i in Af_Letter[0]:
                            Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
                    for i in Af_Letter[1]:
                            Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
                    Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                    if not Af_Info_all==0:
                        Letter_Rec.append(Af_Letter)
                        BP_Rec.append(Af_BP)
                        Af_IL_Penal=Af_Info_all[0]
                        Af_RD_Rec=Af_Info_all[1]
                        Af_DR_Penal=(Af_Info_all[2])**2
                        Af_TB_Penal_a=Af_Info_all[4]
                        Af_TB_Penal=Af_TB_Penal_Less_Important_Caldu(BP_para_dict['num_of_reads'],Af_Info_all,Af_Letter,GC_para_dict['IL_Statistics'],BP_para_dict['ReadLength'],GC_para_dict['GC_Overall_Median_Num'])
                        Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                        for key in Af_Info_all[5].keys():
                            Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_para_dict['GC_Var_Coverage'][chrom_N]/2)
                        #P_IL.append(Af_IL_Penal-IL_GS)
                        #P_RD.append((Af_RD_Penal-RD_GS)*(len(Af_Letter[0])+len(Af_Letter[1])+1))
                        P_IL.append(Af_IL_Penal)
                        P_RD.append(Af_RD_Penal)
                        P_DR.append(Af_DR_Penal/num_of_read_pairs)
                        P_TB.append(Af_TB_Penal)
                if P_IL==[]: return 'Error'
                else:
                    Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                    Regu_IL=[i*K_IL_new for i in Regu_IL]
                    Regu_RD=[P_RD[i]*(1-TB_Weight*P_TB[i]) for i in range(len(P_RD))] 
                    return [Regu_IL,Regu_RD,Letter_Rec,BP_Rec]
            def Af_Rearrange_Info_Collect_2(BP_para_dict,Letter_Candidates):
                P_IL=[]
                P_RD=[]
                P_DR=[]
                P_TB=[]
                Letter_Rec=[]
                BP_Rec=[]
                for Af_Letter in Letter_Candidates:
                    Af_BP=[[BP_para_dict['original_bp_list'][0]],[BP_para_dict['original_bp_list'][0]]]
                    for i in Af_Letter[0]:
                            Af_BP[0].append(Af_BP[0][-1]+Be_BP_Letter[i])
                    for i in Af_Letter[1]:
                            Af_BP[1].append(Af_BP[1][-1]+Be_BP_Letter[i])
                    Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                    if not Af_Info_all==0:
                        Letter_Rec.append(Af_Letter)
                        BP_Rec.append(Af_BP)
                        Af_IL_Penal=Af_Info_all[0]
                        Af_RD_Rec=Af_Info_all[1]
                        Af_DR_Penal=(Af_Info_all[2])**2
                        Af_TB_Penal_a=Af_Info_all[4]
                        Af_TB_Rec=Af_Info_all[3]
                        Af_TB_Penal=-float(Af_TB_Penal_a)/float(BP_para_dict['num_of_reads'])+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                        Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                        for key in Af_Info_all[5].keys():
                                Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_para_dict['GC_Var_Coverage'][chrom_N]/2)
                        P_IL.append(Af_IL_Penal)
                        P_RD.append(Af_RD_Penal)
                        P_DR.append(Af_DR_Penal/num_of_read_pairs)
                        P_TB.append(Af_TB_Penal)
                if P_IL==[]: return 'Error'
                else:
                    Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                    Regu_IL=[i*K_IL_new for i in Regu_IL]
                    Regu_RD=[P_RD[i]*(1-TB_Weight*P_TB[i]) for i in range(len(P_RD))] 
                    return [Regu_IL,Regu_RD,Letter_Rec,BP_Rec]
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
                RDm=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[0],flank))/2]+[0 for j in range(len(Initial_block_RD)-1),Window_Size]
                RDp=[(Initial_block_RD[0]+left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj[1],flank))/2]+[0 for j in range(len(Initial_block_RD)-1),Window_Size]
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
                for j in be_info_2:
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
                    j_m_3a=candidate_QC_Control(j_m_new)
                    if j_m_3a==[]:
                        jMapPenam+=1
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
                    j_p_3a=candidate_QC_Control(j_p_new)
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
            def Be_Info_3_rearrange(BP_para_dict,Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal):
                be_info_3=Be_Info[2]
                for j in be_info_3:
                    j_m_new=[]
                    if j[0] in temp_letter[0] and j[2] in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]]:
                            for kb in Let_BP_Info['m'][j[2]]:
                                temp_single=[ka[0]+j[1],kb[0]+j[3]]
                                if not temp_single[1]-temp_single[0]>BP_para_dict['ReadLength']*1.2 and temp_single[1]-temp_single[0]>0:
                                    j_m_new.append(temp_single)
                    if j[0]+'^' in temp_letter[0] and j[2]+'^' in temp_letter[0]:
                        for ka in Let_BP_Info['m'][j[0]+'^']:
                            for kb in Let_BP_Info['m'][j[2]+'^']:
                                temp_single=[kb[1]-j[3],ka[1]-j[1]]
                                if not temp_single[1]-temp_single[0]>BP_para_dict['ReadLength']*1.2 and temp_single[1]-temp_single[0]>0:
                                    j_m_new.append(temp_single)
                    j_p_new=[]
                    if j[0] in temp_letter[1] and j[2] in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]]:
                            for kb in Let_BP_Info['p'][j[2]]:
                                temp_single=[ka[0]+j[1],kb[0]+j[3]]
                                if not temp_single[1]-temp_single[0]>BP_para_dict['ReadLength']*1.2 and temp_single[1]-temp_single[0]>0:
                                    j_p_new.append(temp_single)
                    if j[0]+'^' in temp_letter[1] and j[2]+'^' in temp_letter[1]:
                        for ka in Let_BP_Info['p'][j[0]+'^']:
                            for kb in Let_BP_Info['p'][j[2]+'^']:
                                temp_single=[kb[1]-j[3],ka[1]-j[1]]
                                if not temp_single[1]-temp_single[0]>BP_para_dict['ReadLength']*1.2 and temp_single[1]-temp_single[0]>0:
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
            def Block_Assign_To_Letters(bp_list,letter_list,flank):
                #Eg of bp_list:[184569179, 184569775, 184571064, 184572009, 184572016]
                #Eg of letter_list:['a', 'b', 'c', 'd']
                #Eg of flank:446
                number_of_blocks=(numpy.max(bp_list)-numpy.min(bp_list)+2*flank)/Window_Size+1
                blocks={}
                bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
                relative_bp_list=[i-numpy.min(bp_list_new) for i in bp_list_new]
                bp_length=[(bp_list_new[i+1]-bp_list_new[i]) for i in range(len(bp_list_new)-1)]
                letter_list_new=['left']+letter_list+['right']
                bp_blocks=[[letter_list_new[j]]+range(relative_bp_list[j]/Window_Size,relative_bp_list[j+1]/Window_Size+1) for j in range(len(relative_bp_list)-1)]
                blocks_bp={}
                for i in range(number_of_blocks):
                    blocks_bp[i+1]=[bp_list_new[0]+i*Window_Size,bp_list_new[0]+i*Window_Size+Window_Size-1]
                    for j in bp_blocks:
                        if i in j:
                            blocks_bp[i+1].append(j[0])
                blocks_bp[0]=[blocks_bp[1][0]-Window_Size,blocks_bp[1][0]-1,'0']
                blocks_bp[number_of_blocks+1]=[blocks_bp[number_of_blocks][1]+1,blocks_bp[number_of_blocks][1]+Window_Size,'0']
                return blocks_bp
            def block_Info_ReadIn(GC_para_dict,BP_para_dict,chr_letter_bp,blocks_read_in,Multi_Dup):
                block_bps={}
                block_rds={}
                for k1 in chr_letter_bp.keys():
                    block_bps[k1]={}
                    block_rds[k1]={}
                    for k2 in chr_letter_bp[k1].keys():
                        if not k2 in Multi_Dup:
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
                        multi_dup_flag=multi_dup_check(k2,Multi_Dup)
                        if multi_dup_flag==0:
                            k2a=[]
                            k2b=[]
                            for k3 in k2:
                                if type(k3)==type(1):
                                    k2a.append(k3)
                                else:
                                    k2b.append(k3)
                            fbam=os.popen(r'''samtools view -F 256 %s %s:%d-%d'''%(Initial_Bam,k1,min(k2a)-BP_para_dict['flank'],max(k2a)+BP_para_dict['flank']))
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
                                #if not int(pbam[4])>QCAlign:continue
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
                                    pos_block_assign(block_bps[k1],read_pos,tolerance_bp)
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
                                    #elif not mate_pos[1]<flank_region[0] and not mate_pos[0]>flank_region[1]:
                                    #   del_flag+=1
                                    if del_flag>0:
                                        del temp_rec[k3]
                                        pos_block_assign(block_bps[k1],read_pos,tolerance_bp)
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
                                    read_pos=[int(temp_rec[k3][0][2]),int(temp_rec[k3][0][2])+cigar2reaadlength(temp_rec[k3][0][4]),int(temp_rec[k3][1][2]),int(temp_rec[k3][1][2])+cigar2reaadlength(temp_rec[k3][1][4])]+Reads_Direction_Detect_flag(temp_rec[k3][0][0])
                                    #print temp_rec[k3]
                                    #if k3 in test2:
                                    #   print read_pos
                                    if read_pos[0]>read_pos[2]:
                                        read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
                                    pos_block_assign(block_bps[k1],read_pos,tolerance_bp)
                                    if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
                                        block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                                        block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
                                    elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
                                        Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
                                    else:
                                        Double_Read_ThroughBP[k1].append(read_pos)
                                    del temp_rec[k3]
                                    #if k3 in test2:
                                    #   print read_pos
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
                            pos_block_assign(block_bps[k1],read_pos,tolerance_bp)
                            if read_pos[-1]==read_pos[-2]:
                                block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                            else:
                                Single_Read_ThroughBP[k1].append(read_pos)
                    elif len(total_rec[k3])==2:
                        read_pos=[int(total_rec[k3][0][2]),int(total_rec[k3][0][2])+cigar2reaadlength(total_rec[k3][0][4]),int(total_rec[k3][1][2]),int(total_rec[k3][1][2])+cigar2reaadlength(total_rec[k3][1][4])]+Reads_Direction_Detect_flag(total_rec[k3][0][0])
                        #print read_pos
                        if read_pos[0]>read_pos[2]:
                            read_pos=read_pos[2:4]+read_pos[:2]+[read_pos[-1],read_pos[-2]]
                        pos_block_assign(block_bps[k1],read_pos,tolerance_bp)
                        if read_pos[6]==read_pos[7]==read_pos[8]==read_pos[9]:
                            block_rds[k1][read_pos[-1]]+=read_pos[1]-read_pos[0]
                            block_rds[k1][read_pos[-1]]+=read_pos[3]-read_pos[2]
                        elif read_pos[8]==read_pos[9] and read_pos[6]==read_pos[7]:
                            Pair_ThroughBP[k1].append(read_pos[:6]+[read_pos[6],read_pos[8]])
                        else:
                            Double_Read_ThroughBP[k1].append(read_pos)
                        del total_rec[k3]
                #print total_rec
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
                            #if -i[2]+block_bp2[i2][i[8]][1]>200 and i[8]=='a':
                                #print i
                            #if i[3]-block_bp2[i2][i[9]][0]>200 and i[9]=='a':
                                #print i
                        elif i[8]==i[9]:
                            block_rd2[i2][i[8]]+=i[3]-i[2]
                            block_rd2[i2][i[6]]+=-i[0]+block_bp2[i2][i[6]][1]
                            block_rd2[i2][i[7]]+=i[1]-block_bp2[i2][i[7]][0]
                            #if -i[0]+block_bp2[i2][i[6]][1]>101:
                                #print i
                            #if i[1]-block_bp2[i2][i[7]][0]>101:
                                #print i
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
                total_rd_calcu(GC_para_dict['GC_Median_Num'],GC_para_dict['GC_Overall_Median_Num'],letter_RD2,letter_GC,chr_letter_bp,block_rd2)
            def block_RD_Calculate_2a(Initial_GCRD_Adj,original_bp_list,flank):
                allele_BP=[0]+[flank+j-original_bp_list[0] for j in original_bp_list]+[2*flank+original_bp_list[-1]-original_bp_list[0]]
                allele_Letter=['left']+[chr(97+i) for i in range(len(original_bp_list)-1)]
                allele_RD=[]
                for k in range(len(allele_Letter)):
                    length=allele_BP[k+1]-allele_BP[k]
                    block=[allele_BP[k],allele_BP[k+1]]
                    temp=[]
                    if not block[0]==block[0]/Window_Size*Window_Size:
                        blf=float((block[0]/Window_Size+1)*Window_Size-block[0])/Window_Size*Initial_GCRD_Adj[block[0]/Window_Size+1][3]
                        temp.append(blf)
                        for m in range(block[0]/Window_Size+2,block[1]/Window_Size+1):
                            temp.append(Initial_GCRD_Adj[m][3])
                        if not block[1]==block[1]/Window_Size*Window_Size:
                            brf=float(block[1]-block[1]/Window_Size*Window_Size)/Window_Size*Initial_GCRD_Adj[block[1]/Window_Size+1][3]
                            temp.append(brf)
                        allele_RD.append(numpy.sum(temp)/length*Window_Size)        
                    elif block[0]==block[0]/Window_Size*Window_Size:
                        for m in range(block[0]/Window_Size+1,block[1]/Window_Size+1):
                            temp.append(Initial_GCRD_Adj[m][3])
                        if not block[1]==block[1]/Window_Size*Window_Size:
                            brf=float(block[1]-block[1]/Window_Size*Window_Size)/Window_Size*Initial_GCRD_Adj[block[1]/Window_Size+1][3]
                            temp.append(brf)
                        allele_RD.append(numpy.sum(temp)/length*Window_Size)        
                return allele_RD
            def copy_num_estimate_calcu(GC_para_dict,BP_para_dict,bps2):
                chr_letter_bp=letter_rearrange(BP_para_dict['flank'],bps2)
                Initial_GCRD_Adj_pre=letter_RD_ReadIn(letter_RD_test_calcu(chr_letter_bp))
                global Initial_GCRD_Adj
                Initial_GCRD_Adj={}
                for k1 in Initial_GCRD_Adj_pre.keys():
                    for k2 in Initial_GCRD_Adj_pre[k1].keys():
                        Initial_GCRD_Adj[k2]=Initial_GCRD_Adj_pre[k1][k2]
                Initial_GCRD_Adj['left']=numpy.mean([GC_para_dict['GC_Mean_Coverage'][key_chr[0]] for key_chr in bps2])
                Initial_GCRD_Adj['right']=numpy.mean([GC_para_dict['GC_Mean_Coverage'][key_chr[0]] for key_chr in bps2])
                Copy_num_estimate={}
                for i in Initial_GCRD_Adj.keys():
                    if not i in ['left','right']:
                        Copy_num_estimate[i]=round(Initial_GCRD_Adj[i]*2/GC_para_dict['GC_Mean_Coverage'][Chr])
                        if Initial_GCRD_Adj[i]<float(GC_para_dict['GC_Mean_Coverage'][Chr])/10.0:
                            Copy_num_estimate[i]=-1
                Copy_num_Check=[]
                for CNE in Copy_num_estimate.keys():
                    if Copy_num_estimate[CNE]>4:
                        Copy_num_Check.append(CNE)
                return [Copy_num_estimate,Copy_num_Check]
            def calcu_chr_letter_bp_left(bps2):
                out={}
                for i in bps2:
                    if not i[0] in out.keys():
                        out[i[0]]={}
                    out[i[0]]['a']=[i[1]-1000,i[1]]
                return out
            def calcu_chr_letter_bp_right(bps2):
                out={}
                for i in bps2:
                    if not i[0] in out.keys():
                        out[i[0]]={}
                    out[i[0]]['a']=[i[-1],i[-1]+1000]
                return out
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
                region_length_new=(region_length/100+1)*100-2*flank
                Number_Of_Blocks=len(Ori_1_Seq)/100
                GC_Content={}
                for i in range(Number_Of_Blocks):
                    GC_Content[i+1]=GC_Content_Calculate(Ori_1_Seq[i*100:(i+1)*100])[0]
                return GC_Content
            def c_Coverage_Calculate_2a(Letter_Single,Letter_Double,Chromo,original_bp_list,original_letters,flank):
                letter_list=original_letters
                bp_list=[i-original_bp_list[0] for i in original_bp_list]
                bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
                coverage=Block_Assign_To_Letters(bp_list,letter_list,flank)
                for key in coverage.keys():
                        coverage[key].append(0)
                for key in Letter_Single.keys():
                        for i in Letter_Single[key]:
                                keynumL=(i[0]+flank)/Window_Size+1
                                keynumR=(i[1]+flank)/Window_Size+1
                                lenL=coverage[keynumL][1]-i[0]
                                lenR=i[1]-coverage[keynumR][0]+1
                                if lenL>lenR:
                                        coverage[keynumL][-1]+=1
                                else:
                                        coverage[keynumR][-1]+=1
                for key in Letter_Double.keys():
                    for i in Letter_Double[key]:
                        keynumL=(i[0]+flank)/Window_Size+1
                        keynumR=(i[1]+flank)/Window_Size+1
                        if keynumL in coverage.keys() and keynumR in coverage.keys():
                            lenL=coverage[keynumL][1]-i[0]
                            lenR=i[1]-coverage[keynumR][0]+1
                            if lenL>lenR:
                                    coverage[keynumL][-1]+=1
                            else:
                                    coverage[keynumR][-1]+=1
                        keynumL=(i[2]+flank)/Window_Size+1
                        keynumR=(i[3]+flank)/Window_Size+1
                        if keynumL in coverage.keys() and keynumR in coverage.keys():
                            lenL=coverage[keynumL][1]-i[0]
                            lenR=i[1]-coverage[keynumR][0]+1
                            if lenL>lenR:
                                    coverage[keynumL][-1]+=1
                            else:
                                    coverage[keynumR][-1]+=1
                return coverage
            def c_Coverage_Calculate_2b(Letter_Through,Chromo,original_bp_list,original_letters,flank):
                #Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
                letter_list=original_letters
                bp_list=[i-original_bp_list[0] for i in bp_MP[0]]
                bp_list_new=[bp_list[0]-flank]+bp_list+[bp_list[-1]+flank]
                coverage=Block_Assign_To_Letters(bp_list,letter_list,flank)
                for key in coverage.keys():
                    coverage[key].append(0)
                for key in Letter_Through.keys():
                    i=Letter_Through[key]
                    keynumL=(i[0]+flank)/Window_Size+1
                    keynumR=(i[1]+flank)/Window_Size+1
                    lenL=coverage[keynumL][1]-i[0]
                    lenR=i[1]-coverage[keynumR][0]+1
                    if lenL>lenR:
                        coverage[keynumL][-1]+=1
                    elif lenL<lenR:
                        coverage[keynumR][-1]+=1
                    elif lenL==lenR:
                        coverage[keynumL][-1]+=0.5
                        coverage[keynumR][-1]+=0.5
                    keynumL=(i[2]+flank)/Window_Size+1
                    keynumR=(i[3]+flank)/Window_Size+1
                    lenL=coverage[keynumL][1]-i[0]
                    lenR=i[1]-coverage[keynumR][0]+1
                    if lenL>lenR:
                        coverage[keynumL][-1]+=1
                    elif lenL<lenR:
                        coverage[keynumR][-1]+=1
                    elif lenL==lenR:
                        coverage[keynumL][-1]+=0.5
                        coverage[keynumR][-1]+=0.5
                return coverage
            def c_Coverage_Calculate_2d(Full_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
                #Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
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
                for key in Full_Info.keys():
                    if not len(Full_Info[key])==8:
                        Halfa=Full_Info[key][:2]+[Full_Info[key][4]]+[Full_Info[key][6]]
                        Halfb=Full_Info[key][2:4]+[Full_Info[key][5]]+[Full_Info[key][6]]
                        for Half in [Halfa,Halfb]:
                            if Half[0]<-flank-Window_Size: continue
                            else:
                                if Half[-1]=='M':
                                    M_coverage[(Half[0]-(M_New_bp[0]))/Window_Size+1][-1]+=1
                                elif Half[-1]=='P':
                                    P_coverage[(Half[0]-(P_New_bp[0]))/Window_Size+1][-1]+=1
                    elif  len(Full_Info[key])==8:
                        Halfa=Full_Info[key][:2]+[Full_Info[key][4]]+[Full_Info[key][6]]
                        Halfb=Full_Info[key][2:4]+[Full_Info[key][5]]+[Full_Info[key][6]]
                        for Half in [Halfa,Halfb]:
                            if Half[0]<-flank-Window_Size: continue
                            else:
                                if Half[-1]=='M':
                                    M_coverage[(Half[0]-(M_New_bp[0]))/Window_Size+1][-1]+=float(Full_Info[key][7])
                                elif Half[-1]=='P':
                                    P_coverage[(Half[0]-(P_New_bp[0]))/Window_Size+1][-1]+=float(Full_Info[key][7])
                return [M_coverage,P_coverage]
            def c_Coverage_Calculate_2e(Af_Info,Chromo,bp_MP,letter_MP,original_bp_list,flank):
                #Eg of RD_Full_Info_of_Reads (a hash list) elements: 'HWI-ST177_136:2:1:7920:85270': [1202, 1302, 1443, 1543, '+', '-']
                hashM={}
                for i in letter_MP[0]:
                    if not i[0] in hashM.keys():
                        hashM[i[0]]=[i[0]]
                        if (letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))>1:
                            hashM[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))[1:]]
                hashP={}
                for i in letter_MP[1]:
                    if not i[0] in hashP.keys():
                        hashP[i[0]]=[i[0]]
                        if (letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))>1:
                            hashP[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))[1:]]    
                hashMPLetterBP={}
                hashMPLetterBP['M']={}
                hashMPLetterBP['P']={}
                for j in range(len(letter_MP[0])):
                    hashMPLetterBP['M'][hashM[letter_MP[0][j][0]][0]]=[bp_MP[0][j],bp_MP[0][j+1]]
                    hashM[letter_MP[0][j][0]].remove(hashM[letter_MP[0][j][0]][0])
                for j in range(len(letter_MP[1])):
                    hashMPLetterBP['P'][hashP[letter_MP[1][j][0]][0]]=[bp_MP[1][j],bp_MP[1][j+1]]
                    hashP[letter_MP[1][j][0]].remove(hashP[letter_MP[1][j][0]][0])
                hashM={}
                hashM['left']=['left']
                hashM['right']=['right']
                for i in letter_MP[0]:
                    if not i[0] in hashM.keys():
                        hashM[i[0]]=[i[0]]
                        if (letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))>1:
                            hashM[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[0].count(i[0])+letter_MP[0].count(i[0]+'^'))[1:]]
                hashP={}
                hashP['left']=['left']
                hashP['right']=['right']
                for i in letter_MP[1]:
                    if not i[0] in hashP.keys():
                        hashP[i[0]]=[i[0]]
                        if (letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))>1:
                            hashP[i[0]]+=[i[0]+'_'+str(j) for j in range(letter_MP[1].count(i[0])+letter_MP[1].count(i[0]+'^'))[1:]]    
                M_Coverage={}
                M_Coverage['left']=0
                for key_1 in hashMPLetterBP['M'].keys():
                    M_Coverage[key_1]=[0 for i in range((hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])/Window_Size)]
                    if ((hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])-(hashMPLetterBP['M'][key_1][1]-hashMPLetterBP['M'][key_1][0])/Window_Size*Window_Size)>30:
                        M_Coverage[key_1].append(0)
                P_Coverage={}
                P_Coverage['left']=0
                for key_1 in hashMPLetterBP['P'].keys():
                    P_Coverage[key_1]=[0 for i in range((hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])/Window_Size)]
                    if ((hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])-(hashMPLetterBP['P'][key_1][1]-hashMPLetterBP['P'][key_1][0])/Window_Size*Window_Size)>30:
                        P_Coverage[key_1].append(0)
                for key in Af_Info.keys():
                    if Af_Info[key][0]==Af_Info[key][1]==Af_Info[key][2]==Af_Info[key][3]==(-flank/2):
                        M_Coverage['left']+=0.5
                        P_Coverage['left']+=0.5
                    else:
                        if key in Letter_Through.keys():
                            if Af_Info[key][6]=='M':
                                lele=hashM[Letter_Through[key][6]]
                                rile=hashM[Letter_Through[key][9]]
                                lebl=Af_Info[key][:2]
                                ribl=Af_Info[key][2:4]
                                for lele1 in lele:
                                    if lele1=='left' or lele1=='right': continue
                                    block=[lele2-bps[0] for lele2 in hashMPLetterBP['M'][lele1]]
                                    if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
                                        lebl=[k-block[0] for k in lebl]
                                        M_Coverage[lele1][lebl[0]/Window_Size]+=float(lebl[0]/Window_Size*Window_Size+Window_Size-lebl[0])/float(lebl[1]-lebl[0])
                                        if lebl[1]/Window_Size<len(M_Coverage[lele1]):
                                            M_Coverage[lele1][lebl[1]/Window_Size]+=float(lebl[1]-lebl[1]/Window_Size*Window_Size)/float(lebl[1]-lebl[0])
                                for rile1 in rile:
                                    if rile1=='left' or rile1=='right':continue
                                    block=[rile2-bps[0] for rile2 in hashMPLetterBP['M'][rile1]]
                                    if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
                                        ribl=[k-block[0] for k in ribl]
                                        M_Coverage[rile1][ribl[0]/Window_Size]+=float(ribl[0]/Window_Size*Window_Size+Window_Size-ribl[0])/float(ribl[1]-ribl[0])
                                        if ribl[1]/Window_Size<len(M_Coverage[rile1]):
                                            M_Coverage[rile1][ribl[1]/Window_Size]+=float(ribl[1]-ribl[1]/Window_Size*Window_Size)/float(ribl[1]-ribl[0])
                            if Af_Info[key][6]=='P':
                                lele=hashP[Letter_Through[key][6]]
                                rile=hashP[Letter_Through[key][9]]
                                lebl=Af_Info[key][:2]
                                ribl=Af_Info[key][2:4]
                                for lele1 in lele:
                                    if lele1=='left' or lele1=='right': continue
                                    block=[lele2-bps[0] for lele2 in hashMPLetterBP['P'][lele1]]
                                    if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
                                        lebl=[k-block[0] for k in lebl]
                                        P_Coverage[lele1][lebl[0]/Window_Size]+=float(lebl[0]/Window_Size*Window_Size+Window_Size-lebl[0])/float(lebl[1]-lebl[0])
                                        if lebl[1]/Window_Size<len(P_Coverage[lele1]):
                                            P_Coverage[lele1][lebl[1]/Window_Size]+=float(lebl[1]-lebl[1]/Window_Size*Window_Size)/float(lebl[1]-lebl[0])
                                for rile1 in rile:
                                    if rile1=='left' or rile1=='right':continue
                                    block=[rile2-bps[0] for rile2 in hashMPLetterBP['P'][rile1]]
                                    if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
                                        ribl=[k-block[0] for k in ribl]
                                        P_Coverage[rile1][ribl[0]/Window_Size]+=float(ribl[0]/Window_Size*Window_Size+Window_Size-ribl[0])/float(ribl[1]-ribl[0])
                                        if ribl[1]/Window_Size<len(P_Coverage[rile1]):
                                            P_Coverage[rile1][ribl[1]/Window_Size]+=float(ribl[1]-ribl[1]/Window_Size*Window_Size)/float(ribl[1]-ribl[0])
                        if not key in Letter_Through.keys():    
                                key2='_'.join(key.split('_')[:-1])
                                if Af_Info[key][6]=='M':
                                    lele=hashM[Letter_Through[key2][6]]
                                    rile=hashM[Letter_Through[key2][9]]
                                    lebl=Af_Info[key][:2]
                                    ribl=Af_Info[key][2:4]
                                    for lele1 in lele:
                                        if lele1=='left' or lele1=='right': continue
                                        block=[lele2-bps[0] for lele2 in hashMPLetterBP['M'][lele1]]
                                        if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
                                            lebl=[k-block[0] for k in lebl]
                                            M_Coverage[lele1][lebl[0]/Window_Size]+=float(lebl[0]/Window_Size*Window_Size+Window_Size-lebl[0])/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
                                            if lebl[1]/Window_Size<len(M_Coverage[lele1]):
                                                M_Coverage[lele1][lebl[1]/Window_Size]+=float(lebl[1]-lebl[1]/Window_Size*Window_Size)/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
                                    for rile1 in rile:
                                        if rile1=='left' or rile1=='right':continue
                                        block=[rile2-bps[0] for rile2 in hashMPLetterBP['M'][rile1]]
                                        if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
                                            ribl=[k-block[0] for k in ribl]
                                            M_Coverage[rile1][ribl[0]/Window_Size]+=float(ribl[0]/Window_Size*Window_Size+Window_Size-ribl[0])/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
                                            if ribl[1]/Window_Size<len(M_Coverage[rile1]):
                                                M_Coverage[rile1][ribl[1]/Window_Size]+=float(ribl[1]-ribl[1]/Window_Size*Window_Size)/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
                                if Af_Info[key][6]=='P':
                                    lele=hashP[Letter_Through[key2][6]]
                                    rile=hashP[Letter_Through[key2][9]]
                                    lebl=Af_Info[key][:2]
                                    ribl=Af_Info[key][2:4]
                                    for lele1 in lele:
                                        if lele1=='left' or lele1=='right': continue
                                        block=[lele2-bps[0] for lele2 in hashMPLetterBP['P'][lele1]]
                                        if numpy.min(lebl)+15>block[0] and numpy.max(lebl)-15<block[1]:
                                            lebl=[k-block[0] for k in lebl]
                                            P_Coverage[lele1][lebl[0]/Window_Size]+=float(lebl[0]/Window_Size*Window_Size+Window_Size-lebl[0])/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
                                            if lebl[1]/Window_Size<len(P_Coverage[lele1]):
                                                P_Coverage[lele1][lebl[1]/Window_Size]+=float(lebl[1]-lebl[1]/Window_Size*Window_Size)/float(lebl[1]-lebl[0])*float(Af_Info[key][7])
                                    for rile1 in rile:
                                        if rile1=='left' or rile1=='right':continue
                                        block=[rile2-bps[0] for rile2 in hashMPLetterBP['P'][rile1]]
                                        if numpy.min(ribl)+15>block[0] and numpy.max(ribl)-15<block[1]:
                                            ribl=[k-block[0] for k in ribl]
                                            P_Coverage[rile1][ribl[0]/Window_Size]+=float(ribl[0]/Window_Size*Window_Size+Window_Size-ribl[0])/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
                                            if ribl[1]/Window_Size<len(P_Coverage[rile1]):
                                                P_Coverage[rile1][ribl[1]/Window_Size]+=float(ribl[1]-ribl[1]/Window_Size*Window_Size)/float(ribl[1]-ribl[0])*float(Af_Info[key][7])
                return [M_Coverage,P_Coverage]
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
                            Qual_Filter_1[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_1]
                            return Qual_Filter_1
                        else:
                            Qual_Filter_2=[]
                            for j2 in Qual_Filter_1:
                                if j2[-2:]==['+','-']:
                                      Qual_Filter_2.append(j2)
                            if not Qual_Filter_2==[]:
                                if len(Qual_Filter_2)==1:
                                    Qual_Filter_2[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                    return Qual_Filter_2
                                else:
                                    Qual_Filter_3=[]
                                    Qual_IL=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                    for jq in range(len(Qual_IL)):
                                        if Qual_IL[jq]==max(Qual_IL) and not Qual_Filter_1[jq] in Qual_Filter_3:
                                            Qual_Filter_3.append(Qual_Filter_1[jq]+[max(Qual_IL)])
                                    return Qual_Filter_3                    
                            else:
                                Qual_Filter_2=Qual_Filter_1
                                if len(Qual_Filter_2)==1:
                                    Qual_Filter_2[0]+=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                                    return Qual_Filter_2
                                else:
                                    Qual_Filter_3=[]
                                    Qual_IL=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
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
                    IL_Qual=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
                    for j in range(len(IL_Qual)):
                        if IL_Qual[j]==max(IL_Qual) and not Qual_Filter_2[j] in Qual_Filter_3:
                            Qual_Filter_3.append(Qual_Filter_2[j])
                else:
                    Qual_Filter_2=Qual_Filter_1
                    Qual_Filter_3=[]
                    IL_Qual=[pdf_calculate(max(j3[:4])-min(j3[:4]),GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero) for j3 in Qual_Filter_2]
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
            def Define_Default_SVPredict():
                global tolerance_bp
                tolerance_bp=10
                global min_resolution
                min_resolution=70
                global Best_IL_Score    
                Best_IL_Score=0
                global Best_RD_Score
                Best_RD_Score=0
                global deterministic_flag
                deterministic_flag=0
                if '--deterministic-flag' in dict_opts.keys():
                    deterministic_flag=int(dict_opts['--deterministic-flag'])
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
                global Local_Minumum_Number
                Local_Minumum_Number=100
                global IL_Weight
                global DR_Weight
                global TB_Weight
                IL_Weight=1
                DR_Weight=5
                TB_Weight=5
            def Full_Info_of_Reads_Product(Initial_Bam,bps,total_bps,total_letters,bamChr,flank,QCAlign,ReadLength,chr_link):
                #   letters=[chr(97+i) for i in range(len(bps)-1)]
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
                fbam=os.popen(r'''samtools view -F 256 %s %s:%d-%d'''%(Initial_Bam,bamChr,bps[0]-flank,bps[-1]+flank))
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
                        block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1)
                        block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2)
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
                                    block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1)
                                    block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2)
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
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double[key][0][1])
                        elif pos1>pos2:
                            pos1=int(Letter_Double[key][1][3])
                            pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                            pos2=int(Letter_Double[key][0][3])
                            pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double[key][1][1])
                        block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1+low_qual_edge)
                        block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2+low_qual_edge)
                        block1b=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1b-low_qual_edge)
                        block2b=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2b-low_qual_edge)
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
                        if Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][7]))==0:
                            if Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
                                BlockCov[Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
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
                    Initial_ILPenal+=[pdf_calculate(j,GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero)/len(Initial_IL)]
                return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov,Letter_Double]
            def Full_Info_of_Reads_Product_3(Initial_Bam,temp_bp,temp_let,bamChr,target_region,Chr_Link):
                Letter_Double={}
                Pair_ThroughBP=[]
                Double_Read_ThroughBP=[]
                Single_Read_ThroughBP=[]
                blackList=[]
                fbam=os.popen(r'''samtools view -F 256 %s %s:%d-%d'''%(Initial_Bam,bamChr,target_region[0]-flank,target_region[-1]+flank))
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
                        block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1)
                        block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2)
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
                                    block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1)
                                    block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2)
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
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double[key][0][1])
                        elif pos1>pos2:
                            pos1=int(Letter_Double[key][1][3])
                            pos1b=pos2+cigar2reaadlength(Letter_Double[key][1][5])
                            pos2=int(Letter_Double[key][0][3])
                            pos2b=pos1+cigar2reaadlength(Letter_Double[key][0][5])
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double[key][1][1])
                        block1=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1+low_qual_edge)
                        block2=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2+low_qual_edge)
                        block1b=Reads_block_assignment_1(flank,temp_bp,temp_let,pos1b-low_qual_edge)
                        block2b=Reads_block_assignment_1(flank,temp_bp,temp_let,pos2b-low_qual_edge)
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
                        if Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][3]))==Reads_block_assignment_1(flank,temp_bp,temp_let,int(Letter_Double[key][0][3])+cigar2reaadlength(Letter_Double[key][0][5])):
                            BlockCov[Reads_block_assignment_1(flank,flank,temp_bp,temp_let,int(Letter_Double[key][0][3]))]+=cigar2reaadlength(Letter_Double[key][0][5])
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
            def Full_Info_of_Reads_Integrate(GC_para_dict,BP_para_dict,bps2):
                bps2_left=[]
                bps2_right=[]
                for x in bps2:
                    bps2_left.append([x[0],x[1]-5000,x[1]])
                    bps2_right.append([x[0],x[-1],x[-1]+5000])
                chr_letter_bp=letter_rearrange(BP_para_dict['flank'],bps2)
                letter_GC=letter_GC_ReadIn(chr_letter_bp)
                letter_RD_test=letter_RD_ReadIn(letter_RD_test_calcu(chr_letter_bp))
                if len(bps2)==1 and len(bps2[0])==3 and letter_RD_test[bps2[0][0]]['a']>GC_para_dict['GC_Overall_Median_Coverage'][bps2[0][0]]*4:
                    return [letter_RD_test[bps2[0][0]],letter_RD_test[bps2[0][0]],0,0,[],[],[],letter_GC[bps2[0][0]]]+original_bp_let_produce(chr_letter_bp,bps2)
                letter_RD=letter_RD_ReadIn(chr_letter_bp)
                Multi_Dup=multi_dup_define(letter_RD,GC_para_dict['GC_Overall_Median_Coverage'])
                global letter_RD_left_control
                letter_RD_left_control=letter_RD_ReadIn(letter_rearrange(BP_para_dict['flank'],bps2_left))
                global letter_RD_right_control
                letter_RD_right_control=letter_RD_ReadIn(letter_rearrange(BP_para_dict['flank'],bps2_right))
                letter_range_report(BP_para_dict['flank'],chr_letter_bp)
                blocks_read_in=block_Read_From_Bam(chr_letter_bp)
                read_info=block_Info_ReadIn(GC_para_dict,BP_para_dict,chr_letter_bp,blocks_read_in,Multi_Dup)
                block_rds=read_info[0]
                block_rd2=read_info[1]
                letter_RD2={}
                for k1 in letter_RD.keys():
                    for k2 in letter_RD[k1].keys():
                        if k2 in Multi_Dup:
                            letter_RD2[k2]=letter_RD[k1][k2]
                            if not k1 in block_rd2.keys():
                                block_rd2[k1]={}
                            if not k2 in block_rd2[k1].keys():
                                block_rd2[k1][k2]=0
                        else:
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
                Initial_RD=total_rd_calcu(GC_para_dict['GC_Median_Num'],GC_para_dict['GC_Overall_Median_Num'],letter_RD2,letter_GC,chr_letter_bp,block_rd2)
                DR_Penal=DR_Penal_Calcu(read_info)
                IL_Penal=IL_Penal_Calcu(read_info,GC_para_dict['IL_Statistics'],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero)
                letter_GC_out={}
                for k1 in letter_GC.keys():
                        for k2 in letter_GC[k1].keys():
                                letter_GC_out[k2]=letter_GC[k1][k2]
                return [letter_RD2,Initial_RD,DR_Penal,numpy.mean(IL_Penal),Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,letter_GC_out]+original_bp_let_produce(chr_letter_bp,bps2)
            def global_para_declaration():
                global chrom_N
                global chrom_X
                global chrom_Y  
                global workdir
                workdir=path_modify(dict_opts['--workdir'])
                global bp_txt_Path
                global BPPath
                global NullPath
                global ref_path
                global ref_file
                global ref_index
                global ref_ppre
                global ref_prefix
                ref_path=workdir+'reference_SVelter/'
                ref_file=ref_path+'genome.fa'
                ref_index=ref_file+'.fai'
                ref_ppre=ref_path
                ref_prefix='.'.join(ref_file.split('.')[:-1])
                global GC_hash
                GC_hash=GC_Index_Readin(ref_prefix+'.GC_Content')
            def letter_GC_ReadIn(chr_letter_bp):
                block_GC_temp={}
                filein=ref_prefix+'.GC_Content'
                block_range={}
                GC_hash_temp={}
                test_flag=0
                for i in chr_letter_bp.keys():
                    if not os.path.isfile(filein):
                        test_flag+=1
                if test_flag==0:
                    for i in chr_letter_bp.keys():
                        GC_hash_temp[i]={}
                        block_range[i]=[]
                        for j in chr_letter_bp[i].keys():
                            block_range[i]+=chr_letter_bp[i][j]
                        block_range[i]=[min(block_range[i]),max(block_range[i])]
                        for xa in GC_hash[i].keys():
                            for xb in GC_hash[i][xa].keys():
                                if not xb<block_range[i][0] and not xa>block_range[i][1]:
                                    GC_hash_temp[i][str(xa)+'-'+str(xb)]=GC_hash[i][xa][xb]
                    for k1 in chr_letter_bp.keys():
                        block_GC_temp[k1]={}
                        for k2 in GC_hash_temp[k1].keys():
                            bl2=[int(k2.split('-')[0]),int(k2.split('-')[1])]
                            for k3 in chr_letter_bp[k1].keys():
                                if min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])<bl2[1]+1:
                                    block_GC_temp[k1][k3]=GC_hash_temp[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]
                                elif min(chr_letter_bp[k1][k3])>bl2[0]-1 and max(chr_letter_bp[k1][k3])>bl2[1]:
                                    if not k3 in block_GC_temp[k1].keys():
                                        block_GC_temp[k1][k3]=GC_hash_temp[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:]
                                    else:
                                        block_GC_temp[k1][k3]+=GC_hash_temp[k1][k2][(min(chr_letter_bp[k1][k3])-bl2[0])/100:]                        
                                elif min(chr_letter_bp[k1][k3])<bl2[0] and max(chr_letter_bp[k1][k3])>bl2[0]-1:
                                    if not k3 in block_GC_temp[k1].keys():
                                        block_GC_temp[k1][k3]=GC_hash_temp[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]
                                    else:
                                        block_GC_temp[k1][k3]+=GC_hash_temp[k1][k2][:(max(chr_letter_bp[k1][k3])-bl2[0])/100+1]                      
                                elif min(chr_letter_bp[k1][k3])<bl2[0]+1 and max(chr_letter_bp[k1][k3])>bl2[1]-1:
                                    if not k3 in block_GC_temp[k1].keys():
                                        block_GC_temp[k1][k3]=GC_hash_temp[k1][k2]
                                    else:
                                        block_GC_temp[k1][k3]+=GC_hash_temp[k1][k2]                     
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
                    filein=NullPath+'RD_Stat/'+BamN+'.'+k1+'.RD.index'
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
                        filein=NullPath+'RD_Stat/'+BamN+'.'+k1+'.RD.index'
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
            def letter_rearrange(flank,bps2):
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
            def letter_RD_test_calcu(chr_letter_bp):
                out={}
                for x in chr_letter_bp.keys():
                    out[x]={}
                    for y in chr_letter_bp[x].keys():
                        if not y in ['left','right']:
                            if len(chr_letter_bp[x][y])==2:
                                out[x][y]=[chr_letter_bp[x][y][0]-500]+chr_letter_bp[x][y]+[chr_letter_bp[x][y][1]+500]
                            else:
                                out[x][y]=chr_letter_bp[x][y]
                return out
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
            def RD_Index_ReadIn(ppre_Path,BamN, chromo, region):
                if not ppre_Path[-1]=='/':
                    ppre_Path+='/'
                path_in=NullPath+'RD_Stat/'
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
            def read_Pair_Single_Read_ThroughBP(chr_letter_bp,Single_Read_ThroughBP):
                out=[]
                for k1 in Single_Read_ThroughBP.keys():
                    for k2 in Single_Read_ThroughBP[k1]:
                        rela=[k2[2],k2[0]-chr_letter_bp[k1][k2[2]][0],
                              k2[3],k2[1]-chr_letter_bp[k1][k2[3]][0]]
                        out.append(rela)
                return out
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
                        rela=[k2[6],k2[0]-chr_letter_bp[k1][k2[6]][0],
                              k2[7],k2[1]-chr_letter_bp[k1][k2[7]][0],
                              k2[8],k2[2]-chr_letter_bp[k1][k2[8]][0],
                              k2[9],k2[3]-chr_letter_bp[k1][k2[9]][0],k2[4],k2[5]]
                        out.append(rela)
                return out
            def Single_Rec_Read_Locate(BP_para_dict,Letter_Double_rec,temp_bp, temp_let):
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
                        fbamtemp=os.popen(r'''samtools view -F 256 %s %s:%d-%d'''%(Initial_Bam,bamChr,pos2,pos2+ReadLength))
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
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double_rec[key][0][1])
                        elif pos1>pos2:
                            pos1=int(Letter_Double_rec[key][1][3])
                            pos1b=pos2+cigar2reaadlength(Letter_Double_rec[key][1][5])
                            pos2=int(Letter_Double_rec[key][0][3])
                            pos2b=pos1+cigar2reaadlength(Letter_Double_rec[key][0][5])
                            direct_temp=Reads_Direction_Detect_flag(Letter_Double_rec[key][1][1])
                        if not pos1<temp_bp[0]-BP_para_dict['flank']+1 and not pos2b>temp_bp[-1]+BP_para_dict['flank']-1:
                            block1=Reads_block_assignment_1(BP_para_dict['flank'],temp_bp,temp_let,pos1+low_qual_edge)
                            block2=Reads_block_assignment_1(BP_para_dict['flank'],temp_bp,temp_let,pos2+low_qual_edge)
                            block1b=Reads_block_assignment_1(BP_para_dict['flank'],temp_bp,temp_let,pos1b-low_qual_edge)
                            block2b=Reads_block_assignment_1(BP_para_dict['flank'],temp_bp,temp_let,pos2b-low_qual_edge)
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
                    Initial_ILPenal+=[pdf_calculate(j,GC_para_dict['IL_Statistics'][4],GC_para_dict['IL_Statistics'][0],GC_para_dict['IL_Statistics'][1],GC_para_dict['IL_Statistics'][2],GC_para_dict['IL_Statistics'][3],BP_para_dict['Cut_Upper'],BP_para_dict['Cut_Lower'],Penalty_For_InsertLengthZero)/len(Initial_IL)]
                return [Initial_DR_Penal,Initial_ILPenal,Pair_ThroughBP,Double_Read_ThroughBP,Single_Read_ThroughBP,BlockCov,Initial_Cov]
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
            def Insert_Seq_Pool_Prod_2(original_bp_list,ori_1_Seq,flank):
                ini_letters=['left']+['I'+chr(97+i) for i in range(len(original_bp_list)-1)]+['right']+['I'+chr(97+i)+'^' for i in range(len(original_bp_list)-1)]
                relative_bps=[0]+[j-original_bp_list[0]+flank for j in original_bp_list]+[original_bp_list[-1]+flank-original_bp_list[0]+flank]
                Insert_Seq_Pool={}
                for k in range(len(original_bp_list)+1):
                    Insert_Seq_Pool[ini_letters[k]]=ori_1_Seq[relative_bps[k]:relative_bps[k+1]]
                for k in range(len(original_bp_list)+1,len(ini_letters)):
                    Insert_Seq_Pool[ini_letters[k]]=complementary(ori_1_Seq[relative_bps[k-len(original_bp_list)]:relative_bps[k+1-len(original_bp_list)]])
                return Insert_Seq_Pool
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
            def penal_calculate(GC_para_dict,BP_para_dict,Map_All,temp_bp, Af_Letter,Af_BP,letters_numbers,NoMapPenal):
                out_rd=[[0 for i in temp_bp[0][:-1]],[0 for i in temp_bp[1][:-1]]]
                IL_Rec={}
                DR_Penal=0
                out_tb=[[0 for i in temp_bp[0]],[0 for i in temp_bp[1]]]
                for i in Map_All:
                    if len(i)>4:
                        if not i[6] in IL_Rec.keys():
                            IL_Rec[i[6]]=i[8]
                        else:
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
                #for i in Num_Read_TB:
                #   for j in i:
                #       if j<Through_BP_Minimum:
                #           TB_Pena_2_out+=1
                Af_Block_Len=[[BP_para_dict['flank']]+[Af_BP[0][i+1]-Af_BP[0][i] for i in range(len(Af_BP[0])-1)]+[BP_para_dict['flank']],[BP_para_dict['flank']]+[Af_BP[1][i+1]-Af_BP[1][i] for i in range(len(Af_BP[1])-1)]+[BP_para_dict['flank']]]
                out_rd=[[out_rd[0][i]/Af_Block_Len[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]/Af_Block_Len[1][i] for i in range(len(out_rd[1]))]]
                out_rd_new=[[(BP_para_dict['RD_within_B']['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[0][0],
                (BP_para_dict['RD_within_B']['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[0][-1]],
                [(BP_para_dict['RD_within_B']['left']-out_rd[0][0]-out_rd[1][0])/2.0+out_rd[1][0],
                (BP_para_dict['RD_within_B']['right']-out_rd[0][-1]-out_rd[1][-1])/2.0+out_rd[1][-1]]]
                out_rd=[[out_rd_new[0][0]]+out_rd[0][1:-1]+[out_rd_new[0][-1]],[out_rd_new[1][0]]+out_rd[1][1:-1]+[out_rd_new[1][-1]]]
                out_rd_within=[[BP_para_dict['RD_within_B'][Af_Letter[0][i]]/letters_numbers[0][i] for i in range(len(Af_Letter[0]))],[BP_para_dict['RD_within_B'][Af_Letter[1][i]]/letters_numbers[1][i] for i in range(len(Af_Letter[1]))]]
                out_rd_within[0]=[0]+out_rd_within[0]+[0]
                out_rd_within[1]=[0]+out_rd_within[1]+[0]   
                cov_bp2=[[out_rd[0][i]+out_rd_within[0][i] for i in range(len(out_rd[0]))],[out_rd[1][i]+out_rd_within[1][i] for i in range(len(out_rd[1]))]]
                Cov_GC=[[BP_para_dict['BlockGC2'][k] for k in Af_Letter[0]],[BP_para_dict['BlockGC2'][k] for k in Af_Letter[1]]]
                adj_cov_bp=[GC_RD_Adj(GC_para_dict['GC_Median_Num'],GC_para_dict['GC_Overall_Median_Num'],chrom_N,Cov_GC[0],cov_bp2[0][1:-1]),GC_RD_Adj(GC_para_dict['GC_Median_Num'],GC_para_dict['GC_Overall_Median_Num'],chrom_N,Cov_GC[1],cov_bp2[1][1:-1])]
                return [IL_Output,adj_cov_bp,DR_Penal,TB_Pena_2_out,Num_total_TB]
            def Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP):
                Total_Cov_For_Pen={}
                for key in BP_para_dict['RD_within_B'].keys():
                    Total_Cov_For_Pen[key]=0
                Map_M=[]
                Map_P=[]
                Map_Both=[]
                Let_BP_Info={}
                Let_BP_Info['m']={}
                Let_BP_Info['p']={}
                temp_letter=[['left']+Af_Letter[0]+['right'],['left']+Af_Letter[1]+['right']]
                temp_bp=[[Af_BP[0][0]-BP_para_dict['flank']]+Af_BP[0]+[Af_BP[0][-1]+BP_para_dict['flank']],[Af_BP[1][0]-BP_para_dict['flank']]+Af_BP[1]+[Af_BP[1][-1]+BP_para_dict['flank']]]
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
                NoMapPenal=Be_Info_3_rearrange(BP_para_dict,Be_Info,temp_letter,Let_BP_Info,Total_Cov_For_Pen,Map_M,Map_P,Map_Both,NoMapPenal)
                best_structure_sign_flag=0
                for key in Total_Cov_For_Pen.keys():
                    if Total_Cov_For_Pen[key]==0:
                        del Total_Cov_For_Pen[key]
                    else:
                        Total_Cov_For_Pen[key]/=float(Be_BP_Letter[key])
                for key in BP_para_dict['RD_within_B'].keys():
                    if not key[-1]=='^' and not key in ['left','right','left^', 'right^']:
                        if not key in Af_Letter[0]+Af_Letter[1] and not key+'^' in Af_Letter[0]+Af_Letter[1]:
                            if not key in Total_Cov_For_Pen.keys():
                                Total_Cov_For_Pen[key]=0
                            Total_Cov_For_Pen[key]+=BP_para_dict['RD_within_B'][key]
                if NoMapPenal>0:
                    best_structure_sign_flag+=1
                for key1 in Total_Cov_For_Pen.keys():
                    if Total_Cov_For_Pen[key1]>2.58*GC_para_dict['GC_Std_Coverage'][chrom_N]:
                        best_structure_sign_flag+=1
                if not Map_M+Map_P+Map_Both==[]:
                    penals=penal_calculate(GC_para_dict,BP_para_dict,Map_M+Map_P+Map_Both,temp_bp,Af_Letter,Af_BP,letters_numbers,NoMapPenal)
                    if penals[2]>0:
                        best_structure_sign_flag+=1
                    return penals[:-1]+[NoMapPenal,Total_Cov_For_Pen,best_structure_sign_flag]+[penals[-1]]
                else:
                    return 0
            def write_best_letter(bps_all,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash,original_letters):
                fo=open(output_Score_File,'a')
                time2=time.time()
                Best_Letter_2=[]
                if not Score_rec_hash=={}:
                    temp1=Best_Let_modify(original_letters,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash)
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
                        if Uniparental_disomy_check(original_letters,bestletter)=='Pass':
                            print >>fo, ' '.join([str(bp_ele) for bp_ele in bps3])
                            print >>fo, '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
                            #print  ' '.join([str(bp_ele) for bp_ele in bps3])
                            #print  '/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
                            print >>fo, 'Theoretical Best Score: '+str(Best_IL_Score+Best_RD_Score+20)
                            if Best_Score_Rec>80:
                                Best_Score_Rec=80
                            print >>fo, 'Current Best Scure: '+str(Best_Score_Rec+20)
                            print>>fo, 'Time Consuming:'+str(datetime.timedelta(seconds=(time2-time1)))
                fo.close()
            def score_rec_hash_Modify_for_short_del(Score_rec_hash):
                Score_rec_hash_new={}
                for x in sorted(Score_rec_hash.keys())[::-1][:1]:
                    Score_rec_hash_new[x]=Score_rec_hash[x]
                for x in sorted(Score_rec_hash.keys())[::-1][1:]:
                    Score_rec_hash_new[x-1.1]=Score_rec_hash[x]
                return Score_rec_hash_new
            def one_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash):
                #Letter_Candidates=[[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']],[['a','a^'], []],[['a^','a'], []],[['a^','a^'], []]]
                Letter_Candidates=[[[],[]],[['a'], []],[['a^'], []],[['a'], ['a']],[['a^'], ['a']],[['a^'], ['a^']]]
                if Ploidy==2:
                    Letter_Candidates=Letter_Candidates
                elif Ploidy==1:
                    Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
                elif Ploidy==0:
                    Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
                IL_RD_Temp_Info=Af_Rearrange_Info_Collect(GC_para_dict,BP_para_dict,Be_BP_Letter,Be_Info,Letter_Candidates)
                if not IL_RD_Temp_Info=='Error':
                    ILTemp=IL_RD_Temp_Info[0]
                    RDTemp=IL_RD_Temp_Info[1]
                    Letter_Rec=IL_RD_Temp_Info[2]
                    BP_Rec=IL_RD_Temp_Info[3]
                    if not ILTemp==[]:
                        DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_para_dict['GC_Var_Coverage'])
                        Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                        Best_Score_Rec=ILTemp[DECISION_Score[0]]+RDTemp[DECISION_Score[0]]
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
                else:
                    return 'Error'
            def two_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash):
                Letter_Candidates=struc_propose_single_block(2)+struc_propose_single_block(3)+struc_propose_single_block(4)+struc_propose_single_block(5)   
                Letter_Candidates=[i for i in Letter_Candidates if not [] in i]
                if [[], ['a', 'a']] in Letter_Candidates:
                    del Letter_Candidates[Letter_Candidates.index([[], ['a', 'a']])]
                if [[], ['a^', 'a^']] in Letter_Candidates:
                    del Letter_Candidates[Letter_Candidates.index([[], ['a^', 'a^']])]
                if Ploidy==2:
                    Letter_Candidates=Letter_Candidates
                elif Ploidy==1:
                    Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
                elif Ploidy==0:
                    Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
                IL_RD_Temp_Info=Af_Rearrange_Info_Collect_2(BP_para_dict,Letter_Candidates)
                if not IL_RD_Temp_Info=='Error':
                    ILTemp=IL_RD_Temp_Info[0]
                    RDTemp=IL_RD_Temp_Info[1]
                    Letter_Rec=IL_RD_Temp_Info[2]
                    BP_Rec=IL_RD_Temp_Info[3]
                    if not ILTemp==[]:
                        DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_para_dict['GC_Var_Coverage'])
                        Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                        Best_Score_Rec=ILTemp[DECISION_Score[0]]+RDTemp[DECISION_Score[0]]
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
                else: return 'Error'
            def few_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash):
                Letter_Candidates=struc_propose_single_block(copy_num_a)+struc_propose_single_block(copy_num_b)
                Letter_Candidates=[i for i in Letter_Candidates if not [] in i]
                if Ploidy==2:
                    Letter_Candidates=Letter_Candidates
                elif Ploidy==1:
                    Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
                elif Ploidy==0:
                    Letter_Candidates=[i for i in Letter_Candidates if ['a'] in i]
                IL_RD_Temp_Info=Af_Rearrange_Info_Collect_2(BP_para_dict,Letter_Candidates)
                if not IL_RD_Temp_Info=='Error':
                    ILTemp=IL_RD_Temp_Info[0]
                    RDTemp=IL_RD_Temp_Info[1]
                    Letter_Rec=IL_RD_Temp_Info[2]
                    BP_Rec=IL_RD_Temp_Info[3]
                    if not ILTemp==[]:
                        DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_para_dict['GC_Var_Coverage'])
                        Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                        Best_Score_Rec=ILTemp[DECISION_Score[0]]+RDTemp[DECISION_Score[0]]
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
                else: return 'Error'
            def two_block_RD_Process(GC_para_dict,BP_para_dict,run_flag):
                Letter_Candidates=struc_produce_two_block(Copy_num_estimate)
                if Ploidy==2:
                    Letter_Candidates=Letter_Candidates
                elif Ploidy==1:
                    Letter_Candidates=[i for i in Letter_Candidates if i[0]==i[1]]
                elif Ploidy==0:
                    Letter_Candidates=[i for i in Letter_Candidates if ['a','b'] in i]
                IL_RD_Temp_Info=Af_Rearrange_Info_Collect_2(BP_para_dict,Letter_Candidates)
                if not IL_RD_Temp_Info=='Error':
                    ILTemp=IL_RD_Temp_Info[0]
                    RDTemp=IL_RD_Temp_Info[1]
                    Letter_Rec=IL_RD_Temp_Info[2]
                    BP_Rec=IL_RD_Temp_Info[3]
                    if not ILTemp==[]:
                        DECISION_Score=Move_Decide_3(ILTemp,RDTemp,GC_para_dict['GC_Var_Coverage'])
                        Best_Letter_Rec=[Letter_Rec[DECISION_Score[0]]]
                        Best_Score_Rec=ILTemp[DECISION_Score[0]]+RDTemp[DECISION_Score[0]]
                        run_flag+=1
                    else:
                        Best_Letter_Rec=[]
                        Best_Score_Rec=100
                        run_flag+=1
                    return([Best_Letter_Rec,Best_Score_Rec,run_flag])
                else: return 'Error'
            Define_Default_SVPredict()
            if not '--workdir' in dict_opts.keys():
                print 'Error: please specify working directory using: --workdir'
            else:
                global_para_declaration()
                if not '--bp-file' in dict_opts.keys():
                    print 'Error: please specify input txt file using : --bp-file'
                else:
                    if not '--out-path' in dict_opts.keys():
                        dict_opts['--out-path']='/'.join(dict_opts['--bp-file'].split('/')[:-1])
                    if not dict_opts['--out-path'][-1]=='/':
                        dict_opts['--out-path']+='/'
                    if not os.path.isfile(ref_file):
                        print 'Error: wrong reference genome provided'
                    else:
                        if not os.path.isfile(ref_index):
                            print 'Error: reference genome not indexed'
                        else:
                            global chromos_all
                            chromos_all=chromos_readin_list(ref_file)
                            if not '--sample' in dict_opts.keys():
                                print 'Error: please specify either input file using --sample'
                            else:
                                time1=time.time()
                                bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                                BamN=dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'')
                                Input_File=dict_opts['--bp-file']
                                bp_txt_Path='/'.join(Input_File.split('/')[:-1])+'/'
                                BPPath='/'.join(bp_txt_Path.split('/')[:-2])+'/'+'.'.join(['BreakPoints']+bp_txt_Path.split('/')[-2].split('.')[1:])+'/'
                                NullPath='/'.join(bp_txt_Path.split('/')[:-2])+'/'+'.'.join(['NullModel']+bp_txt_Path.split('/')[-2].split('.')[1:])+'/'
                                if not bp_txt_Path[0]=='/':
                                    if BPPath[0]=='/':
                                        BPPath=BPPath[1:]
                                    if NullPath[0]=='/':
                                        NullPath=NullPath[1:]
                                Insert_Len_Stat=NullPath+'ILNull.'+BamN+'.'+genome_name+'.Bimodal'
                                RD_NB_Stat=NullPath+'RDNull.'+BamN+'.'+genome_name+'.NegativeBinomial'
                                global RD_Weight
                                RD_Weight=Insert_len_stat_readin(Insert_Len_Stat)/RD_NB_stat_readin(RD_NB_Stat)
                                if not os.path.isfile(Insert_Len_Stat):
                                    print 'Error: cannot access file: '+Insert_Len_Stat
                                else:
                                    ReadLenFin=NullPath+BamN+'.'+genome_name+'.Stats'
                                    if not os.path.isfile(ReadLenFin):
                                        print 'Error: cannot access file: '+ReadLenFin
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
                                    Initial_Bam_Name=BamN+'.'+bam_files_appdix
                                    Initial_Bam=dict_opts['--sample']
                                    flank=cdf_solver_application(Insert_Len_Stat,0.95,model_comp)
                                    Cut_Lower=cdf_solver_application(Insert_Len_Stat,0.005,model_comp)
                                    Cut_Upper=cdf_solver_application(Insert_Len_Stat,0.995,model_comp)
                                    IL_Stat_all=IL_Stat(Insert_Len_Stat)
                                    IL_Statistics=IL_Stat_all[0]
                                    IL_Normal_Stat=IL_Stat_all[1]
                                    fi_test=os.popen(r'''wc -l %s'''%(Input_File))
                                    line_test=fi_test.readline().strip().split()
                                    fi_test.close()
                                    if not line_test[0]=='0':
                                        IL_Estimate=IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5]
                                        IL_SD=((IL_Statistics[2]*IL_Statistics[4])**2+(IL_Statistics[3]*IL_Statistics[5])**2)**(0.5)
                                        IL_Penal_Two_End_Limit=min([pdf_calculate(IL_Estimate-3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero),pdf_calculate(IL_Estimate+3*IL_SD,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)])
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
                                                        print pi
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
                                        output_Score_File=dict_opts['--out-path']+'_'.join(dict_opts['--bp-file'].split('/')[-1].split('.')[:-1])+'.coverge'
                                        file_setup(output_Score_File)
                                        for bpsk1 in sorted(bps_hash.keys()):
                                            for bps2 in bps_hash[bpsk1]:
                                                for i in bps2:
                                                    if len(i)<3:
                                                        i.append(str(int(i[-1])+Window_Size))
                                        GC_Stat_Path=NullPath+'RD_Stat'
                                        Affix_GC_Stat='_MP'+str(QCAlign)+'_GC_Coverage_ReadLength'
                                        [GC_Content_Coverage,Chromosome,Coverage_0]=GC_Stat_ReadIn(BamN,GC_Stat_Path,genome_name,Affix_GC_Stat)
                                        Coverage=[int(k) for k in Coverage_0]
                                        GC_Overall_Median_Coverage={}
                                        GC_Overall_Median_Num=[]
                                        GC_Median_Coverage={}
                                        GC_Median_Num={}
                                        GC_Mean_Coverage={}
                                        GC_Std_Coverage={}
                                        GC_Var_Coverage={}
                                        for a in Chromosome:
                                            if a in GC_Content_Coverage.keys():
                                                GC_Overall_temp=[]
                                                for b in Coverage:
                                                    if not b in GC_Content_Coverage[a].keys(): continue
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
                                        GC_Overall_Median_Num=Median_Pick([i for i in GC_Overall_Median_Num if not i==0])
                                        for a in GC_Median_Num.keys():
                                            if GC_Median_Num[a]==[]:
                                                GC_Median_Num[a]=GC_Overall_Median_Num
                                            else:
                                                GC_Median_Num[a]=Median_Pick(GC_Median_Num[a])
                                        GC_Median_Num=GC_Median_Num_Correct(GC_Median_Num)
                                        ChrN_Median_Coverage={}
                                        for i in GC_Median_Coverage.keys():
                                            for j in GC_Median_Coverage[i].keys():
                                                if not j in ChrN_Median_Coverage.keys():
                                                    ChrN_Median_Coverage[j]=[GC_Median_Coverage[i][j]]
                                                else:
                                                    ChrN_Median_Coverage[j]+=[GC_Median_Coverage[i][j]]
                                        [chrom_N,chrom_X,chrom_Y,GC_Median_Coverage,GC_Overall_Median_Coverage,GC_Var_Coverage,GC_Mean_Coverage,GC_Std_Coverage]=GC_RD_Info_Complete(ref_file,GC_Median_Coverage,ChrN_Median_Coverage,GC_Overall_Median_Coverage,GC_Var_Coverage,GC_Mean_Coverage,GC_Std_Coverage,Chromosome)
                                        GC_para_dict={'IL_Statistics':IL_Statistics,'GC_Overall_Median_Coverage':GC_Overall_Median_Coverage,'GC_Overall_Median_Num':GC_Overall_Median_Num,'GC_Median_Coverage':GC_Median_Coverage,'GC_Median_Num':GC_Median_Num,'GC_Mean_Coverage':GC_Mean_Coverage,'GC_Std_Coverage':GC_Std_Coverage,'GC_Var_Coverage':GC_Var_Coverage,'Coverage':Coverage}
                                        for bpsk1 in sorted(bps_hash.keys()):
                                            if bpsk1>50: continue
                                            for bps2_new in bps_hash[bpsk1]:
                                                bps2_new_2=modify_bps2_new(bps2_new)
                                                bps2=LN_bps2_Modify(bps2_new_2,chromos_all)
                                                if len(bps2)>0 and qual_check_bps2(bps2)=='right':
                                                    Chromo=bps2[0][0]
                                                    if not str(Chromo) in GC_Std_Coverage.keys(): continue
                                                    if not str(Chromo) in GC_Mean_Coverage.keys(): continue
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
                                                    if len(bps2)<1: continue
                                                    original_bps_all=[]
                                                    for obas in bps2:
                                                        original_bps_all+=obas
                                                    original_structure=bp_to_let(original_bps_all,chromos_all)
                                                    chr_letter_tbp=letter_rearrange(flank,bps2)
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
                                                                        bps4[k1].append(bps3[k1][sorted(bps3[k1].keys())[k2+1]])
                                                        bps2=bps4_to_bps2(bps4)
                                                        Chr=bps2[0][0]
                                                        Flank_para_dict={'flank':flank,'Cut_Lower':Cut_Lower,'Cut_Upper':Cut_Upper,'ReadLength':ReadLength}
                                                        Copy_num_esti_info=copy_num_estimate_calcu(GC_para_dict,Flank_para_dict,bps2)
                                                        Copy_num_Check=Copy_num_esti_info[1]
                                                        Copy_num_estimate=Copy_num_esti_info[0]
                                                        if Copy_num_Check==[]:
                                                            Full_Info=Full_Info_of_Reads_Integrate(GC_para_dict,Flank_para_dict,bps2)
                                                            RD_within_B=RD_within_B_calcu(GC_Mean_Coverage,Full_Info,bps2)
                                                            for j in range(Cut_Lower,Cut_Upper+1):
                                                                Single_ILScore=pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)
                                                                Best_IL_Score+=Single_ILScore*exp(Single_ILScore)
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
                                                            #if Copy_num_Check==[]:
                                                            median_CN=GC_Overall_Median_Coverage[chrom_N]/2
                                                            for key in Initial_GCRD_Adj.keys():
                                                                if not key in ['left','right']:
                                                                    Block_CN_Upper[key]=Initial_GCRD_Adj[key]/median_CN+2
                                                            Initial_DR=Full_Info[2]
                                                            Initial_IL=Full_Info[3]
                                                            #Initial_Info=Full_Info[4]+Full_Info[5]+Full_Info[6]
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
                                                            #best_score_rec=[]
                                                            num_of_reads=(original_bp_list[-1]-original_bp_list[0])*GC_Mean_Coverage[Chr]/2/ReadLength
                                                            Best_Score_Rec=0
                                                            Score_rec_hash={}
                                                            break_Iteration_Flag=0
                                                            run_flag=0
                                                            Best_Letter_Rec=[]
                                                            BP_para_dict={'flank':flank,'Cut_Lower':Cut_Lower,'Cut_Upper':Cut_Upper,'ReadLength':ReadLength,'Be_Letter':Be_Letter,'num_of_reads':num_of_reads,'original_letters':original_letters,'BlockGC2':BlockGC2,'BlockGC':BlockGC,'original_bp_list':original_bp_list,'RD_within_B':RD_within_B}
                                                            if len(Full_Info[9])==1:
                                                                if Full_Info[1]['a']<GC_Mean_Coverage[Chr]/4 and Full_Info[2]<3:
                                                                    Run_Result=zero_RD_Process(original_bp_list,run_flag,Best_IL_Score,Best_RD_Score)
                                                                    if Run_Result=='Error': continue
                                                                    Best_Letter_Rec=Run_Result[0]
                                                                    Best_Score_Rec=Run_Result[1]
                                                                    run_flag=Run_Result[2]
                                                                    Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                else:
                                                                    if Full_Info[1]['a']<GC_Mean_Coverage[Chr]:
                                                                        Run_Result=one_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash)
                                                                        if Run_Result=='Error': continue
                                                                        Best_Letter_Rec=Run_Result[0]
                                                                        Best_Score_Rec=Run_Result[1]
                                                                        run_flag=Run_Result[2]
                                                                        Score_rec_hash=Run_Result[3]
                                                                        #Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                    else:
                                                                        if Full_Info[1]['a']<2*GC_Mean_Coverage[Chr]:
                                                                            Run_Result=two_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash)
                                                                            if Run_Result=='Error': continue
                                                                            Best_Letter_Rec=Run_Result[0]
                                                                            Best_Score_Rec=Run_Result[1]
                                                                            run_flag=Run_Result[2]
                                                                            Score_rec_hash=Run_Result[3]
                                                                            #Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                        else:
                                                                            copy_num_a=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))
                                                                            copy_num_b=int(float(Full_Info[1]['a'])/(float(GC_Mean_Coverage[Chr])/2))+1
                                                                            if copy_num_b<4:
                                                                                Run_Result=few_RD_Process(GC_para_dict,BP_para_dict,run_flag,Score_rec_hash)
                                                                                if Run_Result=='Error': continue
                                                                                Best_Letter_Rec=Run_Result[0]
                                                                                Best_Score_Rec=Run_Result[1]
                                                                                run_flag=Run_Result[2]
                                                                                Score_rec_hash=Run_Result[3]
                                                                                #Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                            else:
                                                                                Run_Result=many_RD_Process(copy_num_a,run_flag)
                                                                                if Run_Result=='Error': continue
                                                                                Best_Letter_Rec=Run_Result[0]
                                                                                Best_Score_Rec=Run_Result[1]
                                                                                run_flag=Run_Result[2]
                                                                                Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                            elif len(Full_Info[9])==2 and deterministic_flag==0:
                                                                bl2_flag=0
                                                                for keyCNE in Copy_num_estimate.keys():
                                                                    if not Copy_num_estimate[keyCNE]<2:
                                                                        bl2_flag+=1
                                                                if bl2_flag==0:
                                                                    Run_Result=two_block_RD_Process(GC_para_dict,BP_para_dict,run_flag)
                                                                    if Run_Result=='Error': continue
                                                                    Best_Letter_Rec=Run_Result[0]
                                                                    Best_Score_Rec=Run_Result[1]
                                                                    run_flag=Run_Result[2]
                                                                    Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                            if run_flag==0:
                                                                if deterministic_flag==0:
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
                                                                        if M_Move_Choices=='ERROR!' and P_Move_Choices=='ERROR!':
                                                                                Move_Step-=1
                                                                                continue
                                                                        if not M_Move_Choices=='ERROR!': 
                                                                                P_IL=[]
                                                                                P_RD=[]
                                                                                P_DR=[]
                                                                                P_TB=[]
                                                                                Letter_Rec=[]
                                                                                BP_Rec=[]
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
                                                                                        Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                        if Af_Info_all==0:continue
                                                                                        Letter_Rec.append(Af_Letter)
                                                                                        BP_Rec.append(Af_BP)
                                                                                        Af_IL_Penal=Af_Info_all[0]
                                                                                        Af_RD_Rec=Af_Info_all[1]
                                                                                        Af_DR_Penal=(Af_Info_all[2])**2
                                                                                        Af_TB_Penal_a=Af_Info_all[4]
                                                                                        Af_TB_Rec=Af_Info_all[3]
                                                                                        Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                        Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                        for key in Af_Info_all[5].keys():
                                                                                                Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
                                                                                        P_IL.append(Af_IL_Penal)
                                                                                        P_RD.append(Af_RD_Penal)
                                                                                        P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                        P_TB.append(Af_TB_Penal)
                                                                                if len(P_IL)==0: continue
                                                                                Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                Regulator=1
                                                                                ILTemp=[j/Regulator for j in Regu_IL]
                                                                                RDTemp=[i for i in Regu_RD]
                                                                                DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                        if not P_Move_Choices=='ERROR!':
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
                                                                                        Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                        if Af_Info_all==0:
                                                                                                continue
                                                                                        Letter_Rec.append(Af_Letter)
                                                                                        BP_Rec.append(Af_BP)
                                                                                        Af_IL_Penal=Af_Info_all[0]
                                                                                        Af_RD_Rec=Af_Info_all[1]
                                                                                        Af_DR_Penal=(Af_Info_all[2])**2
                                                                                        Af_TB_Penal_a=Af_Info_all[4]
                                                                                        Af_TB_Rec=Af_Info_all[3]
                                                                                        Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                        Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                        for key in Af_Info_all[5].keys():
                                                                                                Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
                                                                                        P_IL.append(Af_IL_Penal)
                                                                                        P_RD.append(Af_RD_Penal)
                                                                                        P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                        P_TB.append(Af_TB_Penal)
                                                                                if len(P_IL)==0: continue
                                                                                Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                                Regulator=1
                                                                                ILTemp=[j/Regulator for j in Regu_IL]
                                                                                RDTemp=[i for i in Regu_RD]
                                                                                DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                                #best_score_rec.append(Best_Score)
                                                                    t2_sptest=time.time()
                                                                    if t2_sptest-t1_sptest<10 or bpsk1<4:
                                                                        while True:
                                                                                if Move_Step>Trail_Number: break
                                                                                if best_iterations>Local_Minumum_Number: 
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
                                                                                if M_Move_Choices=='ERROR!' and P_Move_Choices=='ERROR!':
                                                                                        Move_Step-=1
                                                                                        continue
                                                                                if not M_Move_Choices=='ERROR!':
                                                                                        P_IL=[]
                                                                                        P_RD=[]
                                                                                        P_DR=[]
                                                                                        P_TB=[]
                                                                                        Letter_Rec=[]
                                                                                        BP_Rec=[]
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
                                                                                                Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                                if Af_Info_all==0:continue
                                                                                                Letter_Rec.append(Af_Letter)
                                                                                                BP_Rec.append(Af_BP)
                                                                                                Af_IL_Penal=Af_Info_all[0]
                                                                                                Af_RD_Rec=Af_Info_all[1]
                                                                                                Af_DR_Penal=(Af_Info_all[2])**2
                                                                                                Af_TB_Penal_a=Af_Info_all[4]
                                                                                                Af_TB_Rec=Af_Info_all[3]
                                                                                                Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                                Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                                for key in Af_Info_all[5].keys():
                                                                                                        Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
                                                                                                P_IL.append(Af_IL_Penal)
                                                                                                P_RD.append(Af_RD_Penal)
                                                                                                P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                                P_TB.append(Af_TB_Penal)
                                                                                        if len(P_IL)==0: continue
                                                                                        Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                        Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                        Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                        Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                        Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                                        Regulator=1
                                                                                        ILTemp=[j/Regulator for j in Regu_IL]
                                                                                        RDTemp=[i for i in Regu_RD]
                                                                                        DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                                        #best_score_rec.append(Best_Score)
                                                                                if not P_Move_Choices=='ERROR!':
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
                                                                                                Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                                if Af_Info_all==0:
                                                                                                        continue
                                                                                                Letter_Rec.append(Af_Letter)
                                                                                                BP_Rec.append(Af_BP)
                                                                                                Af_IL_Penal=Af_Info_all[0]
                                                                                                Af_RD_Rec=Af_Info_all[1]
                                                                                                Af_DR_Penal=(Af_Info_all[2])**2
                                                                                                Af_TB_Penal_a=Af_Info_all[4]
                                                                                                Af_TB_Rec=Af_Info_all[3]
                                                                                                Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                                Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                                for key in Af_Info_all[5].keys():
                                                                                                        Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
                                                                                                P_IL.append(Af_IL_Penal)
                                                                                                P_RD.append(Af_RD_Penal)
                                                                                                P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                                P_TB.append(Af_TB_Penal)
                                                                                        if len(P_IL)==0: continue
                                                                                        Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                        Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                        Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                        Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                        Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                                        Regulator=1
                                                                                        ILTemp=[j/Regulator for j in Regu_IL]
                                                                                        RDTemp=[i for i in Regu_RD]
                                                                                        DECISION_Score=Move_Decide_2(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                                        #best_score_rec.append(Best_Score)
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
                                                                else:
                                                                    speed_test=10
                                                                    t1_sptest=time.time()
                                                                    while True:
                                                                        if Move_Step>speed_test: break
                                                                        Move_Step+=1
                                                                        M_Key_0='_'.join(['Step',str(Move_Step),'M'])
                                                                        P_Key_0='_'.join(['Step',str(Move_Step),'P'])
                                                                        Move_Sample_Pool=['delete','invert','insert']
                                                                        Move_M_P=[Move_Sample_Pool[Move_Step%len(Move_Sample_Pool)],Move_Sample_Pool[Move_Step%len(Move_Sample_Pool)]]
                                                                        M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
                                                                        P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
                                                                        if not M_Move_Choices=='ERROR!': 
                                                                            P_IL=[]
                                                                            P_RD=[]
                                                                            P_DR=[]
                                                                            P_TB=[]
                                                                            Letter_Rec=[]
                                                                            BP_Rec=[]
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
                                                                                    Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                    if Af_Info_all==0:continue
                                                                                    Letter_Rec.append(Af_Letter)
                                                                                    BP_Rec.append(Af_BP)
                                                                                    Af_IL_Penal=Af_Info_all[0]
                                                                                    Af_RD_Rec=Af_Info_all[1]
                                                                                    Af_DR_Penal=(Af_Info_all[2])**2
                                                                                    Af_TB_Penal_a=Af_Info_all[4]
                                                                                    Af_TB_Rec=Af_Info_all[3]
                                                                                    Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                    Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                    for key in Af_Info_all[5].keys():
                                                                                            Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
                                                                                    P_IL.append(Af_IL_Penal)
                                                                                    P_RD.append(Af_RD_Penal)
                                                                                    P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                    P_TB.append(Af_TB_Penal)
                                                                            if len(P_IL)==0: continue
                                                                            Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                            Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                            Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                            Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                            Regulator=1
                                                                            ILTemp=[j/Regulator for j in Regu_IL]
                                                                            RDTemp=[i for i in Regu_RD]
                                                                            DECISION_Score=Move_Decide_deterministic(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                            #best_score_rec.append(Best_Score)
                                                                        if not P_Move_Choices=='ERROR!':
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
                                                                                    Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                    if Af_Info_all==0:
                                                                                            continue
                                                                                    Letter_Rec.append(Af_Letter)
                                                                                    BP_Rec.append(Af_BP)
                                                                                    Af_IL_Penal=Af_Info_all[0]
                                                                                    Af_RD_Rec=Af_Info_all[1]
                                                                                    Af_DR_Penal=(Af_Info_all[2])**2
                                                                                    Af_TB_Penal_a=Af_Info_all[4]
                                                                                    Af_TB_Rec=Af_Info_all[3]
                                                                                    Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                    Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                    for key in Af_Info_all[5].keys():
                                                                                            Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
                                                                                    P_IL.append(Af_IL_Penal)
                                                                                    P_RD.append(Af_RD_Penal)
                                                                                    P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                    P_TB.append(Af_TB_Penal)
                                                                            if len(P_IL)==0: continue
                                                                            Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                            Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                            Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                            Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                            Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                            Regulator=1
                                                                            ILTemp=[j/Regulator for j in Regu_IL]
                                                                            RDTemp=[i for i in Regu_RD]
                                                                            DECISION_Score=Move_Decide_deterministic(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                            #best_score_rec.append(Best_Score)
                                                                    t2_sptest=time.time()
                                                                    if t2_sptest-t1_sptest<10 or bpsk1<4:
                                                                        while True:
                                                                                if Move_Step>Trail_Number: break
                                                                                if best_iterations>Local_Minumum_Number: 
                                                                                    break_Iteration_Flag=1
                                                                                    best_iterations=0
                                                                                    Best_Score_Rec=Best_Score
                                                                                    Best_Letter_Rec=Best_Letter
                                                                                    Score_rec_hash[Best_Score_Rec]=Best_Letter_Rec
                                                                                    Best_BPs_Rec=Best_BPs
                                                                                    Be_Letter=Best_Letter[0]
                                                                                    Be_BP=Best_BPs[0]
                                                                                    break_Iteration_Flag+=1
                                                                                if break_Iteration_Flag>0: break
                                                                                Move_Step+=1
                                                                                M_Key_0='_'.join(['Step',str(Move_Step),'M'])
                                                                                P_Key_0='_'.join(['Step',str(Move_Step),'P'])
                                                                                Move_Sample_Pool=['delete','invert','insert']
                                                                                Move_M_P=[Move_Sample_Pool[Move_Step%len(Move_Sample_Pool)],Move_Sample_Pool[Move_Step%len(Move_Sample_Pool)]]
                                                                                M_Move_Choices=Move_Choice_procedure_2(Move_M_P[0],Be_Letter[0],original_letters,'2m')
                                                                                P_Move_Choices=Move_Choice_procedure_2(Move_M_P[1],Be_Letter[1],original_letters,'2p')
                                                                                if not M_Move_Choices=='ERROR!':
                                                                                    P_IL=[]
                                                                                    P_RD=[]
                                                                                    P_DR=[]
                                                                                    P_TB=[]
                                                                                    Letter_Rec=[]
                                                                                    BP_Rec=[]
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
                                                                                            Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                            if Af_Info_all==0:continue
                                                                                            Letter_Rec.append(Af_Letter)
                                                                                            BP_Rec.append(Af_BP)
                                                                                            Af_IL_Penal=Af_Info_all[0]
                                                                                            Af_RD_Rec=Af_Info_all[1]
                                                                                            Af_DR_Penal=(Af_Info_all[2])**2
                                                                                            Af_TB_Penal_a=Af_Info_all[4]
                                                                                            Af_TB_Rec=Af_Info_all[3]
                                                                                            Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                            Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                            for key in Af_Info_all[5].keys():
                                                                                                    Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N]/2)
                                                                                            P_IL.append(Af_IL_Penal)
                                                                                            P_RD.append(Af_RD_Penal)
                                                                                            P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                            P_TB.append(Af_TB_Penal)
                                                                                    if len(P_IL)==0: continue
                                                                                    Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                    Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                    Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                    Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                    Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                                    Regulator=1
                                                                                    ILTemp=[j/Regulator for j in Regu_IL]
                                                                                    RDTemp=[i for i in Regu_RD]
                                                                                    DECISION_Score=Move_Decide_deterministic(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                                if not P_Move_Choices=='ERROR!':
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
                                                                                            Af_Info_all=Letter_Through_Rearrange_4(GC_para_dict,BP_para_dict,Be_Info,Af_Letter,Af_BP)
                                                                                            if Af_Info_all==0:
                                                                                                    continue
                                                                                            Letter_Rec.append(Af_Letter)
                                                                                            BP_Rec.append(Af_BP)
                                                                                            Af_IL_Penal=Af_Info_all[0]
                                                                                            Af_RD_Rec=Af_Info_all[1]
                                                                                            Af_DR_Penal=(Af_Info_all[2])**2
                                                                                            Af_TB_Penal_a=Af_Info_all[4]
                                                                                            Af_TB_Rec=Af_Info_all[3]
                                                                                            Af_TB_Penal=float(Af_TB_Penal_a)/float(num_of_reads)+float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
                                                                                            Af_RD_Penal=RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chr,Af_RD_Rec,Af_Letter)
                                                                                            for key in Af_Info_all[5].keys():
                                                                                                    Af_RD_Penal+=Prob_Norm(Af_Info_all[5][key],0,GC_Var_Coverage[chrom_N])
                                                                                            P_IL.append(Af_IL_Penal)
                                                                                            P_RD.append(Af_RD_Penal)
                                                                                            P_DR.append(Af_DR_Penal/num_of_read_pairs)
                                                                                            P_TB.append(Af_TB_Penal)
                                                                                    if len(P_IL)==0: continue
                                                                                    Regu_IL=[P_IL[i]*(1+DR_Weight*P_DR[i]) for i in range(len(P_IL))]
                                                                                    Regu_RD=[P_RD[i]+P_TB[i] for i in range(len(P_RD))]
                                                                                    Regu_IL=[(i-IL_GS)*K_IL_new for i in Regu_IL]
                                                                                    Regu_RD=[i-RD_GS for i in Regu_RD]
                                                                                    Regulator=numpy.median(Regu_IL)/numpy.median(Regu_RD)
                                                                                    Regulator=1
                                                                                    ILTemp=[j/Regulator for j in Regu_IL]
                                                                                    RDTemp=[i for i in Regu_RD]
                                                                                    DECISION_Score=Move_Decide_deterministic(ILTemp,RDTemp,GC_Var_Coverage)
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
                                                                                    #best_score_rec.append(Best_Score)
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
                                                            temp_Full_Info=original_bp_let_produce(chr_letter_bp,bps2)
                                                            original_letters=temp_Full_Info[1]
                                                            original_bp_list=temp_Full_Info[0]
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
    if function_name=='SVIntegrate':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','help=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_svintegrate()
        else:
            def add_csv_info(csv1,flag_sex,k1,k2):
                if flag_sex==1:
                    del_let=[csv1[0],[]]
                    inv_let=[csv1[1],[]]
                    dup_let=[csv1[2],[]]
                elif flag_sex==2:
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
                        del_bp.append(bp_to_hash(k3,del_let[0]),chromos)
                    else:
                        del_bp.append([])
                    if not del_let[1]==[]:
                        del_bp.append(bp_to_hash(k3,del_let[1]),chromos)
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
                tempa=bp_to_hash(k3[:-1],del_let[0],chromos)
                tempb=bp_to_hash(k3[:-1],del_let[1],chromos)
                for k1 in tempa:
                    if k1 in tempb:
                        tempc='hom'
                        tempb.remove(k1)
                    else:
                        tempc='heta'
                    if not k1[0] in del1.keys():
                        del1[k1[0]]=[]
                    del1[k1[0]].append(k1[1:]+[tempc,k3[-1],';'.join(k3[:-1]+['S'])])
                for k1 in tempb:
                    if not k1[0] in del1.keys():
                        del1[k1[0]]=[]
                    del1[k1[0]].append(k1[1:]+['hetb',k3[-1],';'.join(k3[:-1]+['S'])])
            def del_csv_info_add(k3,del_let):
                tempa=bp_to_hash(k3[:-1],del_let[0],chromos)
                tempb=bp_to_hash(k3[:-1],del_let[1],chromos)
                for k1 in tempa:
                    if k1 in tempb:
                        tempc='hom'
                        tempb.remove(k1)
                    else:
                        tempc='heta'
                    if not k1[0] in del1.keys():
                        del1[k1[0]]=[]
                    del1[k1[0]].append(k1[1:]+[tempc,k3[-1],';'.join(k3[:-1]+['C'])])
                for k1 in tempb:
                    if not k1[0] in del1.keys():
                        del1[k1[0]]=[]
                    del1[k1[0]].append(k1[1:]+['hetb',k3[-1],';'.join(k3[:-1]+['C'])])
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
                            temp=bp_to_hash(k3[:-1],[i for i in k4[0]],chromos)
                            for k5 in temp:
                                if not k5[0] in disperse_dup.keys():
                                    disperse_dup[k5[0]]=[]
                                if k4[1]>1:
                                    disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['S']),k4[1]])
                        elif dup_subtype_current=='Disperse':
                            temp=bp_to_hash(k3[:-1],[i for i in k4[0]],chromos)
                            for k5 in temp:
                                if not k5[0] in disperse_dup.keys():
                                    disperse_dup[k5[0]]=[]
                                if k4[1]>1:
                                    disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['S']),k4[1]])
            def dup_info_add(k3,dup_let):
                for k2x in dup_let:
                    for k4 in k2x:
                        temp=bp_to_hash(k3[:-1],[i for i in k4],chromos)
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
                        temp=bp_to_hash(k3[:-1],[i for i in k4[0]],chromos)
                        for k5 in temp:
                            if not k5[0] in dup1.keys():
                                dup1[k5[0]]=[]
                            if k4[1]>1:
                                dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['S']),k4[1]])
            def Define_Default_SVIntegrate():
                global score_Cff
                if not '--qc-structure' in dict_opts:
                    score_Cff=-20
                else:
                    score_Cff=int(dict_opts['--qc-structure'])
            def disperse_dup_info_2_add(k3,dup_let):
                temprec=-1
                for k2x in dup_let:
                    temprec+=1
                    hetx=['heta','hetb'][temprec]
                    for k4 in k2x:
                        temp=bp_to_hash(k3[:-1],[i for i in k4[0]],chromos)
                        for k5 in temp:
                            if not k5[0] in disperse_dup.keys():
                                disperse_dup[k5[0]]=[]
                            if k4[1]>1:
                                disperse_dup[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['S']),k4[1]])
            def dup_csv_info_2_add(k3,dup_let):
                temprec=-1
                for k2x in dup_let:
                    temprec+=1
                    hetx=['heta','hetb'][temprec]
                    for k4 in k2x:
                        temp=bp_to_hash(k3[:-1],[i for i in k4[0]],chromos)
                        for k5 in temp:
                            if not k5[0] in dup1.keys():
                                dup1[k5[0]]=[]
                            if k4[1]>1:
                                dup1[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['C']),k4[1]])
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
                    ks1=ka1.split(';')[0]
                    ks2=';'.join(ka1.split(';')[:-2]+[ka1.split(';')[-1]])
                    SV_Score=float(ka1.split(';')[-2])
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
            def inv_csv_info_add(k3,inv_let):
                temprec=-1
                for k2x in inv_let:
                    temprec+=1
                    hetx=['heta','hetb'][temprec]
                    for k4 in k2x:
                        temp=bp_to_hash(k3[:-1],[i for i in k4],chromos)
                        for k5 in temp:
                            if not k5[0] in inv1.keys():
                                inv1[k5[0]]=[]
                            inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['C'])])
            def inv_info_add(k3,inv_let):
                temprec=-1
                for k2x in inv_let:
                    temprec+=1
                    hetx=['heta','hetb'][temprec]
                    for k4 in k2x:
                        temp=bp_to_hash(k3[:-1],[i for i in k4],chromos)
                        for k5 in temp:
                            if not k5[0] in inv1.keys():
                                inv1[k5[0]]=[]
                            inv1[k5[0]].append(k5[1:]+[hetx,k3[-1],';'.join(k3[:-1]+['S'])])
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
            def ROC_writing(filename,hash):
                fo=open(filename,'w')
                for k1 in hash.keys():
                    for k2 in hash[k1].keys():
                        for k3 in chromos:
                            if k3 in hash[k1][k2].keys():
                                print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+hash[k1][k2][k3]])
                fo.close()
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
                        let1=bp_to_let([pin1],chromos)
                        if not let1==0:
                            let2='/'.join(sorted(pin2[0].split('/')))
                            if not let1 in sv_info.keys():
                                sv_info[let1]={}
                            if not let2 in sv_info[let1].keys():
                                sv_info[let1][let2]=[]
                            if not pin1 in sv_info[let1][let2]:
                                sv_info[let1][let2].append(pin1+[float(pin4[-1])-float(pin3[-1])])
                fin.close()
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
            def tra_csv_info_add(k1,k2):
                for k3 in sv_info[k1][k2]:
                    SV_ID=';'.join([str(i) for i in k3]+['C'])
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
                    SV_ID=';'.join([str(i) for i in k3]+['S'])
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
                        if os.path.isdir(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]):
                            InputPath.append(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1])
                            print 'Reading Result from default path: '+workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]
                        else:
                            print 'Error: please specify input path using --input-path'
                    ref_path=workdir+'reference_SVelter/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if '--reference' in dict_opts.keys():
                        ref_file=dict_opts['--reference']
                        ref_path='/'.join(ref_file.split('/')[:-1])+'/'
                        ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed'
                    else:
                        if not '--prefix' in dict_opts.keys():
                            print 'Warning: output file name not specified. output file: '+workdir+'Output.vcf'
                            output_file=workdir+'Output.vcf'
                        else:
                            output_file=dict_opts['--prefix']+'.vcf'
                        time1=time.time()
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
                            write_VCF_header(output_file,time)
                            write_VCF_main(output_file,sv_out,chromos,ref)
                        time2=time.time()
                        print 'SVIntegrate Complete !'
                        print 'Time Consuming: '+str(time2-time1)
    if function_name=='PredefinedBP':
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[2:],'o:h:S:',['deterministic-flag=','help=','input-bed=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters_predefinedbp()
        else:
            import time
            import datetime
            if not '--input-bed' in dict_opts.keys():
                print 'Error: please specify predefined breakpoints using --input-bed'
            else:
                def Code_Files_Define():
                    global input_bed
                    input_bed=dict_opts['--input-bed']
                    global workdir
                    workdir=path_modify(dict_opts['--workdir'])
                    global Code_File
                    global Code0_Function
                    global Code1_Function
                    global Code2_Function
                    global Code2_Predefined_Function
                    global Code3_Function
                    global Code4_Function
                    global Code5_Function
                    global RCode_Path
                    global Code1a_file
                    global Code1d_file
                    global Code1d2_file
                    Code_File=script_name
                    Code0_Function='Setup'
                    Code1_Function='NullModel'
                    Code2_Function='BPSearch'
                    Code2_Predefined_Function='BPSearch_Predefined'
                    Code3_Function='BPIntegrate'
                    Code4_Function='SVPredict'
                    Code5_Function='SVIntegrate'
                    RCode_Path=workdir+'reference_SVelter/'
                    Code1a_file=RCode_Path+'SVelter1.NullModel.Figure.a.r'
                    Code1d_file=RCode_Path+'SVelter1.NullModel.Figure.b.r'
                    Code1d2_file=RCode_Path+'SVelter1.NullModel.Figure.c.r'
                def Define_Default_AllInOne():
                    global deterministic_flag
                    deterministic_flag=0
                    if '--deterministic-flag' in dict_opts.keys():
                        deterministic_flag=int(dict_opts['--deterministic-flag'])
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
                    global Local_Minumum_Number
                    Local_Minumum_Number=100
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
                def run_SVelter1_chrom_predefine(sin_bam_file):
                    os.system(r'''%s %s --keep-temp-files %s --keep-temp-figs %s --null-model %s --workdir %s --sample %s --out-path %s'''%(Code_File,Code1_Function,KeepFile,KeepFigure,model_comp,workdir,sin_bam_file,NullModel_out_folder)) 
                def run_SVelter1_Single_chrom_predefine(sin_bam_file,chromos_single):
                    os.system(r'''%s %s --keep-temp-files %s --keep-temp-figs %s --null-model %s --workdir %s --sample %s --chromosome %s --out-path %s'''%(Code_File,Code1_Function,KeepFile,KeepFigure,model_comp,workdir,sin_bam_file,chromos_single,NullModel_out_folder)) 
                def run_SVelter2_chrom_predefine(chrom_name,sin_bam_file,ILCff_STD_Time):
                    os.system(r'''%s %s --chromosome %s --workdir %s --sample %s --null-model %s -S %s --out-path %s'''%(Code_File,Code2_Predefined_Function,chrom_name,workdir,sin_bam_file,model_comp,ILCff_STD_Time,BPPredict_out_folder))
                def run_SVelter3_chrom_predefine(sin_bam_file,out_folder):
                    os.system(r'''%s %s --batch %s --workdir %s --sample %s --bp-path %s'''%(Code_File,Code3_Function,dict_opts['--batch'],workdir,sin_bam_file,BPPredict_out_folder)) 
                def run_SVelter4_chrom(txt_name,sin_bam_file):
                    os.system(r'''%s %s --workdir %s --bp-file %s --sample %s --num-iteration %s --ploidy %s --null-model %s --deterministic-flag %s'''%(Code_File,Code4_Function,workdir,txt_name,sin_bam_file,str(Trail_Number),str(Ploidy),model_comp,deterministic_flag))
                    print txt_name+' done!'
                def run_SVelter5_chrom(path2,out_vcf):
                    os.system(r'''%s %s --workdir %s --input-path %s --prefix %s'''%(Code_File,Code5_Function,workdir,path2,out_vcf))
                def SamplingPercentage_read_in():
                    if '--null-copyneutral-perc' in dict_opts.keys():
                        SamplingPercentage=float(dict_opts['--null-copyneutral-perc'])
                    else:
                        SamplingPercentage=0.001
                    return SamplingPercentage
                def main():
                    Code_Files_Define()
                    Define_Default_AllInOne()
                    if '--sample' in dict_opts.keys():
                        bam_path='/'.join(dict_opts['--sample'].split('/')[:-1])+'/'
                        bam_files=[dict_opts['--sample']]
                        bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                    else:
                        bam_path=path_modify(dict_opts['--samplePath'])
                        bam_files=[]
                        for file in os.listdir(bam_path):
                            if file.split('.')[-1]==bam_files_appdix:
                                bam_files.append(bam_path+file)
                    ref_path=workdir+'reference_SVelter/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed '
                    else:
                        global whole_genome
                        global len_genome
                        [whole_genome,len_genome]=calculate_len_genome(ref)
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
                            print 'Warning: please make sure the reference file matches the sample file'
                        chr_flag=0
                        if 'chr' in chr_ref_check[0]:
                            chr_flag=1
                        SamplingPercentage=float(SamplingPercentage_read_in())
                        cn2_file=cn2_file_read_in(dict_opts,workdir)
                        ex_file=ex_file_read_in(dict_opts,workdir)
                        cn2_length=int(cn2_length_readin(dict_opts))
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
                            out_vcf=workdir+dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'.vcf')
                            out_svelter=workdir+dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'.svelter')
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
                            global NullModel_out_folder
                            global BPPredict_out_folder
                            global bp_files_out_folder
                            BPPredict_out_folder=workdir+'BreakPoints.'+'.'.join(sin_bam_file.split('/')[-1].split('.')[:-1])+'.predefinedBP.'+'.'.join(dict_opts['--input-bed'].split('/')[-1].split('.')[:-1])+'/'
                            NullModel_out_folder=workdir+'NullModel.'+'.'.join(sin_bam_file.split('/')[-1].split('.')[:-1])+'.predefinedBP.'+'.'.join(dict_opts['--input-bed'].split('/')[-1].split('.')[:-1])+'/'
                            bp_files_out_folder=workdir+'bp_files.'+'.'.join(sin_bam_file.split('/')[-1].split('.')[:-1])+'.predefinedBP.'+'.'.join(dict_opts['--input-bed'].split('/')[-1].split('.')[:-1])+'/'
                            running_time=[]
                            print ' '
                            print 'Step1: Running null parameters for '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')+' ...'
                            time1=time.time()
                            if len(chromos)>1:
                                run_SVelter1_chrom_predefine(sin_bam_file)
                            elif len(chromos)==1:
                                run_SVelter1_Single_chrom_predefine(sin_bam_file,chromos[0])
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Null model built for '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step2: Integrate predefined breakpoitns of sample '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')+' ...'
                            time1=time.time()
                            for x in chromos:
                                print x
                                run_SVelter2_chrom_predefine(x,sin_bam_file,ILCff_STD_Time)
                            if os.path.isfile(input_bed):
                                bed_info=bed_readin(input_bed)
                                path_mkdir(BPPredict_out_folder)
                                bed_write(bed_info,BPPredict_out_folder,sin_bam_file.split('/')[-1],input_bed)
                            else:
                                print 'Error: predefined breakpoints file not exist !'
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Breakpointse set for sample:'+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step3: Integrating breakpoints ... '
                            if not '--batch' in dict_opts.keys():
                                dict_opts['--batch']='0'
                            time1=time.time()
                            run_SVelter3_chrom_predefine(sin_bam_file,BPPredict_out_folder)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Break points cluster done for sample:'+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step4: Resolving structure ... '
                            time1=time.time()
                            for k3 in os.listdir(bp_files_out_folder):
                                if k3.split('.')[-1]=='txt':
                                    run_SVelter4_chrom(bp_files_out_folder+k3,sin_bam_file)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Structure resolved !'
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step5: Integrating results in VCF file: '+out_vcf+' ... '
                            time1=time.time()
                            run_SVelter5_chrom(workdir+bp_files_out_folder,'.'.join(out_vcf.split('.')[:-1]))
                            time2=time.time() 
                            running_time.append(time2-time1)
                            if temp_inter_replace==0:
                                print out_vcf+' completed! '
                                print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print 'Total Running Time:'+' '.join([str(i) for i in running_time])
                        if os.path.isfile(out_vcf):
                            os.system(r'''rm -r %s'''%(NullModel_out_folder))
                            os.system(r'''rm -r %s'''%(BPPredict_out_folder))
                            os.system(r'''rm -r %s'''%(TXTPath))
            main()
    if not function_name in ['BPSearch_Predefined','PredefinedBP','Setup','NullModel','BPSearch','BPIntegrate','SVPredict','SVIntegrate','Clean']:
        import glob
        import getopt
        opts,args=getopt.getopt(sys.argv[1:],'o:h:S:',['deterministic-flag=','help=','long-insert=','prefix=','batch=','sample=','workdir=','reference=','chromosome=','exclude=','copyneutral=','ploidy=','svelter-path=','input-path=','null-model=','null-copyneutral-length=','null-copyneutral-perc=','null-random-length=','null-random-num=','null-random-length=','null-random-num=','qc-align=','qc-split=','qc-structure=','qc-map-tool=','qc-map-file=','split-min-len=','read-length=','keep-temp-files=','keep-temp-figs=','bp-file=','num-iteration=','keep-interval-files='])
        dict_opts=dict(opts)
        if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
            readme.print_default_parameters()
        else:
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
                Code0_Function='Setup'
                Code1_Function='NullModel'
                Code2_Function='BPSearch'
                Code3_Function='BPIntegrate'
                Code4_Function='SVPredict'
                Code5_Function='SVIntegrate'
                RCode_Path=workdir+'reference_SVelter/'
                Code1a_file=RCode_Path+'SVelter1.NullModel.Figure.a.r'
                Code1d_file=RCode_Path+'SVelter1.NullModel.Figure.b.r'
                Code1d2_file=RCode_Path+'SVelter1.NullModel.Figure.c.r'
            def check_scripts(Code_path):
                flag=0
                out=[]
                Code0_file=Code_path+'SVelter0.Ref.Setup.py'
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
                Code0_file=Code_path+'SVelter1.NullModel.Figure.b.r'
                if not os.path.isfile(Code0_file):
                    flag+=1
                    out.append(Code0_file)
                Code0_file=Code_path+'SVelter1.NullModel.Figure.c.r'
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
            def Define_Default_AllInOne():
                global deterministic_flag
                deterministic_flag=0
                if '--deterministic-flag' in dict_opts.keys():
                    deterministic_flag=int(dict_opts['--deterministic-flag'])
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
                os.system(r'''%s %s --workdir %s --bp-file %s --sample %s --num-iteration %s --ploidy %s --null-model %s --deterministic-flag %s'''%(Code_File,Code4_Function,workdir,txt_name,sin_bam_file,str(Trail_Number),str(Ploidy),model_comp,deterministic_flag))
                print txt_name+' done!'
            def run_SVelter5_chrom(path2,out_vcf):
                os.system(r'''%s %s --workdir %s --input-path %s --prefix %s'''%(Code_File,Code5_Function,workdir,path2,out_vcf))
            def SamplingPercentage_read_in():
                global SamplingPercentage
                if '--null-copyneutral-perc' in dict_opts.keys():
                    SamplingPercentage=float(dict_opts['--null-copyneutral-perc'])
                else:
                    SamplingPercentage=0.001
            def clean_path(path):
                if os.path.isdir(path):
                    os.system(r'''rm -r %s'''%(path))
            def global_para_declaration_all():
                    global whole_genome
                    global len_genome
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
            Define_Default_AllInOne()
            global_para_declaration_all()
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
                        bam_files_appdix=dict_opts['--sample'].split('.')[-1]
                    else:
                        bam_path=path_modify(dict_opts['--samplePath'])
                        bam_files=[]
                        for file in os.listdir(bam_path):
                            if file.split('.')[-1]==bam_files_appdix:
                                bam_files.append(bam_path+file)
                    ref_path=workdir+'reference_SVelter/'
                    ref_file=ref_path+'genome.fa'
                    ref_index=ref_file+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed '
                    else:
                        [whole_genome,len_genome]=calculate_len_genome(ref_file)
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
                            print 'Warning: please make sure the reference file matches the sample file'
                        chr_flag=0
                        if 'chr' in chr_ref_check[0]:
                            chr_flag=1
                        SamplingPercentage_read_in()
                        cn2_file=cn2_file_read_in(dict_opts,workdir)
                        ex_file=ex_file_read_in(dict_opts,workdir)
                        cn2_length=int(cn2_length_readin(dict_opts))
                        Gap_Refs=[ex_file]
                        if not os.path.isfile(cn2_file):
                            print 'Error: CN2 file not correctly setup!'
                        if not os.path.isfile(ex_file):
                            random_produce_exclude_region(ex_file,chromos)
                        if '--prefix' in dict_opts.keys():
                            out_vcf=dict_opts['--prefix']+'.vcf'
                            out_svelter=dict_opts['--prefix']+'.svelter'
                        else:
                            out_vcf=workdir+dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'.vcf')
                            out_svelter=workdir+dict_opts['--sample'].split('/')[-1].replace('.'+bam_files_appdix,'.svelter')
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
                            print 'Step1: Running null parameters for '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')+' ...'
                            time1=time.time()
                            if len(chromos)>1:
                                run_SVelter1_chrom(sin_bam_file)
                            elif len(chromos)==1:
                                run_SVelter1_Single_chrom(sin_bam_file,chromos[0])
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Null model built for '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step2: Searching for BreakPoints of sample '+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')+' ...'
                            time1=time.time()
                            for x in chromos:
                                print x
                                run_SVelter2_chrom(x,sin_bam_file,ILCff_STD_Time)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Break points searching done for sample:'+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step3: Integrating breakpoints ... '
                            if not '--batch' in dict_opts.keys():
                                dict_opts['--batch']='0'
                            time1=time.time()
                            run_SVelter3_chrom(sin_bam_file)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Break points cluster done for sample:'+sin_bam_file.split('/')[-1].replace('.'+bam_files_appdix,'')
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step4: Resolving structure ... '
                            time1=time.time()
                            for k1 in os.listdir(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'):
                                if k1.split('.')[-1]=='txt':
                                    run_SVelter4_chrom(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/'+k1,sin_bam_file)
                            time2=time.time()
                            running_time.append(time2-time1)
                            print 'Structure resolved !'
                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            time1=time.time()
                            run_SVelter5_chrom(workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]+'/','.'.join(out_vcf.split('.')[:-1]))
                            time2=time.time() 
                            running_time.append(time2-time1)
                            if temp_inter_replace==0:
                                print out_vcf+' completed! '
                                print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print 'Total Running Time:'+' '.join([str(i) for i in running_time])
                        #if os.path.isfile(out_vcf):
                        NullPath=workdir+'NullModel.'+dict_opts['--sample'].split('/')[-1]
                        BPPath=workdir+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]
                        TXTPath=workdir+'bp_files.'+dict_opts['--sample'].split('/')[-1]
                        if not '--keep-interval-files' in dict_opts.keys():
                            clean_path(NullPath)
                            clean_path(BPPath)
                            clean_path(TXTPath)
                        elif dict_opts['--keep-interval-files']=='FALSE':
                            clean_path(NullPath)
                            clean_path(BPPath)
                            clean_path(TXTPath)                            