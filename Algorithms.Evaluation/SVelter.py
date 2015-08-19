#!/usr/bin/python

#!Python
#Usage:
#python SVelter.py [options]

#--KeepFigure: yes/Y/Yes/y : will keep figures showing distribution of parameters ; otherwise, no/N/No/n ; default: no
#--KeepFile: yes/Y/Yes/y : will keep large detailed useless statistical files ; otherwise, no/N/No/n ; default: no
#--SplitLength: 1-10:minimum length of clipped region that a read can be considered to be soft clipped;50-95%: quantile cutoff; default: 10% quantile cutoff
#--alignQC: minimum mappability of a certain genomic region required; default: 0.0; 
#--SPCff: minimum number of clipped reads required to preliminary define a break point; default:1/10*mean_coverge
#--LNCff: minimum number of reads required to form a cluster to define a break point; default: 2*/10*mean_coverge
#--qc: minimum requirement for mapping quality; default:20
#--cn2length: minimum length of cn2 regions used for default para calculate; default: 2000
#--SamplingPerc: percentage to Sample from cn2 regions or other Ref regions to build null model; default: 10%
#--chr:the certain chromosome to focus on;name of chromosome has to consist with the name in fasta file default: whole genome recored in Ref genome.fa
#-o: output path; defult: workdir/NullModel for code1; workdir/BreakPoints for code1b, code2
#--filter_bed:folder containing other bed files that are used to select regions on genome for running analysis; default: workdir/tools/wgEncodeDacMapabilityConsensusExcludable.bed / workdir/tools/wgEncodeDukeMapabilityRegionsExcludable.bed , both downloaded from ucsc
#--scriptPath: folder containing code; default: workdir/code
#--pbsSample: absolute path of a pbs script used in users' server
#--out_Null: output path of step1 (build null model); default: workdir/NullModel
#--out_Readable:output path of step 1b (produce mappable regions in bed format); default: workdir/BreakPoints 
#--out_BP: output path of step2 (.LN, .SP, file containing breadk points); default: workdir/BreakPoints 
#--out_txt: output path of step3 (integrated bps in txt files); default: workdir/BreakPoints/bp_files
#--out_struc: output path of step4(txt files containing solved structure); default: workdir/BreakPoints/bp_files
#--out_vcf:output path of step5(integrated results in vcf format); default: workdir/vcf
#For debug use only
#command='SVelter.py --ppre /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/ --ScriptPath /mnt/EXT/Mills-data/xuefzhao/SVelter -s /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/BamFiles/NA12878_S1.bam --ref /mnt/EXT/Mills-data/xuefzhao/projects/Pedigree1463.axiom/reference/genome.fa'
#sys.argv=command.split()
import os
import sys
import getopt
import numpy
import re
import random
import pickle
import time
import datetime
from multiprocessing import Pool
pool = Pool(processes=4)
opts,args=getopt.getopt(sys.argv[1:],'o:h:s:x:',['ScriptPath=','NullModel=','cn=','ex=','help=','ppre=','PathSVelter=','chr=','ref=','NullSplitLength=','NullSampleLength=','NullSampleNumber=','ToolMappingQ=','FileMappingQ=','NullSamplePercentage=','SplitLength=','BPSPCff=','BPLNCff=','BPAlignQC=','BPAlignQCFlank=','ReadLen=','SPCluLen=','QCSplit=','QCAlign=','KeepFigure=','KeepFile='])
dict_opts=dict(opts)

if not '--NullModel' in dict_opts.keys():
    model_comp='S'
else:
    if dict_opts['--NullModel'] in ['S','Simple']:
        model_comp='S'
    else:
        model_comp='C'

if '--QCAlign' in dict_opts.keys():
    QCAlign=int(dict_opts['--QCAlign'])
else:
    QCAlign=20

if '--QCSplit' in dict_opts.keys():
    QCSplit=int(dict_opts['--QCSplit'])
else:
    QCSplit=20

if '--NullSplitLength' in dict_opts.keys():
    NullSplitLen_perc=int(dict_opts['--NullSplitLength'])
else:
    NullSplitLen_perc=0.9

if '--KeepFile' in dict_opts.keys():
    KeepFile=dict_opts['--KeepFile']
else:
    KeepFile='No'

if '--KeepFigure' in dict_opts.keys():
    KeepFigure=dict_opts['--KeepFigure']
else:
    KeepFigure='No'

if '--NIteration' in dict_opts.keys():
    Trail_Number=int(dict_opts['--NIteration'])
else:
    Trail_Number=10000

if '--MoveFlag' in dict_opts.keys():
    MoveFlag=int(dict_opts['--MoveFlag'])
else:
    MoveFlag=2

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
    if '--NullSamplePercentage' in dict_opts.keys():
        SamplingPercentage=float(dict_opts['--NullSamplePercentage'])
    else:
        SamplingPercentage=0.001
    return SamplingPercentage

def cn2_file_read_in():
    if '--cn' in dict_opts.keys():
        cn2_file=dict_opts['--cn']
    else:
        cn2_file=workdir+'tools/CN2.bed'
    return cn2_file

def ex_file_read_in():
    if '--ex' in dict_opts.keys():
        ex_file=dict_opts['--ex']
    else:
        ex_file=workdir+'tools/Exclude.bed'
    return ex_file

def cn2_length_read_in():
    if '--NullBedLength' in dict_opts.keys():
        cn2_length=int(dict_opts['--NullBedLength'])
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
    os.system(r'''%s --ppre %s --ref %s --ex %s -s %s --chr %s'''%(Code0_file,workdir,ref_file,ex_file,sin_bam_file,chrom_name))

def run_SVelter2_chrom(chrom_name):
    os.system(r'''%s --chr %s --ppre %s -s %s --ref %s --NullModel %s'''%(Code2_file,chrom_name,workdir,sin_bam_file,ref_file,model_comp))

def run_SVelter4_chrom(txt_name):
    os.system(r'''%s --ppre %s -f %s -s %s --NIteration %s --ref %s --MoveFlag %s --NullModel %s'''%(Code4_file,workdir,txt_name,sin_bam_file,str(Trail_Number),ref_file,str(MoveFlag),model_comp))

if dict_opts=={} or dict_opts.keys()==['-h'] or dict_opts.keys()==['--help']:
    print 'SVelter-0.1          Last Update:2014-10-27'
    print 'Required Parameters:'
    print '--ppre: writable working directory; default: ./(current directory)'
    print '--ScriptPath:folder containing code; default: workdir/code'
    print '--ref: absolute path of Reference genome; default: workdir/Reference/genome.fa; ref path shuold be writable'
    print '-s: absolute path of one sample need to run; default: all Samples under workdir/BamFiles'
    print '--SamplePath:absolute path of multiple samples to run; if -s specified, this parameter would be ignored'
    print ' '
    print 'Optional Parameters:'
    print '--cn2: absolute path of cn2 files; or other Ref bed files that contains regions users are interested in; default: workdir/tools/CN2.bed'
    print '--ex: absolute path of exculde files; default: workdir/tools/Exclude.bed'
    print '--MoveFlag: 0:heterozygous 1:homozygous 2: diploid module'
    print '-o, absolute path of output vcf file. eg: .../SVelter/sample.vcf'
    print '--ex, bed file containing genomic regions to exclude. eg: Excludable.bed'
    print '--cn, bed file containing genomic regions where no SVs have yet been reported, based on which null model were built'
    print '--NullModel, specify which stat model to be fitted on each parameter. if --NullModel==C / Complex, negative bimodal distribution will be fitted to insertlenth; else, normal will be used'
    print '--NullBedLength, minimum requirement for length of regions used to build null model; default: 2000bp'
    print '--NullSamplePercentage, sampling percentage of regions from Nullbed'
    print '--NullSampleLength, if not --NullBed provided, SVelter will randomly Sample regions on genome to build null model; --NullSampleLength specify the length of region picked; default: 5kb'
    print '--NullSampleNumber, if not --NullBed provided, SVelter will randomly Sample regions on genome to build null model; --NullSampleNumber specify the number of region picked; default: 10k'
    print '--QCAlign, minimum alignment quality required for mapped reads in bam file; default: 20'
    print '--QCSplit, minimum alighment of clipped parts of reads considered as a soft clip; default: 20'
    print '--NullSplitLength, the minumum length of cutoff considered as split; default:10% of read length '
    print '--KeepFile, whether to keep interval files during the program; to remove interval file, specify the parameter with "no"/"N"/"No"/"n"; else, interval files will be kept'
    print '--KeepFigure, whether to keep interval figures during the program; to remove interval file, specify the parameter with "no"/"N"/"No"/"n"; else, interval files will be kept'
else:
    if not '--ppre' in dict_opts.keys():
        print 'Error: please specify working directory using: --ppre'
    else:
        workdir=path_modify(dict_opts['--ppre'])
        if not os.path.isdir(workdir):
            print 'Error: working directory does not exit!'
        if not '--ScriptPath' in dict_opts.keys():
            Code_path='/'.join(sys.argv[0].split('/')[:-1])+'/'
        else:
            Code_path=path_modify(dict_opts['--ScriptPath'])
        Code0_file=Code_path+'SVelter0.Ref.Index.py'
        Code1_file=Code_path+'SVelter1.NullModel.py'
        Code1a_file=Code_path+'SVelter1.NullModel.Figure.a.r'
        Code1d_file=Code_path+'SVelter1.NullModel.Figure.d.r'
        Code1d2_file=Code_path+'SVelter1.NullModel.Figure.d2.r'
        Code2_file=Code_path+'SVelter2.BP.Searching.py'
        Code3_file=Code_path+'SVelter3.BPIntegrate.py'
        Code4_file=Code_path+'SVelter4.StructureResolvation.py'
        Code5_file=Code_path+'SVelter5.result.integrate.py'
        script_test=check_scripts(Code_path)
        if not script_test==[]:
            print 'Error: please make sure all required scripts are under working directory'
        else:
            if not '-s' in dict_opts.keys() and not '--SamplePath' in dict_opts.keys():
                print 'Error: please specify input file using -s'
            else:
                if '-s' in dict_opts.keys():
                    bam_path='/'.join(dict_opts['-s'].split('/')[:-1])+'/'
                    bam_files=[dict_opts['-s']]
                else:
                    bam_path=path_modify(dict_opts['--SamplePath'])
                    bam_files=[]
                    for file in os.listdir(bam_path):
                        if file.split('.')[-1]=='bam':
                            bam_files.append(bam_path+file)
                if not '--ref' in dict_opts.keys():
                    print 'Error: please specify refrence genome using --ref'
                else:
                    ref_path='/'.join(dict_opts['--ref'].split('/')[:-1])+'/'
                    ref_file=dict_opts['--ref']
                    ref_index=dict_opts['--ref']+'.fai'
                    if not os.path.isfile(ref_index):
                        print 'Error: reference genome not indexed !'
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
                            if not '--NullSampleLength' in dict_opts.keys():
                                dict_opts['--NullSampleLength']=5000
                            else:
                                dict_opts['--NullSampleLength']=int(dict_opts['--NullSampleLength'])
                            if not '--NullSampleNumber' in dict_opts.keys():
                                dict_opts['--NullSampleNumber']=10000
                            else:
                                dict_opts['--NullSampleNumber']=int(dict_opts['--NullSampleNumber'])
                            cn2_length=dict_opts['--NullSampleLength']-100
                            fo=open(cn2_file,'w')
                            for i in sorted(whole_genome.keys()):
                                num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--NullSampleNumber'])
                                reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--NullSampleLength']) for j in range(num_i)]
                                for j in sorted(reg_i):
                                    print >>fo, ' '.join([i,str(j),str(j+dict_opts['--NullSampleLength']-1)])
                            fo.close()
                            SamplingPercentage=1
                        if not os.path.isfile(ex_file):
                                fo=open(ex_file,'w')
                                for chr_ex in chromos:
                                    print >>fo, ' '.join([chr_ex,'0','0'])
                                fo.close()
                        if '-o' in dict_opts.keys():
                            out_vcf=dict_opts['-o']
                        else:
                            out_vcf=workdir+dict_opts['-s'].split('/')[-1].replace('.bam','.vcf')
                            print 'Warning: output file is not specified'
                            print 'output file: '+out_vcf
                        #temp_inter=workdir+dict_opts['-s'].split('/')[-1].replace('.bam','.temp')
                        #temp_inter_replace=0
                        for sin_bam_file in bam_files:
                            running_time=[]
                            print 'Step1: Index reference genome ...'
                            time1=time.time()
                            results = [pool.apply_async(run_SVelter0_chrom, args=(x,)) for x in chromos]                            
                            time2=time.time()
                            running_time.append(time2-time1)
                            #print 'reference: '+ref_file+' indexed!'
                            #print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step2: Running null parameters for '+sin_bam_file.split('/')[-1].replace('.bam','')+' ...'
                            time1=time.time()
                            os.system(r'''%s --KeepFile %s --KeepFigure %s --NullModel %s --ppre %s -s %s --ref %s --cn %s --ex %s'''%(Code1_file,KeepFile,KeepFigure,model_comp,workdir,sin_bam_file,ref_file,cn2_file,ex_file)) 
                            time2=time.time()
                            running_time.append(time2-time1)
                            #print 'Null model built for '+sin_bam_file.split('/')[-1].replace('.bam','')
                            #print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step3: Searching for BreakPoints of sample '+sin_bam_file.split('/')[-1].replace('.bam','')+' ...'
                            time1=time.time()
                            results = [pool.apply_async(run_SVelter2_chrom, args=(x,)) for x in chromos]
                            time2=time.time()
                            running_time.append(time2-time1)
                            #print 'Break points searching done for sample:'+sin_bam_file.split('/')[-1].replace('.bam','')
                            #print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step4: clustering breakpoints ... '
                            if not '--batch' in dict_opts.keys():
                                dict_opts['--batch']='0'
                            time1=time.time()
                            os.system(r'''%s --batch %s --ppre %s -s %s --ref %s '''%(Code3_file,dict_opts['--batch'],workdir,sin_bam_file,ref_file)) 
                            time2=time.time()
                            running_time.append(time2-time1)
                            #print 'Break points cluster done for sample:'+sin_bam_file.split('/')[-1].replace('.bam','')
                            #print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print ' '
                            print 'Step5: Resolving structure ... '
                            break_flag=0
                            time1=time.time()
                            for k1 in os.listdir(workdir+'bp_files.'+dict_opts['-s'].split('/')[-1]+'/'):
                                if k1==dict_opts['-s'].split('/')[-1].replace('.bam',''):
                                    path1=workdir+'bp_files.'+dict_opts['-s'].split('/')[-1]+'/'+k1+'/'
                                    for k2 in os.listdir(path1):
                                        path2=path1+k2+'/'
                                        all_txt_files=[]
                                        for k3 in os.listdir(path2):
                                            if k3.split('.')[-1]=='txt':
                                                all_txt_files.append(path2+k3)
                                            results = [pool.apply_async(run_SVelter4_chrom, args=(x,)) for x in all_txt_files]                            
                            time2=time.time()
                            running_time.append(time2-time1)
                            if not break_flag==1:
                                print 'Structure resolved !'
                                print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                                print ' '
                                for k1 in os.listdir(workdir+'bp_files.'+dict_opts['-s'].split('/')[-1]+'/'):
                                    path1=workdir+'bp_files.'+dict_opts['-s'].split('/')[-1]+'/'+k1+'/'
                                    for k2 in os.listdir(path1):
                                        path2=path1+k2+'/'
                                        print 'Step6: Producing VCF file: '+out_vcf+' ... '
                                        time1=time.time()
                                        os.system(r'''%s --ppre %s --ref %s --RSPath %s -o %s'''%(Code5_file,workdir,ref_file,path2,out_vcf))
                                        time2=time.time() 
                                        if temp_inter_replace==0:
                                            #os.system(r'''rm %s'''%(temp_inter))
                                            print out_vcf+' completed! '
                                            print 'Time Consuming: '+str(datetime.timedelta(seconds=(time2-time1)))
                            print 'Running Time:'+' '.join([str(i) for i in running_time])
                        if os.path.isfile(out_vcf):
                            NullPath=workdir+'NullModel.'+dict_opts['-s'].split('/')[-1]
                            BPPath=workdir+'BreakPoints.'+dict_opts['-s'].split('/')[-1]
                            TXTPath=workdir+'bp_files.'+dict_opts['-s'].split('/')[-1]
                            os.system(r'''rm -r %s'''%(NullPath))
                            os.system(r'''rm -r %s'''%(BPPath))
                            os.system(r'''rm -r %s'''%(TXTPath))



