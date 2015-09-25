#!/usr/bin/env python

#!python
#command='SV.Output.Process.py vcf --input input.vcf --reference genome.fa'
#sys.argv=command.split()
import os
import sys
script_name=sys.argv[0]
if len(sys.argv)<2:
    print 'SV.Output.Process.py         Last Update:2015-08-20'
    print ''
    print 'this script is used to process output of different algorithms'
    print ''
    print 'Usage:'
    print 'SV.Output.Process.python [options]  <parameters>'
    print ' '
    print 'Options:'
    print 'vcf-to-bed:  extract simple SVs from vcf files and output in separate bed files' 
    print 'bedpe-to-bed:    extract simple SVs from bedpe files and output in separate bed files' 
    print 'Mappable-Control:    remove SVs located outside mappable regions' 
    print 'Size-Control:    filter out SVs of size outside defined range'
    print 'TRA-Control: remove SVs overlap with defined SVs' 
    print ' '
    print 'Parameters for vcf-to-bed:'
    print '--input: input file'
    print ' '
    print 'Parameters for bedpe-to-bed'
    print '--input: input file'
    print '--reference: reference.genome.fa'
    print ' '
    print 'Parameters for Mappable-Control:'
    print '--input: input file'
    print '--ref-prefix: reference.genome.fa'
    print ' '
    print 'Parameters for Size-Control:'
    print '--input: input.bed'
    print '--min-size: reference.Mappable. describing mappable regions'
    print '--max-size: reference.Mappable. describing mappable regions'
    print ' '
    print 'Parameters for TRA-Control:'
    print '--input: input.bed'
    print '--TRA-rec: TRA information kept in .rec files'
    print ' '
else:
    function_name=sys.argv[1]
    import numpy
    import getopt
    if function_name=='vcf-to-bed':
        if len(sys.argv)<3:
            print 'SV.Output.Process.py vcf-to-bed <parameters>'
            print ' '
            print 'Parameters for vcf-to-bed:'
            print '--input: input file'
            print ' '
        else:
            import numpy
            import getopt
            def end_cordi_calcu(pin):
                    #info=pin
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
            def genotype_calcu(pin):
                    gt_index=pin[8].split(':').index('GT')
                    geno=pin[9].split(':')[gt_index]
                    out=0
                    if '/' in geno:
                            if geno.split('/')==['0','0']:
                                    out='Error!'
                            elif sorted(geno.split('/'))==['0','1']:
                                    out='het'
                            elif sorted(geno.split('/'))==['1','1']:
                                    out='homo'
                    elif '|' in geno:
                            if geno.split('|')==['0','0']:
                                    out='Error!'
                            elif sorted(geno.split('|'))==['0','1']:
                                    out='het'
                            elif sorted(geno.split('|'))==['1','1']:
                                    out='homo'
                    return out
            def DEL_add(pin,sv_hash):
                    if not 'DEL' in sv_hash.keys():
                            sv_hash['DEL']={}
                    gt=genotype_calcu(pin)
                    if not gt=='Error!':
                            if not gt in sv_hash['DEL'].keys():
                                    sv_hash['DEL'][gt]=[]
                            pos=end_cordi_calcu(pin)
                            if not pos=='Error!':
                                if not pos in sv_hash['DEL'][gt]:
                                    sv_hash['DEL'][gt].append(pos)
                                    #print pos
            def DUP_add(pin,sv_hash):
                    if not 'DUP' in sv_hash.keys():
                            sv_hash['DUP']={}
                    gt=genotype_calcu(pin)
                    if not gt=='Error!':
                            if not gt in sv_hash['DUP'].keys():
                                    sv_hash['DUP'][gt]=[]
                            pos=end_cordi_calcu(pin)
                            if not pos=='Error!':
                                if not pos in sv_hash['DUP'][gt]:
                                    sv_hash['DUP'][gt].append(pos)
                                    #print pos
            def TANDEMDUP_add(pin,sv_hash):
                    if not 'DUP_TANDEM' in sv_hash.keys():
                            sv_hash['DUP_TANDEM']={}
                    gt=genotype_calcu(pin)
                    if not gt=='Error!':
                            if not gt in sv_hash['DUP_TANDEM'].keys():
                                    sv_hash['DUP_TANDEM'][gt]=[]
                            pos=end_cordi_calcu(pin)
                            if not pos=='Error!':
                                if not pos in sv_hash['DUP_TANDEM'][gt]:
                                    sv_hash['DUP_TANDEM'][gt].append(pos)
                                    #print pos
            def INV_add(pin,sv_hash):
                    if not 'INV' in sv_hash.keys():
                            sv_hash['INV']={}
                    gt=genotype_calcu(pin)
                    if not gt=='Error!':
                            if not gt in sv_hash['INV'].keys():
                                    sv_hash['INV'][gt]=[]
                            pos=end_cordi_calcu(pin)
                            if not pos=='Error!':
                                if not pos in sv_hash['INV'][gt]:
                                    sv_hash['INV'][gt].append(pos)
                                    #print pos
            def vcf_read_in(file_name,sv_hash):
                    fin=open(file_name)
                    for line in fin:
                            pin=line.strip().split()
                            if not pin[0][0]=='#':
                                    if pin[6]=='PASS':
                                            if 'DEL' in pin[4]:
                                                    DEL_add(pin,sv_hash)
                                            if 'DUP' in pin[4] and not 'TANDEM' in pin[4]:
                                                    DUP_add(pin,sv_hash)
                                            if 'INV' in pin[4]:
                                                    INV_add(pin,sv_hash)
                                            if 'DUP:TANDEM' in pin[4]:
                                                    TANDEMDUP_add(pin,sv_hash)
                    fin.close()
            def sv_detect(pin):
                    sv_type=0
                    if 'DEL' in pin[4]:
                            sv_type='DEL'
                    elif 'DUP' in pin[4] and not 'TANDEM' in pin[4]:
                            sv_type='DUP'
                    elif 'INV' in pin[4]:
                            sv_type='INV'
                    elif 'DUP:TANDEM' in pin[4]:
                        sv_type='DUP_TANDEM'
                    else:
                            for x in pin[7].split(';'):
                                   if 'SVTYPE' in x:
                                            sv_type=x.split('=')[1]
                    if sv_type=='DUP:TANDEM':
                        sv_type='DUP_TANDEM'
                    return sv_type
            def vcf_read_in_Pindel(file_name,sv_hash):
                    fin=open(file_name)
                    for line in fin:
                            pin=line.strip().split()
                            if not pin[0][0]=='#':
                                if pin[6]=='PASS':
                                    sv_type=sv_detect(pin)
                                    if not sv_type in sv_hash.keys():
                                            sv_hash[sv_type]={}
                                    gt=genotype_calcu(pin)
                                    if not gt=='Error!':
                                            if not gt in sv_hash[sv_type].keys():
                                                    sv_hash[sv_type][gt]=[]
                                            pos=end_cordi_calcu(pin)
                                            if not pos=='Error!':
                                                if not pos in sv_hash[sv_type][gt]:
                                                    sv_hash[sv_type][gt].append(pos)
                    fin.close()
            def hash_reorder(sv_hash):
                    out={}
                    for k1 in sv_hash.keys():
                            out[k1]={}
                            for k2 in sv_hash[k1].keys():
                                    for k3 in sv_hash[k1][k2]:
                                            if not k3[0] in out[k1].keys():
                                                    out[k1][k3[0]]={}
                                            if not int(k3[1]) in out[k1][k3[0]].keys():
                                                    out[k1][k3[0]][int(k3[1])]={}
                                            if not int(k3[2]) in out[k1][k3[0]][int(k3[1])].keys():
                                                    out[k1][k3[0]][int(k3[1])][int(k3[2])]=[]
                                            if not k2 in out[k1][k3[0]][int(k3[1])][int(k3[2])]:
                                                out[k1][k3[0]][int(k3[1])][int(k3[2])].append(k2)
                    return out
            def bed_writing(file_name,sv_hash):
                    for k1 in sv_hash.keys():
                            fo=open(file_name.replace('.vcf','.'+k1+'.bed'),'w')
                            for k2 in sorted(sv_hash[k1].keys()):
                                    for k3 in sorted(sv_hash[k1][k2].keys()):
                                            for k4 in sorted(sv_hash[k1][k2][k3].keys()):
                                                    tmp=[]
                                                    for k5 in sv_hash[k1][k2][k3][k4]:
                                                            if not k5 in tmp:
                                                                    tmp.append(k5)
                                                    for k5 in tmp:
                                                            print >>fo, ' '.join([str(i) for i in [k2,k3,k4,k5]])
                            fo.close()
            def bed_compare(sv_1_hash,sv_2_hash):
                    report={}
                    diff_hash={}
                    for k1 in sv_1_hash.keys():
                            if k1 in sv_2_hash.keys():
                                    report[k1]={}
                                    diff_hash[k1]={}
                                    for k2 in sv_1_hash[k1].keys():
                                            if k2 in sv_2_hash[k1].keys():
                                                    report[k1][k2]=[]
                                                    diff_hash[k1][k2]={}
                                                    bed1=[]
                                                    bed2=[]
                                                    uni1=[]
                                                    uni2=[]
                                                    for k3 in sorted(sv_1_hash[k1][k2].keys()):
                                                            for k4 in sorted(sv_1_hash[k1][k2][k3].keys()):
                                                                if not [k3,k4] in bed1:
                                                                    bed1.append([k3,k4])
                                                    for k3 in sorted(sv_2_hash[k1][k2].keys()):
                                                            for k4 in sorted(sv_2_hash[k1][k2][k3].keys()):
                                                                if not [k3,k4] in bed2:
                                                                    bed2.append([k3,k4])
                                                    report[k1][k2].append(len(bed1))
                                                    report[k1][k2].append(len(bed2))
                                                    for k3 in bed1:
                                                            flag1=0
                                                            for k4 in bed2:
                                                                    if k4[0]>k3[1]: continue
                                                                    elif k4[1]<k3[0]:continue
                                                                    else:
                                                                            if (sorted(k3+k4)[2]-sorted(k3+k4)[1])/max([k3[1]-k3[0],k4[1]-k4[0]])>0.5:
                                                                                    flag1+=1
                                                            if flag1==0:
                                                                    uni1.append(k3)
                                                    for k4 in bed2:
                                                            flag1=0
                                                            for k3 in bed1:
                                                                    if k4[0]>k3[1]: continue
                                                                    elif k4[1]<k3[0]:continue
                                                                    else:
                                                                            if (sorted(k3+k4)[2]-sorted(k3+k4)[1])/max([k3[1]-k3[0],k4[1]-k4[0]])>0.5:
                                                                                    flag1+=1
                                                            if flag1==0:
                                                                    uni2.append(k4)
                                                    diff_hash[k1][k2]['a']=uni1
                                                    diff_hash[k1][k2]['b']=uni2
            bp_collaps_CI=[-51,51]
            def hash_refilter(sv_1_hash):
                    for k1 in sv_1_hash.keys():
                                    for k2 in sv_1_hash[k1].keys():
                                            rec_temp=[]
                                            for k3 in range(len(sv_1_hash[k1][k2].keys())-1):
                                                    if sorted(sv_1_hash[k1][k2].keys())[k3+1]-sorted(sv_1_hash[k1][k2].keys())[k3]<bp_collaps_CI[1]:
                                                                    for k4 in sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3+1]].keys():
                                                                            if not k4 in sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3]].keys():
                                                                                    sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3]][k4]=sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3+1]][k4]
                                                                            else:
                                                                                    sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3]][k4]+=sv_1_hash[k1][k2][sorted(sv_1_hash[k1][k2].keys())[k3+1]][k4]
                                                                    rec_temp.append(sorted(sv_1_hash[k1][k2].keys())[k3+1])
                                            for k3 in rec_temp:
                                                    del sv_1_hash[k1][k2][k3]
                                    for k2 in sv_1_hash[k1].keys():
                                            for k3 in sv_1_hash[k1][k2].keys():
                                                    if len(sv_1_hash[k1][k2][k3])>1:
                                                            ke1=sorted(sv_1_hash[k1][k2][k3].keys())[::-1]
                                                            ke2=[ke1[0]]
                                                            for k4 in range(len(ke1)-1):
                                                                    if ke1[k4]-ke1[k4+1]<bp_collaps_CI[1]:
                                                                            continue
                                                                    else:
                                                                            ke2.append(ke1[k4+1])
                                                            for k4 in sv_1_hash[k1][k2][k3].keys():
                                                                    if not k4 in ke2:
                                                                            del sv_1_hash[k1][k2][k3][k4]
                                    for k2 in sv_1_hash[k1].keys():
                                            for k3 in sorted(sv_1_hash[k1][k2].keys()):
                                                    for k4 in sorted(sv_1_hash[k1][k2][k3].keys()):
                                                            if len(sv_1_hash[k1][k2][k3][k4])>1:
                                                                    if 'homo' in sv_1_hash[k1][k2][k3][k4]:
                                                                            sv_1_hash[k1][k2][k3][k4]=['homo']
                                                                    else:
                                                                            sv_1_hash[k1][k2][k3][k4]=['het']
                    return sv_1_hash
            opts,args=getopt.getopt(sys.argv[2:],'i:',['input=','reference='])
            dict_opts=dict(opts)
            file_in_1=dict_opts['--input']
            sv_1_hash={}                                                                                                  
            vcf_read_in(file_in_1,sv_1_hash)
            vcf_read_in_Pindel(file_in_1,sv_1_hash)
            sv_1_hash=hash_reorder(sv_1_hash)
            sv_1_hash=hash_refilter(sv_1_hash)
            bed_writing(file_in_1,sv_1_hash)
    elif function_name=='bedpe-to-bed':
        if len(sys.argv)<3:
            print 'SV.Output.Process.py bedpe-to-bed <parameters>'
            print ' '
            print 'Parameters for bedpe-to-bed'
            print '--input: input file'
            print '--reference: reference.genome.fa'
            print ' '
        else:
            import numpy
            import getopt
            opts,args=getopt.getopt(sys.argv[2:],'i:',['input=','reference='])
            dict_opts=dict(opts)
            out_hash={}
            fin=open(dict_opts['--input'])
            for line in fin:
                pin=line.strip().split()
                chrom=pin[0]
                x1=int(numpy.mean([int(pin[1]),int(pin[2])]))
                x2=int(numpy.mean([int(pin[4]),int(pin[5])]))
                sv_key=pin[10].split(':')[1][:3]
                if not sv_key in out_hash.keys():
                    out_hash[sv_key]={}
                if not chrom in out_hash[sv_key].keys():
                    out_hash[sv_key][chrom]={}
                if not x1 in out_hash[sv_key][chrom].keys():
                    out_hash[sv_key][chrom][x1]=[]
                if not x2 in out_hash[sv_key][chrom][x1]:
                    out_hash[sv_key][chrom][x1].append(x2)
            fin.close()
            chromos=[]
            ref=dict_opts['--reference']
            fin=open(ref+'.fai')
            for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
            fin.close()
            for k1 in out_hash.keys():
                fo=open(dict_opts['--input'].replace('bedpe',k1+'.bed'),'w')
                for k2 in chromos:
                    if k2 in out_hash[k1].keys():
                        for k3 in sorted(out_hash[k1][k2].keys()):
                            for k4 in sorted(out_hash[k1][k2][k3]):
                                print >>fo, '\t'.join([str(i) for i in [k2,k3,k4,'het']])
            fo.close()
    elif function_name=='Mappable-Control':
        if len(sys.argv)<3:
            print 'SV.Output.Process.py Mappable-Control <parameters>'
            print ' '
            print 'Parameters for Mappable-Control:'
            print '--input: input.bed'
            print '--ref-prefix: reference.Mappable. describing mappable regions'
            print ' '
        else:
            opts,args=getopt.getopt(sys.argv[2:],'i:',['input=','ref-prefix='])
            dict_opts=dict(opts)
            ref_hash={}
            fin=open(dict_opts['--ref-prefix'])
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in ref_hash.keys():
                    ref_hash[pin[0]]=[]
                ref_hash[pin[0]].append([int(pin[1]),int(pin[2])])
            fin.close()
            #ref_path='/'.join(dict_opts['--ref-prefix'].split('/')[:-1])
            #ref_prefix=dict_opts['--ref-prefix'].split('/')[-1]
            #ref_files=[]
            #for k1 in os.listdir(ref_path):
            #    if ref_prefix in k1 and k1.split('.')[-1]=='bed':
            #        ref_files.append(ref_path+'/'+k1)
            #ref_hash={}
            #for k1 in ref_files:
            #    chrom=k1.split('/')[-1].replace('.bed','').replace(ref_prefix+'.','')
            #    ref_hash[chrom]=[]
            #    fin=open(k1)
            #    for line in fin:
            #        pin=line.strip().split()
            #        ref_hash[chrom].append([int(pin[1]),int(pin[2])])
            #   fin.close()
            file_in=dict_opts['--input']
            file_out=file_in.replace('.bed','.Mappable.bed')
            fin=open(file_in)
            in_hash={}
            chroms=[]
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in in_hash.keys():
                    in_hash[pin[0]]=[]
                    chroms.append(pin[0])
                in_hash[pin[0]].append([int(i) for i in pin[1:3]]+pin[3:])
            fin.close()
            fo=open(file_out,'w')
            for k1 in chroms:
                if k1 in in_hash.keys() and k1 in ref_hash.keys():
                    for k2 in in_hash[k1]:
                        flag=0
                        for k3 in ref_hash[k1]:
                            if k3[0]-1<k2[0] and k3[1]+1>k2[1]:
                                flag+=1
                        if not flag==0:
                            print >>fo, ' '.join([str(i) for i in [k1]+k2])
            fo.close()
    elif function_name=='Size-Control':
        if len(sys.argv)<3:
            print 'SV.Output.Process.py Size-Control <parameters>'
            print ' '
            print 'Parameters for Mappable-Control:'
            print '--input: input.bed'
            print '--min-size: reference.Mappable. describing mappable regions'
            print '--max-size: reference.Mappable. describing mappable regions'
            print ' '
        else:
            opts,args=getopt.getopt(sys.argv[2:],'input:',['input=','min-size=','max-size='])
            dict_opts=dict(opts)
            min=int(dict_opts['--min-size'])
            max=int(dict_opts['--max-size'])
            filein=dict_opts['--input']
            fin=open(filein)
            fo=open(filein.replace('.bed','.min'+str(min)+'.max'+str(max)+'.bed'),'w')
            for line in fin:
                pin=line.strip().split()
                if int(pin[2])-int(pin[1])>min and int(pin[2])-int(pin[1])<max:
                    print >>fo, ' '.join(pin)
            fin.close()
            fo.close()
    elif function_name=='TRA-Control':
        if len(sys.argv)<3:
            print 'SV.Output.Process.py TRA-Control <parameters>'
            print ' '
            print 'Parameters for TRA-Control:'
            print '--input: input.bed'
            print '--TRA-rec: TRA information kept in .rec files'
            print ' '
        else:
            opts,args=getopt.getopt(sys.argv[2:],'i:',['TRA-rec=','input='])
            dict_opts=dict(opts)
            ref_hash={}
            fin=open(dict_opts['--TRA-rec'])
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in ref_hash.keys():
                    ref_hash[pin[0]]={}
                if not int(pin[1]) in ref_hash[pin[0]].keys():
                    ref_hash[pin[0]][int(pin[1])]=[]
                if not [pin[0]]+[int(i) for i in pin[1:4]] in ref_hash[pin[0]][int(pin[1])]:
                    ref_hash[pin[0]][int(pin[1])].append([pin[0]]+[int(i) for i in pin[1:4]])
            fin.close()
            in_hash={}
            chromos=[]
            fin=open(dict_opts['--input'])
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in in_hash.keys():
                    in_hash[pin[0]]={}
                    chromos.append(pin[0])
                if not int(pin[1]) in in_hash[pin[0]].keys():
                    in_hash[pin[0]][int(pin[1])]=[]
                in_hash[pin[0]][int(pin[1])].append([pin[0]]+[int(i) for i in pin[1:-1]])
            fin.close()
            def hash_to_list(in_hash):
                ref_list={}
                for k1 in in_hash.keys():
                    ref_list[k1]=[]
                    for k2 in sorted(in_hash[k1].keys()):
                        for k3 in in_hash[k1][k2]:
                            ref_list[k1].append(k3)
                return ref_list
            def bp_check(list1,list2):
                #eg: list1=['1', 3615657, 3638623]; list2=['1', 3615658, 3638624, 3751340];
                flag=0
                rec=[]
                for x in list1[1:]:
                    for y in list2[1:]:
                        if abs(x-y)<51:
                            flag+=1
                            rec.append(list2.index(y))
                return [flag,rec]
            def compare_list(in_list,ref_list):
                out={}
                out2=[]
                for k1 in in_list.keys():
                    if k1 in ref_list.keys():
                        out[k1]={}
                        for x1 in ref_list[k1]:
                            for x2 in in_list[k1]:
                                if x2[-1]<x1[1]: continue
                                elif x2[1]>x1[-1]:continue
                                else:
                                    bp_result=bp_check(x2,x1)
                                    if bp_result[0]==2:
                                        if bp_result[1]==[1,2]:
                                            block='a'
                                        elif bp_result[1]==[2,3]:
                                            block='b'
                                        elif bp_result[1]==[1,3]:
                                            block='ab'
                                        if not '_'.join([str(i) for i in x1]) in out[k1].keys():
                                            out[k1]['_'.join([str(i) for i in x1])]=[]
                                        out[k1]['_'.join([str(i) for i in x1])].append(x2+[block])
                                        out2.append(x2)
                return [out,out2]       
            in_list=hash_to_list(in_hash)
            ref_list=hash_to_list(ref_hash)
            comp_all=compare_list(in_list,ref_list)
            comp_hash=comp_all[0]
            comp_list=comp_all[1]
            fo=open(dict_opts['--input'].replace('.bed','.TRAFree.bed'),'w')
            for k1 in chromos:
                for k2 in in_list[k1]:
                    if not k2 in comp_list:
                        print >>fo, ' '.join([str(i) for i in k2])
            fo.close()
            #fo2=open(dict_opts['--input'].replace('.bed','.TRARelated.bed'),'w')
            #for k1 in chromos:
            #    for k2 in comp_hash[k1].keys():
            #        for k3 in comp_hash[k1][k2]:
            #            print >>fo2, ' '.join([str(i) for i in [k2]+k3])
            #fo2.close()

