import os
import re
import sys
import numpy
import random
import getopt
import math
import pickle
from math import sqrt,pi,exp
import scipy
from scipy.stats import norm
import time
import datetime
import itertools
import random
from scipy.stats import norm
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
def Af_TB_Penal_Less_Important_Caldu(num_of_reads,Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
    Af_TB_Rec=0
    Af_TB_Rec_list=Af_Info_all[-1]
    Af_TB_Penal_a=Af_Info_all[4]
    for x in Af_TB_Rec_list:
        if x <(IL_Statistics[0]*IL_Statistics[4]+IL_Statistics[1]*IL_Statistics[5])/ReadLength*GC_Overall_Median_Num/float(10):
            Af_TB_Rec+=1
    return -float(Af_TB_Penal_a)/float(num_of_reads)-float(Af_TB_Rec)/float(len(Af_Letter[0]+Af_Letter[1])+2)
def Af_TB_Penal_More_Important_Caldu(num_of_reads,Af_Info_all,Af_Letter,IL_Statistics,ReadLength,GC_Overall_Median_Num):
    Af_TB_Penal=Af_Info_all[3]
    Af_TB_Penal_a=Af_Info_all[4]
    return -float(Af_TB_Penal_a)/float(num_of_reads)-Af_TB_Penal
def bed_hash_short(hash_in,Chromo_Length):
    out={}
    for k1 in hash_in.keys():
        if hash_in[k1]==[[]]:
            if k1 in Chromo_Length.keys():
                out[k1]=[[0,Chromo_Length[k1]]]
        else:
            temp_hash={}
            for k2 in hash_in[k1]:
                if len(k2)>1:
                    if not int(k2[0]) in temp_hash.keys():
                        temp_hash[int(k2[0])]={}
                    if not int(k2[1]) in temp_hash[int(k2[0])].keys():
                        temp_hash[int(k2[0])][int(k2[1])]=[]
            temp_list=[]
            for k2 in sorted(temp_hash.keys()):
                for k3 in sorted(temp_hash[k2].keys()):
                    temp_list.append([k2,k3])
            temp_list2=[temp_list[0]]
            for k2 in temp_list[1:]:
                if not k2[0]-1000>temp_list2[-1][-1]:
                    temp_list2[-1][-1]=k2[1]
                else:
                    temp_list2.append(k2)
            out[k1]=temp_list2
    return out            
def Best_Let_modify(original_letters,Best_Letter_Rec,Best_Score_Rec,Score_rec_hash):
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
def Block_Assign_To_Letters(bp_list,letter_list,flank,Window_Size):
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
    #Eg:Ref_Sequence_Path_File='/home/xuefzhao/References/chr2:184569179-184572016.fa'
    #This function produce a hash list that contains the original sequences
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
    #Eg:Ref_Sequence_Path_File='/home/xuefzhao/References/chr2:184569179-184572016.fa'
    #This function produce a hash list that contains the original sequences
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
def bps4_to_bps2(bps4):
    bps2=[]
    for k1 in bps4.keys():
        for k2 in bps4[k1]:
            bps2.append(k2)
    return bps2
def Indicator_Readin(temp_inter):
    ft=open(temp_inter)
    F_test='ERROR'
    for line in ft:
        pin=line.strip().split()
        if pin[0]=='SVelter':
            F_test=pin[1]
    ft.close()
    return F_test
def bed_readin(input_bed):
    data_hash={}
    fin=open(input_bed)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in data_hash.keys():
            data_hash[pin[0]]=[]
        data_hash[pin[0]].append([int(pin[1]),int(pin[2])])
    fin.close()
    return data_hash
def bed_write(bed_info,out_folder,sample_name,input_bed):
    for k1 in bed_info.keys():
        fout=out_folder+'.'.join(sample_name.split('.')[:-1])+'.'+k1+'.predefinedBP.'+'.'.join(input_bed.split('/')[-1].split('.')[:-1])+'.LNs'
        fo=open(fout,'w')
        for k2 in bed_info[k1]:
            print >>fo, ' '.join([str(i) for i in [k1,k2[0],100,k2[1],100]])
        fo.close()
        fo=open(out_folder+'.'.join(sample_name.split('.')[:-1])+'.'+k1+'.predefinedBP.'+'.'.join(input_bed.split('/')[-1].split('.')[:-1])+'.SPs','w')
        fo.close()
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
def IL_CI_Decide(ILStats,Times_of_stds,model_comp):
    if model_comp=='C':
        return CI_Calculate_Bimodal(ILStats,Times_of_stds)
    if model_comp=='S':
        return CI_Calculate_Normal(ILStats,Times_of_stds)
def IL_Stat_Simp(StatFile):
    #for normal distribution used as IL null model 
    fstat=open(StatFile)
    temp=fstat.readline()
    temp=fstat.readline()
    temp=fstat.readline()
    p1=fstat.readline().strip().split()
    Mean1=float(p1[1])
    STD1=float(p1[2])
    fstat.close()
    return [[Mean1,0,STD1,1,1,0],[1,Mean1,STD1]]
def IL_Stat_readin(StatFile):  
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
def RD_Stat_readin(StatFile):
    #eg of StatFile='/scratch/remills_flux/xuefzhao/SV_discovery_index/download/NullModel.NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram/RDNull.NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.genome.NegativeBinomial'
    fstat=open(StatFile)
    temp=fstat.readline()
    p1=fstat.readline().strip().split()
    fstat.close()
    return [float(i) for i in p1]
def CI_Calculate_Bimodal(ILStats,Times_of_stds):
    alpha1=float(ILStats['bimodal']['Bimodal1'])
    miu1=float(ILStats['bimodal']['Mean1'])
    std1=float(ILStats['bimodal']['STD1'])
    alpha2=float(ILStats['bimodal']['Bimodal2'])
    miu2=float(ILStats['bimodal']['Mean2'])
    std2=float(ILStats['bimodal']['STD2'])
    mean=alpha1*miu1+alpha2*miu2
    theta1=float(ILStats['bimodal']['Mean1'])-mean
    theta2=float(ILStats['bimodal']['Mean2'])-mean
    var=alpha1*(std1**2+theta1**2)+alpha2*(std2**2+theta2**2)
    std=numpy.sqrt(var)
    return [max([0,mean-Times_of_stds*std]), mean+Times_of_stds*std]
def CI_Calculate_Normal(ILStats,Times_of_stds):
    le_CI=float(ILStats['normal']['Mean'])-Times_of_stds*float(ILStats['normal']['STD'])
    ri_CI=float(ILStats['normal']['Mean'])+Times_of_stds*float(ILStats['normal']['STD'])
    return [max([0,le_CI]),ri_CI]
def ClusterLen_Calculation(ILStats,model_comp,ReadLength):
    if model_comp=='C':
        return float(ILStats['normal']['Mean'])+2*float(ILStats['normal']['STD'])-ReadLength
    else:
        alpha1=float(ILStats['bimodal']['Bimodal1'])
        miu1=float(ILStats['bimodal']['Mean1'])
        std1=float(ILStats['bimodal']['STD1'])
        alpha2=float(ILStats['bimodal']['Bimodal2'])
        miu2=float(ILStats['bimodal']['Mean2'])
        std2=float(ILStats['bimodal']['STD2'])
        mean=alpha1*miu1+alpha2*miu2
        theta1=float(ILStats['bimodal']['Mean1'])-mean
        theta2=float(ILStats['bimodal']['Mean2'])-mean
        var=alpha1*(std1**2+theta1**2)+alpha2*(std2**2+theta2**2)
        std=numpy.sqrt(var)
        return mean+2*std-ReadLength
def Compare_LN_list_against_Segdup_list(segdup_list,LN_list,chromo_name):
    in_list=[]
    out_list=[]
    rec1=0
    rec2=0
    while True:
        if not chromo_name in segdup_list.keys(): break
        if rec2==len(segdup_list[chromo_name]) or rec1==len(LN_list): break
        if LN_list[rec1][1]<segdup_list[chromo_name][rec2][0]:
            in_list.append(rec1)
            rec1+=1
        elif LN_list[rec1][0]<segdup_list[chromo_name][rec2][1]+1 and LN_list[rec1][0]>segdup_list[chromo_name][rec2][0]-1:
            out_list.append(rec1)
            rec1+=1
        elif LN_list[rec1][1]<segdup_list[chromo_name][rec2][1]+1 and LN_list[rec1][1]>segdup_list[chromo_name][rec2][0]-1:
            out_list.append(rec1)
            rec1+=1
        elif LN_list[rec1][0]>segdup_list[chromo_name][rec2][1]:
            rec2+=1
        elif LN_list[rec1][0]<segdup_list[chromo_name][rec2][0] and LN_list[rec1][1]>segdup_list[chromo_name][rec2][1]:
            in_list.append(rec1)
            rec1+=1
    if rec1<len(LN_list):
        in_list+=range(len(LN_list))[rec1+1:]
    return [LN_list[i] for i in in_list if i<len(LN_list)]
def Compare_SP_list_against_Segdup_list(segdup_list,SP_list,chromo_name):
    in_list=[]
    out_list=[]
    rec1=0
    rec2=0
    while True:
        if not chromo_name in segdup_list.keys(): break
        if rec1==len(SP_list) or rec2==len(segdup_list[chromo_name]): break
        if SP_list[rec1][0]<segdup_list[chromo_name][rec2][0]:
            in_list.append(rec1)
            rec1+=1
        elif SP_list[rec1][0]>segdup_list[chromo_name][rec2][0]-1 and SP_list[rec1][0]<segdup_list[chromo_name][rec2][1]+1:
            out_list.append(rec1)
            rec1+=1
        elif SP_list[rec1][0]>segdup_list[chromo_name][rec2][1]:
            rec2+=1
    return [SP_list[i] for i in in_list if i<len(SP_list)]
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
def calculate_interval_region(temp3_list,Chromo_Length,chrom):
    out=[] 
    if len(temp3_list)==1:
        if temp3_list[0][0]==0 and temp3_list[0][1]==Chromo_Length[chrom]:
            return out
    if not temp3_list[0][0]==0:
        out.append([0,temp3_list[0][0]-1])
    for x in temp3_list[:-1]:
        out.append([x[1]+1,temp3_list[temp3_list.index(x)+1][0]-1])
    if temp3_list[-1][1]<Chromo_Length[chrom]:
        out.append([temp3_list[-1][1]+1,Chromo_Length[chrom]-1])
    return out
def calculate_len_genome(ref):
    whole_genome=chromos_read_in(ref)
    len_genome=0
    for i in whole_genome.keys():
        len_genome+=whole_genome[i][0]
    return [whole_genome,len_genome]

def cdf_solver_application(Insert_Len_Stat,cdf,model_comp):
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
def chromos_info_readin(ref_index):
    chromos=[]
    Chromo_Length={}
    fin=open(ref_index)
    for line in fin:
        pin=line.strip().split()
        chromos.append(pin[0])
        Chromo_Length[pin[0]]=int(pin[1])
    fin.close()
    return [chromos,Chromo_Length]
def chromos_readin_list(ref):
    fin=open(ref+'.fai')
    chromos=[]
    for line in fin:
            pin=line.strip().split()
            chromos.append(pin[0])
    fin.close()
    return chromos
def chromos_read_in(ref):
    whole_genome={}
    fref=open(ref+'.fai')
    for line in fref:
        pref=line.strip().split()
        whole_genome[pref[0]]=[int(pref[1])]
    fref.close()
    return whole_genome
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
def cigar2cigars(cigar):
    pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
    cigars=[]
    for m in pcigar.finditer(cigar):
        cigars.append((m.groups()[0],m.groups()[1]))
    return cigars
def cigar2split(cigar):
    cigars=cigar2cigars(cigar)
    MapLen=[]
    maplen=0
    if cigars[0][1] in ['S','H']:
        MapLen.append(0)
        for n in cigars[1:]:
            if n[1]=='M' or n[1]=='D' or n[1]=='N':
                maplen+=int(n[0])
            if n[1] in ['S','H']: 
                MapLen.append(maplen-1)
    else:
        for n in cigars:
            if n[1]=='M' or n[1]=='D' or n[1]=='N':
                maplen+=int(n[0])
            if n[1] in ['S','H']: 
               MapLen.append(maplen-1)         
    return MapLen
def cigar2splitlength(cigar):
    cigars=cigar2cigars(cigar)
    splitlen=[]
    for i in cigars:
        if i[1] in ['S','H']:
            splitlen.append(int(i[0]))
    return splitlen
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
def clean_svelter_set(path):
    if os.path.isdir(path):
        os.system(r'''rm -r %s'''%(path))
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
        return [out2,out3]
def clusterNums_svpredict(data, ClusterLen, direction):
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
    data1=clusterNums(hash.keys(),ClusterLen,direction)[0]
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
    #eg: cluSP= clusterNumSP(SP,SPCluMin,SplitMin)
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
    #eg: hash= KeyF3
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
def clusterNums_bpintegrate(Data1, ClusterLen, direction):
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
def clusterQC_bpintegrate(hash, CluNum):
    out=[]
    for i in range(len(hash[1])):
        if not hash[1][i]<CluNum:
            out.append(hash[0][i])
    return out
def clusterSupVis_bpintegrate(DataRou1,DataCen1,ClusterLen):
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
def cn2_file_read_in(dict_opts,workdir):
    if '--copyneutral' in dict_opts.keys():
        cn2_file=dict_opts['--copyneutral']
    else:
        cn2_file=workdir+'reference_SVelter/CN2.bed'
    return cn2_file
def cn2_length_readin(dict_opts):
    if '--null-copyneutral-length' in dict_opts.keys():
        cn2_length=int(dict_opts['--null-copyneutral-length'])
    else:
        cn2_length=2000
    return cn2_length
def cn2_length_read_in(dict_opts):
    if '--null-copyneutral-length' in dict_opts.keys():
        cn2_length=int(dict_opts['--null-copyneutral-length'])
    else:
        cn2_length=2000
    return cn2_length
def cn2_region_write(cn2_file,dict_opts,ref):
    if not '--null-random-length' in dict_opts.keys():
        dict_opts['--null-random-length']=5000
    else:
        dict_opts['--null-random-length']=int(dict_opts['--null-random-length'])
    if not '--null-random-num' in dict_opts.keys():
        dict_opts['--null-random-num']=10000
    else:
        dict_opts['--null-random-num']=int(dict_opts['--null-random-num'])
    cn2_length=dict_opts['--null-random-length']-100
    [whole_genome,len_genome]=calculate_len_genome(ref)
    fo=open(cn2_file,'w')
    for i in whole_genome.keys():
        num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--null-random-num'])
        reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--null-random-length']) for j in range(num_i)]
        for j in sorted(reg_i):
            print >>fo, ' '.join([i,str(j),str(j+dict_opts['--null-random-length']-1)])
    fo.close()
    SamplingPercentage=1
    return [cn2_length,SamplingPercentage,whole_genome,len_genome]            
def complement(direction):
    if direction=='+':
        return '-'
    elif direction=='-':
        return '+'
    else:
        return 'error'
def complementary(seq):
    seq2=[]
    for i in seq:
        if i in 'ATGCN':
            seq2.append('ATGCN'['TACGN'.index(i)])
        elif i in 'atgcn':
            seq2.append('atgcn'['tacgn'.index(i)])
    return ''.join(seq2)        
def Direction_Penalize_2a(Coordinate_Info):     
    #This function is used to calculate the Direction_Penalize for Letter_Double only
    #Eg of Coordinate_Info: hash containing many: 'HWI-ST177_136:2:47:18313:117413': ['1032', '1132', '1407', '1507', '+', '-'], 
    Direction_Penalty=0
    for key in Coordinate_Info.keys():
        LETT=Coordinate_Info[key]
        for lett in LETT:
            if lett[4]==lett[5]:
                Direction_Penalty+=1
            else: continue
    return Direction_Penalty
def Direction_Penalize_2b(Coordinate_Info):     
    #This function is used to calculate the Direction_Penalize for Letter_Double only
    #Eg of Coordinate_Info: hash containing many: 'HWI-ST177_136:2:47:18313:117413': ['1032', '1132', '1407', '1507', '+', '-'], 
    Direction_Penalty=0
    for key in Coordinate_Info.keys():
        lett=Coordinate_Info[key]
        if lett[4]==lett[5]:
            Direction_Penalty+=1
        else: continue
    return Direction_Penalty
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
def ex_file_read_in(dict_opts,workdir):
    if '--exclude' in dict_opts.keys():
        ex_file=dict_opts['--exclude']
    else:
        ex_file=workdir+'reference_SVelter/Exclude.bed'
    return ex_file
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
def Fasta_To_Sequence_nullmodel(fasta_file):
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
def clusterQC(hash, CluNum):
    out=[]
    for i in range(len(hash[1])):
        if not hash[1][i]<CluNum:
            out.append(hash[0][i])
    return out
def createPwd(sList,kList):
    assert len(sList)*len(kList)>0
    upper,sList=len(sList),sList*len(kList)
    for i in range(len(kList)):
        for j in range(upper):
            sList[i*upper+j]+=kList[i]
    return sList
def getPwd(length):
    s=''.join([chr(97+i) for i in range(length)])
    sSet=set(s)
    if len(sSet)<length:return 0
    result=tmp=list(sSet)
    for i in range(length-1):
        result=createPwd(result,tmp)
    return result
def path_modify(path):
    if not path[-1]=='/':
        path+='/'
    return path
def path_mkdir(path):
    if not os.path.isdir(path):
            os.system(r'''mkdir %s'''%(path))
def random_produce_cn2_region(cn2_file,whole_genome,len_genome,dict_opts):
    cn2_path='/'.join(cn2_file.split('/')[:-1])+'/'
    path_mkdir(cn2_path)
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
    return [SamplingPercentage,cn2_length]
def random_produce_exclude_region(ex_file,chromos):
    fo=open(ex_file,'w')
    for chr_ex in chromos:
        print >>fo, ' '.join([chr_ex,'0','0'])
    fo.close()
def random_pick_cn2_region(cn2_file,whole_genome,chromos,len_genome,dict_opts):
    out=[]
    if not '--null-random-length' in dict_opts.keys():
        dict_opts['--null-random-length']=5000
    else:
        dict_opts['--null-random-length']=int(dict_opts['--null-random-length'])
    if not '--null-random-num' in dict_opts.keys():
        dict_opts['--null-random-num']=10000
    else:
        dict_opts['--null-random-num']=int(dict_opts['--null-random-num'])
    cn2_length=dict_opts['--null-random-length']-100
    for i in sorted(whole_genome.keys()):
        if i in chromos:
            num_i=int(float(whole_genome[i][0])/float(len_genome)*dict_opts['--null-random-num'])
            reg_i=[random.randint(1,whole_genome[i][0]-dict_opts['--null-random-length']) for j in range(num_i)]
            for j in sorted(reg_i):
                 out.append([i,str(j),str(j+dict_opts['--null-random-length']-1)])
    SamplingPercentage=1
    return out
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
def file_setup(file):
    fin=open(file,'w')
    fin.close()
def file_initiate(file_name):
    if not os.path.isfile(file_name):
        fo=open(file_name,'w')
        fo.close()
def Gap_Hash_Initiate(chromos):
    Gap_Hash={}
    for chr_ex in chromos:
        Gap_Hash[chr_ex]=[]
    return Gap_Hash
def GC_Content_Calculate_num(seq2):
    NumAT=0
    NumGC=0
    for n in seq2:
        if n=='A' or n=='T' or n=='a' or n=='t':
            NumAT+=1
        elif n=='G' or n=='C' or n=='g' or n=='c':
            NumGC+=1
    return [NumGC,NumAT]            
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
def GC_Stat_ReadIn(BamN,GC_Stat_Path,genome_name,affix):
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
    for k1 in CN2_Region.keys():
        for k2 in CN2_Region[k1].keys():
            if CN2_Region[k1][k2][0]=='@:':
                del CN2_Region[k1][k2]
            else:
                if sum([float(i) for i in CN2_Region[k1][k2][0].split(':')[1].split(',')])==0.0:
                    del CN2_Region[k1][k2]
        if CN2_Region[k1]=={}:
            del CN2_Region[k1]
    return [CN2_Region,Chromos,GC_Content]
def GC_index(ref_prefix,chrbam):
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
def genome_name_readin(dict_opts):
    if not '--NullGenomeName' in dict_opts.keys():
        genome_name='genome'
    else:
        genome_name=dict_opts['--NullGenomeName']
    return genome_name
def GC_Index_Readin(fileGC):
    fin=open(fileGC)
    out={}
    while True:
        pin=fin.readline().strip().split()
        if not pin: break
        pin2=fin.readline().strip().split()
        if not pin[0] in out.keys():
            out[pin[0]]={}
        if not int(pin[1]) in out[pin[0]].keys():
            out[pin[0]][int(pin[1])]={}
        if not int(pin[2]) in out[pin[0]][int(pin[1])].keys():
            out[pin[0]][int(pin[1])][int(pin[2])]=[]
        out[pin[0]][int(pin[1])][int(pin[2])]=pin2
    fin.close()
    return out
def GC_RD_Adj_hash(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
    Coverage_af_Adj={}
    Overall_Median_Coverage=float(GC_Overall_Median_Num)
    for key_1 in GC_Content.keys():
            Be_Adj_RD=Coverage[key_1]
            GC_Con=int(round(GC_Content[key_1]*100))
            if GC_Con in GC_Median_Num.keys():
                Median_Coverage=GC_Median_Num[GC_Con]
            else:
                Median_Coverage=Overall_Median_Coverage
            Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
            Coverage_af_Adj[key_1]=Af_Adj_RD
    return  Coverage_af_Adj
def GC_RD_Info_Complete(ref_file,GC_Median_Coverage,ChrN_Median_Coverage,GC_Overall_Median_Coverage,GC_Var_Coverage,GC_Mean_Coverage,GC_Std_Coverage,Chromosome):
    global refFlag
    refFlag=0
    refgenome=ref_file
    reftest=open(refgenome)
    preftest=reftest.readline().strip().split()
    if len(preftest[0])>3:
        refFlag=1
    reftest.close()
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
    return [chrom_N,chrom_X,chrom_Y,GC_Median_Coverage,GC_Overall_Median_Coverage,GC_Var_Coverage,GC_Mean_Coverage,GC_Std_Coverage]
def hash_Cor_modify(hash_Cor,Chromo_Length,hash_key,chrom):
    out=[]
    if hash_Cor==[[0]]:
        out=[[0,Chromo_Length[chrom]]]
    for hts in hash_Cor:
        if len(hts)==2:
            if hts[1]>Chromo_Length[hash_key]:
                    break
            else:
                out.append(hts)
        else:
            continue
    return out
def IL_Stat_Calculate(InsertLength):
    TotalILNum=0
    for key in InsertLength.keys():
            TotalILNum+=InsertLength[key]
    return TotalILNum
def ILStats_readin(ILStat):
    ILStats={}
    fILS=open(ILStat)
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
    fILS.close() 
    return ILStats
def IL_Prob_Bimodal(IL_Length,ILStats):
    log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std1']))-(float(IL_Length)-float(ILStats['bimodal']['Mean1']))**2/2/float(ILStats['bimodal']['std1'])**2
    log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['bimodal']['std2']))-(float(IL_Length)-float(ILStats['bimodal']['Mean2']))**2/2/float(ILStats['bimodal']['std2'])**2
    return log_y1*float(ILStats['bimodal']['Bimodal1'])+log_y2*float(ILStats['bimodal']['Bimodal2'])
def IL_Prob_Normal(IL_Length,ILStats):
    log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(ILStats['normal']['std']))-(float(IL_Length)-float(ILStats['normal']['Mean']))**2/2/float(ILStats['normal']['std'])**2
    return log_y
def IL_Penal_Calcu(read_info,IL_Statistics,Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero):
    Initial_IL=[]
    for j2 in read_info[2].keys():
        for j in read_info[2][j2]:
            Initial_IL.append(j[3]-j[1])
    for j2 in read_info[3].keys():
        for j in read_info[3][j2]:
            Initial_IL.append(j[3]-j[1])
    Initial_ILPenal=[]
    for j in Initial_IL:
        Initial_ILPenal+=[pdf_calculate(j,IL_Statistics[4],IL_Statistics[0],IL_Statistics[1],IL_Statistics[2],IL_Statistics[3],Cut_Upper,Cut_Lower,Penalty_For_InsertLengthZero)/len(Initial_IL)]
    return Initial_ILPenal
def tag_inv(letter_in):
    #eg of letter_in=[['a^'], []]
    out=0
    for x in letter_in:
        for y in x:
            if '^' in y:    out+=1
    return out
def insert_block_produce(Letter_List,Letter_List_origin):
    Insert_Pool=[]
    for i in Letter_List_origin:
        if not i in Letter_List and not i+'^' in Letter_List:
            Insert_Pool.append(i)
    return Insert_Pool
def insert_BPs_produce(Letter_List,Letter_List_origin):
    Insert_Pool=[]
    for i in Letter_List_origin:
        if not i in Letter_List and not i+'^' in Letter_List:
            Insert_Pool.append(i)   
    Insert_Pool2=[[(ord(j)-96),(ord(j)-95)] for j in Insert_Pool]
    return Insert_Pool2
def invert_block_produce(Letter_List):
    invert_block=[x.upper() for x in Letter_List]
    return invert_block
def invert_BPs_produce(Letter_List):
    invert_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
    return invert_BPs
def letter_range_report(flank,chr_letter_bp):
    ks_range={}
    for k1 in chr_letter_bp.keys():
        ks_range[k1]=[]
        for k2 in chr_letter_bp[k1].keys():
            ks_range[k1]+=chr_letter_bp[k1][k2]
        ks_range[k1]=[min(ks_range[k1]),max(ks_range[k1])]
    for k1 in chr_letter_bp.keys():
        chr_letter_bp[k1]['left']=[ks_range[k1][0]-flank,ks_range[k1][0]]
        chr_letter_bp[k1]['right']=[ks_range[k1][1],ks_range[k1][1]+flank]
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
def LN_Filter(LN_filein,SP_filein,workdir):
    filein=workdir+'reference_SVelter/Segdup.bed'
    if os.path.isfile(filein):
        chromo_name=LN_filein.split('.')[-6]
        segdup_list=Segdup_readin(filein,chromo_name)
        LN_list=LN_file_readin(LN_filein)
        LN_keep_list=Compare_LN_list_against_Segdup_list(segdup_list,LN_list,chromo_name)
        SP_list=SP_file_readin(SP_filein)
        SP_keep_list=Compare_SP_list_against_Segdup_list(segdup_list,SP_list,chromo_name)
        Write_filtered_LN_list(LN_keep_list,LN_filein,chromo_name)
        Write_filtered_SP_list(SP_keep_list,SP_filein,chromo_name)
def LN_file_readin(LN_filein):
    LN_hash={}
    fin=open(LN_filein)
    for line in fin:
        pin=line.strip().split()
        if not int(pin[1]) in LN_hash.keys():
            LN_hash[int(pin[1])]={}
        if not int(pin[3]) in LN_hash[int(pin[1])].keys():
            LN_hash[int(pin[1])][int(pin[3])]=[]
        LN_hash[int(pin[1])][int(pin[3])].append([pin[2],pin[4]])
    fin.close()
    LN_list=[]
    for k1 in sorted(LN_hash.keys()):
        for k2 in sorted(LN_hash[k1].keys()):
            for k3 in LN_hash[k1][k2]:
                LN_list.append([k1,k2]+k3)
    return LN_list
def LNd_hash_readin(chrom_LN,bps_in_path):
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
        clu1=clusterNums_bpintegrate(sorted([int(i) for i in LNd_hash[k1][0].keys()]), 1000, 'f')
        cen1=clusterQC_bpintegrate(clu1, 2)
        if not cen1==[]:
            sup1=clusterSupVis_bpintegrate(sorted([int(i) for i in LNd_hash[k1][0].keys()]),cen1,1000)
            for k3 in sup1.keys():
                temp3=[k1.split('_')[0]]+[str(i) for i in sorted(sup1[k3])]+[k1.split('_')[1]]
                temp4=[]
                for k4 in sup1[k3]:
                    temp4+=LNd_hash[k1][0][str(k4)]
                temp4.sort()
                temp3+=temp4
                if not temp3 in temp2:
                    temp2.append(temp3)
        clu2=clusterNums_bpintegrate(sorted([int(i) for i in LNd_hash[k1][1].keys()]), 1000, 'f')
        cen2=clusterQC_bpintegrate(clu2, 2)
        if not cen2==[]:
            sup2=clusterSupVis_bpintegrate(sorted([int(i) for i in LNd_hash[k1][1].keys()]),cen2,1000)
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
def LNe_hash_Filter(LNe_hash,allchromos):
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
def LNd_hash_pair(LNd_hash,SP_Info):
    out={}
    for k1 in LNd_hash.keys():
        if k1.split('_')[0] in SP_Info.keys() and k1.split('_')[1] in SP_Info.keys():
            out[k1]=[]
            cen1=[int(i) for i in LNd_hash[k1][0].keys()]
            if not cen1==[]:
                clu1=clusterSupVis_bpintegrate(sorted(SP_Info[k1.split('_')[0]].keys()),cen1,10000)
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
                        clu2=clusterSupVis_bpintegrate(sorted(SP_Info[k1.split('_')[0]].keys()),k3,10000)
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
def LN_Info_Correct(LN_hash,unique_SPs):
    #eg of LN_hash: LN_hash=SP_LN_Info[0]
    out=[]
    for k1 in LN_hash:
        out.append([])
        for k2 in k1:
            if not k2 in unique_SPs:
                temp=[i for i in unique_SPs if i<k2+200 and i>k2-200]
                if temp==[]:
                    out[-1].append(k2)
                else:
                    diff=[abs(i-k2) for i in temp]
                    new_k2=temp[diff.index(min(diff))]
                    out[-1].append(new_k2)
            else:
                out[-1].append(k2)
    return out
def unify_list(list):
    out=[]
    for x in list:
        if not x in out:
            out.append(x)
    return out
def order_int_list(list):
    out1=[int(i) for i in list]
    out1.sort()
    out2=[]
    for x in out1:
        if not x in out2:
            out2.append(x)
    return out2
def LN_bps_write(bps_hash,bps_folder,S_Sample,dict_opts,chromos,allchromos,bps_in_path):
    SP_LI_in=SP_info_ReadIn_pre(bps_hash,S_Sample,bps_in_path)
    SP_links=SP_LI_in[0]
    SP_Info=SP_LI_in[1]
    SP_Check=SP_LI_in[2]
    S_Chr=bps_hash[S_Sample].keys()[0] 
    for k1 in os.listdir(bps_in_path):
        if k1.split('.')[-1]=='chromLNs':
            chrom_LN=k1
            LNall_hash=LNd_hash_readin(chrom_LN,bps_in_path)
            LNd_hash=LNall_hash[0]
            LNe_hash=LNe_hash_Filter(LNall_hash[1],allchromos)
            LNd_hash=LNd_hash_pair(LNd_hash,SP_Info)
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
                fout=bps_folder+S_Sample+'.txt'
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
                    fout=bps_folder+S_Sample+'.'+'LN'+'.'+'txt'
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
                    LN_out3=[[[i for i in j if len(i)>1]for j in k]for k in LN_out2]
                    LN_out4=[[j for j in k if len(j)>0]for k in LN_out3] 
                    LN_out2=[] 
                    for x in LN_out4:
                        LN_out2.append([])
                        for y in x:
                            test={}
                            for z in y:
                                test_flag=0
                                for u in z[1:]:
                                    if not u.isdigit():
                                        test_flag+=1
                                if test_flag>0:continue
                                if not z[0] in test.keys():
                                    test[z[0]]=[]
                                test[z[0]]+=z[1:]
                            for z in test.keys():
                                test[z]=order_int_list(test[z])
                            if len(test)>0:
                                LN_out2[-1].append([])
                                for z in test.keys():
                                    LN_out2[-1][-1].append([z]+[str(i) for i in test[z]])
                    rec=0
                    for k1 in LN_out2:
                        rec+=1
                        fout=bps_folder+S_Sample+'.LN.'+str(rec)+'.'+'txt'
                        file_setup(fout)
                        fo=open(fout,'a')
                        for k2 in k1:
                            for k3 in k2:
                                print >>fo, ' '.join(k3)
                            print >>fo, ' '
                        fo.close()
def LN_bps2_Modify(bps2,chromos_all):
    out=[]
    for x in bps2:
        flag_chr_QC=x_qc_check(x,chromos_all)
        if flag_chr_QC>0:
            for y in x:
                if y in chromos_all:
                    out.append([y])
                else:
                    out[-1].append(y)
        else:
            print 'Error: chromosome not found in reference !'
    out2=[]
    for x in out:
        if len(x)>1:
            out2.append(x)
    return out2
def LN_Merge_Final_Check(LN_LN_Merge):
    out=[LN_LN_Merge[0]]
    for x in LN_LN_Merge[1:]:
        if x[0]>out[-1][-1]:
            out.append(x)
        else:
            if max(out[-1])-min(out[-1])<100000:
                out[-1]+=x
            else:
                out.append(x)
        out[-1].sort()
    out2=[]
    for x in out:
        out2.append(list_unify(x))
    return out2
def modify_bps2_new(bps2_new):
    if len(bps2_new)>0:
        out1={}
        for k2 in bps2_new:
            if not k2[0] in out1.keys():
                out1[k2[0]]=[]
            out1[k2[0]]+=k2[1:]
        out=[]
        for x in out1.keys():
            out_temp=sorted([int(i) for i in out1[x]])
            out.append([x]+[str(i) for i in out_temp])
        return out
    else:
        return bps2_new
def list_unify(list):
    out=[]
    list.sort()
    for x in list:
        if not x in out:
            out.append(x)
    return out
def log_factorial(K):
    k=int(K)
    k_list=[0]
    for i in range(k+1)[2:]:
        k_list.append(k_list[i-2]+numpy.log(i))
    return k_list[-1]
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
def left_RD_Calculate_2a(Through_GCRD_Adj,Af_GCRD_Adj,flank,Window_Size):
    left_blocks=range(flank/Window_Size+1)[1:]
    left_Changes=[Af_GCRD_Adj[j][3]-Through_GCRD_Adj[j][3] for j in left_blocks]
    return numpy.mean(left_Changes)
def length_control(out3):
    #control length of each cluster
    out=[]
    for x in out3:
        if sum([len(y) for y in x])==1: continue
        else:
            if sum([len(y) for y in x])>10:
                for y in x:
                    if len(y)>1:
                        out.append(y)
            else:
                out.append([])
                for y in x:
                    out[-1]+=y
    return out
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
def link_SP_link3(SP_link3,SP_Info,min_length):
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
def Median_Pick(number_list):
    #number_list: a list of number
    #Output: median of all numbers in the list
    if 2*(len(number_list)/2)+1==len(number_list):
        return float(sorted(number_list)[len(number_list)/2])
    elif 2*(len(number_list)/2)==len(number_list):
        return float(sorted(number_list)[len(number_list)/2-1]+sorted(number_list)[len(number_list)/2])/float(2)
def Move_Choose(Move_Sample_Pool,Ploidy,Move_probs):
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
    #ChrAllele Eg: 2m, 2p
    #Only delete, invert, insert
    from random import choice
    #Choice_Pool=['delete','invert','CopyPaste','CutPaste','insert']
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
def Move_Decide_deterministic(IL_List,RD_List,GC_Var_Coverage):
    IL_Weight=1
    RD_Weight=1
    #regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
    regulator=1
    T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
    T2_Penal=[math.exp(k) for k in T_Penal]
    Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
    output=Normalized_Penal.index(max(Normalized_Penal))
    return [output,IL_List[output]*IL_Weight+RD_List[output]*RD_Weight]
def Move_Decide_2(IL_List,RD_List,GC_Var_Coverage):
    #here, we add direction as an extra part to IL penalty
    IL_Weight=1
    RD_Weight=1
    #regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
    regulator=1
    T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
    T2_Penal=[math.exp(k) for k in T_Penal]
    if sum(T2_Penal)>0:
        Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
        indicator=float(random.choice(range(1000)))/1000
        for j in range(len(IL_List)):
            cdf=numpy.sum(Normalized_Penal[:j+1])
            if cdf>indicator or cdf==indicator:
                return [j,IL_List[j]*IL_Weight+RD_List[j]*RD_Weight]
    return ''
def Move_Decide_3(IL_List,RD_List,GC_Var_Coverage):
    #pick the highest scored structure as the final
    IL_Weight=1
    RD_Weight=1
    #regulator=-numpy.max([(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight) for j in range(len(IL_List))])/5
    regulator=1
    T_Penal=[(IL_List[j]*IL_Weight+RD_List[j]*RD_Weight)/regulator for j in range(len(IL_List))]
    T_Penal_modify=[]
    T_Penal_Rec=[]
    rec=-1
    for x in T_Penal:
        rec+=1
        if x<100:
            T_Penal_modify.append(x)
            T_Penal_Rec.append(rec)
    T2_Penal=[math.exp(k) for k in T_Penal_modify]
    Normalized_Penal=[l/numpy.sum(T2_Penal) for l in T2_Penal]
    j=Normalized_Penal.index(max(Normalized_Penal))
    k=T_Penal_Rec[j]
    return [k,IL_List[k]*IL_Weight+RD_List[k]*RD_Weight]
def merge_SPs_into_LNs(unique_SPs,multi_removed_LNs):
    unique_SPs.sort()
    rec1=0
    rec2=0
    out=[]
    while True:
        if rec1==len(unique_SPs):
            if rec2<len(multi_removed_LNs):
                out.append(multi_removed_LNs[rec2+1:])
                break
            else:
                break
        elif rec2==len(multi_removed_LNs):
            if rec1<len(unique_SPs):
                out.append(unique_SPs[rec1+1:])
                break
            else:break
        elif multi_removed_LNs[rec2][1]-multi_removed_LNs[rec2][0]>10**5:
            out.append(multi_removed_LNs[rec2])
            rec2+=1
        elif unique_SPs[rec1]>multi_removed_LNs[rec2][0]-1 and unique_SPs[rec1]<multi_removed_LNs[rec2][1]+1:
            multi_removed_LNs[rec2].append([unique_SPs[rec1]])
            rec1+=1
        elif unique_SPs[rec1]<multi_removed_LNs[rec2][0]:
            out.append([unique_SPs[rec1]])
            rec1+=1
        elif unique_SPs[rec1]>multi_removed_LNs[rec2][1]:
            out.append(multi_removed_LNs[rec2])
            rec2+=1
    out2=[unify_list(i) for i in out]
    out3=[i for i in out2 if not len(i)==1]
    return out3
def merge_LNs_overlap_unit(out2):
    #merge LNs  units that overlap with eacg other
    out3=[]
    for k1 in out2:
        if max(k1)-min(k1)<10**6:
            if out3==[]:
                out3.append([k1])
            if k1[0]<max([max(i) for i in out3[-1]]):
                out3[-1].append(k1)
            else:
                out3.append([k1])
        else:
            out3.append([k1])
    return out3            
def merge_LNs_into_LNs(multi_removed_LNs):
    out=[]
    for k1 in multi_removed_LNs:
        if out==[]:
            out.append(k1)
        if k1[0] in out[-1] or k1[1] in out[-1]:
            out[-1]+=k1
        else:
            out.append(k1)
    out2=[unify_list(i) for i in out]
    out3=merge_LNs_overlap_unit(out2)
    out4=length_control(out3)
    return out4
def multi_trans_detect(modified_LNs):
    #remove bps with >10 mate bps, which might be ploy A/ T tails of mobile elements
    number,counts=numpy.unique(modified_LNs,return_counts=True)
    remove=number[counts>10]
    out=[]
    for x in modified_LNs:
        if not x[0] in remove and not x[1] in remove:
            out.append(x)
    return out
def multi_dup_check(k2,Multi_Dup):
    #decide if k2 is a multi_dup block
    flag=0
    for x in k2:
        if x in Multi_Dup:
            flag+=1
    return flag
def multi_dup_define(letter_RD,GC_Overall_Median_Coverage):
    #pick out blocks with RD > 4* mean RD;
    #such blocks are worthless to rearrange
    Multi_Dup=[]
    for x in letter_RD.keys():
        for y in letter_RD[x].keys():
            if letter_RD[x][y]>GC_Overall_Median_Coverage[x]*4:
                Multi_Dup.append(y)
    return Multi_Dup
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
        data_cen=clusterNums_bpintegrate(temp6, 10**6, 'r')[0]
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
def NullILCI_Calculate(InsertLength,TotalILNum,NullILCIs):
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
def NullPath_SetUp(out_path,dict_opts):
    if '--out-path' in dict_opts.keys():
        NullPath=path_modify(dict_opts['--out-path'])
    else:
        NullPath=out_path+'NullModel.'+dict_opts['--sample'].split('/')[-1]+'/'
    if not os.path.isdir(NullPath):
        os.system(r'''mkdir %s'''%(NullPath))
    return NullPath
def Null_Stats_Readin_One(NullPath,bamF,NullSplitLen_perc,genome_name,bam_files_appdix):
    fin=open(NullPath+bamF.split('/')[-1].replace('.'+bam_files_appdix,'.')+genome_name+'.Stats')
    pin=fin.readline().strip().split()
    pin=fin.readline().strip().split()
    pin=fin.readline().strip().split()
    Window_Size=int(float(pin[0])/3)
    sub_loc_size=max(int(float(pin[3]))*100,10**6)
    for line in fin:
        pin=line.strip().split()
    ReadLength=int(pin[-1].split(':')[1])
    if NullSplitLen_perc>1:
        NullSplitLen_perc=float(NullSplitLen_perc)/float(ReadLength)
    fin.close()
    return [Window_Size,ReadLength,sub_loc_size]
def many_RD_Process(copy_num_a,run_flag):
    Best_Letter_Rec=[[['a' for i in range(copy_num_a/2)],['a' for i in range(copy_num_a/2)]]]
    Best_Score_Rec=-100
    run_flag+=1     
    return([Best_Letter_Rec,Best_Score_Rec,run_flag])
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
def oppo_direct(x):
    if x=='+':
        return '-'
    elif x=='-':
        return '+'
def out_pair_modify_check(out_pair_modify,link_limit):
    for i in sorted(out_pair_modify.keys()):
        if len(out_pair_modify[i])>5:
            all_interval=[out_pair_modify[i][x+1]-out_pair_modify[i][x] for x in range(len(out_pair_modify[i])-1)]
            del out_pair_modify[i]
    return out_pair_modify
def out_pair_bp_check(out_pair_bp,link_limit):
    out=[]
    for x in out_pair_bp:
        if not x in out:
            if not x[0]==x[1]:
                out.append(x)
    out_hash={}
    for k1 in out:
        if not k1[0] in out_hash.keys():
            out_hash[k1[0]]=[]
        if not k1[1] in out_hash.keys():
            out_hash[k1[1]]=[]
        out_hash[k1[0]].append(k1[1])
        out_hash[k1[1]].append(k1[0])
    for x in out_hash.keys():
        if len(out_hash[x])>link_limit:
            del out_hash[x]
    out2=[]
    for x in out_hash.keys():
        for y in out_hash[x]:
            if not sorted([x,y]) in out2:
                out2.append(sorted([x,y]))
    return out2
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
def Prob_Possion(Number, Mean):
    #This function is used to calculate the probability of RD of each region
    #a possion model is applied
    #prob is produced in a log scale
    #float is transferred to the closest integer
    if not Mean==0:
        probs=[]
        from scipy.stats import poisson
        lamda=Mean
        K=Number
        K2=numpy.round(K)
        pdf=K2*math.log(lamda)-lamda-math.log(math.factorial(K2))
        return pdf
    elif Mean==0:
        return -Number
def Prob_NB(Number, Mean, Variance):
    if not Mean==0:
        P=1-Mean/Variance
        R=numpy.round(Mean*(1-P)/P)
        K2=numpy.round(Number)
        if P>0 and R>1 and K2+R-1>0 :
            log_pdf=math.log(math.factorial(K2+R-1))-math.log(math.factorial(K2))-math.log(math.factorial(R-1))+R*math.log(1-P)+K2*math.log(P)
        else:
            log_pdf=Prob_Possion(Number, Mean)
        return log_pdf
    elif Mean==0:
        return -Number
def Prob_Norm(Number, Mean, Variance):
    return -0.5*math.log(2*numpy.pi*Variance)-(Number-Mean)**2/2/Variance
def Prob_Standard_Norm(Number,Mean,Variance):
    Number_new=(float(Number)-float(Mean))/float(sqrt(Variance))
    Mean_new=0
    Variance_new=1
    return -0.5*math.log(2*numpy.pi*Variance_new)-(Number_new-Mean_new)**2/2/Variance_new
def pos_define_and_redefine(pbam,ReadLen,real_region,Window_Size,RD_RealRegion):
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
def qual_check_bps2(bps2):
    flag=0
    if len(bps2)==1 and len(bps2[0])==3:
        for x1 in bps2:
            for x2 in range(len(x1)-2):
                if int(x1[x2+2])-int(x1[x2+1])<100: 
                    flag=1
                elif int(x1[x2+2])-int(x1[x2+1])>10**7: 
                    flag=1
    if flag==0:
        return 'right'
    else:
        return 'error'
def RDStats_readin(RDStat):    
    RDStats={}
    fRDS=open(RDStat)
    pRDS1=fRDS.readline().strip().split()
    pRDS2=fRDS.readline().strip().split()
    for i in range(len(pRDS1)):
            RDStats[pRDS1[i]]=pRDS2[i]
    fRDS.close()
    return RDStats
def RD_Prob_NegativeBinomial(ReadDepth,RDStats):
    P=1-float(RDStats['Mean'])/(float(RDStats['STD'])**2)
    R=float(RDStats['Mean'])*(1-P)/P
    log_y=log_factorial(ReadDepth+R-1)-log_factorial(ReadDepth)-log_factorial(R-1)+R*numpy.log(1-P)+ReadDepth*numpy.log(P)
    return log_y
def RD_Prob_Normal(ReadDepth,RDStats):
    log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(RDStats['STD']))-(float(ReadDepth)-float(RDStats['Mean']))**2/2/float(RDStats['STD'])**2
    return log_y
def Segdup_readin(filein,chromo_name):
    data_hash={}
    fin=open(filein)
    pin=fin.readline().strip().split()
    for line in fin:
        pin=line.strip().split()
        if pin[0]==chromo_name:
            if not pin[0] in data_hash.keys():
                data_hash[pin[0]]={}
            if not int(pin[1]) in data_hash[pin[0]].keys():
                data_hash[pin[0]][int(pin[1])]=[]
            if not int(pin[2]) in data_hash[pin[0]][int(pin[1])]:
                data_hash[pin[0]][int(pin[1])].append(int(pin[2]))
    fin.close()
    data_list={}
    for k1 in data_hash.keys():
        data_list[k1]=[]
        for k2 in sorted(data_hash[k1].keys()):
            for k3 in sorted(data_hash[k1][k2]):
                data_list[k1].append([k2,k3])
    return data_list
def SP_file_readin(SP_filein):
    SP_hash={}
    fin=open(SP_filein)
    for line in fin:
        pin=line.strip().split()
        SP_hash[int(pin[1])]=int(pin[2])
    fin.close()
    SP_list=[]
    for k1 in sorted(SP_hash.keys()):
        SP_list.append([k1,SP_hash[k1]])
    return SP_list
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
def SP_Info_Cluster(SP_Info):
    for S_Chr in SP_Info.keys():
        bpCluster=clusterNums_bpintegrate(sorted(SP_Info[S_Chr].keys()), 11, 'f')
        bpMilti=[]
        for k1 in range(len(bpCluster[1])):
            if bpCluster[1][k1]>1:
                bpMilti.append(bpCluster[0][k1])
        if not bpMilti==[]:
            bpClu2=clusterSupVis_bpintegrate(sorted(SP_Info[S_Chr].keys()),sorted(bpMilti),11)
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
def SP_info_ReadIn_pre(bps_hash,S_Sample,bps_in_path):
    single_chromos=bps_hash[S_Sample].keys()
    SP_links={}
    SP_Link1={}
    SP_Info={}
    SP_Link2={}
    SP_Link3={}
    for k1 in bps_hash[S_Sample].keys():
        filein1=bps_hash[S_Sample][k1][0]
        fin=open(bps_in_path+filein1)
        for line in fin:
            pin=line.strip().split()
            if not pin[0] in SP_links.keys():
                SP_links[pin[0]]={}
            if not pin[0] in SP_Link1.keys():
                SP_Link1[pin[0]]=[]
            if not pin[0] in SP_Info.keys():
                SP_Info[pin[0]]={}
            if not int(pin[1]) in SP_Info[pin[0]].keys():
                SP_Info[pin[0]][int(pin[1])]=int(pin[2])
            else:
                SP_Info[pin[0]][int(pin[1])]=max([int(pin[2]),SP_Info[pin[0]][int(pin[1])]])
            if not int(pin[3]) in SP_Info[pin[0]].keys():
                SP_Info[pin[0]][int(pin[3])]=int(pin[4])
            else:
                SP_Info[pin[0]][int(pin[3])]=max(int(pin[4]),SP_Info[pin[0]][int(pin[3])])
            if not pin[1]==pin[2]:
                SP_Link1[pin[0]]+=[pin[1],pin[3]]
        fin.close()
        filein1=bps_hash[S_Sample][k1][1]
        fin=open(bps_in_path+filein1)
        for line in fin:
            pin=line.strip().split()
            if not pin[0] in SP_Info.keys():
                SP_Info[pin[0]]={}
            if not int(pin[1]) in SP_Info[pin[0]].keys():
                SP_Info[pin[0]][int(pin[1])]=int(pin[2])
            else:
                SP_Info[pin[0]][int(pin[1])]=max([int(pin[2]),SP_Info[pin[0]][int(pin[1])]])
        fin.close()
    SP_Info_Cluster(SP_Info)
    for S_Chr in SP_Link1.keys():
        SP_Link2[S_Chr]=[int(i) for i in SP_Link1[S_Chr]]
    for S_Chr in SP_Link1.keys():
        bpCluster=clusterSupVis_bpintegrate(SP_Info[S_Chr].keys(),SP_Link2[S_Chr],10)
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
def SP_Links_Info_Comple(SP_links,SP_link4):
    for x in SP_links.keys():
        temp=[]
        temp2=[]
        for k2 in sorted(SP_links[x].keys()):
            for k3 in sorted(SP_links[x][k2]):
                if not k3==k2:
                    temp.append([k2,k3])
        rec1=0
        rec2=0
        while True:
            if rec1==len(temp): break
            if rec2==len(SP_link4[x]): break
            if SP_link4[x][rec2][-1]<temp[rec1][0]:
                rec2+=1
            elif SP_link4[x][rec2][0]>temp[rec1][1]:
                rec1+=1
            else:
                if temp[rec1]==SP_link4[x][rec2]:
                    rec1+=1
                else:
                    flag_temp=0
                    for x1 in temp[rec1]:
                        if x1 in SP_link4[x][rec2]:
                            flag_temp+=1
                    if flag_temp<2:
                        if temp[rec1][1]-temp[rec1][0]>50:
                            temp2.append(temp[rec1])
                    else:
                        if len(SP_link4[x][rec2])>5:
                            if SP_link4[x][rec2].index(temp[rec1][1])-SP_link4[x][rec2].index(temp[rec1][0])>2:
                                if temp[rec1][1]-temp[rec1][0]>50:
                                    temp2.append(temp[rec1])
                    rec1+=1
        SP_link4[x]+=temp2
def SP_LN_info_ReadIn(LN_list,all_SPs,LN_filein,SP_filein):
    LN_Links={}
    fin=open(LN_filein)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in LN_Links.keys():
            LN_Links[pin[0]]={}
        if not int(pin[1]) in LN_Links[pin[0]].keys():
            LN_Links[pin[0]][int(pin[1])]=[]
        if not int(pin[3]) in LN_Links[pin[0]][int(pin[1])]:
            LN_Links[pin[0]][int(pin[1])].append(int(pin[3]))
        if not pin[0] in all_SPs.keys():
            all_SPs[pin[0]]={}
        if not int(pin[1]) in all_SPs[pin[0]].keys():
            all_SPs[pin[0]][int(pin[1])]=int(pin[2])
        else:
            all_SPs[pin[0]][int(pin[1])]=max([all_SPs[pin[0]][int(pin[1])],int(pin[2])])
        if not int(pin[3]) in all_SPs[pin[0]].keys():
            all_SPs[pin[0]][int(pin[3])]=int(pin[4])
        else:
            all_SPs[pin[0]][int(pin[3])]=max([all_SPs[pin[0]][int(pin[3])],int(pin[2])])
    fin.close()
    fin=open(SP_filein)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in all_SPs.keys():
            all_SPs[pin[0]]={}
        if not int(pin[1]) in all_SPs[pin[0]].keys():
            all_SPs[pin[0]][int(pin[1])]=int(pin[2])
        else:
            all_SPs[pin[0]][int(pin[1])]=max([all_SPs[pin[0]][int(pin[1])],int(pin[2])])
    fin.close()         
    for k1 in LN_Links.keys():
        if not k1 in LN_list.keys():
            LN_list[k1]=[]
        for k2 in sorted(LN_Links[k1].keys()):
            for k3 in sorted(LN_Links[k1][k2]):
                LN_list[k1].append([k2,k3])
    return [LN_list,all_SPs]
def SP_Info_Merge(SP_hash):
    #eg of SP_hash: SP_hash=SP_LN_Info[1]
    SP_pos=sorted(SP_hash.keys())
    SP_step=[SP_pos[i+1]-SP_pos[i] for i in range(len(SP_pos)-1)]
    min_step=min(SP_step)
    while True:
        if min(SP_step)>60: break
        min_index=[i for i in range(len(SP_step)) if SP_step[i]==min_step]
        for x in min_index:
            if SP_pos[x] in SP_hash.keys() and SP_pos[x+1] in SP_hash.keys():
                if SP_hash[SP_pos[x]]>SP_hash[SP_pos[x+1]]:
                    del SP_hash[SP_pos[x+1]]
                else:
                    del SP_hash[SP_pos[x]]
        SP_pos=sorted(SP_hash.keys())
        SP_step=[SP_pos[i+1]-SP_pos[i] for i in range(len(SP_pos)-1)]
        min_step=min(SP_step)
    return SP_pos
def struc_produce_two_block(Copy_num_estimate):
    cpNum_a=[j for j in [Copy_num_estimate['a']+i for i in [-1,0,1]] if j>-1]
    cpNum_b=[j for j in [Copy_num_estimate['b']+i for i in [-1,0,1]] if j>-1]
    all_Strucs=two_block_a(cpNum_a,cpNum_b)
    all_Str2=two_block_b(all_Strucs)
    return all_Str2
def simpler_struc_pick(list_of_structures,original_letters):
    #eg of list_of_structures:[[['b^', 'a^', 'b'], ['a', 'b']],[['b^', 'a', 'b'], ['a', 'b']]]
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
def split_loc_to_loc2(loc,ClusterLen):
    if not loc[1]-loc[0]>2*10**6:
        loc2=[loc]
    else:
        sublocNum=int(float(loc[1]-loc[0])/10**6)
        loc2=[[loc[0],loc[0]+10**6+int(ClusterLen)]]
        for slNum in range(sublocNum)[1:]:
            loc2.append([loc[0]+slNum*10**6-int(ClusterLen),loc[0]+(slNum+1)*10**6+int(ClusterLen)])
        if sublocNum>1:
            loc2.append([loc[0]+(slNum+1)*10**6-int(ClusterLen),loc[1]])
    return loc2
def SPLCff_Calculate(NullSplitLen_perc,SPLenStat,ReadLength):
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
            if SPLCff>ReadLength/20:
                SPLCff=ReadLength/20
        else:
            SPLCff=ReadLength/20
        if SPLCff<5:
            SPLCff=5
        fSPLStat.close()
        return SPLCff
def split_loc_to_subloc(loc,sub_loc_size,ClusterLen2):
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
def tempIL_Info_Add(clu_c_F,LinkCluMin,SPCluLen,abInfF,tempIL):
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
    return [tempIL,clu_c_F]
def tempIL_Filter_One(BPAlignQCFlank,BPAlignQC,ToolMappingQ,FileMappingQ,align_QCflag,tempIL,chrom):
    if not tempIL=={}:
        for aqb in tempIL.keys():
            if aqb<BPAlignQCFlank:
                del tempIL[aqb]
                continue
            if align_QCflag==1:
                if not 'chr' in chrom:
                    chrom_Mappability='chr'+chrom
                LNSPFR_aqb=os.popen(r'''%s %s %s %d %d 1'''%(ToolMappingQ,FileMappingQ,chrom_Mappability,aqb-BPAlignQCFlank,aqb+BPAlignQCFlank))
                tPairF_b=float(LNSPFR_aqb.read().strip()) 
                if tPairF_b<float(BPAlignQC):
                    del tempIL[aqb]
    return tempIL
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
def PathBP_SetUp(NullPath):
    path_BP='/'.join(NullPath.split('/')[:-2])+'/'+'.'.join(['BreakPoints']+NullPath.split('/')[-2].split('.')[1:])+'/'
    #path_BP=out_path+'BreakPoints.'+dict_opts['--sample'].split('/')[-1]+'/'
    if not os.path.isdir(path_BP):
        os.system(r'''mkdir %s'''%(path_BP))
    return path_BP
def pdf_bimodal_func(x,IL_Statistics):
    #eg of IL_Statistics=[mean1,mean2,std1,std2,alpha1,alpha2]
    return -IL_Statistics[4]*scipy.stats.norm.pdf(x,IL_Statistics[0],IL_Statistics[2])-IL_Statistics[5]*scipy.stats.norm.pdf(x,IL_Statistics[1],IL_Statistics[3])
def pdf_bimodal_calcu(x,IL_Statistics):
    #eg of IL_Statistics=[mean1,mean2,std1,std2,alpha1,alpha2]
    return IL_Statistics[4]*scipy.stats.norm.pdf(x,IL_Statistics[0],IL_Statistics[2])+IL_Statistics[5]*scipy.stats.norm.pdf(x,IL_Statistics[1],IL_Statistics[3])
def pdf_normal_calcu(x,Stat_list):
    #eg of Stat_list=[mean,std]
    return scipy.stats.norm.pdf(x,Stat_list[0],Stat_list[1])
def pdf_normal_func(x,Stat_list):
    #eg of Stat_list=[mean,std]
    return -scipy.stats.norm.pdf(x,Stat_list[0],Stat_list[1])
def find_max_bimodal(IL_Statistics):
    #return maximum pdf of IL_ditribution
    if not IL_Statistics[0]<IL_Statistics[1]: IL_new=IL_Statistics
    else:
        IL_new=[IL_Statistics[1],IL_Statistics[0],IL_Statistics[3],IL_Statistics[2],IL_Statistics[5],IL_Statistics[4]]
    return -scipy.optimize.fminbound(pdf_bimodal_func,IL_new[1],IL_new[0],args=(IL_new,),full_output=True)[1]
def find_max_negative_binomial(RD_Statistics):
    return -scipy.optimize.fminbound(pdf_normal_func,RD_Statistics[0]-1,RD_Statistics[0]+1,args=([RD_Statistics[0],RD_Statistics[1]],),full_output=True)[1]
def pdf_calculate_newer(x,IL_Statistics,upper_limit,lower_limit,Penalty_For_InsertLengthZero):
    [Mean1,Mean2,Std1,Std2,Alpha,Beta]=IL_Statistics
    if not Alpha==0:
        if x<upper_limit and x>lower_limit:
            return math.log(pdf_bimodal_calcu(x,IL_Statistics))
        else:
            return Penalty_For_InsertLengthZero
    else:
        if not x==0:
            return  math.log(1-Alpha)-math.log(math.sqrt(2*math.pi*math.pow(Std2,2)))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
        elif x==0:
            return Penalty_For_InsertLengthZero
def pdf_calculate(x,alpha,mean1,mean2,std1,std2,upper_limit,lower_limit,Penalty_For_InsertLengthZero):
    [Alpha,Mean1,Mean2,Std1,Std2]=[alpha,mean1,mean2,std1,std2]
    if not Alpha==0:
        if x<upper_limit and x>lower_limit:
            return math.log(Alpha/math.sqrt(2*math.pi*math.pow(Std1,2))*math.exp(-math.pow((x-Mean1),2)/(2*math.pow(Std1,2)))+(1-Alpha)/math.sqrt(2*math.pi*math.pow(Std2,2))*math.exp(-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))))
        else:
            return Penalty_For_InsertLengthZero
    else:
        if not x==0:
            return math.log(1-Alpha)-math.log(math.sqrt(2*math.pi*math.pow(Std2,2)))-math.pow((x-Mean2),2)/(2*math.pow(Std2,2))
        elif x==0:
            return Penalty_For_InsertLengthZero
def pos_block_assign(block_bps_chr,read_pos,tolerance_bp):
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
def RDNullfigure2_Modify(RDNullfigure2,Window_Size):
    fin=open(RDNullfigure2)
    pin1=fin.readline().strip().split()
    pin2=fin.readline().strip().split()
    fin.close()
    fo=open(RDNullfigure2,'w')
    print >>fo, ' '.join(pin1)
    print >>fo, ' '.join([str(i) for i in [float(j)/Window_Size for j in pin2]])
    fo.close()
def RD_Adj_Penal(GC_para_dict,Initial_GCRD_Adj,Chromo,RD_List,Let_List):
    Coverage_af_Adj=[]
    Letters=[['left']+Let_List[0]+['right'],['left']+Let_List[1]+['right']]
    Overall_Median_Coverage=float(GC_para_dict['GC_Overall_Median_Num'])
    Coverage_af_Adj=RD_List
    Theo_RD=GC_para_dict['GC_Overall_Median_Coverage'][str(Chromo)]
    Theo_Var=GC_para_dict['GC_Var_Coverage'][str(Chromo)]
    Prob_out=[]
    if Let_List==[[], []]:
        for i in Initial_GCRD_Adj.keys():
            if not i in ['left','right']:
                Prob_out.append(Prob_Norm(Initial_GCRD_Adj[i]*2,0,Theo_Var))
    else:
        for i in Coverage_af_Adj:
            for j in i:
                Prob_out.append(Prob_Norm(j*2,Theo_RD,Theo_Var))
    return numpy.mean(Prob_out)
def Read_Block_From_Position(bps,letters,bounce_length,position,end,flank):
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
def Reads_block_assignment(flank,bps,letters,block):
    if numpy.mean(block[:2])<bps[0]:
        return 'left'
    elif not numpy.mean(block[:2])<bps[-1]:
        return 'right'
    else:
        for i in letters:
            if not numpy.mean(block[:2])<bps[letters.index(i)] and numpy.mean(block[:2])<bps[letters.index(i)+1]:
                return i
def Reads_block_assignment_1(flank,bps,letters,position):
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
def Reads_block_assignment_2(bps,letters,block1,block2,flank):
    #In this function, we assign a read to the block where majority of the read falls in
    #'majority' is defined by the relative length. if the larger part of a read fells in block a, then the whole read are in block a 
    #to save time, we use this for now
    #may complex it later
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
def Reads_Direction_Detect_flag(flag):
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
def Region_Coverage_Calculate(sam_file,Number_Of_Windows,Region_Info,Window_Size):
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
def SamplingPercentage_readin(dict_opts):
    if '--null-copyneutral-perc' in dict_opts.keys():
        SamplingPercentage=float(dict_opts['--null-copyneutral-perc'])
    else:
        SamplingPercentage=0.05
    return SamplingPercentage
def SplitLenPNum_Calculate(SplitLength,NullPath,bamF,bam_file_appdix,genome_name,NullSplitLen_perc):
    SplitLengthPath=NullPath
    path_mkdir(SplitLengthPath)
    SplitLengthOutput=SplitLengthPath+'/'+bamF.split('/')[-1].replace('.'+bam_file_appdix,'')+'.'+genome_name+'.SplitLength'
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
    return [TotalILNum,TBNullDensity]
def TBStats_readin(TBStat):
    TBStats={}
    fTBS=open(TBStat)
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
    fTBS.close()
    return TBStats
def TB_Prob_Bimodal(TB_Length,TBStats):
    log_y1=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std1']))-(float(TB_Length)-float(TBStats['bimodal']['Mean1']))**2/2/float(TBStats['bimodal']['std1'])**2
    log_y2=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['bimodal']['std2']))-(float(TB_Length)-float(TBStats['bimodal']['Mean2']))**2/2/float(TBStats['bimodal']['std2'])**2
    return log_y1*float(TBStats['bimodal']['Bimodal1'])+log_y2*float(TBStats['bimodal']['Bimodal2'])
def TB_Prob_Normal(TB_Length,TBStats):
    log_y=-0.5*numpy.log(2*numpy.pi)-numpy.log(float(TBStats['normal']['std']))-(float(TB_Length)-float(TBStats['normal']['Mean']))**2/2/float(TBStats['normal']['std'])**2
    return log_y
def unify_list(list):
    out=[i for i in numpy.unique(list)]
    out.sort()
    return out
    #out=[]
    #for x in list:
    #    if not x in out:
    #        out.append(x)
    #out.sort()
    #return out
def write_bp_1a(LN_LN_Merge,bps_folder,chromo_name,S_Sample):
    fout_rec=[]
    fout=bps_folder+S_Sample+'.txt'
    if not os.path.isfile(fout):
        fo=open(fout,'w')
    else:
        fo=open(fout,'a')
    for x in LN_LN_Merge:
        print >>fo, ' '.join([chromo_name]+[str(i) for i in sorted(numpy.unique(x))])
        print >>fo, ' '
    fo.close()
def write_bp_2a(LN_LN_Merge,bps_folder,chromo_name,S_Sample):
    fout=bps_folder+S_Sample+'.'+str(chromo_name)+'.txt'
    if not os.path.isfile(fout):
        fo=open(fout,'w')
    else:
        fo=open(fout,'a')
    for x in LN_LN_Merge:
        print >>fo, ' '.join([chromo_name]+[str(i) for i in sorted(numpy.unique(x))])
        print >>fo, ' '
    fo.close()
def write_bp_3a(LN_LN_Merge,bps_folder,file_length,chromo_name,S_Sample):
    out_list=[[]]
    file_index=0
    for x in LN_LN_Merge:
        if len(out_list[-1])<file_length:
            out_list[-1].append(x)
        else:
            out_list.append([x])
    for x in out_list:
        file_index+=1
        fout=bps_folder+S_Sample+'.'+str(chromo_name)+'.'+str(file_index)+'.txt'
        fo=open(fout,'w')
        for y in x:
            print >>fo, ' '.join([chromo_name]+[str(i) for i in sorted(numpy.unique(y))])
            print >>fo, ' '
        fo.close()
def Write_filtered_LN_list(LN_keep_list,LN_filein,chromo_name):
    fo=open(LN_filein,'w')
    for k1 in LN_keep_list:
        print >>fo, ' '.join([str(i) for i in [chromo_name,k1[0],k1[2],k1[1],k1[3]]])
    fo.close()
def Write_filtered_SP_list(SP_keep_list,SP_filein,chromo_name):
    fo=open(SP_filein,'w')
    for x in SP_keep_list:
        print >>fo, ' '.join([str(i) for i in [chromo_name]+x])
    fo.close()
def x_qc_check(x,chromos_all):
    flag=0
    for y in x:
        if y in chromos_all:
            flag+=1
    return flag
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
def bp_to_hash(bp_list,sv_let,chromos):
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
def bp_to_let(del_info_unit,chromos):
    flag=0
    for i in del_info_unit[0]:
        if i in chromos:
            flag+=1
    if not flag==0:
        letter=''.join([chr(i+97) for i in range(len(del_info_unit[0])-2*flag)])
        letters='/'.join([letter,letter])
        return letters
    else:
        return 0
def bed_writing(filename,hash,chromos):
    fo=open(filename,'w')
    for k3 in chromos:
        if k3 in hash.keys():
            for k4 in hash[k3]:
                print >>fo, ' '.join([str(i) for i in k4])
    fo.close()
def score_rec_hash_Modify_for_short_del(Score_rec_hash):
    Score_rec_hash_new={}
    for x in sorted(Score_rec_hash.keys())[::-1][:1]:
        Score_rec_hash_new[x]=Score_rec_hash[x]
    for x in sorted(Score_rec_hash.keys())[::-1][1:]:
        Score_rec_hash_new[x-1.1]=Score_rec_hash[x]
    return Score_rec_hash_new
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
def delete_block_produce(Letter_List):
    delete_block=[x.upper() for x in Letter_List]
    return delete_block
def delete_BPs_produce(Letter_List):
    delete_BPs=[(j+1,j+2) for j in range(len(Letter_List))]
    return delete_BPs
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
def dup_type_decide(dup_let,flag_sex,k1,k2):
    if flag_sex==1:
        k1x=k1.split('/')[0]
        k2x=k2.split('/')[0]
    else:
        k1x=k1.split('/')[1]
        k2x=k2.split('/')[1]
    if not k2x=='':
        k1x_temp=[]
        for x_temp in k2x:
            if not x_temp=='^':
                k1x_temp.append(x_temp)
            else:
                k1x_temp[-1]+=x_temp
        k2x_temp=[[]]
        for x_temp in k1x_temp:
            if k2x_temp[-1]==[]:
                k2x_temp[-1].append(x_temp)
            else:
                if not '^' in x_temp and not '^' in k2x_temp[-1][-1] and ord(x_temp[0])-ord(k2x_temp[-1][-1][0])==1:
                        k2x_temp[-1].append(x_temp)
                elif '^' in x_temp and '^' in k2x_temp[-1][-1] and ord(x_temp[0])-ord(k2x_temp[-1][-1][0])==-1:
                        k2x_temp[-1].append(x_temp)
                else:
                    k2x_temp.append([x_temp])
        k3x_temp=[]
        for x_temp in k2x_temp:
            if not '^' in x_temp[0]:
                k3x_temp+=x_temp
            else:
                k3x_temp+=[x.replace('^','') for x in x_temp[::-1]]
        out=[]     
        for x in dup_let:
            out.append([])
            for y in x:
                if y[0]+y[0] in k2x:
                    out[-1].append('Tandem')
                pos_index=[z for z in range(len(k3x_temp)-len(y[0])+1) if ''.join(k3x_temp[z:z+len(y[0])]) == y[0]]
                inter_index=[pos_index[z+1]-pos_index[z] for z in range(len(pos_index)-1)]
                inter_index=[int(float(i)/float(len(y[0]))) for i in inter_index]
                while True:
                    if 1 in inter_index:
                        out[-1].append('Tandem')
                        inter_index.remove(1)
                    else:
                        break
                if not inter_index==[]:
                    out[-1].append('Disperse')
        return out
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
def GC_RD_Adj(GC_Median_Num,GC_Overall_Median_Num,Chromo,GC_Content,Coverage):
    Coverage_af_Adj=[]
    GC_Content=[float(i) for i in GC_Content]
    Coverage=[float(i) for i in Coverage]
    Overall_Median_Coverage=float(GC_Overall_Median_Num)
    for key_1 in range(len(GC_Content)):
            Be_Adj_RD=Coverage[key_1]
            GC_Con=GC_Content[key_1]
            GC_Con=int(round(GC_Con*100))
            if GC_Con in GC_Median_Num.keys():
                    Median_Coverage=GC_Median_Num[GC_Con]
                    Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
            elif not GC_Con in GC_Median_Num.keys():
                    Median_Coverage=Overall_Median_Coverage
                    Af_Adj_RD=Be_Adj_RD*Overall_Median_Coverage/Median_Coverage
            Coverage_af_Adj+=[Af_Adj_RD]
    return Coverage_af_Adj
def GC_RD_Correction(GC_Median_Num,GC_Overall_Median_Num,NullPath,BamN,ref_prefix,chrbam):
    ref_index=ref_prefix+'GC_Content'
    cov_index=NullPath+'RD_Stat/'+BamN+'.'+chrbam+'.RD.index'
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
def GC_index(ref_prefix,chrbam):
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
def GC_Stat_ReadIn(BamN,GC_Stat_Path,genome_name,affix):
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
    for k1 in CN2_Region.keys():
        for k2 in CN2_Region[k1].keys():
            if CN2_Region[k1][k2][0]=='@:':
                del CN2_Region[k1][k2]
            else:
                if sum([float(i) for i in CN2_Region[k1][k2][0].split(':')[1].split(',')])==0.0:
                    del CN2_Region[k1][k2]
        if CN2_Region[k1]=={}:
            del CN2_Region[k1]
    return [CN2_Region,Chromos,GC_Content]
def GC_Median_Num_Correct(GC_Median_Num):
    numbers=[]
    for k1 in sorted(GC_Median_Num.keys()):
        if not GC_Median_Num[k1]==0.0:
            numbers.append(GC_Median_Num[k1])
    for k1 in sorted(GC_Median_Num.keys()):
        if GC_Median_Num[k1]==0.0:
            GC_Median_Num[k1]=Median_Pick(numbers)
    return GC_Median_Num
def inv_flag(k1,k2):
    if '^' in k2:
        return 1
    else:
        return 0
def index_element(k1,ele):
    out=[]
    for x in range(len(k1)-len(ele)+1):
        if k1[x:(x+len(ele))]==ele:
            out.append(x)
    return out
def inv_flag_SA(k1,k2):
    out=0
    if '^' in k2:
        if k2.replace('^','')==k1:
            out+=1
    return out
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
    if len(k2_new[0])>1:
        for x1 in k2_new[0][1:]:
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
    if len(k2_new[1])>1:
        for x1 in k2_new[1][1:]:
            if not '^' in x1 and not '^' in out[1][-1]:
                if ord(x1)-ord(out[1][-1][-1])==1:
                    out[1][-1]+=x1
                else:
                    out[1].append(x1)
            elif '^' in x1 and '^' in out[1][-1]:
                if ord(x1[0])-ord(out[1][-1][0])==-1:
                    out[1][-1]=x1[0]+out[1][-1]
                else:
                    out[1].append(x1)
            else:
                out[1].append(x1)
    return out
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
def ref_base_readin(ref,chromo,pos):
    fref=os.popen(r'''samtools faidx %s %s:%s-%s'''%(ref,chromo,str(pos),str(pos)))
    tre=fref.readline().strip().split()
    REF_AL=fref.readline().strip().split()
    if not REF_AL==[]:
        return REF_AL[0]
    else:
        return 'N'
def RD_NB_stat_readin(RD_NB_Stat):
    if os.path.isfile(RD_NB_Stat):
        fin=open(RD_NB_Stat)
        for line in fin:
            pin=line.strip().split()
        fin.close()
        pdf_exp=norm.pdf(float(pin[1]),float(pin[1]),float(pin[2]))
        return numpy.log(pdf_exp)
    else:
        return 1
def RD_within_B_calcu(GC_Mean_Coverage,Full_Info,bps2):
    RD_within_B=Full_Info[0]
    RD_within_B['left']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
    RD_within_B['right']=numpy.mean([GC_Mean_Coverage[key_chr[0]] for key_chr in bps2])
    for i in RD_within_B.keys():
        RD_within_B[i+'^']=RD_within_B[i]
    return RD_within_B
def samtools_sort_process(mini_fout_N2,mini_fout_N3,mini_fout_N4):
    import subprocess
    try: subprocess.check_output("samtools", stderr=subprocess.STDOUT)
    except Exception, e: x= e.output
    samtools_type_decide=x.split('\n')
    version=samtools_type_decide[2].split(' ')[1].split('.')[0]
    if version=='0':
        os.system(r'''samtools sort %s %s'''%(mini_fout_N2,mini_fout_N3))
    else:
        os.system(r'''samtools sort %s -T %s -o %s'''%(mini_fout_N2,mini_fout_N3,mini_fout_N3+'.bam'))
    os.system(r'''samtools index %s '''%(mini_fout_N4))
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
                if outdup3[i+1][-1]==outdup3[i][-1] and ord(outdup3[i+1][0][0])-ord(outdup3[i][0][-1])==1 and temp2.count(outdup3[i][0]+outdup3[i+1][0])+temp2.count(outdup3[i][0]+outdup3[i+1][0]+'^')>1:
                    outdup4[-1][0]+=outdup3[i+1][0]
                else:
                    outdup4.append(outdup3[i+1])
        else:
            outdup4=outdup3
        return [outdel,outinv,outdup4,outtra]
def simple_flag_SA_svpredict(k1,k2):
    #print [k1,k2]
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
def total_rd_calcu(GC_Median_Num,GC_Overall_Median_Num,letter_RD2,letter_GC,chr_letter_bp,block_rd2):
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
def Uniparental_disomy_check(original_letters,bestletter):
    flag_out=0
    if sorted(original_letters+original_letters)==sorted(bestletter[0]+bestletter[1]):
        for x in bestletter:
            if len(x)>1:
                for y in range(len(x[1:])):
                    if ord(x[y+1])-ord(x[y])<0:
                        flag_out+=1
    else:
        flag_out+=1
    if flag_out==0:
        print 'Uniparental_disomy: '+'/'.join([''.join(bestletter[0]),''.join(bestletter[1])])
        return 'Uniparental_disomy'
    else:
        return 'Pass'
def chromo_contig_prepare(workdir):
    out=[]
    ppre=path_modify(workdir)
    ref_file=workdir+'reference_SVelter/genome.fa'
    out.append('##reference='+ref_file)
    if os.path.isfile(ref_file+'.fai'):
        fin=open(ref_file+'.fai')
        for line in fin:
            pin=line.strip().split()
            if not '_' in pin[0] and not '-' in pin[0]:
                out.append('##contig=<ID='+pin[0]+',length='+pin[1]+'>')
        fin.close()
    return out
def write_VCF_header(output_file,time,workdir):
    fo=open(output_file,'w')
    print>>fo, '##fileformat=VCFv4.1'
    print>>fo,'##fileDate='+time.strftime("%Y%m%d")
    ref_info=chromo_contig_prepare(workdir)
    for x in ref_info:
        print >>fo,x
    #print>>fo,'##reference=hg19'
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
    print>>fo,'##INFO=<ID=SVCLASS,Number=1,Type=String,Description="Classification of structural variant">'
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
def sv_reorganize_reorder(sv_reorganize):
    out={}
    for k1 in sv_reorganize.keys():
        out[k1]={}
        for k2 in sv_reorganize[k1].keys():
            for k3 in sv_reorganize[k1][k2].keys():
                for k4 in sv_reorganize[k1][k2][k3]:
                    if not k4[1] in out[k1].keys():
                        out[k1][k4[1]]={}
                    if not k4[2] in out[k1][k4[1]].keys():
                        out[k1][k4[1]][k4[2]]=[]
                    out[k1][k4[1]][k4[2]].append(k4[:10])
    return out
def write_VCF_main(output_file,sv_out,chromos,ref,sv_type_record):
    fo=open(output_file,'a')
    sv_reorganize={}
    for k1 in sv_out.keys():
        #sv_reorganize[k1]={}
        for k2 in sv_out[k1].keys():
            start=int(k2.split(';')[1])
            #if not start in sv_reorganize[k1].keys():
                #sv_reorganize[k1][start]={}
            SVtemp_a=[]
            SVtemp_b=[]
            for k3 in sv_out[k1][k2]:
                if not k3[:-1] in SVtemp_a:
                    SVtemp_a.append(k3[:-1])
                    SVtemp_b.append([k3[-1]])
                else:
                    SVtemp_b[SVtemp_a.index(k3[:-1])].append(k3[-1])
            SVtemp=[]
            #sv_reorganize[k1][start][k2]=[]
            for k3 in range(len(SVtemp_a)):
                if len(SVtemp_b[k3])==2 and SVtemp_b[k3] in [['0|1', '1|0'],['1|0', '0|1']]:
                    SVtemp_b[k3]=['1|1']
            for k3 in range(len(SVtemp_a)):
                for k4 in SVtemp_b[k3]:
                    SVtemp_info=SVtemp_a[k3]+[k4]
                    if not SVtemp_info[0] in sv_reorganize.keys():
                        sv_reorganize[SVtemp_info[0]]={}
                    if not SVtemp_info[1] in sv_reorganize[SVtemp_info[0]].keys():
                        sv_reorganize[SVtemp_info[0]][SVtemp_info[1]]={}
                    if not SVtemp_info[2] in sv_reorganize[SVtemp_info[0]][SVtemp_info[1]].keys():
                        sv_reorganize[SVtemp_info[0]][SVtemp_info[1]][SVtemp_info[2]]=[]
                    if not SVtemp_info in sv_reorganize[SVtemp_info[0]][SVtemp_info[1]][SVtemp_info[2]]:
                        sv_reorganize[SVtemp_info[0]][SVtemp_info[1]][SVtemp_info[2]].append(SVtemp_info)
                    #sv_reorganize[k1][start][k2].append(SVtemp_a[k3]+[k4])
    for k1 in sv_reorganize.keys():
        for k2 in sv_reorganize[k1].keys():
            for k3 in sv_reorganize[k1][k2].keys():
                for k4 in sv_reorganize[k1][k2][k3]:
                    if k4[5]>0: continue
                    else:
                        k4[5]=0.0
                        k4[6]='LowQual'
    sv_reorganize=sv_reorganize_reorder(sv_reorganize)
    for k1 in chromos:
        if k1 in sv_reorganize.keys():
            for k2 in sorted(sv_reorganize[k1].keys()):
                for k3 in sorted(sv_reorganize[k1][k2].keys()):
                    for k4 in sv_reorganize[k1][k2][k3]:
                        if k4[3]=='N':
                            k4[3]=ref_base_readin(ref,k4[0],k4[1])
                        if ':'.join(k4[2].split(';')[:-1]) in sv_type_record.keys():
                            SV_type=sv_type_record[':'.join(k4[2].split(';')[:-1])]
                        else:
                            SV_type=['unclassified/unclassified']
                        k4[7]+=';SVCLASS='+SV_type[0]
                        print >>fo, '\t'.join([str(i) for i in k4])
    fo.close()
def zero_RD_Process(original_bp_list,run_flag,Best_IL_Score,Best_RD_Score):
    Best_Letter_Rec=[[[], []]]
    Af_Letter=[[],[]]
    Af_BP=[[original_bp_list[0]],[original_bp_list[0]]]
    Best_Score_Rec=Best_IL_Score+Best_RD_Score
    run_flag+=1
    return([Best_Letter_Rec,Best_Score_Rec,run_flag])
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
def Insert_len_stat_readin(Insert_Len_Stat):
    if os.path.isfile(Insert_Len_Stat):
        fin=open(Insert_Len_Stat)
        for line in fin:
            pin=line.strip().split()
        fin.close()
        pdf_exp=norm.pdf(float(pin[1]),float(pin[1]),float(pin[2]))
        return numpy.log(pdf_exp)
    else:
        return 1
def sv_info_qc_score_extract(sv_info):
    sv_info_out={}
    qc_score_out={}
    for k1 in sv_info.keys():
        sv_info_out[k1]={}
        for k2 in sv_info[k1].keys():
            sv_info_out[k1][k2]=[]
            for k3 in sv_info[k1][k2]:
                sv_info_out[k1][k2].append(k3[:-1])
                qc_score_out[':'.join(k3[:-1])]=k3[-1]
    return [sv_info_out,qc_score_out]
def sv_info_score_modify(sv_info):
    all_scores_pos=0
    all_scores_neg=0
    for x in sv_info.keys():
        for y in sv_info[x].keys():
            if y.split('/')[0].count('a')>2 or y.split('/')[1].count('a')>2:
                for z in sv_info[x][y]:
                    z[-1]=0
            elif y=='/':
                for z in sv_info[x][y]:
                    z[-1]=100
            else:
                for z in sv_info[x][y]:
                    if z[-1]>0:
                        if z[-1]>all_scores_pos:
                            all_scores_pos=z[-1]
                    elif z[-1]<0:
                        if z[-1]<all_scores_neg:
                            all_scores_neg=z[-1]
    for x in sv_info.keys():
        for y in sv_info[x].keys():
            if not y=='/':
                for z in sv_info[x][y]:
                    if z[-1]>0:
                        z[-1]=int(float(z[-1])/float(all_scores_pos)*100)
                    elif z[-1]<0:
                        z[-1]=int(float(z[-1])/float(-all_scores_neg)*100)
    return sv_info