def read_in_run_time(log_file):
        fin=os.popen(r'''head -1 %s'''%(log_file))
        pin=fin.readline().strip().split()
        start=float(pin[0])
        fin.close()
        fin=os.popen(r'''tail -1 %s'''%(log_file))
        pin=fin.readline().strip().split()
        if len(pin)>1:
                end='error'
        else:
                end=float(pin[0])
        fin.close()
        if not end=='error':
                return [end-start]
        else:
                return ['*']

def read_in_run_time_SVeltr(log_file):
    fin=os.popen(r'''tail -2 %s'''%(log_file))
    pin=fin.readline().strip().split()
    temp=[float(i) for i in pin[3:]]+[float(pin[2].split(':')[1])]
    fin.close()
    return temp

def chromos_readin(ref):
        fin=open(ref+'.fai')
        chromos=[]
        for line in fin:
                pin=line.strip().split()
                chromos.append(pin[0])
        fin.close()
        return chromos

import os
run_time_hash={}
pbs_path='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/pbs/'

delly_path=pbs_path+'pbs_Delly/'
run_time_hash['delly']={}
for k1 in os.listdir(delly_path):
        if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-2]
                run_time_hash['delly'][chromo_name]=read_in_run_time(delly_path+k1)


lumpy_path=pbs_path+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
run_time_hash['lumpy']={}
for k1 in os.listdir(lumpy_path_1):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['lumpy'][chromo_name]=read_in_run_time(lumpy_path_1+k1)

for k1 in os.listdir(lumpy_path_2):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-2]
                run_time_hash['lumpy'][chromo_name]+=read_in_run_time(lumpy_path_2+k1)


svelter_path=pbs_path+'pbs_SVelter/'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]=read_in_run_time_SVeltr(svelter_path+k1)


svelter_path=pbs_path+'pbs_SVelter'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'1/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]=read_in_run_time(svelter_path+'1/'+k1)

for k1 in os.listdir(svelter_path+'2/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-2]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'2/'+k1)

for k1 in os.listdir(svelter_path+'3/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'3/'+k1)

for k1 in os.listdir(svelter_path+'4/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'4/'+k1)

for k1 in os.listdir(svelter_path+'5/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-2]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'5/'+k1)



pindel_path=pbs_path+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
run_time_hash['pindel']={}
for k1 in os.listdir(pindel_path_1):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-4]
                run_time_hash['pindel'][chromo_name]=read_in_run_time(pindel_path_1+k1)

for k1 in os.listdir(pindel_path_2):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-4]
                if chromo_name in run_time_hash['pindel'].keys():
                        run_time_hash['pindel'][chromo_name]+=read_in_run_time(pindel_path_2+k1)


erds_path=pbs_path+'pbs_erds/'
run_time_hash['erds']={}
for k1 in os.listdir(erds_path):
        if k1.split('.')[-1]=='log' and 'chr' in k1:
                chromo_name=k1.split('.')[-3]
                run_time_hash['erds'][chromo_name]=read_in_run_time(erds_path+k1)


pbs_path='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19_deterministic_2/pbs/'
svelter_path=pbs_path+'pbs_SVelter/'
run_time_hash['svelter_Deter']={}
for k1 in os.listdir(svelter_path):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter_Deter'][chromo_name]=read_in_run_time_SVeltr(svelter_path+k1)


algorithms=run_time_hash.keys()
chromos=chromos_readin('/scratch/remills_flux/xuefzhao/reference/hg19_platinum/hg19_platinum.fa')
fo=open(pbs_path+'Run_Time_assessment.NA12878.txt','w')
print >>fo, ' '.join(['algorithms']+algorithms)
for k1 in chromos:
        test=[]
        for k2 in algorithms:
                if k1 in run_time_hash[k2].keys():
                        test.append(sum(run_time_hash[k2][k1])/3600)
                else:
                        test.append('*')
        print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()


#Calculate PacbioValidation
ppre='/scratch/remills_flux/xuefzhao/NA12878.NGS/hg19/VaLoR_Vali/'
import os
data_hash={}
for k1 in os.listdir(ppre):
    if k1.split('.')[-1]=='PacVal':
        fin=open(k1)
        for line in fin:
            pin=line.strip().split()
            if not pin[-4] in data_hash.keys():
                data_hash[pin[-4]]={}
            if not pin[-3] in data_hash[pin[-4]].keys():
                data_hash[pin[-4]][pin[-3]]={}
            if not pin[0] in data_hash[pin[-4]][pin[-3]].keys():
                data_hash[pin[-4]][pin[-3]][pin[0]]=[]
            data_hash[pin[-4]][pin[-3]][pin[0]].append(float(pin[-1]))
        fin.close()

stat_hash={}
fo=open(ppre+'Integrated.Validated.Result.Cff0.1','w')
for k1 in data_hash.keys():
    stat_hash[k1]={}
    for k2 in data_hash[k1].keys():
        stat_hash[k1][k2]=[]
        test=[]
        for k3 in data_hash[k1][k2].keys():
            test.append([len([i for i in data_hash[k1][k2][k3] if i>0.1]),len(data_hash[k1][k2][k3])])
        stat_hash[k1][k2]+=[sum([i[0] for i in test]),sum([i[1] for i in test])]
        print>>fo, ' '.join([str(i) for i in [k1,k2]+stat_hash[k1][k2]]) 

fo.close()













