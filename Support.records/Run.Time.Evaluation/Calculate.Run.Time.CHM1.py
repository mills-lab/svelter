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

def chromos_readin(ref):
    fin=open(ref+'.fai')
    chromos=[]
    for line in fin:
            pin=line.strip().split()
            chromos.append(pin[0])
    fin.close()
    return chromos

def read_in_run_time_SVeltr(log_file):
    fin=os.popen(r'''tail -2 %s'''%(log_file))
    pin=fin.readline().strip().split()
    temp=[float(i) for i in pin[3:]]+[float(pin[2].split(':')[1])]
    fin.close()
    return temp

import os
ppre='/scratch/remills_flux/xuefzhao/CHM1/IL500/hg19/pbs/'
run_time_hash={}
delly_path=ppre+'pbs_Delly/'
for k1 in os.listdir(delly_path):
    if k1.split('.')[-1]=='log':
            key_1='Delly'
            key_2=k1.split('.')[2]
            if not key_1 in run_time_hash.keys():
                    run_time_hash[key_1]={}
            run_time_hash[key_1][key_2]=read_in_run_time(delly_path+k1)

lumpy_path=ppre+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
for k1 in os.listdir(lumpy_path_1):
        if k1.split('.')[-1]=='log':
                key_1='Lumpy'
                key_2=k1.split('.')[-4]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time(lumpy_path_1+k1)

for k1 in os.listdir(lumpy_path_2):
        if k1.split('.')[-1]=='log':
                key_1='Lumpy'
                key_2=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]+=read_in_run_time(lumpy_path_2+k1)

svelter_path=ppre+'pbs_SVelter/'
for k1 in os.listdir(svelter_path):
        if k1.split('.')[-1]=='log':
                key_1='SVelter'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time_SVeltr(svelter_path+k1)




svelter_path=ppre+'pbs_SVelter'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'1/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-4]
                run_time_hash['svelter'][chromo_name]=read_in_run_time(svelter_path+'1/'+k1)

for k1 in os.listdir(svelter_path+'2/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-2]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'2/'+k1)

for k1 in os.listdir(svelter_path+'3/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-4]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'3/'+k1)

for k1 in os.listdir(svelter_path+'4/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'4/'+k1)

for k1 in os.listdir(svelter_path+'5/'):
         if k1.split('.')[-1]=='log':
                chromo_name=k1.split('.')[-3]
                run_time_hash['svelter'][chromo_name]+=read_in_run_time(svelter_path+'5/'+k1)



svelter_path_deter='/scratch/remills_flux/xuefzhao/CHM1/IL500/hg19_deterministic/pbs/pbs_SVelter/'
for k1 in os.listdir(svelter_path_deter):
        if k1.split('.')[-1]=='log':
                key_1='SVelter_Deter'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time_SVeltr(svelter_path_deter+k1)

pindel_path=ppre+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
for k1 in os.listdir(pindel_path_1):
        if k1.split('.')[-1]=='log':
                key_1='Pindel'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time(pindel_path_1+k1)

for k1 in os.listdir(pindel_path_2):
        if k1.split('.')[-1]=='log':
                key_1='Pindel'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]+=read_in_run_time(pindel_path_2+k1)

erds_path=ppre+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step2_SNPs_calling/'
erds_path_2=erds_path+'pbs_erders_step2_SNPs_calling/'
for k1 in os.listdir(erds_path_1):
        if k1.split('.')[-1]=='log':
                key_1='erds.SNPs'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time(erds_path_1+k1)

for k1 in os.listdir(erds_path_2):
        if k1.split('.')[-1]=='log':
                key_1='erds'
                key_2=k1.split('.')[3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                run_time_hash[key_1][key_2]=read_in_run_time(erds_path_2+k1)

algorithms=sorted(run_time_hash.keys())
chromosomes=run_time_hash['Delly'].keys()
fo=open(ppre+'Run_Time_assessment.CHM1.txt','w')
print >>fo, ' '.join(['algorithms']+algorithms)
for k1 in chromosomes:
        test=[]
        for k2 in algorithms:
                if k1 in run_time_hash[k2].keys():
                        if not '*' in run_time_hash[k2][k1]:
                                test.append(sum(run_time_hash[k2][k1])/3600)
                        else:
                              test.append('*')   
                else:
                        test.append('*')
        print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()
