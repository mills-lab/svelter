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
ppre='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/Simulate_each_chr/pbs/'
run_time_hash={}
delly_path=ppre+'pbs_Delly/'
for k1 in os.listdir(delly_path):
        if k1.split('.')[-1]=='log':
                if not k1.split('.')[0] in run_time_hash.keys():
                        run_time_hash[k1.split('.')[0]]={}
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                if not key_2 in run_time_hash[k1.split('.')[0]].keys():
                        run_time_hash[k1.split('.')[0]][key_2]={}
                key_3=k1.split('.')[-7]
                if not key_3 in run_time_hash[k1.split('.')[0]][key_2].keys():
                        run_time_hash[k1.split('.')[0]][key_2][key_3]={}
                key_4=k1.split('.')[-3]
                if not key_4 in run_time_hash[k1.split('.')[0]][key_2][key_3].keys():
                        run_time_hash[k1.split('.')[0]][key_2][key_3][key_4]=[]
                run_time_hash[k1.split('.')[0]][key_2][key_3][key_4]=read_in_run_time(delly_path+k1)

lumpy_path=ppre+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
for k1 in os.listdir(lumpy_path_1):
        if k1.split('.')[-1]=='log':
                key_1='Lumpy'
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[4:6])
                else:
                        key_2=k1.split('.')[4]
                key_3=k1.split('.')[-8]
                key_4=k1.split('.')[-4]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]=read_in_run_time(lumpy_path_1+k1)

for k1 in os.listdir(lumpy_path_2):
        if k1.split('.')[-1]=='log':
                key_1='Lumpy'
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-7]
                key_4=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]+=read_in_run_time(lumpy_path_2+k1)

svelter_path=ppre+'pbs_SVelter/'
for k1 in os.listdir(svelter_path):
        if k1.split('.')[-1]=='log':
                key_1=k1.split('.')[0]
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-7]
                key_4=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]=read_in_run_time(svelter_path+k1)

pindel_path=ppre+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
for k1 in os.listdir(pindel_path_1):
        if k1.split('.')[-1]=='log':
                key_1=k1.split('.')[0]
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-7]
                key_4=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]=read_in_run_time(pindel_path_1+k1)

for k1 in os.listdir(pindel_path_2):
        if k1.split('.')[-1]=='log':
                key_1='Pindel'
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-7]
                key_4=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]+=read_in_run_time(pindel_path_2+k1)

erds_path=ppre+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step2_SNPs_calling/'
erds_path_2=erds_path+'pbs_erders_step3_erds/'
for k1 in os.listdir(erds_path_1):
        if k1.split('.')[-1]=='log':
                key_1='erds.SNPs'
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-7]
                key_4=k1.split('.')[-3]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]=read_in_run_time(erds_path_1+k1)

for k1 in os.listdir(erds_path_2):
        if k1.split('.')[-1]=='log':
                key_1=k1.split('.')[0]
                if 'comp' in k1:
                        key_2='_'.join(k1.split('.')[1:3])
                else:
                        key_2=k1.split('.')[1]
                key_3=k1.split('.')[-8]
                print key_3
                key_4=k1.split('.')[-4]
                if not key_1 in run_time_hash.keys():
                        run_time_hash[key_1]={}
                if not key_2 in run_time_hash[key_1].keys():
                        run_time_hash[key_1][key_2]={}
                if not key_3 in run_time_hash[key_1][key_2].keys():
                        run_time_hash[key_1][key_2][key_3]={}
                if not key_4 in run_time_hash[key_1][key_2][key_3].keys():
                        run_time_hash[key_1][key_2][key_3][key_4]=[]
                run_time_hash[key_1][key_2][key_3][key_4]+=read_in_run_time(erds_path_2+k1)

algorithms=sorted(run_time_hash.keys())
chromosomes=['1','2']
simus=['het', 'homo','comp_het', 'comp_homo']
RDs=['RD10', 'RD20', 'RD30', 'RD40', 'RD50']

fo=open(ppre+'Run_Time_assessment.Simulate.chr1.chr2','w')
print >>fo, ' '.join(['chromosome','event','RD']+algorithms)
for k2 in simus:
        for k3 in RDs:
                for k1 in chromosomes:
                        test=[]
                        for k4 in algorithms:
                                if not '*'  in run_time_hash[k4][k2][k3][k1]:
                                        test.append(sum(run_time_hash[k4][k2][k3][k1])/3600)
                                else:
                                        test.append('*')
                        print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+test])

fo.close()



