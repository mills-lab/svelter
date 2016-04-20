def read_in_run_time(log_file):
    fin=os.popen(r'''head -1 %s'''%(log_file))
    pin=fin.readline().strip().split()
    start=float(pin[0])
    fin.close()
    fin=os.popen(r'''tail -1 %s'''%(log_file))
    pin=fin.readline().strip().split()
    end=float(pin[0])
    fin.close()
    return [end-start]

def chromos_readin(ref):
    fin=open(ref+'.fai')
    chromos=[]
    for line in fin:
        pin=line.strip().split()
        chromos.append(pin[0])
    fin.close()
    return chromos

import os
ppre='/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.FussyJunc/'

#Het
run_time_hash={}
#Delly
het=ppre+'Simulate.het/pbs/'
delly_path=het+'pbs_Delly/'
run_time_hash['delly']={}
for k1 in os.listdir(delly_path):
    if k1.split('.')[-1]=='log':
        RD=k1.split('.')[1]
        if not RD in run_time_hash['delly'].keys():
            run_time_hash['delly'][RD]=[read_in_run_time(delly_path+k1)]
        else:
            run_time_hash['delly'][RD]+=[read_in_run_time(delly_path+k1)]

#Lumpy
lumpy_path=het+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
run_time_hash['lumpy']={}
for k1 in os.listdir(lumpy_path_1):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-5]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_1+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_1+k1)]

for k1 in os.listdir(lumpy_path_2):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-4]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_2+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_2+k1)]

#SVelter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#SVelter_Deter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter_deterministic']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter_deterministic'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4_deterministic/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4_deterministic/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#pindel
pindel_path=het+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
run_time_hash['pindel']={}
for k1 in os.listdir(pindel_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['pindel'][chromo_name]=[read_in_run_time(pindel_path_1+k1)]

for k1 in os.listdir(pindel_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-2]
        if chromo_name in run_time_hash['pindel'].keys():
            run_time_hash['pindel'][chromo_name]+=[read_in_run_time(pindel_path_2+k1)]

#erds
erds_path=het+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step3_erds/'
run_time_hash['erds']={}
for k1 in os.listdir(erds_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['erds'][chromo_name]=[read_in_run_time(erds_path_1+k1)]

erds_path_2=erds_path+'pbs_erders_step2_SNPs_calling/'
run_time_hash['erds_SNP']={}
for k1 in os.listdir(erds_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['erds_SNP'][chromo_name]=[read_in_run_time(erds_path_2+k1)]

algorithms=run_time_hash.keys()
rds=sorted(run_time_hash['delly'].keys())

fo=open(ppre+'Run_time_Sumi','w')
print >>fo,' '.join(['het']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([sum(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo,' '.join([k1]+[str(i) for i in test])

print >>fo,' '.join(['het_Paralelle']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([max(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()




#Homo
run_time_hash={}

#Delly
het=ppre+'Simulate.homo/pbs/'
delly_path=het+'pbs_Delly/'
run_time_hash['delly']={}
for k1 in os.listdir(delly_path):
    if k1.split('.')[-1]=='log':
        RD=k1.split('.')[1]
        if not RD in run_time_hash['delly'].keys():
            run_time_hash['delly'][RD]=[read_in_run_time(delly_path+k1)]
        else:
            run_time_hash['delly'][RD]+=[read_in_run_time(delly_path+k1)]

#Lumpy
lumpy_path=het+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
run_time_hash['lumpy']={}
for k1 in os.listdir(lumpy_path_1):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-5]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_1+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_1+k1)]

for k1 in os.listdir(lumpy_path_2):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-4]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_2+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_2+k1)]

#SVelter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['svelter'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
        run_time_hash['svelter'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4/'):
    if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#SVelter_Deter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter_deterministic']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter_deterministic'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4_deterministic/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4_deterministic/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#pindel
pindel_path=het+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
run_time_hash['pindel']={}
for k1 in os.listdir(pindel_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['pindel'][chromo_name]=[read_in_run_time(pindel_path_1+k1)]

for k1 in os.listdir(pindel_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-2]
        if chromo_name in run_time_hash['pindel'].keys():
            run_time_hash['pindel'][chromo_name]+=[read_in_run_time(pindel_path_2+k1)]

#erds
erds_path=het+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step3_erds/'
run_time_hash['erds']={}
for k1 in os.listdir(erds_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['erds'][chromo_name]=[read_in_run_time(erds_path_1+k1)]

erds_path_2=erds_path+'pbs_erders_step2_SNPs_calling/'
run_time_hash['erds_SNP']={}
for k1 in os.listdir(erds_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['erds_SNP'][chromo_name]=[read_in_run_time(erds_path_2+k1)]

fo=open(ppre+'Run_time_Sumi','a')
print >>fo,' '.join(['homo']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([sum(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo,' '.join([k1]+[str(i) for i in test])


print >>fo,' '.join(['homo_Paralelle']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([max(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()






#Comp Het
run_time_hash={}
#Delly
het=ppre+'Simulate.comp.het/pbs/'
delly_path=het+'pbs_Delly/'
run_time_hash['delly']={}
for k1 in os.listdir(delly_path):
    if k1.split('.')[-1]=='log':
        RD=k1.split('.')[2]
        if not RD in run_time_hash['delly'].keys():
            run_time_hash['delly'][RD]=[read_in_run_time(delly_path+k1)]
        else:
            run_time_hash['delly'][RD]+=[read_in_run_time(delly_path+k1)]

#Lumpy
lumpy_path=het+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
run_time_hash['lumpy']={}
for k1 in os.listdir(lumpy_path_1):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-5]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_1+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_1+k1)]

for k1 in os.listdir(lumpy_path_2):
     if k1.split('.')[-1]=='log':
        RD=k1.split('.')[-4]
        if not RD in run_time_hash['lumpy'].keys():
            run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_2+k1)]
        else:
            run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_2+k1)]

#SVelter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#SVelter_Deter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter_deterministic']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter_deterministic'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4_deterministic/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4_deterministic/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#pindel
pindel_path=het+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
run_time_hash['pindel']={}
for k1 in os.listdir(pindel_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['pindel'][chromo_name]=[read_in_run_time(pindel_path_1+k1)]

for k1 in os.listdir(pindel_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-2]
        if chromo_name in run_time_hash['pindel'].keys():
            run_time_hash['pindel'][chromo_name]+=[read_in_run_time(pindel_path_2+k1)]

#erds
erds_path=het+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step3_erds/'
run_time_hash['erds']={}
for k1 in os.listdir(erds_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['erds'][chromo_name]=[read_in_run_time(erds_path_1+k1)]

erds_path_2=erds_path+'pbs_erders_step2_SNPs_calling/'
run_time_hash['erds_SNP']={}
for k1 in os.listdir(erds_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['erds_SNP'][chromo_name]=[read_in_run_time(erds_path_2+k1)]

algorithms=run_time_hash.keys()
rds=sorted(run_time_hash['delly'].keys())

fo=open(ppre+'Run_time_Sumi','a')
print >>fo,' '.join(['comp.het']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([sum(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo,' '.join([k1]+[str(i) for i in test])

print >>fo,' '.join(['comp.het_Paralelle']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([max(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()








#Comp Homo
run_time_hash={}
#Delly
het=ppre+'Simulate.comp.homo/pbs/'
delly_path=het+'pbs_Delly/'
run_time_hash['delly']={}
for k1 in os.listdir(delly_path):
    if k1.split('.')[-1]=='log':
        RD=k1.split('.')[3]
        if not RD in run_time_hash['delly'].keys():
            run_time_hash['delly'][RD]=[read_in_run_time(delly_path+k1)]
        else:
            run_time_hash['delly'][RD]+=[read_in_run_time(delly_path+k1)]

#Lumpy
lumpy_path=het+'pbs_Lumpy/'
lumpy_path_1=lumpy_path+'Prepare.Ab.Reads.axiom/'
lumpy_path_2=lumpy_path+'Lumpy2/'
run_time_hash['lumpy']={}
for k1 in os.listdir(lumpy_path_1):
    if k1.split('.')[-1]=='log':
       RD=k1.split('.')[-5]
       if not RD in run_time_hash['lumpy'].keys():
          run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_1+k1)]
       else:
          run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_1+k1)]

for k1 in os.listdir(lumpy_path_2):
    if k1.split('.')[-1]=='log':
       RD=k1.split('.')[-4]
       if not RD in run_time_hash['lumpy'].keys():
          run_time_hash['lumpy'][RD]=[read_in_run_time(lumpy_path_2+k1)]
       else:
          run_time_hash['lumpy'][RD]+=[read_in_run_time(lumpy_path_2+k1)]

#SVelter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#SVelter_Deter
svelter_path=het+'pbs_SVelter/'
run_time_hash['svelter_deterministic']={}
for k1 in os.listdir(svelter_path+'SVelter1/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name]=[read_in_run_time(svelter_path+'SVelter1/'+k1),[]]

for k1 in os.listdir(svelter_path+'SVelter2/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter2/'+k1)

for k1 in os.listdir(svelter_path+'SVelter3/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-4]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter3/'+k1))
           run_time_hash['svelter_deterministic'][chromo_name].append([])

for k1 in os.listdir(svelter_path+'SVelter4_deterministic/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name][-1]+=read_in_run_time(svelter_path+'SVelter4_deterministic/'+k1)

for k1 in os.listdir(svelter_path+'SVelter5/'):
    if k1.split('.')[-1]=='log':
           chromo_name=k1.split('.')[-5]
           run_time_hash['svelter_deterministic'][chromo_name].append(read_in_run_time(svelter_path+'SVelter5/'+k1))

#pindel
pindel_path=het+'pbs_Pindel/'
pindel_path_1=pindel_path+'pbs_Pindel_step1/'
pindel_path_2=pindel_path+'pbs_Pindel_step2/'
run_time_hash['pindel']={}
for k1 in os.listdir(pindel_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['pindel'][chromo_name]=[read_in_run_time(pindel_path_1+k1)]

for k1 in os.listdir(pindel_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-2]
        if chromo_name in run_time_hash['pindel'].keys():
            run_time_hash['pindel'][chromo_name]+=[read_in_run_time(pindel_path_2+k1)]

#erds
erds_path=het+'pbs_erds/'
erds_path_1=erds_path+'pbs_erders_step3_erds/'
run_time_hash['erds']={}
for k1 in os.listdir(erds_path_1):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-5]
        run_time_hash['erds'][chromo_name]=[read_in_run_time(erds_path_1+k1)]

erds_path_2=erds_path+'pbs_erders_step2_SNPs_calling/'
run_time_hash['erds_SNP']={}
for k1 in os.listdir(erds_path_2):
     if k1.split('.')[-1]=='log':
        chromo_name=k1.split('.')[-4]
        run_time_hash['erds_SNP'][chromo_name]=[read_in_run_time(erds_path_2+k1)]

algorithms=run_time_hash.keys()
rds=sorted(run_time_hash['delly'].keys())

fo=open(ppre+'Run_time_Sumi','a')
print >>fo,' '.join(['comp.homo']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([sum(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo,' '.join([k1]+[str(i) for i in test])

print >>fo,' '.join(['comp.homo_Paralelle']+algorithms)
for k1 in rds:
    test=[]
    for k2 in algorithms:
        if k1 in run_time_hash[k2].keys():
            test.append(sum([max(i) for i in run_time_hash[k2][k1]])/3600)
        else:
            test.append('*')
    print >>fo, ' '.join([k1]+[str(i) for i in test])

fo.close()

