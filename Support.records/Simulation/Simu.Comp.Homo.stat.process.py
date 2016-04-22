import os
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='stat' and not k1=='BPLink.Integrated.Categorized.stat':
		filein=k1
		fileout=k1+'.new'
		fin=open(filein)
		fo=open(fileout,'w')
		for line in fin:
			pin=line.strip().split()
			pin+=pin[4].split('_')
			print >>fo, ' '.join(pin)
		fin.close()
		fo.close()

import os
CSV_type_1=['ab/ab_aba^/aba^','ab/ab_b^ab/b^ab','ab/ab_b^a^b/b^a^b']
CSV_type_2=['abc/abc_b^/b^','abc/abc_ac^/ac^','abc/abc_a^c/a^c','abc/abc_a^c^/a^c^','abc/abc_c^a^/c^a^']
CSV_type_3=['abc/abc_aba/aba','abc/abc_cbc/cbc']
CSV_type_4=['abc/abc_aba^/aba^','abc/abc_c^bc/c^bc','abc/abc_aca^/aca^']
for k1 in os.listdir('./'):
	if k1.split('.')[-1]=='stat' and not k1=='BPLink.Integrated.Categorized.stat':
		filein=k1
		fileout1=fileout=k1+'.SV1'
		fileout2=fileout=k1+'.SV2'
		fileout3=fileout=k1+'.SV3'
		fileout4=fileout=k1+'.SV4'
		fin=open(filein)
		fo1=open(fileout1,'w')
		fo2=open(fileout2,'w')
		fo3=open(fileout3,'w')
		fo4=open(fileout4,'w')
		for line in fin:
			pin=line.strip().split()
			if not pin==[]:
				if '_'.join(pin[4].split('_')[:2]) in CSV_type_1:
					print >>fo1, ','.join(pin)
				elif '_'.join(pin[4].split('_')[:2]) in CSV_type_2:
					print >>fo2, ','.join(pin)
				elif '_'.join(pin[4].split('_')[:2]) in CSV_type_3:
					print >>fo3, ','.join(pin)
				elif '_'.join(pin[4].split('_')[:2]) in CSV_type_4:
					print >>fo4, ','.join(pin)
		fin.close()
		fo1.close()
		fo2.close()
		fo3.close()
		fo4.close()

data_hash={}
for k1 in os.listdir('./'):
	if 'SV' in k1.split('.')[-1]:
		key_algorithm=k1.split('.')[0]
		key_RD=k1.split('.')[1]
		key_SV=k1.split('.')[-1]
		if not key_algorithm in data_hash.keys():
			data_hash[key_algorithm]={}
		if not key_RD in data_hash[key_algorithm].keys():
			data_hash[key_algorithm][key_RD]={}
		if not key_SV in data_hash[key_algorithm][key_RD].keys():
			data_hash[key_algorithm][key_RD][key_SV]=[]
		fin=open(k1)
		FP=0
		TotalP=0
		TotalPredict=0
		for line in fin:
			pin=line.strip().split(',')
			FP+=int(pin[0])
			TotalP+=int(pin[1])
			TotalPredict+=int(pin[2])
		fin.close()
		data_hash[key_algorithm][key_RD][key_SV]=[float(FP)/float(TotalP),float(FP)/float(TotalPredict),FP,TotalP,TotalPredict]

fo=open('Integrated.CSV.stat','w')
for k1 in sorted(data_hash.keys()):
	for k2 in sorted(data_hash[k1].keys()):
		for k3 in sorted(data_hash[k1][k2].keys()):
			print >>fo, ' '.join([str(i) for i in [k1,k2,k3]+data_hash[k1][k2][k3]])
