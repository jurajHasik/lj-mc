#!/usr/bin/python3
import sys
import math
import random

print('Opening'+sys.argv[1])
f = open(sys.argv[1],'r') #First arg is filename
eql = 0
prl = 0

rT = 0.0
eAvg = 0
pAvg = 0
Elst = []
Plst = []
for line in f:
	fields = line.split()
	if fields[0]=='rT':
		print('rT: '+fields[1])
		rT = float(fields[1])
	if fields[0]=='rDens':
		print('rDens: '+fields[1])
	if fields[0]=='E_lrc:':
		print(line)
	if fields[0]=='EQ':
		eql=eql+1
	if fields[0]=='PR':
		prl = prl+1
		eAvg = eAvg + float(fields[4])
		Elst.append(float(fields[4]))
		pAvg = pAvg + float(fields[6]) 		
		Plst.append(float(fields[6]))
print('Num of samples: '+str(prl))
eAvg = eAvg / prl
eVar = 0.0
for i in range(prl):
	eVar = eVar + (Elst[i]-eAvg)**2.0
print('full sample CV: {0:.5f}'.format((eVar/len(Elst))/(rT*rT)))	

# Draw randomly (lets say) thousand values and compute
# C_V, repeat (let say) thousand times and take average
# of C_V sample + error of mean of the C_V sample
CV = []
CV_avg = 0.0
CV_var = 0.0
CV_avg_err = 0.0
bin_s = int(sys.argv[2])
cv_smpl = 1000
for i in range(cv_smpl):
	E_draw = []
	E_avg = 0.0
	E_var = 0.0
	for j in range(bin_s):
		ind = random.randint(0,len(Elst)-1)
		E_draw.append(Elst[ind])
		E_avg = E_avg + Elst[ind]
	E_avg = E_avg / bin_s
	for j in range(bin_s):
		E_var = E_var + (E_draw[j]-E_avg)**2.0
	E_var = E_var / bin_s
	CV.append(E_var/(rT*rT))
	CV_avg = CV_avg + E_var/(rT*rT)
	print('C_V: {0:.5f} +/- {1:.5f}'.format(E_var,0.0))
CV_avg = CV_avg / cv_smpl
for i in range(cv_smpl):
	CV_var = CV_var + (CV[i]-CV_avg)**2.0
CV_avg_err = math.sqrt(CV_var / cv_smpl)
print('CV_var: {0:.5f}'.format(CV_var))
print('CV_avg: {0:.5f} +/- {1:.5f}'.format(CV_avg,CV_avg_err))