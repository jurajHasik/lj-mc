#!/usr/bin/python3
import sys
import math

print('Opening'+sys.argv[1])
f = open(sys.argv[1],'r') #First arg is filename
eql = 0
prl = 0

eAvg = 0
pAvg = 0
rhoAvg = 0
Elst = []
Plst = []
rholst = []
for line in f:
	fields = line.split()
	if fields[0]=='rT':
		print('rT: '+fields[1])
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
		rhoAvg = rhoAvg + float(fields[8])
		rholst.append(float(fields[8]))
print('Num of samples: '+str(prl))
eAvg = eAvg / prl
pAvg = pAvg / prl
rhoAvg = rhoAvg / prl
print('E_smpl_avg: {0:.5f}'.format(eAvg))
print('P_smpl_avg: {0:.5f}'.format(pAvg))
print('rho_smpl_avg: {0:.5f}'.format(rhoAvg))

#Start binning
print('Bin_size E_AvgEr P_AvgEr rho_AvgEr'
	'Bin_size^1/2 E_AvgEr P_AvgEr rho_AvgEr')
bb = int(sys.argv[2]) # Largest bin size
for bin_size in range(1,bb+1):
	# We will loop only through whole number bb-tuples
	# discarding the rest
	E_b = []
	P_b = []
	rho_b = []
	d_point = iter(Elst)
	p_point = iter(Plst)
	rho_point = iter(rholst)
	bin_num = prl//bin_size
	for bin_no in range(bin_num):
		E_bin_avg = 0.0
		P_bin_avg = 0.0
		rho_bin_avg = 0.0
		for d in range(bin_size):
			E_bin_avg += next(d_point)
			P_bin_avg += next(p_point)
			rho_bin_avg += next(rho_point)
		E_bin_avg = E_bin_avg/bin_size
		P_bin_avg = P_bin_avg/bin_size
		rho_bin_avg = rho_bin_avg/bin_size
		E_b.append(E_bin_avg)
		P_b.append(P_bin_avg)
		rho_b.append(rho_bin_avg)
	E_tot_var = 0.0
	P_tot_var = 0.0
	rho_tot_var = 0.0
	for bin_no in range(bin_num):
		E_tot_var += (E_b[bin_no]-eAvg)**2.0
		P_tot_var += (P_b[bin_no]-pAvg)**2.0
		rho_tot_var += (rho_b[bin_no]-rhoAvg)**2.0
	E_tot_var = E_tot_var/bin_num
	P_tot_var = P_tot_var/bin_num
	rho_tot_var = rho_tot_var/bin_num
	if(bin_size == 1):
		E_smpl_var = E_tot_var
		P_smpl_var = P_tot_var
		rho_smpl_var = rho_tot_var
	E_mean_var = math.sqrt(E_tot_var/bin_num)
	P_mean_var = math.sqrt(P_tot_var/bin_num)
	rho_mean_var = math.sqrt(rho_tot_var/bin_num)
	print('{0:05d} {1:.5f} {2:.5f} {3:.5f} {4:.5f} {5:.5f}'
		' {6:.5f} {7:.5f}'.format(bin_size, E_mean_var, P_mean_var, \
		rho_mean_var, math.sqrt(bin_size), bin_size*E_tot_var/E_smpl_var, \
		bin_size*P_tot_var/P_smpl_var, bin_size*rho_tot_var/rho_smpl_var))


