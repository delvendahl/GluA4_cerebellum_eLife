# Get population sparseness, variance & covariance

import numpy as np
import os


def get_spar(x):
	N = x.shape[0]
	T = x.shape[1]
	sptemp = np.zeros(T)
	sptemp2 = np.zeros(T)
	for t in range(T):
		sptemp[t] = (N-np.sum(x[:, t])**2. / np.sum(x[:, t]**2.)) / (N - 1.)
		sptemp2[t] = len(np.where(x[:, t] > 0)[0]) * 1. / N
	spar = np.nanmean(sptemp)
	active = np.nanmean(sptemp2)
	return spar, active


def get_var_cov(x):
	N = x.shape[0]
	L, V = np.linalg.eig(np.cov(x))
	L = np.real(np.sqrt(L+0J))
	var_x = np.sum(L**2)
	cov_x = (np.max(L)/np.sum(L) - 1./N)/(1.-1./N)
	return var_x, cov_x


if __name__ == '__main__':
	# basedir = 'results/_freq/data_150_r20_orig/r20_10ms'
	desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')
	basedir = desktop + '/simulation november 2020/RESULTS/Heterogeneous/r30_dat_75Hz' + '/'

	N_syn = 4
	f_mf = np.linspace(0.1, 0.9, 5)

	spar_mf, spar_grc, active_mf, active_grc = (np.zeros(len(f_mf), float) for i in range(4))
	var_mf, cov_mf, var_grc, cov_grc = (np.zeros(len(f_mf), float) for i in range(4))

	print(N_syn)
	for k in range(len(f_mf)):
		samples_mf = np.loadtxt(basedir + 'MF_samples_{:.0f}_{:.2f}.txt'.format(N_syn, f_mf[k]))
		samples_grc = np.loadtxt(basedir + 'GrC_samples_{:.0f}_{:.2f}.txt'.format(N_syn, f_mf[k]))
		spar_mf[k], active_mf[k] = get_spar(samples_mf)
		spar_grc[k], active_grc[k] = get_spar(samples_grc)
		var_mf[k], cov_mf[k] = get_var_cov(samples_mf)
		var_grc[k], cov_grc[k] = get_var_cov(samples_grc)

	filename1 = 'grc_spar_biophys_' + basedir.split('_')[-1][:-1] + '.txt'
	np.savetxt(basedir + filename1, np.transpose([spar_mf, spar_grc, active_mf, active_grc]), delimiter='\t')

	filename2 = 'grc_cov_biophys_' + basedir.split('_')[-1][:-1] + '.txt'
	np.savetxt(basedir + filename2, np.transpose([var_mf, var_grc, cov_mf, cov_grc]), delimiter='\t')

	print(np.transpose([spar_mf, spar_grc, var_mf, var_grc, cov_mf, cov_grc]))
