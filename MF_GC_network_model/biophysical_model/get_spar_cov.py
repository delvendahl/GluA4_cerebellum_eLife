# Get population sparseness, variance & covariance

import numpy as np
import sys


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

	basedir = '/' + sys.argv[1] + '/'

	N_syn = 4
	f_mf = np.linspace(0.1, 0.9, 9)

	spar_mf, spar_gc, active_mf, active_gc = (np.zeros(len(f_mf), float) for i in range(4))
	var_mf, cov_mf, var_gc, cov_gc = (np.zeros(len(f_mf), float) for i in range(4))

	print(N_syn)
	for k in range(len(f_mf)):
		samples_mf = np.loadtxt(basedir + 'MF_samples_{:.0f}_{:.2f}.txt'.format(N_syn, f_mf[k]))
		samples_gc = np.loadtxt(basedir + 'GC_samples_{:.0f}_{:.2f}.txt'.format(N_syn, f_mf[k]))
		spar_mf[k], active_mf[k] = get_spar(samples_mf)
		spar_gc[k], active_gc[k] = get_spar(samples_gc)
		var_mf[k], cov_mf[k] = get_var_cov(samples_mf)
		var_gc[k], cov_gc[k] = get_var_cov(samples_gc)

	filename1 = 'gc_spar_biophys_' + basedir.split('_')[-1][:-1] + '.txt'
	np.savetxt(basedir + filename1, np.transpose([spar_mf, spar_gc, active_mf, active_gc]), delimiter='\t')

	filename2 = 'gc_cov_biophys_' + basedir.split('_')[-1][:-1] + '.txt'
	np.savetxt(basedir + filename2, np.transpose([var_mf, var_gc, cov_mf, cov_gc]), delimiter='\t')

	print(np.transpose([spar_mf, spar_gc, var_mf, var_gc, cov_mf, cov_gc]))
