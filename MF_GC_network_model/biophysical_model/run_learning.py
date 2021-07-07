# Backpropagation learning

import numpy as np
from datetime import datetime
from multiprocessing import Pool
from functools import partial
import sys


def get_learning_speed(loss_curve):
    cutoff = 0.2   # cutoff value for learning speed estimation
    temp = np.where(np.array(loss_curve) <= cutoff)
    speed = 1. / (temp[0][0] + 1) if len(temp[0]) > 0 else 0

    return speed


def backprop_step_nohid(w_out, gamma, training_pattern, target_pattern):
    # Dynamics of units
    s = lambda x: 1. / (1. + np.exp(-x))
    ds = lambda s: s * (1. - s)  # sigmoidal
    #####################################
    # First step: feedforward propagation
    o_in = training_pattern
    o_out = s(np.dot(np.append(o_in, 1), w_out))
    # Second step: output layer backpropagation
    diag = np.diag(ds(o_out))
    err = o_out - target_pattern
    err_d = np.prod((target_pattern == target_pattern.max()) == (o_out == o_out.max()))
    delta_out = np.dot(diag, err)
    dw_out = - gamma * np.outer(np.append(o_in, 1), delta_out)
    # Third step: update weights
    w_out = w_out + dw_out
    
    return err, err_d, w_out


def backprop_nohid(training_set, target, n_epochs, gamma):
    
    n_in = training_set.shape[0]
    n_out = target.shape[0]
    w_out = np.random.uniform(-1., 1., size=(n_in + 1, n_out)) * 1. / (n_in + 1)
    # Shuffle order of training set
    temp = range(training_set.shape[1])
    np.random.shuffle(temp)
    training_set = training_set[:, temp]
    target = target[:, temp]
    
    errors_rms = np.zeros(n_epochs, float)
    errors_discrim = np.zeros(n_epochs, float)
    for ep in range(n_epochs):
        errors_temp = np.zeros(target.shape[1], float)
        errors_d_temp = np.zeros(target.shape[1], float)
        for k in range(target.shape[1]):
            # Backpropagation step
            err, err_d, w_out = backprop_step_nohid(w_out, gamma, training_set[:, k], target[:, k])
            # Record errors
            errors_temp[k] = np.sqrt(np.mean(err ** 2))  # RMS error
            errors_d_temp[k] = err_d  # Discrimination error
        # Record average error for the epoch
        errors_rms[ep] = errors_temp.mean()
        errors_discrim[ep] = errors_d_temp.mean()
        # Reshuffle order of training data
        temp = range(training_set.shape[1])
        np.random.shuffle(temp)
        training_set = training_set[:, temp]
        target = target[:, temp]
    
    return errors_rms, errors_discrim, w_out


def analyse_learning(basedir, f_mf):
    
    # Network parameters
    n_syn = 4
    # Backprop parameters
    n_epochs = 5000
    gamma = 0.01
    c = 10
    num_patterns = 64 * c
    
    err_mf = np.zeros(n_epochs, float)
    err_rms_mf = np.zeros(n_epochs, float)
    err_grc = np.zeros(n_epochs, float)
    err_rms_grc = np.zeros(n_epochs, float)
    
    samples_mf = np.loadtxt(basedir + '/MF_samples_' + str(n_syn) + '_' + '{:.2f}'.format(f_mf) + '.txt')
    samples_grc = np.loadtxt(basedir + '/GrC_samples_' + str(n_syn) + '_' + '{:.2f}'.format(f_mf) + '.txt')
    # Get pattern classifications
    target = np.zeros((c, num_patterns))
    for k in range(num_patterns):
        target[np.random.choice(c), k] = 1
    # Single layer backpropagation
    err_rms_mf[:], err_mf[:], w_mf = backprop_nohid(samples_mf, target, n_epochs, gamma)
    err_rms_grc[:], err_grc[:], w_grc = backprop_nohid(samples_grc, target, n_epochs, gamma)
    
    # File name to save
    filename = 'grc_bp_biophys_' + '{:.2f}'.format(f_mf) + '_' + basedir.split('_')[-1][:-1] + '.txt'
    # Save results
    np.savetxt(basedir + filename, np.transpose([err_rms_mf, err_rms_grc, err_mf, err_grc]), delimiter='\t')
    
    mf_learning_speed = get_learning_speed(err_rms_mf)
    grc_learning_speed = get_learning_speed(err_rms_grc)
    
    print(basedir, f_mf, mf_learning_speed, grc_learning_speed, err_rms_mf[-1], err_rms_grc[-1])
    
    return mf_learning_speed, grc_learning_speed, err_rms_mf[-1], err_rms_grc[-1]


if __name__ == '__main__':
    
    startTime = datetime.now()

    basedir = sys.argv[1]
    print(basedir)
    
    f_mf = np.linspace(0.1, 0.9, 9)

    pool = Pool()
    
    fun = partial(analyse_learning,
                  basedir)
    results = pool.map(fun, f_mf)
    
    pool.close()
    pool.join()

    print('')
    np.savetxt(basedir + '/learning_results.txt', results, delimiter='\t')
    print(datetime.now() - startTime)
