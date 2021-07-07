import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pyspike as spk
import numpy as np
import itertools
import os


def spike_sync(inputfilename):
    """
    Calculate average spike synchronicity using pyspike
    """
    spike_trains = spk.load_spike_trains_from_txt(inputfilename, separator='\t', edges=(100, 1100))
    spike_sync = spk.spike_sync_matrix(spike_trains, interval=(100, 1100))

    num_trains = np.array(spike_trains).shape[0]
    sync_vals = np.zeros(num_trains/2)

    for i, j in zip(xrange(0, num_trains, 2), xrange(num_trains/2)):
        sync_vals[j] = spike_sync[i, i+1]

    poisson_trains = np.array(xrange(0, num_trains, 2))
    combos = list(itertools.combinations(poisson_trains, 2))
    poisson_sync = np.zeros(len(combos))

    for i in xrange(len(combos)):
        poisson_sync[i] = spike_sync[combos[i][0], combos[i][1]]

    return np.average(sync_vals), np.average(poisson_sync)


def spike_dist(inputfilename):
    """
    Calculate average spike distance using pyspike
    """
    spike_trains = spk.load_spike_trains_from_txt(inputfilename, separator='\t', edges=(100, 1100))
    spike_dist = spk.spike_distance_matrix(spike_trains, interval=(100, 1100))

    num_trains = np.array(spike_trains).shape[0]
    dist_vals = np.zeros(num_trains/2)

    for i, j in zip(xrange(0, num_trains, 2), xrange(num_trains/2)):
        dist_vals[j] = spike_dist[i, i+1]

    poisson_trains = np.array(xrange(0, num_trains, 2))
    combos = list(itertools.combinations(poisson_trains, 2))
    poisson_dist = np.zeros(len(combos))

    for i in xrange(len(combos)):
        poisson_dist[i] = spike_dist[combos[i][0], combos[i][1]]

    return np.average(dist_vals), np.average(poisson_dist)


def get_field_sub(x): return x.split('_')[0]


if __name__ == '__main__':

    files = list()
    for filename in os.listdir('./SpikeTimes'):
        if filename.endswith('.txt'):
            files.append(filename)
            continue
        else:
            continue

    files.sort()
    mylist = sorted(files, key=get_field_sub)
    grouped_files = [list(y) for x, y in itertools.groupby(mylist, get_field_sub)]

    print [x[0:1] for x in grouped_files]
    num_freq = len(grouped_files[0])
    num_conditions = len(grouped_files)

    spk_sync_results = np.zeros((num_conditions, num_freq))
    poisson_sync_results = np.zeros((num_conditions, num_freq))
    spk_dist_results = np.zeros((num_conditions, num_freq))
    poisson_dist_results = np.zeros((num_conditions, num_freq))

    for i in xrange(num_freq):
        for j in xrange(num_conditions):
            spk_sync_results[j, i], poisson_sync_results[j, i] = spike_sync('./SpikeTimes/' + grouped_files[j][i])
            spk_dist_results[j, i], poisson_dist_results[j, i] = spike_dist('./SpikeTimes/' + grouped_files[j][i])

    frequencies = np.arange(1, num_freq+1) * 50

    fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(6, 5))

    for i in xrange(num_conditions):
        axs[0, 0].plot(frequencies, spk_sync_results[i])
    axs[0, 0].set_xlabel('Hz')
    axs[0, 0].set_ylabel('spike_sync')

    for i in xrange(num_conditions):
        axs[0, 1].plot(frequencies, poisson_sync_results[i])
    axs[0, 1].set_xlabel('Hz')
    axs[0, 1].set_ylabel('spike_sync')

    for i in xrange(num_conditions):
        axs[1, 0].plot(frequencies, spk_dist_results[i])
    axs[1, 0].set_xlabel('Hz')
    axs[1, 0].set_ylabel('spike_dist')

    for i in xrange(num_conditions):
        axs[1, 1].plot(frequencies, poisson_dist_results[i])
    axs[1, 1].set_xlabel('Hz')
    axs[1, 1].set_ylabel('spike_dist')

    plt.show()

    results = np.concatenate(([frequencies], spk_sync_results, poisson_sync_results,
                              spk_dist_results, poisson_dist_results), axis=0)
    results = np.transpose(results)
    np.savetxt('results_all.txt', results, fmt='%1.8e', delimiter='\t')
