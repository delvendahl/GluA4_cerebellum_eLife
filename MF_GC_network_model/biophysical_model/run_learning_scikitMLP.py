import numpy as np
from sklearn.neural_network import MLPClassifier
from datetime import datetime
import sys


def get_learning_speed(loss_curve):
    cutoff = 1.5    # cutoff value for learning speed estimation
    temp = np.where(np.array(loss_curve) <= cutoff)
    speed = 1. / (temp[0][0] + 1) if len(temp[0]) > 0 else 0

    return speed


def run_learning(samples, target):

    clf = MLPClassifier(solver='sgd',
                        alpha=0.0001,
                        random_state=1,
                        max_iter=5000,
                        learning_rate='adaptive',
                        early_stopping=False,
                        hidden_layer_sizes=(samples.shape[1],))

    clf.fit(samples, target)
    score = clf.score(samples, target)
    learning_speed = get_learning_speed(clf.loss_curve_) if hasattr(clf, "loss_curve_") else np.nan

    return score, clf.loss_, learning_speed, clf.n_iter_


if __name__ == '__main__':

    startTime = datetime.now()
    np.random.seed(42)

    basedir = sys.argv[1]
    print(basedir)

    n_syn = 4
    n_classes = 10
    f_mf = np.linspace(0.1, 0.9, 9)

    mf_results = np.zeros((len(f_mf), 4))
    gc_results = np.zeros((len(f_mf), 4))

    for i, fraction in enumerate(f_mf):
        mf_samples = np.loadtxt(basedir + '/MF_samples_{}_{:.2f}.txt'.format(n_syn, fraction)).T
        gc_samples = np.loadtxt(basedir + '/GrC_samples_{}_{:.2f}.txt'.format(n_syn, fraction)).T

        y = np.random.choice(n_classes, mf_samples.shape[0], replace=True)

        mf_results[i] = run_learning(samples=mf_samples, target=y)
        gc_results[i] = run_learning(samples=gc_samples, target=y)

    print('MF patterns:')
    print(mf_results)
    print('')
    print('GC patterns:')
    print(gc_results)

    np.savetxt(basedir + '/sklearn_result.txt', np.hstack((mf_results, gc_results)), delimiter='\t')
    print(datetime.now() - startTime)
