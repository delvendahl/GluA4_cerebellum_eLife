import numpy as np
import pickle as pkl
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


file = open('GCLconnectivity_4.pkl')
p = pkl.load(file)
file.close()
conn_mat = p['conn_mat']
grc_pos = p['grc_pos']
glom_pos = p['glom_pos']
N_mf, N_grc = conn_mat.shape
assert (np.all(conn_mat.sum(axis=0) == 4)), 'Connectivity matrix is incorrect.'


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k1 in range(N_mf):
    for k2 in range(N_grc):
        if conn_mat[k1, k2] == 1:
            ax.plot([glom_pos[k1, 0], grc_pos[k2, 0]],
                    [glom_pos[k1, 1], grc_pos[k2, 1]],
                    [glom_pos[k1, 2], grc_pos[k2, 2]],
                    color='grey', alpha=0.6)

ax.plot(grc_pos[:, 0], grc_pos[:, 1], grc_pos[:, 2], c='b', marker='o', linestyle='')
ax.plot(glom_pos[:, 0], glom_pos[:, 1], glom_pos[:, 2], c='r', marker='o', linestyle='')

plt.show()

