#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.colors
from matplotlib import cm
import sys

input_file = sys.argv[1]
forces = []
tensions = []

with open(input_file) as fh:
    for line in fh:
        split_line = line.split()
        if len(split_line) == 5:
            # This is a member with tension
            tt = tuple([float(s) for s in split_line])
            tensions.append(tt)
        elif len(split_line) == 4:
            # This is a reaction force to the truss
            tt= tuple([float(s) for s in split_line])
            forces.append(tt)
        else:
            raise Exception('error reading output file')

# pos1x, pos1y, pos2x, pos2y, tension
# posx, posy, f_x, f_y
min_force = min([x[4] for x in tensions])
max_force = max([x[4] for x in tensions])
norm = matplotlib.colors.DivergingNorm(vmin=min_force, vcenter=0, vmax=max_force, )
cmap = cm.get_cmap('seismic', 256)

for tension in tensions:
    plt.plot([tension[0], tension[2]], [tension[1],tension[3]], c=cmap(norm(tension[4])))
    x_mid = (tension[0]+ tension[2])/2
    y_mid = (tension[1]+ tension[3])/2
    try:
        angle = np.arctan((tension[3]-tension[1])/(tension[2]-tension[0]))
        angle *= 180 / np.pi
        if angle > 180:
            angle -= 360
    except:
        angle =0
    plt.text(x_mid, y_mid, '%5.2f'%tension[4], horizontalalignment='center', verticalalignment='bottom', rotation=angle, rotation_mode='anchor')
plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
plt.savefig(sys.argv[1]+'.png')
# plt.show()
