import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos
from scipy import loadtxt, optimize
import os

M = 1.3669

plt.figure(figsize=(10,7), dpi=80)
ax = plt.axes()

t_min = 2000

t_plot = []
horizon_plot = []

filename_list = sorted(os.listdir("rays"), key=lambda a: float(a[:-4][4:]))
for filename in filename_list:
    dat = loadtxt("./rays/{}".format(filename), skiprows=0, delimiter="\t")
    t = dat.transpose()[0]
    rays = dat.transpose()[1:]

    #horizon_idx = max(np.argmax(np.diff(rays[:, -1]))-1, 0)
    #horizon_idx_2 = np.argmax(np.diff(rays[:, -1])) + 2
    #horizon_plot.append(rays[:,0][horizon_idx])
    for i in range(0, len(rays)):
    #for i in range(max(horizon_idx-4, 0), horizon_idx_2):
        ax.plot(rays[i], t)

    idx = rays[:,-1] < 25
    if idx.any():
        horizon_plot.append((rays[:,0][rays[:,-1]<25]).max())
        t_plot.append(t.min())        
    t_min = min(t_min, t.min())

r_max = rays.max()

t_plot = np.array(t_plot)
horizon_plot = np.array(horizon_plot)


ax.axis([0, horizon_plot.max() * 1.2, t_min, t_plot.max() * 1.1])
ax.plot(horizon_plot / M, t_plot, label="Event Horizon", linewidth=4)

plot_name = r"Event Horizon"
ax.set_title(plot_name)
ax.set_xlabel("r")
ax.set_ylabel("t")
ax.legend()


plt.show()

