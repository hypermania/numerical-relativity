import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos
from scipy import loadtxt, optimize
import os

M = 1.3669

plt.figure(figsize=(10,7), dpi=80)
ax = plt.axes()

dat = loadtxt("./mat_tracers/mat_tracers.tsv", skiprows=0, delimiter="\t")
t = dat.transpose()[0]
tracers = dat.transpose()[1:]

r_max = tracers.max()

ax.axis([0, r_max * 1.1 , 0, t.max() * 1.0])


for i in range(0, len(tracers)):
    ax.plot(tracers[i] / M , t, label = r"m/M = {}".format((i+1)*0.2), linestyle='--')

ax.plot(horizon_plot / M, t_plot, label="Event Horizon")
    
plot_name = r"Lagrangian matter tracers"
ax.set_title(plot_name)
ax.set_xlabel("r/M")
ax.set_ylabel("t")
ax.legend()

#plt.savefig("Tracer.jpg")
plt.show()


