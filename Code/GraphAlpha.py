import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos
from scipy import loadtxt, optimize
import os

M = 1.3669

plt.figure(figsize=(10,7), dpi=80)
ax = plt.axes()

alpha_max = 0
alpha_min = 1

alpha_filename_list = sorted(os.listdir("alpha"), key=lambda a: float(a[:-4]))
#f_filename_list = sorted(os.listdir("solution"), key=lambda a: float(a[:-4]))

# Compute ADM mass of system
r, f = loadtxt("./solution/{}".format(alpha_filename_list[0]), unpack=True, skiprows=0, delimiter="\t")
M = (f[int(len(f)*0.1)] - 1) * 2 * r[int(len(f)*0.1)]

for filename in alpha_filename_list:
    path_alpha = "./alpha/{}".format(filename)
    path_f = "./solution/{}".format(filename)
    r, alpha = loadtxt(path_alpha, unpack=True, skiprows=0, delimiter="\t")
    r, f = loadtxt(path_f, unpack=True, skiprows=0, delimiter="\t")
    t = filename[:-4]
    ax.plot(f*f*r/M, alpha, label = 't = {}'.format(t))
    alpha_max = max(alpha.max(), alpha_max)
    alpha_min = min(alpha.min(), alpha_min)

ax.axis([0.1 / M, r.max() * 1.1 / M, alpha_min * 0.5, alpha_max * 1.1])

plot_name = r"Time lapse $\alpha$ versus areal radius"

ax.set_title(plot_name)
ax.set_xlabel("R/M")
ax.set_ylabel(r'$\alpha$')
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

plt.savefig("alpha.jpg")
plt.show()

