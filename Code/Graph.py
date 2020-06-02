import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos
from scipy import loadtxt, optimize
import os

filename_list = sorted(os.listdir("solution"), key=lambda a: float(a[:-4]))

# Compute ADM mass of system
r, f = loadtxt("./solution/{}".format(filename_list[0]), unpack=True, skiprows=0, delimiter="\t")
M = (f[int(len(f)*0.15)] - 1) * 2 * r[int(len(f)*0.15)]
#M_adm = 0.5 * f[-2]**2 * r[-2] * (1 - (1 + (2 * r[-2] / f[-2]) * (f[-1] - f[-3]) / (r[-1] - r[-3]))**2)

plt.figure(figsize=(10,7), dpi=80)
ax = plt.axes()

f_max = 0
for filename in filename_list:
    path = "./solution/{}".format(filename)
    r, f = loadtxt(path, unpack=True, skiprows=0, delimiter="\t")
    t = filename[:-4]
    ax.plot(r / M * f**2, f**2, label = 't = {}'.format(t))
    #ax.plot(((M+2*r)**2 / (4*r))[40:] / M, (f**2)[40:], label = 't = {}'.format(t))
    f_max = max(f.max(), f_max)


ax.axis([0.1, r.max() * 1.1 / M, 1, (f_max**2) * 1.1])

plot_name = r"Metric component $A$ versus areal radius"

ax.set_title(plot_name)
ax.set_xlabel("R/M")
ax.set_ylabel("A")
ax.set_xscale("log")
ax.legend()
#ax.plot(r, f, label = 't = {}'.format(0))


#r_vacuum = np.linspace(0.1, r.max() * 5, 1000)
#f_vacuum = 1 + M / (2 * r_vacuum)

#ax.plot(r_vacuum / M, f_vacuum**2, label = 'vacuum')

#save_name = "plot_3.pdf"
#plt.savefig(save_name)

#plt.savefig("A.jpg")
plt.show()
