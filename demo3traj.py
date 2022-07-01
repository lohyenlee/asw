#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

if os.system("make main") != 0:
  sys.exit("Failed to compile main.cc!")

L = 64
T = 2.27
h = 0.003
steps = 10000

with open("pars.dat", "w") as f:
    f.write ("xmax        = {}\n".format(L))
    f.write ("ymax        = {}\n".format(L))
    f.write ("temperature = {}\n".format(T))
    f.write ("field       = {}\n".format(h))
    f.write ("steps       = {}\n".format(steps))

os.system ("./main")

data = np.loadtxt("history.dat")
M = data[:,0]
U = data[:,1]
 
plt.plot(M/L**2, label="Magnetization $M/N$")  
plt.plot(U/L**2, label="Energy per site $U/N$")  
plt.xlabel("Simulation time t (MCS")
plt.title("Monte Carlo trajectories")
plt.legend(loc="upper right")  
plt.show()
