#!/usr/bin/env python3
import os,sys
import numpy as np
import matplotlib.pyplot as plt

if os.system("make main") != 0:
  sys.exit("Failed to compile main.cc!")

L = 64
T = 2.26919
steps = 10000
hlist = np.concatenate([ 
  0.0001*np.arange(1,10), 
  0.001*np.arange(1,10),
  0.01*np.arange(1,10),
  0.1*np.arange(1,10),
  [1.0]
  ])
hlist = np.round(hlist, 12)

codeDir = os.getcwd()

for h in hlist:
  dataDir = "DATA/h{}".format(h, 'g')
  print ("Creating directory ", dataDir)
  os.chdir (codeDir)
  os.makedirs (dataDir, mode=0o755, exist_ok=True)
  os.chdir (dataDir)
  
  with open("pars.dat", "w") as f:
    f.write ("xmax        = {}\n".format(L))
    f.write ("ymax        = {}\n".format(L))
    f.write ("temperature = {}\n".format(T))
    f.write ("field       = {}\n".format(h))
    f.write ("steps       = {}\n".format(steps))

  os.system (codeDir+"/main")
