#!/usr/bin/env python3
import os,sys,re
import numpy as np
import matplotlib.pyplot as plt

#======== ESTIMATE ERROR FROM CORRELATED DATA
def estimateError (timeSeries, numBlocks=10):
  timeSeries = np.array (timeSeries)
  n = timeSeries.size
  if (n < numBlocks*2):
    return "NOT ENOUGH DATA POINTS TO PARTITION INTO numBlocks BLOCKS!"
  blocks = np.split(timeSeries, np.floor(np.arange(n/numBlocks,n,n/numBlocks)).astype(int))
  blockVars = []
  for block in blocks:
    blockMean = np.mean (block)
    blockVar = np.var (block, ddof=1)
    blockVars += [blockVar]
  meanBlockVar = np.mean (blockVars)
  return np.sqrt(meanBlockVar/numBlocks)

#======== READ MONTE CARLO PARAMETERS AND TRAJECTORIES FROM DATA DIRECTORIES
codeDir = os.getcwd()
dataDirList = sorted ( os.listdir("DATA") )
data = np.empty ([0,5])       # prepare to add stats lines of 5 elements each
for dataDir in dataDirList:
  dataPath = "DATA/" + dataDir

  #-------- READ PARAMETERS h, xmax, ymax
  with open(dataPath + "/pars.dat") as f:
    s = f.read()
  match = re.search ('field\s*=\s*(.*)', s); h = float (match.group(1))
  match = re.search ('xmax\s*=\s*(.*)', s);  xmax = float (match.group(1))
  match = re.search ('ymax\s*=\s*(.*)', s);  ymax = float (match.group(1))
  
  #-------- READ TRAJECTORIES M(t), U(t)
  trajs = np.loadtxt(dataPath + "/history.dat")
  M = trajs[5000:,0]   # discard first 1000 MCS
  Mavg = np.mean(M).item() / (xmax*ymax)
  Merr = estimateError (M) / (xmax*ymax)
  
  U = trajs[5000:,1]
  Uavg = np.mean(U).item() / (xmax*ymax)
  Uerr = estimateError (U) / (xmax*ymax)
  
  print ("dataDir={:8s}  h={:8g}   Mavg={:8.4f}   Merr={:8.4f}".
          format (dataDir, h, Mavg, Merr, Uavg, Uerr))

  data = np.append (data, [[h, Mavg, Merr, Uavg, Uerr]], axis=0)
      
print ("data.shape = ", data.shape)
[h,Mavg,Merr,Uavg,Uerr] = np.transpose(data)

#======== PLOT DATA
plt.errorbar (h,Mavg,Merr,linewidth=1,capsize=3,label="Magnetic moment M")
plt.errorbar (h,Uavg,Uerr,linewidth=1,capsize=3,label="Energy U")
plt.xscale("log")
plt.xlabel("Field h")
plt.legend(loc="center right")
plt.show()
