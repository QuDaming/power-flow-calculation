#coding=utf-8
import numpy as np

def getLineData(myDataFile):
    return np.array(np.loadtxt(myDataFile,skiprows=1))
def getNodeData(myNodeFile):
    return np.array(np.loadtxt(myNodeFile,skiprows=1))