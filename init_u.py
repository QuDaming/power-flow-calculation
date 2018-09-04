#coding=utf-8
import numpy as np

def init_U(nodeData):       
    num_node = nodeData.shape[0]
    for i in range(num_node):
        if nodeData[i,1] == 1:   #pq节点
            nodeData[i,6],nodeData[i,7] = 1,0
        elif nodeData[i,1] == 2: #pv节点
            nodeData[i,7] = 0