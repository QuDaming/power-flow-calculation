#coding=utf-8
#得到Y阵
import numpy as np

def get_Y(lineData,nodeData):
    num_node = nodeData.shape[0]                      
    num_line = lineData.shape[0]
    Y = np.zeros((num_node,num_node),dtype=complex)
    #print(Y)
    for i in range(num_line):
        node1 = int(lineData[i,0] - 1)
        node2 = int(lineData[i,1] - 1)
        R = lineData[i,2]
        X = lineData[i,3]
        G = lineData[i,4]
        B = lineData[i,5]
        K = lineData[i,6]
        Z = R+1j*X
        if K == 1:    #普通的线路
            Y[node1][node1] = Y[node1][node1] + 1/Z + (G+1j*B)
            Y[node2][node2] = Y[node2][node2] + 1/Z + (G+1j*B)
            Y[node1][node2] = Y[node1][node2] - 1/Z
            Y[node2][node1] = Y[node2][node1] - 1/Z
        elif K == 0:  #接地支路 
            Y[node1][node2] = Y[node1][node2] + (G+1j*B)
        elif K != 1:  #变压器支路
            Y[node1,node1] = Y[node1,node1] + (K-1)/(K*Z) + 1/(K*Z) 
            Y[node2,node2] = Y[node2,node2] + (1-K)/(K**2*Z) + 1/(K*Z)
            Y[node1,node2] = Y[node1,node2] - 1/(K*Z) 
            Y[node2,node1] = Y[node2,node1] - 1/(K*Z)
    return Y
