#coding=utf-8
import numpy as np

ERROR = 0.001
def get_numPQ(nodeData):
    num_node = nodeData.shape[0]    #得到节点数
    num_PQ = 0
    for i in range(num_node):
        if nodeData[i,1] == 1:
            num_PQ += 1
    return num_PQ
def get_numPV(nodeData):
    num_node = nodeData.shape[0]    #得到节点数
    num_PV = 0
    for i in range(num_node):
        if nodeData[i,1] == 2:
            num_PV += 1
    return num_PV

def get_difPower(Y,nodeData):
    """
    @func: 得到功率误差\n
    @param：节点导纳阵(Y)，节点参数信息(nodeData)\n
    @return:功率误差 dif_power
    """
    ans_P = equ(Y,nodeData)[0,:]    #计算有功方程
    ans_Q = equ(Y,nodeData)[1,:]    #计算无功方程
    num_node = nodeData.shape[0]    #得到节点数
    dif = np.zeros((2,num_node))    #定义差值的数组
    dif_P = dif[0,:]                #差值P
    dif_Q = dif[1,:]                #差值Q
    # 分解nodedata信息
    num_PQ= 0                       #PQ节点计数
    ntype = nodeData[:,1]
    Pg,Qg,Pl,Ql = nodeData[:,2],nodeData[:,3],nodeData[:,4],nodeData[:,5]
    U,theta = nodeData[:,6],nodeData[:,7]    
    # 计算dif并存储
    for i in range(num_node):
        if ntype[i] == 1:           #pq
            dif_P[i] = Pg[i]-Pl[i]-ans_P[i]
            dif_Q[num_PQ] = Qg[i]-Ql[i]-ans_Q[i]
           # solve[i] = theta[i] 
           # solve[num_node-1+num_PQ] = U[i]
            num_PQ = num_PQ + 1
        elif ntype[i] == 2:         #pv
            dif_P[i] = Pg[i]-Pl[i]-ans_P[i]
         #   solve[i] = theta[i] 
    #修正量
    dif_power = np.hstack([dif_P[0:num_node-1],dif_Q[0:num_PQ]]).reshape((-1,1))  #修正量
    return dif_power

def get_Jacbo(Y,nodeData):
    #得到Jacbo矩阵
    num_node = nodeData.shape[0]    #得到节点数
    num_PQ = get_numPQ(nodeData)
    Pg,Qg,Pl,Ql = nodeData[:,2],nodeData[:,3],nodeData[:,4],nodeData[:,5]
    U,theta = nodeData[:,6],nodeData[:,7]  
    #计算H　N　M　L
    H = np.zeros((num_node-1,num_node-1))   #  H(n-1,n-1)
    N = np.zeros((num_node-1,num_PQ))       #  N(n-1,m)
    M = np.zeros((num_PQ,num_node-1))       #  M(m.n-1)
    L = np.zeros((num_PQ,num_PQ))           #  L(m,m)
    # H 的非对角元
    for i in range(num_node-1):
        for j in range(num_node-1):
            if i != j:
                H[i,j] = -U[i]*U[j]*(Y[i,j].real*np.sin(theta[i]-theta[j]) \
                                    -Y[i,j].imag*np.cos(theta[i]-theta[j]))                   
    # N 的非对角元
    for i in range(num_node-1):
        for j in range(num_PQ):
            if i != j:
                N[i,j] = -U[i]*U[j]*(Y[i,j].real*np.cos(theta[i]-theta[j]) \
                                    +Y[i,j].imag*np.sin(theta[i]-theta[j]))
    # M 的非对角元
    for i in range(num_PQ):
        for j in range(num_node-1):
            if i != j:
                M[i,j] = U[i]*U[j]*(Y[i,j].real*np.cos(theta[i]-theta[j]) \
                                    +Y[i,j].imag*np.sin(theta[i]-theta[j]))
    # L 的非对角元
    for i in range(num_PQ):
        for j in range(num_PQ):
            if i != j:
                L[i,j] = -U[i]*U[j]*(Y[i,j].real*np.sin(theta[i]-theta[j]) \
                                    -Y[i,j].imag*np.cos(theta[i]-theta[j]))
    # H 的对角元
    for i in range(num_node-1):
        for j in range(num_node):
            if i != j:
                H[i,i] = H[i,i] + U[j]*(Y[i,j].real*np.sin(theta[i]-theta[j]) \
                                      -Y[i,j].imag*np.cos(theta[i]-theta[j]))
        H[i,i] = H[i,i]*U[i]   
    # N 的dui]角元
    for i in range(np.min([num_node-1,num_PQ])):
        for j in range(num_node):
            if i != j:
                N[i,i] = N[i,i] + U[j]*(Y[i,j].real*np.cos(theta[i]-theta[j]) \
                                      +Y[i,j].imag*np.sin(theta[i]-theta[j]))
        N[i,i] = -1*(N[i,i]*U[i] + 2*U[i]**2*Y[i,i].real) 
    # M 的对角元
    for i in range(np.min([num_node-1,num_PQ])):
        for j in range(num_node):
            if i != j:
                M[i,i] = M[i,i] + U[j]*(Y[i,j].real*np.cos(theta[i]-theta[j]) \
                             +Y[i,j].imag*np.sin(theta[i]-theta[j]))
        M[i,i] = - M[i,i]*U[i]  
    # L 的对角元
    for i in range(num_PQ):
        for j in range(num_node):
            if i != j:
                L[i,i] = L[i,i] + U[j]*(Y[i,j].real*np.sin(theta[i]-theta[j]) \
                             -Y[i,j].imag*np.cos(theta[i]-theta[j]))
        L[i,i] = - L[i,i]*U[i] + 2*U[i]**2*Y[i,i].imag 

    hn = np.hstack([H,N])
    ml = np.hstack([M,L])
    Jacbo = np.vstack([hn,ml])
    return Jacbo

def get_solve(Jacbo,dif_power,nodeData,solve):
    #求解修正方程
    num_node = nodeData.shape[0]    #得到节点数
    iter_U_theta = -np.dot(np.linalg.inv(Jacbo),dif_power)
    solve[0:num_node-1] =  solve[0:num_node-1] + iter_U_theta[0:num_node-1]
    solve[num_node-1:] = solve[num_node-1:] + solve[num_node-1:]*iter_U_theta[num_node-1:]
    return solve

def correct_nodeData(nodeData,solve):
    num_node = nodeData.shape[0]   
    for i in range(num_node):
        if nodeData[i,1] == 1:
            nodeData[i,7] = solve[i]
            nodeData[i,6] = solve[num_node-1+i]
        elif nodeData[i,1] == 2:
            nodeData[i,7] = solve[i]
    return nodeData

def equ(Y,nodeData):
    """
    计算功率方程的函数\n
    参数：节点导纳阵Y，节点数据\n
    返回值：功率方程的数值：第一维有功；第二维无功
    """
    num_node = nodeData.shape[0]
    U,theta = nodeData[:,6],nodeData[:,7]
    ans = np.zeros((2,4))
    for i in range(num_node):
        for j in range(num_node):
            ans[0,i] = ans[0,i]+U[j]*(Y[i,j].real*np.cos(theta[i]-theta[j])+Y[i,j].imag*np.sin(theta[i]-theta[j]))
            ans[1,i] = ans[1,i]+U[j]*(Y[i,j].real*np.sin(theta[i]-theta[j])-Y[i,j].imag*np.cos(theta[i]-theta[j]))
        ans[0,i] = ans[0,i]*U[i]
        ans[1,i] = ans[1,i]*U[i]
    return ans

def solve_shape(nodeData):
    num_PQ = get_numPQ(nodeData)
    num_PV = get_numPV(nodeData) 
    num_node = nodeData.shape[0]    #得到节点数
    solve = np.zeros(num_PQ*2+num_PV)
    solve[num_node-1:] = np.full((1,num_PQ),1)
    solve = solve.reshape(-1,1)
    return solve