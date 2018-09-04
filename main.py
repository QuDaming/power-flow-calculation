#coding=utf-8
import numpy as np
from init_u import init_U
from readfile import *
from get_Y import get_Y
from iteration import *
# 1、传递文件进程序
myLineFile = "F:\\myPython\\power flow caculate\\linedata3.txt"
myNodeFile = "F:\\myPython\\power flow caculate\\nodedata3.txt"
# 2、形成节点导纳矩阵
myNodeData = getNodeData(myNodeFile)
myLineData = getLineData(myLineFile)
Y = get_Y(myLineData,myNodeData)
solve = solve_shape(myNodeData)
# 3、初始化电压幅值和相位
init_U(myNodeData)
# 4、得到修正矩阵（u，theta）  计算功率误差初值
dif_power = get_difPower(Y,myNodeData)
# 5、计算功率误差判断是否循环  (deta_p,deta_q)
num_iter = 0
while np.max(np.fabs(dif_power)) > ERROR:
# 6、计算雅可比矩阵元素 (jacbo)
    Jacbo = get_Jacbo(Y,myNodeData)
# 7、求解修正方程式   (deta_u,deta_theta)
    solve = get_solve(Jacbo,dif_power,myNodeData,solve)
# 8、修正nodeData   (u,theta)
    mynodeData = correct_nodeData(myNodeData,solve)
    dif_power = get_difPower(Y,myNodeData)
    num_iter = num_iter + 1
print(mynodeData[:,6:])
print(num_iter)
