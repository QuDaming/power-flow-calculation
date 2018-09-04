# power-flow-calculation
Python N-Lmethod
用牛拉法进行潮流计算
数据格式：
文件1：linedata
i j   R    X   G   B   k
1 2 0.04 0.250 0 0.25 1.00
1 3 0.10 0.350 0 0.00 1.00
2 3 0.08 0.300 0 0.25 1.00
4 2 0.00 0.015 0 0.00 1.05
5 3 0.00 0.030 0 0.00 1.05
文件2：nodedata
i	 type	Pg	Qg	Pl	 Ql 	Um	theta       （pq=1  pv=2）
1	 1	    0	  0	  0	   0	  1	    0
2	 1	    0	  0	  0.5	 0.3	1	    0
3	 2	    0.2	0	  0	   0	  1.05	0
4  0	    0	  0	  0.15 0.1	1.05	0	
