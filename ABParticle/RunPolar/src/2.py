import numpy as np
import matplotlib.pyplot as plt 

# 空间单位 ：L（um） =  sqrt(K0/a) = sqrt(25 * 10-12 m2 ) = 5
# 时间单位 ：T (sec) = 1/Gamma0/a = 2.5
# 应力单位 ：sigma （N /m2） =  0.4
# 力单位 ：K_0 (pN) = 10



K0 = 10 #10-12N
Er0 = 3.75 
Gamma0 = 1 #m sec /kg   
eta0 = 0.5 # kg/m/sec
xian0 = 0.5 #1/sec
h0 = 20 # um

a0 = 0.4#N/m2
b0 = 0.8#N/m2

l0 = 5 #um
V0 = 15 #um/sec
tau = 200 #sec
Dc0 = 200 # um2/sec
AA0 = 187# 10-9 N um


