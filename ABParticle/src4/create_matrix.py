import numpy as np
savedir = "../data/test13/"
import os
if os.path.exists(savedir)==0:
    os.makedirs(savedir)
# create matrixC 
def CreateMatrix1(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    for i in range(lenth_y):
        matrix1x[i,:] = np.cos(i/L*np.pi)
        matrix1y[i,:] = -np.sin(i/L*np.pi)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T



def CreateMatrix2(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    for i in range(lenth_y):
        matrix1x[:,i] = np.cos(i/L*np.pi)
        matrix1y[:,i] = np.sin(i/L*np.pi)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T

def CreateMatrixXdirction(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))+0.5
    matrix1y = np.zeros((lenth_x,lenth_y))

    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T
def CreateZeros(lenth_x,lenth_y):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))

    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T
# lenth_x = 256
# lenth_y = 256
# L = 64
lenth_x = 128
lenth_y = 128
L = 32
# anxx,anxy = CreateMatrix1(lenth_x,lenth_y,L)
anxx,anxy = CreateZeros(lenth_x,lenth_y)

np.savetxt(savedir+'anchx_0.dat', anxx, fmt='%f')
np.savetxt(savedir+'anchy_0.dat', anxy, fmt='%f')