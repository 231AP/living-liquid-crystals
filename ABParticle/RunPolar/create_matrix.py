
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

savenames = ["test06",]

for savename in savenames:
    datadir = "../data/"+savename+"/"
    savedir = "../photo_video/"+savename+"/"
    if not os.path.exists(savedir):
        os.makedirs(savedir,exist_ok=True)
    if not os.path.exists(datadir):
        os.makedirs(datadir,exist_ok=True)


# create matrixC 
def CreateZeros(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T
def CreateMatrixC(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    for i in range(lenth_y):
        matrix1x[i,:] = np.cos(i/L*np.pi+np.pi)
        matrix1y[i,:] = np.sin(i/L*np.pi+np.pi)
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

def CreateMatrixF(lenth_x,lenth_y,L):
    matrix1x = -np.zeros((lenth_x,lenth_y))
    matrix1y = -np.zeros((lenth_x,lenth_y))
    for i in range(lenth_y):
        matrix1x[i,:] = np.cos(i/L*np.pi+np.pi)
        matrix1y[i,:] = -np.sin(i/2/L*np.pi+np.pi/2)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T

def CreateMatrixS1(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    # for theta in np.arange(0,2*np.pi,0.1):
    for i in range(lenth_y):
        matrix1x[i,:] = np.cos(i/L*np.pi+np.pi)
        matrix1y[i,:] = -np.sin(i/2/L*np.pi+np.pi/2)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T


def CreateMatrixS0(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    # for theta in np.arange(0,2*np.pi,0.1):
    for i in range(lenth_y):
            matrix1x[i,:] = np.cos(i/L*np.pi+np.pi)
            matrix1y[i,:] = np.sin(i/L*np.pi+(-1)**(np.sin(i*np.pi)+2**1.5)*np.pi)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T
def CreateMatrixS1(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    # for theta in np.arange(0,2*np.pi,0.1):
    for i in range(lenth_y):
        for j in range(lenth_x):
            theta = np.arctan2(i-lenth_x/2,j-lenth_x/2)*2+np.pi
            R = np.sqrt((i-lenth_x/2)**2+(j-lenth_x/2)**2)
            matrix1x[i,j] = np.cos(R/L*np.pi+theta+np.pi)
            matrix1y[i,j] = np.sin(R/L*np.pi+theta+np.pi)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T
def CreateMatrixS1(lenth_x,lenth_y,L):
    matrix1x = np.zeros((lenth_x,lenth_y))
    matrix1y = np.zeros((lenth_x,lenth_y))
    for i in range(lenth_y):
        for j in range(lenth_x):
            theta = np.arctan2(i-lenth_x/2,j-lenth_x/2)*2+np.pi
            R = np.sqrt((i-lenth_x/2)**2+(j-lenth_x/2)**2)
            if R<=lenth_x/2:
                matrix1x[i,j] = np.cos(R/L*np.pi+theta+np.pi)
                matrix1y[i,j] = np.sin(R/L*np.pi+theta+np.pi)
            else:
                matrix1x[i,j] = np.cos(theta-np.pi)
                matrix1y[i,j] = np.sin(theta-np.pi)
    return matrix1x.reshape(1,-1).T,matrix1y.reshape(1,-1).T


# lenth_x = 128
# lenth_y = 128
# L = 32
steph = 1
size = 128 * steph
lenth_x = 128
lenth_y = 128

L = 32
anxx,anxy = CreateMatrixF(lenth_x,lenth_y,L)
# anxx,anxy = CreateZeros(lenth_x,lenth_y,L)
# anxx,anxy = CreateMatrixS1(lenth_x,lenth_y,L)

# np.savetxt(savedir+'anchx_0.dat', anxx, fmt='%f')
# np.savetxt(savedir+'anchy_0.dat', anxy, fmt='%f')
np.savetxt(datadir+'anchx_0.dat', anxx, fmt='%f')
np.savetxt(datadir+'anchy_0.dat', anxy, fmt='%f')
np.savetxt(datadir+'abParticle.Qxx_0.dat', anxx, fmt='%f')
np.savetxt(datadir+'abParticle.Qxy_0.dat', anxy, fmt='%f')




fig = plt.figure(figsize = [20,20])
ax = fig.add_subplot(1,1,1)
anchx_df = pd.read_csv(datadir + "anchx_0.dat",header = None)
anchy_df = pd.read_csv(datadir + "anchy_0.dat",header = None)
anchx = np.array(anchx_df).reshape(size,-1,order = "F").T
anchy = np.array(anchy_df).reshape(size,-1,order = "F").T
A2 = 4*anchx**2+4*anchy**2
X,Y = np.meshgrid(np.arange(size),np.arange(size))
cut_size = 4
theta2 = 0.5 * np.arctan2(anchy,anchx)
a1 = np.cos(theta2)
a2 = np.sin(theta2)
a1cut = a1[::cut_size,::cut_size]
a2cut = a2[::cut_size,::cut_size]
im = ax.imshow(A2,vmax =0.9,vmin= 0,norm=None,origin="lower")
ax.quiver(X[::cut_size,::cut_size],Y[::cut_size,::cut_size],a1cut,a2cut,color = "red" ,angles='xy', scale_units='xy', scale=1/6, headwidth=0, headlength=0 ,headaxislength=0, pivot='middle')
plt.savefig(savedir+'anchor.png')
