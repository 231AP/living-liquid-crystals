import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import VisualVideo
import os
def d1xO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[:, n-2:n], u, u[:, 0:2]), axis=1)
    dxdu = (1 / (12 * h)) * (u1[:, 0:n] - 8 * u1[:, 1:n+1] + 8 * u1[:, 3:n+3] - u1[:, 4:n+4])
    return dxdu

def d1yO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[m-2:m, :], u, u[0:2, :]), axis=0)
    dxdu = (1 / (12 * h)) * (u1[0:m, :] - 8 * u1[1:m+1, ] + 8 * u1[3:m+3, :] - u1[4:m+4, :])
    return dxdu
def d1xF(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[:, n-2:n], u, u[:, 0:2]), axis=1)
    dxdu = (1 / h) * (u1[:, 0:n] - u1[:, 1:n+1] )
    return dxdu

def rotate_image(arr, angle):
   height, width = arr.shape[:2]
   # get the image centers
   image_center = (width/2, height/2)

   rotation_arr = cv2.getRotationMatrix2D(image_center, angle, 1)

   abs_cos = abs(rotation_arr[0,0])
   abs_sin = abs(rotation_arr[0,1])

   bound_w = int(height * abs_sin + width * abs_cos)
   bound_h = int(height * abs_cos + width * abs_sin)

   rotation_arr[0, 2] += bound_w/2 - image_center[0]
   rotation_arr[1, 2] += bound_h/2 - image_center[1]

   rotated_mat = cv2.warpAffine(arr, rotation_arr, (bound_w, bound_h))
 
   return rotated_mat

def d1xO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[:, n-2:n], u, u[:, 0:2]), axis=1)
    dxdu = (1 / (12 * h)) * (u1[:, 0:n] - 8 * u1[:, 1:n+1] + 8 * u1[:, 3:n+3] - u1[:, 4:n+4])
    return dxdu

def d1yO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[m-2:m, :], u, u[0:2, :]), axis=0)
    dxdu = (1 / (12 * h)) * (u1[0:m, :] - 8 * u1[1:m+1, ] + 8 * u1[3:m+3, :] - u1[4:m+4, :])
    return dxdu

def d2xO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[:, n-2:n], u, u[:, 0:2]), axis=1)
    d2xdu = (1 / (12 * h**2)) * (-u1[:, 0:n] + 16 * u1[:, 1:n+1] - 30 * u1[:, 2:n+2] + 16 * u1[:, 3:n+3] - u1[:, 4:n+4])
    return d2xdu

def d2yO4(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[m-2:m, :], u, u[0:2, :]), axis=0)
    d2ydu = (1 / (12 * h**2)) * (-u1[0:m, :] + 16 * u1[1:m+1, :] - 30 * u1[2:m+2, :] + 16 * u1[3:m+3, :] - u1[4:m+4, :])
    return d2ydu

def LaplO4(u, h):
    return d2xO4(u,h) + d2yO4(u,h)
def diff_ccw(a, b):
    c = b - a

    # There are 3 cases
    # First, C is greater than pi/2 then we have to invert the angle of b to get the smallest difference
    c[c >= np.pi / 2] = c[c >= np.pi / 2] - np.pi

    # Second case is when C is less than -pi/2. We have to invert the angle of a
    c[c <= -np.pi / 2] = c[c <= -np.pi / 2] + np.pi

    delta = c
    return delta
def detect_defect(theta):
    # Consider the following matrix where each entry represents the angle of the particle.
    # theta = [ A | B ]
    #         [ C | D ]

    # Calculate angles around a point D by using circshift
    # theta_left = [ B | A ]
    #              [ D | C ]
    #              --------> shifted right so that element (i, j) of theta_left
    #                        represents a grid to the left of the original theta
    theta = np.mod(theta, np.pi)
    theta_left = np.roll(theta, shift=(0, 1), axis=(0, 1))
    theta_up = np.roll(theta, shift=(1, 0), axis=(0, 1))
    theta_left_up = np.roll(theta, shift=(1, 1), axis=(0, 1))

    # Calculate the difference between theta of a loop in space
    up_diff = diff_ccw(theta, theta_up)
    up_left_diff = diff_ccw(theta_up, theta_left_up)
    left_diff = diff_ccw(theta_left_up, theta_left)
    diff = diff_ccw(theta_left, theta)

    # Check if angle is pi/2 or -pi/2 or 0.
    tol = 0.3
    strength = up_diff + up_left_diff + left_diff + diff
    loc12 = np.abs(strength + np.pi) / np.pi < tol
    loc_12 = np.abs(strength - np.pi) / np.pi < tol

    return loc12, loc_12
def CreatePeriodicMatrix(matrix,delta_r):
    lenth = matrix.shape[0]
    width = matrix.shape[1]
    matrix_bigger = np.zeros((lenth+2*delta_r,width+2*delta_r))
   
    matrix_bigger[delta_r:lenth+delta_r , delta_r:width+delta_r] = matrix



    matrix_bigger[0:delta_r,0:delta_r] = matrix[lenth-delta_r:lenth,width-delta_r:]
    matrix_bigger[0:delta_r,delta_r:lenth+delta_r] = matrix[lenth-delta_r:lenth,:]
    matrix_bigger[0:delta_r,lenth+delta_r:] = matrix[lenth-delta_r:lenth,0:delta_r]


    matrix_bigger[lenth+delta_r:,0:delta_r] = matrix[0 :delta_r,width-delta_r:]
    matrix_bigger[lenth+delta_r:,delta_r:lenth+delta_r] = matrix[0:delta_r,:]
    matrix_bigger[lenth+delta_r:,lenth+delta_r:] = matrix[0:delta_r,0:delta_r]


    matrix_bigger[delta_r:lenth+delta_r,0:delta_r] = matrix[:,width-delta_r:]
    matrix_bigger[delta_r:lenth+delta_r,width+delta_r:] = matrix[:,0:delta_r]

    return matrix_bigger

def defect_orie(xD, k, Qxx, Qxy, Loop,h):
    dn = Loop
    m, n = Qxx.shape

    dxQxx = d1xO4(Qxx, h)
    dxQxy = d1xO4(Qxy, h)
    dyQxx = d1yO4(Qxx, h)
    dyQxy = d1yO4(Qxy, h)

    signk = np.sign(k)

    top = signk * dxQxy - dyQxx
    bottom = dxQxx + signk * dyQxy

    top1 = np.block([[top[m - dn:m, n - dn:n], top[m - dn:m, :], top[m - dn:m, 0:dn]],
                        [top[:, n - dn:n], top, top[:, 0:dn]],
                        [top[0:dn, n - dn:n], top[0:dn, :], top[0:dn, 0:dn]]])

    bottom1 = np.block([[bottom[m - dn:m, n - dn:n], bottom[m - dn:m, :],
                            bottom[m - dn:m, 0:dn]],
                        [bottom[:, n - dn:n], bottom, bottom[:, 0:dn]],
                        [bottom[0:dn, n - dn:n], bottom[0:dn, :], bottom[0:dn, 0:dn]]])

    topa = np.zeros((m, n))
    bottoma = np.zeros((m, n))

    for j in range(-Loop, Loop):
        for i in range(-Loop, Loop):
            topa += top1[dn + i:dn + m + i, dn + j:dn + n + j]
            bottoma += bottom1[dn + i:dn + m + i, dn + j:dn + n + j]

    pD = k / (1 - k) * xD * np.arctan2(topa, bottoma)

    return pD

def ConvertToPassive(angle):
    return (angle%(np.pi)+np.pi)%(np.pi)

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def DrawPhase(datadirs, savedir, lenth_fig, hight_fig,tStart=1, tStop=90, size=128, V0=7.5):
    
    font = FontProperties()
    font.set_size("xx-large")
    id_data = 0
    os.makedirs(savedir,exist_ok=True)
    h = 1
    for time in np.arange(tStart, tStop):
        fig, axes = plt.subplots(lenth_fig, hight_fig, figsize=(4*hight_fig, 4*lenth_fig))
        id_data = 0
        for datadir in datadirs:
            i = id_data % lenth_fig
            j = int(id_data/lenth_fig)
        # Loop through the time steps and process data
        
            omega_df = pd.read_csv(datadir + "incompFlow.omega_%d.dat"%time,header=None)
            vx_df = pd.read_csv(datadir + "incompFlow.vx_%d.dat"%time,header=None)
            vy_df = pd.read_csv(datadir + "incompFlow.vy_%d.dat"%time,header=None)
            Qxx_df = pd.read_csv(datadir + "abParticle.Qxx_%d.dat"%time,header = None)
            Qxy_df = pd.read_csv(datadir + "abParticle.Qxy_%d.dat"%time,header = None)
            Pxx_df = pd.read_csv(datadir + "abParticle.Pxx_%d.dat"%time,header = None)
            Pxy_df = pd.read_csv(datadir + "abParticle.Pxy_%d.dat"%time,header = None)
            anchx_df = pd.read_csv(datadir + "anchx_%d.dat"%time,header = None)
            anchy_df = pd.read_csv(datadir + "anchy_%d.dat"%time,header = None)
            C_df = pd.read_csv(datadir + "abParticle.Concentration_%d.dat"%time,header=None)
            kk = pd.read_csv(datadir+"conf_%d"%time+".dat",sep=" ",header=None)
            pts = np.array(kk)
            Pxx = np.array(Pxx_df).reshape(size,-1,order = "F").T
            Pxy = np.array(Pxy_df).reshape(size,-1,order = "F").T
            anchx = np.array(anchx_df).reshape(size,-1,order = "F").T
            anchy = np.array(anchy_df).reshape(size,-1,order = "F").T
            omega = np.array(omega_df).reshape(size,size,order = "F").T
            vx = np.array(vx_df).reshape(size,size,order = "F").T
            vy = np.array(vy_df).reshape(size,size,order = "F").T
            Qxx = np.array(Qxx_df).reshape(size,size,order = "F").T
            Qxy = np.array(Qxy_df).reshape(size,size,order = "F").T
            C = np.array(C_df).reshape(size,size,order = "F").T
            theta = 0.5 * np.arctan2(Qxy,Qxx)
            cbsize = 16
            titlesize = 20
            cut_size = 4
            cut_size1 = 1
            S2 = 4*Qxx**2+4*Qxy**2
            A2 = 4*anchx**2+4*anchy**2
            P2 = 4*Pxx**2+4*Pxy**2
            v_abs = np.sqrt(vx **2 + vy **2)
            n1 = np.cos(theta)
            n2 = np.sin(theta)
            n1cut = n1[::cut_size,::cut_size]
            n2cut = n2[::cut_size,::cut_size]
            theta1 = 0.5 * np.arctan2(Pxy,Pxx)
            p11 = np.cos(theta1)
            p22 = np.sin(theta1)
            p11cut = p11[::cut_size,::cut_size]
            p22cut = p22[::cut_size,::cut_size]
            axes[i, j].imshow(omega,vmin = -3,vmax = 3,origin='lower', animated=True, cmap='jet')
            axes[i, j].quiver(pts[:,0],pts[:,1],100*pts[:,2],100*pts[:,3],0.01)
            axes[i, j].set_xlim(0,size)
            axes[i, j].set_ylim(0,size)
            axes[i,j].axis("off")
            print(i,j)
            print(datadir)
            id_data+=1
            # if id_data == len(datadirs)-1:
            # print(id_data)
         
            if id_data == len(datadirs):
                print(time)
                plt.savefig(savedir + "DrawPhase%04d.png"%time)
            
        plt.close()
        
def DrawPhaseS(datadirs, savedir, lenth_fig, hight_fig,tStart=1, tStop=90, size=128, V0=7.5):

    font = FontProperties()
    font.set_size("xx-large")
    id_data = 0
    os.makedirs(savedir,exist_ok=True)
    h = 1
    steph = h
    for time in np.arange(tStart, tStop):
        fig, axes = plt.subplots(lenth_fig, hight_fig, figsize=(4*hight_fig, 4*lenth_fig))
        id_data = 0
        for datadir in datadirs:
            i = id_data % lenth_fig
            j = int(id_data/lenth_fig)
        # Loop through the time steps and process data
        
            omega_df = pd.read_csv(datadir + "incompFlow.omega_%d.dat"%time,header=None)
            vx_df = pd.read_csv(datadir + "incompFlow.vx_%d.dat"%time,header=None)
            vy_df = pd.read_csv(datadir + "incompFlow.vy_%d.dat"%time,header=None)
            Qxx_df = pd.read_csv(datadir + "abParticle.Qxx_%d.dat"%time,header = None)
            Qxy_df = pd.read_csv(datadir + "abParticle.Qxy_%d.dat"%time,header = None)
            Pxx_df = pd.read_csv(datadir + "abParticle.Pxx_%d.dat"%time,header = None)
            Pxy_df = pd.read_csv(datadir + "abParticle.Pxy_%d.dat"%time,header = None)
            anchx_df = pd.read_csv(datadir + "anchx_%d.dat"%time,header = None)
            anchy_df = pd.read_csv(datadir + "anchy_%d.dat"%time,header = None)
            C_df = pd.read_csv(datadir + "abParticle.Concentration_%d.dat"%time,header=None)
            kk = pd.read_csv(datadir+"conf_%d"%time+".dat",sep=" ",header=None)
            pts = np.array(kk)
            Pxx = np.array(Pxx_df).reshape(size,-1,order = "F").T
            Pxy = np.array(Pxy_df).reshape(size,-1,order = "F").T
            anchx = np.array(anchx_df).reshape(size,-1,order = "F").T
            anchy = np.array(anchy_df).reshape(size,-1,order = "F").T
            omega = np.array(omega_df).reshape(size,size,order = "F").T
            vx = np.array(vx_df).reshape(size,size,order = "F").T
            vy = np.array(vy_df).reshape(size,size,order = "F").T
            Qxx = np.array(Qxx_df).reshape(size,size,order = "F").T
            Qxy = np.array(Qxy_df).reshape(size,size,order = "F").T
            C = np.array(C_df).reshape(size,size,order = "F").T
            theta = 0.5 * np.arctan2(Qxy,Qxx)
            cbsize = 16
            titlesize = 20
            cut_size = 4
            cut_size1 = 1
            S2 = 4*Qxx**2+4*Qxy**2
            A2 = 4*anchx**2+4*anchy**2
            P2 = 4*Pxx**2+4*Pxy**2
            v_abs = np.sqrt(vx **2 + vy **2)
            n1 = np.cos(theta)
            n2 = np.sin(theta)
            n1cut = n1[::cut_size,::cut_size]
            n2cut = n2[::cut_size,::cut_size]
            theta1 = 0.5 * np.arctan2(Pxy,Pxx)
            p11 = np.cos(theta1)
            p22 = np.sin(theta1)
            p11cut = p11[::cut_size,::cut_size]
            p22cut = p22[::cut_size,::cut_size]
            
            
            X,Y = np.meshgrid(np.arange(size),np.arange(size))
            xDP, xDM = detect_defect(theta)
            pDP, pDM = detect_defect(theta1)
            theta_xDP = defect_orie(xDP,0.5,Qxx,Qxy,2,steph)
            theta_pDP = defect_orie(pDP,0.5,Pxx,Pxy,2,steph)
            pyDP1, pxDP1 = np.where(pDP > 0.5)
            pyDM1, pxDM1 = np.where(pDM > 0.5)
            yDP1, xDP1 = np.where(xDP > 0.5)
            yDM1, xDM1 = np.where(xDM > 0.5)
        
            axes[i, j].imshow(S2,vmax =1,vmin= 0,norm=None,origin="lower")
       
            axes[i,j].quiver(X[::cut_size,::cut_size],Y[::cut_size,::cut_size],n1cut,n2cut,color = "red" ,angles='xy', scale_units='xy', scale=1/6, headwidth=0, headlength=0 ,headaxislength=0, pivot='middle')
            axes[i, j].quiver(pts[:,0],pts[:,1],100*pts[:,2],100*pts[:,3],0.01)
            axes[i, j].set_xlim(0,size)
            axes[i, j].set_ylim(0,size)
            axes[i,j].axis("off")
            # print(i,j)
            # print(datadir)
            id_data+=1
            # if id_data == len(datadirs)-1:
            # print(id_data)
         
            if id_data == len(datadirs):
                print(time)
                plt.savefig(savedir + "DrawPhase%04d.png"%time)
            
        plt.close()
if __name__=="__main__":
    datadir = "../data/numParticles400-gammaB3-alpha-8-xi0.9-convection0-Dp0.01-kBT0.05-MatrixC/"
    datadir1 = "../data/numParticles4000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir2 = "../data/numParticles6000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir3 = "../data/numParticles400-gammaB3-alpha-10-xi0.9-convection0-Dp0.01-kBT0.05-MatrixC/"
    arr_numParticles = [2000,4000,6000, 8000]
    arr_alpha = [-12,-15,-20]
    arr_k = [0.5,1,2,4]
    arr_gammaB = [0.5,1,3,10]
    datadirs = []

    
    for alpha in arr_alpha:
        for numParticles in arr_numParticles:
            datadirs.append(
                f"../data/numParticles{numParticles}-gammaB3-alpha{alpha}-xi0.9-convection1-Dp0-kBT0.5/"
            )

    # Print the results for verification
    # for path in datadirs:
    # #     print(path)


    
    # datadirs = []

    
    # for k in arr_k:
    #     for gammaB in arr_gammaB:
    #         datadirs.append(
    #             f"../data/numParticles400-gammaB{gammaB}-alpha-8-xi0.9-convection1-Dp0.01-kBT0.05-K{k}-MatrixC/"
    #         )
    # for k in k_arr:
    #     for 

    # datadirs[i] = ""
    # datadirs = [datadir,datadir1,datadir2,datadir3]
    # datadirs = [datadir,datadir3]
    savedir = "../photo_video/analysis/"+ "example5"+"alpha_num/"
    # DrawPhaseS(datadirs,savedir,lenth_fig=4, hight_fig=3,tStart=1, tStop=300, size=128, V0=7.5)
    VisualVideo.png_to_mp4(savedir, "../photo_video/analysis/"+"example5"+"alpha_num", 10)


