import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import cv2

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

def MSD_defect(datadir,savedir,savename,tStart,tStop,size,steph):
    
    for time in range(tStart,tStop):
        Qxx_df = pd.read_csv(datadir + "abParticle.Qxx_%d.dat"%time,header = None)
        Qxy_df = pd.read_csv(datadir + "abParticle.Qxy_%d.dat"%time,header = None)
        Pxx_df = pd.read_csv(datadir + "abParticle.Pxx_%d.dat"%time,header = None)
        Pxy_df = pd.read_csv(datadir + "abParticle.Pxy_%d.dat"%time,header = None)
      
        Qxx = np.array(Qxx_df).reshape(size,size,order = "F").T
        Qxy = np.array(Qxy_df).reshape(size,size,order = "F").T

        Pxx = np.array(Pxx_df).reshape(size,size,order = "F").T
        Pxy = np.array(Pxy_df).reshape(size,size,order = "F").T
        theta = 0.5 * np.arctan2(Qxy,Qxx)
        theta1 = 0.5 * np.arctan2(Pxy,Pxx)
      
        X,Y = np.meshgrid(np.arange(size),np.arange(size))
        xDP, xDM = detect_defect(theta)
        pDP, pDM = detect_defect(theta1)
        theta_xDP = defect_orie(xDP,0.5,Qxx,Qxy,2,steph)
        theta_pDP = defect_orie(pDP,0.5,Pxx,Pxy,2,steph)
        pyDP1, pxDP1 = np.where(pDP > 0.5)
        pyDM1, pxDM1 = np.where(pDM > 0.5)
        yDP1, xDP1 = np.where(xDP > 0.5)
        yDM1, xDM1 = np.where(xDM > 0.5)
        
        arrDP1 = []
        print(arrDP1)


    
    
if __name__ == "__main__":

    AA = -187
    DC = 10
    V0=10
    vmax1 = 160

    # "N500cMean00005CoAA-3xian1",
    savenames = ["numParticles2000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5-Video20000.0","numParticles2000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5-Video20000.0"]
    time_range = np.arange(1,500)
    steph = 1
    size = 128 * steph
    for savename in savenames:

        # datadir = "/home/lyuan/share/LLC/T_phase/data/heat_Dc10/4NL128/"+savename+"/"
        # savedir = "/home/lyuan/share/LLC/T_phase/photo_video/heat_Dc10/4NL128/"+savename+"/"
        datadir = "../data/"+savename+"/"
        savedir = "../photo_video/analysis/"+savename+"/"
        
        if not os.path.exists(savedir):
            os.makedirs(savedir,exist_ok=True)
        time = 90
        MSD_defect(datadir,savedir,savename,10000,10010,size,steph)
        