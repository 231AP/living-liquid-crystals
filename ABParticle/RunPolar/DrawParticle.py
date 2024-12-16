import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


savename = "test11"
datadir = "../data/" + savename +"/"
savedir = "../photo_video/" + savename + "/"
import os 
os.makedirs(savedir,exist_ok=True)
time_range = np.arange(1,100)
for time in time_range:
    kk = pd.read_csv(datadir+"conf_%d"%time+".dat",sep=" ",header=None)
    # print(kk.shape)
    pts = np.array(kk)
    # print(pts[1,:])
    fig = plt.figure(figsize=(20,20))
    
    ax = fig.add_subplot(1,1,1)
    # print(pts.shape)
    ax.quiver(pts[:,0],pts[:,1],100*pts[:,2],100*pts[:,3],0.01)
    
    ax.set_xlim(0,1*128)
    ax.set_ylim(0,1*128)
    plt.savefig(savedir+"%04d"%time+".png")
    plt.close()
    
    
    
    