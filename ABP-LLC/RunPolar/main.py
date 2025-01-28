import os
import subprocess
import numpy as np 
import CreateMatrix
import DrawPolar
import LLC_ABP
import PDFall
import Ek
import VisualVideo

if __name__ == "__main__":
    llc = LLC_ABP.LLC_ABP()
    arr_numParticles = [500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000]
    arr_alpha = [-23,-20,-18,-16]
    arr_gammaB = [3]
    arr_xi = [0.9]
    arr_convection = [1]
    arr_kBT = [0.5]
    llc.tExpo = 0.1
    llc.convection = 1
    llc.xiAn= 0
    llc.Dp =0.01
    photoNum = llc.tStop/llc.tExpo
    draw_photo = 0
    draw_video = 0
    for xi in arr_xi:
        for gammaB in arr_gammaB:
            for alpha in arr_alpha:
                for numParticles in arr_numParticles:
                    for convection in arr_convection:
                        for kBT in arr_kBT:
                            llc.kBT = kBT
                            llc.numParticles = numParticles
                            llc.alpha = alpha
                            llc.gammaB = gammaB
                            llc.xi = xi
                            llc.convection = convection
                            llc.savename = f"numParticles{llc.numParticles}-gammaB{llc.gammaB}-alpha{llc.alpha}-xi{llc.xi}-convection{llc.convection}-Dp{llc.Dp}-kBT{llc.kBT}"
                            llc.datadir = f"../data/{llc.savename}/"
                            llc.savedir = f"../photo_video/{llc.savename}/"
                            llc.run()
                            # print(llc.xiAn)
                            
                            CreateMatrix.CreateMatrix(-1,llc.datadir,llc.numGridX,llc.numGridY)
                            #1 为C型锚定，0为圆形锚定，-1为无锚定
                            
                            # subprocess.run(['bash','bash.sh'])
                            os.system(f"cp -r ../RunPolar {llc.datadir}RunPolar")
                            os.system(f"cp -r ../src {llc.datadir}src")
                            subprocess.run(['perl',"run_cuda.pl" ])
                            if draw_photo:
                                for time in np.arange(llc.tStart+1,photoNum):
                                    DrawPolar.DrawPolar(llc.datadir,llc.savedir,llc.savename,time,llc.numGridX,llc.dx)
                                savedir = llc.savedir+"../"+llc.savename
                                PDFall.PDF(llc.datadir,savedir,size = llc.numGridX)
                                Ek.Ek(llc.datadir,savedir,size= llc.numGridX)
                            if draw_video:
                                outputFile = savedir + ".mp4"
                                VisualVideo.png_to_mp4(llc.savedir,savedir,fps = 10)
                                

            

