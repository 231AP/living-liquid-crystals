import os
import subprocess
import numpy as np 
import CreateMatrix
import DrawPolar
import LLC_ABP
import PDF
import Ek
import VisualVideo

if __name__ == "__main__":
    llc = LLC_ABP.LLC_ABP()
    arr_numParticles = [1000,2000,5000,8000]
    arr_alpha = [-1,-3,-5,-10,-20]
    arr_gammaB = [10,1,3]
    arr_xi = [0.9,0.3,0]
    # arr_convection = [1,0]
    llc.convection = 0
    llc.xiAn= 0
    draw_photo = 1
    draw_video = 1
    for xi in arr_xi:
        for gammaB in arr_gammaB:
            for alpha in arr_alpha:
                for numParticles in arr_numParticles:
                    llc.numParticles = numParticles
                    llc.alpha = alpha
                    llc.gammaB = gammaB
                    llc.xi = xi
                    llc.run()    
                    CreateMatrix.CreateMatrix(-1,llc.datadir,llc.numGridX,llc.numGridY)
                    #1 为C型锚定，0为圆形锚定，-1为无锚定
                    os.system(f"cp -r ../RunPolar {llc.datadir}RunPolar")
                    os.system(f"cp -r ../src {llc.datadir}src")
                    subprocess.run(['perl',"run_cuda.pl" ])
                    if draw_photo:
                        for time in np.arange(llc.tStart+1,llc.tStop,llc.tExpo):
                            DrawPolar.DrawPolar(llc.datadir,llc.savedir,llc.savename,time,llc.numGridX,llc.dx)
                        savedir = llc.savedir+"../"+llc.savename
                        PDF.PDF(llc.datadir,savedir,size = llc.numGridX)
                        Ek.Ek(llc.datadir,savedir,size= llc.numGridX)
                    if draw_video:
                        outputFile = savedir + ".mp4"
                        VisualVideo.png_to_mp4(llc.savedir,savedir,fps = 10)
                        

    
