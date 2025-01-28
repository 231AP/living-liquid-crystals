import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Ek(datadirs, savedir, tStart=50, time_range=40, size=128, num_bins=50):
    id_data = 0
    for datadir in datadirs:
        
        L = size
        dx = 1
        dy = 1
        u = 0
        v = 0
        
        # 计算频率空间的坐标
        kx = np.fft.fftfreq(L, dx) * 2 * np.pi  # 对应的 kx 频率
        ky = np.fft.fftfreq(L, dy) * 2 * np.pi  # 对应的 ky 频率
        KX, KY = np.meshgrid(kx, ky)

        # 计算频率模 k = sqrt(kx^2 + ky^2)
        k = np.sqrt(KX**2 + KY**2)
        k_max = np.max(k)

        # Logarithmic binning for the radial energy spectrum
        k_bins = np.logspace(np.log10(0.05), np.log10(k_max), num_bins)
        E_radial = np.zeros(num_bins - 1)

        # Loop through the time range
        for time in np.arange(tStart, tStart + time_range):

            vx_df = pd.read_csv(datadir + f"incompFlow.vx_{time}.dat", header=None)
            vy_df = pd.read_csv(datadir + f"incompFlow.vy_{time}.dat", header=None)

            vx = np.array(vx_df).reshape(size, -1, order="F").T
            vy = np.array(vy_df).reshape(size, -1, order="F").T
            
            # Compute mean velocity components
            u_0 = np.mean(vx)
            v_0 = np.mean(vy)
            u = vx
            v = vy

            # Fourier transform of velocity components
            u_hat = np.fft.fft2(u)
            v_hat = np.fft.fft2(v)

            # Compute energy spectrum (squared amplitude)
            E = 0.5 * np.abs(u_hat)**2 + np.abs(v_hat)**2

            # Plot 2D energy spectrum (log scale)
            # plt.figure(figsize=(8, 6))
            # plt.imshow(np.log10(E), extent=[kx.min(), kx.max(), ky.min(), ky.max()], origin='lower', aspect='auto')
            # plt.colorbar(label='Log10 of Energy Spectrum')
            # plt.xlabel('kx')
            # plt.ylabel('ky')
            # plt.title('2D Energy Spectrum (Log Scale)')
            # plt.savefig(savedir + "Ek.png")
            # plt.close()

            # Radial averaging of energy spectrum
            for i in range(1, num_bins):
                k_min = k_bins[i - 1]
                k_max = k_bins[i]
                mask = (k >= k_min) & (k < k_max)
                E_radial[i - 1] += np.mean(E[mask])

            # Prepare for reference lines (y = x^1 and y = x^-6)
            X1 = k_bins[:-30]
            Y1 = X1**1
            X2 = k_bins[-20:]
            Y2 = X2**-4
            if time == tStart+time_range-1:
            # Plot radial average of energy spectrum
                if id_data ==0 :
                    plt.figure(figsize=(8, 6))
                    id_data+=1
                
                plt.loglog(k_bins[:-10], E_radial[:-9], marker='o',label = f"Ek{id_data*2000}")
                plt.loglog(X1, 10**10*Y1, label="y = x^1")
                plt.loglog(X2, 10**8*Y2, label="y = x^-4")
                # plt.xlim(0,1)
                plt.legend()
                plt.xlabel('k')
                plt.ylabel('Energy Spectrum E(k)')
                plt.title('Radially Averaged 2D Energy Spectrum')
                plt.grid(True)
                plt.savefig(savedir + "Ek-k.png")
                id_data+=1
            # plt.close()

if __name__=="__main__":
    
    
    datadir = "../data/numParticles2000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir1 = "../data/numParticles4000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir2 = "../data/numParticles6000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir3 = "../data/numParticles8000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadirs = [datadir,datadir1,datadir2,datadir3]
    savedir = "../photo_video/analysis/"+ "example1"+"gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5"
    Ek(datadirs,savedir)
