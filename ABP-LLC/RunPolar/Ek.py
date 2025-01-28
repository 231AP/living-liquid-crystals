import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Ek(datadir, savedir, start_time=50, time_range=40, size=128, num_bins=50):
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
    for time in np.arange(start_time, start_time + time_range):

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
        # u_hat = np.fft.fftshift(np.fft.fft2(u))
        u_hat = np.fft.fft2(u)

        v_hat = np.fft.fft2(v)
        

        # Compute energy spectrum (squared amplitude)
        E = 0.5 * np.abs(u_hat)**2 + np.abs(v_hat)**2

        # Plot 2D energy spectrum (log scale)
        plt.figure(figsize=(8, 6))
        plt.imshow(np.log10(E), extent=[kx.min(), kx.max(), ky.min(), ky.max()], origin='lower', aspect='auto')
        plt.colorbar(label='Log10 of Energy Spectrum')
        plt.xlabel('kx')
        plt.ylabel('ky')
        plt.title('2D Energy Spectrum (Log Scale)')
        plt.savefig(savedir + "Ek.png")
        plt.close()

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

        # Plot radial average of energy spectrum
        plt.figure(figsize=(8, 6))
        plt.loglog(k_bins[:-1], E_radial, marker='o')
        plt.loglog(X1, Y1, label="y = x^1")
        plt.loglog(X2, Y2, label="y = x^-4")
        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Energy Spectrum E(k)')
        plt.title('Radially Averaged 2D Energy Spectrum')
        plt.grid(True)
        plt.savefig(savedir + "Ek-k.png")
        plt.close()

