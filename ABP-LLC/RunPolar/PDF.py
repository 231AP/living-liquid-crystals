import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def d1xF(u, h):
    m, n = u.shape
    u1 = np.concatenate((u[:, n-2:n], u, u[:, 0:2]), axis=1)
    dxdu = (1 / h) * (u1[:, 0:n] - u1[:, 1:n+1] )
    return dxdu

def PDF(datadir, savedir,tStart=50,tStop=90,max_velocity=80, max_vorticity=80, max_C=2, size=128):
    unit_velocity = max_velocity / 100
    unit_vorticity = max_vorticity / 100
    unit_C = max_C / 100
    
    # Initialize PDF arrays
    PDF_velocity = np.zeros(100)
    PDF_vx = np.zeros(100)
    PDF_vy = np.zeros(100)
    PDF_vorticity = np.zeros(100)
    PDF_C = np.zeros(100)
    
    # Loop through the time steps and process data
    for time in np.arange(tStart, tStop):
        omega_df = pd.read_csv(datadir + f"incompFlow.omega_{time}.dat", header=None)
        vx_df = pd.read_csv(datadir + f"incompFlow.vx_{time}.dat", header=None)
        vy_df = pd.read_csv(datadir + f"incompFlow.vy_{time}.dat", header=None)
        C_df = pd.read_csv(datadir + f"abParticle.Concentration_{time}.dat", header=None)
        
        omega = np.array(omega_df).reshape(size, size)
        vx = np.array(vx_df).reshape(size, size)
        vy = np.array(vy_df).reshape(size, size)
        C = np.array(C_df).reshape(size, size)
        
        v = np.sqrt(vx**2 + vy**2)
        ALLv4 = 0
        ALLv2 = 0
        ALLdv4 = 0
        ALLdv2 = 0
        dvx = d1xF(vx,1)
        ALLv4 += np.mean(vx**4)
        ALLv2 += np.mean(vx**2)**2
        kurtosis = ALLv4/(ALLv2)
        print("kurtosis_vx = %f"%kurtosis)
        
        ALLdv4 += np.mean(dvx**4)
        ALLdv2 += np.mean(dvx**2)**2
        kurtosis_dvx = ALLdv4/(ALLdv2)
        print("kurtosis_dvx = %f"%kurtosis_dvx)
        for j in range(100):
            i = j - 50
            PDF_velocity[j] += np.sum((v < i * unit_velocity) * (v > (i - 1) * unit_velocity))
            PDF_vx[j] += np.sum((vx < i * unit_velocity) * (vx > (i - 1) * unit_velocity))
            PDF_vy[j] += np.sum((vy < i * unit_velocity) * (vy > (i - 1) * unit_velocity))
            PDF_vorticity[j] += np.sum((omega < i * unit_vorticity) * (omega > (i - 1) * unit_vorticity))
            PDF_C[j] += np.sum((C < i * unit_C) * (C > (i - 1) * unit_C))
        
        normal_velocity = np.sum(PDF_velocity)
        normal_vx = np.sum(PDF_vx)
        normal_vy = np.sum(PDF_vy)
        normal_vorticity = np.sum(PDF_vorticity)
        normal_C = np.sum(PDF_C)
        
        # Normalize the PDFs
        PDF_velocity /= normal_velocity
        PDF_vx /= normal_vx
        PDF_vy /= normal_vy
        PDF_vorticity /= normal_vorticity
        PDF_C /= normal_C

        # Set font properties
        font = FontProperties()
        font.set_size("xx-large")

        # Plot Velocity PDF
        fig = plt.figure(figsize=(10, 9))
        plt.plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vx, label="u")
        plt.plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vy, label="v")
        plt.legend(loc=1, prop=font)
        plt.yscale("log")
        plt.xlabel("Velocity", fontsize=24)
        plt.title(f"kurtosis = {kurtosis},kurtosis_dvx = {kurtosis_dvx}")
        plt.ylabel("PDF", fontsize=24)
        plt.rc("xtick", labelsize=24)
        plt.rc("ytick", labelsize=24)
        plt.ylim((10**-5, 10**-0.1))
        plt.savefig(savedir + "Velocity.png")
        plt.close()

        # Plot Vorticity PDF
        fig = plt.figure(figsize=(10, 9))
        plt.plot(np.arange(-50, 50) / 100 * max_vorticity, PDF_vorticity, label="vorticity")
        plt.legend(loc=1, prop=font)
        plt.yscale("log")
        plt.xlabel("Vorticity", fontsize=24)
        plt.ylabel("PDF", fontsize=24)
        plt.rc("xtick", labelsize=24)
        plt.rc("ytick", labelsize=24)
        plt.ylim((10**-5, 10**-0.1))
        plt.savefig(savedir + "Vorticity.png")
        plt.close()

        # Plot Concentration PDF
        fig = plt.figure(figsize=(10, 9))
        plt.plot(np.arange(-50, 50) / 100 * max_C, PDF_C, label="Concentration")
        plt.legend(loc=1, prop=font)
        plt.xlim((0, 0.5))
        plt.yscale("log")
        plt.xlabel("Concentration", fontsize=24)
        plt.ylabel("PDF", fontsize=24)
        plt.rc("xtick", labelsize=24)
        plt.rc("ytick", labelsize=24)
        plt.ylim((10**-5, 10**-0.1))
        plt.savefig(savedir + "Concentration.png")
        print(savedir + "Concentration.png")
        plt.close()

        # Revert PDFs to pre-normalized values
        PDF_velocity *= normal_velocity
        PDF_vx *= normal_vx
        PDF_vy *= normal_vy
        PDF_vorticity *= normal_vorticity
        PDF_C *= normal_C

if __name__=="__main__":
    datadir = "../data/numParticles1000-gammaB10-alpha-1-xi0.9-convection1/"
    savedir = "../photo_video/numParticles1000-gammaB10-alpha-1-xi0.9-convection1"
    PDF(datadir, savedir)
