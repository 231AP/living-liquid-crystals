import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
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

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def PDF(datadirs, savedir, tStart=50, tStop=90, max_F =8,max_velocity=80, max_vorticity=80, max_C=2, size=128, V0=7.5):
    id_data = -4000
    os.makedirs(savedir,exist_ok=True)
    h = 1
    for datadir in datadirs:
        unit_velocity = max_velocity / 100
        unit_vorticity = max_vorticity / 100
        unit_C = max_C / 100
        unit_F = max_F / 100
        
        # Initialize PDF arrays
        PDF_velocity = np.zeros(100)
        PDF_vx = np.zeros(100)
        PDF_vy = np.zeros(100)
        PDF_vorticity = np.zeros(100)
        PDF_C = np.zeros(100)
        PDF_vxB = np.zeros(100)
        PDF_vyB = np.zeros(100)
        PDF_vB = np.zeros(100)
        PDF_test = np.zeros(100)
        PDF_test1 = np.zeros(100)
        PDF_deltaTheta = np.zeros(100)
        PDF_actForce = np.zeros(100)
        ALLv2B = 0
        ALLv4B = 0
        ALLv4 = 0
        ALLv2 = 0
        ALLdv4 = 0
        ALLdv2 = 0
        
        # Loop through the time steps and process data
        for time in np.arange(tStart, tStop):
            Qxx_df = pd.read_csv(datadir + f"abParticle.Qxx_{time}.dat",header = None)
            Qxy_df = pd.read_csv(datadir + f"abParticle.Qxy_{time}.dat",header = None)
            Pxx_df = pd.read_csv(datadir + f"abParticle.Pxx_{time}.dat",header = None)
            Pxy_df = pd.read_csv(datadir + f"abParticle.Pxy_{time}.dat",header = None)
            omega_df = pd.read_csv(datadir + f"incompFlow.omega_{time}.dat", header=None)
            vx_df = pd.read_csv(datadir + f"incompFlow.vx_{time}.dat", header=None)
            vy_df = pd.read_csv(datadir + f"incompFlow.vy_{time}.dat", header=None)
            C_df = pd.read_csv(datadir + f"abParticle.Concentration_{time}.dat", header=None)
            kk = pd.read_csv(datadir + f"conf_{time}.dat", sep=" ", header=None)
            pts = np.array(kk)
            
            Qxx = np.array(Qxx_df).reshape(size,size)
            Qxy = np.array(Qxy_df).reshape(size,size)
            Pxx = np.array(Pxx_df).reshape(size,size)
            Pxy = np.array(Pxy_df).reshape(size,size)
            
            omega = np.array(omega_df).reshape(size, size)
            vx = np.array(vx_df).reshape(size, size)
            vy = np.array(vy_df).reshape(size, size)
            C = np.array(C_df).reshape(size, size)
            # print(np.random.normal(len(pts[:,2])))
            mean = 1
            std_dev = 0.2
            num_samples = len(pts[:,2])
# Generate random numbers
            random_numbers = np.random.normal(loc=mean, scale=std_dev, size=num_samples)
            vx_bacteria = pts[:, 2] * V0 
            vy_bacteria = pts[:, 3] * V0 
            
            for idB in np.arange(len(vx_bacteria)):
                vx_bacteria[idB] += vx[int(pts[idB, 0]), int(pts[idB, 1])]
                vy_bacteria[idB] += vy[int(pts[idB, 0]), int(pts[idB, 1])]
            vB = np.sqrt(vx_bacteria**2 + vy_bacteria**2)
            v = np.sqrt(vx**2 + vy**2)
            # actFx = d1xO4(C * Pxx,h) + d1yO4(C*Pxy,h)
            # actFy = d1xO4(C*Pxy,h) - d1yO4(C*Pxx,h)
            # actF = np.sqrt(actFx**2+actFy**2)
            actFx = d1xO4(C*Pxx,h) + d1yO4(C*Pxy,h)
            actFy = d1xO4(C*Pxy,h) - d1yO4(C*Pxx,h)
            actF = np.sqrt(actFx**2+actFy**2)
            ALLv4B += np.mean(vx_bacteria**4)
            ALLv2B += np.mean(vx_bacteria**2)**2
            kurtosisB = ALLv4B / ALLv2B
            print("kurtosis_vxB = %f" % kurtosisB)
            
            dvx = d1xF(vx, 1)
            ALLv4 += np.mean(vx**4)
            ALLv2 += np.mean(vx**2)**2
            kurtosis = ALLv4 / ALLv2
            print("kurtosis_vx = %f" % kurtosis)
            
            ALLdv4 += np.mean(dvx**4)
            ALLdv2 += np.mean(dvx**2)**2
            kurtosis_dvx = ALLdv4 / ALLdv2
            print("kurtosis_dvx = %f" % kurtosis_dvx)
            
            # Calculate PDF values
            for j in range(100):
                i = j - 50
                PDF_velocity[j] += np.sum((v < i * unit_velocity) * (v > (i - 1) * unit_velocity))
                PDF_vx[j] += np.sum((vx < i * unit_velocity) * (vx > (i - 1) * unit_velocity))
                PDF_vy[j] += np.sum((vy < i * unit_velocity) * (vy > (i - 1) * unit_velocity))
                PDF_vB[j] += np.sum((vB < i * unit_velocity) * (vB > (i - 1) * unit_velocity))
                PDF_vxB[j] += np.sum((vx_bacteria < i * unit_velocity) * (vx_bacteria > (i - 1) * unit_velocity))
                PDF_vyB[j] += np.sum((vy_bacteria < i * unit_velocity) * (vy_bacteria > (i - 1) * unit_velocity))
                PDF_test[j] += np.sum((pts[:, 2] * V0 < i * unit_velocity) * (pts[:, 2] * V0 > (i - 1) * unit_velocity))
                PDF_test1[j] += np.sum((pts[:, 3] * V0 < i * unit_velocity) * (pts[:, 3] * V0 > (i - 1) * unit_velocity))
                PDF_vorticity[j] += np.sum((omega < i * unit_vorticity) * (omega > (i - 1) * unit_vorticity))
                PDF_C[j] += np.sum((C < i * unit_C) * (C > (i - 1) * unit_C))
                # PDF_actForceQ[j] += np.sum((actFQ < i*unit_F) * (actFQ>(i-1)*unit_F))
                PDF_actForce[j] += np.sum((actF < i*unit_F) * (actF>(i-1)*unit_F))

            # Normalize the PDFs
            normal_velocity = np.sum(PDF_velocity)
            normal_vx = np.sum(PDF_vx)
            normal_vy = np.sum(PDF_vy)
            normal_vB = np.sum(PDF_vB)
            normal_vxB = np.sum(PDF_vxB)
            normal_vyB = np.sum(PDF_vyB)
            normal_test = np.sum(PDF_test)
            normal_test1 = np.sum(PDF_test1)
            normal_vorticity = np.sum(PDF_vorticity)
            normal_C = np.sum(PDF_C)
            normal_actForce = np.sum(PDF_actForce)

            PDF_velocity /= normal_velocity
            PDF_vx /= normal_vx
            PDF_vy /= normal_vy
            PDF_vB /= normal_vB
            PDF_vxB /= normal_vxB
            PDF_vyB /= normal_vyB
            PDF_test /= normal_test
            PDF_test1 /= normal_test1
            PDF_vorticity /= normal_vorticity
            PDF_C /= normal_C
            PDF_actForce /= normal_actForce

            # Set font properties
            font = FontProperties()
            font.set_size("xx-large")
            if time == tStop-1:
                if id_data == -4000:
            # Create a single figure with multiple subplots
                    fig, axes = plt.subplots(2, 3, figsize=(18, 12))  # 2 rows and 3 columns of subplots
            # Plot Velocity PDF (fluid vs. bacteria vs. velocity)
                id_data+=6000
                print(id_data)
                axes[0, 0].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vx, label=f"u{id_data}")
                axes[0, 0].scatter(np.arange(-50, 50) / 100 * max_velocity, PDF_vx)
                axes[0, 0].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vy, label=f"v{id_data}")
                axes[0, 0].scatter(np.arange(-50, 50) / 100 * max_velocity, PDF_vy)
                axes[0, 0].legend(loc=1, prop=font)
                axes[0, 0].set_yscale("log")
                axes[0, 0].set_xlabel("Velocity", fontsize=24)
                axes[0, 0].set_ylabel("PDF", fontsize=24)
                axes[0, 0].set_ylim((10**-5, 10**-0.1))

                # Plot Velocity PDF (bacteria vs. test)
                axes[0, 1].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vxB, label=f"u{id_data}")
                axes[0, 1].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vyB, label=f"v{id_data}")
                axes[0, 1].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_test, label=f"test{id_data}")
                axes[0, 1].legend(loc=1, prop=font)
                axes[0, 1].set_yscale("log")
                axes[0, 1].set_xlabel("Velocity", fontsize=24)
                axes[0, 1].set_ylabel("PDF", fontsize=24)
                axes[0, 1].set_ylim((10**-5, 10**-0.1))

                # Plot Velocity PDF (fluid vs. bacteria)
                axes[0, 2].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_velocity, label=f"fluid{id_data}")
                axes[0, 2].plot(np.arange(-50, 50) / 100 * max_velocity, PDF_vB, label=f"bacteria{id_data}")
                axes[0, 2].legend(loc=1, prop=font)
                axes[0, 2].set_yscale("log")
                axes[0, 2].set_xlabel("Velocity", fontsize=24)
                axes[0, 2].set_ylabel("PDF", fontsize=24)
                axes[0, 2].set_ylim((10**-5, 10**-0.1))
                axes[0, 2].set_xlim(0, 30)

                # Plot Vorticity PDF
                axes[1, 0].plot(np.arange(-50, 50) / 100 * max_vorticity, PDF_vorticity, label=f"vorticity{id_data}")
                axes[1, 0].legend(loc=1, prop=font)
                axes[1, 0].set_yscale("log")
                axes[1, 0].set_xlabel("Vorticity", fontsize=24)
                axes[1, 0].set_ylabel("PDF", fontsize=24)
                axes[1, 0].set_ylim((10**-5, 10**-0.1))

                # Plot Concentration PDF
                axes[1, 1].plot(np.arange(-50, 50) / 100 * max_C, PDF_C, label=f"Concentration{id_data}")
                axes[1, 1].scatter(np.arange(-50, 50) / 100 * max_C, PDF_C, label=f"Concentration{id_data}")
                axes[1, 1].legend(loc=1, prop=font)
                axes[1, 1].set_yscale("log")
                axes[1, 1].set_xlim((0, max_C/2))
                axes[1, 1].set_xlabel("Concentration", fontsize=24)
                axes[1, 1].set_ylabel("PDF", fontsize=24)
                axes[1, 1].set_ylim((10**-5, 10**-0.1))

                
                axes[1, 2].plot(np.arange(-50, 50) / 100 * max_F, PDF_actForce, label=f"activeForce{id_data}")
                axes[1, 2].scatter(np.arange(-50, 50) / 100 * max_F, PDF_actForce, label=f"activeForce{id_data}")
                axes[1, 2].legend(loc=1, prop=font)
                axes[1, 2].set_yscale("log")
                axes[1, 2].set_xlim((0, max_F/2))
                axes[1, 2].set_xlabel("ActiveForce", fontsize=24)
                axes[1, 2].set_ylabel("PDF", fontsize=24)
                axes[1, 2].set_ylim((10**-5, 10**-0.1))
                
                
                # axes[0,0].plot(np.arange(-50, 50) / 100 * max_velocity*2, PDF_actForce, label=f"activeForce{id_data}")
                # axes[0, 0].scatter(np.arange(-50, 50) / 100 * max_velocity*2, PDF_actForce, label=f"activeForce{id_data}")
                # axes[0, 0].legend(loc=1, prop=font)
                # axes[0, 0].set_yscale("log")
                # # axes[1, 2].set_xlim((0, max_F/2))
                # axes[0, 0].set_xlabel("ActiveForce", fontsize=24)
                # axes[0, 0].set_ylabel("PDF", fontsize=24)
                # axes[0,0 ].set_ylim((10**-5, 10**-0.1))
                # Hide the last empty subplot (since we have 5 plots in total)
                # axes[2, 0].axis('off')
                # axes[2, 1].axis('off')
                # axes[2, 2].axis('off')
                plt.tight_layout()  # Automatically adjust subplot layout
                print(1)
                plt.savefig(savedir + "All_PDFs.png")
            # plt.close()

            # Revert PDFs to pre-normalized values
            PDF_velocity *= normal_velocity
            PDF_vx *= normal_vx
            PDF_vy *= normal_vy
            PDF_vorticity *= normal_vorticity
            PDF_C *= normal_C
            PDF_vB *= normal_vB
            PDF_vxB *= normal_vxB
            PDF_vyB *= normal_vyB
            PDF_test *= normal_test
            PDF_test1 *= normal_test1
            PDF_actForce *= normal_actForce
if __name__=="__main__":
    datadir = "../data/numParticles2000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir1 = "../data/numParticles4000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir2 = "../data/numParticles6000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    datadir3 = "../data/numParticles8000-gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5/"
    
    datadirs = [datadir,datadir1,datadir2,datadir3]
    datadirs = [datadir,datadir3]
    savedir = "../photo_video/analysis/"+ "example3"+"gammaB3-alpha-20-xi0.9-convection1-Dp0-kBT0.5"
    PDF(datadirs,savedir,max_velocity=40,max_vorticity=40,max_C=10)


