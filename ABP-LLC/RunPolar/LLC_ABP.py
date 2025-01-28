import os
import numpy as np
class LLC_ABP:
    def __init__(self):
        # General Parameters (Values provided in your example)
        self.numParticles = 231
        self.gammaB = 10
        self.alpha = -15
        self.xi = 0.9
        self.Lx = 128
        self.Ly = 128
        self.convection = 1
        # Other parameters
        self.parameter_file = "./parameter.dat"
        self.dire_export = "../photo_video/"
        self.tStart = 0
        self.tStop = 100
        self.tStepField = 0.001
        self.tStepParticle = 0.0001
        self.tExpo = 1

        # Field parameters
        self.numGridX = 128
        self.numGridY = 128
        self.numGridBounX = 5
        self.numGridBounY = 5
        self.blockNumFields = 128
        self.threadNumFields = 128
        self.dx = 1
        self.dy = 1
        self.initCondFields = 1
        self.hLDGa = 1
        self.hLDGb = 2
        self.hLDGk = 1
        self.xiAn = 0.15
        self.eta = 1
        self.h = 5
        self.zta = 12*self.eta/self.h/self.h
        self.gammaV = 1
        self.Gamma = 1 

        # Particle parameters
        self.blockNumParticles = 128
        self.threadNumParticles = 64
        self.cellNumX = 32
        self.cellNumY = 32
        self.cellSizeX = 4
        self.cellSizeY = 4
        self.initCondParticles = 1
        self.ABParticle = 1
        self.V0 = 10
        self.boundarysize = 5
        self.kBT = 0.01
        self.TEffectiveBase = 1
        self.maxParticlePerCell = 200
        self.maxParticlePerGrid = 800
        self.maxParticleNeighbor = 800
        self.rd = 2
        self.minDistance = 1
        self.r0 = 1
        self.epsilon = 1
        self.Dc = 0.1
        self.Dp = 0.1
        self.Dparall = 5
        self.Dverti = 2
        self.seed = 231
        self.rUpdateCellList = 1


        # Generate savename, datadir, and savedir dynamically
        self.savename = f"numParticles{self.numParticles}-gammaB{self.gammaB}-alpha{self.alpha}-xi{self.xi}-convection{self.convection}-Dp{self.Dp}-kBT{self.kBT}"
        self.datadir = f"../data/{self.savename}/"
        self.savedir = f"../photo_video/{self.savename}/"

    def create_directories(self):
        # self.savename = f"numParticles{self.numParticles}-gammaB{self.gammaB}-alpha{self.alpha}-xi{self.xi}-convection{self.convection}-Dp{self.Dp}-kBT{self.kBT}"
        # self.datadir = f"../data/{self.savename}/"
        # self.savedir = f"../photo_video/{self.savename}/"
        # Create the directories if they don't exist
        os.makedirs(self.datadir, exist_ok=True)
        os.makedirs(self.savedir, exist_ok=True)

    def save_parameters(self):
        # Write parameters to file
        with open(self.parameter_file, "w") as f:
            # General Parameters
            f.write(f"parameter_file = {self.parameter_file}\n")
            f.write(f"dire_export = {self.dire_export}\n")
            f.write(f"Lx = {self.Lx}\n")
            f.write(f"Ly = {self.Ly}\n")
            f.write(f"convection = {self.convection}\n")

            # Time Parameters
            f.write(f"tStart = {self.tStart}\n")
            f.write(f"tStop = {self.tStop}\n")
            f.write(f"tStepField = {self.tStepField}\n")
            f.write(f"tStepParticle = {self.tStepParticle}\n")
            f.write(f"tExpo = {self.tExpo}\n")

            # Field Parameters
            f.write(f"numGridX = {self.numGridX}\n")
            f.write(f"numGridY = {self.numGridY}\n")
            f.write(f"numGridBounX = {self.numGridBounX}\n")
            f.write(f"numGridBounY = {self.numGridBounY}\n")
            f.write(f"blockNumFields = {self.blockNumFields}\n")
            f.write(f"threadNumFields = {self.threadNumFields}\n")
            f.write(f"dx = {self.dx}\n")
            f.write(f"dy = {self.dy}\n")
            f.write(f"initCondFields = {self.initCondFields}\n")
            f.write(f"hLDGa = {self.hLDGa}\n")
            f.write(f"hLDGb = {self.hLDGb}\n")
            f.write(f"hLDGk = {self.hLDGk}\n")
            f.write(f"xi = {self.xi}\n")
            f.write(f"xiAn = {self.xiAn}\n")
            f.write(f"alpha = {self.alpha}\n")
            f.write(f"eta = {self.eta}\n")
            f.write(f"h = {self.h}\n")
            f.write(f"zta = {self.zta}\n")
            f.write(f"gammaV = {self.gammaV}\n")
            f.write(f"Gamma = {self.Gamma}\n")
            f.write(f"gammaB = {self.gammaB}\n")
            # Particle Parameters
            f.write(f"numParticles = {self.numParticles}\n")
            f.write(f"blockNumParticles = {self.blockNumParticles}\n")
            f.write(f"threadNumParticles = {self.threadNumParticles}\n")
            f.write(f"cellNumX = {self.cellNumX}\n")
            f.write(f"cellNumY = {self.cellNumY}\n")
            f.write(f"cellSizeX = {self.cellSizeX}\n")
            f.write(f"cellSizeY = {self.cellSizeY}\n")
            f.write(f"initCondParticles = {self.initCondParticles}\n")
            f.write(f"ABParticle = {self.ABParticle}\n")
            f.write(f"V0 = {self.V0}\n")
            f.write(f"boundarysize = {self.boundarysize}\n")
            f.write(f"kBT = {self.kBT}\n")
            f.write(f"TEffectiveBase = {self.TEffectiveBase}\n")
            f.write(f"maxParticlePerCell = {self.maxParticlePerCell}\n")
            f.write(f"maxParticlePerGrid = {self.maxParticlePerGrid}\n")
            f.write(f"maxParticleNeighbor = {self.maxParticleNeighbor}\n")
            f.write(f"rd = {self.rd}\n")
            f.write(f"minDistance = {self.minDistance}\n")
            f.write(f"r0 = {self.r0}\n")
            f.write(f"epsilon = {self.epsilon}\n")
            f.write(f"Dc = {self.Dc}\n")
            f.write(f"Dp = {self.Dp}\n")
            f.write(f"Dparall = {self.Dparall}\n")
            f.write(f"Dverti = {self.Dverti}\n")
            f.write(f"savename = {self.savename}\n")
            f.write(f"datadir = {self.datadir}\n")
            f.write(f"savedir = {self.savedir}\n")
            
            f.write(f"seed = {self.seed}\n")
            f.write(f"rUpdateCellList = {self.rUpdateCellList}\n")

      

        print(f"All parameters saved to {self.parameter_file}")

    def run(self):
        # Create necessary directories
        self.create_directories()

        # Save the parameters
        self.save_parameters()
