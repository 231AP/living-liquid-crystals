#ifndef ABPARTICLECLASSGPU_CU
#define ABPARTICLECLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "ABParticleClass.h"

using namespace std; 

//================================================================
__global__ void getABParticlePxPy(Particle PT, Parameter PM, int Nx,int Ny) {
    int i = blockIdx.x;
    int j = threadIdx.x; 
    int cellId = blockDim.x*i + j;
    // printf("wrong = %d",i);
    // printf("%d",i);
    // printf("%d",j);
    if (i<Ny && j<Nx) {
        PT.cellPx[cellId] = 0;
        PT.cellPy[cellId] = 0;
        PT.cellPxx[cellId] = 0;
        PT.cellPxy[cellId] = 0;


        for (int num=1; num<=PT.cellOffsetsCL[cellId];num++){
            PT.cellPx[cellId] += abs(PT.px[PT.cellList[cellId * PM.maxParticlePerCell + num]]);
            PT.cellPy[cellId] += abs(PT.py[PT.cellList[cellId * PM.maxParticlePerCell + num]]);

        }

        if(PT.cellOffsetsCL[i*blockDim.x+j]){
        PT.cellPxx[cellId] = pow(PT.cellPx[cellId]/PT.cellOffsetsCL[cellId],2)-0.5;
        PT.cellPxy[cellId] = \
                            (PT.cellPy[cellId]/PT.cellOffsetsCL[cellId])\
                                 *(PT.cellPx[cellId]/PT.cellOffsetsCL[cellId]);
        }
        
    }
}





__global__ void updateConcentration(Particle PT,Parameter PM,double* Pxx, double* Pxy, double* Concentration, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
   
    
    if (i<Ny && j<Nx) { 
        // printf("1");

        Concentration[idx] = 0.6*PT.cellOffsetsCL[i*blockDim.x+j];
        Pxx[idx] = 0.6*PT.cellOffsetsCL[i*blockDim.x+j]*PT.cellPxx[i*blockDim.x +j];
        Pxy[idx] = 0.6*PT.cellOffsetsCL[i*blockDim.x+j]*PT.cellPxy[i*blockDim.x +j];//这里

        for (int x = -1;x <= 1;x++) {
        for (int y = -1;y <= 1;y++) {
            int kk1 =(i+x+PM.cellNumX)%PM.cellNumX;
            int kk2 = (j+y+PM.cellNumY)%PM.cellNumY;
            Concentration[idx] += 0.2*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += 0.2*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += 0.2*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];

            
        }
        }
        for (int x = -2;x <= 2;x++) {
        for (int y = -2;y <= 2;y++) {
            int kk1 =(i+x+PM.cellNumX)%PM.cellNumX;
            int kk2 = (j+y+PM.cellNumY)%PM.cellNumY;
            Concentration[idx] += 0.05*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += 0.05*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += 0.05*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
        }
        }

        for (int x = -3;x <= 3;x++) {
        for (int y = -3;y <= 3;y++) {
            int kk1 =(i+x+PM.cellNumX)%PM.cellNumX;
            int kk2 = (j+y+PM.cellNumY)%PM.cellNumY;
            Concentration[idx] += 0.05*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += 0.05*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += 0.05*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
        }
        }

        for (int x = -4;x <= 4;x++) {
        for (int y = -4;y <= 4;y++) {
            int kk1 =(i+x+PM.cellNumX)%PM.cellNumX;
            int kk2 = (j+y+PM.cellNumY)%PM.cellNumY;
            Concentration[idx] += 0.06*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += 0.06*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += 0.06*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
            
        }
        }
        if (Concentration[idx]>0){
            Pxx[idx]/= Concentration[idx];
            Pxy[idx]/= Concentration[idx];
        }
        


        // Pxx[idx]=0;
        // Pxy[idx]=0;


        // if(PT.cellOffsetsCL[i*blockDim.x+j]){
        // Pxx[idx] = pow(PT.cellPx[i*blockDim.x+j]/PT.cellOffsetsCL[i*blockDim.x+j],2)-0.5;
        // Pxy[idx] = pow(PT.cellPy[i*blockDim.x+j]/PT.cellOffsetsCL[i*blockDim.x+j],2);    
        // };
        
        // printf("%f",PT.cellPx[i*blockDim.x+j]);
        

        PT.cellPx[i*blockDim.x+j] = 0;
        PT.cellPy[i*blockDim.x+j] = 0;


  
    }
};
//=================================================================================
__global__ void initState(curandState* state,unsigned long long seed, int particleNum) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= particleNum)return;
    curand_init(seed, id, 0, &state[id]);
};
__device__ int sign(real x) {
    return -(x < 0.f) + (x > 0.f);
};

__device__ int sign01(real x) {
    return (sign(x) + 1) / 2;
}
__global__ void getCellList(Particle PT, Parameter PM) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= PM.particleNum)return;
    PT.cellX[id] = std::floor(PT.x[id] / PM.cellSizeX);
    PT.cellY[id] = std::floor(PT.y[id] / PM.cellSizeY);
    int cellId = PT.cellY[id] * PM.cellNumX + PT.cellX[id];
    int offsetsCL = atomicAdd(&PT.cellOffsetsCL[cellId], 1);
    if (offsetsCL < PM.maxParticlePerCell) {
        PT.cellList[cellId * PM.maxParticlePerCell + offsetsCL] = id;
    }
    else {
        
        // printf("wrong");//append cout error later
        printf("wrong = %d",offsetsCL);
    };
};

//===========================================================================================







__device__ int getNeighborListTry(real x0, real y0, real x1, real y1, Parameter PM) {
    real dx = sign(x1 - x0) * (x1 - x0);
    real dy = sign(y1 - y0) * (y1 - y0);
    dx = sign01(0.5 * PM.boxX - dx) * dx + sign01(dx - 0.5 * PM.boxX) * (PM.boxX - dx);
    dy = sign01(0.5 * PM.boxY - dy) * dy + sign01(dy - 0.5 * PM.boxY) * (PM.boxY - dy);
    real dr2 = dx * dx + dy * dy;
    if (dr2 < PM.rd * PM.rd)return 1;
    else return 0;
}

//==========================================================================================================
__global__ void getAroundCellParticleId(Particle PT, Parameter PM) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= PM.particleNum)return;
    int offsetPAI = 0;//particleAroundId put particleId in PAI
    int periodicBoundaryFlagX, periodicBoundaryFlagY;
    int cellXAround, cellYAround;
    int cellAroundId;
    // printf("333");
    for (int x = -1;x <= 1;x++) {
        for (int y = -1;y <= 1;y++) {
            if (PT.cellX[id] + x == -1) {
                cellXAround = PM.cellNumX - 1;
                periodicBoundaryFlagX = 1;
            }
            else if (PT.cellX[id] + x == PM.cellNumX) {
                cellXAround = 0;
                periodicBoundaryFlagX = -1;

            }
            else {
                cellXAround = PT.cellX[id] + x;
                periodicBoundaryFlagX = 0;
            }
            if (PT.cellY[id] + y == -1) {
                cellYAround = PM.cellNumY - 1;
                periodicBoundaryFlagY = 1;
            }
            else if (PT.cellY[id] + y == PM.cellNumY) {
                cellYAround = 0;
                periodicBoundaryFlagY = -1;
            }
            else {
                cellYAround = PT.cellY[id] + y;
                periodicBoundaryFlagY = 0;
            }
            int cellAroundId = cellYAround * PM.cellNumX + cellXAround;
            // printf("%d",PT.cellOffsetsCL[cellAroundId]);
            for (int i = 0;i < PT.cellOffsetsCL[cellAroundId];i++) {
                if (PT.cellList[cellAroundId * PM.maxParticlePerCell + i] == id)continue;
                int ifNeighbor = getNeighborListTry(PT.x[id], PT.y[id], PT.x[PT.cellList[cellAroundId * PM.maxParticlePerCell + i]]\
                    , PT.y[PT.cellList[cellAroundId * PM.maxParticlePerCell + i]], PM);
                if (ifNeighbor) {
                    PT.NeighborList[id * PM.maxParticlePerCell + PT.offsetsNL[id]] = PT.cellList[cellAroundId * PM.maxParticlePerCell + i];
                    PT.NeighborListFlagX[id * PM.maxParticlePerCell + PT.offsetsNL[id]] = periodicBoundaryFlagX;
                    PT.NeighborListFlagY[id * PM.maxParticlePerCell + PT.offsetsNL[id]] = periodicBoundaryFlagY;//nodebug
                    atomicAdd(&PT.offsetsNL[id], 1);
                }
	    }
        }
    }
}

//===========================================================================
__global__ void saveXY0ToUpdateHybridList(Particle PT, Parameter PM) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= PM.particleNum)return;
    PT.x0ToUpdateHybridList[id] = PT.x[id];
    PT.y0ToUpdateHybridList[id] = PT.y[id];
}
//===========================================================================

__global__ void checkUpdate(Particle PT, Parameter PM) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= PM.particleNum)return;
    real x1 = PT.x[id], x0 = PT.x0ToUpdateHybridList[id];
    real y1 = PT.y[id], y0 = PT.y0ToUpdateHybridList[id];
    real dx = sign(x1 - x0) * (x1 - x0);
    real dy = sign(y1 - y0) * (y1 - y0);
    dx = sign01(0.5 * PM.boxX - dx) * dx + sign01(dx - 0.5 * PM.boxX) * (PM.boxX - dx);
    dy = sign01(0.5 * PM.boxY - dy) * dy + sign01(dy - 0.5 * PM.boxY) * (PM.boxY - dy);
    if ((dx * dx + dy * dy) > PM.rOutUpdateList * PM.rOutUpdateList) atomicExch(&updateListFlag, 1);
}

//=========================================================================
__global__ void getForce (Particle PT, Parameter PM, double* vx, double* vy, double* Pxx, double* Pxy, double* Qxx, double* Qxy,
                         int Nx, int Ny, int Nbx, int Nby) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // printf("%d,%d,%d",id,id,PM.particleNum);
    if (id >= PM.particleNum)return;


    real x0, y0, x1, y1, dx, dy, dr, f12;


    
    PT.fx[id] = 0;
    PT.fy[id] = 0;
    int idNL;
    for (idNL = 0;idNL < PT.offsetsNL[id];idNL++) {
        x0 = PT.x[id];
        y0 = PT.y[id];
        x1 = PT.x[PT.NeighborList[id * PM.maxParticlePerCell + idNL]];
        y1 = PT.y[PT.NeighborList[id * PM.maxParticlePerCell + idNL]];
        dx = sign01(0.5 * PM.boxX - x0 + x1) * sign01(0.5 * PM.boxX + x0 - x1) * (x0 - x1) + \
            sign01(sign(x0 - x1) * (x0 - x1) - 0.5 * PM.boxX) * -sign(x0 - x1) * (PM.boxX - sign(x0 - x1) * (x0 - x1));
        dy = sign01(0.5 * PM.boxY - y0 + y1) * sign01(0.5 * PM.boxY + y0 - y1) * (y0 - y1) + \
            sign01(sign(y0 - y1) * (y0 - y1) - 0.5 * PM.boxY) * -sign(y0 - y1) * (PM.boxY - sign(y0 - y1) * (y0 - y1));
        
        dr = sqrt(dx * dx + dy * dy);



        // f12 = 1;
// f12 = 0.01/pow(dr,6);
        // f12 = 24 * PM.epsilon * pow(PM.r0, 6) * (2 * pow(PM.r0, 6) - pow(dr, 6)) / pow(dr, 14);
        f12 = 0.01/pow(dr,2);
        PT.fx[id] += f12 * dx ;
        PT.fy[id] += f12 * dy ;
        if (PT.fx[id] > 10000 || PT.fx[id] < -10000 || PT.fy[id] > 10000 || PT.fy[id] < -10000) {
            break;
        }
    }


        // int cellId = PT.cellY[id] * PM.cellNumX + PT.cellX[id];//此时未考虑周期边界条件
        int idx =  (PT.cellY[id]+Nby) * (PM.cellNumX+2*Nbx)+Nbx+PT.cellX[id];
        // printf(",,%f,,%f,,",vx[idx],PT.fx[id]);
        // printf(",,,,");
        PT.fx[id] += vx[idx];
        PT.fy[id] += vy[idx];
        double SQ = 2*sqrt(Qxx[idx]*Qxx[idx]+Qxy[idx]*Qxy[idx]);
        double theta = 0.5* atan2(Qxy[idx],Qxx[idx]);
        double theta1 = atan2(PT.py[id],PT.px[id]);
        
        // printf("%f",theta);
        theta1 += SQ*sin(2*theta - 2*theta1);
        PT.px[id] = cos(theta1);
        PT.py[id] = sin(theta1);
// printf("%f",PT.px[id]);







    if (PT.fx[id] > 10000 || PT.fx[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,x0:%f,x1:%f,y0:%f,y1:%f,NLFX:%d\n", id,PT.fx[id], PT.fy[id], dx,dy,x0,x1, y0, y1, PT.NeighborListFlagX[id * PM.maxParticlePerCell + idNL]);
        wrongFlag = 1;
    }
    if (PT.fy[id] > 10000 || PT.fy[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,y0:%f,y1:%f,x0:%f,x1:%f,NLFY:%d\n", id, PT.fx[id], PT.fy[id],dx,dy,y0,y1, x0, x1, PT.NeighborListFlagY[id * PM.maxParticlePerCell + idNL]);
        wrongFlag = 1;
    }
}

//====================================================================
__device__ real generateNormal(curandState* state) {
    return curand_normal(&(*state));
}

//==========================================================================================================
__global__ void updatePosition(Particle PT, Parameter PM) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= PM.particleNum)return;
    real fT = sqrt(2 * PM.kBT * PM.gammaValue * PM.tStep);
    real FRx = generateNormal(&PT.state[id]);
    real FRy = generateNormal(&PT.state[id]);
    PT.x[id] = fmod(PT.x[id] + (PT.fx[id] * PM.tStep  + fT * FRx + PM.V0*PT.px[id]*PM.tStep) / PM.gammaValue + PM.boxX, PM.boxX);
    PT.y[id] = fmod(PT.y[id] + (PT.fy[id] * PM.tStep + fT * FRy + PM.V0*PT.py[id]*PM.tStep)/ PM.gammaValue + PM.boxY, PM.boxY);
    // int cellId = PT.cellY[id] * PM.cellNumX + PT.cellX[id];
    
    // PT.cellPx[cellId] += PT.px[id];
    // PT.cellPy[cellId] += PT.py[id];
    // printf("%f",PT.cellPx[cellId]);

}

//=========================================================================
#endif
