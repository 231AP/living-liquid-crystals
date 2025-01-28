#ifndef ABPARTICLECLASSGPU_CU
#define ABPARTICLECLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "ABParticleClass.h"

using namespace std; 

//================================================================

__global__ void test11() {

    printf("1");
   
}
// =====
__global__ void getABParticlePxPy(Particle PT, LLCParams &params, int Nx,int Ny) {
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


        for (int num=0; num<PT.cellOffsetsCL[cellId];num++){
            PT.cellPx[cellId] += PT.px[PT.cellList[cellId * params.maxParticlePerGrid + num]];
            PT.cellPy[cellId] += PT.py[PT.cellList[cellId * params.maxParticlePerGrid + num]];

        }

        if(PT.cellOffsetsCL[i*blockDim.x+j]){
        PT.cellPxx[cellId] = pow(PT.cellPx[cellId]/PT.cellOffsetsCL[cellId],2)-0.5;
        PT.cellPxy[cellId] = \
                            (PT.cellPy[cellId]/PT.cellOffsetsCL[cellId])\
                                 *(PT.cellPx[cellId]/PT.cellOffsetsCL[cellId]);
        }
        
    }
}


__global__ void findConcentration(Particle PT,LLCParams &params,double* C1, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
   
    
    if (i<Ny && j<Nx) { 

        float a0 = 1;
        C1[idx] += a0*PT.cellOffsetsCL[i*blockDim.x+j];
  
    }
};





__global__ void updateConcentration(Particle PT,LLCParams &params,double* Pxx, double* Pxy,double* C1, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
   
    
    if (i<Ny && j<Nx) { 

        float a0 = 0.8;
        float a1 = 0.1;
        float a2 = 0.025;
        float a3 = 0.012;
        float a4 = 0.007;
        float a5 = 0.004;
        // printf("1");

        C1[idx] = a0*PT.cellOffsetsCL[i*blockDim.x+j];
        Pxx[idx] = PT.cellPxx[i*blockDim.x +j];
        Pxy[idx] = PT.cellPxy[i*blockDim.x +j];//这里

        for (int x = -1;x <= 1;x++) {
        for (int y = -1;y <= 1;y++) {
            int kk1 =(i+x+params.numGridX)%params.numGridX;
            int kk2 = (j+y+params.numGridY)%params.numGridY;
            C1[idx] += a1*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += a1*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += a1*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
        }
        }
        for (int x = -2;x <= 2;x++) {
        for (int y = -2;y <= 2;y++) {
            int kk1 =(i+x+params.numGridX)%params.numGridX;
            int kk2 = (j+y+params.numGridY)%params.numGridY;
            C1[idx] += a2*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += a2*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += a2*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
        }
        }

        for (int x = -3;x <= 3;x++) {
        for (int y = -3;y <= 3;y++) {
            int kk1 =(i+x+params.numGridX)%params.numGridX;
            int kk2 = (j+y+params.numGridY)%params.numGridY;
            C1[idx] += a3*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += a3*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += a3*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
        }
        }

        for (int x = -4;x <= 4;x++) {
        for (int y = -4;y <= 4;y++) {
            int kk1 =(i+x+params.numGridX)%params.numGridX;
            int kk2 = (j+y+params.numGridY)%params.numGridY;
            C1[idx] += a4*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += a4*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += a4*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
            
        }
        }

        for (int x = -5;x <= 5;x++) {
        for (int y = -5;y <= 5;y++) {
            int kk1 =(i+x+params.numGridX)%params.numGridX;
            int kk2 = (j+y+params.numGridY)%params.numGridY;
            C1[idx] += a5*PT.cellOffsetsCL[(kk1)*blockDim.x+(kk2)];
            Pxx[idx] += a5*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxx[(kk1)*blockDim.x+(kk2)];
            Pxy[idx] += a5*PT.cellOffsetsCL[kk1*blockDim.x+kk2]*PT.cellPxy[(kk1)*blockDim.x+(kk2)];
            
        }
        }

       
        if (C1[idx]>0){
            Pxx[idx]/= C1[idx];
            Pxy[idx]/= C1[idx];
        }
        PT.cellPx[i*blockDim.x+j] = 0;
        PT.cellPy[i*blockDim.x+j] = 0;


  
    }
};



// =======================================================================================



__device__ real LaplO4I(real * u, int di, int idx) {
    // int idx=(p.Nx+2*p.Nb)*(i+p.Nb)+j+p.Nb;
    int dj=1;

    return 1.0/pow(1,2)*( -21.0/5.0*u[idx]
    +13.0/15.0*( u[idx+di] + u[idx-di] + u[idx+dj] + u[idx-dj] )
    +4.0/15.0*( u[idx+di+dj] + u[idx-di+dj] + u[idx+di-dj] + u[idx-di-dj] )
    -1.0/60.0*( u[idx+2*di] + u[idx-2*di] + u[idx+2*dj] + u[idx-2*dj] )
    -1.0/30.0*( u[idx+di+2*dj] + u[idx-di+2*dj] + u[idx+di-2*dj] + u[idx-di-2*dj]
    + u[idx+2*di+dj] + u[idx-2*di+dj] + u[idx+2*di-dj] + u[idx-2*di-dj] )
    );

}


__global__ void getConcentration(double* Concentration,double* C1, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;  
    if (i<Ny && j<Nx) { 
        Concentration[idx] = C1[idx];
    }
};


__global__ void smoothConcentration(LLCParams &params,double* C1,double* Pxx,double* Pxy, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
    int di = blockDim.x+2*Nbx;


    if (i<Ny && j<Nx) { 
        C1[idx] +=   params.Dc * LaplO4I(C1,di,idx);
        Pxx[idx] +=  params.Dp * LaplO4I(Pxx,di,idx);
        Pxy[idx] +=  params.Dp * LaplO4I(Pxy,di,idx);
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
__global__ void getCellList(Particle PT, LLCParams params) {
    // printf("1");
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= params.numParticles)return;
    // printf("1");
    PT.cellX1[id] = std::floor(PT.x[id] / params.cellSizeX);
    PT.cellY1[id] = std::floor(PT.y[id] / params.cellSizeY);
    int cellId1 = PT.cellY1[id] * params.cellNumX + PT.cellX1[id];
    int offsetsCL1 = atomicAdd(&PT.cellOffsetsCL1[cellId1], 1);
    if (offsetsCL1 < params.maxParticlePerCell) {
        PT.cellList1[cellId1 * params.maxParticlePerCell + offsetsCL1] = id;
    }
    else {
        
        // printf("wrong");//append cout error later
        printf("wrong = %d",offsetsCL1);
    }

    PT.cellX[id] = std::floor(PT.x[id] / params.dx);
    PT.cellY[id] = std::floor(PT.y[id] / params.dy);

    
    int cellId = PT.cellY[id] * params.numGridX + PT.cellX[id];
    int offsetsCL = atomicAdd(&PT.cellOffsetsCL[cellId], 1);
    if (offsetsCL < params.maxParticlePerGrid) {
        PT.cellList[cellId * params.maxParticlePerGrid + offsetsCL] = id;
    }
    else {
        
        // printf("wrong");//append cout error later
        printf("wrong = %d",offsetsCL);
    }
};

//===========================================================================================



__device__ int getNeighborListTry(real x0, real y0, real x1, real y1, LLCParams &params) {
    real dx = sign(x1 - x0) * (x1 - x0);
    real dy = sign(y1 - y0) * (y1 - y0);
    dx = sign01(0.5 * params.Lx - dx) * dx + sign01(dx - 0.5 * params.Lx) * (params.Lx - dx);
    dy = sign01(0.5 * params.Ly - dy) * dy + sign01(dy - 0.5 * params.Ly) * (params.Ly - dy);
    real dr2 = dx * dx + dy * dy;
    if (dr2 < (params.rd+2*params.rUpdateCellList) * (params.rd+2*params.rUpdateCellList))return 1;
    else return 0;
}




//==========================================================================================================
__global__ void getAroundCellParticleId(Particle PT, LLCParams &params) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= params.numParticles)return;
    int offsetPAI = 0;//particleAroundId put particleId in PAI
    int periodicBoundaryFlagX, periodicBoundaryFlagY;
    int cellXAround, cellYAround;
    int cellAroundId;
    // printf("333");
    for (int x = -1;x <= 1;x++) {
        for (int y = -1;y <= 1;y++) {
            if (PT.cellX[id] + x == -1) {
                // printf("1");
                cellXAround = params.numGridX - 1;
                periodicBoundaryFlagX = 1;
            }
            else if (PT.cellX[id] + x == params.numGridX) {
                cellXAround = 0;
                periodicBoundaryFlagX = -1;

            }
            else {
                cellXAround = PT.cellX[id] + x;
                periodicBoundaryFlagX = 0;
            }
            if (PT.cellY[id] + y == -1) {
                cellYAround = params.numGridY - 1;
                periodicBoundaryFlagY = 1;
            }
            else if (PT.cellY[id] + y == params.numGridY) {
                cellYAround = 0;
                periodicBoundaryFlagY = -1;
            }
            else {
                cellYAround = PT.cellY[id] + y;
                periodicBoundaryFlagY = 0;
            }
            int cellAroundId = cellYAround * params.numGridX + cellXAround;
            // printf("%d",PT.cellOffsetsCL[cellAroundId]);
            for (int i = 0;i < PT.cellOffsetsCL[cellAroundId];i++) {
                if (PT.cellList[cellAroundId * params.maxParticlePerGrid + i] == id)continue;
                int ifNeighbor = getNeighborListTry(PT.x[id], PT.y[id], PT.x[PT.cellList[cellAroundId * params.maxParticlePerGrid + i]]\
                    , PT.y[PT.cellList[cellAroundId * params.maxParticlePerGrid + i]], params);
                if (ifNeighbor) {
                    PT.NeighborList[id * params.maxParticlePerGrid + PT.offsetsNL[id]] = PT.cellList[cellAroundId * params.maxParticlePerGrid + i];
                    PT.NeighborListFlagX[id * params.maxParticlePerGrid + PT.offsetsNL[id]] = periodicBoundaryFlagX;
                    PT.NeighborListFlagY[id * params.maxParticlePerGrid + PT.offsetsNL[id]] = periodicBoundaryFlagY;//nodebug
                    atomicAdd(&PT.offsetsNL[id], 1);
                }
	    }
        }
    }

int offsetPAI1 = 0;//particleAroundId put particleId in PAI
    int periodicBoundaryFlagX1, periodicBoundaryFlagY1;
    int cellXAround1, cellYAround1;
    int cellAroundId1;
    // printf("333");
    for (int x = -1;x <= 1;x++) {
        for (int y = -1;y <= 1;y++) {
            if (PT.cellX1[id] + x == -1) {
                cellXAround1 = params.cellNumX - 1;
                periodicBoundaryFlagX1 = 1;
            }
            else if (PT.cellX1[id] + x == params.cellNumX) {
                cellXAround1 = 0;
                periodicBoundaryFlagX1 = -1;

            }
            else {
                cellXAround1 = PT.cellX1[id] + x;
                periodicBoundaryFlagX1 = 0;
            }
            if (PT.cellY1[id] + y == -1) {
                cellYAround1 = params.cellNumY - 1;
                periodicBoundaryFlagY1 = 1;
            }
            else if (PT.cellY1[id] + y == params.cellNumY) {
                cellYAround1 = 0;
                periodicBoundaryFlagY1 = -1;
            }
            else {
                cellYAround1 = PT.cellY1[id] + y;
                periodicBoundaryFlagY1 = 0;
            }
            int cellAroundId1 = cellYAround1 * params.cellNumX + cellXAround1;
            // printf("%d",PT.cellOffsetsCL[cellAroundId]);
            for (int i = 0;i < PT.cellOffsetsCL1[cellAroundId1];i++) {
                if (PT.cellList1[cellAroundId1 * params.maxParticlePerCell + i] == id)continue;
                int ifNeighbor1 = getNeighborListTry(PT.x[id], PT.y[id], PT.x[PT.cellList1[cellAroundId1 * params.maxParticlePerCell + i]]\
                    , PT.y[PT.cellList1[cellAroundId1 * params.maxParticlePerCell + i]], params);
                if (ifNeighbor1) {
                    PT.NeighborList1[id * params.maxParticlePerCell + PT.offsetsNL1[id]] = PT.cellList1[cellAroundId1 * params.maxParticlePerCell + i];
                    PT.NeighborListFlagX1[id * params.maxParticlePerCell + PT.offsetsNL1[id]] = periodicBoundaryFlagX1;
                    PT.NeighborListFlagY1[id * params.maxParticlePerCell + PT.offsetsNL1[id]] = periodicBoundaryFlagY1;//nodebug
                    atomicAdd(&PT.offsetsNL1[id], 1);
                }
	    }
        }
    }



}

//===========================================================================
__global__ void saveXY0ToUpdateHybridList(Particle PT, LLCParams &params) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= params.numParticles)return;
    PT.x0ToUpdateHybridList[id] = PT.x[id];
    PT.y0ToUpdateHybridList[id] = PT.y[id];

    PT.x0ToUpdateHybridList1[id] = PT.x[id];
    PT.y0ToUpdateHybridList1[id] = PT.y[id];


}
//===========================================================================

__global__ void checkUpdate(Particle PT, LLCParams &params) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= params.numParticles)return;
    real x1 = PT.x[id], x0 = PT.x0ToUpdateHybridList[id];
    real y1 = PT.y[id], y0 = PT.y0ToUpdateHybridList[id];
    real dx = sign(x1 - x0) * (x1 - x0);
    real dy = sign(y1 - y0) * (y1 - y0);
    dx = sign01(0.5 * params.Lx - dx) * dx + sign01(dx - 0.5 * params.Lx) * (params.Lx - dx);
    dy = sign01(0.5 * params.Ly - dy) * dy + sign01(dy - 0.5 * params.Ly) * (params.Ly - dy);
    if ((dx * dx + dy * dy) > params.rUpdateCellList * params.rUpdateCellList) atomicExch(&updateListFlag, 1);

}
//====================================================================
__device__ real generateNormal(curandState* state) {
    return curand_normal(&(*state));
}

//=========================================================================
__global__ void getForce (Particle PT, LLCParams &params, double* vx, double* vy, double* Pxx, double* Pxy, double* Qxx, double* Qxy,
                         int Nx, int Ny, int Nbx, int Nby) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    
    printf("%d,%d,%d",id,id,params.numParticles);
    if (id >= params.numParticles)return;


    real x0, y0, x1, y1, dx, dy, dr, f12;

    
    PT.fx[id] = 0;
    PT.fy[id] = 0;
    int idNL;
    for (idNL = 0;idNL < PT.offsetsNL1[id];idNL++) {
        x0 = PT.x[id];
        y0 = PT.y[id];
        x1 = PT.x[PT.NeighborList1[id * params.maxParticlePerCell + idNL]];
        y1 = PT.y[PT.NeighborList1[id * params.maxParticlePerCell + idNL]];
        dx = sign01(0.5 * params.Lx - x0 + x1) * sign01(0.5 * params.Lx + x0 - x1) * (x0 - x1) + \
            sign01(sign(x0 - x1) * (x0 - x1) - 0.5 * params.Lx) * -sign(x0 - x1) * (params.Lx - sign(x0 - x1) * (x0 - x1));
        dy = sign01(0.5 * params.Ly - y0 + y1) * sign01(0.5 * params.Ly + y0 - y1) * (y0 - y1) + \
            sign01(sign(y0 - y1) * (y0 - y1) - 0.5 * params.Ly) * -sign(y0 - y1) * (params.Ly - sign(y0 - y1) * (y0 - y1));
        // printf("1");
        dr = sqrt(dx * dx + dy * dy);
        if(dr<params.rd){
            f12 = 1*(params.rd-dr);
            // f12 = 0;
            PT.fx[id] += f12 * dx;
            PT.fy[id] += f12 * dy;
        }else f12 = 0;

        if (PT.fx[id] > 10000 || PT.fx[id] < -10000 || PT.fy[id] > 10000 || PT.fy[id] < -10000) {
            break;
        }
    }

        int idx =  (PT.cellY[id]+Nby) * (params.numGridX+2*Nbx)+Nbx+PT.cellX[id];
        PT.fx[id] += vx[idx];
        PT.fy[id] += vy[idx];
        double SQ = 2*sqrt(Qxx[idx]*Qxx[idx]+Qxy[idx]*Qxy[idx]);
        double theta = 0.5* atan2(Qxy[idx],Qxx[idx]);
        double theta1 = atan2(PT.py[id],PT.px[id]);


        theta1 += params.gammaB*SQ*sin(2*theta - 2*theta1)*params.tStepParticle;
        PT.px[id] = cos(theta1);
        PT.py[id] = sin(theta1);
        // if (generateNormal(&PT.state[id])>8){
        //     // printf("%f",generateNormal(&PT.state[id]));
        //     PT.px[id] *= -1;
        //     PT.py[id] *= -1;
        // }
    if (PT.fx[id] > 10000 || PT.fx[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,x0:%f,x1:%f,y0:%f,y1:%f,NLFX:%d\n", id,PT.fx[id], PT.fy[id], dx,dy,x0,x1, y0, y1, PT.NeighborListFlagX[id * params.maxParticlePerGrid + idNL]);
        wrongFlag = 1;
    }
    if (PT.fy[id] > 10000 || PT.fy[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,y0:%f,y1:%f,x0:%f,x1:%f,NLFY:%d\n", id, PT.fx[id], PT.fy[id],dx,dy,y0,y1, x0, x1, PT.NeighborListFlagY[id * params.maxParticlePerGrid + idNL]);
        wrongFlag = 1;
    }
}


//==========================================================================================================
__global__ void updatePosition(Particle PT, LLCParams &params) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= params.numParticles)return;
    real fT = sqrt(2 * params.kBT *  params.gammaV* params.tStepParticle);
    real fT1 = 0.2*sqrt(2 * params.kBT * params.gammaV * params.tStepParticle);
    real FRx = generateNormal(&PT.state[id]);
    real FRy = generateNormal(&PT.state[id]);
    real PRx = generateNormal(&PT.state[id]);
    real PRy = generateNormal(&PT.state[id]);
    // printf("123");
    PT.px[id] += fT1*PRx;
    // printf("%f",PT.px[id]);
    PT.py[id] += fT1*PRy;
    PT.x[id] = fmod(PT.x[id] + (PT.fx[id] * params.tStepParticle  + fT * FRx + params.V0*PT.px[id]*params.tStepParticle) / params.gammaV + params.Lx, params.Lx);
    // printf("fx%f,FRx%f,px%f",PT.fx[id],fT,PT.px[id]);
    PT.y[id] = fmod(PT.y[id] + (PT.fy[id] * params.tStepParticle + fT * FRy + params.V0*PT.py[id]*params.tStepParticle)/ params.gammaV + params.Ly, params.Ly);
    // int cellId = PT.cellY[id] * params.numGridX + PT.cellX[id];
    
    // PT.cellPx[cellId] += PT.px[id];
    // PT.cellPy[cellId] += PT.py[id];
    // printf("%f",PT.cellPx[cellId]);

}

// ------------------------------------------------------------------



//=========================================================================
#endif
