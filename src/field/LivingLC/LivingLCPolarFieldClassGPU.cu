#ifndef LIVINGLCPOLARFIELDCLASSGPU_CU
#define LIVINGLCPOLARFIELDCLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "LivingLCPolarFieldClass.h"

using namespace std; 

//================================================================
__global__ void getLivingLCPxPyThetaGPU(double* flip,double* Pxx, double* Pxy, double* px, double* py, double* theta, double* theta_old, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
    
    if (i<Ny && j<Nx) {
      // flip[idx] = px[idx];
      theta_old[idx] = theta[idx];
      // double p = sqrt(px[idx]*px[idx]+py[idx]*py[idx]);

      // Pxx[idx] = px[idx]*px[idx]- p/2;
      // Pxy[idx] = py[idx]*py[idx];
      theta[idx] =0.5* atan2(Pxy[idx],Pxx[idx]);
   
      double p = 2*sqrt(Pxx[idx]*Pxx[idx]+Pxy[idx]*Pxy[idx]);
      px[idx] = p*cos(theta[idx]);
      py[idx] = p*sin(theta[idx]);
      
    };
};


//================================================================
__global__ void getLivingLCFlipGPU(double* Omega,double* cplus, double* cminus,double* theta_old, double* theta,double* flip, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
    int di = (blockDim.x+2*Nbx);
    if (i<Ny && j<Nx) {
      if (abs(theta[idx+di]-theta[idx-di])>2.5) {
          flip[idx] = 1;
          // theta_old[idx-di] = cplus[idx-di];
          // cplus[idx-di] = cminus[idx-di];
          // cminus[idx-di] = theta_old[idx-di];
          Omega[idx-di] = -Omega[idx-di];
       
          
      };

    };

};


#endif
