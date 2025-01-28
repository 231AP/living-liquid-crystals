#ifndef USERDEFINEDFUNCTION_H
#define USERDEFINEDFUNCTION_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "../src/field/fieldFunction.h"

// ===============================================================
// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double getPressure(double * f, int idx, FFuncArgs f_func_args) {
    double func=0;
    if (f[idx]>f_func_args.f_func_arg2) {
        func=f_func_args.f_func_arg1*(f[idx]-f_func_args.f_func_arg2);
    }
    return func;
};
__device__ FFuncType getPressure_dev=getPressure;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double oneOverPhi(double * f, int idx, FFuncArgs f_func_args) {
    double func=f[idx];
    if (f[idx]>f_func_args.f_func_arg1) {
        func=f_func_args.f_func_arg1;
    }
    return func;
};
__device__ FFuncType oneOverPhi_dev=oneOverPhi;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double changeF(double * f, int idx, FFuncArgs f_func_args) {
    
    return  (double(3.14/2.0) < abs(f[idx])) - (abs(f[idx]) < double(3.14/2.0));

};
__device__ FFuncType changeF_dev=changeF;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double atan2F(double * f, int idx, FFuncArgs f_func_args) {
   
    /* return atan2(f[idx],f_func_args.field2[idx]); */
  return 1+600*f[idx]*exp(-20*f[idx]);
};
__device__ FFuncType atan2F_dev=atan2F;


// ---------------------------------------------------------
CUDA_CALLABLE_MEMBER double getAlpha(double * f, int idx, FFuncArgs f_func_args) {
   
  /* return 1.0+1.0/(f[idx]+0.1); */
  return 2;
};
__device__ FFuncType getAlpha_dev=getAlpha;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double gausF(double * f, int idx, FFuncArgs f_func_args) {
    double e1 = 1;
    return exp(-f[idx]*f[idx]*e1);
};
__device__ FFuncType gausF_dev=gausF;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1yRL1(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    /* int dj=f_func_args.dj; */
    double kk=1;

    if ( abs(f_func_args.field2[idx+di]-f_func_args.field2[idx-di]) > 2.9) {
        kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
    };
    return 1.0/(2.0*dy)*(f[idx+di] - kk*f[idx-di]);
};




// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1xCO2RLI2D(double * f, int idx, FFuncArgs f_func_args) {
    /* di-y, dj-x */
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    // double kk = 1
    double kk_1 = 1;
    double kk_11 = 1;

// if ( abs(f_func_args.field2[idx+di]-f_func_args.field2[idx-di]) > 2.9) {
//         /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
//         // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
//       kk=-1;
//     };

if ( abs(f_func_args.field2[idx+dj+di]-f_func_args.field2[idx+dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_1=-1;
    };

if ( abs(f_func_args.field2[idx-dj+di]-f_func_args.field2[idx-dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_11=-1;
    };

    return 1.0/(12.0*dx)*(
        (f[idx+di+dj] - f[idx+di-dj])
        +4*(f[idx+dj] - f[idx-dj])
        +(kk_1*f[idx-di+dj] - kk_11*f[idx-di-dj]));
};

__device__ FFuncType d1xCO2RLI2D_dev=d1xCO2RLI2D;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1yRL(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    /* int dj=f_func_args.dj; */
    double kk=1;

    if ( abs(f_func_args.field2[idx+di]-f_func_args.field2[idx-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk=-1;
    };
    return 1.0/(2.0*dy)*(f[idx+di] - kk*f[idx-di]);
};
__device__ FFuncType d1yRL_dev=d1yRL;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1yCO2RLI2D(double * f, int idx, FFuncArgs f_func_args) {
    /* di-y, dj-x */
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    double kk = 1;
    double kk_1 = 1;
    double kk_11 = 1;

if ( abs(f_func_args.field2[idx+di]-f_func_args.field2[idx-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk=-1;
    };

if ( abs(f_func_args.field2[idx+dj+di]-f_func_args.field2[idx+dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_1=-1;
    };

if ( abs(f_func_args.field2[idx-dj+di]-f_func_args.field2[idx-dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_11=-1;
    };

    return 1.0/(12.0*dx)*(
        (f[idx+di+dj] - kk_1*f[idx+dj-di])
        +4*(f[idx+di] - kk*f[idx-di])
        +(f[idx+di-dj] - kk_11*f[idx-di-dj]));
};

__device__ FFuncType d1yCO2RLI2D_dev=d1yCO2RLI2D;

// ---------------------------------------------------------------


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double laplaceCO2RLI2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;

      double kk = 1;
    double kk_1 = 1;
    double kk_11 = 1;

if ( abs(f_func_args.field2[idx+di]-f_func_args.field2[idx-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk=-1;
    };

if ( abs(f_func_args.field2[idx+dj+di]-f_func_args.field2[idx+dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_1=-1;
    };

if ( abs(f_func_args.field2[idx-dj+di]-f_func_args.field2[idx-dj-di]) > 2.9) {
        /* kk=f_func_args.field4[idx]/f_func_args.field3[idx]; */
        // kk=f_func_args.field4[idx-di]/f_func_args.field3[idx-di];
      kk_11=-1;
    };

    return 1.0/(dx*dy)*( -10.0/3.0*f[idx]
    +2.0/3.0*( f[idx+di] + kk*f[idx-di]
    +f[idx+dj] + f[idx-dj] )
    +1.0/6.0*( f[idx+di+dj] + kk_1*f[idx-di+dj]
    + f[idx+di-dj] + kk_11*f[idx-di-dj] ));
};

__device__ FFuncType laplaceCO2RLI2D_dev=laplaceCO2RLI2D;




void addUserDefinedFuncs () {
    f_func_map_all[{"gausF",""}]=gausF;
    f_func_map_all_dev[{"gausF",""}] = getFFuncDevPtr(&gausF_dev);
    f_func_map_all[{"atan2F",""}]=atan2F;
    f_func_map_all[{"getAlpha",""}]=getAlpha;
    f_func_map_all[{"change",""}]=changeF;
    f_func_map_all_dev[{"atan2F",""}]=getFFuncDevPtr(&atan2F_dev);
    f_func_map_all_dev[{"change",""}]=getFFuncDevPtr(&changeF_dev);
    f_func_map_all[{"getPressure",""}]=getPressure;
    f_func_map_all[{"oneOverPhi",""}]=oneOverPhi;
    f_func_map_all_dev[{"getPressure",""}]=getFFuncDevPtr(&getPressure_dev);
    f_func_map_all_dev[{"oneOverPhi",""}]=getFFuncDevPtr(&oneOverPhi_dev);
    f_func_map_all_dev[{"d1yRL",""}]=getFFuncDevPtr(&d1yRL_dev);
    f_func_map_all_dev[{"d1xRLI",""}]=getFFuncDevPtr(&d1xCO2RLI2D_dev);
    f_func_map_all_dev[{"d1yRLI",""}]=getFFuncDevPtr(&d1yCO2RLI2D_dev);
    f_func_map_all_dev[{"laplaceRL",""}]=getFFuncDevPtr(&laplaceCO2RLI2D_dev);
};


#endif
