#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "cuda.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
#include "userDefinedFunction.h"
#include "../src/evolver/evolverclass.cu"
// #include "LLCParams.cu"

using namespace std;


// ======================================================================
int main() {
    // LLCParams params;
    // Read the parameter file
// test11 <<< 20,10 >>>();
    readParams("./parameter.dat", params);
    string direExpo= params.datadir ;
    // string direExpo='data/'
    string device="gpu";
    string FDMScheme="CentralDifferenceO2Iso2D";
    string timeScheme="EulerForward";

    double dt=params.tStepField;
    double T=params.tStop;
    double dtExpo=params.tExpo;
    int    NGrid=params.numGridX;
    double L=params.Lx;//NGrid;
    // Model parameters
        
    double a =params.hLDGa;
    // double AA =  -1; 
    double AA =  params.alpha; 
    
    double xian = params.xiAn;
    // double xian = 0;

    double K = params.hLDGk;
    double Gamma = params.Gamma;
    double eta = params.eta;
    double h = params.h;
    double xi = params.xi;  
    double b = params.hLDGb;
    double convection = params.convection;

    
    double zta = params.zta;

    

    // Generating system and mesh
    System mySys(direExpo);
    Mesh mesh(2);
    mesh.setGridNumber(NGrid,NGrid,1);
    mesh.setBoxSize(L,L,1);
    mySys.mesh_ptr=&mesh;

    // Declare fields
    // Nematics
    Field Anxx(&mesh,"anchx",0,"sin");
    Field Anxy(&mesh,"anchy",0,"sin");
    // Field Qxx(&mesh, "Qxx",0,"sin");
    // Field Qxy(&mesh, "Qxy",0,"sin","periodic","on");
    Field S2(&mesh, "S2",1);
    Field Hxx(&mesh,"Hxx",1,"sin","periodic","off");
    Field Hxy(&mesh,"Hxy",1,"sin","periodic","off");
    Field trQH(&mesh,"trQH",1,"sin","periodic","off");
    Field trQW(&mesh,"trQW",1,"sin","periodic","off");

    // Velocity
    IncompFlow incompFlow(&mesh, "incompFlow",1);
    Field dxvx(&mesh,"dxvx",1,"sin","periodic","off");
    Field dxvy(&mesh,"dxvy",1,"sin","periodic","off");
    Field dyvx(&mesh,"dyvx",1,"sin","periodic","off");    
    // Symmetric and anti-symmetric parts of velocity gradient
    Field Axx(&mesh,"Axx",1,"sin","periodic","off");
    Field Axy(&mesh,"Axy",1,"sin","periodic","off");
    Field Omega(&mesh,"Omega",1,"sin","periodic","off");
    Field sigmaA(&mesh,"sigmaA",1,"sin","periodic","off");
    Field sigmaSxx(&mesh,"sigmaSxx",1,"sin","periodic","off");
    Field sigmaSxy(&mesh,"sigmaSxy",1,"sin","periodic","off");
    Field sigmaKxx(&mesh,"sigmaKxx",1,"sin","periodic","off");
    Field sigmaKxy(&mesh,"sigmaKxy",1,"sin","periodic","off");    
    Field sigmaKyy(&mesh,"sigmaKyy",1,"sin","periodic","off");

    // Bacteria
    // Field cplus(&mesh,"cplus",0);
    // Field cminus(&mesh,"cminus",0);
    // LivingLC livingLC(&mesh, "livingLC", 0);  
    ABParticle abParticle(&mesh,"abParticle", 0);     
    
    
    Field P2(&mesh,"P2",1);       
    // Field actx(&mesh,"actx",1);
    // Field acty(&mesh,"acty",1);
    // Field actF(&mesh,"actF",1);    
    
    // Set equations
    // Nematics

    Anxx.setRhsTerms({
        {0,{{&abParticle.Qxx}}}
    });

    Anxy.setRhsTerms({
        {0,{{&abParticle.Qxy}}}
    });

    abParticle.Qxx.setRhsTerms({
        {Gamma,{{&Hxx}}},
        {-1*convection,{{"d1x",&abParticle.Qxx},{&incompFlow.vx}}},
        {-1*convection,{{"d1y",&abParticle.Qxx},{&incompFlow.vy}}},
        {2,{{&abParticle.Qxy},{&Omega}}},
        {xi,{{&Axx}}},
        {2*xi,{{&Axx},{&abParticle.Qxx}}},
        {2*xi,{{&Axy},{&abParticle.Qxy}}},
        {-xi,{{&trQW}}},
        {-2*xi,{{&trQW},{&abParticle.Qxx}}},
        {xian,{{&Anxx}}},       
    });
    abParticle.Qxy.setRhsTerms({        
        {Gamma,{{&Hxy}}},
        {-1*convection,{{"d1x",&abParticle.Qxy},{&incompFlow.vx}}},
        {-1*convection,{{"d1y",&abParticle.Qxy},{&incompFlow.vy}}},
        {-2,{{&abParticle.Qxx},{&Omega}}},
        {xi,{{&Axy}}},
        {-2*xi,{{&abParticle.Qxy},{&trQW}}},
        {xian,{{&Anxy}}},
    });
    
    S2.setRhsTerms({
        {4,{{&abParticle.Qxx},{&abParticle.Qxx}}},
        {4,{{&abParticle.Qxy},{&abParticle.Qxy}}}
    });
    Hxx.setRhsTerms({
        {a,{{&abParticle.Qxx}}},{-b/2,{{&S2},{&abParticle.Qxx}}},
        {K,{{"laplace",&abParticle.Qxx}}}
    });
    Hxy.setRhsTerms({
        {a,{{&abParticle.Qxy}}},{-b/2,{{&S2},{&abParticle.Qxy}}},
        {K,{{"laplace",&abParticle.Qxy}}}
    });
    trQW.setRhsTerms({
        {1,{{&abParticle.Qxy},{&dxvy}}},
        {1,{{&abParticle.Qxy},{&dyvx}}},
        {2,{{&abParticle.Qxx},{&dxvx}}}
    });
    trQH.setRhsTerms({
        {2,{{&Hxx},{&abParticle.Qxx}}},
        {2,{{&Hxy},{&abParticle.Qxy}}},
    });    

    // Velocity
    incompFlow.omega.setRhsTerms({
        {1,{{"laplace",&sigmaA}}},
        {1,{{"d2x",&sigmaSxy}}},
        {-1,{{"d2y",&sigmaSxy}}},
        {-2,{{"d1x1y",&sigmaSxx}}},
        {1,{{"d2x",&sigmaKxy}}},
        {-1,{{"d2y",&sigmaKxy}}},
        {-1,{{"d1x1y",&sigmaKxx}}},
        {1,{{"d1x1y",&sigmaKyy}}}
    });
    // incompFlow.omega.setRhsTerms({

    //      {-1,{{"d1x",&incompFlow.omega},{&incompFlow.vx}}},
    //     {-1,{{"d1y",&incompFlow.omega},{&incompFlow.vy}}},
    //     {1,{{"laplace",&sigmaA}}},
    //     {1,{{"d2x",&sigmaSxy}}},
    //     {-1,{{"d2y",&sigmaSxy}}},
    //     {-2,{{"d1x1y",&sigmaSxx}}},
    //     {1,{{"d2x",&sigmaKxy}}},
    //     {-1,{{"d2y",&sigmaKxy}}},
    //     {-1,{{"d1x1y",&sigmaKxx}}},
    //     {1,{{"d1x1y",&sigmaKyy}}},
    //     // {-gamma,{{&incompFlow.omega}}}
    // });

    abParticle.vx.setRhsTerms({
        {1*convection,{{&incompFlow.vx}}}
    });
    abParticle.vy.setRhsTerms({
        {1*convection,{{&incompFlow.vy}}}
    });



    dxvy.setRhsTerms({ {1,{{"d1x",&incompFlow.vy}}} });
    dyvx.setRhsTerms({ {1,{{"d1y",&incompFlow.vx}}} });
    dxvx.setRhsTerms({ {1,{{"d1x",&incompFlow.vx}}} });
    Axx.setRhsTerms({ {1,{{&dxvx}}} });
    Axy.setRhsTerms({ {0.5,{{&dyvx}}}, {0.5,{{&dxvy}}} });
    Omega.setRhsTerms({ {-0.5,{{&dxvy}}}, {0.5,{{&dyvx}}} });    
    // Anti-symmetric part of viscous stress from nematics
    sigmaA.setRhsTerms({
        {2,{{&abParticle.Qxy},{&Hxx}}},
        {-2,{{&abParticle.Qxx},{&Hxy}}},
    });

    //     actx.setRhsTerms({
    //     {AA,{{&abParticle.Concentration},{&abParticle.Pxx}}},
    //     // {AA,{{&abParticle.Concentrationminus},{&abParticle.Pxx}}}
    // });
    // acty.setRhsTerms({
    //     {AA,{{&abParticle.Concentration},{&abParticle.Pxy}}},
    //     // {AA,{{&abParticle.Concentrationminus},{&abParticle.Pxy}}}
    // });
    // actF.setRhsTerms({
    //     {1,{{"d1x",&actx}}},
    //     {1,{{"d1y",&acty}}},
    // });
    // Symmetric part of viscous stress from nematic
    sigmaSxx.setRhsTerms({
        {-2*xi,{{&Hxy},{&abParticle.Qxy}}},
        {-2*xi,{{&Hxx},{&abParticle.Qxx}}},
        {-xi,{{&Hxx}}},
        {2*xi,{{&abParticle.Qxx},{&trQH}}},
        {xi,{{&trQH}}},

        // {AA,{{&abParticle.Qxx}}},
        {AA,{{&abParticle.Concentration},{&abParticle.Pxx}}},
        // {AA,{{&abParticle.Concentration},{&abParticle.Qxx}}},
      
    });
    sigmaSxy.setRhsTerms({
        {-xi,{{&Hxy}}},
        {2*xi,{{&abParticle.Qxy},{&trQH}}},    
        // {AA,{{&abParticle.Qxy}}},   
        {AA,{{&abParticle.Concentration},{&abParticle.Pxy}}},
        // {AA,{{&abParticle.Concentration},{&abParticle.Qxy}}},

    });
    // Elastic stress from nematics
    sigmaKxx.setRhsTerms({
        {-2*K,{{"d1x",&abParticle.Qxx},{"d1x",&abParticle.Qxx}}},
        {-2*K,{{"d1x",&abParticle.Qxy},{"d1x",&abParticle.Qxy}}}
    });
    sigmaKxy.setRhsTerms({
        {-2*K,{{"d1x",&abParticle.Qxx},{"d1y",&abParticle.Qxx}}},
        {-2*K,{{"d1x",&abParticle.Qxy},{"d1y",&abParticle.Qxy}}}
    });
    sigmaKyy.setRhsTerms({
        {-2*K,{{"d1y",&abParticle.Qxx},{"d1y",&abParticle.Qxx}}},
        {-2*K,{{"d1y",&abParticle.Qxy},{"d1y",&abParticle.Qxy}}}
    });
    // Bacteria

    abParticle.Concentration.setRhsTerms({
        {0,{{&abParticle.Concentration}}}
    });

    abParticle.Pxx.setRhsTerms({
        {0,{{&abParticle.Pxx}}}
    });

    abParticle.Pxy.setRhsTerms({
        {0,{{&abParticle.Pxy}}}
    });
   
    P2.setRhsTerms({
        {1,{{&abParticle.Pxx},{&abParticle.Pxx}}},
        {1,{{&abParticle.Pxy},{&abParticle.Pxy}}}
    });

    
    // Adding fields to systems
    // Priority 1
    mySys.addField(&Anxx);
    mySys.addField(&Anxy);
    mySys.addField(&S2);
    mySys.addField(&Hxx);
    mySys.addField(&Hxy);
    mySys.addField(&trQH);
    // mySys.addField(&actx);
    // mySys.addField(&acty);
    // mySys.addField(&actF);
    mySys.addField(&sigmaSxx);
    mySys.addField(&sigmaSxy);
    mySys.addField(&sigmaA);
    mySys.addField(&sigmaKxx);
    mySys.addField(&sigmaKxy);
    mySys.addField(&sigmaKyy);
    mySys.addIncompFlow(&incompFlow);
    mySys.addField(&dxvx);
    mySys.addField(&dxvy);
    mySys.addField(&dyvx);
    mySys.addField(&Axx);
    mySys.addField(&Axy);
    mySys.addField(&Omega);    
    mySys.addField(&trQW);        
    mySys.addField(&P2);
    // Priority 0
    // mySys.addField(&Qxx);
    // mySys.addField(&Qxy);
    // mySys.addLivingLC(&livingLC);
    mySys.addABParticle(&abParticle);
// test11 <<< 20,10 >>>();
    
    // Print system information.
    mySys.printSysInfo();

    // Set initial conditions

    Anxx.initFieldImport();
    Anxy.initFieldImport();

    abParticle.Qxx.initFieldConst(0,0.01);
    abParticle.Qxy.initFieldConst(0,0.01);
    // abParticle.Pxx.initFieldConst(0,0.01);
    // abParticle.Pxy.initFieldConst(0,0.01);
    abParticle.Qxx.initFieldImport();
    abParticle.Qxy.initFieldImport();


    incompFlow.setOmegaEq(2,zta,eta);


// test11 <<< 20,10 >>>();
    MemAlloc(params);
    Init_Coords(1, pt, params);
    InitOffset(params);
    // test11 <<< 20,10 >>>();
    initState << <params.blockNumParticles, params.threadNumParticles >> > (PT.state, params.seed, params.numParticles);
    cudaDeviceSynchronize();
    // test11 <<< 20,10 >>>();

    HostUpdataToDevice(params);
    // test11 <<< 20,10 >>>();
    listUpdate(PT, params);
    // test11 <<< 20,10 >>>();
    cudaDeviceSynchronize();
// test11 <<< 20,10 >>>();
    // Creating an evolver:   
    //    printf("1");
    // test11 <<< 20,10 >>>();
    Evolver evolver(&mySys,pt,PT,params,0,T,dt,dtExpo,device,timeScheme,FDMScheme);
    // test11 <<< 20,10 >>>();
    printParams(params);

// test11 <<< 20,10 >>>();
    // Running simulations
    evolver.run();
    
    MemFree();//

    // err = cudaGetLastError();
    // if (err != cudaSuccess) {
    //     printf("CUDA error: %s\n", cudaGetErrorString(err));
    // }

    return 0; // 
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};
