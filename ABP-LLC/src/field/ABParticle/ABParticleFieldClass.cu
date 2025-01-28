#ifndef ABPARTICLEFIELDCLASS_CU
#define ABPARTICLEFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "ABParticleFieldClass.h"
#include "ABParticleClassGPU.cu"


using namespace std;

// =============================================================
// Constructors
// -------------------------------------------------------------
ABParticleField::ABParticleField (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=0;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.location = "both";
    traits_host.expo_data = "on";    
    initPolarField ();
};


// -------------------------------------------------------------
ABParticleField::ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initPolarField ();
};


// -------------------------------------------------------------
ABParticleField::ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initPolarField ();
};


// -------------------------------------------------------------
ABParticleField::ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initPolarField ();
};


// -------------------------------------------------------------
void ABParticleField::initPolarField() {
    
};


// -------------------------------------------------------------
// =============================================================

// void ExpoConf(const std::string& str_t,LLCParams &params) {
//     std::ofstream ConfFile;
//     int PrecData = 8;
//     std::string ConfFileName = params.datadir;
//     ConfFile.open(ConfFileName.c_str());

//     if (!ConfFile.is_open()) {
//         std::cerr << "Can't Open File: " << ConfFileName << std::endl;
//         return;
//     }
//     for (int idx = 0; idx < params.numParticles; idx++) {
//         ConfFile << std::fixed << std::setprecision(PrecData)
//             << pt.x[idx] << ' '
//             << pt.y[idx] << ' '
//             << pt.px[idx] << ' '
//             << pt.py[idx] << ' ';
//             // << pt.cellPxx[idx] << ' '
//             // << pt.cellPxy[idx] << ' ';
//         ConfFile << std::endl; 
//     }

//     ConfFile.close();
// }

void ExpoConf(const std::string& str_t, LLCParams &params) {
    std::ofstream ConfFile;
    int PrecData = 8; // Precision for data
    // std::cout << params.datadir <<std::endl;
    std::string ConfFileName = params.datadir+"conf_" + str_t + ".dat"; // Assuming this holds the path
    ConfFile.open(ConfFileName.c_str());

    if (!ConfFile.is_open()) {
        std::cerr << "Failed to open file: " << ConfFileName
                  << " (Check permissions or directory existence)" << std::endl;
        return;
    }

    // Assuming pt is part of params, and you want to write data for each particle
    for (int idx = 0; idx < params.numParticles; idx++) {
        ConfFile << std::fixed << std::setprecision(PrecData)
                 << pt.x[idx] << ' '
                 << pt.y[idx] << ' '
                 << pt.px[idx] << ' '
                 << pt.py[idx] << ' ';
                 // Uncomment and add more parameters as needed
                 // << params.pt.cellPxx[idx] << ' '
                 // << params.pt.cellPxy[idx] << ' ';

        ConfFile << std::endl;
    }

    ConfFile.close();
}



// =============================================================


void MemFree() {
    // Free host memory
    delete[] pt.x;
    delete[] pt.y;

    delete[] pt.px;
    delete[] pt.py;

    delete[] pt.cellPx;
    delete[] pt.cellPy;

    

    delete[] pt.cellList;
    delete[] pt.cellOffsetsCL;
    
    delete[] pt.offsetsNL;
    delete[] pt.offsetsAL;
    delete[] pt.NeighborList;
    delete[] pt.NeighborListFlagX;
    delete[] pt.NeighborListFlagY;
    delete[] pt.fx;
    delete[] pt.fy;
    delete[] pt.x0ToUpdateHybridList;
    delete[] pt.y0ToUpdateHybridList;
    delete[] pt.state;


    delete[] pt.cellList1;
    delete[] pt.cellOffsetsCL1;
   
    delete[] pt.offsetsNL1;
    delete[] pt.offsetsAL1;
    delete[] pt.NeighborList1;
    delete[] pt.NeighborListFlagX1;
    delete[] pt.NeighborListFlagY1;

    delete[] pt.x0ToUpdateHybridList1;
    delete[] pt.y0ToUpdateHybridList1;





    // Free device memory
    cudaFree(PT.x);
    cudaFree(PT.y);

    cudaFree(PT.px);
    cudaFree(PT.py);
   


    cudaFree(PT.cellX);
    cudaFree(PT.cellY);
    cudaFree(PT.cellX1);
    cudaFree(PT.cellY1);

    cudaFree(PT.cellPx);
    cudaFree(PT.cellPy);




    cudaFree(PT.cellList);
    cudaFree(PT.cellOffsetsCL);
    
    cudaFree(PT.offsetsAL);
    cudaFree(PT.offsetsNL);
    cudaFree(PT.NeighborList);
    cudaFree(PT.NeighborListFlagX);
    cudaFree(PT.NeighborListFlagY);
    cudaFree(PT.fx);
    cudaFree(PT.fy);
    cudaFree(PT.x0ToUpdateHybridList);
    cudaFree(PT.y0ToUpdateHybridList);
    cudaFree(PT.state);





    cudaFree(PT.cellList1);
    cudaFree(PT.cellOffsetsCL1);
    
    cudaFree(PT.offsetsAL1);
    cudaFree(PT.offsetsNL1);
    cudaFree(PT.NeighborList1);
    cudaFree(PT.NeighborListFlagX1);
    cudaFree(PT.NeighborListFlagY1);

    cudaFree(PT.x0ToUpdateHybridList1);
    cudaFree(PT.y0ToUpdateHybridList1);
}




// =============================================================

//=================================================================================





void MemAlloc(LLCParams &params) {
    // Allocate particle mem in host memory.
    pt.x = new real[params.numParticles];
    pt.y = new real[params.numParticles];


    pt.px = new real[params.numParticles];
    pt.py = new real[params.numParticles];



    pt.cellList = new int[params.numGridX * params.numGridY * params.maxParticlePerGrid];
    pt.cellOffsetsCL = new int[params.numGridX * params.numGridY * params.maxParticlePerGrid];//删了*params.maxParticlePerGrid会有问题

    pt.cellList1 = new int[params.cellNumX * params.cellNumY * params.maxParticlePerCell];
    pt.cellOffsetsCL1 = new int[params.cellNumX * params.cellNumY * params.maxParticlePerCell];//删了*params.maxParticlePerGrid会有问题

    pt.cellPx= new real[params.numGridX * params.numGridY * params.maxParticlePerGrid];
    pt.cellPy = new real[params.numGridX * params.numGridY * params.maxParticlePerGrid];
    pt.cellPxx = new real[params.numGridX * params.numGridY * params.maxParticlePerGrid];
    pt.cellPxy = new real[params.numGridX * params.numGridY * params.maxParticlePerGrid];


  
    pt.offsetsNL = new int[params.numParticles];
    pt.NeighborList = new int[params.numParticles * params.maxParticlePerGrid];
    pt.NeighborListFlagX = new int[params.numParticles];
    pt.NeighborListFlagY = new int[params.numParticles];
    pt.fx = new real[params.numParticles];
    pt.fy = new real[params.numParticles];
    pt.x0ToUpdateHybridList = new real[params.numParticles];
    pt.y0ToUpdateHybridList = new real[params.numParticles];


   
    pt.offsetsNL1 = new int[params.numParticles];
    pt.NeighborList1 = new int[params.numParticles * params.maxParticlePerCell];
    pt.NeighborListFlagX1 = new int[params.numParticles];
    pt.NeighborListFlagY1 = new int[params.numParticles];

    pt.x0ToUpdateHybridList1 = new real[params.numParticles];
    pt.y0ToUpdateHybridList1 = new real[params.numParticles];
    pt.state = new curandState[params.numParticles];


    // Allocate memory of fields in device.
    cudaMalloc((void**)&PT.x, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.y, params.numParticles * sizeof(real));



    cudaMalloc((void**)&PT.px,params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.py,params.numParticles * sizeof(real));



    cudaMalloc((void**)&PT.cellX, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.cellY, params.numParticles * sizeof(int));
 cudaMalloc((void**)&PT.cellX1, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.cellY1, params.numParticles * sizeof(int));

    cudaMalloc((void**)&PT.cellList, params.numGridX * params.numGridY * params.maxParticlePerGrid * sizeof(int));
    cudaMalloc((void**)&PT.cellOffsetsCL, params.numGridX * params.numGridY   * params.maxParticlePerGrid* sizeof(int));
cudaMalloc((void**)&PT.cellList1, params.cellNumX * params.cellNumY * params.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.cellOffsetsCL1, params.cellNumX * params.cellNumY   * params.maxParticlePerCell* sizeof(int));


    cudaMalloc((void**)&PT.cellPx, params.numGridX * params.numGridY  * params.maxParticlePerGrid* sizeof(real));
    cudaMalloc((void**)&PT.cellPy, params.numGridX * params.numGridY  * params.maxParticlePerGrid* sizeof(real));

    cudaMalloc((void**)&PT.cellPxx, params.numGridX * params.numGridY  * params.maxParticlePerGrid* sizeof(real));
    cudaMalloc((void**)&PT.cellPxy, params.numGridX * params.numGridY  * params.maxParticlePerGrid* sizeof(real));


    
    cudaMalloc((void**)&PT.offsetsAL, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.offsetsNL, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.NeighborList, params.numParticles * params.maxParticlePerGrid * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagX, params.numParticles * params.maxParticlePerGrid * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagY, params.maxParticlePerGrid * params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.fx, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.fy, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.x0ToUpdateHybridList, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.y0ToUpdateHybridList, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.state, params.numParticles * sizeof(curandState));


    
    cudaMalloc((void**)&PT.offsetsAL1, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.offsetsNL1, params.numParticles * sizeof(int));
    cudaMalloc((void**)&PT.NeighborList1, params.numParticles * params.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagX1, params.numParticles * params.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagY1, params.maxParticlePerCell * params.numParticles * sizeof(int));

    cudaMalloc((void**)&PT.x0ToUpdateHybridList1, params.numParticles * sizeof(real));
    cudaMalloc((void**)&PT.y0ToUpdateHybridList1, params.numParticles * sizeof(real));
}

//===============================================================================

//=============================================================================



// //==============================================
void Init_Coords(int flag, Particle pt, LLCParams &params) {
    memset(pt.cellList, 0, params.numGridX * params.numGridY * params.maxParticlePerGrid * sizeof(int));
    memset(pt.cellOffsetsCL, 0, params.numGridX * params.numGridY * params.maxParticlePerGrid * sizeof(int));
    int N=params.numParticles;
    real xBox = params.Lx;
    real yBox = params.Ly;
    std::default_random_engine e;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    e.seed(time(0));
    real x0, y0, dx, dy,px0,py0;
    int cellX,cellY,cellX1,cellY1,cellAround;
    int xFlag,yFlag;
    int wrongFlag=0;
    for(int i=0;i<N;i++){
        while(1){
            x0=u(e)*xBox;
            y0=u(e)*yBox;
            px0 = cos(3140*u(e));
            py0 = sin(3140*u(e));
            //printf("id:%d,x:%f,y:%f\n",i,x0,y0);
            wrongFlag=0;
            cellX=std::floor(x0/params.dx);
            cellY=std::floor(y0/params.dy);
            for(int x=-1;x<=1;x++){
                for(int y=-1;y<=1;y++){
                    if(cellX+x==-1){
                        cellX1=params.numGridX-1;
                        xFlag=1;
                    }else if(cellX+x==params.numGridX){
                        cellX1=0;
                        xFlag=-1;
                    }else{
                        cellX1=cellX+x;
                        xFlag=0;
                    }
                    if(cellY+y==-1){
                        cellY1=params.numGridY-1;
                        yFlag=1;
                    }else if(cellY+y==params.numGridY){                    
                        cellY1=0;
                        yFlag=-1;
                    }else{
                        cellY1=cellY+y;
                        yFlag=0;
                    }
                    cellAround=cellX1+cellY1*params.numGridX;
                    for(int j=0;j<pt.cellOffsetsCL[cellAround];j++){
                        //printf("cell:%d,cellAround:%d,j:%d,x0:%f,y0:%f,x:%f,y:%f\n",cellX+cellY*params.numGridX,cellAround,j,x0,y0,pt.x[pt.cellList[cellAround*params.maxParticlePerGrid+j]],pt.y[pt.cellList[cellAround*params.maxParticlePerGrid+j]]);
                        dx=(x0-pt.x[pt.cellList[cellAround*params.maxParticlePerGrid+j]])+xFlag*params.Lx;
                        dy=(y0-pt.y[pt.cellList[cellAround*params.maxParticlePerGrid+j]])+yFlag*params.Ly;
                        if(dx*dx+dy*dy<params.r0*params.r0){
                            wrongFlag=1;
                            break;
                        }
                    }
                    if(wrongFlag==1){
                        break;
                    }
                }
                if(wrongFlag==1){ 
                    break;
                }
            }
            if(wrongFlag==0){
                break;
            }else continue;
        }        
        pt.x[i]=x0;
        pt.y[i]=y0;
        pt.px[i] = px0;
        pt.py[i] = py0;
        //printf("id:%d,x:%f,y:%f,cell:%d,%d,cellSize:%f,%f,xbox:%f,ybox:%f\n",\
        i,pt.x[i],pt.y[i],cellX,cellY,params.dx,params.dy,xBox,yBox);
        pt.cellList[(cellX+cellY*params.numGridX)*params.maxParticlePerGrid+pt.cellOffsetsCL[cellX+cellY*params.numGridX]]=i;
        pt.cellOffsetsCL[cellX+cellY*params.numGridX]++;
        //printf("id:%d,cell:%d,%d,cellList:%d,cellOffsetsCL:%d\n",i,cellX,cellY,pt.cellList[cellX+cellY*params.numGridX],pt.cellOffsetsCL[cellX+cellY*params.numGridX]);
    }
}



//上传=============================================================================================
void HostUpdataToDevice(LLCParams &params) {
    cudaMemcpy(PT.x, pt.x, params.numParticles * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.y, pt.y, params.numParticles * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.px, pt.px, params.numParticles * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.py, pt.py, params.numParticles * sizeof(real), cudaMemcpyHostToDevice);
}

//下载=============================================================================================
void DeviceUpdataToHost(LLCParams &params) {
    cudaMemcpy(pt.x, PT.x, params.numParticles * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.y, PT.y, params.numParticles * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.px, PT.px, params.numParticles * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.py, PT.py, params.numParticles * sizeof(real), cudaMemcpyDeviceToHost);
}


//===========================================================================

void InitOffset(LLCParams &params) {
    cudaMemset(PT.cellOffsetsCL, 0, sizeof(int) * params.numGridX * params.numGridY * params.maxParticlePerGrid);
    cudaMemset(PT.offsetsNL, 0, sizeof(int) * params.numParticles);
    cudaMemset(PT.offsetsAL, 0, sizeof(int) * params.numParticles);
    cudaMemset(PT.cellOffsetsCL1, 0, sizeof(int) * params.cellNumX * params.cellNumY * params.maxParticlePerCell);
    cudaMemset(PT.offsetsNL1, 0, sizeof(int) * params.numParticles);
    cudaMemset(PT.offsetsAL1, 0, sizeof(int) * params.numParticles);
}

void listUpdate(Particle PT,LLCParams &params) {
    // printf("444");
    InitOffset(params);
    // test11 <<< 20,10 >>>();
    getCellList << <params.blockNumParticles, params.threadNumParticles >> > (PT,params);
    // test11 <<< 20,10 >>>();
    cudaDeviceSynchronize();
    getAroundCellParticleId << <params.blockNumParticles, params.threadNumParticles >> > (PT, params);
        // test11 <<< 20,10 >>>();

    cudaDeviceSynchronize();
    saveXY0ToUpdateHybridList << <params.blockNumParticles, params.threadNumParticles >> > (PT, params);
    cudaDeviceSynchronize();
}



//=========================================================================
void ABParticleField::forceAndPositionUpdate(Particle PT, LLCParams &params,int i_field) {
 int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    // printf("2");
    getForce << <params.blockNumParticles, params.threadNumParticles>> > (PT, params,(*ptr_vx).f[i_field],(*ptr_vy).f[i_field],(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],(*ptr_Qxx).f[i_field],(*ptr_Qxy).f[i_field], Nx, Ny, Nbx, Nby);
    cudaDeviceSynchronize();
    updatePosition << <params.blockNumParticles, params.threadNumParticles >> > (PT, params);
    cudaDeviceSynchronize();
}



//=========================================================================
void ABParticleField::iterate(Particle PT,LLCParams &params,int i_field) {
    forceAndPositionUpdate(PT,params,i_field);
    checkUpdate << <params.blockNumParticles, params.threadNumParticles >> > (PT, params);
    cudaMemcpyFromSymbol(&updateListFlagHost, updateListFlag, sizeof(int));



    if (updateListFlagHost){
        listUpdate(PT, params);
        updateListFlagHost = 0;
        cudaMemcpyToSymbol(updateListFlag, &updateListFlagHost, sizeof(int));
    };

};

void ABParticleField::ParticleToField(int i_field,LLCParams &params) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;

    for (int t1 = 0; t1 <params.tStepField/params.tStepParticle; t1++){
        iterate(PT,params,i_field);
        // printf("1");
    };
    getABParticlePxPy<<<Ny,Nx>>>(PT,params,Nx,Ny);
    updateConcentration << <Ny,Nx>> > (PT, params,(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],(*ptr_C1).f[i_field],Nx, Ny, Nbx, Nby);
    // printf("1");
     for (int t2 = 0; t2 <30; t2++){
        (*ptr_C1).applyBounCondPeriGPU((*ptr_C1).f[i_field]);
        (*ptr_Pxx).applyBounCondPeriGPU((*ptr_Pxx).f[i_field]);
        (*ptr_Pxy).applyBounCondPeriGPU((*ptr_Pxy).f[i_field]);
        smoothConcentration << <Ny,Nx>> > (params,(*ptr_C1).f[i_field],(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],Nx,Ny,Nbx,Nby);
    };

    getConcentration << <Ny,Nx>> > (f[i_field],(*ptr_C1).f[i_field],Nx,Ny,Nbx,Nby);
    
    (*ptr_Pxx).applyBounCondPeriGPU((*ptr_Pxx).f[i_field]);
    // (*ptr_Concentration).applyBounCondPeriGPU((*ptr_Concentration).f[i_field]);
    (*ptr_Pxy).applyBounCondPeriGPU((*ptr_Pxy).f[i_field]);
    (*ptr_Qxx).applyBounCondPeriGPU((*ptr_Qxx).f[i_field]);
    (*ptr_Qxy).applyBounCondPeriGPU((*ptr_Qxy).f[i_field]);
    (*ptr_vx).applyBounCondPeriGPU((*ptr_vx).f[i_field]);
    (*ptr_vy).applyBounCondPeriGPU((*ptr_vy).f[i_field]);
    
    // };
    
};


// =============================================================

#endif
