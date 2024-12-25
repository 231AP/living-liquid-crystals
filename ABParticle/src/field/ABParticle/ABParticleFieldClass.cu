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

void ExpoConf(const std::string& str_t) {
    std::ofstream ConfFile;
    int PrecData = 8;
    std::string ConfFileName = "../data/test54/conf_" + str_t + ".dat";
    ConfFile.open(ConfFileName.c_str());

    if (!ConfFile.is_open()) {
        std::cerr << "Can't Open File: " << ConfFileName << std::endl;
        return;
    }
    for (int idx = 0; idx < PM.particleNum; idx++) {
        ConfFile << std::fixed << std::setprecision(PrecData)
            << pt.x[idx] << ' '
            << pt.y[idx] << ' '
            << pt.px[idx] << ' '
            << pt.py[idx] << ' ';
            // << pt.cellPxx[idx] << ' '
            // << pt.cellPxy[idx] << ' ';
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
    delete[] pt.particleAroundId;
    delete[] pt.particleAroundFlagX;
    delete[] pt.particleAroundFlagY;
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
    delete[] pt.particleAroundId1;
    delete[] pt.particleAroundFlagX1;
    delete[] pt.particleAroundFlagY1;
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
    cudaFree(PT.particleAroundId);
    cudaFree(PT.particleAroundFlagX);
    cudaFree(PT.particleAroundFlagY);
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
    cudaFree(PT.particleAroundId1);
    cudaFree(PT.particleAroundFlagX1);
    cudaFree(PT.particleAroundFlagY1);
    cudaFree(PT.offsetsAL1);
    cudaFree(PT.offsetsNL1);
    cudaFree(PT.NeighborList1);
    cudaFree(PT.NeighborListFlagX1);
    cudaFree(PT.NeighborListFlagY1);

    cudaFree(PT.x0ToUpdateHybridList1);
    cudaFree(PT.y0ToUpdateHybridList1);
}




// =============================================================
void getInput() {
    std::ifstream InputFile("input.dat");

    if (!InputFile.is_open()) {
        std::cerr << "Error opening input file!" << std::endl;
        return; // 退出函数
    }

    std::string line;
    int lineCount = 0;

    while (std::getline(InputFile, line)) {
        // 检查是否为注释行
        if (line.empty() || line.find('#') != std::string::npos) {
            continue; // 
        }

        std::istringstream iss(line);
        switch (lineCount) {
        case 0: iss >> PM.boxX; break;
        case 1: iss >> PM.boxY; break;
        case 2: iss >> PM.cellSizeX; break;
        case 3: iss >> PM.cellSizeY; break;
        case 4: iss >> PM.cellNumX; break;
        case 5: iss >> PM.cellNumY; break;
        case 6: iss >> PM.rho; break;
        case 7: iss >> PM.maxParticlePerCell; break;
        case 8: iss >> PM.rd; break;
        case 9: iss >> PM.mask0; break;
        case 10: iss >> PM.mask1; break;
        case 11: iss >> PM.miniInstanceBetweenParticle; break;
        case 12: iss >> PM.r0; break;
        case 13: iss >> PM.epsilon; break;
        case 14: iss >> PM.kBT; break;
        case 15: iss >> PM.gammaValue; break;
        case 16: iss >> PM.rOutUpdateList; break;
        case 17: iss >> PM.particleNum; break;
        case 18: iss >> PM.tStart; break;
        case 19: iss >> PM.tStop; break;
        case 20: iss >> PM.tStep; break;
        case 21: iss >> PM.tExpo; break;
        case 22: iss >> PM.V0; break;
        case 23: iss >> PM.ABParticle; break;


        case 24: iss >> PM.cellSizeX1; break;
        case 25: iss >> PM.cellSizeY1; break;
        case 26: iss >> PM.cellNumX1; break;
        case 27: iss >> PM.cellNumY1; break;
        case 28: iss >> PM.maxParticlePerCell1; break;
        case 29: iss >> PM.rOutUpdateList1; break;
        default: break; // 
        }
        lineCount++;
    }

    InputFile.close();
}







//=================================================================================





void MemAlloc() {
    // Allocate particle mem in host memory.
    pt.x = new real[PM.particleNum];
    pt.y = new real[PM.particleNum];


    pt.px = new real[PM.particleNum];
    pt.py = new real[PM.particleNum];



    pt.cellList = new int[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
    pt.cellOffsetsCL = new int[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];//删了*PM.maxParticlePerCell会有问题

    pt.cellList1 = new int[PM.cellNumX1 * PM.cellNumY1 * PM.maxParticlePerCell1];
    pt.cellOffsetsCL1 = new int[PM.cellNumX1 * PM.cellNumY1 * PM.maxParticlePerCell1];//删了*PM.maxParticlePerCell会有问题

    pt.cellPx= new real[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
    pt.cellPy = new real[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
    pt.cellPxx = new real[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
    pt.cellPxy = new real[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];


    pt.particleAroundId = new int[9 * PM.particleNum * PM.maxParticlePerCell];
    pt.particleAroundFlagX = new int[9 * PM.particleNum * PM.maxParticlePerCell];
    pt.particleAroundFlagY = new int[9 * PM.particleNum * PM.maxParticlePerCell];
    pt.offsetsNL = new int[PM.particleNum];
    pt.NeighborList = new int[PM.particleNum * PM.maxParticlePerCell];
    pt.NeighborListFlagX = new int[PM.particleNum];
    pt.NeighborListFlagY = new int[PM.particleNum];
    pt.fx = new real[PM.particleNum];
    pt.fy = new real[PM.particleNum];
    pt.x0ToUpdateHybridList = new real[PM.particleNum];
    pt.y0ToUpdateHybridList = new real[PM.particleNum];


    pt.particleAroundId1 = new int[9 * PM.particleNum * PM.maxParticlePerCell1];
    pt.particleAroundFlagX1 = new int[9 * PM.particleNum * PM.maxParticlePerCell1];
    pt.particleAroundFlagY1 = new int[9 * PM.particleNum * PM.maxParticlePerCell1];
    pt.offsetsNL1 = new int[PM.particleNum];
    pt.NeighborList1 = new int[PM.particleNum * PM.maxParticlePerCell1];
    pt.NeighborListFlagX1 = new int[PM.particleNum];
    pt.NeighborListFlagY1 = new int[PM.particleNum];

    pt.x0ToUpdateHybridList1 = new real[PM.particleNum];
    pt.y0ToUpdateHybridList1 = new real[PM.particleNum];
    pt.state = new curandState[PM.particleNum];


    // Allocate memory of fields in device.
    cudaMalloc((void**)&PT.x, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.y, PM.particleNum * sizeof(real));



    cudaMalloc((void**)&PT.px,PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.py,PM.particleNum * sizeof(real));



    cudaMalloc((void**)&PT.cellX, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.cellY, PM.particleNum * sizeof(int));
 cudaMalloc((void**)&PT.cellX1, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.cellY1, PM.particleNum * sizeof(int));

    cudaMalloc((void**)&PT.cellList, PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.cellOffsetsCL, PM.cellNumX * PM.cellNumY   * PM.maxParticlePerCell* sizeof(int));
cudaMalloc((void**)&PT.cellList1, PM.cellNumX1 * PM.cellNumY1 * PM.maxParticlePerCell1 * sizeof(int));
    cudaMalloc((void**)&PT.cellOffsetsCL1, PM.cellNumX1 * PM.cellNumY1   * PM.maxParticlePerCell1* sizeof(int));


    cudaMalloc((void**)&PT.cellPx, PM.cellNumX * PM.cellNumY  * PM.maxParticlePerCell* sizeof(real));
    cudaMalloc((void**)&PT.cellPy, PM.cellNumX * PM.cellNumY  * PM.maxParticlePerCell* sizeof(real));

    cudaMalloc((void**)&PT.cellPxx, PM.cellNumX * PM.cellNumY  * PM.maxParticlePerCell* sizeof(real));
    cudaMalloc((void**)&PT.cellPxy, PM.cellNumX * PM.cellNumY  * PM.maxParticlePerCell* sizeof(real));


    cudaMalloc((void**)&PT.particleAroundId, 9 * PM.maxParticlePerCell * PM.particleNum * sizeof(int));  //这里以后可以把9改成作用力范围
    cudaMalloc((void**)&PT.particleAroundFlagX, 9 * PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.particleAroundFlagY, 9 * PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.offsetsAL, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.offsetsNL, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.NeighborList, PM.particleNum * PM.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagX, PM.particleNum * PM.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagY, PM.maxParticlePerCell * PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.fx, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.fy, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.x0ToUpdateHybridList, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.y0ToUpdateHybridList, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.state, PM.particleNum * sizeof(curandState));


     cudaMalloc((void**)&PT.particleAroundId1, 9 * PM.maxParticlePerCell1 * PM.particleNum * sizeof(int));  //这里以后可以把9改成作用力范围
    cudaMalloc((void**)&PT.particleAroundFlagX1, 9 * PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.particleAroundFlagY1, 9 * PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.offsetsAL1, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.offsetsNL1, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.NeighborList1, PM.particleNum * PM.maxParticlePerCell1 * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagX1, PM.particleNum * PM.maxParticlePerCell1 * sizeof(int));
    cudaMalloc((void**)&PT.NeighborListFlagY1, PM.maxParticlePerCell1 * PM.particleNum * sizeof(int));

    cudaMalloc((void**)&PT.x0ToUpdateHybridList1, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.y0ToUpdateHybridList1, PM.particleNum * sizeof(real));
}

//=================================================================================







void printInput() {
    std::cout << "Box X: " << PM.boxX << std::endl;
    std::cout << "Box Y: " << PM.boxY << std::endl;
    std::cout << "Cell size X: " << PM.cellSizeX << std::endl;
    std::cout << "Cell size Y: " << PM.cellSizeY << std::endl;
    std::cout << "Cell num X: " << PM.cellNumX << std::endl;
    std::cout << "Cell num Y: " << PM.cellNumY << std::endl;
    std::cout << "Density: " << PM.rho << std::endl;
    std::cout << "Max particle per cell: " << PM.maxParticlePerCell << std::endl;
    std::cout << "Deadline distance: " << PM.rd << std::endl;
    std::cout << "Mask 0: " << PM.mask0 << std::endl;
    std::cout << "Mask 1: " << PM.mask1 << std::endl;
    std::cout << "Mini instance between particle: " << PM.miniInstanceBetweenParticle << std::endl;
    std::cout << "Equilibrium position: " << PM.r0 << std::endl;
    std::cout << "Epsilon: " << PM.epsilon << std::endl;
    std::cout << "kBT: " << PM.kBT << std::endl;
    std::cout << "Gamma value: " << PM.gammaValue << std::endl;
    std::cout << "Update list distance threshold: " << PM.rOutUpdateList << std::endl;
    std::cout << "Particle num: " << PM.particleNum << std::endl;
    std::cout << "Start time: " << PM.tStart << std::endl;
    std::cout << "Stop time: " << PM.tStop << std::endl;
    std::cout << "Time step: " << PM.tStep << std::endl;
    std::cout << "TExpo: " << PM.tExpo << std::endl;
    std::cout << "V0: " << PM.V0 << std::endl;
    printf("1");
    std::cout << "ABParticle: " << PM.ABParticle << std::endl;
    printf("2");
     std::cout << "ABParticle: " << PM.ABParticle << std::endl;
         std::cout << "Cell size X1: " << PM.cellSizeX1 << std::endl;
    std::cout << "Cell size Y1: " << PM.cellSizeY1 << std::endl;
    std::cout << "Cell num X1: " << PM.cellNumX1 << std::endl;
    std::cout << "Cell num Y1: " << PM.cellNumY1 << std::endl;
    std::cout << "Density: " << PM.rho << std::endl;
    std::cout << "Max particle per cell1: " << PM.maxParticlePerCell1 << std::endl;
};

//=============================================================================

void Init_Coords(int flag, Particle pt, Parameter PM) {
   
    if (flag == 0) {

    //    printf("12312313");
  
        int N = PM.particleNum;
        real rho = PM.rho;
        real L = sqrt(N / rho);
   
        real xBox = L;
        real yBox = L;
        PM.boxX = xBox;
        PM.boxY = yBox;
        int initUcell = sqrt(N); 
        real d_lattice = L / sqrt(N); 

        int n, nx, ny;
        n = 0;
        for (ny = 0;ny < initUcell; ny++) {
            for (nx = 0;nx < initUcell; nx++) {
                pt.x[n] = nx * d_lattice;
                pt.y[n] = ny * d_lattice;
                n++;
            }
        }
    }

    else if (flag == 1) {
        // printf("12312313");
        std::default_random_engine e;
        std::uniform_real_distribution<double> u(0.0, 1.0);
        // std::uniform_real_distribution<double> u1(0.0, 1.0);
        e.seed(time(0));
        for (int n = 0; n < PM.particleNum; n++) {
            int flag = 0;
            pt.x[n] = u(e) * PM.boxX;
            pt.y[n] = u(e) * PM.boxY;
            pt.px[n] = cos(3140*u(e));
            pt.py[n] = sin(3140*u(e));
            while (1) {
                for (int m = 0; m < n; m++) {

                    float dx = fmod((pt.x[n] - pt.x[m] + PM.boxX), PM.boxX);
                    float dy = fmod((pt.y[n] - pt.y[m] + PM.boxY), PM.boxY);


                    if (dx > PM.boxX / 2) dx -= PM.boxX;
                    if (dy > PM.boxY / 2) dy -= PM.boxY;

   
                    float dist2 = dx * dx + dy * dy;
                    if (dist2 < PM.r0 * PM.r0) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 1) {
                    pt.x[n] = u(e) * PM.boxX;
                    pt.y[n] = u(e) * PM.boxY;
                    flag = 0;
                }
                else {
                    break; 
                }
            }
        }
    }
    else if (flag == 2) {
 
        int n = 0;
        int Ln = sqrt(PM.particleNum);

        real dx = PM.boxX / (Ln - 1);
        real dy = PM.boxY / (Ln - 1);

        for (int i = 0; i < Ln; i++) {
            for (int j = 0; j < Ln; j++) {
                real x = j * dx; 
                real y = i * dy; 
                pt.x[n] = x;
                pt.y[n] = y;
                n++;
            }
        }
    }
}

//==============================================

void parameterInit() {
    PM.lenBit = 0;
    real boxToIntreal = PM.boxX / PM.miniInstanceBetweenParticle;
    while (++PM.lenBit) {//ignore boxX very small
        if (boxToIntreal < (1 << PM.lenBit)) break;
    }
    PM.lenBit++;
    PM.mask0 = (1 << PM.lenBit) + (1 << 2 * PM.lenBit + 1);
    int bitRd = ceil(log(PM.rd / PM.miniInstanceBetweenParticle) / log(2.0f));
    PM.mask1 = (((1 << (PM.lenBit - bitRd)) - 1) << bitRd) + (((1 << (PM.lenBit - bitRd)) - 1) << (bitRd + PM.lenBit + 1));
}

//上传=============================================================================================
void HostUpdataToDevice() {
    cudaMemcpy(PT.x, pt.x, PM.particleNum * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.y, pt.y, PM.particleNum * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.px, pt.px, PM.particleNum * sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(PT.py, pt.py, PM.particleNum * sizeof(real), cudaMemcpyHostToDevice);
}

//下载=============================================================================================
void DeviceUpdataToHost() {
    cudaMemcpy(pt.x, PT.x, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.y, PT.y, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.px, PT.px, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.py, PT.py, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
}


//===========================================================================

void InitOffset() {
    cudaMemset(PT.cellOffsetsCL, 0, sizeof(int) * PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell);
    cudaMemset(PT.offsetsNL, 0, sizeof(int) * PM.particleNum);
    cudaMemset(PT.offsetsAL, 0, sizeof(int) * PM.particleNum);
    cudaMemset(PT.cellOffsetsCL1, 0, sizeof(int) * PM.cellNumX1 * PM.cellNumY1 * PM.maxParticlePerCell1);
    cudaMemset(PT.offsetsNL1, 0, sizeof(int) * PM.particleNum);
    cudaMemset(PT.offsetsAL1, 0, sizeof(int) * PM.particleNum);
}

void listUpdate(Particle PT,Parameter PM) {
    // printf("444");
    InitOffset();
    getCellList << <PM.blockNum, PM.threadNum >> > (PT, PM);
    
    cudaDeviceSynchronize();
    getAroundCellParticleId << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
    saveXY0ToUpdateHybridList << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
}



//===================================================================
void initBlockAndThreadNum() {
    PM.threadNum = 256;
    PM.blockNum = (PM.particleNum + PM.threadNum - 1) / PM.threadNum;
    printf("blockNum:%d,threadNum:%d\n", PM.blockNum, PM.threadNum);
}

//==========================================================================================================
// void showProgress(real tNow, real tStart, real tStop, clock_t clockNow, clock_t clockStart) {
//     real progress = (tNow - tStart) / (tStop - tStart);
//     real tUsed = double(clockNow - clockStart) / CLOCKS_PER_SEC;
//     real tUsePrediction = (tStop - tNow) * tUsed / (tNow - tStart);
//     printf("%.8f,%.8f  ", pt.x[0], pt.y[0]);
//     printf("  Progress:%.4f%%,%Prediction:%.1f\r", progress*100, tUsePrediction);
//     fflush(stdout);
// }

//=======================================================

//=========================================================================
void ABParticleField::forceAndPositionUpdate(Particle PT, Parameter PM,int i_field) {
 int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    getForce << <PM.blockNum, PM.threadNum>> > (PT, PM,(*ptr_vx).f[i_field],(*ptr_vy).f[i_field],(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],(*ptr_Qxx).f[i_field],(*ptr_Qxy).f[i_field], Nx, Ny, Nbx, Nby);
    cudaDeviceSynchronize();
    updatePosition << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
}

//=========================================================================

// void ABParticleField::CellPUpdate(Particle PT,Parameter PM,int i_field){
    
//     int Nx=gridNumber().x;
//     int Ny=gridNumber().y;
//     int Nbx=gridNumberBoun().x;
//     int Nby=gridNumberBoun().y;
//     double dx=gridSize().x;
//     double dy=gridSize().y;
//     // printf("%d",PM.threadNum);
//     // cout <<"ABParticle detected3."<<endl;
    
//     updateConcentration << <Ny,Nx>> > (PT, PM,(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],f[i_field], Nx, Ny, Nbx, Nby);
// // cout <<"ABParticle detected4."<<endl;
// };




//=========================================================================
void ABParticleField::iterate(Particle PT,Parameter PM,int i_field) {
    forceAndPositionUpdate(PT,PM,i_field);



    checkUpdate << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaMemcpyFromSymbol(&updateListFlagHost, updateListFlag, sizeof(int));



    if (updateListFlagHost){
        // printf("update");
        listUpdate(PT, PM);
        updateListFlagHost = 0;
        // int updatefiled =1;
        cudaMemcpyToSymbol(updateListFlag, &updateListFlagHost, sizeof(int));
        // printf("1");
    };
    // printf("1");
    
    // printf("2");
};
// void ABParticleField::getPxxPxy(Particle PT, Parameter PM) {
//     int Nx = 

// }
void ABParticleField::ParticleToField(int i_field) {
// void ABParticleField::getConcentration(int i_field) {
    // Get velocity from vorticity field
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;

    for (int t1 = 0; t1 <PM.tExpo/PM.tStep; t1++){
        iterate(PT,PM,i_field);
        // printf("1");
    };
    getABParticlePxPy<<<Ny,Nx>>>(PT,PM,Nx,Ny);
    updateConcentration << <Ny,Nx>> > (PT, PM,(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],(*ptr_C1).f[i_field],Nx, Ny, Nbx, Nby);
    // printf("1");
     for (int t2 = 0; t2 <30; t2++){
        (*ptr_C1).applyBounCondPeriGPU((*ptr_C1).f[i_field]);
        smoothConcentration << <Ny,Nx>> > (PM,(*ptr_C1).f[i_field],(*ptr_Pxx).f[i_field],(*ptr_Pxy).f[i_field],Nx,Ny,Nbx,Nby);
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
