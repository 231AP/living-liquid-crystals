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
    //设置输出精度
    int PrecData = 8;

    // 文件名
    std::string ConfFileName = "conf_" + str_t + ".dat";
    ConfFile.open(ConfFileName.c_str());

    if (!ConfFile.is_open()) {
        std::cerr << "无法打开文件: " << ConfFileName << std::endl;
        return;
    }
    for (int idx = 0; idx < PM.particleNum; idx++) {
        // 使用固定格式和精度输出数据
        ConfFile << std::fixed << std::setprecision(PrecData)
            << pt.x[idx] << ' '
            << pt.y[idx];
        ConfFile << std::endl; // 换行
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

    // Free device memory
    cudaFree(PT.x);
    cudaFree(PT.y);

    cudaFree(PT.px);
    cudaFree(PT.py);
   


    cudaFree(PT.cellX);
    cudaFree(PT.cellY);

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
            continue; // 跳过空行和注释行
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

        default: break; // 超过预期行数时不处理
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
    pt.cellOffsetsCL = new int[PM.cellNumX * PM.cellNumY *2];
    pt.cellPx= new real[PM.cellNumX * PM.cellNumY * 2];
    pt.cellPy = new real[PM.cellNumX * PM.cellNumY * 2];


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
    pt.state = new curandState[PM.particleNum];


    // Allocate memory of fields in device.
    cudaMalloc((void**)&PT.x, PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.y, PM.particleNum * sizeof(real));



    cudaMalloc((void**)&PT.px,PM.particleNum * sizeof(real));
    cudaMalloc((void**)&PT.py,PM.particleNum * sizeof(real));



    cudaMalloc((void**)&PT.cellX, PM.particleNum * sizeof(int));
    cudaMalloc((void**)&PT.cellY, PM.particleNum * sizeof(int));


    cudaMalloc((void**)&PT.cellList, PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell * sizeof(int));
    cudaMalloc((void**)&PT.cellOffsetsCL, PM.cellNumX * PM.cellNumY   *2* sizeof(int));


    cudaMalloc((void**)&PT.cellPx, PM.cellNumX * PM.cellNumY  * 2* sizeof(real));
    cudaMalloc((void**)&PT.cellPy, PM.cellNumX * PM.cellNumY  * 2* sizeof(real));


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
    std::cout << "ABParticle: " << PM.ABParticle << std::endl;
}

//=============================================================================

void Init_Coords(int flag, Particle pt, Parameter PM) {
    /*
    flag代表系统的初始化方式，flag=0代表均匀分布，flag=1代表随机分布
    当按照均匀分布时，需给定粒子密度，会同时按照初始粒子数目,初始系统的周期盒大小；
    当按照随机分布时，需给定粒子数目，随机生成粒子坐标
    */

    if (flag == 0) {
        //初始周期盒长度
        int N = PM.particleNum;
        real rho = PM.rho;
        real L = sqrt(N / rho);
        //考虑正方形盒子
        real xBox = L;
        real yBox = L;
        PM.boxX = xBox;
        PM.boxY = yBox;
        int initUcell = sqrt(N); //初始x,y,方向粒子数目
        real d_lattice = L / sqrt(N); //晶格间距
        //均匀分布 系统以原点为中心
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
    //随机分布 均匀分布的随机数生成器
    else if (flag == 1) {
        std::default_random_engine e;
        std::uniform_real_distribution<double> u(0.0, 1.0);
        e.seed(time(0));
        for (int n = 0; n < PM.particleNum; n++) {
            int flag = 0;
            pt.x[n] = u(e) * PM.boxX;
            pt.y[n] = u(e) * PM.boxY;

            pt.px[n] = u(e);
            pt.py[n] = u(e);



            while (1) {
                for (int m = 0; m < n; m++) {
                    // 计算两个粒子之间的距离，考虑周期性边界条件
                    float dx = fmod((pt.x[n] - pt.x[m] + PM.boxX), PM.boxX);
                    float dy = fmod((pt.y[n] - pt.y[m] + PM.boxY), PM.boxY);





                    // 若计算结果为负数，则调整到正值
                    if (dx > PM.boxX / 2) dx -= PM.boxX;
                    if (dy > PM.boxY / 2) dy -= PM.boxY;

                    // 计算距离的平方
                    float dist2 = dx * dx + dy * dy;

                    // 如果距离小于某个阈值（如 r0/2），则重新生成位置
                    if (dist2 < PM.r0 * PM.r0) {
                        flag = 1;
                        break;
                    }
                }

                // 如果找到距离太近的粒子，重新生成位置
                if (flag == 1) {
                    pt.x[n] = u(e) * PM.boxX;
                    pt.y[n] = u(e) * PM.boxY;
                    flag = 0;
                }
                else {
                    break;  // 如果所有的粒子距离都合适，退出循环
                }
            }

            //cout << u(e)<<"," << PM.boxX <<"," << pt.x[n] << endl;
        }
    }
    else if (flag == 2) {
        //计算粒子数
        int n = 0;
        int Ln = sqrt(PM.particleNum);
        //计算间距
        real dx = PM.boxX / (Ln - 1);
        real dy = PM.boxY / (Ln - 1);
        // 生成二维晶格的格点
        for (int i = 0; i < Ln; i++) {
            for (int j = 0; j < Ln; j++) {
                real x = j * dx; // 计算x坐标
                real y = i * dy; // 计算y坐标
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
}

//下载=============================================================================================
void DeviceUpdataToHost() {
    cudaMemcpy(pt.x, PT.x, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(pt.y, PT.y, PM.particleNum * sizeof(real), cudaMemcpyDeviceToHost);
}


//===========================================================================

void InitOffset() {
    cudaMemset(PT.cellOffsetsCL, 0, sizeof(int) * PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell);
    cudaMemset(PT.offsetsNL, 0, sizeof(int) * PM.particleNum);
    cudaMemset(PT.offsetsAL, 0, sizeof(int) * PM.particleNum);
}

void listUpdate(Particle PT,Parameter PM) {
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
void ABParticleField::iterate(Particle PT,Parameter PM,int i_field) {
    forceAndPositionUpdate(PT,PM,i_field);
    checkUpdate << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaMemcpyFromSymbol(&updateListFlagHost, updateListFlag, sizeof(int));
    if (updateListFlagHost){
        listUpdate(PT, PM);
        updateListFlagHost = 0;
        cudaMemcpyToSymbol(updateListFlag, &updateListFlagHost, sizeof(int));
        
    };
    
    // CellPUpdata(PT,PM);

};

void ABParticleField::getConcentration(int i_field) {
    // Get velocity from vorticity field
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    

    if(i_field==0){
        iterate(PT,PM,i_field);
    // getLivingLCPxPyThetaGPU<<<Ny,Nx>>>((*ptr_flip).f[i_field],(*ptr_Pxx).f[i_field], f[i_field], (*ptr_px).f[i_field], (*ptr_py).f[i_field], (*ptr_theta).f[i_field], (*ptr_theta_old).f[i_field], Nx, Ny, Nbx, Nby);
    // getLivingLCFlipGPU<<<Ny,Nx>>>((*ptr_Omega).f[i_field],(*ptr_cplus).f[i_field],(*ptr_cminus).f[i_field],(*ptr_theta_old).f[i_field], (*ptr_theta).f[i_field],(*ptr_flip).f[i_field], Nx, Ny, Nbx, Nby);
    
    
    
    (*ptr_Pxx).applyBounCondPeriGPU((*ptr_Pxx).f[i_field]);
    // (*ptr_Concentration).applyBounCondPeriGPU((*ptr_Concentration).f[i_field]);
    (*ptr_Pxy).applyBounCondPeriGPU((*ptr_Pxy).f[i_field]);
    (*ptr_Qxx).applyBounCondPeriGPU((*ptr_Qxx).f[i_field]);
    (*ptr_Qxy).applyBounCondPeriGPU((*ptr_Qxy).f[i_field]);
     (*ptr_vx).applyBounCondPeriGPU((*ptr_vx).f[i_field]);
    (*ptr_vy).applyBounCondPeriGPU((*ptr_vy).f[i_field]);
    };
    
};


// =============================================================

#endif
