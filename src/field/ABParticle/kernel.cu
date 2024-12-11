#include <cufftXt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <time.h>
#include <math.h>
#include <random> // The header for the generators.
#include <iomanip>
#include "cuda_runtime.h"
#include <cuda.h>
#include <curand_kernel.h>

//if no run in a long time ,maybe beacause box is too small
//there exist a bug that position of one particle will change abruptly 
//change 387 PM.lenBit-- to lenBit++

//Definitions=======================================================================
// Define the precision of real numbers, could be real/double.
#define real double
#define Pi 3.1415926535897932384626433832795
#define Zero 0
//typedef double2 Complex;
using namespace std;

struct Particle {
    real* x;//save x position in GPU
    real* y;//save y position in GPU

    
    real* px;//save x dirction in GPU  保存方向x
    real* py;//save y dirction in GPU  保存方向y
    real* fx;//force on the x direction
    real* fy;//force on the y direction
    real* x0ToUpdateHybridList;//save xGpu[id] to judge whether update hybrid list 
    real* y0ToUpdateHybridList;//save yGpu[id] to judge whether update hybrid list
 
    int* cellX;//save xth cell of nth particle
    int* cellY;//save yth cell of nth particle

   

    int* cellList;//cell particle id for all particle, as [maxParticlePerCell*id + offsetsCL]?????
    int* cellOffsetsCL;//offset of every cell list to save particle number in this cell 
    int* particleAroundId;//save ids around this on particle, use rd to judge wether is "around"????????
    int* particleAroundFlagX;//mask whether cell of idth particle at the edge of box
    int* particleAroundFlagY;//mask whether cell of idth particle at the edge of box
    int* offsetsAL;//offset of every particle's around list
    int* offsetsNL;//offset of every particle's neighbor list to save neighbor particle id
    int* NeighborList;//neighbor list
    int* NeighborListFlagX;//translate from particleAroundFlagX
    int* NeighborListFlagY;//translate from particleAroundFlagY

    curandState* state;
} PT, pt;

struct Parameter {
    real boxX;//box size X
    real boxY;//box size Y
    real cellSizeX;//cell size in the x direction
    real cellSizeY;//cell size in the y direction
    int cellNumX;//num of cell in the x direction
    int cellNumY;//num of cell in the y direction
 

    int useForce;
    int ABParticle;
    // real Fx;
    // real Fy; //外力
    // int boundarysize; //这个是设置的边界的宽度，如果我们取9个cell做计算的话，那么只要1个cell的boundary就可以了
    

    real rho;		//密度
    int maxParticlePerCell;//theory maxmum particle number in one cell
    real rd;//deadline distance to get NeighborList
    int mask0;//use for bit calculate
    int mask1;//use for bit calculate 	
    real miniInstanceBetweenParticle;//theory minimum distance from two particle
    real r0;//balance position
    real epsilon;//coefficient of force
    float kBT;//kB*T
    real gammaValue;//Viscosity coefficien
    real rOutUpdateList;//update hybrid list when any one particle move a distance greater than rOutUpdateList
    int particleNum; //粒子数目
    real tStart;
    real tStop;
    real tStep;
    real tExpo;
    int lenBit;
    unsigned long long seed;
    int blockNum;
    int threadNum;
} PM;

__device__ int updateListFlag = 0;
__device__ int wrongFlag = 0;
int updateListFlagHost = 0;
int wrongFlagHost = 0;


//output===========================================================================================

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

//===========================================================================
void MemFree() {
    // Free host memory
    delete[] pt.x;
    delete[] pt.y;

    delete[] pt.px;
    delete[] pt.py;



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

void getInput(Parameter PM) {
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
        case 22: iss >> PM.useForce; break;
        case 23: iss >> PM.ABParticle; break;

        default: break; // 超过预期行数时不处理
        }
        lineCount++;
    }

    InputFile.close();
}

//mem ============================================================================================================

void MemAlloc() {
    // Allocate particle mem in host memory.
    pt.x = new real[PM.particleNum];
    pt.y = new real[PM.particleNum];


    pt.px = new real[PM.particleNum];
    pt.py = new real[PM.particleNum];



    pt.cellList = new int[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
    pt.cellOffsetsCL = new int[PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell];
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
    cudaMalloc((void**)&PT.cellOffsetsCL, PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell * sizeof(int));
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
__global__ void initState(curandState* state,unsigned long long seed, int particleNum) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= particleNum)return;
    curand_init(seed, id, 0, &state[id]);
}

//==============================================================================
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
    std::cout << "useForce: " << PM.useForce << std::endl;
    std::cout << "ABParticle: " << PM.ABParticle << std::endl;
}

//=============================================================================
__device__ int sign(real x) {
    return -(x < 0.f) + (x > 0.f);
}

__device__ int sign01(real x) {
    return (sign(x) + 1) / 2;
}

//===================================================================================
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
void InitOffset() {
    cudaMemset(PT.cellOffsetsCL, 0, sizeof(int) * PM.cellNumX * PM.cellNumY * PM.maxParticlePerCell);
    cudaMemset(PT.offsetsNL, 0, sizeof(int) * PM.particleNum);
    cudaMemset(PT.offsetsAL, 0, sizeof(int) * PM.particleNum);
}


//===============================================================
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

//==========================================================================================================
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
        printf("wrong");//append cout error later
    }
}

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
void listUpdate(Particle PT,Parameter PM) {
    InitOffset();

    getCellList << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
    getAroundCellParticleId << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
    saveXY0ToUpdateHybridList << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
}

//==================================================================================
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
__global__ void getForce (Particle PT, Parameter PM, *ABParticle) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;


    if (id >= PM.particleNum)return;
    real x0, y0, x1, y1, dx, dy, dr, f12;


    
    PT.fx[id] = 0;
    PT.fy[id] = 0;
    int i;
    for (i = 0;i < PT.offsetsNL[id];i++) {
        x0 = PT.x[id];
        y0 = PT.y[id];
        x1 = PT.x[PT.NeighborList[id * PM.maxParticlePerCell + i]];
        y1 = PT.y[PT.NeighborList[id * PM.maxParticlePerCell + i]];
        dx = sign01(0.5 * PM.boxX - x0 + x1) * sign01(0.5 * PM.boxX + x0 - x1) * (x0 - x1) + \
            sign01(sign(x0 - x1) * (x0 - x1) - 0.5 * PM.boxX) * -sign(x0 - x1) * (PM.boxX - sign(x0 - x1) * (x0 - x1));
        dy = sign01(0.5 * PM.boxY - y0 + y1) * sign01(0.5 * PM.boxY + y0 - y1) * (y0 - y1) + \
            sign01(sign(y0 - y1) * (y0 - y1) - 0.5 * PM.boxY) * -sign(y0 - y1) * (PM.boxY - sign(y0 - y1) * (y0 - y1));
        dr = sqrt(dx * dx + dy * dy);





        f12 = 24 * PM.epsilon * pow(PM.r0, 6) * (2 * pow(PM.r0, 6) - pow(dr, 6)) / pow(dr, 14);
        PT.fx[id] += f12 * dx ;
        PT.fy[id] += f12 * dy ;
        if (PT.fx[id] > 10000 || PT.fx[id] < -10000 || PT.fy[id] > 10000 || PT.fy[id] < -10000) {
            break;
        }
    }

    if (PM.useForce ==1) {
        int cellId = PT.cellY[id] * PM.cellNumX + PT.cellX[id];//此时未考虑周期边界条件
        real forceX,forceY;
        forceX = ABParticle.Qxx[cellID];
        forceY = ABParticle.Qxy[cellID];
        PT.fx[id] += forceX;
        PT.fy[id] += forceY;




    };




    if (PT.fx[id] > 10000 || PT.fx[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,x0:%f,x1:%f,y0:%f,y1:%f,NLFX:%d\n", id,PT.fx[id], PT.fy[id], dx,dy,x0,x1, y0, y1, PT.NeighborListFlagX[id * PM.maxParticlePerCell + i]);
        wrongFlag = 1;
    }
    if (PT.fy[id] > 10000 || PT.fy[id] < -10000) {
        printf("wrong!!!!!!!!!id:%d,fx:%f,fy:%f,dx:%f,dy:%f,y0:%f,y1:%f,x0:%f,x1:%f,NLFY:%d\n", id, PT.fx[id], PT.fy[id],dx,dy,y0,y1, x0, x1, PT.NeighborListFlagY[id * PM.maxParticlePerCell + i]);
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
    PT.x[id] = fmod(PT.x[id] + (PT.fx[id] * PM.tStep + fT * FRx) / PM.gammaValue + PM.boxX, PM.boxX);
    PT.y[id] = fmod(PT.y[id] + (PT.fy[id] * PM.tStep + fT * FRy) / PM.gammaValue + PM.boxY, PM.boxY);
}

//=========================================================================
void forceAndPositionUpdate(Particle PT, Parameter PM) {
    getForce << <PM.blockNum, PM.threadNum >> > (PT, PM);
    
    cudaDeviceSynchronize();
    updatePosition << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaDeviceSynchronize();
}

//=========================================================================
void iterate(Particle PT,Parameter PM) {
    forceAndPositionUpdate(PT,PM);
    checkUpdate << <PM.blockNum, PM.threadNum >> > (PT, PM);
    cudaMemcpyFromSymbol(&updateListFlagHost, updateListFlag, sizeof(int));
    if (updateListFlagHost){
        listUpdate(PT, PM);
        updateListFlagHost = 0;
        cudaMemcpyToSymbol(updateListFlag, &updateListFlagHost, sizeof(int));
    }
}

//===================================================================
void initBlockAndThreadNum() {
    PM.threadNum = 256;
    PM.blockNum = (PM.particleNum + PM.threadNum - 1) / PM.threadNum;
    printf("blockNum:%d,threadNum:%d\n", PM.blockNum, PM.threadNum);
}

//==========================================================================================================
void showProgress(real tNow, real tStart, real tStop, clock_t clockNow, clock_t clockStart) {
    real progress = (tNow - tStart) / (tStop - tStart);
    real tUsed = double(clockNow - clockStart) / CLOCKS_PER_SEC;
    real tUsePrediction = (tStop - tNow) * tUsed / (tNow - tStart);
    printf("%.8f,%.8f  ", pt.x[0], pt.y[0]);
    printf("  Progress:%.4f%%,%Prediction:%.1f\r", progress*100, tUsePrediction);
    fflush(stdout);
}

//=======================================================
int main()
{
    real tNow = PM.tStart;
    getInput();
    MemAlloc();
    printInput();
    Init_Coords(1, pt, PM);
    initBlockAndThreadNum();
    InitOffset();
    initState << <PM.blockNum, PM.threadNum >> > (PT.state, PM.seed, PM.particleNum);


    
    cudaDeviceSynchronize();
    parameterInit();
    HostUpdataToDevice();
    PM.seed = static_cast<unsigned long long>(time(0));
    printf("seed:%d\n", PM.seed);
    ExpoConf("0");
       

    listUpdate(PT, PM);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
        return 0;
    }

    clock_t clockStart = clock();
    for (tNow = PM.tStart;tNow < PM.tStop;tNow += PM.tStep) {
        iterate(PT, PM);

        cudaMemcpyFromSymbol(&wrongFlagHost, wrongFlag, sizeof(int));
        if (wrongFlagHost == 1)return;

        //printf("---------------------------------------------tNow:%f\n",tNow);
        if (floor(tNow / PM.tExpo) > floor((tNow - PM.tStep) / PM.tExpo)) {
            showProgress(tNow, PM.tStart, PM.tStop, clock(), clockStart);
            DeviceUpdataToHost();//下载数据到主机
            int te = floor(tNow / PM.tExpo) + 1;
            string str_t = to_string(te);
            ExpoConf(str_t);
        }
    }
    MemFree();//释放内存

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }

    return 0; // 返回成功状态
}
