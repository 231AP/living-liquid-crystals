#ifndef ABPARTICLECLASS_H
#define ABPARTICLECLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "ABParticleFieldClass.cu"


using namespace std; 


// ===============================================================
class ABParticle {
    
    public:

    string name;
    int priority=0;
    string boun_cond = "periodic";
    string init_cond = "sin";
    int rank=0;
    string location = "both";
    string expo_data = "on";
    string equation = "";
    Mesh* mesh_ptr=NULL;
    string FDMScheme;

    
    ABParticleField Concentration;
    Field Qxx;
    Field Qxy;
    Field Pxx;
    Field Pxy;
    Field vx;
    Field vy;
    Field C1;


    // Field;
  
    
    



    
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    ABParticle (Mesh* mesh_ptr_t, string name_t);
    ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t);
    ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);

    // Fields
    void initFields ();
    void setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t);
 
    // ===========================================================
};







// ===========================================================







// struct Particle {
//     real* x;//save x position in GPU
//     real* y;//save y position in GPU

    
//     real* px;//save x dirction in GPU  保存方向x
//     real* py;//save y dirction in GPU  保存方向y
//     real* fx;//force on the x direction
//     real* fy;//force on the y direction

//     real* Fx;
//     real* Fy; //force which from field

//     real* x0ToUpdateHybridList;//save xGpu[id] to judge whether update hybrid list 
//     real* y0ToUpdateHybridList;//save yGpu[id] to judge whether update hybrid list
 
//     int* cellX;//save xth cell of nth particle
//     int* cellY;//save yth cell of nth particle

   

//     int* cellList;//cell particle id for all particle, as [maxParticlePerCell*id + offsetsCL]?????
//     int* cellOffsetsCL;//offset of every cell list to save particle number in this cell 
//     int* particleAroundId;//save ids around this on particle, use rd to judge wether is "around"????????
//     int* particleAroundFlagX;//mask whether cell of idth particle at the edge of box
//     int* particleAroundFlagY;//mask whether cell of idth particle at the edge of box
//     int* offsetsAL;//offset of every particle's around list
//     int* offsetsNL;//offset of every particle's neighbor list to save neighbor particle id
//     int* NeighborList;//neighbor list
//     int* NeighborListFlagX;//translate from particleAroundFlagX
//     int* NeighborListFlagY;//translate from particleAroundFlagY

//     curandState* state;


    






// }PT,pt;






// // ===========================================================


// struct Parameter {
//     real boxX;//box size X
//     real boxY;//box size Y
//     real cellSizeX;//cell size in the x direction
//     real cellSizeY;//cell size in the y direction
//     int cellNumX;//num of cell in the x direction
//     int cellNumY;//num of cell in the y direction
 

//     int V0;
//     int ABParticle;
//     // int boundarysize; //这个是设置的边界的宽度，如果我们取9个cell做计算的话，那么只要1个cell的boundary就可以了
    

//     real rho;		//密度
//     int maxParticlePerCell;//theory maxmum particle number in one cell
//     real rd;//deadline distance to get NeighborList
//     int mask0;//use for bit calculate
//     int mask1;//use for bit calculate 	
//     real miniInstanceBetweenParticle;//theory minimum distance from two particle
//     real r0;//balance position
//     real epsilon;//coefficient of force
//     float kBT;//kB*T
//     real gammaValue;//Viscosity coefficien
//     real rOutUpdateList;//update hybrid list when any one particle move a distance greater than rOutUpdateList
//     int particleNum; //粒子数目
//     real tStart;
//     real tStop;
//     real tStep;
//     real tExpo;
//     int lenBit;
//     unsigned long long seed;
//     int blockNum;
//     int threadNum;
// }PM;

// __device__ int updateListFlag = 0;
// __device__ int wrongFlag = 0;
// int updateListFlagHost = 0;
// int wrongFlagHost = 0;




// // =======================================================================


// void getInput();


// // methods
//     void ExpoConf(const std::string& str_t);
//     void MemFree();
   
//     void MemAlloc();
//     void printInput();
//     void MemAlloc();
//     void printInput();
//     void Init_Coords(int flag, Particle pt, Parameter PM);
//     void InitOffset();
//     void parameterInit();
//     void HostUpdataToDevice();
//     void listUpdate(Particle PT,Parameter PM);
//     void forceAndPositionUpdate(Particle PT, Parameter PM);
//     void iterate(Particle PT,Parameter PM);
//     void initBlockAndThreadNum();
//     void showProgress(real tNow, real tStart, real tStop, clock_t clockNow, clock_t clockStart);
    
#endif
