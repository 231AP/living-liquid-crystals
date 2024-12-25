#ifndef ABPARTICLEFIELDCLASS_H
#define ABPARTICLEFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 

struct Particle {
    real* x;//save x position in GPU
    real* y;//save y position in GPU


    real* px;//save x dirction in GPU  保存方向x
    real* py;//save y dirction in GPU  保存方向y
    real* fx;//force on the x direction
    real* fy;//force on the y direction

    real* Fx;
    real* Fy; //force which from field

    real* x0ToUpdateHybridList;//save xGpu[id] to judge whether update hybrid list 
    real* y0ToUpdateHybridList;//save yGpu[id] to judge whether update hybrid list
    real* x0ToUpdateHybridList1;//save xGpu[id] to judge whether update hybrid list 
    real* y0ToUpdateHybridList1;//save yGpu[id] to judge whether update hybrid list

    real* cellPx;
    real* cellPy;
    real* cellPxx;
    real* cellPxy;
    
 
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


    int* cellX1;//save xth cell of nth particle
    int* cellY1;//save yth cell of nth particle
    int* cellList1;//cell particle id for all particle, as [maxParticlePerCell*id + offsetsCL]?????
    int* cellOffsetsCL1;//offset of every cell list to save particle number in this cell 
    int* particleAroundId1;//save ids around this on particle, use rd to judge wether is "around"????????
    int* particleAroundFlagX1;//mask whether cell of idth particle at the edge of box
    int* particleAroundFlagY1;//mask whether cell of idth particle at the edge of box
    int* offsetsAL1;//offset of every particle's around list
    int* offsetsNL1;//offset of every particle's neighbor list to save neighbor particle id
    int* NeighborList1;//neighbor list
    int* NeighborListFlagX1;//translate from particleAroundFlagX
    int* NeighborListFlagY1;//translate from particleAroundFlagY
    curandState* state;
}PT,pt;






// ===========================================================


struct Parameter {
    real boxX;//box size X
    real boxY;//box size Y
    real cellSizeX;//cell size in the x direction
    real cellSizeY;//cell size in the y direction
        real cellSizeX1;//cell size in the x direction
    real cellSizeY1;//cell size in the y direction
    int cellNumX;//num of cell in the x direction
    int cellNumY;//num of cell in the y direction
    int cellNumX1;//num of cell in the x direction
    int cellNumY1;//num of cell in the y direction
    
    int V0;
    int ABParticle;
    // int boundarysize; //这个是设置的边界的宽度，如果我们取9个cell做计算的话，那么只要1个cell的boundary就可以了
    

    real rho;		//密度
    int maxParticlePerCell;//theory maxmum particle number in one cell
    int maxParticlePerCell1;//theory maxmum particle number in one cell
    real rd;//deadline distance to get NeighborList
    int mask0;//use for bit calculate
    int mask1;//use for bit calculate 	
    real miniInstanceBetweenParticle;//theory minimum distance from two particle
    real r0;//balance position
    real epsilon;//coefficient of force
    float kBT;//kB*T
    real gammaValue;//Viscosity coefficien
    real rOutUpdateList;//update hybrid list when any one particle move a distance greater than rOutUpdateList
    real rOutUpdateList1;//update hybrid list when any one particle move a distance greater than rOutUpdateList
    int particleNum; //粒子数目
    real tStart;
    real tStop;
    real tStep;
    real tExpo;
    int lenBit;
    unsigned long long seed;
    int blockNum;
    int threadNum;
}PM;
// ===============================================================
class ABParticleField :public Field {
    
public:
          
    // ===========================================================
    Field* ptr_Pxx;
    Field* ptr_Pxy;
    // Field* ptr_Concentration;
    Field* ptr_Qxx;
    Field* ptr_Qxy;   
    Field* ptr_vx;
    Field* ptr_vy;
    Field* ptr_C1;

   
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    ABParticleField () {};
    ABParticleField (Mesh* mesh_ptr_t, string name_t);
    ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    ABParticleField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Get velocity field
    void initPolarField();
    void ParticleToField(int i_field);
    // void getConcentration(int i_field);
    void forceAndPositionUpdate(Particle PT, Parameter PM,int i_field);
    void iterate(Particle PT, Parameter PM,int i_field);
    void CellPUpdate(Particle PT,Parameter PM, int i_field);
    
    // ===========================================================
};


__device__ int updateListFlag = 0;
__device__ int wrongFlag = 0;
int updateListFlagHost = 0;
int wrongFlagHost = 0;




// =======================================================================


   


// methods
    void getInput();
    void ExpoConf(const std::string& str_t);
    void MemFree();
   
    void MemAlloc();
    void printInput();

    
    void Init_Coords(int flag, Particle pt, Parameter PM);
    void InitOffset();
    void parameterInit();
    void HostUpdataToDevice();
    void listUpdate(Particle PT,Parameter PM);
    // void forceAndPositionUpdate(Particle PT, Parameter PM);
    // void iterate(Particle PT,Parameter PM);
    void initBlockAndThreadNum();
    void showProgress(real tNow, real tStart, real tStop, clock_t clockNow, clock_t clockStart);
    


#endif
