#ifndef ABPARTICLEFIELDCLASS_H
#define ABPARTICLEFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "LLCParams.cu"


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

    int* offsetsAL;//offset of every particle's around list
    int* offsetsNL;//offset of every particle's neighbor list to save neighbor particle id
    int* NeighborList;//neighbor list
    int* NeighborListFlagX;//translate from particleAroundFlagX
    int* NeighborListFlagY;//translate from particleAroundFlagY


    int* cellX1;//save xth cell of nth particle
    int* cellY1;//save yth cell of nth particle
    int* cellList1;//cell particle id for all particle, as [maxParticlePerCell*id + offsetsCL]?????
    int* cellOffsetsCL1;//offset of every cell list to save particle number in this cell 

    int* offsetsAL1;//offset of every particle's around list
    int* offsetsNL1;//offset of every particle's neighbor list to save neighbor particle id
    int* NeighborList1;//neighbor list
    int* NeighborListFlagX1;//translate from particleAroundFlagX
    int* NeighborListFlagY1;//translate from particleAroundFlagY
    curandState* state;
}PT,pt;






// ===========================================================



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
    void ParticleToField(int i_field,LLCParams &params);
    // void getConcentration(int i_field);
    void forceAndPositionUpdate(Particle PT, LLCParams &params ,int i_field);
    void iterate(Particle PT, LLCParams &params,int i_field);
    void CellPUpdate(Particle PT,LLCParams &params, int i_field);
    
    // ===========================================================
};


__device__ int updateListFlag = 0;
__device__ int wrongFlag = 0;
int updateListFlagHost = 0;
int wrongFlagHost = 0;




// =======================================================================


   


// methods
    void ExpoConf(const std::string& str_t, LLCParams &params);
    void MemFree();
    void MemAlloc();

    
    void Init_Coords(int flag, Particle pt, LLCParams &params);
    void InitOffset();
    void parameterInit();
    void HostUpdataToDevice(LLCParams &params);
    void DeviceUpdataToHost(LLCParams &params);
    void listUpdate(Particle PT,LLCParams &params);

    


#endif
