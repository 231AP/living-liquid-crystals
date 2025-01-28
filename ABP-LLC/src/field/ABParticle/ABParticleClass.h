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

#endif
