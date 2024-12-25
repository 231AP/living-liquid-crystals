#ifndef ABPARTICLECLASS_CU
#define ABPARTICLECLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "ABParticleClass.h"
#include "ABParticleClassGPU.cu"


using namespace std;

// =============================================================
// Constructors
ABParticle::ABParticle (Mesh* mesh_ptr_t, string name_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=0;
    boun_cond = "periodic";
    init_cond = "sin";
    location = "both";
    expo_data = "on";
    initFields ();
};


// -------------------------------------------------------------
ABParticle::ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = "sin";
    expo_data = "on";
    initFields ();    
};


// -------------------------------------------------------------
ABParticle::ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = init_cond_t;
    expo_data = "on";
    initFields ();
};

// -------------------------------------------------------------
ABParticle::ABParticle (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = boun_cond_t;
    init_cond = init_cond_t;
    expo_data = expo_data_t;
    initFields ();
};


// -------------------------------------------------------------
void ABParticle::initFields () {
    setFieldProperties(&Concentration, name+".Concentration", priority , init_cond);
    setFieldProperties(&Pxx, name+".Pxx", -1, init_cond);
    setFieldProperties(&Pxy, name+".Pxy", -1, init_cond);
    setFieldProperties(&C1, name+".C1", -1, init_cond);



    setFieldProperties(&Qxx,name+".Qxx",priority, init_cond);
    setFieldProperties(&Qxy,name+".Qxy",priority, init_cond);
    
    setFieldProperties(&vx,name+".vx",priority,init_cond);
    setFieldProperties(&vy,name+".vy",priority,init_cond);
    

    
    Concentration.specialty="ABParticle";
    Concentration.ptr_Pxx=&Pxx;
    Concentration.ptr_Pxy=&Pxy;

    Concentration.ptr_Qxx=&Qxx;
    Concentration.ptr_Qxy=&Qxy;
    Concentration.ptr_vx = &vx;
    Concentration.ptr_vy = &vy;
    Concentration.ptr_C1 = &C1;




    
    // No RHS for vx, vy, phi, as they need special functions to get values
    // Qxx.setRhsTerms({});
    // Qxy.setRhsTerms({});
    // vx.setRhsTerms({});
    // vy.setRhsTerms({});
    // Concentration.setRhsTerms({});
    // Pxx.setRhsTerms({});
    // Pxy.setRhsTerms({});
   
    
    
    
};


// -------------------------------------------------------------
void ABParticle::setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t) {
    (*field_ptr).traits_host.mesh_ptr=mesh_ptr;
    (*field_ptr).traits_host.name=field_name;
    (*field_ptr).traits_host.priority=priority_t;
    (*field_ptr).traits_host.boun_cond = boun_cond;
    (*field_ptr).traits_host.init_cond = init_cond_t;
    (*field_ptr).traits_host.expo_data = expo_data;
    (*field_ptr).num_f_funcs=0;
    
    (*field_ptr).initFieldAddi();
    for (int i=0; i<200; i++) {
        (*field_ptr).f_funcs_host[i]=NULL;
    };
};


// // =============================================================

#endif
