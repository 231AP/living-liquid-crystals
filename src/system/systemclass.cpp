#ifndef SYSTEMCLASS_CPP
#define SYSTEMCLASS_CPP

#include "systemclass.h"
#include <iostream> 
#include <vector>

using namespace std; 


// ======================================================================
void System::printSysInfo () {
    cout << "System's info:" <<endl;
    // Printing grid info.
    if (mesh_ptr != NULL) {
        cout <<"This is a "<<(*mesh_ptr).host.space_dim<<" dimensional system.\n";
        switch ((*mesh_ptr).host.space_dim) {
            case 1:
                cout <<"The box size is: "<<(*mesh_ptr).host.box_size.x<<".\n";
                cout <<"There are "<<(*mesh_ptr).host.grid_number.x<<" grid points, ";
                cout <<"with "<<(*mesh_ptr).host.grid_number_boun.x<<" layers of boundary grids.\n";
                break;
            case 2:
                cout <<"The box size is "<<(*mesh_ptr).host.box_size.x;
                cout <<"*"<<(*mesh_ptr).host.box_size.y<<".\n";
                cout <<"There are "<<(*mesh_ptr).host.grid_number.x<<"*";
                cout <<(*mesh_ptr).host.grid_number.y<<" grid points, ";
                cout <<"and "<<(*mesh_ptr).host.grid_number_boun.x<<"*";
                cout <<(*mesh_ptr).host.grid_number_boun.y<<" layers of boundary grids.\n";
                break;
            case 3:
                cout <<"The box size is "<<(*mesh_ptr).host.box_size.x;
                cout <<"*"<<(*mesh_ptr).host.box_size.y;
                cout <<"*"<<(*mesh_ptr).host.box_size.z<<".\n";
                cout <<"There are "<<(*mesh_ptr).host.grid_number.x<<"*";
                cout <<(*mesh_ptr).host.grid_number.y<<"*";
                cout <<(*mesh_ptr).host.grid_number.z<<" grid points, ";
                cout <<"and "<<(*mesh_ptr).host.grid_number_boun.x<<"*";
                cout <<(*mesh_ptr).host.grid_number_boun.y<<"*";
                cout <<(*mesh_ptr).host.grid_number_boun.z<<" layers of boundary grids.\n";
        };
        // cout <<(*mesh_ptr).host.
    } else {
        cout <<"The system contains NO mesh!" <<endl;
    };
    
    // Printing field info.
    if (field_ptrs.size()==0) {
        cout << "There are no field in this system."<<endl;
    } else {
        cout << "There are " <<field_ptrs.size()<<" fields: ";
        for (auto iter = field_ptrs.begin(); iter != field_ptrs.end(); ) {        
            std::cout << (*iter)->name();
            if ( ++iter == field_ptrs.end() ) {
                cout<<"." <<endl;
            } else {
                cout <<", ";
            };
        };
    };
    
};


// ------------------------------------------------------------
void System::addField (Field* f_ptr_t) {    
    field_ptrs.push_back(f_ptr_t);
    (*f_ptr_t).traits_host.dire_expo=dire_expo;
};


// ------------------------------------------------------------
void System::addIncompFlow (IncompFlow* incompFlow_ptr_t) {
    incompFlow_ptr=incompFlow_ptr_t;
    addField(&(*incompFlow_ptr).omega);
    addField(&(*incompFlow_ptr).phi);
    addField(&(*incompFlow_ptr).vx);
    addField(&(*incompFlow_ptr).vy);
};

// ------------------------------------------------------------
void System::addConstantField (ConstantField* constF_ptr_t) {
    constF_ptr=constF_ptr_t;
    addField(&(*constF_ptr).one);
    addField(&(*constF_ptr).zero);
    addField(&(*constF_ptr).pi);
};

// ------------------------------------------------------------
void System::addLivingLC (LivingLC* LivingLC_ptr_t) {
    addField(&(*LivingLC_ptr_t).Pxx);
    addField(&(*LivingLC_ptr_t).Pxy);
    addField(&(*LivingLC_ptr_t).theta);
    addField(&(*LivingLC_ptr_t).theta_old);
    addField(&(*LivingLC_ptr_t).px);
    addField(&(*LivingLC_ptr_t).py);
    addField(&(*LivingLC_ptr_t).flip);
    addField(&(*LivingLC_ptr_t).cplus);
    addField(&(*LivingLC_ptr_t).cminus);
    addField(&(*LivingLC_ptr_t).C);
    addField(&(*LivingLC_ptr_t).Omega);
    
};
// =======================================================================

void System::addABParticle (ABParticle* ABParticle_ptr_t) {
    addField(&(*ABParticle_ptr_t).vx);
    addField(&(*ABParticle_ptr_t).vy);
    addField(&(*ABParticle_ptr_t).Qxx);
    addField(&(*ABParticle_ptr_t).Qxy);
    addField(&(*ABParticle_ptr_t).Concentration);
    addField(&(*ABParticle_ptr_t).Pxx);
    addField(&(*ABParticle_ptr_t).Pxy);
    
};
// =======================================================================
#endif
