//
//  problem.h
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef problem_h
#define problem_h

#include <string>
#include "DSC.h"
#include "define.h"

class problem{
public:
    problem(){}; 
    ~problem(){};
    
// SPH parametter
    double m_deltap; // distance between particle
    double m_influenceRadius; // Radius for searching neighbor
    double m_slength; // Smooth length, for computing the kernel
    int m_nb_phases;
    
// related function
    DSC::DeformableSimplicialComplex<> * init_dsc_domain(double);
    virtual DSC::DeformableSimplicialComplex<> * init_dsc(double scale = 1)=0; // Pure virtual function
    virtual vec3 domain_size()=0;
    
    void init(std::string sumary_file_path);
};

class two_phase_fluid:public problem
{
public:
    two_phase_fluid(){};
    ~two_phase_fluid(){};
    
    virtual DSC::DeformableSimplicialComplex<> * init_dsc(double scale = 1);
    virtual vec3 domain_size(){return vec3(0.15, 0.15, 0.15);};
};

class dam_break_fluid:public problem
{
public:
    dam_break_fluid(){};
    ~dam_break_fluid(){};
    
    virtual DSC::DeformableSimplicialComplex<> * init_dsc(double scale = 1);
    virtual vec3 domain_size(){return vec3(1.6, 0.67, 0.6);};
};


#endif /* problem_h */
