//
//  fluid_motion.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef fluid_motion_hpp
#define fluid_motion_hpp

#include <stdio.h>
#include "file_load.hpp"
#include "DSC.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "dam_break.h"
#include "bubble.h"
#include "dambreak2.h"
#include "dambreak_high_res.h"

#include "problem.h"
#include "particle_manager.hpp"

#define DT_NORM 0.1

class fluid_motion
{
public:
    fluid_motion(); // dangerous to use with static variable. Use init() after the static variables are set
    ~fluid_motion(){};
    
    void init();
    
    void draw();
    
    static std::string m_data_path;
    
private:
    std::vector<std::string> m_out_path; // Write surfaces for different phases
    
public:
    DSC::DeformableSimplicialComplex<>* s_dsc; // shared DSC
    std::unique_ptr<problem> m_problem; // Will be casted to specific problem
    std::vector<std::shared_ptr<particle_manager>> m_particles;
    
    void load_next_particle();
    int subdivide_time_step();
//    dam_break m_file_load;
    
    void deform();
    void project_interface();
    
    void log_dsc_surface(int idx);
};



#endif /* fluid_motion_hpp */
