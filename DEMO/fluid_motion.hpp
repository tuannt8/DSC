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

#include "anisotrpic_kernel.h"

#define DT_NORM 0.2



class fluid_motion
{
public:
    fluid_motion(); // dangerous to use with static variable. Use init() after the static variables are set
    ~fluid_motion(){};
    
    void load_configuration();
    void init(DSC::DeformableSimplicialComplex<> *dsc);
    
    void draw();
    
    static std::string m_data_path;
    
private:
    std::vector<std::string> m_out_path; // Write surfaces for different phases
    
    double m_max_dsc_displacement;
    
    double m_max_displacement_projection; // for binary search
public:
    DSC::DeformableSimplicialComplex<>* s_dsc; // shared DSC
    std::unique_ptr<problem> m_problem; // Will be casted to specific problem
    std::vector<std::shared_ptr<particle_manager>> m_particles;
    anisotropic_kernel m_share_aniso_kernel;
    
    int m_cur_global_idx = 0;
    double t = 0, dt = 1;
    
    void load_first_particle();
    void load_next_particle(); // return true if loading new file
    void add_ghost_particles();
    
    void compute_advection(std::vector<vec3> & vertex_dis);
    
    void deform();
    void project_interface_itteratively();
    void project_interface_one_iter();
    int project_interface( );
    void project_interface_test();
    void snapp_boundary_vertices();
    void advect_velocity();
    
    void log_dsc();
    void log_dsc_surface();
    void extract_surface_phase(int phase, std::string path);
    
    std::vector<bool> is_vertices_boundary;
    void update_vertex_boundary();
    bool is_boundary_work_around(is_mesh::FaceKey fkey);
    
    void reset_projected_flag();
    
    void build_anisotropic_kernel();
    void draw_anisotropic_kernel_plane();
};



#endif /* fluid_motion_hpp */
