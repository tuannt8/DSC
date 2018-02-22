//
//  particle.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef particle_hpp
#define particle_hpp

#include <stdio.h>

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include <fstream>
#include <memory>
#include "KDTree.h"

#include "DSC.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#ifdef _WIN32 // WINDOWS
#include <GL/glut.h>
#include <GL/glew.h>
#elif defined(__APPLE__) // IOS
#include <OpenGL/gl3.h>
#include <GLUT/glut.h>
#else // LINUX
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include "anisotrpic_kernel.h"


class particle_manager{
public:
    particle_manager(){};
    ~particle_manager(){};
    
    // Data that need input
    std::string m_data_path;
    double m_influence_radius; // influent radius
    double m_deltap; // Spacing distance
    double m_slength; // Smoothing radius
    
    // Particles
    std::vector<particle> m_current_particles; // Current file
    std::vector<particle> m_next_particles; // Next file
    
    void init_first(int idx);
    void load_next(int idx, double t);
    
    int m_cache_idx;
    std::vector<particle> m_cache_particles; // Current file
    std::vector<particle> m_cache_particles_next; // Next file
    
//    int m_cur_sub_step, m_max_step; // estimated by max displacement
//    std::vector<particle> m_sub_step_particles; // When subdivide to sub time steps
//    std::vector<vec3> m_sub_step_vel; // Linear approximation
    
    // accelerating finding neighbor
    Geometry::KDTree<vec3, int> m_vtree;
    void build_kd_tree();
    void rebuild_density();
    
    // Aniso tropic kernel
    anisotropic_kernel m_aniso_kernel;
    void build_anisotropic_kernel();
    bool get_projection(vec3 pos, vec3 direction, bool &bInside, double &t);//Using anisotropic kernel
    
    // Load next time step
    void load_first_time(int idx);
    void load_time_step(int idx);
    void interpolate(int sub_idx, int sub_count);
    void load_time_step();
    void load(int idx, std::vector<particle> & par);
    
    // Estimate displacements
    bool get_displacement(vec3 pos, vec3 & dis);
    bool get_displacement_avg(vec3 pos, vec3 & dis);
    bool get_displacement_weighted_avg(vec3 pos, vec3 & dis);
    bool get_displacement_sph_kernel(vec3 pos, vec3 & dis);
    vec3 get_displacement_closet_point(vec3 pos);
    vec3 get_displacement_cubic_kernel(vec3 pos);
    bool get_displacement_WENLAND_kernel(vec3 pos, vec3 & dis);
    bool get_displacement_MLS_kernel(vec3 pos, vec3 & dis);
    
    enum weight_type
    {
        POLY_4 = 0
    };
    double weight_function(double r, int type = POLY_4);
    
    // Get max displacement between two time steps
    double get_max_displacement();
    
    void draw(double yunder = -INFINITY, double ylimit = INFINITY);
    void draw_intermediate_vel();
    
    std::vector<vec3> pos;
    std::vector<double> phi;
    void draw_anisotropic_kernel(vec3 domain_size, vec3 c);
    void draw_orientation_anisotropic();
    void draw_anisotropic_kernel(double yunder, double yupper);
};

#endif /* particle_hpp */
