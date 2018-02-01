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
    double m_influence_radius;
    double m_deltap;
    double m_slength;
    
    // Particles
    int m_cur_idx = 0;
    std::vector<particle> m_current_particles; // Current file
    std::vector<particle> m_next_particles; // Next file
    
    int m_cur_sub_step, m_max_step; // estimated by max displacement
    std::vector<particle> m_sub_step_particles; // When subdivide to sub time steps
    std::vector<vec3> m_sub_step_vel; // Linear approximation
    
    // accelerating finding neighbor
    Geometry::KDTree<vec3, int> m_vtree;
    void build_kd_tree();
    
    // Aniso tropic kernel
    anisotropic_kernel m_aniso_kernel;
    void build_anisotropic_kernel();
    bool get_projection(vec3 pos, vec3 direction, bool &bInside, double &t);//Using anisotropic kernel
    
    // Load next time step
    void load_time_step(int idx);
    void interpolate(int sub_idx, int sub_count);
    void load_time_step();
    void load(int idx, std::vector<particle> & par);
    
    // Estimate displacements
    bool get_displacement(vec3 pos, vec3 & dis);
    bool get_displacement_avg(vec3 pos, vec3 & dis);
    vec3 get_displacement_closet_point(vec3 pos);
    vec3 get_displacement_cubic_kernel(vec3 pos);
    bool get_displacement_WENLAND_kernel(vec3 pos, vec3 & dis);
    
    // Some other
    double get_max_displacement();
    void draw();
};

#endif /* particle_hpp */
