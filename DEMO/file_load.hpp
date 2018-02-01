//
//  vtkWraper.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef vtkWraper_hpp
#define vtkWraper_hpp

#include <stdio.h>
#include "define.h"

#define MAX_ITER 300

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

#include "particle_manager.hpp"
#include "anisotrpic_kernel.h"
#include "glut_menu.hpp"

//#define DEMO_INTERFACE

class hash3
{
public:
    hash3(vec3 domain_bound, double cell_size);
    ~hash3(){};
    
    void insert_point(vec3 pos, int index);
    std::vector<long> get_close_point(double x, double y, double z, double radius);
    
    void draw();
private:
    std::map<long, std::vector<long>> m_bins;
    vec3i m_dimension;
    double m_cell_size;
    
    inline int get_idx_cell(vec3 & pos);
    inline vec3i get_idx_cell3(vec3& pos)
    {
        return vec3i(floor(pos[0]/m_cell_size),
                    floor(pos[1]/m_cell_size),
                    floor(pos[2]/m_cell_size) );
    }
    inline int idx_int(vec3i idx_h)
    {
        return idx_h[2]*m_dimension[0]*m_dimension[1] + idx_h[1]*m_dimension[0] + idx_h[0];
    }
};

// This class loads obj object and does some processing
class fluid_interface
{
    std::vector<vec3> m_points;
    std::vector<int> m_faces;
public:
    fluid_interface(){};
    ~fluid_interface(){};
    
    void load_surface(int idx)
    {
        m_points.clear();
        m_faces.clear();
        
        std::stringstream file_path;
        file_path << "../Large_data/DamBreak3D/mesh/dam_break_1_obstacle/dsc_" << std::setfill('0') << std::setw(5) << idx << ".obj";
        is_mesh::import_surface_mesh(file_path.str(), m_points, m_faces);
    };
    
    void project_boundary(vec3 ld, vec3 ru, double epsilon)
    {
        for(auto & p : m_points)
        {
            for (int i = 0; i < 3; i++)
            {
                if(p[i] < ld[i])
                    p[i] = ld[i];
                
                if(p[i] > ru[i])
                    p[i] = ru[i];
                
                if (p[i] - ld[i] < epsilon)
                {
                    p[i] = ld[i];
                }
                if (ru[i] - p[i] < epsilon)
                {
                    p[i] = ru[i];
                }
            }
        }
    }
    
    void load_surface(std::string path, int idx)
    {
        m_points.clear();
        m_faces.clear();
        
        std::stringstream file_path;
        file_path << path << std::setfill('0') << std::setw(5) << idx << ".obj";
        is_mesh::import_surface_mesh(file_path.str(), m_points, m_faces);
    };
    
    void write_surface(std::string path, int idx)
    {
        std::stringstream file_path;
        file_path << path << std::setfill('0') << std::setw(5) << idx << ".obj";
        
        is_mesh::export_surface_mesh(file_path.str(), m_points, m_faces);
    }
    
    void draw();
};

class file_load
{
public:
    file_load(){};
    ~file_load();
    
    std::vector<particle> m_current_particles;
    std::vector<particle> m_next_particles;
    
    std::vector<particle> m_sub_step_particles;
    std::vector<vec3> m_sub_step_vel;
    
    int m_cur_idx = 0;
    int m_cur_sub_step = 0;
#ifdef __APPLE__
    int m_max_step=1;
#else
    int m_max_step=1;
#endif
    
    std::shared_ptr<hash3> m_hashTable;
    Geometry::KDTree<vec3, int> m_vtree;
    void build_hash();
    
    anisotropic_kernel m_aniso_kernel;
    void build_anisotropic_kernel();
    
    void load_time_step();
    void load(int idx, std::vector<particle> & par);
    
    void draw();
    
    bool get_displacement(vec3 pos, vec3 & dis);
    bool get_displacement_avg(vec3 pos, vec3 & dis);
    vec3 get_displacement_closet_point(vec3 pos);
    vec3 get_displacement_cubic_kernel(vec3 pos);
    bool get_displacement_WENLAND_kernel(vec3 pos, vec3 & dis);
    
    void fix_output_boundary();
    
    bool get_projection(vec3 pos, vec3 direction, bool &bInside, double &t);//Using anisotropic kernel
    
    std::string m_data_path;
    
    virtual void init_dsc(DSC::DeformableSimplicialComplex<> * dsc)=0;//{};
    virtual vec3 get_domain_dimension(){return vec3(0.0);};
    virtual double get_influence_radius(){return 0;};
    virtual double get_spacing_distance()=0;
    virtual void personal_draw(){};
    
public:
    fluid_interface m_interface;
};


#endif /* vtkWraper_hpp */
