//
//  anisotrpic_kernel.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 04/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef anisotrpic_kernel_h
#define anisotrpic_kernel_h

#include "KDTree.h"
#include "define.h"
#include "util.h"

class particle
{
public:
    particle(){}
    ~particle(){}
    
    double pressure;
    double density;
    double mass;
    int type;
    int flag;
    int object; // fluid number?
    double volume;
    double sigma;
    vec3 vel, pos;
    
    void draw()
    {
//        static std::vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
//
//        glBegin(GL_POINTS);
//        glColor3dv(_color[type].get());
//        glVertex3dv(pos.get());
//        glEnd();
    }
};


class anisotropic_kernel
{
public:
    Geometry::KDTree<vec3, int> m_vtree;
    std::vector<int> m_connected_component_label;
    int m_nb_component;
    double m_h;// Smoothing radius. Used for connected component
    double m_ra; // Spacing distance betwwen particle
    
    double m_r;// Use for everything.
    
    std::vector<mat3x3d> m_G;
    std::vector<double> m_det_G;
    std::vector<bool> m_b_kernel_computed;
    std::vector<particle> *m_shared_particles;
    
    
    const mat3x3d & get_transform_mat(int idx);
public:
    anisotropic_kernel(){};
    ~anisotropic_kernel(){};
    
    void build();
    double weight_func(int i, int j, double radius);
    double weight_func(int i, vec3 posi, int j, double);
    void Taubin_smooth();
    void compute_kd_tree();
    void build_connected_component();
    
    double get_value(vec3 pos);
    
    void draw_connected_component();
    void compute_tranformation_mat_for_particle(int);
    
    std::vector<int> neighbor_search(vec3 pos, double radius);
    
public:
    // Debugging
    std::vector<mat3x3d> m_U;
    std::vector<mat3x3d> m_S;
    std::vector<mat3x3d> m_C;
};

#endif /* anisotrpic_kernel_h */
