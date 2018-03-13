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
    
    double density;
    double mass;
    long fluid;
    vec3 vel, pos;
};


class anisotropic_kernel
{
public:
    Geometry::KDTree<vec3, int> m_vtree;
    std::vector<int> m_connected_component_label;
    int m_nb_component;
    double m_h;// Smoothing radius. Used for connected component
//    double m_ra; // Spacing distance betwwen particle
    
    double m_r;// Use for everything. = 2*m_h
    
    
    
    std::vector<mat3x3d> m_G;
    std::vector<double> m_det_G;
    std::vector<bool> m_b_kernel_computed;
    std::vector<particle> m_particles;
    
    
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
    
    // 0.0 means outside
    //  Inside is larger than 0
    double get_value(vec3 pos);
    bool get_projection(vec3 pos, vec3 direction, bool &bInside, vec3& projected_point);
    vec3 get_displacement_projection(vec3 pos, vec3 norm, double max_displace, bool & bLast);
    
    vec3 get_displacement_projection(vec3 pos, vec3 norm, double max_displace);
    bool is_inside(vec3 pos, std::vector<int> & neighbor);
    bool is_inside(vec3 pos);
    
    vec3 get_displacement_projection(vec3 pos, vec3 norm, int phase, double max_displace, bool & bLast);
    bool is_inside(vec3 pos, int phase);
    
    void draw_connected_component();
    void compute_tranformation_mat_for_particle(int);
    
    double get_coeff(vec3 pos, int idx);
    
    std::vector<int> neighbor_search(vec3 pos, double radius);
    
    vec3 estimate_norm(vec3 pos);
public:
    // Debugging
    std::vector<mat3x3d> m_U;
    std::vector<mat3x3d> m_S;
    std::vector<mat3x3d> m_C;
};

#endif /* anisotrpic_kernel_h */
