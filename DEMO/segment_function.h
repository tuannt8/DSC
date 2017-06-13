//
//  segment_function.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef segment_function_hpp
#define segment_function_hpp

#include <stdio.h>
#include "image3d.h"
#include "DSC.h"

#include "define.h"

extern std::bitset<4> X_direction, Y_direction, Z_direction;


struct intersect_pt
{
    intersect_pt(int z_, bool b_in_z_){
        z = z_;
        b_in = b_in_z_;
    }
    int z;
    bool b_in;
};

struct ray_z
{
    int x, y;
    std::vector<intersect_pt> intersects;
};

class segment_function
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
    
public:
    segment_function(){};
    ~segment_function(){};
    
    void init(); // Load image
    void initialze_segmentation(); // Label initialization
    void random_initialization(); // Label initialization
    void initialization_discrete_opt(); // Optimize the initialization automatically
    
    void segment();
    
public:// Configuration parametters
    int NB_PHASE;
    double VARIATION_THRES_FOR_RELABELING;
    double ALPHA;
    double QALPHA = 0.2;
    double _dt = 1;
    std::string _directory_path;
    
    // Face plit
    double ratio_signed_and_mag_mean = 0.2;
    
public:// Variables
    image3d _img; // Store crossection -> voxel
    dsc_class *_dsc; // Shared dsc
    
public:
    std::vector<double> _mean_intensities;
    std::vector<double> _total_intensities; // To update mean intensity during relabeling
    std::vector<double> _phase_volume;
    void update_average_intensity();
    void update_average_intensity1();
    
    std::vector<int> _vertex_stability_map;
    std::vector<vec3> _forces;
    std::vector<vec3> _internal_forces;
    std::vector<vec3> _quality_control_forces; // curvature based force
    std::vector<vec3> _quality_angle_forces; // angle based
    
    void compute_surface_curvature();
    void compute_external_force();
    void compute_internal_force();
    void compute_mesh_quality_control_force();
    
    void work_around_on_boundary_vertices();
    void update_vertex_stability();
    
    
    double get_energy_tetrahedron(is_mesh::TetrahedronKey tkey, int );
    void relabel_tetrahedra();
    
    void face_split();
    void adapt_tetrahedra();
public:
    std::vector<ray_z> _d_rayz;
    
    std::vector<std::vector<vec3>> _mean_curvature_of_each_hat;
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::TetrahedronKey>>> _tets_in_hat;
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::NodeKey>>> _node_in_hat;
    
    // For debuging
    std::vector<vec3> boundary_vertices_displacements;
    
    std::vector<unsigned int> d_is_image_boundary;
    std::vector<std::bitset<4>> d_direction_state;
    
    std::vector<std::vector<vec3>> _curvature_force;
    std::vector<std::vector<vec3>> _area_force;
};

#endif /* segment_function_hpp */
