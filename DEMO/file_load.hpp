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
    int object;
    vec3 vel, pos;
    
    void draw();
};

class file_load
{
public:
    file_load();
    ~file_load();
    
    std::vector<particle> m_current_particles;;
    std::vector<particle> m_next_particles;
    int m_cur_idx = 0;
    
    std::shared_ptr<hash3> m_hashTable;
    Geometry::KDTree<vec3, int> m_vtree;
    void build_hash();
    
    void load_time_step();
    void load(int idx, std::vector<particle> & par);
    
    void draw();
    
    vec3 get_displacement(vec3 pos);
    vec3 get_displacement_avg(vec3 pos);
    vec3 get_displacement_closet_point(vec3 pos);
    vec3 get_displacement_cubic_kernel(vec3 pos);
    
    virtual void init_dsc(DSC::DeformableSimplicialComplex<> * dsc){};
    virtual vec3 get_domain_dimension(){return vec3(0.0);};
    virtual double get_influence_radius(){return 0;};
};


#endif /* vtkWraper_hpp */
