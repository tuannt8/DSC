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

class hash3
{
public:
    hash3(vec3 domain_bound, double cell_size);
    ~hash3(){};
    
    void insert_point(vec3 pos, int index);
    std::vector<long> get_close_point(double x, double y, double z, double radius);
    
private:
    std::map<long, std::vector<long>> m_bins;
    vec3i m_dimension;
    double m_cell_size;
    
    inline long vec3i_to_index(vec3i const & v){
        return v[0]*m_dimension[1]*m_dimension[2] + v[1]*m_dimension[2] + v[2];
    };
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
    void build_hash();
    
    void load_time_step();
    void load(int idx, std::vector<particle> & par);
    
    void draw();
    
    vec3 get_displacement(vec3 pos);
};


#endif /* vtkWraper_hpp */
