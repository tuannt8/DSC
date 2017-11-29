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

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellIterator.h>
#include <vtkPointData.h>

#define MAX_ITER 300

#include <map>

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

class vtkWrapper
{
public:
    vtkWrapper();
    ~vtkWrapper(){};

private:
    vtkSmartPointer<vtkUnstructuredGrid> m_current_grid;
    vtkSmartPointer<vtkUnstructuredGrid> m_next_grid;
    std::vector<vec3> m_current_pos;
    std::vector<vec3> m_next_pos;
    
    void convert(vtkUnstructuredGrid* grid, std::vector<vec3> &pos);
    
    int m_cur_idx = 0;
    std::shared_ptr<hash3> m_hashTable;
    
    void build_hash();
    vtkSmartPointer<vtkUnstructuredGrid> load_file(int idx);
    
    vec3 m_bound = vec3(0.4, 0.67, 1.6);
public:
    void draw();
    vec3 get_dispacement(vec3 pos);
    
    void load_next_grid();
    
    vec3 get_bound_size();
    vec3 get_bound_left_down();
    
    // DamBreak initial condition
    void dam_break_bound(vec3& left_down, vec3& right_up);
    
    void test();
};

#endif /* vtkWraper_hpp */
