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

class vtkWrapper
{
public:
    vtkWrapper();
    ~vtkWrapper(){};

private:
    vtkSmartPointer<vtkUnstructuredGrid> m_current_grid;
    vtkSmartPointer<vtkUnstructuredGrid> m_next_grid;
    
    int m_cur_idx = 0;
    
    vtkSmartPointer<vtkUnstructuredGrid> load_file(int idx);
public:
    void draw();
    vec3 get_displace_ment(vec3 pos);
    
    void load_next_grid();
    
    vec3 get_bound_size();
    vec3 get_bound_left_down();
    
    // DamBreak initial condition
    void dam_break_bound(vec3& left_down, vec3& right_up);
};

#endif /* vtkWraper_hpp */
