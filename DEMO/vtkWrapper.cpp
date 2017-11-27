//
//  vtkWraper.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "vtkWrapper.hpp"

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

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

vtkWrapper::vtkWrapper()
{
    m_current_grid = load_file(0);
    m_next_grid = load_file(1);
    m_cur_idx = 1;
}

void vtkWrapper::load_next_grid()
{
    if(m_cur_idx < MAX_ITER - 1)
    {
        m_current_grid = m_next_grid;
        m_next_grid = load_file(++m_cur_idx);
    }
}

vec3 vtkWrapper::get_bound_size()
{
    double *bound = m_current_grid->GetBounds();
    return vec3(bound[1]-bound[0], bound[3] - bound[2], bound[5] - bound[4]);
}

vtkSmartPointer<vtkUnstructuredGrid> vtkWrapper::load_file(int idx)
{
    std::stringstream s;
    s << "../Large_data/DamBreak3D/data/PART_" << setfill('0') << setw(5) << idx << ".vtu";

    
    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(s.str().c_str());
    reader->Update();
    
    return vtkSmartPointer<vtkUnstructuredGrid>(reader->GetOutput());
}

void vtkWrapper::draw()
{
    if (m_current_grid)
    {
        glBegin(GL_POINTS);
        
        auto pos = m_current_grid->GetPointData()->GetVectors();
        vtkPointData * trait = m_current_grid->GetPointData();
        
        std::vector<vec3> _color ={vec3(1,0,0), vec3(0,1,0), vec3(0,0,1), vec3(1,1,0)};
        
        for (int i = 0; i < m_current_grid->GetNumberOfPoints(); i++)
        {
            auto id = trait->GetArray(3)->GetTuple(i);
            
            glColor3dv(_color[(int)(*id)].get());
            
            double * pi = m_current_grid->GetPoint(i);
            glVertex3dv(pi);
        }
        
        glEnd();
    }
}
