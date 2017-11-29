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

#include <unordered_map>
#include <cassert>


hash3::hash3(vec3 domain_bound, double cell_size)
{
    m_cell_size = cell_size;
    m_dimension = vec3i(domain_bound[0]/cell_size + 1, domain_bound[1]/cell_size + 1, domain_bound[2]/cell_size + 1);
}

void hash3::insert_point(vec3 pos, int index)
{
    vec3i idx_h(pos[0]/m_cell_size, pos[1]/m_cell_size, pos[2]/m_cell_size);
    
    m_bins[vec3i_to_index(idx_h)].push_back(index);
}

std::vector<long> hash3::get_close_point(double x, double y, double z, double radius)
{
    vec3 ld(x-radius, y-radius, z-radius);
    vec3 ru(x+radius, y+radius, z+radius);
    
    vec3i ldi(ld/m_cell_size);
    vec3i rui(ru/m_cell_size);
    std::vector<long> out;
    for (int i = ldi[0]; i<=rui[0]; i++)
    {
        for (int j = ldi[1]; j <= rui[1]; j++)
        {
            for (int k = ldi[2]; k <= rui[2]; k++)
            {
                auto & list = m_bins[vec3i_to_index(vec3i(i,j,k))];
                out.insert(out.end(), list.begin(), list.end());
            }
        }
    }
    
    return  out;
}


void vtkWrapper::build_hash()
{
    convert(m_current_grid, m_current_pos);
    convert(m_next_grid, m_next_pos);
    
    double cell_size = 0.3;
    m_hashTable = std::shared_ptr<hash3>(new hash3(m_bound, cell_size));
    
    for (int i = 0; i < m_current_pos.size(); i++)
    {
        m_hashTable->insert_point(m_current_pos[i], i);
    }
}
void vtkWrapper::test()
{
    std::stringstream s;
    s << "../Large_data/DamBreak3D/data/PART_" << setfill('0') << setw(5) << 0 << ".vtu";
    
    
    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(s.str().c_str());
    reader->Update();
    
    std::cout << reader->GetNumberOfPoints() << " particles\n";
    // Print info
    for (int i =0; i < reader->GetNumberOfPointArrays(); i++)
    {
        std::cout << reader->GetPointArrayName(i) << "; " ;
    }
    std::cout<<std::endl;
    
}

vec3 vtkWrapper::get_dispacement(vec3 pos)
{
    double radius = 0.052;
    auto list = m_hashTable->get_close_point(pos[0],pos[1], pos[2], radius);
    
    
    double denum = 0;
    double epsilon = 0.000001;
    vec3 dis_v(0.0);
    
    double max_particle_move = -INFINITY;
    
    for (auto idx : list)
    {
        vec3 cur_pos = m_current_pos[idx];
        vec3 next_pos = m_next_pos[idx];
        
        auto dis = (pos - cur_pos).length() + epsilon;
        
        if(dis < radius)
        {
            dis_v += (next_pos - cur_pos)*dis;
            denum += dis;
            
//            max_particle_move = std::max( (next_pos - cur_pos).length(), max_particle_move);
        }
    }
    
    if (denum > epsilon)
    {
        dis_v /= denum;
    }
    
//    std::cout << "Particle max: " << max_particle_move << "; vertex dis: " << dis_v.length() << std::endl;
    
    return dis_v;
}

vtkWrapper::vtkWrapper()
{
    test();
    
    m_current_grid = load_file(0);
    m_next_grid = load_file(1);
    m_cur_idx = 1;
    
    build_hash();
}

void vtkWrapper::load_next_grid()
{
    if(m_cur_idx < MAX_ITER - 1)
    {
        m_current_grid = m_next_grid;
        m_next_grid = load_file(++m_cur_idx);
        
        build_hash();
    }
}

void vtkWrapper::convert(vtkUnstructuredGrid* grid, std::vector<vec3> & pos)
{
    pos.resize(grid->GetNumberOfPoints());
    
    for (int i = 0; i < grid->GetNumberOfPoints(); i++)
    {
        auto p = grid->GetPoint(i);
        auto idx = grid->GetPointData()->GetArray(6)->GetTuple(i);
        
        pos[*idx] = vec3(p[0], p[1], p[2]);
    }
}
void vtkWrapper::dam_break_bound(vec3& left_down, vec3& right_up)
{
    
}

vec3 vtkWrapper::get_bound_left_down()
{
    double *bound = m_current_grid->GetBounds();
    return vec3(bound[0], bound[2], bound[4]);
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
        glPointSize(2.0);
        
        glBegin(GL_POINTS);
        
        vtkPointData * trait = m_current_grid->GetPointData();
        
        std::vector<vec3> _color ={vec3(1,0,0), vec3(0,1,0), vec3(0,0,1), vec3(1,1,0)};
        
        for (int i = 0; i < m_current_grid->GetNumberOfPoints(); i++)
        {
            
            auto id = trait->GetArray(3)->GetTuple(i);
            
            glColor3dv(_color[(int)(*id)].get());
            
//            if (*(trait->GetArray(6)->GetTuple(i)) > 10000
//                && *(trait->GetArray(6)->GetTuple(i)) < 11000)
//            {
//                glColor3d(0, 0, 1);
//            }
//            
            double * pi = m_current_grid->GetPoint(i);
            glVertex3dv(pi);
        }
        
        glEnd();
    }
}
