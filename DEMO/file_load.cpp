//
//  vtkWraper.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "file_load.hpp"

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

using namespace std;

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


file_load::file_load()
{
    load_time_step();
}

vec3 file_load::get_displacement(vec3 pos)
{
    double r = 0.052;
    auto list = m_hashTable->get_close_point(pos[0], pos[1], pos[2], r);
    
    vec3 sum_vec(0.0);
    double sum_dis = 0;
    double epsilon = 1e-8;
    for (auto p : list)
    {
        vec3 cur_pos = m_current_particles[p].pos;
        vec3 nex_pos = m_next_particles[p].pos;
        
        auto cur_dis = (cur_pos - pos).length() + epsilon;
        if (cur_dis < r)
        {
            sum_vec += (nex_pos - cur_pos)*cur_dis;
            sum_dis += cur_dis;
        }
    }
    
    if (sum_dis > epsilon)
    {
        sum_vec /= sum_dis;
    }
    
    return sum_vec;
}
file_load::~file_load()
{
    
}

void file_load::load_time_step()
{
    if(m_cur_idx==0)
        load(0, m_current_particles);
    else
        m_current_particles = m_next_particles;
    
    load(++m_cur_idx, m_next_particles);
    
    build_hash();
}

inline std::istream& operator>> (std::istream&is, particle& p)
{
    is >> p.pos[0] >> p.pos[1] >> p.pos[2]
    >> p.pressure
    >> p.density
    >> p.mass
    >> p.type
    >> p.flag
    >> p.object
    >> p.vel[0] >> p.vel[1] >> p.vel[2];
    
    return is;
}

void particle::draw()
{
    static vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
    
    glBegin(GL_POINTS);
    glColor3dv(_color[type].get());
    glVertex3dv(pos.get());
    glEnd();
}

void file_load::draw()
{
    glPointSize(2.5);
    for (auto &p : m_current_particles)
    {
        p.draw();
    }
}

void file_load::build_hash()
{
    m_hashTable = shared_ptr<hash3>( new hash3(vec3(0.4, 0.67, 1.6), 0.1));
    int idx = 0;
    for (auto &p : m_current_particles)
    {
        if (p.type == 0)
        {
            m_hashTable->insert_point(p.pos, idx++);
        }
    }
}

void file_load::load(int idx, std::vector<particle> & par)
{
    try
    {
        stringstream s;
#if defined(__APPLE__)
        s << "../Large_data/DamBreak3D/my_format/iter_" << setfill('0') << setw(5) << idx << ".particle";
#else
        s << "../../Large_data/DamBreak3D/my_format/iter_" << setfill('0') << setw(5) << idx << ".particle";
#endif
        
        std::ifstream f(s.str());
        if(f.is_open())
        {
            int num_points;
            f >> num_points;
            string comment;
            getline(f, comment);
            
            par.resize(num_points);
            for (int i = 0; i < num_points; i++)
            {
                f >> par[i];
            }
        }
        else{
        	cout << s.str();
            throw "Fail to load particle";
        }
    }
    catch (exception e)
    {
        cout << "Error fail to load file" <<  e.what();
    }

    
}
