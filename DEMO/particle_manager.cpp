//
//  particle.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#include "particle_manager.hpp"

#include <unordered_map>
#include <cassert>

using namespace std;

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

inline std::ostream& operator<<(std::ostream&os, particle& p)
{
    os << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
    << p.pressure << " "
    << p.density << " " << endl;
    
    return os;
}


vec3 particle_manager::get_displacement_closet_point(vec3 pos)
{
    double r = m_influence_radius;
    int idx;
    vec3 pp;
    if(m_vtree.closest_point(pos, r, pp, idx))
    {
        return m_sub_step_vel[idx];
    }
    else
    {
        return vec3(0.0);
    }
}

bool particle_manager::get_displacement_WENLAND_kernel(vec3 pos, vec3 &dis)
{
    double h = m_influence_radius;
    double r = h;
    vector<int> pt_in_sphere;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
    
    if(pt_in_sphere.size()==0)
        return false;
    
    vec3 sum_vec(0);
    for (auto key : pt_in_sphere)
    {
        vec3 pos_key = m_sub_step_particles[key].pos;
        vec3 vel = m_sub_step_vel[key];
        auto cur_dis = (pos_key - pos).length();
        
        double q = cur_dis/h;
        double contribute = 21.0/16.0/3.14159/pow(h,3)*pow(1-q/2.0,4)*(1 +2*q);
        
        sum_vec += vel*(contribute * m_current_particles[key].mass/m_current_particles[key].density);
    }
    
    dis = sum_vec;
    return true;
    
}

vec3 particle_manager::get_displacement_cubic_kernel(vec3 pos)
{
    double h = m_influence_radius;
    double r = h*2;
    vector<int> pt_in_sphere;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
    
    vec3 sum_vec(0);
    for (auto key : pt_in_sphere)
    {
        vec3 pos_key = m_sub_step_particles[key].pos;
        vec3 vel = m_sub_step_vel[key];
        auto cur_dis = (pos_key - pos).length();
        
        double q = cur_dis/h;
        double contribute;
        if (cur_dis < 1)
        {
            contribute = 1 - 1.5*q*q + 0.75*q*q*q;
        }
        else
        {
            contribute = 0.25*pow(2-q, 3);
        }
        sum_vec += vel*(contribute/3.1415/(h*h*h) * m_current_particles[key].mass/m_current_particles[key].density);
    }
    
    return sum_vec;
}

bool particle_manager::get_displacement_avg(vec3 pos, vec3 & dis)
{
    double r = m_influence_radius;
    vector<int> list;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, list);
    
    if (list.size() == 0) // found nothing
    {
        return false;
    }
    
    vec3 sum_vec(0.0);
    double sum_dis = 0;
    for (auto p : list)
    {
        vec3 pos_key = m_sub_step_particles[p].pos;
        vec3 vel = m_sub_step_vel[p];
        auto cur_dis = (pos_key - pos).length();
        
        sum_vec += vel*cur_dis;
        sum_dis += cur_dis;
    }
    
    sum_vec /= sum_dis;
    
    dis = sum_vec;
    return true;
}

void particle_manager::draw()
{
    glDisable(GL_LIGHTING);
    glPointSize(6);
    glBegin(GL_POINTS);
    //            int idx = 0;
    //            for (int idx : idx_list)
    for(int idx = 0; idx < m_sub_step_particles.size(); idx++)
    {
        auto &p = m_sub_step_particles[idx];
        
        if (p.type != 0)
        {
            continue;
        }
        
//        static vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
//        auto c = _color[p.type];
//        glColor3f(c[0], c[1], c[2]);
        glVertex3dv(p.pos.get());
    }
    glEnd();
}

bool particle_manager::get_displacement(vec3 pos, vec3 & dis)
{
    
    //    return get_displacement_avg(pos, dis);
    return get_displacement_WENLAND_kernel(pos, dis);
    //    return get_displacement_cubic_kernel(pos);
}

void particle_manager::load_time_step(int idx)
{
    if (idx == -1)
    {
        load(0, m_next_particles);
    }
    else
    {
        m_current_particles = m_next_particles;
        load(idx+1, m_next_particles);
    }
}
void particle_manager::interpolate(int sub_idx, int sub_count)
{
    // Load sub step by linear interpolation
    m_sub_step_particles = m_current_particles;

    for (int i = 0; i<m_current_particles.size(); i++)
    {
        auto pre_pos = m_current_particles[i].pos;
        auto next_pos = m_next_particles[i].pos;
        
        m_sub_step_particles[i].pos = pre_pos + (next_pos - pre_pos)*(sub_idx/(double)sub_count);
    }
    
    if(sub_idx == 0)//velocity is unchanged hence only need one computation
    {
        m_sub_step_vel.resize(m_current_particles.size());
        for (int i = 0; i<m_current_particles.size(); i++)
        {
            auto pre_pos = m_current_particles[i].pos;
            auto next_pos = m_next_particles[i].pos;
            
            m_sub_step_vel[i] = (next_pos - pre_pos)/(double)sub_count;
        }
    }
}

double particle_manager::get_max_displacement()
{
    double max_displacement = -INFINITY;
    for (int i = 0; i<m_current_particles.size(); i++)
    {
        auto pre_pos = m_current_particles[i].pos;
        auto next_pos = m_next_particles[i].pos;
        
        max_displacement = std::max(max_displacement, (next_pos - pre_pos).length());
    }
    
    return max_displacement;
}

void particle_manager::load_time_step()
{
    static bool first_load = true;
    if (first_load)
    {
        first_load = false;
        m_cur_idx = 0;
        m_cur_sub_step = 0;
        load(0, m_current_particles);
        load(++m_cur_idx, m_next_particles);
    }
    else
    {
        if(m_cur_sub_step == m_max_step)
        {
            m_current_particles = m_next_particles;
            load(++m_cur_idx, m_next_particles);
            m_cur_sub_step = 0;
        }
    }
    
    m_sub_step_particles = m_current_particles;
    m_sub_step_vel.resize(m_current_particles.size());
    for (int i = 0; i<m_current_particles.size(); i++)
    {
        auto pre_pos = m_current_particles[i].pos;
        auto next_pos = m_next_particles[i].pos;
        
        m_sub_step_vel[i] = (next_pos - pre_pos)/m_max_step;
        m_sub_step_particles[i].pos = pre_pos + (next_pos - pre_pos)*(m_cur_sub_step/(double)m_max_step);
    }
    
    m_cur_sub_step++;
    
    build_kd_tree(); // should not do it every iteration
}

bool particle_manager::get_projection(vec3 pos, vec3 direction, bool &bInside, double &t)
{
    double ra = m_influence_radius;
    double phi_pre = m_aniso_kernel.get_value(pos);
    double epsilon = -ra*(phi_pre>0? -1:1);
    
    bInside = phi_pre > 0;
    
    double phi_epsilon = m_aniso_kernel.get_value(pos + direction*epsilon);
    
    if (phi_pre * phi_epsilon < 0) //the point is somewhere in between
    {
        double ep1 = epsilon;
        double ep2 = 0;
        double phi1 = phi_epsilon;
        double phi2 = phi_pre;
        
        for (int i = 0; i < 4; i++)
        {
            double phi_middle = m_aniso_kernel.get_value(pos + direction*(ep1+ep2)/2);
            if (std::abs(phi_middle) < 0.0001)
            {
                break;
            }
            
            if (phi_middle*phi1 < 0)
            {
                ep2 = (ep1 + ep2)/2;
            }
            else
            {
                ep1 = (ep1 + ep2)/2;
                phi1 = phi_middle;
            }
            
            
        }
        
        t = (ep1 + ep2)/2;
        
        return true;
    }
    
    //    t = epsilon;
    // fail to project
    return false;
}

void particle_manager::build_anisotropic_kernel()
{
    m_aniso_kernel.m_shared_tree = &m_vtree;
    m_aniso_kernel.m_shared_particles = m_sub_step_particles;
    m_aniso_kernel.m_h = m_influence_radius;
    m_aniso_kernel.m_ra = m_deltap;
    
    m_aniso_kernel.build();
}

void particle_manager::build_kd_tree()
{
    m_vtree = Geometry::KDTree<vec3, int>();
    
    int idx = 0;
    for (auto &p : m_sub_step_particles)
    {
        if (p.type == 0)
        {
            m_vtree.insert(p.pos, idx);
        }
        idx++;
    }
    
    m_vtree.build();
    
}

void particle_manager::load(int idx, std::vector<particle> & par)
{
    stringstream s;
    s << m_data_path << "/iter_" << setfill('0') << setw(5) << idx << ".particle";
    
    
    std::ifstream f(s.str());
    if(f.is_open())
    {
        int num_points;
        f >> num_points;
        
        par.resize(num_points);
        for (int i = 0; i < num_points; i++)
        {
            f >> par[i];
        }
    }
    else{
        cout << "No more file to load. Last file: " <<  s.str() << endl;
        exit(EXIT_SUCCESS);
    }
}
