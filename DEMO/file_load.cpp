//
//  vtkWraper.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "file_load.hpp"




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

void hash3::draw()
{
    
}

hash3::hash3(vec3 domain_bound, double cell_size)
{
    m_cell_size = cell_size;
    m_dimension = vec3i(ceil(domain_bound[0]/cell_size), ceil(domain_bound[1]/cell_size), ceil(domain_bound[2]/cell_size));
}

inline int hash3::get_idx_cell(vec3 & pos)
{
    return idx_int(get_idx_cell3(pos));
    
}

void hash3::insert_point(vec3 pos, int index)
{
    m_bins[get_idx_cell(pos)].push_back(index);
}

void fluid_interface::draw()
{
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < m_faces.size(); )
    {
        vector<vec3> p;
        p.push_back(m_points[m_faces[i++]]);
        p.push_back(m_points[m_faces[i++]]);
        p.push_back(m_points[m_faces[i++]]);
        
        auto norm = Util::normal_direction(p[0], p[1], p[2]);
        glNormal3dv(norm.get());
        
        for (auto & pt : p)
        {
            glVertex3dv(pt.get());
        }
        
    }
    glEnd();
}

std::vector<long> hash3::get_close_point(double x, double y, double z, double radius)
{
    vec3 ld(x-radius, y-radius, z-radius);
    vec3 ru(x+radius, y+radius, z+radius);
    
    vec3i ldi = get_idx_cell3(ld);
    vec3i rui = get_idx_cell3(ru);
    
    ldi = vec3i(std::max(ldi[0], 0), std::max(ldi[1], 0), std::max(ldi[2], 0));
    rui = vec3i(std::min(rui[0], m_dimension[0]), std::min(rui[1], m_dimension[1]), std::min(rui[2], m_dimension[2]));
    
    std::vector<long> out;
    for (int i = ldi[0]; i<=rui[0]; i++)
    {
        for (int j = ldi[1]; j <= rui[1]; j++)
        {
            for (int k = ldi[2]; k <= rui[2]; k++)
            {
                auto & list = m_bins[idx_int(vec3i(i,j,k))];
                out.insert(out.end(), list.begin(), list.end());
            }
        }
    }
    
    return  out;
}


//file_load::file_load()
//{
//    load_time_step();
//}

vec3 file_load::get_displacement_closet_point(vec3 pos)
{
    
    double r = get_influence_radius()*2;
    int idx;
    vec3 pp;
    if(m_vtree.closest_point(pos, r, pp, idx))
    {
        vec3 cur_pos = m_current_particles[idx].pos;
        vec3 nex_pos = m_next_particles[idx].pos;
        
        return nex_pos - cur_pos;
    }
    else
    {
        return vec3(0.0);
    }
}

vec3 file_load::get_displacement_WENLAND_kernel(vec3 pos)
{
    double h = get_influence_radius();
    double r = h*2;
    vector<int> pt_in_sphere;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
    
    vec3 sum_vec(0);
    for (auto key : pt_in_sphere)
    {
        vec3 cur_pos = m_current_particles[key].pos;
        vec3 nex_pos = m_next_particles[key].pos;
        auto cur_dis = (cur_pos - pos).length();
        
        double q = cur_dis/h;
        double contribute = 21.0/16.0/3.14159/pow(h,3)*pow(1-q/2.0,4)*(1 +2*q);
        
        sum_vec += (nex_pos - cur_pos)*(contribute * m_current_particles[key].mass/m_current_particles[key].density);
    }
    
    return sum_vec;
    
}

vec3 file_load::get_displacement_cubic_kernel(vec3 pos)
{
    double h = get_influence_radius();
    double r = h*2;
    vector<int> pt_in_sphere;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
    
    vec3 sum_vec(0);
    for (auto key : pt_in_sphere)
    {
        vec3 cur_pos = m_current_particles[key].pos;
        vec3 nex_pos = m_next_particles[key].pos;
        auto cur_dis = (cur_pos - pos).length();
        
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
        
        sum_vec += (nex_pos - cur_pos)*(contribute/3.1415/(h*h*h) * m_current_particles[key].mass/m_current_particles[key].density);
    }
    
//    sum_vec += closest_v*(1 - pow((pp-pos).length()/r, 3));
    
    return sum_vec;
}

vec3 file_load::get_displacement_avg(vec3 pos)
{
    double r = get_influence_radius();
    vector<int> list;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, r, pos_in_sphere, list);
    
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
    
    if (sum_dis < epsilon) // found nothing
    {
        //        assert(0);
        // Displace to closset point
        if(list.size() == 0)
        {
            // Tuan: Must be optimized later
            list.resize(m_current_particles.size());
            for (int i = 0; i < list.size(); i++)
            {
                list[i] = i;
            }
            
        }
        
        double min_dis = INFINITY;
        int min_pt = -1;
        for (auto p : list)
        {
            vec3 cur_pos = m_next_particles[p].pos;
            
            auto cur_dis = (cur_pos - pos).length() + epsilon;
            if (cur_dis < min_dis)
            {
                min_dis = cur_dis;
                min_pt = p;
            }
        }
        
        sum_vec = m_next_particles[min_pt].pos - pos;
        sum_dis = 1;
    }
    
    sum_vec /= sum_dis;
    
    return sum_vec;
}

vec3 file_load::get_displacement(vec3 pos)
{
    return get_displacement_avg(pos);
//    return get_displacement_WENLAND_kernel(pos);
//    return get_displacement_cubic_kernel(pos);
}
file_load::~file_load()
{
    
}

void file_load::load_time_step()
{
#ifdef DEMO_INTERFACE
    m_interface.load_surface(m_cur_idx);
#endif
    
    if(m_cur_idx==0)
        load(0, m_current_particles);
    else
        m_current_particles = m_next_particles;
    
    load(++m_cur_idx, m_next_particles);
    
    build_hash();
    build_anisotropic_kernel();
}

bool file_load::get_projection(vec3 pos, vec3 direction, bool &bInside, double &t)
{
    double max_search = get_influence_radius()*2;
    double phi0 = m_aniso_kernel.get_value(pos);
    
    static double epsilon = 1e-8;
    bInside = phi0 > epsilon;

    if (phi0 < epsilon) // outside
    {
        t = -max_search;
    }
    else{
        t = max_search;
    }
    
    double phi_t0 = m_aniso_kernel.get_value(pos + direction*t);
    bool btInside0 = phi_t0 > epsilon;
    if (!(btInside0 ^ bInside)) // both inside or outside
    {
        return false;
    }
    
    double t2 = 0;
    for (int iter = 0; iter < 3; iter++)
    {
        double phi_t = m_aniso_kernel.get_value(pos + direction*t);
        bool btInside = phi_t > epsilon;
        double phi_t2 = m_aniso_kernel.get_value(pos + direction*t2);
        bool btInside2 = phi_t2 > epsilon;
        
        if (btInside2 ^ bInside) // t2 and pos
        {
            t = t2 + (t-t2)/2;
        }
        else
        {
            t2 = t2 + (t-t2)/2;
        }
    }
    
    t = (t+t2)/2;
    
    return true;
}

void file_load::build_anisotropic_kernel()
{
    m_aniso_kernel.m_shared_tree = &m_vtree;
    m_aniso_kernel.m_shared_particles = m_current_particles;
    m_aniso_kernel.m_h = get_influence_radius();
    
    m_aniso_kernel.build();
}

void file_load::draw()
{
    personal_draw();

#ifdef DEMO_INTERFACE
    glColor3f(0, 0, 1);
    m_interface.draw();
#endif
}

void file_load::build_hash()
{
    m_vtree = Geometry::KDTree<vec3, int>();
    
    int idx = 0;
    for (auto &p : m_current_particles)
    {
        if (p.type == 0)
        {
            m_vtree.insert(p.pos, idx);
        }
        idx++;
    }
    
    m_vtree.build();
    
}

void file_load::load(int idx, std::vector<particle> & par)
{
    stringstream s;
    s << m_data_path << setfill('0') << setw(5) << idx << ".particle";

    
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
        cout << "Error: " <<  s.str() << endl;
    }
}
