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
        return m_sub_step_vel[idx];
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
        vec3 pos_key = m_sub_step_particles[key].pos;
        vec3 vel = m_sub_step_vel[key];
        auto cur_dis = (pos_key - pos).length();
        
        double q = cur_dis/h;
        double contribute = 21.0/16.0/3.14159/pow(h,3)*pow(1-q/2.0,4)*(1 +2*q);
        
        sum_vec += vel*(contribute * m_current_particles[key].mass/m_current_particles[key].density);
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

bool file_load::get_displacement_avg(vec3 pos, vec3 & dis)
{
    double r = get_influence_radius();
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

bool file_load::get_displacement(vec3 pos, vec3 & dis)
{
    
    return get_displacement_avg(pos, dis);
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

    build_hash();
    build_anisotropic_kernel();
}

bool file_load::get_projection(vec3 pos, vec3 direction, bool &bInside, double &t)
{
    double ra = get_spacing_distance();
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

void file_load::build_anisotropic_kernel()
{
    m_aniso_kernel.m_shared_tree = &m_vtree;
    m_aniso_kernel.m_shared_particles = m_sub_step_particles;
    m_aniso_kernel.m_h = get_influence_radius();
    m_aniso_kernel.m_ra = get_spacing_distance();
    
    m_aniso_kernel.build();
}

void file_load::draw()
{
//    personal_draw();
    
    

    if(glut_menu::get_state("particle points", 0))
    {
        glColor3d(1,0,0);
        glDisable(GL_LIGHTING);
        glPointSize(6);
        glBegin(GL_POINTS);
        //            int idx = 0;
        //            for (int idx : idx_list)
        for(int idx = 0; idx < m_sub_step_particles.size(); idx++)
        {
            auto &p = m_sub_step_particles[idx];
            
            static vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
            auto c = _color[p.type];
            glColor3f(c[0], c[1], c[2]);
            glVertex3dv(p.pos.get());
        }
        glEnd();
    }
    if(glut_menu::get_state("Principle component axis", 0))
    {
        std::vector<vec3> axis = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
        
        double ra = get_spacing_distance();
        //            int idx = 0;
        //            for (int idx : idx_list)
        for(int idx = 0; idx < m_sub_step_particles.size(); idx++)
        {
            auto &p = m_sub_step_particles[idx];
            
            auto pos = p.pos;
            
            glPushMatrix();
            glTranslated(pos[0], pos[1], pos[2]);
            
            auto G = m_aniso_kernel.m_G[idx];
            double m[16] = {G[0][0], G[0][1], G[0][2], 0,
                G[1][0], G[1][1], G[1][2], 0,
                G[2][0], G[2][1], G[2][2], 0,
                0, 0, 0, 1
            };
            glMultMatrixd(m);
            
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES);
            for (auto a : axis)
            {
                glColor3f(a[0], a[1], a[2]);
                glVertex3f(0, 0, 0);
                glVertex3dv((a*ra).get());
            }
            glEnd();
            
            
            glPopMatrix();
        }
    }
    
    if(glut_menu::get_state("Iso surface field", 1))
    {
        static vector<vec3> pos;
        static vector<double> phi;
        if (pos.size()==0)
        {
            int N = 200;
            vec3 delta = get_domain_dimension()/N;
            for (int i=0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
//                    for (int k = 0; k < N; k++)
                    int k = 80;
                    {
                        vec3 cur_p(i*delta[0], j*delta[1], k*delta[2]);
                        pos.push_back(cur_p);
                        phi.push_back(m_aniso_kernel.get_value(cur_p));
                    }
                }
            }
            
            // Log
            ofstream f("test.txt");;
            for (auto phi_i : phi)
            {
                f << phi_i << " ";
            }
            f.close();
            
            double max = *std::max_element(phi.begin(), phi.end());
            double min = *std::min_element(phi.begin(), phi.end());
            cout << "\nMax field: " << max << "; min: " << min << endl;
            for (auto & pp : phi)
            {
                if (pp < 0)
                {
                    pp /= abs(min);
                }else{
                    pp /= abs(max);
                }
            }
        }
        
        glDisable(GL_LIGHTING);
        glPointSize(4);
        glBegin(GL_POINTS);

        vec3 RED(1,0,0);
        vec3 BLUE(0,0,1);
        for (int i = 0; i < pos.size(); i++)
        {
            if(phi[i] == 0)
                continue;
            
            double dis_red = abs(phi[i] + 1)/2;
            double dis_blue = abs(1 - phi[i])/2;
            vec3 c = RED*dis_red + BLUE*dis_blue;
//            c.normalize();
            glColor3f(c[0], c[1], c[2]);
            glVertex3dv(pos[i].get());
        }
        glEnd();
    }
    
    if(glut_menu::get_state("Principle component", 1))
    {
        double ra = get_spacing_distance();
        //            int idx = 0;
        //            for (int idx : idx_list)
        for(int idx = 0; idx < m_current_particles.size(); idx++)
        {
            if (idx > 0)
            {
                break;
            }
            
            auto &p = m_current_particles[idx];
            auto pos = p.pos;
            
            glPushMatrix();
            glTranslated(pos[0], pos[1], pos[2]);
            
            //                auto G = m_aniso_kernel.m_principle[idx];
            auto G = m_aniso_kernel.m_G[idx]*0.2;
            double m[16] = {G[0][0], G[0][1], G[0][2], 0,
                G[1][0], G[1][1], G[1][2], 0,
                G[2][0], G[2][1], G[2][2], 0,
                0, 0, 0, 1
            };
            glMultMatrixd(m);
            
            glEnable(GL_LIGHTING);
            glEnable(GL_COLOR_MATERIAL);
            glColor3f(1, 0, 0);
            glutSolidSphere(ra, 10, 10);
            
            
            glPopMatrix();
        }
    }
    
#ifdef DEMO_INTERFACE
    glColor3f(0, 0, 1);
    m_interface.draw();
#endif
}

void file_load::build_hash()
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
