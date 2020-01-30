//
//  particle.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#include "particle_manager.hpp"

#include "eigen_wrapper.hpp"

//#include <unordered_map>
//#include <cassert>

#include "glut_menu.hpp"


#include "debugger.h"

double wendland(double q, double h)
{
    static double alpha = 21./(16*3.14159);
    return alpha/(h*h*h) * std::pow(1 - q/2, 4) * (1 + 2*q);
}

using namespace std;

inline std::istream& operator>> (std::istream&is, particle& p)
{
    is >> p.pos[0] >> p.pos[1] >> p.pos[2]
//    >> p.pressure
    >> p.density
    >> p.mass
    >> p.fluid;
//    >> p.vel[0] >> p.vel[1] >> p.vel[2];
    
    return is;
}

inline std::ostream& operator<<(std::ostream&os, particle& p)
{
    os << p.pos[0] << " " << p.pos[1] << " " << p.pos[2] << " "
    << p.mass << " "
    << p.density << " " << endl;
    
    return os;
}


vec3 particle_manager::get_displacement_closet_point(vec3 pos)
{
//    double r = m_influence_radius;
//    int idx;
//    vec3 pp;
//    if(m_vtree.closest_point(pos, r, pp, idx))
//    {
//        return m_sub_step_vel[idx];
//    }
//    else
//    {
//        return vec3(0.0);
//    }
}

double particle_manager::weight_function(double r,  int type)
{
    if (type == POLY_4)
    {
        return  1 - 6*r*r + 8*r*r*r - 3*r*r*r*r;
    }
    else{
        throw std::invalid_argument("Weight function type unknown");
    }
}
bool particle_manager::get_displacement_MLS_kernel(vec3 pos, vec3 & dis)
{
//    double support_radius = m_deltap*3;
//    
//    vector<int> pt_in_sphere;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, support_radius, pos_in_sphere, pt_in_sphere);
//    
//    if(pt_in_sphere.size() <= 3) //ignore if have only few neighbors
//        return false;
//    
//    // 1. Compute moment matrix
//    mat4x4d A(0);
//    for (auto j : pt_in_sphere)
//    {
//        auto p_n = m_current_particles[j].pos;
//        vec4 poly_n(1, p_n[0], p_n[1], p_n[2]); // Polynominal function
//
//        mat4x4d pxpt;
//        CGLA::outer_product(poly_n, poly_n, pxpt);
//
//        double r = (pos - p_n).length() / support_radius;
//            // As we search for points in sphere then r < 1
//        A += pxpt * weight_function(r);
//    }
//
//    // 1.1. Invert of A
//    mat4x4d A_i;
//    invert4x4(A.get(), A_i.get());
//    
//    // 2. Compute shape functions
//    vector<double> phi(pt_in_sphere.size(), 0);
//    
//    dis = vec3(0);
//    
//    double sum = 0;
//    for (int j = 0; j < pt_in_sphere.size(); j++)
//    {
//        int par_idx = pt_in_sphere[j];
//        auto p_n = m_current_particles[par_idx].pos;
//        auto v_n = m_sub_step_vel[par_idx];
//        vec4 poly(1, p_n[0], p_n[1], p_n[2]);
//        double r = (pos - p_n).length() / support_radius;
//        
//        auto mat_o = (A_i*poly);
//        phi[j] = CGLA::dot(poly, mat_o) *weight_function(r);
//
//        sum += phi[j];
//        dis += v_n*phi[j];
//        
//        
////        assert(v_n.length() < support_radius);
//    }
//    
//    assert(dis.length() < support_radius);
    
//    cout << sum << " ; ";
    
    return true;
}

//bool particle_manager::get_displacement_MLS_kernel(vec3 pos, vec3 & dis)
//{
////    using namespace Eigen;
//
//    double support_radius = m_slength*2;
//
//    vector<int> pt_in_sphere;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, support_radius, pos_in_sphere, pt_in_sphere);
//
//    if(pt_in_sphere.size() < 1) //ignore if have only few neighbors
//        return false;
//
////    // 1. Compute moment matrix
////    Matrix4d A = Matrix4d::Zero();
////    for (auto j : pt_in_sphere)
////    {
////        auto p_n = m_current_particles[j].pos;
////        Vector4d poly_n(1, p_n[0], p_n[1], p_n[2]); // Polynominal function
////
////        Matrix4d pxpt = poly_n * poly_n.transpose();
////
////        double r = (pos - p_n).length() / support_radius;
////            // As we search for points in sphere then r < 1
////        A += pxpt * weight_function(r);
////    }
////
////    // 1.1. Invert of A
////    auto A_i = A.inverse();
//
//    // 2. Compute shape functions
//    vector<double> phi(pt_in_sphere.size(), 0);
//
//    dis = vec3(0);
//
//    for (int j = 0; j < pt_in_sphere.size(); j++)
//    {
////        int par_idx = pt_in_sphere[j];
////        auto p_n = m_current_particles[par_idx].pos;
////        auto v_n = m_current_particles[par_idx].vel;
////        Vector4d poly(1, p_n[0], p_n[1], p_n[2]);
////        double r = (pos - p_n).length() / support_radius;
////        auto mat_o = poly.transpose()*(A_i*poly);
////        phi[j] = mat_o(0,0)*weight_function(r);
//
////        dis += v_n*phi[j];
//    }
//    return true;
//}

bool particle_manager::get_displacement_WENLAND_kernel(vec3 pos, vec3 &dis)
{
//    double h = m_slength;
//    double r = 2*h;
//    vector<int> pt_in_sphere;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
//
//    if(pt_in_sphere.size()==0)
//        return false;
//
//    vec3 sum_vec(0);
//    for (auto key : pt_in_sphere)
//    {
//        vec3 pos_key = m_sub_step_particles[key].pos;
//        vec3 vel = m_sub_step_vel[key];
//        auto cur_dis = (pos_key - pos).length();
//
//        double q = cur_dis/h;
//        double contribute = 21.0/16.0/3.14159/pow(h,3)*pow(1-q/2.0,4)*(1 +2*q);
//
//        sum_vec += vel*(contribute * m_current_particles[key].mass/m_current_particles[key].density);
//    }
//
//    dis = sum_vec;
//    return true;
    
}

vec3 particle_manager::get_displacement_cubic_kernel(vec3 pos)
{
//    double h = m_influence_radius;
//    double r = h*2;
//    vector<int> pt_in_sphere;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, r, pos_in_sphere, pt_in_sphere);
//
//    vec3 sum_vec(0);
//    for (auto key : pt_in_sphere)
//    {
//        vec3 pos_key = m_sub_step_particles[key].pos;
//        vec3 vel = m_sub_step_vel[key];
//        auto cur_dis = (pos_key - pos).length();
//
//        double q = cur_dis/h;
//        double contribute;
//        if (cur_dis < 1)
//        {
//            contribute = 1 - 1.5*q*q + 0.75*q*q*q;
//        }
//        else
//        {
//            contribute = 0.25*pow(2-q, 3);
//        }
//        sum_vec += vel*(contribute/3.1415/(h*h*h) * m_current_particles[key].mass/m_current_particles[key].density);
//    }
//
//    return sum_vec;
}

bool particle_manager::get_displacement_sph_kernel(vec3 pos, vec3 & dis)
{
    double h = m_slength;
    
    vector<int> list;
    vector<vec3> pos_in_sphere;
    m_vtree.in_sphere(pos, 2*h, pos_in_sphere, list);
    
    if (list.size() == 0) // found nothing. Cannot project
    {
        dis = vec3(0.0);
        return false;
    }

    dis = vec3(0);
    double sum = 0;
    for (auto p : list)
    {
        auto part = m_current_particles[p];
        auto part_next = m_next_particles[p];
//        auto vel = m_current_particles[p].vel;
        
        
        auto pdis = part_next.pos - part.pos;
        
        double r_h = (part.pos - pos).length()/h;
        
        double coeff = wendland(r_h, h);
        
        dis += pdis * coeff;
        sum += coeff;
        
//        dis += vel * (
//                      part.mass / part.density
//                      * wendland(r_h, h)
//                      );
        
        
    }
    
    if(sum > 0)
        dis /= sum;
    
    return  true;
}

bool particle_manager::get_displacement_weighted_avg(vec3 pos, vec3 & dis)
{
    assert(0);
//    double r = m_slength;
//    vector<int> list;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, r, pos_in_sphere, list);
//
//    if (list.size() == 0) // found nothing
//    {
//        return false;
//    }
//
//    vec3 sum_vec(0.0);
//    double sum_dis = 0;
//    for (auto p : list)
//    {
//        vec3 pos_key = m_sub_step_particles[p].pos;
//        vec3 vel = m_sub_step_vel[p];
//        auto cur_dis = (pos_key - pos).length();
//
//
//        auto weight = weight_function(cur_dis/r);
//
//        sum_vec += vel*weight;
//        sum_dis += weight;
//    }
//
//    sum_vec /= sum_dis;
//
//    dis = sum_vec;
//    return true;
}

bool particle_manager::get_displacement_avg(vec3 pos, vec3 & dis)
{
    assert(0);
//    double r = m_deltap*3;
//    vector<int> list;
//    vector<vec3> pos_in_sphere;
//    m_vtree.in_sphere(pos, r, pos_in_sphere, list);
//
//    if (list.size() == 0) // found nothing
//    {
//        return false;
//    }
//
//    vec3 sum_vec(0.0);
//    double sum_dis = 0;
//    for (auto p : list)
//    {
//        vec3 pos_key = m_sub_step_particles[p].pos;
//        vec3 vel = m_sub_step_vel[p];
//        auto cur_dis = (pos_key - pos).length();
//
//        sum_vec += vel*cur_dis;
//        sum_dis += cur_dis;
//    }
//
//    sum_vec /= sum_dis;
//
//    dis = sum_vec;
//    return true;
}

void particle_manager::draw_intermediate_vel()
{
//    glDisable(GL_LIGHTING);
//    glBegin(GL_LINES);
//    for (int i = 0; i < m_sub_step_vel.size(); i++)
//    {
//        auto p = m_sub_step_particles[i].pos;
//        glVertex3dv(p.get());
//        glVertex3dv((p + m_sub_step_vel[i]).get());
//    }
//    glEnd();
}


void particle_manager::draw_anisotropic_kernel(double yunder, double yupper)
{
    for (int i = 0; i < m_current_particles.size(); i++)
    {
//        {
//            int i = debugger<>::get_int("particle idx", 0);
            auto &p = m_current_particles[i];
            auto pos = p.pos;
            
            if (pos[1] < yunder || pos[1] > yupper)
            {
                continue;
            }
            
            // Draw kernel
            if(glut_menu::get_state("isotropic sphere", 1))
            {
                auto G = m_aniso_kernel.get_transform_mat(i)*m_slength;
                G = m_aniso_kernel.m_C[i];
                
                glPushMatrix();
                glTranslated(pos[0], pos[1], pos[2]);
                
                double m[16] = {G[0][0], G[0][1], G[0][2], 0,
                    G[1][0], G[1][1], G[1][2], 0,
                    G[2][0], G[2][1], G[2][2], 0,
                    0, 0, 0, 1
                };
                glMultMatrixd(m);
                
                glEnable(GL_LIGHTING);
                glEnable(GL_COLOR_MATERIAL);
                glColor3f(0.6, 0.6, 0.6);
                glutSolidSphere(0.3, 10, 10);
                
                glPopMatrix();
            }

            
            // neighbor
            if(glut_menu::get_state("Neighbor", 0))
            {
                auto neighbor = m_aniso_kernel.neighbor_search(pos, m_slength*2);
                glDisable(GL_LIGHTING);
                glColor3f(1, 0, 0);
                glBegin(GL_POINTS);
                for (auto p : neighbor)
                {
                    glVertex3dv(m_current_particles[p].pos.get());
                }
                glEnd();
            }
            
            // coord
            if(glut_menu::get_state("PCA coords", 0))
            {
                auto U = m_aniso_kernel.m_U[i];
                auto S = m_aniso_kernel.m_S[i];
                glBegin(GL_LINES);
                for (int i = 0; i < 3; i++)
                {
    //                vec3 cc(U[i][0], U[i][1], U[i][2]);
                    vec3 cc(U[0][i], U[1][i], U[2][i]);
                    cc *= S[i][i]*5000;
                    glVertex3dv(pos.get());
                    glVertex3dv((pos + cc).get());
                }
                glEnd();
            }
//        }
    }
}

void particle_manager::draw_orientation_anisotropic()
{
    
}

void particle_manager::draw_anisotropic_kernel(vec3 domain_size, vec3 c)
{

    
    // Build
    if (pos.size() == 0)
    {
        build_anisotropic_kernel();
        
        int N = 300;
        int gap = N/20;
        vec3 delta = domain_size / N;
        pos.clear(); phi.clear();
        
        for (int i=-gap; i < N+gap; i++)
        {
            //                for (int j = 0; j < N; j++)
            int j = 10;
            {
                for (int k = -gap; k < N+gap; k++)
                    //                    int k = 1;
                {
                    vec3 cur_p(i*delta[0], j*delta[1], k*delta[2]);
                    pos.push_back(cur_p);
                    phi.push_back(m_aniso_kernel.is_inside(cur_p));
                }
            }
        }
    }
    
    
    // Draw
    glDisable(GL_LIGHTING);
    glPointSize(4);
    glBegin(GL_POINTS);
    vec3 RED(1,0,0);
    vec3 BLUE(0,0,1);
    for (int i = 0; i < pos.size(); i++)
    {
        if (phi[i] > 0)
        {
            glColor3f(c[0], c[1], c[2]);
            glVertex3dv(pos[i].get());
        }

    }
    glEnd();
    
}

void particle_manager::draw(double y_under, double y_limit)
{
    glDisable(GL_LIGHTING);
    glPointSize(2);
    glBegin(GL_POINTS);
    //            int idx = 0;
    //            for (int idx : idx_list)
    for(int idx = 0; idx < m_current_particles.size(); idx++)
    {
        auto &p = m_current_particles[idx];
        if (p.pos[1] < y_limit
            && p.pos[1] > y_under)
        {
            
            glVertex3dv(p.pos.get());
        }
    }
    glEnd();
}

bool particle_manager::get_displacement(vec3 pos, vec3 & dis)
{
    return get_displacement_sph_kernel(pos, dis);
}

void particle_manager::load_first_time(int idx){
    load(idx, m_current_particles);
    load(idx+1, m_next_particles);
}

void particle_manager::load_time_step(int idx)
{
    m_current_particles = m_next_particles;
    load(idx+1, m_next_particles);
}
void particle_manager::interpolate(int sub_idx, int sub_count)
{
//    // Load sub step by linear interpolation
//    m_sub_step_particles = m_current_particles;
//
//    for (int i = 0; i<m_current_particles.size(); i++)
//    {
//        auto pre_pos = m_current_particles[i].pos;
//        auto next_pos = m_next_particles[i].pos;
//
//        m_sub_step_particles[i].pos = pre_pos + (next_pos - pre_pos)*(sub_idx/(double)sub_count);
//    }
//
//
//    if(sub_idx == 0)//velocity is unchanged hence only need one computation
//    {
//        double max_vel = 0;
//        m_sub_step_vel.resize(m_current_particles.size());
//        for (int i = 0; i<m_current_particles.size(); i++)
//        {
//            auto pre_pos = m_current_particles[i].pos;
//            auto next_pos = m_next_particles[i].pos;
//
//            m_sub_step_vel[i] = (next_pos - pre_pos)/(double)sub_count;
//
//            max_vel = max(max_vel, m_sub_step_vel[i].length());
//        }
////        cout << "Max sub vel: " << max_vel <<endl;
//    }
//
//
//
//    build_kd_tree(); // may not need to be build every time step
//    rebuild_density(); // Because the density output is different, and I dont know why.
//    build_anisotropic_kernel();
}

double particle_manager::get_max_displacement()
{
    double max_displacement = -INFINITY;
    double avg_displace = 0;
    for (int i = 0; i<m_current_particles.size(); i++)
    {
        auto pre_pos = m_current_particles[i].pos;
        auto next_pos = m_next_particles[i].pos;
        
        max_displacement = std::max(max_displacement, (next_pos - pre_pos).length());
        avg_displace += (next_pos - pre_pos).length();
    }
    avg_displace /= m_current_particles.size();
    
    return max_displacement;
    return avg_displace;
}

void particle_manager::load_time_step()
{

}

void particle_manager::init_first(int idx)
{
    load(idx, m_cache_particles);
    load(idx+1, m_cache_particles_next);
    m_cache_idx = idx;
    
    m_next_particles = m_cache_particles;
    m_current_particles = m_cache_particles;
}

// Only load the particle from file
// No anisotropic, kdtree or rebuild density
void particle_manager::load_next(int idx, double t)
{
    assert(t >= 0 && t < 1);
    
    m_current_particles = m_next_particles;
    load(idx, m_next_particles);
    
//    // Interpolate
//    if(idx == m_cache_idx+1)
//    {
//        m_cache_idx = idx;
//        m_cache_particles = m_cache_particles_next;
//        load(m_cache_idx, m_cache_particles_next);
//    }
//
//    double max_vel = 0;
//    for (int i = 0; i < m_next_particles.size(); i++)
//    {
//        m_next_particles[i].pos = m_cache_particles[i].pos * (1-t) + m_cache_particles_next[i].pos * t;
//        m_current_particles[i].vel = m_next_particles[i].pos - m_current_particles[i].pos;
//
//        max_vel = max(max_vel, m_current_particles[i].vel.length());
//    }
    
//    cout << "Max vel: " << max_vel << endl;
}

void particle_manager::add_ghost_particles(vec3 domain_size, double dis)
{
    build_kd_tree(); // may not need to be build every time step
//    rebuild_density(); // Because the density output is different, and I dont know why.
    
//    // Add ghost
//
//    for(int direct = 0; direct < 3; direct++)
//    {
//        auto current_size = m_current_particles.size();
//        for (int i = 0; i < current_size; i++)
//        {
//            auto & p = m_current_particles[i];
//            if (p.pos[direct] < dis)
//            {
//                particle new_p = p;
//                new_p.pos[direct] = -new_p.pos[direct];
//                new_p.vel[direct] = -new_p.vel[direct];
//                m_current_particles.push_back(new_p);
//            }else if(p.pos[direct] > domain_size[direct]-dis)
//            {
//                particle new_p = p;
//                new_p.pos[direct] = 2*domain_size[direct] - new_p.pos[direct];
//                new_p.vel[direct] = -new_p.vel[direct];
//                m_current_particles.push_back(new_p);
//            }
//        }
//    }
//
//    build_kd_tree(); // may not need to be build every time step
}


bool particle_manager::get_projection(vec3 pos, vec3 direction, bool &bInside, double &t)
{
//    static double eps = std::numeric_limits<double>::min();
//    double phi_pre = m_aniso_kernel.get_value(pos);
//
//    if (phi_pre >= 0)
//    {
//
//    }
//    else // No particle nearby.
//    {
//        return false;
//    }
//
//    double ra = m_deltap; // Important parameter
//
//    
//
//    double epsilon = -ra*(phi_pre>eps? -1:1);
//
//    bInside = phi_pre > 0;
//
//    double phi_epsilon = m_aniso_kernel.get_value(pos + direction*epsilon);
//
//    if (phi_pre * phi_epsilon < 0) //the point is somewhere in between
//    {
//        double ep1 = epsilon;
//        double ep2 = 0;
//        double phi1 = phi_epsilon;
//        double phi2 = phi_pre;
//
//        for (int i = 0; i < 4; i++)
//        {
//            double phi_middle = m_aniso_kernel.get_value(pos + direction*(ep1+ep2)/2);
//            if (std::abs(phi_middle) < 0.0001)
//            {
//                break;
//            }
//
//            if (phi_middle*phi1 < 0)
//            {
//                ep2 = (ep1 + ep2)/2;
//            }
//            else
//            {
//                ep1 = (ep1 + ep2)/2;
//                phi1 = phi_middle;
//            }
//        }
//
//        t = (ep1 + ep2)/2;
//
//        return true;
//    }
    
    //    t = epsilon;
    // fail to project
    return false;
}

void particle_manager::build_anisotropic_kernel()
{
    m_aniso_kernel.m_particles = m_current_particles;
    m_aniso_kernel.m_h = m_slength;
//    m_aniso_kernel.m_ra = m_deltap;
    
    m_aniso_kernel.build();
}
std::shared_ptr<Geometry::KDTree<vec3, int>> particle_manager::build_kd_tree(std::vector<particle> *)
{
    std::shared_ptr<Geometry::KDTree<vec3, int>> vtree = std::shared_ptr<Geometry::KDTree<vec3, int>>(new Geometry::KDTree<vec3, int>());
    
    int idx = 0;
    for (auto &p : m_current_particles)
    {
        vtree->insert(p.pos, idx);
        idx++;
    }
    
    vtree->build();
    
    return vtree;
}
void particle_manager::rebuild_density()
{
    for (auto & part : m_current_particles)
    {
        vector<vec3> neighbor_pos;
        vector<int> neighbor_idx;
        
        double h = m_slength;
        double r = 2*h;
        
        m_vtree.in_sphere(part.pos, r, neighbor_pos, neighbor_idx);

        
        assert(neighbor_idx.size() > 0);
        
        double rho = 0;
        for (auto & p : neighbor_idx)
        {
            double q = (m_current_particles[p].pos - part.pos).length() / h;
            rho += m_current_particles[p].mass * wendland(q, h);
        }
        part.density = rho;
    }
}

void particle_manager::build_kd_tree()
{
    m_vtree = Geometry::KDTree<vec3, int>();
    
    int idx = 0;
    for (auto &p : m_current_particles)
    {
       m_vtree.insert(p.pos, idx);
        idx++;
    }
    
    m_vtree.build();
}

void particle_manager::load(int idx, std::vector<particle> & par)
{
    par.clear();
    
    stringstream s;
    s << m_data_path << "/iter_" << setfill('0') << setw(5) << idx << ".particle";
    
    
    std::ifstream f(s.str());
    if(f.is_open())
    {
        int num_points;
        f >> num_points;
        
        for (int i = 0; i < num_points; i++)
        {
            particle p;
            f >> p;

            par.push_back(p);
        }
    }
    else{
        cout << "No more file to load. Last file: " <<  s.str() << endl;
        exit(EXIT_SUCCESS);
    }
}
