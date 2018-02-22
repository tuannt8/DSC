//
//  anisotrpic_kernel.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 05/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include <stdio.h>

#include "anisotrpic_kernel.h"
#include "eigen_wrapper.hpp"

#include "eigensolution.h"
#include <queue>


using namespace std;

vec3 anisotropic_kernel::estimate_norm(vec3 pos)
{
    auto neighbor = neighbor_search(pos, m_r);
    if (neighbor.size() == 0)
    {
        return vec3(0);
    }
    
    vec3 norm(0);
    for(auto n : neighbor){
        norm += (pos - m_particles[n].pos);
    }
    norm.normalize();
    return norm;
}

vec3 anisotropic_kernel::get_displacement_projection(vec3 pos, vec3 norm, double max_displace){
    
    bool last;
    return get_displacement_projection(pos, norm, max_displace, last);
}

vec3 anisotropic_kernel::get_displacement_projection(vec3 pos, vec3 norm, double max_displace, bool & bLast){
    
//    norm = estimate_norm(pos);
    
    // point 1
    vec3 pos1 = pos;
    bool bInside1 = is_inside(pos1);
    
    // point 2
    vec3 pos2 = pos + norm * (max_displace * (bInside1? 1 : -1));
    bool bInside2 = is_inside(pos2);
    
    if (bInside1 == bInside2){
        bLast = false;
        return pos2 - pos1;
    }
    else{
        for (int i = 0; i < 4; i++)
        {
            vec3 mid = (pos1 + pos2)*0.5;
            bool bInside = is_inside(mid);
            
            if (bInside == bInside1)
            {
                pos1 = mid;
            }else
                pos2 = mid;
        }
        
        auto dis = (pos1 + pos2)*0.5 - pos;
        if (dis.length() > 1.001*max_displace)
        {
            cout << pos;
            cout << dis;
            cout << max_displace;
            assert(0);
        }
        
        bLast = true; // Binary search means this is the last projection
        return dis;
    }
}

bool anisotropic_kernel::get_projection(vec3 pos, vec3 direction, bool &bInside, vec3& projected_point)
{
    double max_displace = m_h;
    bInside = is_inside(pos);
    
    vec3 pos1 = pos;
    vec3 pos2 = pos + direction * max_displace;
    bool bInside_2 = is_inside(pos2);
    bool bInside_1 = bInside;
    
    // Two sampling points are both outside or inside
    if (bInside_1 == bInside_2)
    {
        pos2 = pos - direction * max_displace;
        bInside_2 = is_inside(pos2);
    }
    
    if (bInside_1 == bInside_2) // The point is completely inside or outside
    {
        return false;
    }else{
        for (int i = 0; i < 4; i++)
        {
            auto mid = (pos1 + pos2)*0.5;
            bool bInside_mid = is_inside(mid);
            if (bInside_1 == bInside_mid)
            {
                pos1 = mid;
            }else{
                pos2 = mid;
            }
        }
        
        projected_point = (pos1 + pos2)*0.5;
        
        return true;
    }
}

double anisotropic_kernel::get_value(vec3 pos){
    // dam-break hard code test
    static double kernel_sigma = 315.0 / (64 * 3.14159);
    
    std::vector<int> close_particles = neighbor_search(pos, m_r);
    
    if(close_particles.size() == 0)
    {
        return 0; // OUTSIDE
    }
    
    double phi = 0.0;
    for (auto n_p : close_particles)
    {
        auto part = m_particles.at(n_p);
        auto & G = get_transform_mat(n_p);

        vec3 ra_h = G*(pos-part.pos);
        
        
        if (ra_h.length() < 1)
        {
            return 1;
            phi += part.mass/part.density * kernel_sigma *m_det_G[n_p]
                * std::pow(1 - Util::dot(ra_h, ra_h), 3);
        }
    }
    
    return phi;
};

bool anisotropic_kernel::is_inside(vec3 pos)
{
    if (m_b_kernel_computed.size()==0)
    {
        return false;
    }
    
    auto neighbor = neighbor_search(pos, m_r);
    return is_inside(pos, neighbor);
}

bool anisotropic_kernel::is_inside(vec3 pos, std::vector<int> & neighbor)
{
    for (auto & n_p : neighbor)
    {
        auto part = m_particles.at(n_p);
        auto & G = get_transform_mat(n_p);
        
        vec3 ra_h = G*(pos-part.pos);

        if (ra_h.length() < 1)
        {
            return true;
        }
    }
    
    return false;
    
    
}

void anisotropic_kernel::compute_kd_tree()
{
    m_vtree = Geometry::KDTree<vec3, int>();
    

    for (int idx = 0; idx < m_particles.size(); idx++)
    {
        auto & p = m_particles.at(idx);
        m_vtree.insert(p.pos, idx);
    }

    m_vtree.build();
}
std::vector<int> anisotropic_kernel::neighbor_search(vec3 pos, double radius)
{
    std::vector<int> close_particles;
    std::vector<vec3> close_particles_pos;
    m_vtree.in_sphere(pos, radius, close_particles_pos, close_particles);
    
    return close_particles;
    
}
void anisotropic_kernel::build_connected_component()
{
    m_connected_component_label = std::vector<int>(m_particles.size(), -1);
    
    // m_h or 1.1 m_ra
    double connected_radius = 1.1*m_ra;
    
    int cc_count = 0;
    int num_fluid_particles = 0;
    for (int i = 0; i < m_particles.size(); i++)
    {
        num_fluid_particles++;
        
        if (m_connected_component_label[i] == -1)
        {
            m_connected_component_label[i] = cc_count;
            
            std::queue<int> neighbor;
            neighbor.push(i);
            
            while (!neighbor.empty())
            {
                int cur_idx = neighbor.front();
                neighbor.pop();
                
                auto & p_cur = (m_particles)[cur_idx];
                
                // growing the neighbor
                // Search for neighbor particles
                std::vector<int> close_particles = neighbor_search(p_cur.pos, connected_radius);;
                
                for (auto nidx : close_particles)
                {
                    if (m_connected_component_label[nidx] == -1)
                    {
                        neighbor.push(nidx);
                    }
                    m_connected_component_label[nidx] = cc_count;
                }
            }
            
            cc_count++;
        }
    }
    
    m_nb_component = cc_count;
//    cout << cc_count << " connected components of " << num_fluid_particles << " particles" << endl;
}

void draw_connected_component()
{
    
}

void anisotropic_kernel::build(){
    
    m_r = 2*m_h;// Two times the smoothing radius
    
    // 1. Connected component
    compute_kd_tree();
    build_connected_component();
//    // 2. First smooth the particle
//    //  Smooth makes better fit, but in case of two phase, there is gap between phases
//    Taubin_smooth();
//    compute_kd_tree();

    // 3. Build transformation matrix
    m_G.resize(m_particles.size());
    m_det_G.resize(m_particles.size());
    m_b_kernel_computed = vector<bool>(m_particles.size(), false);
    
    m_U.resize(m_particles.size());
    m_S.resize(m_particles.size());
    m_C.resize(m_particles.size());
};

const mat3x3d & anisotropic_kernel::get_transform_mat(int idx)
{
    if (!m_b_kernel_computed[idx])
    {
        m_b_kernel_computed[idx] = true;
        compute_tranformation_mat_for_particle(idx);
    }
    
    return m_G[idx];
}

void anisotropic_kernel::compute_tranformation_mat_for_particle(int i)
{
    auto & pi = m_particles.at(i);
    
    std::vector<int> close_particles = neighbor_search(pi.pos, m_r);
    
    // Weighted mean position
    vec3 x_w(0.0);
    double sum_omega = 0.0;
    for (auto idx : close_particles)
    {
        auto neighbor_p = m_particles.at(idx);
        double omega = weight_func(i, idx, m_r);
        x_w += neighbor_p.pos*omega;
        sum_omega += omega;
    }
    x_w /= sum_omega;
    
    // Orientation matrix
    mat3x3d C(0.0);
    for (auto idx : close_particles)
    {
        auto neighbor_p = m_particles.at(idx);
        double omega = weight_func(i, x_w, idx, m_r);
        
        mat3x3d oo;
        vec3 r_w = neighbor_p.pos - x_w;
        
        // oo = rw*rw'
        CGLA::outer_product(r_w, r_w, oo);
        C = C + oo*omega;
    }
    C /= sum_omega;
    
    // Solve svd using eigen
    mat3x3d Q(0.0);
    mat3x3d L(0.0);
    // C = Q*L*Q'
    svd_solve(C.get(), Q.get(), L.get());
    
    // Modify the strecth matrix
    mat3x3d Sigma(0.0);
    double kr = 4.0;
    double ks = 160000;//cbrt(1 / CGLA::determinant(C));// May optimize latter
//    cout << ks << endl;
    double kn = 0.5;
    for (int d = 0; d < 3; d++)
    {
        if(close_particles.size() > 25)
        {
            Sigma[d][d] = std::max(L[d][d], L[0][0] / kr) * ks;
        }
        else
        {
            Sigma[d][d] = 1.0 * kn;
        }
    }
    
    mat3x3d Sigma_inv(0.0);
    Sigma_inv[0][0] = 1.0 / Sigma[0][0];
    Sigma_inv[1][1] = 1.0 / Sigma[1][1];
    Sigma_inv[2][2] = 1.0 / Sigma[2][2];

    m_G[i] = Q*(Sigma_inv*(CGLA::transpose(Q) *(1/m_h)));
    m_det_G[i] = CGLA::determinant(m_G[i]);
    
    m_U[i] = Q;
    m_S[i] = L;
    m_C[i] = Q*Sigma*CGLA::transpose(Q)*m_h;
    // det(G) is equivilant to 1/h^3
    // Gr is equivilant to r/h
}


double anisotropic_kernel::get_coeff(vec3 pos, int idx)
{
    static double kernel_sigma = 315.0 / (64 * 3.14159);
    
    auto & part = (m_particles)[idx];
    auto & G = get_transform_mat(idx);
    
    vec3 ra_h = G*(pos-part.pos);
    
    if (ra_h.length() < 1)
    {
        return part.mass/part.density * kernel_sigma *m_det_G[idx]
        * std::pow(1 - Util::dot(ra_h, ra_h), 3);
    }
    else
        return 0;
}

void anisotropic_kernel::Taubin_smooth()
{
    double lamda = 0.93;
    std::vector<particle> smoothed_particles = m_particles;
    
    for (int i = 0; i < m_particles.size(); i++)
    {
        auto cur_particle = (m_particles)[i];
        
        // Search for neighbor particles
        std::vector<int> close_particles = neighbor_search(cur_particle.pos, m_r);
        
        if (close_particles.size() > 0)
        {
            vec3 new_pos(0.0);
            double sum_omega = 0;
            for (auto pidx : close_particles)
            {
                auto & pj = (m_particles)[pidx];
                
                double omega = weight_func(i, pidx, m_r);
                sum_omega += omega;
                new_pos += pj.pos*omega;
            }
            new_pos /= sum_omega;
            new_pos = new_pos*lamda + cur_particle.pos*(1-lamda);
            
            smoothed_particles[i].pos = new_pos;
        }
    }
    
    m_particles = smoothed_particles;
}

double anisotropic_kernel::weight_func(int i, vec3 posi, int j, double radius)
{
    if (m_connected_component_label[i] != m_connected_component_label[j])
    {
        return 0;
    }
    
    auto const & pj = m_particles.at(j);
    double r = (posi - pj.pos).length() / radius;
    if (r > 1)
    {
        return 0;
    }
    
    return 1 - r*r*r;
}

double anisotropic_kernel::weight_func(int i, int j, double h)
{
    if (m_connected_component_label[i] != m_connected_component_label[j])
    {
        return 0;
    }
    
    auto const & pi = m_particles.at(i);
    auto const & pj = m_particles.at(j);
    double r = (pi.pos - pj.pos).length() / h;
    if (r > 1)
    {
        return 0;
    }
    
    return 1 - r*r*r;
}
