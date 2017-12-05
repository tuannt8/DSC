//
//  anisotrpic_kernel.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 04/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef anisotrpic_kernel_h
#define anisotrpic_kernel_h

#include "KDTree.h"
#include "define.h"
#include "particle.h"


class anisotropic_kernel
{
public:
    Geometry::KDTree<vec3, int> *m_shared_tree;
    double m_h;//influence radius of the particles
    std::vector<mat3x3d> m_G;
    std::vector<double> m_det_G;
    std::vector<particle> m_shared_particles;
    
public:
    anisotropic_kernel(){};
    ~anisotropic_kernel(){};
    
    void build();
    
    double get_value(vec3 pos){
        // dam-break hard code test
        
        double h = 0.02;
        double r = h;
        
        static double kernel_sigma = 315.0 / (64 * 3.14159 * std::pow(m_h, 6));
        
        std::vector<int> close_particles;
        std::vector<vec3> close_particles_pos;
        m_shared_tree->in_sphere(pos, r, close_particles_pos, close_particles);
        
        if(close_particles.size() == 0)
        {
            return 0;
        }
        double phi = 0.0;
        for (auto n_p : close_particles)
        {
            auto part = m_shared_particles.at(n_p);
            auto & G = m_G[n_p];
            vec3 ro = pos-part.pos;
            vec3 ra = G*(pos-part.pos);
            // Use WENLAND kernel
            double contribute = kernel_sigma/std::pow(m_h,3)*std::pow(h*h - Util::dot(ra, ra), 3)*m_det_G[n_p];
            
            phi += part.mass/part.density * contribute;
        }
        
        return phi;
    };
};

#endif /* anisotrpic_kernel_h */
