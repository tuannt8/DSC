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

void anisotropic_kernel::build(){
    // We build connected coponent later

    m_G.resize(m_shared_particles.size());
    m_det_G.resize(m_shared_particles.size());
    //        // dam-break hard code test
    //        double ra = 0.02;
    //        double h = 2*ra;
    //        double r = 2*h;
    //        m_h = h;

    // Bubble hardcode test
    double ra = 0.002;
    double h = 2*ra;
    double r = 2*h;
    m_h = h;

    for (int i = 0; i < m_shared_particles.size(); i++)
    {
        if(i==98)
        {

        }

        auto & pi = m_shared_particles.at(i);

        if(pi.type != 0)
            continue;

        std::vector<int> close_particles;
        std::vector<vec3> close_particles_pos;
        m_shared_tree->in_sphere(pi.pos, h, close_particles_pos, close_particles);

        if(close_particles.size() == 0)
        {
            assert(0);
        }

        // Mean position
        vec3 x_w(0.0);
        double sum_omega = 0.0;
        for (auto idx : close_particles)
        {
            auto neighbor_p = m_shared_particles.at(idx);
            vec3 r_v = neighbor_p.pos - pi.pos;
            double omega = 1 - std::pow(r_v.length()/r, 3);
            x_w += neighbor_p.pos*omega;
            sum_omega += omega;
        }
        x_w /= sum_omega;

        // Orientation matrix
        mat3x3d C(0.0);
        for (auto idx : close_particles)
        {
            auto neighbor_p = m_shared_particles.at(idx);
            vec3 r_v = neighbor_p.pos - pi.pos;
            double omega = 1 - std::pow(r_v.length()/r, 3);

            mat3x3d oo;
            vec3 r_w = neighbor_p.pos - x_w;
            for (int ri = 0; ri < 3; ri++)
                for (int ci = 0; ci < 3; ci++)
                    oo[ri][ci] = r_w[ri]*r_w[ci];
            C = C + oo*omega;
        }

        C = C / sum_omega;

        // Solve svd using eigen
        mat3x3d Q(0.0);
        mat3x3d L(0.0);
        svd_solve(C.get(), Q.get(), L.get());
        
        
        
        // Singular value decomposition
        mat3x3d Q1;
        mat3x3d L1;
        int nb_singular = CGLA::power_eigensolution<mat3x3d>(C, Q1, L1);

        mat3x3d G;

        // Modify the strech matrix
        mat3x3d Sigma(0.0);
        double kr = 1.0/4.0;
        double ks = 1400;
        double kn = 0.5;
        for (int d = 0; d < 3; d++)
        {
            if(close_particles.size() > 25)
            {
                Sigma[d][d] = (std::max(L[d][d], L[0][0] * kr)*ks);
            }
            else
            {
                Sigma[d][d] = 1.0 * kn;
            }
        }

        auto CCC = CGLA::transpose(Q)*L*Q;

        // Finally
        //            Sigma[0][0] = 1; Sigma[1][1] = 2; Sigma[2][2] = 1;
        Sigma = L;
        double max = std::max(Sigma[0][0], std::max(Sigma[1][1], Sigma[2][2]));
        Sigma[0][0] /= max; Sigma[1][1] /= max; Sigma[2][2] /= max;

        G = CGLA::transpose(Q)*Sigma*Q;// Q*Sigma*CGLA::transpose(Q) *(1.0/h);


        m_G[i] = G;
        m_det_G[i] = CGLA::determinant(G);
    }
};

