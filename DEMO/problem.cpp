//
//  problem.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 01/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#include <stdio.h>
#include "problem.h"
#include "config_file.h"

using namespace std;

double ks; // Param for isotropic

void problem::init(std::string sumary_file_path)
{
    cout << sumary_file_path << endl;
    
    config_file c_f(sumary_file_path);
    
    m_deltap = c_f.get_double("deltap");
    m_influenceRadius = c_f.get_double("influenceRadius");
    m_slength = c_f.get_double("slength");
    m_nb_phases = c_f.get_int("numFluids");
    
    ks = get_ks();
}

DSC::DeformableSimplicialComplex<> * problem::init_dsc_domain(double scale)
{
    double edge_length = m_deltap*scale;
    
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;
    
    double delta = edge_length;
    
    vec3 _dsc_dim = domain_size() + vec3(delta)*2;
    
    cout << "delta " << delta << endl;
    cout << _dsc_dim[0] << " " << _dsc_dim[1] << " " <<  _dsc_dim[2] << " " ;
    
    int NX = round(_dsc_dim[0] / delta) + 1; // number of vertices
    int NY = round(_dsc_dim[1] / delta) + 1;
    int NZ = round(_dsc_dim[2] / delta) + 1;
    
    cout << "Compute point" << NX << " " << NY << " " << NZ << "\n";
    
    double deltax = _dsc_dim[0]/(double)(NX-1);
    double deltay = _dsc_dim[1]/(double)(NY - 1);
    double deltaz = _dsc_dim[2]/(double)(NZ - 1);
    
    
    // points. Push it back
    for (int iz = 0; iz < NZ; iz++)
    {
        for (int iy = 0; iy < NY; iy++)
        {
            for (int ix = 0; ix < NX; ix++)
            {
                points.push_back(vec3(ix*deltax, iy*deltay, iz*deltaz) - vec3(delta));
            }
        }
    }
    
    cout << "Compute tets\n";
    
    // tets
    for (int iz = 0; iz < NZ - 1; iz++)
    {
        for (int iy = 0; iy < NY - 1; iy++)
        {
            for (int ix = 0; ix < NX - 1; ix++)
            {
                // 8 vertices
                int vertices[] = {
                    index_cube(ix, iy, iz),
                    index_cube(ix+1, iy, iz),
                    index_cube(ix+1, iy+1, iz),
                    index_cube(ix, iy+1, iz),
                    index_cube(ix, iy, iz + 1),
                    index_cube(ix+1, iy, iz + 1),
                    index_cube(ix+1, iy+1, iz + 1),
                    index_cube(ix, iy+1, iz + 1)
                };
                
                int tetras[] = {
                    0, 4, 5, 7,
                    0, 7, 5, 1,
                    0, 1, 3, 7,
                    1, 5, 6, 7,
                    1, 6, 7, 3,
                    1, 2, 6, 3
                };
                
                for(int i = 0; i < 6*4; i++)
                {
                    tets.push_back(vertices[tetras[i]]);
                }
            }
        }
    }
    
    long nbTet = tets.size()/4;
    tet_labels = std::vector<int>(nbTet, 0);
    
    cout << "Init DSC from point\n";
    
    auto dsc = new DSC::DeformableSimplicialComplex<>(points, tets, tet_labels);
    dsc->set_avg_edge_length(edge_length);
    
    return dsc;
}
bool two_phase_fluid::is_one_bellow_z(DSC::DeformableSimplicialComplex<> * dsc, is_mesh::TetrahedronKey tkey, double z)
{
    auto node_poses = dsc->get_pos(dsc->get_nodes(tkey));
    
    for (auto n : node_poses)
    {
        if (n[2] <z)
        {
            return true;
        }
    }
    return false;
}
DSC::DeformableSimplicialComplex<> * two_phase_fluid::init_dsc(double scale)
{
    auto dsc = init_dsc_domain(scale);
    
    double z1 = 0.054, z2 = 0.125;
    
    // Init fluid
    auto fluid_size = domain_size();
    for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
    {
        auto center = dsc->barycenter(tit.key());
        
        if (center[0] < 0 || center[0] > fluid_size[0]
            || center[1] < 0 || center[1] > fluid_size[1]
            || center[2] < 0 || center[2] > fluid_size[2])
        {
            continue;
        }
        
        if (is_one_bellow_z(dsc, tit.key(), z1))
        {
            dsc->set_label(tit.key(), 1);
        }else if(is_one_bellow_z(dsc, tit.key(), z2))
        {
            dsc->set_label(tit.key(), 2);
        }
    }
    
    auto avg_edge = dsc->get_avg_edge_length()*1.3;
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto pos = nit->get_pos();
            if (std::abs(pos[2] - z1) < avg_edge)
            {
                pos[2] = z1;
            }
            if (std::abs(pos[2] - z2) < avg_edge)
            {
                pos[2] = z2;
            }

            dsc->set_destination(nit.key(), pos);
        }
    }

    dsc->deform();
    
    return dsc;
    
    
}

DSC::DeformableSimplicialComplex<> * bubble_fluid ::init_dsc(double scale)
{
    auto dsc = init_dsc_domain(scale);
    
    // Label the bubble
    vec3 center(0.0848462, 0.0848462, 0.05);
    double R = 0.025;
    for (auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
    {
        bool bInside = true;
        for (auto p : dsc->get_pos(dsc->get_nodes(tit.key())))
        {
            if ((p - center).length() > R)
            {
                bInside = false;
                break;
            }
        }
        if (bInside)
        {
            dsc->set_label(tit.key(), 1);
        }
    }
    // Project the surface
    for(auto nit = dsc->nodes_begin(); nit!= dsc->nodes_end(); nit++)
    {
        if(nit->is_interface())
        {
            auto pos = nit->get_pos();
            auto r_pos = pos - center;

            r_pos.normalize();
            vec3 new_pos = center + r_pos*R;

            dsc->set_destination(nit.key(), new_pos);
        }
    }

    dsc->deform();

    return dsc;
}

DSC::DeformableSimplicialComplex<> * dam_break_fluid::init_dsc(double scale)
{
    auto dsc = init_dsc_domain(scale);
    
    vec3 fluid_ld(0.0);
    vec3 fluid_ur(0.4, 0.67, 0.4);
    is_mesh::Cube c((fluid_ld+fluid_ur)/2.0, fluid_ur - fluid_ld);
    
    dsc->set_labels(c, 1);
    
//    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
//    {
//        auto pos = nit->get_pos();
//        bool is_inside = false;
//
//        if (c.is_inside(nit->get_pos()))
//        {
//            for(auto t : dsc->get_tets(nit.key()))
//            {
//                if (dsc->get_label(t) == 0)
//                {
//                    dsc->set_label(t, 1);
//                }
//            }
//        }
//    }
    // Make wu
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if (nit->is_boundary())
        {
            for(auto t : dsc->get_tets(nit.key()))
                dsc->set_label(t, 0);
        }
    }
    
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto pos = nit->get_pos();
            for (int i = 0; i<3; i++)
            {
                pos[i] = max(pos[i], fluid_ld[i]);
                pos[i] = min(pos[i], fluid_ur[i]);
            }
        }
    }

    dsc->deform(20);
    
    return dsc;
}
