//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"
#include "tet_dis_coord.hpp"

#include "profile.h"

#include <GLUT/GLUT.h>

using namespace std;

void segment_function::init()
{
    //_img.load("data/sphere_drill");
//    _img.load("../Large_data/hamster");
    _img.load("../Large_data/fuel_cells_small");

}

void segment_function::initialze_segmentation()
{
//    /**
//     Hamster sample
//     */
//    // Initilization by thresholding
//    double thres = 0.6;
//    // initialize by thresholding
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        // test
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//        
//        assert(avgI < 1.01);
//        //assert(total < 1000);
//        
//        if (avgI > thres)
//        {
//            _dsc->set_label(tit.key(), 1);
//        }
//    }
    
    /**
     Fuel cells
     */
    // Initialization by thresholding
    double thres[] = {0.31, 0.57, 0.7};
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        // test
        double total, volume;
        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
        
        double var = _img.get_variation(pts, avgI);
        if (var > 0.1)
        {
            continue;
        }

        assert(avgI < 1.01);
        
        // find closest
        int idx = 0;
        if (std::abs(avgI - thres[1]) < 0.1 )
        {
            idx = 1;
        }else if(std::abs(avgI - thres[2]) < 0.1 )
        {
            idx = 2;
        }

        if (idx != 0)
        {
            _dsc->set_label(tit.key(), idx);
        }
    }
}

void segment_function::update_vertex_stability()
{
    _vertex_stability_map = std::vector<int>(30000, 0); // suppose we have less than 10000 vertices
    
    for (auto vid = _dsc->nodes_begin(); vid != _dsc->nodes_end(); vid++)
    {
        if (_forces[(long)vid.key()].length() < 0.1) // stable
        {
            _vertex_stability_map[vid.key()] = 1;
        }
    }
}

void segment_function::compute_external_force()
{
    auto c = _mean_intensities;
    // Compute forces
    std::vector<vec3> forces = std::vector<vec3>(30000, vec3(0.0)); // suppose we have less than 10000 vertices
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            // Normal
            auto tets = _dsc->get_tets(fid.key());
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            double c0 = c[_dsc->get_label(tets[0])];
            double c1 = c[_dsc->get_label(tets[1])];
            
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            size_t n = std::ceil( sqrt(area) );
            if (n >= tri_coord_size.size())
            {
                n = tri_coord_size.size() - 1;
            }
            
            auto a = tri_dis_coord[n - 1];
            for (auto coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);
                
                //                auto f = - Norm* ((c1 - c0)*(2*g - c0 - c1) / area);
                auto f = - Norm* ((2*g - c0 - c1) / (c1-c0) / area);
                
                // distribute
                forces[verts[0]] += f*coord[0];
                forces[verts[1]] += f*coord[1];
                forces[verts[2]] += f*coord[2];
            }
        }
    }

    _forces = forces;
}

void segment_function::face_split()
{
    auto c = _mean_intensities;
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            // check if it is stable first
            auto nodes = _dsc->get_nodes(fid.key());
            if (!(_vertex_stability_map[nodes[0]] == 1
                and _vertex_stability_map[nodes[1]] == 1
              and _vertex_stability_map[nodes[2]] == 1))
            {
                continue;
            }
            
            double var = 0.0;
            // Normal
            auto tets = _dsc->get_tets(fid.key());
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            double c0 = c[_dsc->get_label(tets[0])];
            double c1 = c[_dsc->get_label(tets[1])];
            
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            size_t n = std::ceil( sqrt(area) );
            if (n >= tri_coord_size.size())
            {
                n = tri_coord_size.size() - 1;
            }
            
            auto a = tri_dis_coord[n - 1];
            for (auto coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);
                
                //                auto f = - Norm* ((c1 - c0)*(2*g - c0 - c1) / area);
                auto f = - Norm* ((2*g - c0 - c1) / (c1-c0) / area);
                
                // add
                var += f.length();
            }
            
            var = var / a.size();
            
            if (var > 0.1) // The face want to move
            {
                cout << "Split face " << fid.key() << ", var = " << var << endl;
                _dsc->split(fid.key());
            }
        }
    }

}


/**
 * Bounding box of point list
 */
void bounding_box(const std::vector<vec3> & pts, vec3 & ld, vec3 & ru)
{
    ld = vec3(INFINITY);
    ru = vec3(-INFINITY);
    for(auto & v : pts)
    {
        ld[0] = std::min(v[0], ld[0]);
        ld[1] = std::min(v[1], ld[1]);
        ld[2] = std::min(v[2], ld[2]);
        
        ru[0] = std::max(v[0], ru[0]);
        ru[1] = std::max(v[1], ru[1]);
        ru[2] = std::max(v[2], ru[2]);
    }
}

bool sort_intersect(intersect_pt p1, intersect_pt p2)
{
    {
        return p1.z < p2.z;
    }
}

void segment_function::update_average_intensity()
{
    cout << "Computing average intensity" << endl;
    
    int nb_phase = 3;
    
    // 1. Init the buffer for intersection
    auto dim = _img.dimension();

    std::vector<ray_z> init_rayz(dim[0] * dim[1]);
    for (int y = 0; y < dim[1]; y++)
    {
        for (int x = 0; x < dim[0]; x ++)
        {
            int idx = y*dim[0] + x;
            init_rayz[idx].x = x;
            init_rayz[idx].y = y;
        }
    }
    
    vector<std::vector<ray_z>> ray_intersect(nb_phase, init_rayz);

    cout << "Find intersection" << endl;
    // 2. Find intersection with interface
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            auto tet = _dsc->get_tets(fid.key());
            auto phase0 = _dsc->get_label(tet[0]);
            auto phase1 = _dsc->get_label(tet[1]);
            
            // check all z-ray that intersect this triangle
            auto pts3 = _dsc->get_pos(_dsc->get_nodes(fid.key()));
            auto pts = pts3;
            pts[0][2] = 0; pts[1][2] = 0; pts[2][2] = 0;
            
            auto n = _dsc->get_normal(fid.key(), tet[0]);
            bool in_1 = Util::dot(n, vec3(0,0,1)) > 0;
            
            vec3 ld, ru;
            bounding_box(pts, ld, ru);
            for (int x = std::floor(ld[0]); x < std::round(ru[0]); x++)
            {
                for (int y = std::floor(ld[1]); y < std::round(ru[1]); y++)
                {
                    try
                    {
                        auto bc = Util::barycentric_coords<double>(vec3(x+0.5, y+0.5, 0), pts[0], pts[1], pts[2]);
                        
                        if (bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                        { // inside
                            auto p = pts3[0]*bc[0] + pts3[1]*bc[1] + pts3[2]*bc[2];
                            
                            assert(p[2] >= 0 and p[2] < 9999);
                            
                            ray_intersect[phase0][y*dim[0] + x].intersects.push_back(intersect_pt(std::floor(p[2]), !in_1));
                            ray_intersect[phase1][y*dim[0] + x].intersects.push_back(intersect_pt(std::floor(p[2]), in_1));
                        }
                    }
                    catch (std::exception e)
                    {
                    
                    }
                }
            }
        }
    }
    
    cout << "Count intersection" << endl;
    // 3. Compute integral
    _d_rayz.clear();
    
    std::fill(_mean_intensities.begin(), _mean_intensities.end(), 0.0);
    vector<double> area(_mean_intensities.size(), 0.0);
    for (int i = 1; i < nb_phase; i++) // we dont compute the background
    {
        int count = 0;
        for (auto r : ray_intersect[i])
        {
            std::vector<intersect_pt> intersect_ps = r.intersects;
            if (intersect_ps.size() > 1)
            {
                std::sort(intersect_ps.begin(), intersect_ps.end(), sort_intersect);
                
                // remove identical intersections
                for (auto p = intersect_ps.begin()+1; p != intersect_ps.end(); p++)
                {
                    auto pre = p -1;
                    if (p->z == pre->z and
                        p->b_in == pre->b_in)
                    {
                        p = intersect_ps.erase(p);
                        
                        if(p == intersect_ps.end())
                        {
                            break;
                        }
                    }
                }
                // Now count
                if (intersect_ps.size() % 2 == 0) // Should be even
                {
                    count++;
                    auto newR = r;
                    newR.intersects = intersect_ps;
                    _d_rayz.push_back(newR);
                    
                    for (int j = 0; j < intersect_ps.size()/2; j++)
                    {
                        int z1 = intersect_ps[2*j].z;
                        int z2 = intersect_ps[2*j + 1].z;
                        
                        area[i] += z2 - z1;
                        _mean_intensities[i] += _img.sum_line_z(r.x, r.y, z1, z2);
                    }
                }
            }
        }
        cout << count << "Intersected rays" << endl;
    }
    
    cout << "Update phase 0" << endl;
    double s0 = _img.sum_area(dim[0]-1, dim[1]-1, dim[2]-1);
    double v0 = dim[0]*dim[1]*dim[2];
    for (int i = 1; i < nb_phase; i++) // we dont compute the background
    {
        s0 -= _mean_intensities[i];
        v0 -= area[i];
    }
    _mean_intensities[0] = s0;
    area[0] = v0;
    
    for (int i  = 0; i < nb_phase; i++)
    {
        _mean_intensities[i] /= area[i];
    }
    
    
    cout << "Done mean intensity" << endl;
    for (int i = 0; i < nb_phase; i++)
    {
        cout << _mean_intensities[i] << " -- ";
    }
    cout << endl;

}

void segment_function::segment()
{
    profile t("Average intensity");
    
    // Compute average intensity
    int num_phases = 3;
//    std::vector<double> c = {0.0 , 0.0};
//    std::vector<double> vols = {0, 0};
//    
//    for (auto tet = _dsc->tetrahedra_begin(); tet != _dsc->tetrahedra_end(); tet++)
//    {
//        double total = 0, volume = 0;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tet.key()));
//        double avg = _img.get_tetra_intensity(pts, &total, &volume);
//        assert( avg < 1.01);
//        assert(total != NAN);
//        int idx = _dsc->get_label(tet.key());
//        c[idx] += total;
//        vols[idx] += volume;
//    }
//    
//    for (int i = 0; i < num_phases; i++)
//    {
//        c[i] = c[i] / vols[i];
//    }
//    
//    t.done();
    
    _mean_intensities.resize(num_phases);
    // This function does not return correct result (in compare to above algorithm).
    // Need further investigation
    update_average_intensity();
    
    /**
     RELABEL TETRAHEDRAL
     */
    static int mesh_opt_counter = 0;
    if (mesh_opt_counter > 20)
    {
        // Perform relabeling
        for (auto tid = _dsc->tetrahedra_begin(); tid != _dsc->tetrahedra_end(); tid++)
        {
            // Evaluate the tetrahedral
            auto pts = _dsc->get_pos(_dsc->get_nodes(tid.key()));
            double meanc = _img.get_tetra_intensity(pts);
            
            // Find new phase
            int phase = -1;
            double min_dis = INFINITY;
            for (int i = 0; i < _mean_intensities.size(); i++)
            {
                if (min_dis > std::abs(_mean_intensities[i] - meanc))
                {
                    min_dis = std::abs(_mean_intensities[i] - meanc);
                    phase = i;
                }
            }
            
            if (phase != _dsc->get_label(tid.key()))
            {
                double var = _img.get_variation(pts, meanc);
                double thres = 0.2*std::abs(_mean_intensities[_dsc->get_label(tid.key())]
                                            - _mean_intensities[phase]);
                
                cout << "thres hold: " << thres << endl;
                if (var > thres) // Large variation, split
                {
                    cout << "Relabel " << tid.key() << " from " << _dsc->get_label(tid.key())
                    << "to" << phase << endl;
                    _dsc->set_label(tid.key(), phase);
                }
                else // If they are stable, split the tetrahedral
                {
                    
                }
            }
        }
        
        // Face spliting
        compute_external_force();
        update_vertex_stability();
        face_split();
    }
    
    compute_external_force();
    
    
    t.done();
    t = profile("Apply force");
    
    double largest = 0;
    for(auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {

        if ( (nid->is_interface() or nid->is_crossing())
            and _dsc->is_movable(nid.key()))
        {
            auto dis = _forces[nid.key()]*_dt;
            //cout << "Node " << nid.key() << " : " << dis << endl;
            _dsc->set_destination(nid.key(), nid->get_pos() + dis);
            if (largest < dis.length())
            {
                largest = dis.length();
            }
        }
    }
    
    t.done();
    t = profile("Deform");

    cout << "--------------------------------Max displacement: " << largest << endl;
    _dsc->deform();
}
