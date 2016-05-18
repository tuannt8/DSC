//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"
#include "tet_dis_coord.hpp"

#include <GLUT/GLUT.h>

using namespace std;

void segment_function::init()
{
    //_img.load("data/sphere_drill");
    _img.load("../Large_data/hamster");

}

void segment_function::initialze_segmentation()
{
//    printf(" { \n");
//    for (int i = 0; i < 15; i++)
//    {
//        _img.generate_sample_point_tri(i);
//    }
//    printf("};\n");
//    
//    
//    return;
    
    // Initilization by thresholding
    double thres = 0.6;
    // initialize by thresholding
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        // test
        double total, volume;
        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
        
        assert(avgI < 1.01);
        //assert(total < 1000);
        
        if (avgI > thres)
        {
            _dsc->set_label(tit.key(), 1);
        }
    }
}

void segment_function::segment()
{
    // Compute average intensity
    int num_phases = 2;
    std::vector<double> c = {0.0 , 0.0};
    std::vector<double> vols = {0, 0};
    
    for (auto tet = _dsc->tetrahedra_begin(); tet != _dsc->tetrahedra_end(); tet++)
    {
        double total = 0, volume = 0;
        auto pts = _dsc->get_pos(_dsc->get_nodes(tet.key()));
        double avg = _img.get_tetra_intensity(pts, &total, &volume);
        assert( avg < 1.01);
        assert(total != NAN);
        int idx = _dsc->get_label(tet.key());
        c[idx] += total;
        vols[idx] += volume;
    }
    
    for (int i = 0; i < num_phases; i++)
    {
        c[i] = c[i] / vols[i];
    }
    
    std::vector<vec3> forces = std::vector<vec3>(30000, vec3(0.0)); // supose we have less than 10000 vertices
    
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
            
            auto & a = tri_dis_coord[n - 1];
            for (auto & coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);
                
                auto f = - Norm* ((c1 - c0)*(2*g - c0 - c1) / area);
                // distribute
                forces[verts[0]] += f*coord[0];
                forces[verts[1]] += f*coord[1];
                forces[verts[2]] += f*coord[2];
            }
        }
    }
    
    double largest = 0;
    for(auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {

        if ( (nid->is_interface() or nid->is_crossing())
            and _dsc->is_movable(nid.key()))
        {
            auto dis = forces[nid.key()]*2;
            //cout << "Node " << nid.key() << " : " << dis << endl;
            _dsc->set_destination(nid.key(), nid->get_pos() + dis);
            if (largest < dis.length())
            {
                largest = dis.length();
            }
        }
    }

    cout << "--------------------------------Max displacement: " << largest << endl;
    _dsc->deform();
}
