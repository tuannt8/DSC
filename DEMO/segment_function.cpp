//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"

using namespace std;

void segment_function::init()
{
    _img.load("data/sphere_drill");
}

void segment_function::initialze_segmentation()
{
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
        assert(total < 1000);
        
        if (avgI > thres)
        {
            _dsc->set_label(tit.key(), 1);
        }
     //   cout << "Inten: " << avgI << endl;
    }
}

void segment_function::segment()
{
//    // Compute average intensity
//    int num_phases = 2;
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
    // Compute force
    for (auto vit = _dsc->nodes_begin(); vit != _dsc->nodes_end(); vit++)
    {
        
//        if (vit->is_interface()
//            and !vit->is_crossing())
//        {
//            vec3 f = _dsc->get_normal(vit.key());
//            auto pos = vit->get_pos();
//            double g = _img.get_value_f(pos);
//            
//            f = - f * ((c[0] - c[1])*(2*g - c[0] - c[1])) * 0.0;
//            
//            std::cout << vit.key() <<  " displace " << f.length() << std::endl;
//            
//            if (_dsc->is_movable(vit.key()))
//            {
//                _dsc->set_destination(vit.key(), vit->get_pos() + f);
//            }
//        }
        
//        if (vit->is_interface())
        {
            _dsc->set_destination(vit.key(), vit->get_pos());
        }
    }
    _dsc->deform();
}
