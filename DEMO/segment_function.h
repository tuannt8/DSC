//
//  segment_function.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef segment_function_hpp
#define segment_function_hpp

#include <stdio.h>
#include "image3d.h"
#include "DSC.h"

struct intersect_pt
{
    intersect_pt(int z_, bool b_in_z_){
        z = z_;
        b_in = b_in_z_;
    }
    int z;
    bool b_in;
};

struct ray_z
{
    int x, y;
    std::vector<intersect_pt> intersects;
};

class segment_function
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
    
public:
    segment_function(){};
    ~segment_function(){};
    
    void init(); //
    void initialze_segmentation();
    
    void segment();
    
public:

    image3d _img; // Store crossection -> voxel
    dsc_class *_dsc; // Shared dsc
    
private:
    std::vector<double> _mean_intensities;
    void update_average_intensity();
    
public:
    std::vector<ray_z> _d_rayz;
};

#endif /* segment_function_hpp */
