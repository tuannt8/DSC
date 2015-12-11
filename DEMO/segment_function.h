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

class segment_function
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
    
public:
    segment_function(){};
    ~segment_function(){};
    
    void init(); //
    void initialze_segmentation();
    
public:
    image3d _img; // Store crossection -> voxel
    dsc_class *_dsc; // Shared dsc
};

#endif /* segment_function_hpp */
