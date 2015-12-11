//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"

void segment_function::init()
{
    _img.load("data/sphere_drill");
}

void segment_function::initialze_segmentation()
{
    // initialize by thresholding
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        
    }
}