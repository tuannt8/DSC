//
//  cache.cpp
//  DSC_parallel
//
//  Created by Tuan Nguyen Trung on 8/1/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#include "cache.hpp"

dsc_cache::dsc_cache()
{
    vertex_flag_clean.reset();
    edge_flag_clean.reset();
    face_flag_clean.reset();
}