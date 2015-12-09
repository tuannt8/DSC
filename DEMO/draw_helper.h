//
//  draw_helper.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef draw_helper_hpp
#define draw_helper_hpp

#include "image3d.h"
#include <stdio.h>

namespace draw_helper
{
    /*
     GL function
     */
    
    /*
     Draw the 3d image
     */
    void draw_image_slice(image3d & im);
}

#endif /* draw_helper_hpp */
