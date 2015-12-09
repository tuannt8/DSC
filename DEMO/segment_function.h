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

class segment_function
{
public:
    segment_function(){};
    ~segment_function(){};
    
    void init();
    image3d & get_image(){return _img;};
private:
    image3d _img;
};

#endif /* segment_function_hpp */
