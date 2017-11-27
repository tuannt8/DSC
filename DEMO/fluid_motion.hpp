//
//  fluid_motion.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef fluid_motion_hpp
#define fluid_motion_hpp

#include <stdio.h>
#include "vtkWrapper.hpp"

class fluid_motion
{
public:
    fluid_motion();
    ~fluid_motion(){};
    
    void draw();

public:
    vtkWrapper m_vtkWrapper;
};


#endif /* fluid_motion_hpp */
