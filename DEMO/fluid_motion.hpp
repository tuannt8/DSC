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
#include "DSC.h"

class fluid_motion
{
public:
    fluid_motion();
    ~fluid_motion(){};
    
    void draw();

public:
    DSC::DeformableSimplicialComplex<>* s_dsc;
    vtkWrapper m_vtkWrapper;
    
    void deform();
    
    void log_dsc_surface(int idx);
};


#endif /* fluid_motion_hpp */
