//
//  fluid_function.h
//  DSC
//
//  Created by Tuan Nguyen Trung on 12/22/14.
//  Copyright (c) 2014 Asger Nyman Christiansen. All rights reserved.
//

#ifndef __DSC__fluid_function__
#define __DSC__fluid_function__

#include "velocity_function.h"
#include "define_fluid.h"

struct dynamic_parameters
{
    float ro;
    float gama;
    float viscosity;
};

class fluid_function: public DSC::VelocityFunc<>{
    
public:
    /**
     Creates a velocity function which moves the interface vertices in the
     normal direction.
     */
    fluid_function(real velocity, real accuracy, int max_time_steps = 200)
    :VelocityFunc<>(velocity, accuracy, max_time_steps)
    {
        
    }
    
    /**
     Returns the name of the velocity function.
     */
    virtual std::string get_name() const
    {
        return std::string("FLUID_MOTION");
    }
    
    /**
     Computes the motion of each interface vertex and stores the new position in
     new_pos in the simplicial complex class.
     */
    virtual void deform(DSC::DeformableSimplicialComplex<>& dsc);
    
private:
    dynamic_parameters param_;
    
    
};


#endif /* defined(__DSC__fluid_function__) */
