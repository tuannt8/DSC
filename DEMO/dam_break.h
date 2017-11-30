//
//  dam_break.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 30/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef dam_break_h
#define dam_break_h

#include "file_load.hpp"

class dam_break:public file_load
{
public:
    dam_break():file_load()
    {
        
    }
    ~dam_break()
    {
        
    }
    
public:
    virtual void  init_dsc(DSC::DeformableSimplicialComplex<> * dsc)
    {
        vec3 dam_bound(0.4, 0.67, 0.4);
        is_mesh::Cube bound_cube(dam_bound/2.0, dam_bound);
        
        dsc->set_labels(bound_cube, 1);
        
        // Project interface to the Cube
        double thres = dsc->get_avg_edge_length()/1.3;
        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                auto p = nit->get_pos();
                for (int i = 0; i < 3; i++)
                {
                    if (std::abs(p[i]) < thres)
                    {
                        p[i] = 0.0;
                    }
                    
                    if (std::abs(dam_bound[i] - p[i]) < thres)
                    {
                        p[i] = dam_bound[i];
                    }
                }
                
                dsc->set_destination(nit.key(), p);
            }
        }
        
        dsc->deform();
    }
    
    virtual vec3 get_domain_dimension()
    {
        return vec3(1.6, 0.67, 0.6);
    }
    
    virtual double get_influence_radius()
    {
        return 0.052;
    }
};

#endif /* dam_break_h */
