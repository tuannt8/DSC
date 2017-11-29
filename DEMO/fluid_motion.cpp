//
//  fluid_motion.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "fluid_motion.hpp"

fluid_motion::fluid_motion()
{
    
}

void fluid_motion::draw()
{
    m_vtkWrapper.draw();
}

void fluid_motion::deform()
{
    // 1. Interpolate the displacement
    static int idx = 0;
    std::cout << "Iter: " << idx << std::endl;
    
    double max_dis = -INFINITY;
    
    profile *t = new profile("compute displacement");
    
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto pos = nit->get_pos();
            vec3 dis = m_vtkWrapper.get_dispacement(pos);
            
            s_dsc->set_destination(nit.key(), pos + dis);
            
            if (max_dis < dis.length())
            {
                max_dis = dis.length();
            }
        }
    }
    
    
    std::cout << "Max displacement: " << max_dis << std::endl;
    
    t->change("Load next grid");
    m_vtkWrapper.load_next_grid();
    
    t->change("displace DSC");
    s_dsc->deform();
    
    delete t;
    
    profile::close();
    
    log_dsc_surface(idx);
    
    if (idx++ > 300)
    {
        exit(0);
    }
}

void fluid_motion::log_dsc_surface(int idx)
{
    try
    {
        std::stringstream s;
        s << "LOG/dsc_" << setfill('0') << setw(5) << idx << ".suf";
        
        std::ofstream f(s.str());
        if (f.is_open())
        {
            for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
            {
                if (fit->is_interface())
                {
                    auto pos = s_dsc->get_pos(s_dsc->get_nodes(fit.key()));
                    f << pos[0] << std::endl << pos[1] << std::endl << pos[2] << std::endl;
                }
            }
        }else{
            throw "Can not open file";
        }
    }
    catch (std::exception e)
    {
        std::cout << "Error " << e.what();
    }

}
