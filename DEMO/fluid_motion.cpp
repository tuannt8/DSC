//
//  fluid_motion.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "fluid_motion.hpp"
#include <iomanip>

using namespace std;

fluid_motion::fluid_motion()
{
    
}

void fluid_motion::draw()
{
    m_file_load.draw();
}

void fluid_motion::deform()
{
//    project_interface();
//    return;
//
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
            vec3 dis;
            if(!m_file_load.get_displacement(pos, dis))
            {
                // invert normal
                auto norm = s_dsc->get_normal(nit.key());
                dis = norm*m_file_load.get_spacing_distance()*DT_NORM;
            }

            s_dsc->set_destination(nit.key(), pos + dis);
            
            if (max_dis < dis.length())
            {
                max_dis = dis.length();
            }
        }
    }


    std::cout << "Max displacement: " << max_dis << std::endl;



    t->change("displace DSC");
    s_dsc->deform();
    
    
    ///////////////////////////////////////////
    //  2. Project interface
    t->change("Project DSC");
    project_interface();
    
    t->change("Load next grid");
    m_file_load.load_time_step();
    
    delete t;

    profile::close();

    log_dsc_surface(idx);
}

void fluid_motion::project_interface()
{
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto norm = s_dsc->get_normal(nit.key());
            bool bInside;
            double t = 0;
            if(!m_file_load.get_projection(nit->get_pos(), norm, bInside, t))
            {
                 t = m_file_load.get_spacing_distance()*(bInside? 1:-1);
            }
            vec3 new_pos = nit->get_pos() + norm*t;
            s_dsc->set_destination(nit.key(), new_pos);
        }
    }
    
    s_dsc->deform();
}

void fluid_motion::log_dsc_surface(int idx)
{
    try
    {
        std::stringstream s;
        s << "LOG/dsc_" << setfill('0') << setw(5) << idx << ".obj";
        
        std::vector<vec3> points;
        std::vector<int> faces;
        s_dsc->extract_surface_mesh(points, faces);
        is_mesh::export_surface_mesh(s.str(), points, faces);
    }
    catch (std::exception e)
    {
        std::cout << "Error " << e.what();
    }

}
