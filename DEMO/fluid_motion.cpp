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
            vec3 dis = m_file_load.get_displacement(pos);

            s_dsc->set_destination(nit.key(), pos + dis);

            if (max_dis < dis.length())
            {
                max_dis = dis.length();
            }
        }
    }


    std::cout << "Max displacement: " << max_dis << std::endl;

    t->change("Load next grid");
    m_file_load.load_time_step();

    t->change("displace DSC");
    s_dsc->deform();
    
    
    ///////////////////////////////////////////
    //  2. Project interface
    t->change("Project DSC");
    project_interface();
    delete t;

    profile::close();

    log_dsc_surface(idx);

    if (idx++ > 300)
    {
        exit(0);
    }
}

void fluid_motion::project_interface()
{
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto norm = s_dsc->get_normal(nit.key());
            bool bInside;
            double t;
            if(m_file_load.get_projection(nit->get_pos(), norm, bInside, t))
            {
                // project ok
                vec3 new_pos = nit->get_pos() + norm*t;
                s_dsc->set_destination(nit.key(), new_pos);
                
//                std::cout << "Project a vertex ---" << std::endl;
            }
            else
            {
                // can not project
//                std::cout << "Can not project a vertex" << std::endl;
//                assert(0);
                // Apply normal displacement
                
            }
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
