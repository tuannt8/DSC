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

string find_name(string input)
{
    auto pos = input.find_last_of("/\\");
    
    return input.substr(pos+1);
}

void fluid_motion::init()
{
    // Make sure the data path dont have '/' at last
    if(m_data_path.back() == '/')
        m_data_path.pop_back();
    
    
    // 1. Init problem name
    string problem_name = find_name(m_data_path);
    if (problem_name.compare("two_phase_fluid")==0)
    {
        m_problem = std::unique_ptr<problem>(new two_phase_fluid);
    }
    else
    {
        throw invalid_argument("Problem name unmatched");
    }
    
    // 2. Load problem parameters
    m_problem->init(m_data_path + "/summary.txt");
    
    // 3. Init particles for all phases
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        auto new_particle_mng = std::shared_ptr<particle_manager>(new particle_manager);
        // Set data
        new_particle_mng->m_data_path = m_data_path + "/my_format_" + std::to_string(i);
        new_particle_mng->m_influence_radius = m_problem->m_influenceRadius;
        new_particle_mng->m_deltap = m_problem->m_deltap;
        new_particle_mng->m_slength = m_problem->m_slength;
        
        m_particles.push_back(new_particle_mng);
    }
    
    // 4. Set up output
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_out_path.push_back(m_data_path + "/surface_" + std::to_string(i));
    }
}

int fluid_motion::subdivide_time_step()
{
    double max_dsc_displacement = m_problem->m_deltap; // Beware of this parametter
    double max_displace = -INFINITY;
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        max_displace = std::max(max_displace, m_particles[i]->get_max_displacement());
    }
    
    return ceil(max_displace / max_dsc_displacement);
}

void fluid_motion::load_next_particle()
{
    static int cur_global_idx = -1;
    static int sub_step_idx = 0;
    static int sub_step_count = 0;
    
    if (cur_global_idx == -1)//first load
    {
        for (int i = 0; i < m_problem->m_nb_phases; i++)
        {
            m_particles[i]->load_time_step(-1);
        }
    }
    
    if(sub_step_idx == sub_step_count)
    {
        cur_global_idx++;
        sub_step_idx = 0;
        
        for (int i = 0; i < m_problem->m_nb_phases; i++)
        {
            m_particles[i]->load_time_step(cur_global_idx);
        }
        
        sub_step_count = subdivide_time_step();
    }
    
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->interpolate(sub_step_idx, sub_step_count);
    }
    sub_step_idx++;
}

void fluid_motion::draw()
{
    static vector<vec3> color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
    int i = 0;
    for (auto pp : m_particles)
    {
        auto & c = color[i++];
        glColor3f(c[0], c[1], c[2]);
        pp->draw();
    }
}

void fluid_motion::deform()
{

//    // 1. Interpolate the displacement
//    static int idx = 0;
//    std::cout << "Iter: " << idx << std::endl;
//
//    double max_dis = -INFINITY;
//
//    profile *t = new profile("compute displacement");
//
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_interface())
//        {
//            auto pos = nit->get_pos();
//            vec3 dis(0.0);
//            if(!m_file_load.get_displacement(pos, dis))
//            {
//                // cannot find displacement, then shrink the surface
//                // invert normal
//                auto norm = s_dsc->get_normal(nit.key());
//                dis = norm*m_file_load.get_spacing_distance()*DT_NORM;
//            }
//
//            s_dsc->set_destination(nit.key(), pos + dis);
//
//            if (max_dis < dis.length())
//            {
//                max_dis = dis.length();
//            }
//        }
//    }
//
//
//    std::cout << "Max displacement: " << max_dis << std::endl;
//
//
//
//    t->change("displace DSC");
//    s_dsc->deform();
//
//
//    ///////////////////////////////////////////
//    //  2. Project interface
//    if(idx%2==0)
//    {
//        t->change("Project DSC");
//        project_interface();
//    }
//    t->change("Load next grid");
//    m_file_load.load_time_step();
//
//    delete t;
//
//    profile::close();
//
//    log_dsc_surface(idx);
//
//    idx++;
}

void fluid_motion::project_interface()
{
//    std::cout << "Project interface \n";
//    double max_projection = 0;
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_interface())
//        {
//            auto norm = s_dsc->get_normal(nit.key());
//            bool bInside;
//            double t = 0;
//            if(!m_file_load.get_projection(nit->get_pos(), norm, bInside, t))
//            {
//                 t = m_file_load.get_spacing_distance()*(bInside? 1:-1);
//            }
//            vec3 new_pos = nit->get_pos() + norm*t;
//            s_dsc->set_destination(nit.key(), new_pos);
//
//            if(max_projection < std::abs(t))
//            {
//                max_projection = t;
//            }
//        }
//    }
//
//    std::cout << "Max projection: " << max_projection << std::endl;
//
//    s_dsc->deform();
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
