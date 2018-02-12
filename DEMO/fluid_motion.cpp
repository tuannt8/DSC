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
    double max_dsc_displacement = s_dsc->get_avg_edge_length(); // Beware of this parametter
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
//    static vector<vec3> color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
//    int i = 0;
//    for (auto pp : m_particles)
//    {
//        auto & c = color[i++];
//        glColor3f(c[0], c[1], c[2]);
////        pp->draw();
//
//        pp->draw_anisotropic_kernel();
//    }
    
    m_particles[1]->draw_anisotropic_kernel();
    
//    m_particles[1]->draw_anisotropic_kernel(m_problem->domain_size());
}

void fluid_motion::deform()
{
    
    // 1. Interpolate the displacement
    static int idx = 0;
    std::cout << "Iter: " << idx << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    // Advect velocity does not work
    // Because of fluid convection and advection, internal particles contribute
    //  false velocities
    ///////////////////////////////////////////////////////////////////////////
//    double max_dis = -INFINITY;
//
//    profile *t = new profile("compute displacement");
//
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_interface() || nit->is_crossing())
//        {
//            auto pos = nit->get_pos();
//            vec3 dis(0.0);
//            // get displacement from multiple phase, then average
//            for (int i = 0; i < m_problem->m_nb_phases; i++)
//            {
//                vec3 phase_dis(0.0);
//                if (!m_particles[i]->get_displacement(pos, phase_dis))
//                {
//                    // if can not find closed particle from any phase
//                    //      should be collaped
//                    // TODO
//                }
//                dis += phase_dis;
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
    
    
    
    /////////////////////////////////////////////////////
    // Try pure projection
    /////////////////////////////////////////////////////
    project_interface();
    // work around DSC boundary
    //  Displacement is capped, so vertices close to boundary will be snapped
    snapp_boundary_vertices();
    s_dsc->deform();
    
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
    
    load_next_particle();
    idx++;
}

void fluid_motion::snapp_boundary_vertices()
{
    vec3 dim = m_problem->domain_size();
    vec3 origin(0.0);
    double thres = m_problem->m_deltap*2; // Should be larger than max projection
    
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface() || nit->is_crossing())
        {
            auto p = nit->get_destination();
            for (int i = 0; i < 3; i++)
            {
                if (abs(p[i] - dim[i]) < thres)
                {
                    p[i] = dim[i];
                }
                if (abs(p[i] - origin[i]) < thres)
                {
                    p[i] = origin[i];
                }
            }// Snap to boundary
            s_dsc->set_destination(nit.key(), p);
        } // If interface
    }// For all nodes
}

#define get_barry_pos(b, pos) pos[0]*b[0] + pos[1]*b[1] + pos[2]*b[2]

void fluid_motion::project_interface()
{
    // Project the DSC interface to particles field
    static const double shrink_vel = m_problem->m_deltap;
    
    vector<vec3> vertex_dis(s_dsc->get_no_nodes(), vec3(0.0));
    
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            static const vector<vec3> sampling_point = {vec3(0.33, 0.33, 0.33)};
            auto node_pts = s_dsc->get_nodes(fit.key());
            auto node_pos = s_dsc->get_pos(node_pts);
            auto tet_keys = s_dsc->get_tets(fit.key());
            
            auto norm_global = s_dsc->get_normal(fit.key()); // From large label to smaller label
            
            for (auto const &sp : sampling_point)
            {
                // Find the point displacement
                vec3 sample_pos = get_barry_pos(sp, node_pos);
                
                vec3 vDisplace(0.0);

                // Have to project on both particles
                //  tuannt8: HARDCODE for 2 phases
                if(s_dsc->get_label(tet_keys[0]) == 0
                   || s_dsc->get_label(tet_keys[1]) == 0)
                {
                    int label = s_dsc->get_label(tet_keys[0]) == 0? s_dsc->get_label(tet_keys[1]) : s_dsc->get_label(tet_keys[0]);
                    
                    double t = 0;
                    bool bInside;
                    if(m_particles[label-1]->get_projection(sample_pos, norm_global, bInside, t))
                    {
                        vDisplace = norm_global*t;
                    }
                    else{
                        vDisplace += norm_global*(shrink_vel)*(bInside? 1:-1);
                    }
                }
                else
                {
                   double count = 0;
                    
                    for (int i = 0; i < m_particles.size(); i++)
                    {
                        auto norm = norm_global *( i==0? -1:1);
                        double t = 0;
                        bool bInside;
                        
                        if(m_particles[i]->get_projection(sample_pos, norm, bInside, t))
                        {
                            vDisplace += norm*t;
                            count += 1;
                        }
                    }
                    
                    if (count > 0)
                    {
                        vDisplace /= count;
                    }
                }
                
                // Distribute
                for (int i =0; i < 3; i++)
                {
                    vertex_dis[node_pts[i]] += vDisplace*sp[i];
                }
            } // If interface
        } // For all triangles
    }
    
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface() || nit->is_crossing())
        {
            s_dsc->set_destination(nit.key(), nit->get_pos() + vertex_dis[nit.key()]);
        }
    }
    
//    double shrink_vel = m_problem->m_deltap;
//
//    std::cout << "Project interface \n";
//    double max_projection = 0;
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_interface() || nit->is_crossing())
//        {
//            auto norm_global = s_dsc->get_normal(nit.key()); // From large label to smaller label
//
//            vec3 vDisplace(0.0);
//            int iCount = 0;
//            bool inside = true;
//
//            auto tid = s_dsc->get_tets(<#const is_mesh::FaceKey &fid#>)
//
//            for (int i = 0; i < m_particles.size(); i++)
//            {
//                auto norm = norm_global;
//
//                double t = 0;
//                bool bInside;
//                if(m_particles[i]->get_projection(nit->get_pos(), norm, bInside, t))
//                {
//                    iCount++;
//                    vDisplace += norm*t;
//                }
//                else{
//                    vDisplace += norm*(shrink_vel)*(inside? 1:-1);
//                }
//            }
//
//            max_projection = max(max_projection, vDisplace.length());
//
//            s_dsc->set_destination(nit.key(), nit->get_pos() + vDisplace);
//        }
//    }
//    std::cout << "Max projection: " << max_projection << std::endl;
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
