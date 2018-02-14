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
    
    m_max_dsc_displacement = s_dsc->get_avg_edge_length()*0.4;
    m_max_displacement_projection = min(m_problem->m_deltap, m_max_dsc_displacement);
    m_threshold_projection = m_max_displacement_projection*0.3;
    
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
    double max_displace = -INFINITY;
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        max_displace = std::max(max_displace, m_particles[i]->get_max_displacement());
    }
    
    cout << "Max particle displace: " << max_displace << " and dsc " << m_max_dsc_displacement << endl;
    return ceil(max_displace / m_max_dsc_displacement);
}

bool fluid_motion::load_next_particle()
{
    bool load_new = false;
    
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
        load_new = true;
        
        cur_global_idx++;
        sub_step_idx = 0;
        
        for (int i = 0; i < m_problem->m_nb_phases; i++)
        {
            m_particles[i]->load_time_step(cur_global_idx);
        }
        
        sub_step_count = subdivide_time_step();
        
        std::cout << "=====================================" << endl
        << "Load particle " << cur_global_idx
        << " and subdivide to " << sub_step_count << " sub steps." << endl;
    }
    
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->interpolate(sub_step_idx, sub_step_count);
    }
    sub_step_idx++;
    
    return load_new;
}

void fluid_motion::draw()
{
    if(glut_menu::get_state("Particles point 0", 0))
    {
        glColor3f(1, 0, 0);
        m_particles[0]->draw();
    }
    
    if(glut_menu::get_state("Particles point 1", 0))
    {
        glColor3f(0, 1, 0);
        m_particles[1]->draw();
    }
    
    if(glut_menu::get_state("anisotropic", 0))
    {
        m_particles[0]->draw_anisotropic_kernel(m_problem->domain_size());
    }
    
    if ((glut_menu::get_state("Isotropic sphere", 0)))
    {
        m_particles[1]->draw_anisotropic_kernel();
    }
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
    double max_dis = -INFINITY;

    profile *t = new profile("compute displacement");

    static double dt = m_problem->m_deltap;

    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface())
        {
            auto pos = nit->get_pos();
            vec3 dis(0.0);
            bool bFound = false;
            // get displacement from multiple phase, then average
            for (int i = 0; i < m_problem->m_nb_phases; i++)
            {
                vec3 phase_dis(0.0);
                if (m_particles[i]->get_displacement(pos, phase_dis))
                {// found
                    if (dis.length() < phase_dis.length())
                    {
                        dis = phase_dis;
                    }
                    bFound = true;
                }
            }

            if (!bFound)
            {
                // Shrink
                auto norm = s_dsc->get_normal(nit.key());
                dis = norm*(-dt);
            }

            s_dsc->set_destination(nit.key(), pos + dis);

            if (max_dis < dis.length())
            {
                max_dis = dis.length();
            }
        }
    }


    std::cout << "Max displacement from velocity projection: " << max_dis << std::endl;

    snapp_boundary_vertices();

    t->change("displace DSC");
    s_dsc->deform(20);

#ifdef __APPLE__
#else
    log_dsc_surface(idx);
#endif
    
    

    
    if(load_next_particle())
    {
        /////////////////////////////////////////////////////
        // Try pure projection
        /////////////////////////////////////////////////////
        project_interface_itteratively();
    }
    idx++;
    
    cout << "DSC statistic (nodes, edges, faces, tets): " << s_dsc->get_no_nodes() << " " << s_dsc->get_no_edges() << " " << s_dsc->get_no_faces() << " " << s_dsc->get_no_tets() << endl;
}

void fluid_motion::project_interface_itteratively(){
    
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->build_anisotropic_kernel();
    }
    
    while (1)
    {
        double max_displace = project_interface(m_threshold_projection*0.8);
        if (max_displace < m_threshold_projection )
        {
            break;
        }
    }
}

bool fluid_motion::is_boundary_work_around(is_mesh::FaceKey fkey)
{
    double thres = m_problem->m_deltap; // Should be larger than max projection
    
    static vec3 dim = m_problem->domain_size();
    static vec3 origin(0.0);
    
    auto pts = s_dsc->get_pos(s_dsc->get_nodes(fkey));
    for (int idx = 0; idx < 3; idx++)
    {
        bool is_bound_vertex = false;
        
        auto p = pts[idx];
        for (int i = 0; i < 3; i++)
        {
            if (abs(p[i] - dim[i]) < thres)
            {
                is_bound_vertex = true;
            }
            if (abs(p[i] - origin[i]) < thres)
            {
                is_bound_vertex = true;
            }
        }// Snap to boundary
        
        if (!is_bound_vertex)
        {
            return false;
        }
    }
    
    return true;
}

void fluid_motion::snapp_boundary_vertices()
{
    double thres = m_problem->m_deltap; // Should be larger than max projection
    
    double max_displace = 0;
    
    vec3 dim = m_problem->domain_size();
    vec3 origin(0.0);

    
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
            
            max_displace = std::max(max_displace, (p - nit->get_pos()).length() );
        } // If interface
    }// For all nodes
    
    cout << "Max displacement: " << max_displace << endl;
}

#define get_barry_pos(b, pos) pos[0]*b[0] + pos[1]*b[1] + pos[2]*b[2]

double fluid_motion::project_interface(double min_displace){
    vector<vec3> vertex_dis(s_dsc->get_no_nodes(), vec3(0.0));
    vector<double> contribution(s_dsc->get_no_nodes(), 0.0);
    
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            // Ignore the face on the work around boundary
            if (is_boundary_work_around(fit.key()))
            {
                continue;
            }
            
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
                   || s_dsc->get_label(tet_keys[1]) == 0) { // Single interface
                    int label = s_dsc->get_label(tet_keys[0]) == 0? s_dsc->get_label(tet_keys[1]) : s_dsc->get_label(tet_keys[0]);
                    
                    vDisplace = m_particles[label - 1]->m_aniso_kernel.get_displacement_projection(sample_pos, norm_global, m_max_displacement_projection);
                }
                else{ // Sharing interface
                    for (int i = 0; i < m_particles.size(); i++)
                    {
                        auto norm = norm_global *( i==0? -1:1);
                        vDisplace += m_particles[i]->m_aniso_kernel.get_displacement_projection(sample_pos, norm, m_max_displacement_projection);
                    }

                    vDisplace *= 0.5;
                }
                
                // Distribute
                for (int i =0; i < 3; i++)
                {
                    vertex_dis[node_pts[i]] += vDisplace*sp[i];
                    contribution[node_pts[i]] += sp[i];
                }
            } // If interface
        } // For all triangles
    }
    
    // Normalize the displacement
    for (int i = 0; i < vertex_dis.size(); i++){
        if (contribution[i] > 0){
            vertex_dis[i] /= contribution[i];
        }
    }
    
    // Set destination and deform the DSC
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface()
            && vertex_dis[nit.key()].length() > min_displace)
        {
            s_dsc->set_destination(nit.key(), nit->get_pos() + vertex_dis[nit.key()]);
        }
    }
    
    snapp_boundary_vertices();
    
    double max_displace = 0;
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++){
        if (nit->is_interface()){
            max_displace = max(max_displace,
                               (nit->get_destination() -nit->get_pos()).length());
        }
    }
    
    s_dsc->deform();
    
    cout << "Max displace during projeciton: " << max_displace << endl;
    return max_displace;
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
