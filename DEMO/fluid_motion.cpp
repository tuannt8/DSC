//
//  fluid_motion.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "fluid_motion.hpp"
#include <iomanip>
#include "eigen_wrapper.hpp"

using namespace std;

extern string g_out_path;

fluid_motion::fluid_motion()
{
    
}

string find_name(string input)
{
    auto pos = input.find_last_of("/\\");
    
    return input.substr(pos+1);
}

void fluid_motion::init(DSC::DeformableSimplicialComplex<> *dsc){
    s_dsc = dsc;
    m_max_dsc_displacement = //m_problem->m_deltap;
                                 max(s_dsc->get_avg_edge_length(), m_problem->m_deltap);
    
    m_max_displacement_projection = m_problem->m_deltap;
    m_threshold_projection = m_problem->m_deltap*0.5;
    
    cout << "\n\n+++++++++++++++++++++++++++++++++++++++"
    << "\n spacing distance: " << m_problem->m_deltap
    << "\n smoothing length: " << m_problem->m_slength
    << "\n average lenght: " << s_dsc->get_avg_edge_length()
    << "\n Initialize paratmeters"
    << "\n Max DSC displacement: " << m_max_dsc_displacement
    << "\n Projection search: " << m_max_displacement_projection
    << "\n thres hold for projection: " << m_threshold_projection
    << "\n++++++++++++++++++++++++++++++++++++++++++\n\n";
}
void fluid_motion::load_configuration()
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
    else if(problem_name.compare("DamBreak3D")==0){
        m_problem = std::unique_ptr<problem>(new dam_break_fluid);
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
        m_out_path.push_back(m_data_path + "/" + g_out_path + "_" + std::to_string(i));
        create_directory(m_out_path[i].c_str());
        cout << "Output files: " << m_out_path[i] << endl;
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

void fluid_motion::load_first_particle()
{
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->load_first_time(m_cur_global_idx);
    }
    
    m_sub_step_count = subdivide_time_step();
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->interpolate(m_sub_step_idx, m_sub_step_count);
    }
    
    m_cur_global_idx ++;
}

void fluid_motion::load_next_particle()
{
    cout << "===========================\nparticle + sub/sum: " << m_cur_global_idx << " + " << m_sub_step_idx << "/" << m_sub_step_count << endl;
    
    if(m_sub_step_idx == m_sub_step_count)
    {
        for (int i = 0; i < m_problem->m_nb_phases; i++)
        {
            m_particles[i]->load_time_step(m_cur_global_idx);
        }
        
        m_sub_step_count = subdivide_time_step();
        
        m_cur_global_idx++;
        m_sub_step_idx = 0;
    }
    
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->interpolate(m_sub_step_idx, m_sub_step_count);
    }
    m_sub_step_idx++;
}

void fluid_motion::draw()
{
    if(glut_menu::get_state("Particles point 0", 0))
    {
        glColor3f(1, 0, 0);
        m_particles[0]->draw(m_problem->domain_size()[1]*0.32, m_problem->domain_size()[1]*0.35);
    }
    
    if(glut_menu::get_state("Particles point 1", 0))
    {
        glColor3f(0, 1, 0);
        m_particles[1]->draw(m_problem->domain_size()[1]*0.32, m_problem->domain_size()[1]*0.35);
    }
    
    if(glut_menu::get_state("anisotropic", 0))
    {
        m_particles[0]->draw_anisotropic_kernel(m_problem->domain_size(), vec3(1,0,0));
        m_particles[1]->draw_anisotropic_kernel(m_problem->domain_size(), vec3(0,1,1));
    }
    
    if ((glut_menu::get_state("Isotropic sphere", 0)))
    {
        m_particles[1]->draw_anisotropic_kernel(m_problem->domain_size()[1]*0.32, m_problem->domain_size()[1]*0.35);
    }
}

void fluid_motion:: advect_velocity()
{
    // 1. Interpolate the displacement
    update_vertex_boundary();
    
    
    ///////////////////////////////////////////////////////////////////////////
    // Advect velocity does not work
    // Because of fluid convection and advection, internal particles contribute
    //  false velocities
    ///////////////////////////////////////////////////////////////////////////
    double max_dis = -INFINITY;
    
    profile *t = new profile("compute displacement");
    
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
            
            auto norm = s_dsc->get_normal(nit.key());
            if (!bFound)
            {
                // Shrink
                dis = norm*(-m_threshold_projection);
            }
            
            s_dsc->set_destination(nit.key(), pos + dis);
            
            if (max_dis < dis.length())
            {
                max_dis = dis.length();
            }
        }
    }
    
    std::cout << "Max advect displacement: " << max_dis << std::endl;
    
    snapp_boundary_vertices();
    
    t->change("displace DSC");
    s_dsc->deform(20);
    

}
void fluid_motion::deform()
{
//    if(load_next_particle())
    {
//        project_interface_itteratively();
        
//        s_dsc->print_mesh_info();
    }
    
    advect_velocity();
    project_interface_one_iter();
    
//    project_interface_itteratively();
    
    
    load_next_particle();
    
    log_dsc_surface();
    
    assert(s_dsc->get_no_faces_buffer() < 10000000);// assert cache
}

void fluid_motion::project_interface_one_iter()
{
    cout << "\n------ Project interface ----------- \n";
    
    update_vertex_boundary();
    
    // Build anisotropic kernel
    for (int i = 0; i < m_problem->m_nb_phases; i++)
    {
        m_particles[i]->build_anisotropic_kernel();
    }
    
    project_interface();
    
//    // Mark shared nodes
//    //  Will be project to smaller index fluid
//    vector<int> correspond_fluid(s_dsc->get_no_nodes_buffer(), -1);
//    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
//    {
//        if(fit->is_interface())
//        {
//            auto tets = s_dsc->get_tets(fit.key());
//            auto l0 = s_dsc->get_label(tets[0]);
//            auto l1 = s_dsc->get_label(tets[0]);
//            auto higher_label = l0>l1? l0 : l1;
//
//            for(auto n : s_dsc->get_nodes(fit.key())){
//                correspond_fluid[n] = higher_label - 1;
//            }
//        }
//    }
//
//    double max_projection = 0;
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_interface()
//            && correspond_fluid[nit.key()] != -1
//            )
//        {
//            if (is_vertices_boundary[nit.key()] && !nit->is_crossing())
//            {
//                continue;
//            }
//
//            auto pos = nit->get_pos();
//            auto norm = s_dsc->get_normal(nit.key()); // from high to low,
//            auto dis = m_particles[correspond_fluid[nit.key()]]->m_aniso_kernel.get_displacement_projection(pos, norm, m_max_displacement_projection);
//
//            s_dsc->set_destination(nit.key(), pos + dis);
//
//            max_projection = max(max_projection, dis.length());
//        }
//    }
//
//    cout << "Max projection: " << max_projection << endl;
//
//    snapp_boundary_vertices();
//
//    s_dsc->deform(20);
}

void fluid_motion::reset_projected_flag()
{
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        fit->set_projected(false);
    }
}

void fluid_motion::project_interface_itteratively(){
    
    reset_projected_flag();
    
    cout << "\n--------------------------------\n Start projecting surface "<<endl;
    
//    for (int i = 0; i < m_problem->m_nb_phases; i++)
//    {
//        m_particles[i]->build_anisotropic_kernel();
//    }
    
    // Becasue we already restrict the particle displacement
    // DSC should not move so many iterations
    for (int idx = 0; idx < 2; idx++)
    {
//        update_vertex_boundary();
//        double max_displace = project_interface(m_threshold_projection*0.5);
//
//        if (max_displace < m_threshold_projection )
//        {
//            break;
//        }
        
        project_interface_one_iter();
        
//        cout << "idx " << idx << " max displacement: " << max_displace <<"/" << m_threshold_projection << endl << endl;
    }
}

bool fluid_motion::is_boundary_work_around(is_mesh::FaceKey fkey)
{
//    double thres = m_problem->m_deltap; // Should be larger than max projection
//
//    static vec3 dim = m_problem->domain_size();
//    static vec3 origin(0.0);
    
    auto nodes = s_dsc->get_nodes(fkey);
    for (int idx = 0; idx < 3; idx++)
    {
        if (!is_vertices_boundary[nodes[idx]])
        {
            return false;
        }
    }
    
    return true;
}

void fluid_motion::snapp_boundary_vertices()
{
    double thres = m_max_dsc_displacement; // Should be larger than max projection
    
    double max_displace = 0;
    
    vec3 dim = m_problem->domain_size();
    vec3 origin(0.0);
    
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (is_vertices_boundary[nit.key()])
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
}

#define get_barry_pos(b, pos) pos[0]*b[0] + pos[1]*b[1] + pos[2]*b[2]

double fluid_motion::project_interface()
{
    vector<vec3> vertex_dis(s_dsc->get_no_nodes_buffer(), vec3(0.0));
    vector<double> contribution(s_dsc->get_no_nodes_buffer(), 0.0);
    
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_interface()
            && !fit->is_projected()
            && !is_boundary_work_around(fit.key())
            )
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
                if(s_dsc->get_label(tet_keys[0]) == 0
                   || s_dsc->get_label(tet_keys[1]) == 0) { // Single interface
                    int label = s_dsc->get_label(tet_keys[0]) == 0? s_dsc->get_label(tet_keys[1]) : s_dsc->get_label(tet_keys[0]);
                    
                    bool bLast = false;
                    vDisplace = m_particles[label - 1]->m_aniso_kernel.get_displacement_projection(sample_pos, norm_global, m_max_displacement_projection, bLast);
                    
                    fit->set_projected(bLast);
                }
                else{ // Sharing interface
//                    for (int i = 0; i < m_particles.size(); i++)
//                    {
//                        //  tuannt8: HARDCODE for 2 phases
//                        auto norm = norm_global *( i==0? -1:1);
//                        vDisplace += m_particles[i]->m_aniso_kernel.get_displacement_projection(sample_pos, norm, m_max_displacement_projection);
//                    }
                    
                    // prior smaller label
                    int i = s_dsc->get_label(tet_keys[0]) < s_dsc->get_label(tet_keys[1])? s_dsc->get_label(tet_keys[0]) : s_dsc->get_label(tet_keys[1]);
                    
                    auto norm = norm_global *( i==0? -1:1);
                    bool bLast = false;
                    vDisplace += m_particles[i]->m_aniso_kernel.get_displacement_projection(sample_pos, norm, m_max_displacement_projection, bLast);
                    fit->set_projected(bLast);
//                    vDisplace *= 0.5;
                    
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
    int nb_move = 0;
    for (int i = 0; i < vertex_dis.size(); i++){
        if (contribution[i] > 0){
            vertex_dis[i] /= contribution[i];
//            if (vertex_dis[i].length() > min_displace)
            {
                nb_move++;
            }
        }
    }
    
    
    int nb_projected = 0;
    for(auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_projected())
        {
            nb_projected++;
        }
    }
    cout << nb_projected << " face projected" << endl;
//    if (nb_move < 10)
//    {
//        return 0;// Ignore if there is only few node displacement
//        // There particles may be the turbulence of the fluid
//    }
    
    // Set destination and deform the DSC
    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
    {
        if (nit->is_interface()
//            && vertex_dis[nit.key()].length() > min_displace
            )
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
    
    s_dsc->deform(10);

    return max_displace;
}

void fluid_motion::log_dsc_surface()
{
    static int idx = 0;
    try
    {
        cout << "Log surface " << idx << endl;
        for (int i = 0; i < m_problem->m_nb_phases; i++)
        {
            std::stringstream s;
            s << m_out_path[i] << "/" << setfill('0') << setw(5) << idx << ".obj";
            extract_surface_phase(i+1, s.str());
        }
//
//        std::stringstream s;
//        s << "LOG/dsc_" << setfill('0') << setw(5) << idx << ".obj";
//
//        std::vector<vec3> points;
//        std::vector<int> faces;
//        s_dsc->extract_surface_mesh(points, faces);
//        is_mesh::export_surface_mesh(s.str(), points, faces);
    }
    catch (std::exception e)
    {
        std::cout << "Error " << e.what();
    }

    idx ++;
}

inline bool is_under_out(vec3 p, vec3 bound)
{
    return p[0] < bound[0]
        || p[1] < bound[1]
    || p[2] < bound[2];
}

inline bool is_upper_out(vec3 p, vec3 bound)
{
    return p[0] > bound[0]
    || p[1] > bound[1]
    || p[2] > bound[2];
}

inline bool is_bound_point(vec3 & p, vec3 & origin, vec3 & domain_size)
{
    return (is_under_out(p, origin) || is_upper_out(p, domain_size));
}

void fluid_motion::update_vertex_boundary()
{
    is_vertices_boundary = vector<bool>(s_dsc->get_no_nodes_buffer(), false);
//    for (auto nit = s_dsc->nodes_begin(); nit != s_dsc->nodes_end(); nit++)
//    {
//        if (nit->is_boundary())
//        {
//            auto neighbor = s_dsc->get_nodes(s_dsc->get_edges(nit.key()));
//            for (int i =0; i < neighbor.size(); i++)
//            {
//                if (s_dsc->get(neighbor[i]).is_interface())
//                {
//                    is_vertices_boundary[neighbor[i]] = true;
//                }
//            }
//        }
//    }
    
    static vec3 origin(0);
    static vec3 domain_size = m_problem->domain_size();
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            auto nodes = s_dsc->get_nodes(fit.key());
            auto tets = s_dsc->get_tets(fit.key());
            auto p0 = s_dsc->get_pos(s_dsc->get_nodes(tets[0]) - nodes);
            auto p1 = s_dsc->get_pos(s_dsc->get_nodes(tets[1]) - nodes);
            
            bool b_is_bound = is_bound_point(p0[0], origin, domain_size) || is_bound_point(p1[0], origin, domain_size);
            
            if (b_is_bound)
            {
                for (auto n : nodes)
                {
                    is_vertices_boundary[n] = true;
                }
            }
        }
    }
}

void fluid_motion::extract_surface_phase(int phase, std::string path)
{
    vector<int> indices_map(s_dsc->get_no_nodes_buffer(), -1);
    int idx = 0;
    
    // Write face first
    stringstream vertices_write, faces_write;
    for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
    {
        if (fit->is_interface())
        {
            auto tets = s_dsc->get_tets(fit.key());
            if(s_dsc->get_label(tets[0]) == phase
               || s_dsc->get_label(tets[1]) == phase)
            {
                auto tid = s_dsc->get_label(tets[0]) == phase? tets[0]:tets[1];
                
                auto nodes = s_dsc->get_sorted_nodes(fit.key(), tid);
                faces_write << "f ";
                for (int i = 0; i < 3; i++)
                {
                    auto n = nodes[i];
                    if (indices_map[n] == -1)
                    {
                        indices_map[n] = idx++;
                    
                        auto pos = s_dsc->get(n).get_pos();
                        vertices_write << "v " << pos[0] << " " << pos[1] << " " << pos[2] << endl;
                    }
                    
                    faces_write << indices_map[n] + 1 << " ";
                }
                faces_write << endl;
            }
        }
    }
    
    ofstream of(path);
    of << vertices_write.str();
    of << faces_write.str();
    of.close();
}
