//
//  Deformabel Simplicial Complex (DSC) method
//  Copyright (C) 2013  Technical University of Denmark
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  See licence.txt for a copy of the GNU General Public License.

#include "user_interface.h"
#include <iostream>
#include "InputParser.h"
#include "fluid_motion.hpp"
#include "debugger.h"
#include "define.h"
#include "eigen_wrapper.hpp"

#include "DSC.h"

#include "fix_DSC.h"

using namespace std;

string fluid_motion::m_data_path;

#ifdef __APPLE__
string data_path = "../Large_data";
#else
string data_path = "../../Large_data";
#endif

string problem = "two_phase_fluid";
//string problem = "DamBreak3D";

string g_out_path; // TO write the surface
double g_res; // Affect DSC resolution

void extract_surface_phase(int phase, std::string path, DSC::DeformableSimplicialComplex<> * s_dsc)
{
    bool shift = false;
    // If not shift, then the shared interface will be removed
    
    if (shift && phase == 2)
    {
        // Shrink the share interface

        double shrink = 0.001*s_dsc->get_avg_edge_length();
        vector<vec3> vertex_shrink(s_dsc->get_no_nodes(), vec3(0.0));
        vector<double> contribute(s_dsc->get_no_nodes(), (0.0));


        for (auto fit = s_dsc->faces_begin(); fit != s_dsc->faces_end(); fit++)
        {
            auto tets = s_dsc->get_tets(fit.key());
            if (s_dsc->get_label(tets[0]) != 0
                && s_dsc->get_label(tets[1]) != 0)
            {
                auto norm = s_dsc->get_normal(fit.key());
                for(auto n : s_dsc->get_nodes(fit.key()))
                {
                    vertex_shrink[n] += -norm*shrink;
                    contribute[n] += 1;
                }
            }
        }
        for (int i = 0; i < vertex_shrink.size(); i++)
        {
            if (contribute[i] > 0)
            {
                vec3 dis = vertex_shrink[i]/contribute[i];
                is_mesh::NodeKey nkey(i);
                vec3 pos = s_dsc->get_pos(nkey) + dis;
                s_dsc->set_pos(nkey, pos);
            }
        }
    }
    
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
                if (!shift
                    && phase == 2
                    && s_dsc->get_label(tets[0]) != 0
                    && s_dsc->get_label(tets[1]) != 0)
                {
                    continue;
                }
                
                auto tid = (s_dsc->get_label(tets[0]) == phase? tets[0]:tets[1]);
                
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
    
    cout << "Write to: " << path << endl;
}

void extract_2_phase_surface(string path)
{

    
    string directory = path.substr(0, path.find_last_of("\\/"));
    string name = path.substr(path.find_last_of("\\/") + 1);
    string phase[2] =  {directory + "/" + name + "_0.obj", directory +"/" + name + "_1.obj"};
    
    // Load
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(path, points, tets, tet_labels);
    
    DSC::DeformableSimplicialComplex<> dsc(points, tets, tet_labels);
    
    // Write
    extract_surface_phase(1, phase[0], &dsc);
    extract_surface_phase(2, phase[1], &dsc);
}

int main(int argc, char** argv)
{
    fix_DSC("/Users/tuannt8/Desktop/iter.dsc", 0.0029, vec3(0.15, 0.15, 0.15));
    return 0;
    
//    ///////////////////////////////////////////////
//    // Extract surface
//    string path(argv[1]);
//    extract_2_phase_surface(path);
//    return 0;
//    ///////////////////////////////////////////////
    
    InputParser input(argc, argv);
    
    fluid_motion::m_data_path = data_path + "/" + input.getCmdOption("-name", problem);
    g_out_path = input.getCmdOption("-out_path", "surface");
    
    bool nodisplay = input.cmdOptionExists("-no_display");
    g_res = atof(input.getCmdOption("-res", "3.0").c_str());
    
    input.print();
    
    if (input.cmdOptionExists("-help"))
    {
        return 0;
    }
    
    if (nodisplay)
    {
        // No display
        UI ui;
        while (1)
        {
            ui.m_fluid.deform();
        }
    }
    else
    {
    
        UI ui(argc, argv);

        glutMainLoop();
    }
    return 0;
}
