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

#include "export.hpp"

using namespace std;

string fluid_motion::m_data_path;

double smooth_ratio;


double scale_anisotropic = 1.0;

#ifdef __APPLE__
string data_path = "../Large_data";
#else
string data_path = "../../Large_data";
#endif

string problem = "two_phase_fluid";
//string problem = "DamBreak3D";

string g_out_path; // TO write the surface
double g_res; // Affect DSC resolution


int main(int argc, char** argv)
{
    #ifdef __APPLE__
//    fix_DSC("/Users/tuannt8/Desktop/iter.dsc", 0.0029, vec3(0.15, 0.15, 0.15));
//    return 0;
#endif
    
//    ///////////////////////////////////////////////
//    // Extract surface
//    string path(argv[1]);
//    extract_2_phase_surface(path);
//    return 0;
//    ///////////////////////////////////////////////
    
    InputParser input(argc, argv);
    
    /////////////////////////////////
    // Export surface
    bool export_mesh = input.cmdOptionExists("-export");
    if(export_mesh)
    {
        string path = input.getCmdOption("-export");
        if (input.cmdOptionExists("-two-phases"))
        {
            dsc_export::export_two_phase_fluid(path);
        }else
            dsc_export::export_surface(path);
        return 0;
    }
    ///////////////////////////////////////////////
    
    fluid_motion::m_data_path = data_path + "/" + input.getCmdOption("-name", problem);
    g_out_path = input.getCmdOption("-out_path", "surface");
    smooth_ratio = atof(input.getCmdOption("-smooth", "0.3").c_str());
    scale_anisotropic = atof(input.getCmdOption("aniso-r", "1.0").c_str());
    
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
//        ui.export_dam_break();
        glutMainLoop();
    }
    return 0;
}
