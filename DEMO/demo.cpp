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

int main(int argc, char** argv)
{
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
