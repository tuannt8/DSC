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

using namespace std;

string fluid_motion::m_data_path;

#ifdef __APPLE__
string data_path = "../Large_data";
#else
string data_path = "../../Large_data";
#endif

string problem = "two_phase_fluid";
//string problem = "DamBreak3D";

int main(int argc, char** argv)
{
    
    InputParser input(argc, argv);
    
    fluid_motion::m_data_path = data_path + "/" + input.getCmdOption("-path", problem);
    
    if (argc > 1)
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
