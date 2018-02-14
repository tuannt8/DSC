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

int main(int argc, char** argv)
{
    
//    {
//        mat3x3d Q1(vec3(0,0,1),
//                   vec3(1,0,0),
//                   vec3(0,1,0));
//        cout << "Q1: " << Q1;
//        mat3x3d sigma(vec3(9,0,0),
//                      vec3(0,9,0),
//                      vec3(0,0,22));
//        mat3x3d Q2(vec3(0,0,1),
//                   vec3(0,1,0),
//                   vec3(1,0,0));
//        cout << "Q2: " << Q2;
//        auto G1 = Q1*sigma*CGLA::transpose(Q1);
//        auto G2 = Q2*sigma*CGLA::transpose(Q2);
//        cout << "G1: " << G1;
//        cout << "G2: " << G2;
//        
//        cout << G1 - G2;
//    }
//    return 0;
    
    InputParser input(argc, argv);
    
    fluid_motion::m_data_path = data_path + "/" + input.getCmdOption("-path", "two_phase_fluid");
    
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
