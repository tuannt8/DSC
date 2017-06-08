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
using namespace std;

// Will display OpenGL window for visualization
bool arg_b_display = true;
string arg_b_no_display_command = "-no_display";

// Max iteration to run, in case we dont have a display
int arg_i_num_iters = 1000;
string arg_i_num_iters_command = "-max_iter";

// Flag: Use original DSC in build table function
bool arg_b_build_table_origin = true;
string arg_b_build_table_origin_command = "-build_table";

string config_file = "fuel_cells_smaller.properties";
//string config_file = "square_arg.properties";
//string config_file = "square_round_arg.properties";
string config_file_command = "-config_file";


void parse_arg(int argc, char** argv)
{
    try
    {
        for (int i = 1; i < argc; i++)
        {
            if (strcmp(argv[i], arg_b_no_display_command.c_str())==0)
            {
                arg_b_display = false;
            }
            else if (strcmp(argv[i], arg_i_num_iters_command.c_str())==0)
            {
                arg_i_num_iters = atoi(argv[++i]);
            }
            else if (strcmp(argv[i], arg_b_build_table_origin_command.c_str())==0)
            {
                arg_b_build_table_origin = false;
            }
            else if (strcmp(argv[i], config_file_command.c_str())==0)
            {
                config_file = atoi(argv[++i]);;
            }
        }
    }
    catch (exception e)
    {
        cout << "Error in input arguments";
        exit(1);
    }
}

int main(int argc, char** argv)
{
 
    parse_arg(argc, argv);
    
    if(arg_b_display)
    {
        UI ui(argc, argv);
        glutMainLoop();

    }
    else{
        UI ui;
        for(int i = 0; i < arg_i_num_iters; i++)
        {
            ui._seg.segment();
        }
    }
    
    return 0;
}
