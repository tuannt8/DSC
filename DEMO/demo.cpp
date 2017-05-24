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

// 0: No interface, for linux cluster; 1: Use interface
#define _USE_INTERFACE 1


#if _USE_INTERFACE == 0
int num_iters = 1000;
int main(int argc, char** argv)
{
    UI ui;
    for(int i = 0; i < num_iters; i++)
    {
        ui._seg.segment();
    }
}
#else
int main(int argc, char** argv)
{
    UI ui(argc, argv);
//    printMenu();
    glutMainLoop();
    return 0;
}
#endif
