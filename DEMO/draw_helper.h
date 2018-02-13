//
//  draw_helper.hpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef draw_helper_hpp
#define draw_helper_hpp

#ifdef _WIN32 // WINDOWS
#include <GL/glut.h>
#include <GL/glew.h>
#elif defined(__APPLE__) // IOS
#include <OpenGL/gl3.h>
#include <GLUT/glut.h>
#else // LINUX
#include <GL/glew.h>
#include <GL/glut.h>
#endif


#include <stdio.h>
#include "DSC.h"

#define DISPLAY_LIM  0.05

class draw_helper
{
    typedef DSC::DeformableSimplicialComplex<> dsc_class;
public:
    /*
     GL function
     */

    /*
     Draw DSC
     */
    static void dsc_draw_edge(dsc_class & dsc);
    static void dsc_draw_domain(dsc_class & dsc);
    static void dsc_draw_interface(dsc_class & dsc, std::vector<double> * color = nullptr);
    static void dsc_draw_interface_edge(dsc_class & dsc);
    static void dsc_draw_face_norm(dsc_class & dsc);
    static void dsc_draw_shared_interface(dsc_class & dsc);
    
    static void dsc_draw_edges_colors(dsc_class & dsc);
    static void dsc_draw_node_color(dsc_class & dsc);
    
    static void save_painting(int WIDTH, int HEIGHT, std::string folder = std::string("LOG"));
    
    /*
     Draw 3 axis X - Y - Z
     */
    static void draw_coord(float length);
    
private:
    CGLA::Vec3i _cur_cross_poss = CGLA::Vec3i(0,0,0);
    
    /*
     For singleton - start
     */
public:
    static draw_helper & get_instance()
    {
        static draw_helper instance;
        return instance;
    }
    
    draw_helper(draw_helper const&) = delete;
    void operator = (draw_helper const&) = delete;
private:
    draw_helper(){};
    /* For singleton - end */
};

#endif /* draw_helper_hpp */
