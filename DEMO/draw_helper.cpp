//
//  draw_helper.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "draw_helper.h"
#include "DSC.h"
#include <iostream>
#include <string>
#include <SOIL/SOIL.h>
#include <fstream>



void draw_helper::dsc_draw_edge(dsc_class &dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        auto pos = dsc.get_pos(dsc.get_nodes(eit.key()));
        glVertex3dv(pos[0].get());
        glVertex3dv(pos[1].get());
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

std::vector<vec3> get_random_color(int N = 5)
{
    std::vector<vec3> cc = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1), vec3(1, 1, 0), vec3(1,0,1), vec3(0,1,1),};
    size_t Nc = cc.size();
    double dl = 1.0/(double)N;
    for (int i = 0; i < N; i++)
    {
        double co = (i+1)*dl;
        for (int j = 0; j < Nc; j++)
        {
            cc.push_back(cc[j]*co);
        }
    }
    
    return cc;
}

void draw_helper::dsc_draw_edges_colors(dsc_class & dsc)
{
    static std::vector<vec3> colors = get_random_color(10);
    glPointSize(2.5);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    int max_color = 0;
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        auto pts = dsc.get_pos(dsc.get_nodes(eit.key()));
        if(//!eit->is_interface() &&
           (pts[0][2] < DISPLAY_LIM || pts[1][2] < DISPLAY_LIM))
           continue;
        
        if (dsc.get_color_edge(eit.key()) >= 0)
        {
            glColor3dv(colors[dsc.get_color_edge(eit.key())].get());
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
            
            if (max_color < dsc.get_color_edge(eit.key()))
            {
                max_color = dsc.get_color_edge(eit.key());
            }
        }
    }
    
//    std::cout << max_color << std::endl;
    
    glEnd();
    glEnable(GL_LIGHTING);
}

void draw_helper::dsc_draw_node_color(dsc_class & dsc)
{
    static std::vector<vec3> colors = get_random_color();
    glPointSize(2.5);
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); ++nit)
    {
        
        if(!nit->is_interface() && dsc.get_pos(nit.key())[2] < DISPLAY_LIM)
            continue;
        
        if (dsc.get_color_node(nit.key()) >= 0)
        {
            glColor3dv(colors[dsc.get_color_node(nit.key())].get());
            glVertex3dv(dsc.get_pos(nit.key()).get());
        }
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

void draw_helper::dsc_draw_face_norm(dsc_class & dsc)
{
    glColor3f(0, 0, 1);
    for (auto fid = dsc.faces_begin(); fid != dsc.faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            auto pts = dsc.get_pos(dsc.get_nodes(fid.key()));
            auto center = (pts[0] + pts[1] + pts[2]) / 3.0;

            vec3 Norm = dsc.get_normal(fid.key());

            
            glBegin(GL_LINES);
            glVertex3dv(center.get());
            glVertex3dv((center + Norm*5).get());
            glEnd();
        }
    }
}

void draw_helper::dsc_draw_interface_edge(dsc_class & dsc)
{

    glColor3f(0, 0, 1);
    glBegin(GL_LINES);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        if (eit->is_interface())
        {
            auto pos = dsc.get_pos(dsc.get_nodes(eit.key()));
            glVertex3dv(pos[0].get());
            glVertex3dv(pos[1].get());
        }
    }
    
    glEnd();
}

void draw_helper::save_painting(int WIDTH, int HEIGHT, std::string folder)
{
    std::ostringstream s;
    s << folder << "/scr";
    int i = 0;
    while (1)
    {
        std::ostringstream name;
        name << s.str() << "_" << i << ".png";
        std::ifstream file(name.str().c_str());
        if (!file)
        {
            // could not open
            s << "_" << i << ".png";
            break;
        }
        file.close();
        i++;
    }
    
    int success = SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, WIDTH, HEIGHT);
    if(!success)
    {
        std::cout << "ERROR: Failed to take screen shot: " << s.str().c_str() << std::endl;
        return;
    }
    else{
        std::cout <<"Screen shot: " << s.str().c_str() << std::endl;
    }
}

void draw_helper::dsc_draw_interface(dsc_class & dsc, std::vector<double> * color)
{
    if(color)
    {
        vec3 red(1,0,0);
        vec3 blue(0,0,1);
        
        glBegin(GL_TRIANGLES);
        for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
        {
            if (f->is_interface())
            {
                auto nodes = dsc.get_nodes(f.key());
                auto pts = dsc.get_pos(dsc.get_nodes(f.key()));
                
                

                for (int i = 0; i < 3; i++)
                {
                    auto v = pts[i];
                    double curvature = color->at((int)nodes[i]);
                    // normalize
                    curvature = curvature/30;
//                    curvature = curvature>0? curvature*2-1 : curvature*2 + 1;
                    auto c = blue*(curvature + 1)/2 + red*(1 - curvature)/2;
                    
                    auto norm = -dsc.get_normal(nodes[i]);
                    glColor3d(c[0], c[1], c[2]);
                    glNormal3dv(norm.get());
                    glVertex3dv(v.get());
                }

            }
        }
        glEnd();
    }
    else
    {
        for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
        {
            if (f->is_interface())
            {
                auto pts = dsc.get_pos(dsc.get_nodes(f.key()));
                //auto norm = Util::normal_direction(pts[0], pts[1], pts[2]);
                auto norm = -dsc.get_normal(f.key());
                
//                glColor3f(0.0, 0.9, 1.0);
                glBegin(GL_TRIANGLES);
                for (auto v : pts)
                {
                    glNormal3dv(norm.get());
                    glVertex3dv(v.get());
                }
                glEnd();
                
//
//                glDisable(GL_LIGHTING);
//                glColor3f(0, 0, 1);
//                glBegin(GL_LINES);
//
//                auto edges = dsc.get_edges(f.key());
//
//                for (int i = 0; i < 3; i++)
//                {
//                    glVertex3dv(pts[i].get());
//                    glVertex3dv(pts[(i+1)%3].get());
//                }
//                glEnd();
//                glEnable(GL_LIGHTING);
                
            }
        }
    }
}

void draw_helper::dsc_draw_domain(dsc_class & dsc)
{
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glBegin(GL_TRIANGLES);
    for (auto fit = dsc.faces_begin(); fit != dsc.faces_end(); fit++)
    {
        if (fit->is_boundary())
        {
            auto verts = dsc.get_pos(dsc.get_sorted_nodes(fit.key()));
            vec3 normal = Util::normal_direction(verts[0], verts[1], verts[2]);
            for (auto v : verts)
            {
                glNormal3dv(normal.get());
                glVertex3dv(v.get());
            }
            
        }
    }
    
    
    glEnd();
    glDisable(GL_CULL_FACE);
}
