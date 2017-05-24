//
//  draw_helper.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "draw_helper.h"
#ifdef __APPLE__
#include <GLUT/GLUT.h>
#endif
#include "DSC.h"
#include <iostream>
#include <string>
#include <SOIL/SOIL.h>
#include <fstream>
#include "define.h"

void draw_helper::draw_image_slice(const image3d & im)
{
    glDisable(GL_LIGHTING);
    // Draw bounding box
    vec3 dim = im.dimension_v();
    
    glPushMatrix();
    
    glTranslatef(dim[0]/2.0, dim[1]/2.0, dim[2]/2.0);
    glScalef(dim[0], dim[1], dim[2]);
    
    glColor3f(1, 0, 0);
    glutWireCube(1.0);
    
    glPopMatrix();
    
    // slice
    vec3 ld(0,0, get_instance()._cur_cross_poss[2]);
    vec3 ru(dim[0], dim[1], get_instance()._cur_cross_poss[2]);
    
    glEnable(GL_TEXTURE_2D);
    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);
    
    glTexCoord3f(0., 0., 0.);
    glVertex3f(ld[0], ld[1], ld[2]);

    glTexCoord3f(1., 0., 0.);
    glVertex3f(ru[0], ld[1], ld[2]);
    
    glTexCoord3f(1., 1., 0.);
    glVertex3f(ru[0], ru[1], ld[2]);
    
    glTexCoord3f(0., 1., 0.);
    glVertex3f(ld[0], ru[1], ld[2]);
    
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    // border
    glColor3f(1, 0, 0);
    glBegin(GL_LINE_LOOP);
    glVertex3f(ld[0], ld[1], ld[2]);
    glVertex3f(ru[0], ld[1], ld[2]);
    glVertex3f(ru[0], ru[1], ld[2]);
    glVertex3f(ld[0], ru[1], ld[2]);
    glEnd();
}

void draw_helper::draw_coord(float length)
{
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex3f(length, 0., 0.);
    glVertex3f(0., 0., 0.);
    
    glColor3f(0, 1, 0);
    glVertex3f(0., length, 0.);
    glVertex3f(0., 0., 0.);
    
    glColor3f(0, 0, 1.);
    glVertex3f(0., 0., length);
    glVertex3f(0., 0., 0.);
    glEnd();
}

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

enum interface_type
{
    type_0_2 = 0,
    type_0_1,
    type_1_2
};

std::vector<vec3> draw_helper::node_normal_vector;
void draw_helper::update_normal_vector_interface(dsc_class & dsc, int phase, vec3 eye_pos)
{
#ifdef _DSC_ORIGIN_
    std::vector<int> neighbor_faces_count(MAX_NUM_ELEMENT_MESH, 0);
#else
    node_normal_vector = std::vector<vec3>(dsc.get_no_nodes_buffer(), vec3(0.0));
    std::vector<int> neighbor_faces_count(dsc.get_no_nodes_buffer(), 0);
#endif
    
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase))
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes = *dsc.get_nodes_cache(f.key());
            auto nodes_pos = dsc.get_pos(nodes);
#else
            auto nodes = dsc.get_nodes(f.key());
            auto nodes_pos = dsc.get_pos(nodes);
#endif
            
            auto norm = Util::normal_direction(nodes_pos[0], nodes_pos[1], nodes_pos[2]);
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[1] : tets[0];
            auto other_node = dsc.get_nodes(other_tet) - nodes;
            vec3 direct = dsc.get_pos(other_node[0]) - nodes_pos[0];
//            auto direct = eye_pos - nodes_pos[0];direct.normalize();
            norm = norm*Util::dot(norm, direct);
            
            for(auto n : nodes)
            {
                node_normal_vector[n] += norm;
                neighbor_faces_count[n]++;
            }
        }
    }
    
    for(int i = 0; i < node_normal_vector.size(); i++)
    {
        if(neighbor_faces_count[i] > 0)
            node_normal_vector[i] /= (double)neighbor_faces_count[i];
    }
}



void draw_helper::dsc_draw_one_interface(dsc_class & dsc, int phase)
{
    
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {
            auto tets = dsc.get_tets(f.key());
            if (!(dsc.get_label(tets[0]) == phase
                  or dsc.get_label(tets[1]) == phase))
            {
                continue;
            }
#ifdef DSC_CACHE
            auto nodes=*dsc.get_nodes_cache(f.key());
#else
            auto nodes=dsc.get_nodes(f.key());
#endif
            auto pts = dsc.get_pos(nodes);
            auto norm = dsc.get_normal(f.key());
            
            // normalize the normal to the eye
            is_mesh::TetrahedronKey other_tet = (dsc.get_label(tets[0]) == phase)? tets[1] : tets[0];
            auto other_node = dsc.get_nodes(other_tet) - nodes;
            vec3 direct = dsc.get_pos(other_node[0]) - pts[0];
            //            auto direct = eye_pos - nodes_pos[0];direct.normalize();
            norm = norm*Util::dot(norm, direct);
            
//            // Draw edges
//            glDisable(GL_LIGHTING);
//            glColor3f(0, 0, 1);
//            glBegin(GL_LINES);
//            
//            auto edges = dsc.get_edges(f.key());
//            
//            for (int i = 0; i < 3; i++)
//            {
//                glVertex3dv(pts[i].get());
//                glVertex3dv(pts[(i+1)%3].get());
//            }
//            glEnd();
//            glEnable(GL_LIGHTING);
            
            // Draw triangle
            glColor3f(0.7, 0.7, 0.7);
            glBegin(GL_TRIANGLES);
            for (int i =0; i < 3; i++)
            {
                auto v = pts[i];
                auto n = nodes[i];
                vec3 real_norm = norm;
                if((int)n < node_normal_vector.size() && node_normal_vector[n].length() > 0.1)
                    real_norm = node_normal_vector[n];
                
                glNormal3dv(real_norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
        
            
        }
    }
}

void draw_helper::dsc_draw_interface(dsc_class & dsc)
{
    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface() && !f->is_boundary())
        {

            auto pts = dsc.get_pos(dsc.get_nodes(f.key()));
            //auto norm = Util::normal_direction(pts[0], pts[1], pts[2]);
            auto norm = -dsc.get_normal(f.key());
            
            glColor3f(0.7, 0.0, 0);
            glBegin(GL_TRIANGLES);
            for (auto v : pts)
            {
                glNormal3dv(norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
            
//            
//            glDisable(GL_LIGHTING);
//            glColor3f(0, 0, 0);
//            glBegin(GL_LINES);
//            
//            auto edges = dsc.get_edges(f.key());
//            
//            for (int i = 0; i < 3; i++)
//            {
//                
//             //   glNormal3dv(norm.get());
//                glVertex3dv(pts[i].get());
//                
//             //   glNormal3dv(norm.get());
//                glVertex3dv(pts[(i+1)%3].get());
//            }
//            glEnd();
//            glEnable(GL_LIGHTING);
            
        }
    }
}

#define P_NONE  0x0000
#define P_ZERO  0x0001
#define P_ONE   0x0010
#define P_TWO   0x0100

#define P_ALL   0x0111

void draw_helper::dsc_draw_triple_edge(dsc_class & dsc)
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for (auto eid = dsc.edges_begin(); eid != dsc.edges_end(); eid++)
    {
        auto tets = dsc.get_tets(eid.key());
        // separate more than 3 phase
        int p = P_NONE;
        int pl[] = {P_ZERO, P_ONE, P_TWO};
        for (auto t : tets)
        {
            int label = dsc.get_label(t);
            p = p & pl[label];
        }
        
        if (p == P_ALL)
        {
            auto pts = dsc.get_pos(dsc.get_nodes(eid.key()));
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
        }
    }
    glEnd();
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

void draw_helper::update_texture(const image3d & im,
                    int const &off_x, int const & off_y, int const & off_z)
{
    
    
    auto dim = im.dimension_v();
    
    get_instance()._cur_cross_poss += CGLA::Vec3i(off_x, off_y, off_z);
    for (int i = 0; i < 3; i++)
    {
        if (get_instance()._cur_cross_poss[i] < 0)
        {
            get_instance()._cur_cross_poss[i] = 0;
        }
        if (get_instance()._cur_cross_poss[i] > dim[i] - 1)
        {
            get_instance()._cur_cross_poss[i] = dim[i] - 1;
        }
    }

    int z = get_instance()._cur_cross_poss[2];
    int width = dim[0];
    int height = dim[1];
    uint8_t * data = (uint8_t *)malloc(dim[0] * dim[1] * 3 * sizeof(uint8_t));
    
    uint8_t *ptr = data;
    for (int j = 0; j < dim[1]; j++)
    {
        for (int i = 0; i < dim[0]; i++)
        {
            uint8_t v = (uint8_t)( im.get_value_f(i,j,z) * 255 );
            *(ptr++) = v;
            *(ptr++) = v;
            *(ptr++) = v;
        }
    }
    
    
    
    static GLuint tex_ID = 0;
    if (tex_ID != 0)
    {
        glDeleteTextures(1, &tex_ID);
        tex_ID = 0;
    }

    
    glGenTextures(1, &tex_ID);
    
    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    glPixelStorei ( GL_UNPACK_ALIGNMENT,   1 );
    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    delete [] data;
}
