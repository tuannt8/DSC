//
//  draw_helper.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "draw_helper.h"
#include <GLUT/GLUT.h>

void draw_helper::draw_image_slice(const image3d & im)
{
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
    glBegin(GL_LINES);
    for (auto eit = dsc.edges_begin(); eit != dsc.edges_end(); eit++)
    {
        auto pos = dsc.get_pos(dsc.get_nodes(eit.key()));
        glVertex3dv(pos[0].get());
        glVertex3dv(pos[1].get());
    }
    
    glEnd();
}

void draw_helper::dsc_draw_interface(dsc_class & dsc)
{

    for (auto f = dsc.faces_begin(); f != dsc.faces_end(); f++)
    {
        if (f->is_interface())
        {
            glColor3f(1.0, 0.0, 0.0);
            glBegin(GL_TRIANGLES);
            auto pts = dsc.get_pos(dsc.get_nodes(f.key()));
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
            glVertex3dv(pts[2].get());
            glEnd();
            
            glColor3f(0, 0, 1);
            glBegin(GL_LINES);
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[1].get());
            
            glVertex3dv(pts[0].get());
            glVertex3dv(pts[2].get());
            
            glVertex3dv(pts[1].get());
            glVertex3dv(pts[2].get());
            glEnd();
        
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
    
    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    
    glBindTexture(GL_TEXTURE_2D, tex_ID);
    
    delete [] data;
}