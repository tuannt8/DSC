//
//  draw_helper.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "draw_helper.h"
#include "DSC.h"
#include <GLUT/GLUT.h>

namespace draw_helper
{
    void draw_image_slice(const image3d & im)
    {
        // Draw bounding box
        vec3 dim = im.dimension_v();
        glPushMatrix();
        glTranslatef(dim[0], dim[1], dim[2]);
        glScalef(1./dim[0], 1./dim[1], 1./dim[2]);
        glColor3f(1, 0, 0);
        glutWireCube(1.0);
        
        glPopMatrix();
        
        glutSolidTeapot(100.0);
    }
    
    void draw_coord(float length)
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
}