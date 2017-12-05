//
//  particle.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 04/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef particle_h
#define particle_h

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

class particle
{
public:
    particle(){}
    ~particle(){}
    
    double pressure;
    double density;
    double mass;
    int type;
    int flag;
    int object;
    vec3 vel, pos;
    
    void draw()
    {
        static std::vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
        
        glBegin(GL_POINTS);
        glColor3dv(_color[type].get());
        glVertex3dv(pos.get());
        glEnd();
    }
};

#endif /* particle_h */
