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
#include "rotate_function.h"
#include "average_function.h"
#include "normal_function.h"
#include "draw_helper.h"


#include <math.h>       /* for cos(), sin(), and sqrt() */

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace DSC;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void reshape_(int width, int height){
    UI::get_instance()->reshape(width, height);
}

void visible_(int v){
    UI::get_instance()->visible(v);
}

void animate_(){
    UI::get_instance()->animate();
}

void motion_(int x, int y)
{
    UI::get_instance()->motion(x, y);
}

void mouse_(int button, int state, int x, int y){
    UI::get_instance()->mouse(button, state, x, y);
}

void UI::mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            moving = 1;
            startx = x;
            starty = y;
        }
        if (state == GLUT_UP) {
            moving = 0;
        }
    }
}

void UI::motion(int x, int y)
{
    if (moving) {
        angle = angle + (x - startx);
        angle2 = angle2 + (y - starty);
        startx = x;
        starty = y;
        glutPostRedisplay();
    }
}

UI* UI::instance = NULL;

UI::UI(int &argc, char** argv)
{
    instance = this;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
    
    glutCreateWindow("Shadowy Leapin' Lizards");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
	glutIgnoreKeyRepeat(true);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    glutMotionFunc(motion_);
    glutMouseFunc(mouse_);


//    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);
    glLineWidth(3.0);
    
    glMatrixMode(GL_PROJECTION);
    gluPerspective( /* field of view in degree */ 40.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 20.0, /* Z far */ 100.0);
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(0.0, 8.0, 60.0,  /* eye is at (0,8,60) */
              0.0, 0.0, 0.0,      /* center is at (0,8,0) */
              0.0, 1.0, 0.);      /* up is in postivie Y direction */
    
    // Lighting
    const float amb = 2.0;
    const float LightAmbient[][4]  = {  { amb, amb, amb, 1.0f },
        { amb, amb, amb, 1.0f }
    };
    const float LightDiffuse[] [4] = {  { 1.0f, 1.0f, 1.0f, 1.0f },
        { 1.0f, 1.0f, 1.0f, 1.0f }
    };
    const float LightPosition[][4] = {  { 1.0f,  4.0f, 2.0f, 0.0f },
        { 0.0f, 10.0f, 0.0f, 1.0f }
    };
    
    glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient[0]);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse[0]);
    glLightfv(GL_LIGHT0, GL_POSITION, LightPosition[0]);
    
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient[1]);
    // etc., snip -- no LIGHT1 for this round
    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
    
    // Read input
    std::string motion = "";
    real discretization = 2.5;
    real velocity = 5.;
    real accuracy = 0.25;
    
    if(argc == 2)
    {
        model_file_name = std::string(argv[1]);
    }
    else if(argc > 2)
    {
        for(int i = 0; i < argc; ++i)
        {
            std::string str(argv[i]);
            if (str == "nu") {
                velocity = std::atof(argv[i+1]);
            }
            else if (str == "delta") {
                discretization = std::atof(argv[i+1]);
            }
            else if (str == "alpha") {
                accuracy = std::atof(argv[i+1]);
            }
            else if (str == "model") {
                model_file_name = argv[i+1];
            }
            else if (str == "motion") {
                motion = argv[i+1];
            }
        }
    }
//    painter = std::unique_ptr<Painter>(new Painter(light_pos));
//    load_model(model_file_name, discretization);
    
    
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
    
    _seg.init();
    auto dims = _seg.get_image().dimension_v();
    gl_dis_max = fmax(dims[0], fmax(dims[1], dims[2]));
    // Update view
    
}

void UI::load_model(const std::string& file_name, real discretization)
{
    std::cout << "\nLoading " << obj_path + file_name + ".dsc" << std::endl;
    dsc = nullptr;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(obj_path + file_name + ".dsc", points, tets, tet_labels);
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    
    vec3 p_min(INFINITY), p_max(-INFINITY);
    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++) {
        for (int i = 0; i < 3; i++) {
            p_min[i] = Util::min(nit->get_pos()[i], p_min[i]);
            p_max[i] = Util::max(nit->get_pos()[i], p_max[i]);
        }
    }
    
    vec3 size = p_max - p_min;
    real var = Util::max(Util::max(size[0], size[1]), size[2]);
    real dist = 1.2*var;
    eye_pos = {dist, var, dist};
    camera_pos = {var, var, -dist};
    light_pos = {0., 0., dist};
    
//    painter->update(*dsc);
    std::cout << "Loading done" << std::endl << std::endl;
}

void UI::update_title()
{
    std::ostringstream oss;
    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
    std::string str(oss.str());
    glutSetWindowTitle(str.c_str());
}

void UI::update_gl()
{
//    glLoadIdentity();
//
    
//    glMatrixMode(GL_PROJECTION);
//    gluPerspective( /* field of view in degree */ 40.0,
//                   /* aspect ratio */ 1.0,
//                   /* Z near */ 0.0, /* Z far */ gl_dis_max*10);
//    glMatrixMode(GL_MODELVIEW);
//    gluLookAt(0.0, 0.0, gl_dis_max,  /* eye is at (0,8,60) */
//              0.0, 0.0, 0.0,      /* center is at (0,8,0) */
//              0.0, 1.0, 0.);      /* up is in postivie Y direction */
    
    glViewport(0, 0, WIN_SIZE_X, WIN_SIZE_Y);
    
//    // Lighting
//    const float amb = 2.0;
//    const float LightAmbient[][4]  = {  { amb, amb, amb, 1.0f },
//        { amb, amb, amb, 1.0f }
//    };
//    const float LightDiffuse[] [4] = {  { 1.0f, 1.0f, 1.0f, 1.0f },
//        { 1.0f, 1.0f, 1.0f, 1.0f }
//    };
//    const float LightPosition[][4] = {  { 1.0f,  4.0f, 2.0f, 0.0f },
//        { 0.0f, 10.0f, 0.0f, 1.0f }
//    };
//    
//    glLightfv(GL_LIGHT0, GL_AMBIENT, LightAmbient[0]);
//    glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse[0]);
//    glLightfv(GL_LIGHT0, GL_POSITION, LightPosition[0]);
//    
//    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient[1]);
//    // etc., snip -- no LIGHT1 for this round
//    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
//    glEnable(GL_LIGHT0);
//    //glEnable(GL_LIGHT1);
//    glEnable(GL_LIGHTING);
}

void UI::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
 //   update_gl();
    
    glPushMatrix();
    /* Perform scene rotations based on user mouse input. */
    glRotatef(angle2, 1.0, 0.0, 0.0);
    glRotatef(angle, 0.0, 1.0, 0.0);
    
    glutSolidTeapot(gl_dis_max/10.0);
//    draw_helper::draw_image_slice(_seg.get_image());
    draw_helper::draw_coord(gl_dis_max);
    
    glPopMatrix();
    
    glutSwapBuffers();

    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    
    update_gl();
    
//    painter->reshape(width, height);
    glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::animate()
{
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {

    }
}

void idle(void)
{
    glutPostRedisplay();
}

void UI::visible(int v)
{
    if (v == GLUT_VISIBLE) {
        if (animation)
            glutIdleFunc(idle);
    } else {
        if (!animation)
            glutIdleFunc(NULL);
    }
}

void UI::stop()
{
    if(RECORD && basic_log)
    {
        basic_log->write_message("MOTION STOPPED");
        basic_log->write_log(*dsc);
        basic_log->write_log(*vel_fun);
        basic_log->write_timings(*vel_fun);
        
        std::vector<vec3> points;
        std::vector<int> faces;
        std::vector<int> tets;
        std::vector<int> tet_labels;
        dsc->extract_tet_mesh(points, tets, tet_labels);
        is_mesh::export_tet_mesh(basic_log->get_path() + std::string("/mesh.dsc"), points, tets, tet_labels);
        points.clear();
        dsc->extract_surface_mesh(points, faces);
        is_mesh::export_surface_mesh(basic_log->get_path() + std::string("/mesh.obj"), points, faces);
        basic_log = nullptr;
    }
    
    CONTINUOUS = false;
//    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    if(RECORD)
    {
//        painter->set_view_position(camera_pos);

    }
    
//    painter->update(*dsc);
//    update_title();
    glutPostRedisplay();
}
