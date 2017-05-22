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
#include "tetralizer.h"
#include "glut_menu.h"


#include <math.h>       /* for cos(), sin(), and sqrt() */

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <mutex>

#include "profile.h"

using namespace DSC;
using namespace std;

int m_iters = 0;

void display_(){
    UI::get_instance()->display();
}

void keyboard_(unsigned char key, int x, int y){
    UI::get_instance()->keyboard(key, x, y);
}

void keyboard_special_(int key, int x, int y){
    UI::get_instance()->keyboard((unsigned char)key, x, y);
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
        angle2= angle2 - (x - startx)*0.03;
        angle = angle + (y - starty)*0.03;

        startx = x;
        starty = y;
        glutPostRedisplay();
    }
}

UI* UI::instance = NULL;

void UI::setup_light()
{
    vec3 center = _obj_dim/ 2.0;
    vec3 eye = center + vec3(gl_dis_max*2.0*cos(angle)*cos(angle2),
                             gl_dis_max*2.0*cos(angle)*sin(angle2),
                             gl_dis_max*2.0*sin(angle));
    
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { -(GLfloat)eye[0], -(GLfloat)eye[1], -(GLfloat)eye[2], 0.0 };
//    glClearColor (0.0, 0.0, 0.0, 0.0);
    glShadeModel (GL_SMOOTH);
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_COLOR_MATERIAL);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
}

void UI::update_draw_list()
{
}

UI::UI(int &argc, char** argv)
{
    instance = this;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
    
    glutCreateWindow("Shadowy Leapin' Lizards");
    
    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
    glutSpecialFunc(keyboard_special_);
    glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    glutMotionFunc(motion_);
    glutMouseFunc(mouse_);


    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.0);
    
    setup_light();
    
	glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();
    
    // Load cross sections
    _seg.init();
    _obj_dim = _seg._img.dimension_v();
    _dsc_dim = _obj_dim + vec3(2*m_edge_length);
    
    cout << "Image dimension" << _obj_dim[0] << " " << _obj_dim[1] << " " <<  _obj_dim[2] << " " ;

    gl_dis_max = fmax(_obj_dim[0], fmax(_obj_dim[1], _obj_dim[2]));

    // Update texture draw
    cout << "Update texture \n";
    draw_helper::update_texture(_seg._img, 0,0,0);

    // Generate DSC
    cout << "Init DSC\n";
    init_dsc();
    cout << "Boundary layer\n";
    set_dsc_boundary_layer();

    std::cout << "Mesh initialized: " << dsc->get_no_nodes() << " nodes; "
                << dsc->get_no_tets() << " tets" << endl;

    _seg._dsc = &*dsc;
//    _seg.initialze_segmentation();
//    _seg.random_initialization();
    _seg.initialization_discrete_opt();
    
}

// Label the gap between DSC boundary and image boundary to BOUND_LABEL (999)
void UI::set_dsc_boundary_layer()
{
    for (auto nit =dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
    {
        if(nit->is_boundary())
        {
            auto tets = dsc->get_tets(nit.key());
            for (auto t : tets)
            {
                dsc->set_label(t, BOUND_LABEL);
            }
        }
    }
}

#define index_cube(x,y,z) ((z)*NX*NY + (y)*NX + (x))
void UI::init_dsc()
{
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;
    

    double delta = m_edge_length;
    
    cout << "delta " << delta << endl;
    cout << _dsc_dim[0] << " " << _dsc_dim[1] << " " <<  _dsc_dim[2] << " " ;
    
    int NX = round(_dsc_dim[0] / delta) + 1; // number of vertices
    int NY = round(_dsc_dim[1] / delta) + 1;
    int NZ = round(_dsc_dim[2] / delta) + 1;
    
        cout << "Compute point" << NX << " " << NY << " " << NZ << "\n";

    // points
    for (int iz = 0; iz < NZ; iz++)
    {
        for (int iy = 0; iy < NY; iy++)
        {
            for (int ix = 0; ix < NX; ix++)
            {
                points.push_back(vec3(ix, iy, iz)*delta - vec3(m_edge_length));
            }
        }
    }

    cout << "Compute tets\n";
    
    // tets
    for (int iz = 0; iz < NZ - 1; iz++)
    {
        for (int iy = 0; iy < NY - 1; iy++)
        {
            for (int ix = 0; ix < NX - 1; ix++)
            {
                // 8 vertices
                int vertices[] = {
                    index_cube(ix, iy, iz),
                    index_cube(ix+1, iy, iz),
                    index_cube(ix+1, iy+1, iz),
                    index_cube(ix, iy+1, iz),
                    index_cube(ix, iy, iz + 1),
                    index_cube(ix+1, iy, iz + 1),
                    index_cube(ix+1, iy+1, iz + 1),
                    index_cube(ix, iy+1, iz + 1)
                };
                
                int tetras[] = {
                    0, 4, 5, 7,
                    0, 7, 5, 1,
                    0, 1, 3, 7,
                    1, 5, 6, 7,
                    1, 6, 7, 3,
                    1, 2, 6, 3
                };
                
                for(int i = 0; i < 6*4; i++)
                {
                    tets.push_back(vertices[tetras[i]]);
                }
            }
        }
    }
    
    long nbTet = tets.size()/4;
    tet_labels = std::vector<int>(nbTet, 0);
    
    cout << "Init DSC from point\n";
    
    dsc = std::unique_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    dsc->set_avg_edge_length(delta);
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
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( /* field of view in degree */ 40.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 10.0, /* Z far */ gl_dis_max*5.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    vec3 center = _obj_dim/ 2.0;
    vec3 eye = center + vec3(gl_dis_max*3.0*cos(angle)*cos(angle2),
                             gl_dis_max*3.0*cos(angle)*sin(angle2),
                             gl_dis_max*3.0*sin(angle));
    vec3 head = vec3(-sin(angle)*cos(angle2),
                     -sin(angle)*sin(angle2),
                     cos(angle));
    gluLookAt(eye[0], eye[1], eye[2], /* eye is at (0,8,60) */
              center[0], center[1], center[2],      /* center is at (0,8,0) */
              head[0], head[1], head[2]);      /* up is in postivie Y direction */
    
    int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
    glViewport((WIN_SIZE_X-size)/2.0, (WIN_SIZE_Y-size)/2.0, size, size);
    
    glClearColor(1.0, 1.0, 1.0, 1.0);
}

void UI::display()
{
    if (CONTINUOUS)
    {
        static int iter = 0;
        if(iter++ > 300)
        {
            CONTINUOUS = false;
            
        }
        _seg.segment();
        m_iters++;
    }
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    update_gl();
    setup_light();
    
    static double total_time = 100;
    static auto init_time = std::chrono::system_clock::now();
    std::chrono::duration<real> t = std::chrono::system_clock::now() - init_time;
    total_time += t.count();
    init_time = std::chrono::system_clock::now();

    //
    if (glut_menu::get_state("Ray line", 0))
    {
        glDisable(GL_LIGHTING);
        glPointSize(2.0);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
        for (auto r : _seg._d_rayz)
        {
            for (int i = 0; i < r.intersects.size()/2; i++)
            {
                for (int j = r.intersects[2*i].z; j < r.intersects[2*i + 1].z; j++)
                {
                    glVertex3f(r.x, r.y, j);
                }
            }

        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    if (glut_menu::get_state("Ray cross section", 0))
    {
        glDisable(GL_LIGHTING);
        glPointSize(2.0);
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
        auto zz = draw_helper::get_instance()._cur_cross_poss[2];
        for (auto r : _seg._d_rayz)
        {
            for (int i = 0; i < r.intersects.size()/2; i++)
            {
                if (zz > r.intersects[2*i].z and zz < r.intersects[2*i + 1].z)
                {
                    glVertex3f(r.x, r.y, zz);
                }
            }

        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    if (glut_menu::get_state("Draw DSC single interface", 1))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
        draw_helper::dsc_draw_one_interface(*dsc, phase_draw);
    }

    if (glut_menu::get_state("Draw tripple edge", 0))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
        glColor3f(1, 0, 0);
        draw_helper::dsc_draw_triple_edge(*dsc);
    }


    if (glut_menu::get_state("Draw DSC edges", 0))
    {
        glColor3f(0, 0, 0);
        draw_helper::dsc_draw_edge(*dsc);
    }

    if (glut_menu::get_state("Draw DSC domain", 1))
    {
        glColor3f(0.3, 0.3, 0.3);
        draw_helper::dsc_draw_domain(*dsc);
    }

    if (glut_menu::get_state("Draw DSC interface edge", 0))
    {
        glDisable(GL_LIGHTING);
        glColor3f(0.0, 0.0, 1.0);
        draw_helper::dsc_draw_interface_edge(*dsc);
    }

    if (glut_menu::get_state("Draw DSC interface", 0))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_LIGHTING);
        draw_helper::dsc_draw_interface(*dsc);
    }


    if (glut_menu::get_state("Draw DSC face normal", 0))
    {
        draw_helper::dsc_draw_face_norm(*dsc);
    }


    if (glut_menu::get_state("Draw Image slide", 1))
    {
        draw_helper::draw_image_slice(_seg._img);
    }
    
    glutSwapBuffers();
    
    std::ostringstream os;
    os << m_iters;
    glutSetWindowTitle(os.str().c_str());

    check_gl_error();
}

void UI::reshape(int width, int height)
{
    WIN_SIZE_X = width;
    WIN_SIZE_Y = height;
    
    update_gl();
    
    glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
}

void UI::animate()
{
}

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case GLUT_KEY_UP:
            draw_helper::update_texture(_seg._img, 0,0,1);
            break;
        case GLUT_KEY_DOWN:
            draw_helper::update_texture(_seg._img, 0,0,-1);
            break;
        case ' ':
            //_seg.segment();
            CONTINUOUS = !CONTINUOUS;
            break;
        case '\t':
            draw_helper::save_painting(WIN_SIZE_X, WIN_SIZE_Y);
            break;
        case 'p':
            profile::close();
            break;
        case 'v':// Change surface type
            phase_draw = (phase_draw+1) % 3;
            break;
        default:
            break;
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
    update_title();
    glutPostRedisplay();
}

void UI::start(const std::string& log_folder_name)
{
    glutPostRedisplay();
}
