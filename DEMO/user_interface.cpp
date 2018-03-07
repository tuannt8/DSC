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

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

#include "profile.h"

#include "draw_helper.h"
#include "glut_menu.hpp"

#include "debugger.h"

#include <sys/types.h>
#include <dirent.h>

extern double g_res;

using namespace DSC;
using namespace std;

bool mouse_press = 0;
int _dx = 0; int _dy = 0;
int _x = 0; int _y=0;

int mode = 0;

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

void mouse_move_(int x, int y)
{
    if (mouse_press)
    {
        UI::get_instance()->angle += _dy/30.;
        UI::get_instance()->angle2 += _dx/30.;
        glutPostRedisplay();
    }
    _dx = x - _x;
    _dy = y - _y;
    _x = x;
    _y = y;
}

void mouse_down_(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON)
    {
        _x = x;
        _y = y;
        
        if (state == GLUT_DOWN)
        {
            mouse_press = 1;
        }
        if (state == GLUT_UP)
        {
            mouse_press = 0;
        }
    }
}

UI* UI::instance = NULL;

void UI::setup_light()
{
    vec3 center = _obj_dim / 2.0;
    vec3 eye = center + vec3(gl_dis_max*2.0*cos(angle)*cos(angle2),
                             gl_dis_max*2.0*sin(angle),
                             gl_dis_max*2.0*cos(angle)*sin(angle2)
                             );
    
    GLfloat light_position[] = { (GLfloat)eye[0], (GLfloat)eye[1], (GLfloat)eye[2], 0.0 };
    glShadeModel (GL_SMOOTH);
    
    GLfloat amb = 0.2, diff = 1., spec = 1.;

    GLfloat light_ambient[] = { amb,amb,amb, 1.0 };
    GLfloat light_diffuse[] = {diff, diff, diff, 1.0 };
    GLfloat light_specular[] = {spec, spec, spec, 1.0 };

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
//
//    GLfloat g_amb = 1.0;
//    GLfloat global_ambient[] = {g_amb, g_amb, g_amb, 0.1};
//    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
//
//
//    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//    GLfloat mat_shininess[] = { 5.0 };
//
//    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
//
//    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
//    glEnable(GL_COLOR_MATERIAL);
}
inline vec3 min_vec(vec3 v1, vec3 v2)
{
    vec3 out;
    out[0] = std::min(v1[0], v2[0]);
    out[1] = std::min(v1[1], v2[1]);
    out[2] = std::min(v1[2], v2[2]);
    return out;
}

inline vec3 max_vec(vec3 v1, vec3 v2)
{
    vec3 out;
    out[0] = std::max(v1[0], v2[0]);
    out[1] = std::max(v1[1], v2[1]);
    out[2] = std::max(v1[2], v2[2]);
    return out;
}

void UI::update_gl()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( /* field of view in degree */ 40.0,
                   /* aspect ratio */ 1.0,
                   /* Z near */ 0.1, /* Z far */ gl_dis_max*10.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    vec3 center = _obj_dim / 2.0;
    double dis = 1.5;
    
    vec3 eye = center + vec3(gl_dis_max*dis*cos(angle)*cos(angle2),
                             gl_dis_max*dis*sin(angle),
                             gl_dis_max*dis*cos(angle)*sin(angle2)
                             );
    vec3 head = vec3(-sin(angle)*cos(angle2),
                     cos(angle),
                     -sin(angle)*sin(angle2)
                     );
    
    gluLookAt(eye[0], eye[1], eye[2], /* eye is at (0,8,60) */
              center[0], center[1], center[2],      /* center is at (0,8,0) */
              head[0], head[1], head[2]);      /* up is in postivie Y direction */
    
    int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
    glViewport((WIN_SIZE_X - size) / 2.0, (WIN_SIZE_Y - size) / 2.0, size, size);
    
    glClearColor(0.1, 0.1, 0.1, 1.0);
}

UI::UI()
{
    init_data();
}

void UI::init_data()
{
#ifdef __APPLE__
    bool test = true;
#else
    bool test = false;
#endif
    
    
    m_fluid.load_configuration();// Load configuration
    
    // init DSC
    _obj_dim = m_fluid.m_problem->domain_size();
    gl_dis_max = std::max(std::max(_obj_dim[0], _obj_dim[1]), _obj_dim[2])*1.7;
    
    if(!test)
        dsc = std::shared_ptr<DeformableSimplicialComplex<>>(m_fluid.m_problem->init_dsc(g_res));
    else
    {
        dsc = load_model("/Users/tuannt8/Desktop/iter.dsc");
//        dsc->validity_check();
    }
//    {
//        std::vector<vec3> points;
//        std::vector<int> faces;
//        dsc->extract_surface_mesh(points, faces);
//        is_mesh::export_surface_mesh("/Users/tuannt8/Desktop/iter.obj", points, faces);
//    }
//    dsc->validity_check();
//    return;
    
    m_fluid.init(&*dsc);
    
    // Load first particle
    m_fluid.load_first_particle();
    
    if(!test)
    {
//        m_fluid.init_mesh();
    }
    
    dsc->print_mesh_info();
}

void UI::export_surface(const std::string& dsc_path)
{
    dsc = load_model(dsc_path);
    std::vector<vec3> points;
    std::vector<int> faces;
    dsc->extract_surface_mesh(points, faces);
    
    // snapp
//    double epsilon = m_fluid.m_problem->m_deltap;
//    vec3 origin(0.0);
//    vec3 bound = m_fluid.m_problem->domain_size();
//    for(auto & p : points)
//    {
//        for (int i = 0; i < 3; i++)
//        {
//            if (p[i] - origin[i] < epsilon)
//            {
//                p[i] = origin[i];
//            }
//            if (bound[i] - p[i] < epsilon)
//            {
//                p[i] = bound[i];
//            }
//        }
//    }
    
    string directory = dsc_path.substr(0, dsc_path.find_last_of("\\/"));
    string name = dsc_path.substr(dsc_path.find_last_of("\\/") + 1);
    string phase =  directory + "/"  + "sur_" + name + ".obj";
    is_mesh::export_surface_mesh(phase, points, faces);
    
    cout << "Export to " << phase << endl;
}

void UI::export_dam_break()
{
    string dir_path = "../Large_data/dam_break_fluid/v_2_dam_0.1_0";
    
    DIR *dp;
    struct dirent *ep;
    dp = opendir ("../Large_data/dam_break_fluid/v_2_dam_0.1_0");
    
    if (dp != NULL)
    {
        while (ep = readdir (dp))
        {
            string name(ep->d_name);
            if(name.length() > 5)
            {
                if (strcmp(name.substr(name.size() - 4).c_str(), ".dsc") == 0)
                {
                    export_surface(dir_path + "/" + name);
                }
            }
        }
        
        (void) closedir (dp);
    }
    else
        perror ("Couldn't open the directory");
}

UI::UI(int &argc, char** argv)
{
    instance = this;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
    
    glutCreateWindow("Shadowy Leapin' Lizards");

    glutDisplayFunc(display_);
    glutKeyboardFunc(keyboard_);
    glutSetKeyRepeat(GLUT_KEY_REPEAT_ON);
    glutVisibilityFunc(visible_);
    glutReshapeFunc(reshape_);
    glutMotionFunc(mouse_move_);
    glutMouseFunc(mouse_down_);
    
 
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.0);
   
    setup_light();
    
    glutReshapeWindow(WIN_SIZE_X, WIN_SIZE_Y);
    check_gl_error();

    init_data();

}

void UI::init_dam_break()
{

}



void UI::init_dsc()
{
    std::vector<vec3> points;
    std::vector<int> tets;
    std::vector<int> tet_labels;
    
#ifdef __APPLE__
    int DISCRETIZATION = 50;
#else
    int DISCRETIZATION = 200;
#endif
    double delta = _obj_dim[0]/(double)DISCRETIZATION;
    
    vec3 _dsc_dim = _obj_dim + vec3(delta)*2;
    
    cout << "delta " << delta << endl;
    cout << _dsc_dim[0] << " " << _dsc_dim[1] << " " <<  _dsc_dim[2] << " " ;
    
    int NX = round(_dsc_dim[0] / delta) + 1; // number of vertices
    int NY = round(_dsc_dim[1] / delta) + 1;
    int NZ = round(_dsc_dim[2] / delta) + 1;
    
    cout << "Compute point" << NX << " " << NY << " " << NZ << "\n";
    
    double deltax = _dsc_dim[0]/(double)(NX-1);
    double deltay = _dsc_dim[1]/(double)(NY - 1);
    double deltaz = _dsc_dim[2]/(double)(NZ - 1);
    
    
    // points. Push it back
    for (int iz = 0; iz < NZ; iz++)
    {
        for (int iy = 0; iy < NY; iy++)
        {
            for (int ix = 0; ix < NX; ix++)
            {
                points.push_back(vec3(ix*deltax, iy*deltay, iz*deltaz) - vec3(delta));
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
    
    dsc = std::shared_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));
    dsc->set_avg_edge_length(delta);
    
    
}

std::shared_ptr<DeformableSimplicialComplex<>> UI::load_model(const std::string& file_name)
{
    std::cout << "\nLoading " <<file_name << std::endl;
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
    
    auto l_dsc = std::shared_ptr<DeformableSimplicialComplex<>>(new DeformableSimplicialComplex<>(points, tets, tet_labels));

    std::cout << "Loading done" << std::endl << std::endl;
    
    return l_dsc;
}

void UI::update_title()
{
//    std::ostringstream oss;
//    oss << "3D DSC\t" << vel_fun->get_name() << ", Time step " << vel_fun->get_time_step();
//    oss << " (Nu = " << vel_fun->get_velocity() << ", Delta = " << dsc->get_avg_edge_length() << ", Alpha = " << vel_fun->get_accuracy() << ")";
//    std::string str(oss.str());
//    glutSetWindowTitle(str.c_str());
}

void UI::draw_dsc_layer(double y_lim)
{
    vec3 look(0,0,1);
    
    for (auto f = dsc->faces_begin(); f != dsc->faces_end(); f++)
    {
        auto pts = dsc->get_pos(dsc->get_nodes(f.key()));
        bool bDraw = true;
        for (auto p : pts)
        {
            bDraw = bDraw & (p[2] > y_lim);
        }
        
        if (bDraw)
        {
            

            //auto norm = Util::normal_direction(pts[0], pts[1], pts[2]);
            auto norm = dsc->get_normal(f.key());
            norm = norm*(Util::dot(norm, look));
            
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

void UI::draw_dsc_layer_1(double y_lim)
{
    vec3 look(0,0,1);
    
    for (auto f = dsc->faces_begin(); f != dsc->faces_end(); f++)
    {
        auto pts = dsc->get_pos(dsc->get_nodes(f.key()));
        bool bDraw = true;
        for (auto p : pts)
        {
            bDraw = bDraw & (p[2] > y_lim);
        }
        
        auto tets = dsc->get_tets(f.key());
        if(tets.size()!=2)
            glColor3f(0.5, 0.5, 0.5);
        else if ( (dsc->get_label(tets[0]) == 1) || (dsc->get_label(tets[1]) == 1) )
            glColor3f(0.0, 0.9, 1.0);
        else
            glColor3f(0.5, 0.5, 0.5);
            
        
        
        if (bDraw)
        {
            
            
            //auto norm = Util::normal_direction(pts[0], pts[1], pts[2]);
            auto norm = dsc->get_normal(f.key());
            norm = norm*(Util::dot(norm, look));
            
            glBegin(GL_TRIANGLES);
            for (auto v : pts)
            {
                glNormal3dv(norm.get());
                glVertex3dv(v.get());
            }
            glEnd();
            
            
            
        }
    }
}

void UI::display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    update_gl();
    setup_light();
    

    if(glut_menu::get_state("Particles", 1))
        m_fluid.draw();
    
    if(dsc)
    {
        if(glut_menu::get_state("DSC domain", 1))
       {
           glEnable(GL_LIGHTING);
            glColor3f(0.8, 0.8, 0.8);
            draw_helper::dsc_draw_domain(*dsc);
       }
        
        if(glut_menu::get_state("DSC interface", 0))
        {
            glEnable(GL_LIGHTING);
            glColor3f(0.0, 0.9, 1.0);
            draw_helper::dsc_draw_interface(*dsc);
        }
        
        if(glut_menu::get_state("inverted tets", 0))
        {
            glEnable(GL_LIGHTING);
            glColor3f(0.0, 0.9, 1.0);
            draw_helper::dsc_draw_inverted_tets(*dsc);
        }
        
        if(glut_menu::get_state("DSC shared interface", 0))
        {
            glEnable(GL_LIGHTING);
            glColor3f(0.0, 0.9, 1.0);
            draw_helper::dsc_draw_shared_interface(*dsc);
        }
        
        if(glut_menu::get_state("DSC interface edges", 1))
        {
            glDisable(GL_LIGHTING);
            glColor3f(0.0, 0.0, 1.0);
            draw_helper::dsc_draw_interface_edge(*dsc);
        }
        
//        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
//        {
//            if (nit->is_interface())
//            {
//                glPushMatrix();
//                auto p = nit->get_pos();
//                glTranslated(p[0], p[1], p[2]);
//                glutSolidSphere(0.0012, 20, 20);
//                glPopMatrix();
//            }
//        }
    }
    

//    if(mode == 0)
//    {
//    glColor3f(0.0, 0.9, 1.0);
//        std::vector<double> * cur = _node_curvature.size()>0? &_node_curvature : nullptr;
//
//    draw_helper::dsc_draw_interface(*dsc, cur);
//
//    glColor3f(0.8, 0.8, 0.8);
//    draw_helper::dsc_draw_domain(*dsc);
//
////    glColor3f(0.5, 0.5, 0.5);
////    draw_dsc_layer(DISPLAY_LIM);
//    }
//    else if(mode == 1)
//    {
//        std::vector<double> * cur = _node_curvature.size()>0? &_node_curvature : nullptr;
//
//        glColor3f(0.0, 0.9, 1.0);
//        draw_helper::dsc_draw_interface(*dsc, cur);
//
//        glColor3f(0.8, 0.8, 0.8);
//        draw_helper::dsc_draw_domain(*dsc);
//    }
//    else if(mode == 2)
//    {
//        glColor3f(0.8, 0.8, 0.8);
//        draw_helper::dsc_draw_domain(*dsc);
//
//        glColor3f(0.5, 0.5, 0.5);
//        draw_dsc_layer_1(DISPLAY_LIM);
//    }
    
    glutSwapBuffers();
//    update_title();
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
    if(CONTINUOUS)
    {
        m_fluid.deform();
    }
    glutPostRedisplay();
}


#ifdef DSC_CACHE

#endif

is_mesh::SimplexSet<is_mesh::EdgeKey> * get_edge_around(DSC::DeformableSimplicialComplex<> *dsc, is_mesh::NodeKey n)
{
    is_mesh::SimplexSet<is_mesh::EdgeKey> * edge_around = new is_mesh::SimplexSet<is_mesh::EdgeKey>;
    
    auto faces = dsc->get_faces(n);
    for (auto f : faces)
    {
        if (dsc->get(f).is_interface())
        {
            auto edges = dsc->get_edges(f);
            for (auto e : edges)
            {
                auto ns = dsc->get_nodes(e);
                if (ns[0] != n && ns[1] != n)
                {
                    edge_around->push_back(e);
                    break;
                }
            }
        }
    }
    
    return edge_around;
}

#ifdef DSC_CACHE

#endif

void UI::keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case 'p':
            profile::close();
            break;
        case 'v':
            mode=(mode+1)%3;
            break;
        case '\033':
            stop();
            exit(0);
            break;
        case '0':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::shared_ptr<VelocityFunc<>>(new VelocityFunc<>(vel_fun->get_velocity(), vel_fun->get_accuracy(), 500));
            start("");
            break;
        case '1':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::shared_ptr<VelocityFunc<>>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("rotate");
            break;
        case '2':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::shared_ptr<VelocityFunc<>>(new AverageFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("smooth");
            break;
        case '3':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::shared_ptr<VelocityFunc<>>(new NormalFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("expand");
            break;
        case ' ':
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            vel_fun->take_time_step(*dsc);
            painter->update(*dsc);
            break;
        case 'r':
//            m_fluid.m_file_load.fix_output_boundary();
            break;
        case 't':
            std::cout << "TEST VELOCITY FUNCTION" << std::endl;
            vel_fun->test(*dsc);
            painter->update(*dsc);
            break;
        case '\t':
            painter->switch_display_type();
            painter->update(*dsc);
            break;
        case 's':
        {
            std::cout << "TAKING SCREEN SHOT" << std::endl;
//            painter->set_view_position(camera_pos);
//            painter->save_painting("LOG");
        }
            break;
        case 'e':
        {
            std::cout << "EXPORTING MESH" << std::endl;
            std::string filename("data/mesh.dsc");
            std::vector<vec3> points;
            std::vector<int> tets;
            std::vector<int> tet_labels;
            dsc->extract_tet_mesh(points, tets, tet_labels);
            is_mesh::export_tet_mesh(filename, points, tets, tet_labels);
        }
            break;
        case 'i':
        {
            std::cout << "EXPORTING SURFACE MESH" << std::endl;
            std::string filename("data/mesh.obj");
            std::vector<vec3> points;
            std::vector<int> faces;
            dsc->extract_surface_mesh(points, faces);
            is_mesh::export_surface_mesh(filename, points, faces);
        }
            break;
        case '+':
        {
            gl_dis_max *= 1.05;
        }
            break;
        case '-':
        {
            gl_dis_max /= 1.05;
        }
            break;
        case '.':
        {
            debugger<>::get_int("particle idx", 0) ++;
        }
            break;
        case ',':
        {
            debugger<>::get_int("particle idx", 0) --;
        }
            break;
        case '<':
        {
            real accuracy = std::min(vel_fun->get_accuracy() + 1., 100.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
        case 'n':
        {
            m_fluid.load_next_particle();
        }
            break;
        case 'l':
        {
            profile::close();
        }
            break;
        case 'b':
        {

            dsc->print_mesh_info();
        }
            break;
    }
}

void UI::visible(int v)
{
    if(v==GLUT_VISIBLE)
        glutIdleFunc(animate_);
    else
        glutIdleFunc(0);
    
}

//#define CACHE_GRADIENT

double UI::mass(const is_mesh::SimplexSet<is_mesh::FaceKey> & face_around, vec3 pos)
{
    double mass = 0;
    for (auto f : face_around)
    {
#ifdef CACHE_GRADIENT
        auto pts = dsc->get_pos(*dsc->get_nodes_cache(f));
#else
        auto pts = dsc->get_pos(dsc->get_nodes(f));
#endif
        auto v = Util::volume<real>(pos, pts[0], pts[1], pts[2]);
        mass += v * 1.0;
    }
    
    return mass;
}

void UI::build_gradient_mass()
{
    profile t("build gradient");
    
    for (auto nk = dsc->nodes_begin(); nk != dsc->nodes_end(); nk++)
    {
        if (!nk->is_interface())
        {
#ifdef CACHE_GRADIENT
            auto tets = *dsc->get_tets_cache(nk.key());
#else
            auto tets = dsc->get_tets(nk.key());
#endif
            if (dsc->get_label(tets[0]) == 0)
            {
                continue; // This node is outside of the object
            }
            
#ifdef CACHE_GRADIENT
            auto faces = *dsc->get_link(nk.key());
#else
            auto faces = dsc->get_faces(tets) - dsc->get_faces(nk.key());
#endif
            
            auto p = dsc->get_pos(nk.key());
            auto mass0 = mass(faces, p);
            vec3 grad(0.0);
            for (int i = 0; i < 3; i++)
            {
                auto mass1 = mass(faces, p + vec3(0.01, 0, 0));
                grad[i] = (mass1 - mass0);
            }
        }
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
        basic_log = std::shared_ptr<Log>(new Log(log_path + log_folder_name));
//        painter->set_view_position(camera_pos);
        painter->save_painting(log_path, vel_fun->get_time_step());
        basic_log->write_message(vel_fun->get_name().c_str());
        basic_log->write_log(*vel_fun);
        basic_log->write_log(*dsc);
    }
    
//    painter->update(*dsc);
    update_title();
    glutPostRedisplay();
}
