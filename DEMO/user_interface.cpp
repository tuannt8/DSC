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


using namespace DSC;


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
    
//    vec3 eye = center + vec3(gl_dis_max*2.0*cos(angle)*cos(angle2),
//                             gl_dis_max*2.0*cos(angle)*sin(angle2),
//                             gl_dis_max*2.0*sin(angle));
    
//    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//    GLfloat mat_shininess[] = { 50.0 };
//    GLfloat light_position[] = { -(GLfloat)eye[0], -(GLfloat)eye[1], -(GLfloat)eye[2], 0.0 };
//    glClearColor(1.0, 1.1, 1.0, 1.0);
    
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 100.0 };
    GLfloat light_position[] = { -(GLfloat)eye[0], -(GLfloat)eye[1], -(GLfloat)eye[2], 0.0 };
    glClearColor(1.0, 1.1, 1.0, 1.0);
    
    glShadeModel(GL_SMOOTH);
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    
    glEnable(GL_COLOR_MATERIAL);
    
//    {
//        GLfloat white[] = {0.8f, 0.8f, 0.8f, 1.0f};
//        GLfloat cyan1[] = {0.f, .8f, .8f, 1.0f};
//        GLfloat cyan[] = {0.f, .8f, .8f, 0.8f};
//        glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan); // other
//        glMaterialfv(GL_FRONT, GL_SPECULAR, cyan1); // shinny
//        GLfloat shininess[] = {2};
//        glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
//        
//        
//        glEnable(GL_LIGHTING);
//        glEnable(GL_LIGHT0);
//        glEnable(GL_DEPTH_TEST);
//        
//    }

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
    
    vec3 center = vec3(0.0);//_obj_dim / 2.0;
    double dis = 1.5;
    
    vec3 eye = center + vec3(gl_dis_max*dis*cos(angle)*cos(angle2),
                             gl_dis_max*dis*sin(angle),
                             gl_dis_max*dis*cos(angle)*sin(angle2)
                             );
    vec3 head = vec3(-sin(angle)*cos(angle2),
                     cos(angle),
                     -sin(angle)*sin(angle2)
                     );
    
//    vec3 eye = center + vec3(gl_dis_max*dis*cos(angle)*cos(angle2),
//                             gl_dis_max*dis*cos(angle)*sin(angle2),
//                             gl_dis_max*dis*sin(angle));
//    vec3 head = vec3(-sin(angle)*cos(angle2),
//                     -sin(angle)*sin(angle2),
//                     cos(angle));
    gluLookAt(eye[0], eye[1], eye[2], /* eye is at (0,8,60) */
              center[0], center[1], center[2],      /* center is at (0,8,0) */
              head[0], head[1], head[2]);      /* up is in postivie Y direction */
    
    int size = std::min(WIN_SIZE_Y, WIN_SIZE_X);
    glViewport((WIN_SIZE_X - size) / 2.0, (WIN_SIZE_Y - size) / 2.0, size, size);
    
    glClearColor(0.1, 0.1, 0.1, 1.0);
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
    
    
    
    load_model("armadillo", 2.5);
    
    real velocity = 5.;
    real accuracy = 0.25;
    vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(velocity, accuracy, 500));
    start("");
    
    //
    vec3 minc(INFINITY);
    vec3 maxc(-INFINITY);
    for (auto vit = dsc->nodes_begin(); vit != dsc->nodes_end(); vit++)
    {
        vec3 pt = dsc->get_pos(vit.key());
        
        minc = min_vec(minc, pt);
        maxc = max_vec(maxc, pt);
    }
    _obj_dim = maxc - minc;
    gl_dis_max = std::max(std::max(_obj_dim[0], _obj_dim[1]), _obj_dim[2])*2;
    
//    build_node_curvature();
    
//    load_model("armadillo", 2.5);
//    real velocity = 5.;
//    real accuracy = 0.25;
//    vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(velocity, accuracy, 500));
//    stop();
//    QUIT_ON_COMPLETION = false;
//    RECORD = false;
//    vel_fun = std::unique_ptr<VelocityFunc<>>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
//    start("smooth");
//    
//    
//    
//    for (int i = 0; i < 20; i++)
//    {
//        std::cout << "Iter " << i << std::endl;
//        vel_fun->take_time_step(*dsc);
//    }
//    
//    profile::close();
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
    
//    vec3 p_min(INFINITY), p_max(-INFINITY);
//    for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++) {
//        for (int i = 0; i < 3; i++) {
//            p_min[i] = Util::min(nit->get_pos()[i], p_min[i]);
//            p_max[i] = Util::max(nit->get_pos()[i], p_max[i]);
//        }
//    }
//    
//    vec3 size = p_max - p_min;
//    real var = Util::max(Util::max(size[0], size[1]), size[2]);
//    real dist = 1.2*var;
//    eye_pos = {dist, var, dist};
//    camera_pos = {var, var, -dist};
//    light_pos = {0., 0., dist};
    
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
    


    if(mode == 0)
    {
    glColor3f(0.0, 0.9, 1.0);
        std::vector<double> * cur = _node_curvature.size()>0? &_node_curvature : nullptr;
    
    draw_helper::dsc_draw_interface(*dsc, cur);
    
    glColor3f(0.8, 0.8, 0.8);
    draw_helper::dsc_draw_domain(*dsc);
    
//    glColor3f(0.5, 0.5, 0.5);
//    draw_dsc_layer(DISPLAY_LIM);
    }
    else if(mode == 1)
    {
        std::vector<double> * cur = _node_curvature.size()>0? &_node_curvature : nullptr;
        
        glColor3f(0.0, 0.9, 1.0);
        draw_helper::dsc_draw_interface(*dsc, cur);
        
        glColor3f(0.8, 0.8, 0.8);
        draw_helper::dsc_draw_domain(*dsc);
    }
    else if(mode == 2)
    {
        glColor3f(0.8, 0.8, 0.8);
        draw_helper::dsc_draw_domain(*dsc);
        
        glColor3f(0.5, 0.5, 0.5);
        draw_dsc_layer_1(DISPLAY_LIM);
    }
    
    //    draw_helper::dsc_draw_node_color(*dsc);
//    draw_helper::dsc_draw_edges_colors(*dsc);
    
//    // test. Get color the first time
//    static bool bRun = false;
//    if (!bRun)
//    {
//        profile t("get color serial");
//        
//        bRun = true;
//        int max_color = 0;
//        for (auto eit = dsc->edges_begin(); eit != dsc->edges_end(); eit++)
//        {
//            if (max_color < dsc->get_color_edge(eit.key()))
//            {
//                max_color = dsc->get_color_edge(eit.key());
//            }
//        }
//        std::cout << max_color << std::endl;
//        
//        
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
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() + 1 <<  " START*************\n" << std::endl;
        vel_fun->take_time_step(*dsc);

//        compute_volume();
//        build_node_curvature();
//        compute_volume_gradient();
//        build_gradient_mass();
        
#ifdef __APPLE__
        static int log_count = 0;
        std::ostringstream s;
        s << "LOG/scr_" << log_count++ << ".png";
        if(!SOIL_save_screenshot(s.str().c_str(), SOIL_SAVE_TYPE_PNG, 0, 0, WIN_SIZE_X, WIN_SIZE_Y))
        {
            std::cout << "Screen shot fail \n";
        }
#endif
        
        std::cout << "\n***************TIME STEP " << vel_fun->get_time_step() <<  " STOP*************\n" << std::endl;
    }
    glutPostRedisplay();
}



void UI::build_node_curvature()
{
    std::vector<double> node_cur(100000, 0);
    
    profile t("curvature");

    for(auto n = dsc->nodes_begin(); n!= dsc->nodes_end(); n++)
    {
        if(n->is_interface())
        {
            // 1. Build nodes around

            
            auto node_around = dsc->node_on_one_ring_cache(n.key());
            
            // 2. Compute the curvature
            auto p = dsc->get_pos(n.key());
            // 2a. Mixed area
//            vec3 normal(0.0);
            double area_mixed = 0;
            for (int i = 0; i < node_around->size(); i++)
            {
                vec3 v0 = p;
                vec3 v1 = dsc->get_pos((*node_around)[i]);
                vec3 v2 = dsc->get_pos( (*node_around)[(i+1)% (node_around->size()) ] );
                
                double f_area = Util::area<real>(v0, v1, v2);
                
                double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
                double a1 = acos(dot(v2-v1, v0-v1)/(length(v2-v1)*length(v0-v1)));
                double a2 = acos(dot(v0-v2, v1-v2)/(length(v0-v2)*length(v1-v2)));
                
                if(a0>(M_PI/2.0) && a1>(M_PI/2.0) && a2>(M_PI/2.0)) // f is non-obtuse
                {
                    // Add Voronoi formula (see Section 3.3)
                    area_mixed += (1.0/8) *
                    ((1.0/tan(a1)) * sqr_length(v2-v0) +
                     (1.0/tan(a2)) * sqr_length(v1-v0));
                }
                else // Voronoi inappropriate
                {
                    // Add either area(f)/4 or area(f)/2
                    area_mixed += f_area/3;
                }
            }
            // 2b. unnormalize curvature normal
            auto curv_normal = vec3(0);
            for (int i = 0; i < node_around->size(); i++)
            {
                auto nbr = dsc->get_pos((*node_around)[i]);
                auto right = dsc->get_pos((*node_around)[(i+1)%node_around->size()]);
                auto left = dsc->get_pos((*node_around)[( i-1 + node_around->size() ) %node_around->size()]);
                
                double d_left = Util::dot(Util::normalize(nbr-left), Util::normalize(p-left));
                double d_right = Util::dot(Util::normalize(nbr-right), Util::normalize(p-right));
                double cos_left = std::min(1.0, std::max(-1.0, d_left));
                double cos_right = std::min(1.0, std::max(-1.0, d_right));
                double sin_left = sqrt(1 - cos_left*cos_left);
                double sin_right = sqrt(1 - cos_right*cos_right);
                
                double w = (sin_left*cos_right + sin_right*cos_left)/(1e-300 + sin_left*sin_right);
                
//                double a_left  = acos(std::min(1.0, std::max(-1.0, d_left)));
//                double a_right = acos(std::min(1.0, std::max(-1.0, d_right)));
//                
//                double w = sin(a_left + a_right) / (1e-300 + sin(a_left)*sin(a_right));
                
                curv_normal += w * (nbr-p);
            }
            // 2c. The curvature
            auto mean_curvature_norm = curv_normal / (4*area_mixed);
//            auto normn = dsc->get_normal(n.key()); // also cache here

            node_cur[n.key()] = Util::dot(mean_curvature_norm, mean_curvature_norm);

            
//            delete node_around;
        }
    }
    

    _node_curvature = node_cur;
}

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

void UI::compute_volume_gradient()
{
    profile t("volume gradient");
    for (auto nk = dsc->nodes_begin(); nk != dsc->nodes_end(); nk++)
    {
        if (nk->is_interface())
        {
            auto edge_around = dsc->get_edge_around_cache(nk.key());
            
            auto p = dsc->get_pos(nk.key());
            vec3 gradient(0.0);
            for (auto e : *edge_around)
            {
                auto pts = dsc->get_pos(dsc->get_nodes(e));
                
                auto p12 = pts[1] - pts[0];
                auto l = p12.length();
                auto p1p = p - pts[0];
                
                auto cc = Util::cross(p12, p1p);
                auto alpha = cc.length()/l;
                
                auto H = p12*(alpha/l);
                
                gradient += (p-H)*l;
            }
        }
    }
}

void UI::compute_volume()
{
    profile t("Compute volume");
    double v = 0;
    for (auto tk = dsc->tetrahedra_begin(); tk != dsc->tetrahedra_end(); tk++)
    {
        if(dsc->get_label(tk.key()) == 1)
        {
            auto pts = dsc->get_pos(*dsc->get_nodes_cache(tk.key()));
            v += Util::volume<real>(pts[0], pts[1], pts[2], pts[3]);
        }
    }
}

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
            vel_fun = std::unique_ptr<VelocityFunc<>>(new VelocityFunc<>(vel_fun->get_velocity(), vel_fun->get_accuracy(), 500));
            start("");
            break;
        case '1':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new RotateFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("rotate");
            break;
        case '2':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new AverageFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("smooth");
            break;
        case '3':
            stop();
            QUIT_ON_COMPLETION = false;
            RECORD = false;
            vel_fun = std::unique_ptr<VelocityFunc<>>(new NormalFunc(vel_fun->get_velocity(), vel_fun->get_accuracy()));
            start("expand");
            break;
        case ' ':
            if(!CONTINUOUS)
            {
                std::cout << "MOTION STARTED" << std::endl;
//                if(RECORD && basic_log)
//                {
//                    painter->set_view_position(camera_pos);
//                    painter->save_painting(basic_log->get_path(), vel_fun->get_time_step());
//                }
            }
            else {
                std::cout << "MOTION PAUSED" << std::endl;
            }
            CONTINUOUS = !CONTINUOUS;
            break;
        case 'm':
            std::cout << "MOVE" << std::endl;
            vel_fun->take_time_step(*dsc);
            painter->update(*dsc);
            break;
        case 'r':
            std::cout << "RELOAD MODEL" << std::endl;
            load_model(model_file_name, dsc->get_avg_edge_length());
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
            real velocity = std::min(vel_fun->get_velocity() + 1., 100.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '-':
        {
            real velocity = std::max(vel_fun->get_velocity() - 1., 0.);
            vel_fun->set_velocity(velocity);
            update_title();
        }
            break;
        case '.':
        {
            real discretization = std::min(dsc->get_avg_edge_length() + 0.5, 100.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case ',':
        {
            real discretization = std::max(dsc->get_avg_edge_length() - 0.5, 1.);
            dsc->set_avg_edge_length(discretization);
            update_title();
        }
            break;
        case '<':
        {
            real accuracy = std::min(vel_fun->get_accuracy() + 1., 100.);
            vel_fun->set_accuracy(accuracy);
            update_title();
        }
            break;
        case '>':
        {
            real accuracy = std::max(vel_fun->get_accuracy() - 1., 1.);
            vel_fun->set_accuracy(accuracy);
            update_title();
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
        basic_log = std::unique_ptr<Log>(new Log(log_path + log_folder_name));
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
