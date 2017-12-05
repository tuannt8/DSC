//
//  bubble.h
//  DEMO
//
//  Created by Tuan Nguyen Trung on 30/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef bubble_h
#define bubble_h

#include "file_load.hpp"

class bubble: public file_load
{
public:
    bubble(){
#ifdef __APPLE__
        m_data_path = "../Large_data/Bubble/my_format/iter_";
#else
        m_data_path = "../../Large_data/Bubble/my_format/iter_";
#endif
        load_time_step();
    }
    ~bubble(){};
    
    virtual void  init_dsc(DSC::DeformableSimplicialComplex<> * dsc)
    {
        vec3 center(0.0848462, 0.0848462, 0.05);
        double R = 0.025;
//        for(auto tit = dsc->tetrahedra_begin(); tit!=dsc->tetrahedra_end(); tit++)
//        {
//            auto center_pos = dsc->barycenter(tit.key());
//            if((center - center_pos).length() < R)
//            {
//                dsc->set_label(tit.key(), 1);
//            }
//        }
        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
        {
            auto node_pos = nit->get_pos();
            if( (node_pos-center).length() < R )
            {
                auto tits = dsc->get_tets(nit.key());
                for(int i = 0; i < tits.size(); i++)
                {
                    dsc->set_label(tits[i], 1);
                }
            }
        }
        
//        // Shake
//        double ll = dsc->get_avg_edge_length()*0.3;
//        for(auto nit = dsc->nodes_begin(); nit!= dsc->nodes_end(); nit++)
//        {
//            if(nit->is_interface())
//            {
//                auto pos = nit->get_pos();
//
//                auto new_pos = pos + vec3(0.5 + (double)rand()/RAND_MAX, 0.5 + (double)rand()/RAND_MAX, 0.5 + (double)rand()/RAND_MAX)*ll;
//                dsc->set_destination(nit.key(), new_pos);
//            }
//        }
//
//        dsc->deform();
//
//        // Project point
//        for(auto nit = dsc->nodes_begin(); nit!= dsc->nodes_end(); nit++)
//        {
//            if(nit->is_interface())
//            {
//                auto pos = nit->get_pos();
//                auto r_pos = pos - center;
//
//                r_pos.normalize();
//                vec3 new_pos = center + r_pos*R;
//
//                dsc->set_destination(nit.key(), new_pos);
//            }
//        }
//
//        dsc->deform();
    }
    
    virtual vec3 get_domain_dimension()
    {
//        return vec3(0.17, 0.17, 0.27);
         return vec3(0.10, 0.10, 0.10);
    }
    
    virtual double get_influence_radius()
    {
        return 0.0064;
    }
    
   virtual void personal_draw()
    {
        std::vector<int> idx_list = {98};
        
        if(glut_menu::get_state("particle points", 1))
        {
            glColor3d(1,0,0);
            glDisable(GL_LIGHTING);
            glPointSize(6);
            glBegin(GL_POINTS);
//            int idx = 0;
//            for (int idx : idx_list)
            for(int idx = 0; idx < m_current_particles.size(); idx++)
            {
                auto &p = m_current_particles[idx];
                
                glVertex3dv(p.pos.get());
//                idx++;
//                if(idx>100)
//                    break;
            }
            glEnd();
        }
        if(glut_menu::get_state("Principle component axis", 1))
        {
            std::vector<vec3> axis = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
            
            double ra = 0.0025*0.4;
//            int idx = 0;
            //            for (int idx : idx_list)
            for(int idx = 0; idx < m_current_particles.size(); idx++)
            {
                auto &p = m_current_particles[idx];
                
                auto pos = p.pos;
                
                glPushMatrix();
                glTranslated(pos[0], pos[1], pos[2]);
                
                auto G = m_aniso_kernel.m_G[idx];
                double m[16] = {G[0][0], G[0][1], G[0][2], 0,
                    G[1][0], G[1][1], G[1][2], 0,
                    G[2][0], G[2][1], G[2][2], 0,
                    0, 0, 0, 1
                };
                glMultMatrixd(m);
                
                glDisable(GL_LIGHTING);
                glBegin(GL_LINES);
                for (auto a : axis)
                {
                    glColor3f(a[0], a[1], a[2]);
                    glVertex3f(0, 0, 0);
                    glVertex3dv((a*ra).get());
                }
                glEnd();
                
                
                glPopMatrix();
//
//                idx++;
//                if(idx>100)
//                    break;
            }
        }
        if(glut_menu::get_state("Principle component", 1))
        {
            double ra = 0.0025*0.4;
//            int idx = 0;
            //            for (int idx : idx_list)
            for(int idx = 0; idx < m_current_particles.size(); idx++)
            {
                auto &p = m_current_particles[idx];
                auto pos = p.pos;
                
                glPushMatrix();
                glTranslated(pos[0], pos[1], pos[2]);
                
                auto G = m_aniso_kernel.m_G[idx];
                double m[16] = {G[0][0], G[0][1], G[0][2], 0,
                    G[1][0], G[1][1], G[1][2], 0,
                    G[2][0], G[2][1], G[2][2], 0,
                    0, 0, 0, 1
                };
                glMultMatrixd(m);
                
                glEnable(GL_LIGHTING);
                glColor3f(1, 0, 0);
                glutSolidSphere(ra, 10, 10);
                
                
                glPopMatrix();
//
//                idx++;
//                if(idx>10)
//                    break;
            }
        }
    }
};
#endif /* bubble_h */
