//
//  dam_break.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 30/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef dam_break_h
#define dam_break_h

#include "file_load.hpp"
#include <math.h>

class dam_break:public file_load
{
public:
    dam_break()
    {
#ifdef __APPLE__
        m_data_path = "../Large_data/DamBreak3D/my_format/iter_";
#else
        m_data_path = "../../Large_data/DamBreak3D/my_format/iter_";
#endif
        load_time_step();

    }
    ~dam_break()
    {
        
    }
    
public:
    virtual void  init_dsc(DSC::DeformableSimplicialComplex<> * dsc)
    {
        vec3 dam_bound(0.4, 0.67, 0.4);
        is_mesh::Cube bound_cube(dam_bound/2.0, dam_bound);
        
        dsc->set_labels(bound_cube, 1);
        
        // Project interface to the Cube
        double thres = dsc->get_avg_edge_length()/1.3;
        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                auto p = nit->get_pos();
                for (int i = 0; i < 3; i++)
                {
                    if (std::abs(p[i]) < thres)
                    {
                        p[i] = 0.0;
                    }
                    
                    if (std::abs(dam_bound[i] - p[i]) < thres)
                    {
                        p[i] = dam_bound[i];
                    }
                }
                
                dsc->set_destination(nit.key(), p);
            }
        }
        
//        // Fix boundary
//        // 4 planes
//        vec3 center(0.9, 0.67/2.0, 0);
//        double size = 0.06;
//        vec3 point[4];
//        vec3 norm[4];
//        for(int i = 0; i < 4; i++)
//        {
//            double alpha = i*M_PI_2 + M_PI_4;
//            norm[i] = vec3(cos(alpha), sin(alpha), 0.0);
//            point[i] = center + norm[i]*size;
//        }
//
//        std::vector<int> interface_boundary(dsc->get_no_nodes(), 0);
//        for (auto tit = dsc->tetrahedra_begin(); tit!=dsc->tetrahedra_end(); tit++)
//        {
//            auto center_p = dsc->barycenter(tit.key());
//            bool is_outside = false;
//            for (int i = 0; i < 4; i++)
//            {
//                double sign = Util::dot(norm[i], center_p - point[i]);
//                if(sign >= 0)
//                    is_outside = true;
//            }
//            if (!is_outside)
//            {
//                dsc->set_label(tit.key(), 2);
//                auto nodes = dsc->get_nodes(tit.key());
//                for (int idx = 0; idx < nodes.size(); idx++)
//                {
//                    interface_boundary[nodes[idx]] = 1;
//                }
//            }
//        }
//        // Project
//        double gap_thres = dsc->get_avg_edge_length()/2;
//        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
//        {
//            if (nit->is_interface() && interface_boundary[nit.key()]==1)
//            {
//                // rotate it
//
//            }
//        }
        
        dsc->deform();
    }
    
    virtual double get_spacing_distance()
    {
        return 0.025;
//        return 0.052;
    }
    
    virtual vec3 get_domain_dimension()
    {
        return vec3(1.6, 0.67, 0.6);
    }
    
    virtual double get_influence_radius()
    {
        return 0.052;
    }
    
    virtual void personal_draw()
    {
//        static std::vector<vec3> _color = {vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)};
//        
//        
//        glDisable(GL_LIGHTING);
//        glPointSize(2.5);
//        glBegin(GL_POINTS);
//        for (auto &p : m_current_particles)
//        {
//            glColor3dv(_color[p.type].get());
//            glVertex3dv(p.pos.get());
//        }
//        glEnd();
    }
};

#endif /* dam_break_h */
