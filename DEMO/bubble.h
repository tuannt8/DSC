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
        for(auto tit = dsc->tetrahedra_begin(); tit!=dsc->tetrahedra_end(); tit++)
        {
            auto center_pos = dsc->barycenter(tit.key());
            if((center - center_pos).length() < R)
            {
                dsc->set_label(tit.key(), 1);
            }
        }
        
        // Project point
        for(auto nit = dsc->nodes_begin(); nit!= dsc->nodes_end(); nit++)
        {
            if(nit->is_interface())
            {
                auto pos = nit->get_pos();
                auto r_pos = pos - center;
                r_pos.normalize();
                vec3 new_pos = center + r_pos*R;
                dsc->set_destination(nit.key(), new_pos);
            }
        }

        
        dsc->deform();
    }
    
    virtual vec3 get_domain_dimension()
    {
        return vec3(0.17, 0.17, 0.27);
    }
    
    virtual double get_influence_radius()
    {
        return 0.0064;
    }
    
   virtual void personal_draw()
    {
        glColor3d(1,0,0);
        glDisable(GL_LIGHTING);
        glPointSize(2.5);
        glBegin(GL_POINTS);
        for (auto &p : m_current_particles)
        {
            glVertex3dv(p.pos.get());
        }
        glEnd();
    }
};
#endif /* bubble_h */
