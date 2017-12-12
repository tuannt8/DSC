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
        // Project point
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
         return vec3(0.15, 0.15, 0.15);
    }
    
    virtual double get_influence_radius()
    {
        return 0.0064;
//        return 0.0034;
    }
    
   virtual void personal_draw()
    {

    }
    
    virtual double get_spacing_distance()
    {
        return 0.003;
    }
};
#endif /* bubble_h */
