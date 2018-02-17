//
//  dambreak2.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 01/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef dambreak2_h
#define dambreak2_h
class dam_break2:public file_load
{
public:
    dam_break2()
    {
#ifdef __APPLE__
        m_data_path = "../Large_data/DamBreak3D/my_format_no_obstacle/iter_";
#else
        m_data_path = "../../Large_data/DamBreak3D/my_format_no_obstacle/iter_";
#endif
        load_time_step();
    }
    ~dam_break2()
    {
        
    }
    
public:
    virtual void  init_dsc(DSC::DeformableSimplicialComplex<> * dsc)
    {
        // Init first cube
        vec3 dam_bound(0.4, 0.67, 0.4);
        is_mesh::Cube bound_cube(dam_bound/2.0, dam_bound);
        
        dsc->set_labels(bound_cube, 1);
        
        // Project interface to the Cube
        std::vector<int> first_bound(dsc->get_no_nodes(), 0);
        double thres = dsc->get_avg_edge_length()/1.3;
        for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
        {
            if (nit->is_interface())
            {
                first_bound[(int)nit.key()] = 1;
                
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
        
        // Second bound
        {
            vec3 dam_bound(0.4, 0.67, 0.4);
            vec3 origin(0.6, 0.0, 0.0);
            is_mesh::Cube bound_cube(dam_bound/2.0 + origin, dam_bound);
            
            dsc->set_labels(bound_cube, 1);
            
            // Project interface to the Cube
            double thres = dsc->get_avg_edge_length()/1.3;
            
            for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
            {
                if (nit->is_interface())
                {
                    if(first_bound[nit.key()] == 1)
                        continue;
                    
                    auto p = nit->get_pos();
                    for (int i = 0; i < 3; i++)
                    {
                        if (std::abs(p[i] - origin[i]) < thres)
                        {
                            p[i] = origin[i];
                        }
                        
                        if (std::abs(dam_bound[i] + origin[i] - p[i]) < thres)
                        {
                            p[i] = dam_bound[i] + origin[i] ;
                        }
                    }
                    
                    dsc->set_destination(nit.key(), p);
                }
            }
        }
        
        dsc->deform();
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
    
    virtual double get_spacing_distance()
    {
        return 0.025;
    }
};

#endif /* dambreak2_h */
