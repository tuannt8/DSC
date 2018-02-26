//
//  fix_DSC.h
//  DEMO
//
//  Created by Tuan Nguyen Trung on 24/02/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#ifndef fix_DSC_h
#define fix_DSC_h
#include "DSC.h"
#include <string>
#include "tetralizer.h"
#include "KDTree.h"

void fix_DSC(std::string file_name, double edge_length, vec3 domain_size, bool gap_bound = true)
{
    // Load DSC
    std::cout << "\nLoading " <<file_name << std::endl;
    
    std::vector<vec3> points;
    std::vector<int>  tets;
    std::vector<int>  tet_labels;
    is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
    
    DSC::DeformableSimplicialComplex<> dsc(points, tets, tet_labels);
    dsc.validity_check();
    
    // Export face
    // Should we make sure there is no face intersection?
    std::vector<vec3> f_points;
    std::vector<int> f_faces;
    dsc.extract_surface_mesh(f_points, f_faces);
    
    vec3 gap(edge_length);
    if(gap_bound)
    {
        for(auto & p : f_points)
            p += gap;
        domain_size = domain_size + 2*gap;
    }
    
    // Fix face intersection
    
    
    // Gen mesh
    std::cout<<"Start generation " << std::endl;
    std::vector<vec3> new_points;
    std::vector<int> new_tets;
    std::vector<int> new_labels;
    Tetralizer::tetralize(domain_size, edge_length, f_points, f_faces, new_points, new_tets, new_labels);
    if(gap_bound)
    {
        for(auto & p : new_points)
            p -= gap;
    }
    
    // fix the label
    DSC::DeformableSimplicialComplex<> dsc_new(new_points, new_tets, new_labels);
    Geometry::KDTree<vec3, int> m_vtree;
    for (auto nit = dsc.nodes_begin(); nit != dsc.nodes_end(); nit++)
    {
        m_vtree.insert(nit->get_pos(), nit.key());
    }
    m_vtree.build();
    
    for(auto tit = dsc_new.tetrahedra_begin(); tit != dsc_new.tetrahedra_end(); tit++)
    {
        auto mid_pos = dsc_new.barycenter(tit.key());
        int closest_point;
        vec3 closest_pos;
        double dist = edge_length*2;
        m_vtree.closest_point(mid_pos, dist, closest_pos, closest_point);
        
        is_mesh::NodeKey cp(closest_point);
        bool found = false;
        for(auto t : dsc.get_tets(cp))
        {
            auto nodes_t = dsc.get_pos(dsc.get_nodes(t));
            if (Util::is_inside(mid_pos, nodes_t[0], nodes_t[1], nodes_t[2], nodes_t[3]))
            {
                found = true;
                dsc_new.set_label(tit.key(), dsc.get_label(t));
                break;
            }
        }
        
        assert(found);
    }
    
    // Write new DSC
    std::string new_name = file_name.substr(file_name.size()-4) + "_new.dsc";
    {
        std::cout << " export to " << new_name << std::endl;
        std::vector<vec3> points;
        std::vector<int> tets;
        std::vector<int> tet_labels;
        dsc_new.extract_tet_mesh(points, tets, tet_labels);
        is_mesh::export_tet_mesh(new_name, points, tets, tet_labels);
    }
}

#endif /* fix_DSC_h */
