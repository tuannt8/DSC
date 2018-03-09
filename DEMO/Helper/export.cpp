//
//  export.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 09/03/2018.
//  Copyright © 2018 Asger Nyman Christiansen. All rights reserved.
//

#include "export.hpp"
#include "define.h"
#include <fstream>

namespace dsc_export
{
    using namespace std;
    using namespace DSC;
    
    dsc_ptr load_dsc(std::string file_name)
    {
        std::vector<vec3> points;
        std::vector<int>  tets;
        std::vector<int>  tet_labels;
        is_mesh::import_tet_mesh(file_name, points, tets, tet_labels);
        
        auto l_dsc = dsc_ptr(new DeformableSimplicialComplex<>(points, tets, tet_labels));

        return l_dsc;
    }
    
    void export_surface(std::string path)
    {
        auto dsc = load_dsc(path);

        vector<bool> phase(11,0);
        for(auto tit = dsc->tetrahedra_begin(); tit != dsc->tetrahedra_end(); tit++)
        {
            phase[min(10, dsc->get_label(tit.key()))] = true;
        }
        
        for (int i = 0; i < 10; i++)
        {
            if(phase[i])
            {
                string path_i = path.substr(0, path.length() - 4) + "_" + std::to_string(i) + ".obj";
                export_surface(dsc, i, path_i);
            }
        }
        
    }
    
    void export_surface(dsc_ptr dsc, int phase, std::string path)
    {
        std::vector<int> points_maps(dsc->get_no_nodes(), -1); // points_maps[dsc_idx] = new_idx
        std::vector<vec3> points;
        std::vector<vec3i> faces;
        
        int cur_point_idx = 0;
        for (auto fit = dsc->faces_begin(); fit != dsc->faces_end(); fit++)
        {
            if (!fit->is_interface())
                continue;
            
            auto tets = dsc->get_tets(fit.key());
            int labels[2] = {dsc->get_label(tets[0]), dsc->get_label(tets[1])};
            
            if (labels[0] != phase && labels[1] != phase)
                continue;
            
            int idx_other = labels[0] == phase ? 0 : 1;
            
            auto nodes = dsc->get_sorted_nodes(fit.key(), tets[idx_other]);
            vec3i face_map;
            for(int i = 0; i < 3; i++)
            {
                auto n = nodes[i];
                if(points_maps[n] == -1)
                {
                    points.push_back(dsc->get_pos(n));
                    points_maps[n] = cur_point_idx++;
                }
                face_map[i] = points_maps[n] + 1;
            }
            faces.push_back(face_map);
        }
        
        ofstream file(path);
        
        for(auto p : points)
            file << "v " << p[0] << " " << p[1] << " " << p[2] << endl;
        for (auto f : faces)
        {
            file << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
        }
    }
}
