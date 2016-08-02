//
//  cache.hpp
//  DSC_parallel
//
//  Created by Tuan Nguyen Trung on 8/1/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#ifndef cache_hpp
#define cache_hpp

#include <stdio.h>
#include <vector>
#include <bitset>
#include "util.h"
#include "is_mesh.h"

class dsc_cache
{
private:
    std::bitset<10000000> vertex_flag_clean;
    std::bitset<10000000> edge_flag_clean;
    std::bitset<10000000> face_flag_clean;
    std::bitset<10000000> tet_flag_clean;
    
public:
    dsc_cache();
    ~dsc_cache(){};
    
    void mark_dirty(is_mesh::NodeKey nk, bool dirty)
    {
        vertex_flag_clean[nk] = dirty;
    }
    
    void mark_dirty(is_mesh::EdgeKey ek, bool dirty)
    {
        edge_flag_clean[ek] = dirty;
    }
    
    void mark_dirty(is_mesh::FaceKey fk, bool dirty)
    {
        face_flag_clean[fk] = dirty;
    }
    
    void mark_dirty(is_mesh::TetrahedronKey tk, bool dirty)
    {
        tet_flag_clean[tk] = dirty;
    }
    
    bool is_dirty(is_mesh::NodeKey nk)
    {
        return !vertex_flag_clean[nk];
    }
    
    bool is_dirty(is_mesh::EdgeKey ek)
    {
        return !edge_flag_clean[ek];
    }
    
    bool is_dirty(is_mesh::FaceKey fk)
    {
        return !face_flag_clean[fk];
    }
    
    bool is_dirty(is_mesh::TetrahedronKey tk)
    {
        return tet_flag_clean[tk];
    }
};

#endif /* cache_hpp */
