//
//  DSC_parallel.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/22/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#include <stdio.h>
#include "DSC.h"

#define NUM_THREADS 6

using namespace std;
typedef DSC::DeformableSimplicialComplex<> dsc_class;



#define FIND_MIN_COLOR \
for (int i = 0; i < colors.size(); i++) \
{   \
    if (fc < colors[i]) \
    {   \
        break; \
    }   \
    else if (fc == colors[i]) \
    {   \
        fc++;   \
    }   \
}

template<> int dsc_class::get_free_color(node_key nk)
{
#ifdef DSC_CACHE
    auto ring_nodes = get_nodes(*get_link(nk));
    std::vector<int> colors = get_colors_cache(ring_nodes);
    
    int fc = 0;
    FIND_MIN_COLOR
    
    if (fc >= MAX_COLORS)
    {
        assert(0);
    }
    return fc;
#else
    return 0;
#endif
}

template<> int dsc_class::get_free_color(edge_key ek)
{
    auto nids = get_nodes(ek);
    
//    auto colors = get_colors_cache(get_edges(get_tets(ek)) + get_edges(nids[0]) + get_edges(nids[1]));
    
    auto colors = get_colors_cache(get_edges(get_tets(nids[0]) + get_tets(nids[1])));
    
    int fc = 0;
    FIND_MIN_COLOR
    
    if (fc >= MAX_COLORS_TET)
    {
        assert(0);
    }
    return fc;
}

template<> int dsc_class::get_free_color(tet_key nk)
{
    auto colors = get_colors(get_tets(get_faces(get_tets(get_nodes(nk)))));
    
    int fc = 0;
    FIND_MIN_COLOR
    
    if (fc >= MAX_COLORS_TET)
    {
        assert(0);
    }
    return fc;
}


template<> void dsc_class::init_vertices_color()
{
    for (auto nk = nodes_begin(); nk != nodes_end(); nk++)
    {
        int c = get_free_color(nk.key());
        set_color(nk.key(), get_free_color(nk.key()));
    }
}

void min_quality_worker(DSC::DeformableSimplicialComplex<> * dsc, const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, real * q_old, real * q_new, int idx_start, int idx_stop)
{
    for (int i = idx_start; i < idx_stop; i++)
    {
        auto f = fids[i];
        
        is_mesh::SimplexSet<is_mesh::NodeKey> nids = dsc->get_nodes(f);
        if(Util::sign(Util::signed_volume<real>(dsc->get_pos(nids[0]), dsc->get_pos(nids[1]), dsc->get_pos(nids[2]), pos_old)) !=
           Util::sign(Util::signed_volume<real>(dsc->get_pos(nids[0]), dsc->get_pos(nids[1]), dsc->get_pos(nids[2]), pos_new)))
        {
            q_old[i] = INFINITY;
            q_new[i] = -INFINITY;
            break;
        }
        
        q_old[i] = std::abs(Util::quality<real>(dsc->get_pos(nids[0]), dsc->get_pos(nids[1]), dsc->get_pos(nids[2]), pos_old));
        q_new[i] = std::abs(Util::quality<real>(dsc->get_pos(nids[0]), dsc->get_pos(nids[1]), dsc->get_pos(nids[2]), pos_new));
    }
}

void min_quality_parallel(DSC::DeformableSimplicialComplex<> * dsc, const is_mesh::SimplexSet<is_mesh::FaceKey>& fids, const vec3& pos_old, const vec3& pos_new, real& min_q_old, real& min_q_new)
{
    min_q_old = INFINITY;
    min_q_new = INFINITY;
    
    real q_old[fids.size()];
    real q_new[fids.size()];
    
    int stride = fids.size()/NUM_THREAD_QUALITY + 1;
    std::thread th[NUM_THREAD_QUALITY];
    
    for (int i = 0; i < NUM_THREAD_QUALITY; i++)
    {
        int start = i*stride;
        int stop = MIN_Q((i+1)*stride, fids.size());
        th[i] = std::thread(min_quality_worker, dsc, fids, pos_old, pos_new, &q_old[0], &q_new[0], start, stop);
    }
    
    for (int i = 0; i < NUM_THREAD_QUALITY; i++)
    {
        th[i].join();
    }
    
    for (int i = 0; i < fids.size(); i++)
    {
        min_q_old = Util::min(min_q_old, q_old[i]);
        min_q_new = Util::min(min_q_new, q_new[i]);
    }
    
    if (min_q_new == -INFINITY)
    {
        min_q_old = INFINITY;
    }
}

template<> bool dsc_class::smart_laplacian(const node_key& nid, real alpha)
{
#ifdef DSC_CACHE // laplacian smooth
//    profile t("smooth get cache 1");
    auto fids = get_link(nid);

//    t.change("get pos cache");
    
    vec3 old_pos = get_pos(nid);
    
    vec3 avg_pos = get_barycenter(*get_nodes_cache(nid));
    vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);

    real q_old, q_new;
    
//    for (auto ff : *fids)
//    {
//        is_mesh::SimplexSet<node_key> nids = *get_nodes_cache(ff); // Funny we need this line of code
//        
//        if(nids.size() != 3)
//        {
//            if(ff.is_valid())
//            {
//                printf("Valid face: %d\n", (int)ff);
//                
//                auto nids1 = get_nodes(ff);
//                assert(nids1.size()==3);
//            }
//            
//        }
//    }
    
//    t.change("quality");
    min_quality(*fids, old_pos, new_pos, q_old, q_new);
       
#else
    
    

    is_mesh::SimplexSet<tet_key> tids1 = get_tets(nid);
    is_mesh::SimplexSet<face_key> fids1 = get_faces(tids1) - get_faces(nid);

    
    vec3 old_pos = get_pos(nid);

    vec3 avg_pos = get_barycenter(get_nodes(fids1));
    vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);
    
    
    real q_old, q_new;
    
    min_quality(fids1, old_pos, new_pos, q_old, q_new);
#endif
    
    if(q_new > pars.MIN_TET_QUALITY || q_new > q_old)
    {
        set_pos(nid, new_pos);
        return true;
    }
    return false;

}


//
//
//template<> void dsc_class::topological_face_removal_worker(dsc_class *dsc, std::vector<dsc_class::tet_key> *tet_list, int start_idx, int stop_idx, Barrier & bar)
//{
//    std::vector<tet_key> tets;
//    std::vector<is_mesh::SimplexSet<tet_key> > neighbor;
//    for (int i = start_idx; i < stop_idx; i++)
//    {
//        auto cur_tet = tet_list->at(i);
//        
//        if (dsc->is_unsafe_editable(cur_tet)
//            && dsc->quality(cur_tet) < dsc->pars.MIN_TET_QUALITY)
//        {
//            tets.push_back(cur_tet);
//            // Neightbor
//            auto nids = dsc->get_nodes(cur_tet);
//            is_mesh::SimplexSet<dsc_class::tet_key> neighbor_c;
////            for (auto n : nids)
////            {
////                auto tets_n = dsc->get_tets(n);
////                for (auto tt : tets_n)
////                {
////                    neighbor_c.push_back(tt);
////                }
////            }
//            
//            is_mesh::SimplexSet<dsc_class::node_key> neighbor_node;
//            for (auto n : nids)
//            {
//                auto tets_n = dsc->get_tets(n);
//                for (auto tt : tets_n)
//                {
//                    neighbor_node += dsc->get_nodes(tt);
//                }
//            }
//            for (auto nn : neighbor_node)
//            {
//                neighbor_c += dsc->get_tets(nn);
//            }
//            
//            
//            neighbor.push_back(neighbor_c);
//        }
//    }
//    
//    bar.Wait();
//
//    for (int i = 0; i < tets.size(); i++)
//    {
//        auto t = tets[i];
//        dsc_lock.lock();
//        if (!_dirty[t])
//        {
//            auto & ns = neighbor[i];
//            for (auto & nn : ns)
//            {
//                _dirty[nn] = 1;
//            }
//            dsc_lock.unlock();
//            
//            for (auto f : dsc->get_faces(t))
//            {
//                if (dsc->is_safe_editable(f))
//                {
//                    auto apices = dsc->get_nodes(dsc->get_tets(f)) - dsc->get_nodes(f);
//                    if(dsc->topological_face_removal(apices[0], apices[1]))
//                    {
//                        processed ++;
//                        break;
//                    }
//                }
//            }
//        }
//        else
//            dsc_lock.unlock();
//    }
//}
//
//
//
//template<> void DSC::DeformableSimplicialComplex<>:: topological_face_removal_parallel()
//{
//    _dirty.reset();
//    total = 0;
//    processed = 0;
//    
//    // Divide the work
//    std::vector<tet_key> face_list;
//    for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
//    {
//        face_list.push_back(tit.key());
//    }
//    auto num_faces = face_list.size();
//    
//    // Launch threads
//    std::thread ths[NUM_THREADS];
//    int work_share = (int)(num_faces / NUM_THREADS) + 1;
//    Barrier bar(NUM_THREADS);
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i] = std::thread(topological_face_removal_worker, this, &face_list, i*work_share, p_min(i*work_share+work_share, num_faces), std::ref(bar) );
//    }
//    
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i].join();
//    }
//    
//    cout << "Topo face remove: " << processed << endl;
//}
//
//is_mesh::SimplexSet<dsc_class::edge_key> affected_neighbor_topo_edge_remove(dsc_class *dsc, dsc_class::edge_key eid)
//{
//    is_mesh::SimplexSet<dsc_class::edge_key> neighbors;
//    
//    auto nids = dsc->get_nodes(eid);
//    nids += dsc->get_polygons(eid).front();
//    neighbors += dsc->get_edges(nids);
//    
//    return neighbors;
//}
//

template<> void dsc_class::topological_edge_removal_worker(dsc_class *dsc, is_mesh::SimplexSet<dsc_class::tet_key> *tet_list, int start_idx, int stop_idx)
{
    printf("Edge rm, %d to %d\n", start_idx, stop_idx);
    
    for (int i = start_idx; i < stop_idx; i++)
    {
        tet_key t = (*tet_list)[i];
        
        if (!dsc->exists(t))
        {
            continue;
        }
        
        if (dsc->quality(t) < dsc->pars.MIN_TET_QUALITY
            && dsc->is_unsafe_editable(t))
        {
            auto ets = dsc->get_edges(t);
            for (auto e : ets)
            {
                if(dsc->is_safe_editable(e))
                {
                    if(dsc->topological_edge_removal(e))
                    {
                        break;
                    }
                }
                else
                {
                    bool bf = dsc->is_flippable(e);
                    
                    if(dsc->exists(e) && (dsc->get(e).is_interface() || dsc->get(e).is_boundary()) && bf)
                    {
                        if(dsc->topological_boundary_edge_removal(e))
                        {
                            break;
                        }
                    }
                }
            }
        }
    }
}

template<> void dsc_class::topological_edge_removal_parallel()
{
    std::vector<is_mesh::SimplexSet<tet_key>> colored_list(MAX_COLORS_TET, is_mesh::SimplexSet<tet_key>());
    for(auto n =tetrahedra_begin(); n!=tetrahedra_end(); n++)
    {
        if (quality(n.key()) < pars.MIN_TET_QUALITY)
        {
            int c = get_color_tet(n.key());
            colored_list[c].push_back(n.key());
        }
    }
    
    // Start thread
    // Simple solution first
    int idx = 0;
    for (auto & cc : colored_list)
    {
        if (cc.size() > 0)
        {
            int num_thread = NUM_THREADS;
            if (cc.size() < 100)
            {
                num_thread = 1;
            }
            
            printf("Start threads %d, %d element\n", idx++ , cc.size());
            
            std::thread th[NUM_THREADS];
            int stride = (cc.size())/num_thread + 1;
            for (int i = 0; i < num_thread; i++)
            {
                int start_idx = std::min(i*stride, (int)cc.size()-1);
                int stop_idx = std::min((int)(i+1)*stride-1, (int)cc.size()-1);
                
//                topological_edge_removal_worker(this, &cc, start_idx, stop_idx);
                th[i] = std::thread(topological_edge_removal_worker, this, &cc, start_idx, stop_idx);
            }
            
            for (int i = 0; i < num_thread; i++)
            {
                th[i].join();
            }
            
            printf("Done threads\n");
        }
    }
    
    garbage_collect();
}

template<> void dsc_class::topological_edge_removal_worker1(DeformableSimplicialComplex<> *dsc, is_mesh::SimplexSet<edge_key> *tet_list, int start_idx, int stop_idx)
{
    for (int i = start_idx; i < stop_idx; i++)
    {
        auto e = (*tet_list)[i];
        if (!dsc->exists(e))
        {
            continue;
        }
        
        if(dsc->is_safe_editable(e))
        {
            if(dsc->topological_edge_removal(e))
            {
                break;
            }
        }
        else
        {
            bool bf = dsc->is_flippable(e);
            
            if(dsc->exists(e) && (dsc->get(e).is_interface() || dsc->get(e).is_boundary()) && bf)
            {
                if(dsc->topological_boundary_edge_removal(e))
                {
                    break;
                }
            }
        }
    }
}


template<> void dsc_class::topological_edge_removal_parallel1()
{
    std::vector<is_mesh::SimplexSet<edge_key>> colored_list(MAX_COLORS_TET, is_mesh::SimplexSet<edge_key>());
    
    std::vector<tet_key> tets;
    for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
    {
        if (quality(tit.key()) < pars.MIN_TET_QUALITY)
        {
            tets.push_back(tit.key());
        }
    }
    
    for (auto &t : tets)
    {
        auto ets = get_edges(t);
        for (auto e : ets)
        {
            int cc = get_color_edge(e);
            colored_list[cc] += e;
        }
    }
    
//    topological_edge_removal();
//    return;
    
    int idx = 0;
    for (auto cc : colored_list)
    {
        if (cc.size() > 0)
        {
            int num_thread = NUM_THREADS;
            if (cc.size() < 100)
            {
                num_thread = 1;
            }
            
//                        printf("Start threads %d, %d element\n", idx++ , cc.size());
            
            std::thread th[NUM_THREADS];
            int stride = (cc.size())/num_thread + 1;
            for (int i = 0; i < num_thread; i++)
            {
                int start_idx = std::min(i*stride, (int)cc.size()-1);
                int stop_idx = std::min((int)(i+1)*stride-1, (int)cc.size()-1);
                
                th[i] = std::thread(topological_edge_removal_worker1, this, &cc, start_idx, stop_idx);
            }
            
            for (int i = 0; i < num_thread; i++)
            {
                th[i].join();
            }
            
//                        printf("Done threads\n");
        }
    }
    

}


template<> void dsc_class::smooth_worker(DeformableSimplicialComplex<> *dsc, is_mesh::SimplexSet<node_key> *node_list, int start_idx, int stop_idx)
{
//        printf("Thread launched %d to %d \n", start_idx, stop_idx);
    
    // smooth it
    for (int i = start_idx; i < stop_idx; i++)
    {
        dsc->smart_laplacian((*node_list)[i]);
    }
}

template<> void dsc_class:: smooth_parallel()
{

//Separate them into group
    std::vector<is_mesh::SimplexSet<node_key>> colored_list(MAX_COLORS, is_mesh::SimplexSet<node_key>());
    for(auto n =nodes_begin(); n!=nodes_end(); n++)
    {
        if (is_safe_editable(n.key()))
        {
            int c = get_color_node(n.key());
            colored_list[c].push_back(n.key());
        }
    }

// Start thread
// Simple solution first
   
    int idx = 0;
    for (auto cc : colored_list)
    {
        if (cc.size() > 0)
        {
            int num_thread = NUM_THREADS;
            if (cc.size() < 100)
            {
                num_thread = 1;
            }
            
//            printf("Start threads %d, %d element\n", idx++ , cc.size());
            
            std::thread th[NUM_THREADS];
            int stride = (cc.size())/num_thread + 1;
            for (int i = 0; i < num_thread; i++)
            {
                int start_idx = std::min(i*stride, (int)cc.size()-1);
                int stop_idx = std::min((int)(i+1)*stride-1, (int)cc.size()-1);
                
                th[i] = std::thread(smooth_worker, this, &cc, start_idx, stop_idx);
            }
            
            for (int i = 0; i < num_thread; i++)
            {
                th[i].join();
            }
            
//            printf("Done threads\n");
        }
    }

//    validity_check();
}