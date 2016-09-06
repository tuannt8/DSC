//
//  DSC_parallel.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 2/22/16.
//  Copyright © 2016 Asger Nyman Christiansen. All rights reserved.
//

#include <stdio.h>
#include "DSC.h"

using namespace std;
typedef DSC::DeformableSimplicialComplex<> dsc_class;


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
    
//    is_mesh::SimplexSet<tet_key> * tids = get_tets_cache(nid);
//    is_mesh::SimplexSet<face_key> fids = get_faces(*tids) - *get_faces_cache(nid);
    
    auto fids = get_link(nid);

//    t.change("get pos cache");
    
    vec3 old_pos = get_pos(nid);
    
    vec3 avg_pos = get_barycenter(*get_nodes_cache(nid));
    vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);

    real q_old, q_new;
    
//    t.change("quality");
    min_quality(*fids, old_pos, new_pos, q_old, q_new);
       
#else
    
    
//    profile t("smooth get");
    is_mesh::SimplexSet<tet_key> tids1 = get_tets(nid);
    is_mesh::SimplexSet<face_key> fids1 = get_faces(tids1) - get_faces(nid);

//    t.change("get pos");
    
    vec3 old_pos = get_pos(nid);

    vec3 avg_pos = get_barycenter(get_nodes(fids1));
    vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);
    
    
    real q_old, q_new;
    
//    t.change("quality");
    
    min_quality(fids1, old_pos, new_pos, q_old, q_new);
        
    
#endif
    
//    t.change("set pos");
    if(q_new > pars.MIN_TET_QUALITY || q_new > q_old)
    {
        set_pos(nid, new_pos);
        return true;
    }
    return false;

}

//
//
//#define NUM_THREADS 8
//#define p_min(a, b) (a<b?a:b)
//std::bitset<10000000> _dirty;
//std::mutex dsc_lock;
//
//std::atomic<int> total(0), processed(0);
//
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
//template<> void dsc_class::topological_edge_removal_worker(dsc_class *dsc, std::vector<dsc_class::tet_key> *tet_list, int start_idx, int stop_idx, Barrier & bar)
//{
//    
//    std::vector<tet_key> tets;
//    for (int i = start_idx; i < stop_idx; i++)
//    {
//        auto cur_tet = tet_list->at(i);
//        
//        if ( dsc->quality(cur_tet) < dsc->pars.MIN_TET_QUALITY)
//        {
//            tets.push_back(cur_tet);
//        }
//    }
//    
//    vector<edge_key> edge_topo_remove;
//    vector<is_mesh::SimplexSet<edge_key>> edge_topo_remove_neighbor;
//    vector<edge_key> edge_topo_boundary_remove;
//    vector<is_mesh::SimplexSet<edge_key>> edge_topo_boundary_remove_neighbor;
//    for (auto &t:tets)
//    {
//        if (dsc->is_unsafe_editable(t) && dsc->quality(t) < dsc->pars.MIN_TET_QUALITY)
//        {
//            for (auto e : dsc->get_edges(t))
//            {
//                if(dsc->is_safe_editable(e))
//                {
//                    edge_topo_remove.push_back(e);
//                    
//                    is_mesh::SimplexSet<edge_key> neightbor = affected_neighbor_topo_edge_remove(dsc, e);
//                    edge_topo_remove_neighbor.push_back(neightbor);
//                }
//                else if(dsc->exists(e)
//                        && (dsc->get(e).is_interface() || dsc->get(e).is_boundary())
//                        && dsc->is_flippable(e))
//                {
//                    edge_topo_boundary_remove.push_back(e);
//                    is_mesh::SimplexSet<edge_key> neightbor = affected_neighbor_topo_edge_remove(dsc, e);
//                    edge_topo_boundary_remove_neighbor.push_back(neightbor);
//                }
//            }
//        }
//    }
//    
//    bar.Wait();
//    
//    for (int i = 0; i < edge_topo_remove.size(); i++)
//    {
//        auto &ekey = edge_topo_remove[i];
//        dsc_lock.lock();
//        if (!_dirty[ekey])
//        {
//            auto & ns = edge_topo_remove_neighbor[i];
//            for (auto ee : ns)
//            {
//                _dirty[ee] = 1;
//            }
//            dsc_lock.unlock();
//            
//            dsc->topological_edge_removal(ekey);
//        }
//        else
//            dsc_lock.unlock();
//    }
//    
//    for (int i = 0; i < edge_topo_boundary_remove.size(); i++)
//    {
//        auto &ekey = edge_topo_boundary_remove[i];
//        dsc_lock.lock();
//        if (!_dirty[ekey])
//        {
//            auto & ns = edge_topo_boundary_remove_neighbor[i];
//            for (auto ee : ns)
//            {
//                _dirty[ee] = 1;
//            }
//            dsc_lock.unlock();
//            
//            dsc->topological_boundary_edge_removal(ekey);
//        }
//        else
//            dsc_lock.unlock();
//    }
//}
//
//template<> void dsc_class::topological_edge_removal_parallel()
//{
//    _dirty.reset();
//    total = 0;
//    processed = 0;
//    
//    // Divide the work
//    std::vector<tet_key> tetra_list;
//    for (auto tit = tetrahedra_begin(); tit != tetrahedra_end(); tit++)
//    {
//        tetra_list.push_back(tit.key());
//    }
//    auto num_faces = tetra_list.size();
//    
//    // Launch threads
//    std::thread ths[NUM_THREADS];
//    int work_share = (int)(num_faces / NUM_THREADS) + 1;
//    Barrier bar(NUM_THREADS);
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i] = std::thread(topological_edge_removal_worker, this, &tetra_list, i*work_share, p_min(i*work_share+work_share, num_faces), std::ref(bar) );
//    }
//    
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i].join();
//    }
//}
//
//template<> void dsc_class::smooth_worker(dsc_class *dsc, std::vector<dsc_class::node_key> *node_list, int start_idx, int stop_idx, Barrier & bar)
//{
//    double alpha = 1.0;
//    
//    vector<vec3> node_new_pos;
//    vector<node_key> nodes;
//    for (int i = start_idx; i < stop_idx; i++)
//    {
//        auto nid = node_list->at(i);
//        
//        if (dsc->is_safe_editable(nid) )
//        {
//            
//            is_mesh::SimplexSet<tet_key> tids = dsc->get_tets(nid);
//            is_mesh::SimplexSet<face_key> fids = dsc->get_faces(tids) - dsc->get_faces(nid);
//            
//            vec3 old_pos = dsc->get_pos(nid);
//            vec3 avg_pos = dsc->get_barycenter(dsc->get_nodes(fids));
//            vec3 new_pos = old_pos + alpha * (avg_pos - old_pos);
//            
//            real q_old, q_new;
//            dsc->min_quality(fids, old_pos, new_pos, q_old, q_new);
//            if (q_new > dsc->pars.MIN_TET_QUALITY || q_new > q_old)
//            {
//                nodes.push_back(nid);
//                node_new_pos.push_back(new_pos);
//            }
//        }
//    }
//    
//    bar.Wait();
//    
//    for (int i = 0; i < nodes.size(); i++)
//    {
//        dsc->set_pos(nodes[i], node_new_pos[i]);
//    }
//}
//
//template<> void dsc_class:: smooth_parallel()
//{
//    _dirty.reset();
//    
//    vector<node_key> nids;
//    for (auto nid = nodes_begin(); nid != nodes_end(); nid++)
//    {
//        nids.push_back(nid.key());
//    }
//    
//    // Launch threads
//    std::thread ths[NUM_THREADS];
//    int work_share = (int)(nids.size() / NUM_THREADS) + 1;
//    Barrier bar(NUM_THREADS);
//    
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i] = std::thread(smooth_worker, this, &nids, i*work_share, p_min(i+work_share+work_share, nids.size()), std::ref(bar) );
//    }
//    
//    for (int i = 0; i < NUM_THREADS; i++)
//    {
//        ths[i].join();
//    }
//}