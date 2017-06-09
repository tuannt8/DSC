//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"
#include "tet_dis_coord.hpp"
#include <queue>

#include "profile.h"

#ifdef _WIN32 // WINDOWS
#include <GL/glut.h>
#include <GL/glew.h>
#elif defined(__APPLE__) // IOS
#include <OpenGL/gl3.h>
#include <GLUT/glut.h>
#else // LINUX
#include <GL/glew.h>
#include <GL/glut.h>
#endif


using namespace std;


std::bitset<4> X_direction("0001");
std::bitset<4> Y_direction("0010");
std::bitset<4> Z_direction("0100");

void segment_function::init()
{
#if defined(__APPLE__) || defined(_WIN32)
    
#else
    _directory_path = std::string("../") + _directory_path;
#endif
    _img.load(_directory_path);
    
    cout << "Done loading " << _directory_path << endl;
}

void segment_function::random_initialization()
{
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if (_dsc->get_label(tit.key()) != BOUND_LABEL)
        {
            int new_label = (rand() % 3);
            _dsc->set_label(tit.key(), new_label);
        }
    }
}

// Threshold the histogram in a range to 3 regions
// Multilevel Thresholding for Image Segmentation through a Fast Statistical Recursive Algorithm
// Arora
void Obstu_threshold(std::vector<int> & thres_T1, std::vector<int> & thres_T2, std::vector<int> & input, int iter)
{
    if (iter == 0)
    {
        return;
    }
    
    int range_min = thres_T1[thres_T1.size()-1];
    int range_max = thres_T2[0];
    
    double thres_cofficient = 1;
    // Compute mean
    int mean = 0, count=0;
    for (int i = range_min; i < range_max; i++)
    {
        mean += input[i]*i;
        count += input[i];
    }
    mean = mean/count;
    
    // Compute deviation
    double deviation = 0;
    for (int i = range_min; i < range_max; i++)
    {
        deviation += input[i]*std::abs(i - mean);
    }
    deviation = (deviation) / count;
    
    int K1 = mean - thres_cofficient*deviation;
    int K2 = mean + thres_cofficient*deviation;
    
    assert(K1 > 0 && K2 < range_max);
    
    thres_T1.push_back(K1);
    thres_T2.insert(thres_T2.begin(), K2);
}

std::vector<int> obstu_recursive(std::vector<int> input, int nb_phase)
{
    std::vector<int> thres_T1; thres_T1.push_back(0);
    std::vector<int> thres_T2; thres_T2.push_back(255);
    int nb_iter = int((nb_phase-1)/2);
    if (nb_iter == 0)
    {
        nb_iter = 1;
    }
    Obstu_threshold(thres_T1, thres_T2, input, nb_iter);
    
    if (nb_phase %2 == 0)
    {
        // remove one threshold
        thres_T2.erase(thres_T2.begin());
    }
    
    thres_T1.insert(thres_T1.end(), thres_T2.begin(), thres_T2.end());
    return thres_T1;
}

void segment_function::initialization_discrete_opt()
{
    cout << "Relabeling " << endl;
    
    
    double portion_keep_for_opt = 0.7; // Depend how sparse the segmenting
    
    
    // Optimize the labels of the tetrahedral
#ifdef _DSC_ORIGIN_
    int no_tets = MAX_NUM_ELEMENT_MESH;
#else
    int no_tets = _dsc->get_no_tets_buffer();
#endif
    std::vector<double> total_intensity_per_tet(no_tets, -1.0);
    std::vector<double> volume_per_tet(no_tets, -1.0);
    std::vector<double> mean_inten_per_tet(no_tets, -1.0);
    std::vector<double> variation_inten_per_tet(no_tets, -1.0);
    std::vector<int> labels(no_tets, -1);
    
    double max_variation = -INFINITY, min_variation = INFINITY;
    // 1. Compute mean intensity and variation of each tetrahedron
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if(_dsc->get_label(tit.key()) == BOUND_LABEL)
            continue;
#ifdef DSC_CACHE
        auto tet_nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tit.key()));
#else
        auto tet_nodes_pos = _dsc->get_pos(_dsc->get_nodes(tit.key()));
#endif
        double total_inten, volume;
        auto mean_inten = _img.get_tetra_intensity(tet_nodes_pos, &total_inten, &volume);
        
        auto variation = _img.get_variation(tet_nodes_pos, mean_inten);
        
        mean_inten_per_tet[tit.key()] = mean_inten;
        total_intensity_per_tet[tit.key()] = total_inten;
        volume_per_tet[tit.key()] = volume;
        variation_inten_per_tet[tit.key()] = variation;
        
        if (max_variation < variation)
        {
            max_variation = variation;
        }
        if (min_variation > variation)
        {
            min_variation = variation;
        }
    }
    
    // 2. Remove tetrahedra that have high variationsimage. More sparse, more tetrahedra to be optimized
    // Use bin size 100 for histogram count
    int bin_size = 200;
    max_variation *= 1.01;
    min_variation *= 0.99;
    double his_step = (max_variation - min_variation)/(double)bin_size;
    std::vector<int> his_count(bin_size, 0);
    long total_count = 0;
    for (auto & v : variation_inten_per_tet)
    {
        if(v > 0)
        {
            int step = (int)((v - min_variation)/his_step);
            his_count[step] ++;
            total_count++;
        }
    }
    
    // Find the threshold
    int index_for_thres = 0;
    long count_cur = 0;
    for (; index_for_thres < his_count.size(); index_for_thres++)
    {
        count_cur += his_count[index_for_thres];
        if (count_cur > total_count*portion_keep_for_opt)
        {
            break;
        }
    }
    double thres_hold = (index_for_thres+1)*his_step; // thres hold to remove high variation tets
    
    // make new list, tets with low variation
    std::vector<int> histogram_for_thresholding(256,0);
    for (int i = 0; i < variation_inten_per_tet.size(); i++)
    {
        if (variation_inten_per_tet[i] > 0
            && variation_inten_per_tet[i] < thres_hold)
        {
            // This tetrahedron is considered for relabeling
            int idx = (int)(mean_inten_per_tet[i]*256); // convert from [0, 1] to [0, 256] image
            
            if(idx > 255) idx = 255;
            histogram_for_thresholding[idx] ++;
        }
    }
    
    // 3. Optimize the label
    // 3.1. Random initialize
    int nb_phases = NB_PHASE;
    vector<int> thres_hold_array = obstu_recursive(histogram_for_thresholding, nb_phases);
    
    // debuging
    cout << "Thresholding with: " << thres_hold_array[1] << "; "
            << thres_hold_array[2] << endl;
    
    // Initialize the label
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if (_dsc->get_label(tit.key()) == BOUND_LABEL
            || variation_inten_per_tet[tit.key()] > thres_hold)
        {
            continue;
        }
        
        int mean_inten_tet = (int)(mean_inten_per_tet[tit.key()]*255);
        if(mean_inten_tet >= 255) mean_inten_tet = 254;
        
        int label = 0;
        for (; label < thres_hold_array.size(); label++)
        {
            if (mean_inten_tet < thres_hold_array[label])
            {
                break;
            }
        }
        assert(label < NB_PHASE+1);
        _dsc->set_label(tit.key(), label-1); // -1 because the thres_array start from threshold = 0
    }
}

void segment_function::initialze_segmentation()
{
//    /**
//     Hamster sample
//     */
//    // Initilization by thresholding
//    double thres = 0.6;
//    // initialize by thresholding
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        // test
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//        
//        assert(avgI < 1.01);
//        //assert(total < 1000);
//        
//        if (avgI > thres)
//        {
//            _dsc->set_label(tit.key(), 1);
//        }
//    }
    
    
//    // Analyzing
//    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
//    {
//        double total, volume;
//        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
//        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
//    }
    
    /**
     Fuel cells
     */
    // Initialization by thresholding
    double thres[] = {0.31, 0.57, 0.7};
    for (auto tit = _dsc->tetrahedra_begin(); tit != _dsc->tetrahedra_end(); tit++)
    {
        if(_dsc->get_label(tit.key()) == BOUND_LABEL)
            continue;
        
        // test
        double total, volume;
        auto pts = _dsc->get_pos(_dsc->get_nodes(tit.key()));
        double avgI = _img.get_tetra_intensity(pts, &total, &volume);
        
        double var = _img.get_variation(pts, avgI);
        if (var > 0.1)
        {
            continue;
        }

        assert(avgI < 1.01);
        
        // find closest
        int idx = 0;
        if (std::abs(avgI - thres[1]) < 0.1 )
        {
            idx = 1;
        }else if(std::abs(avgI - thres[2]) < 0.1 )
        {
            idx = 2;
        }

        if (idx != 0)
        {
            _dsc->set_label(tit.key(), idx);
        }
    }
}

void segment_function::update_vertex_stability()
{
    _vertex_stability_map = std::vector<int>(30000, 0); // suppose we have less than 10000 vertices
    
    for (auto vid = _dsc->nodes_begin(); vid != _dsc->nodes_end(); vid++)
    {
        if (_forces[(long)vid.key()].length() < 0.1) // stable
        {
            _vertex_stability_map[vid.key()] = 1;
        }
    }
}



inline std::bitset<4> get_direction(vec3 a)
{
    std::bitset<4> d("0000");
    if (std::abs(a[0]) > 0.8)
    {
        d = d | X_direction;
    }
    if (std::abs(a[1]) > 0.8)
    {
        d = d | Y_direction;
    }
    if (std::abs(a[2]) > 0.8)
    {
        d = d | Z_direction;
    }
    
    return d;
}

double segment_function::get_energy_tetrahedron(is_mesh::TetrahedronKey tkey, int assumed_label)
{
#ifdef DSC_CACHE
    auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tkey));
#else
    auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tkey));
#endif
    auto old_label = _dsc->get_label(tkey);
    auto energy = _img.get_variation(nodes_pos, _mean_intensities[assumed_label]);
    for (auto fid : _dsc->get_faces(tkey))
    {
        auto cobound_tets = _dsc->get_tets(fid); // We assume that this face is not DSC boundary, as there is a gap between DSC boundary and the image domain
        auto label0 = _dsc->get_label(cobound_tets[0]);
        auto label1 = _dsc->get_label(cobound_tets[1]);
        
        label0 = label0==old_label? label1 : label0;
        
        if (label0 != assumed_label)
        {
            energy += _dsc->area(fid)*ALPHA;
        }
    }
    
    return energy;
}

void segment_function::adapt_tetrahedra()
{
    // Split tetrahedron if there is something inside
    
}

void segment_function::relabel_tetrahedra()
{
    cout << "Relabeling --" << endl;
    
    for (auto tid = _dsc->tetrahedra_begin(); tid != _dsc->tetrahedra_end(); tid++)
    {
        if (_dsc->get_label(tid.key()) == BOUND_LABEL)
        {
            continue;
        }
#ifdef DSC_CACHE
        auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tid.key()));
#else
        auto nodes_pos = _dsc->get_pos(_dsc->get_nodes(tid.key()));
#endif
        double volume, total_inten;
        auto mean_inten_tetra = _img.get_tetra_intensity(nodes_pos, &total_inten, &volume);
        
        
        
        // We check if it is worth changing the label to the phase with closest mean intensity
        double smallest_gap = INFINITY;
        int label_of_closest_phase = -1;
        for (int i = 0; i < _mean_intensities.size(); i++)
        {
            if (smallest_gap > mean_inten_tetra - _mean_intensities[i])
            {
                smallest_gap = mean_inten_tetra - _mean_intensities[i];
                label_of_closest_phase = i;
            }
        }
        
        if (label_of_closest_phase != _dsc->get_label(tid.key()))
        {
            // Check if we reduce the energy
            auto old_energy = get_energy_tetrahedron(tid.key(), _dsc->get_label(tid.key()));
            auto new_energy = get_energy_tetrahedron(tid.key(), label_of_closest_phase);
            
            if (new_energy < old_energy) //Should we a factor here to make sure that the benifit is enough?
            {
                _dsc->set_label(tid.key(), label_of_closest_phase);
                
                // update mean intensity
                int old_label = _dsc->get_label(tid.key());
                _total_intensities[old_label] -= total_inten;
                _total_intensities[label_of_closest_phase] += total_inten;
                _phase_volume[old_label] -= volume;
                _phase_volume[label_of_closest_phase] += volume;
                
                for (int i = 0; i < _mean_intensities.size(); i++)
                {
                    _mean_intensities[i] =_total_intensities[i] / _phase_volume[i];
                }
            }
        }
    }
}

void segment_function::work_around_on_boundary_vertices()
{

    
    // 1. Find boundary vertices
#ifdef _DSC_ORIGIN_
    int node_mem_size = MAX_NUM_ELEMENT_MESH;
#else
    auto node_mem_size = _dsc->get_no_nodes_buffer();
#endif
    std::vector<unsigned int> is_bound_vertex(node_mem_size,0);
    std::vector<std::bitset<4>> direction_state(node_mem_size,std::bitset<4>("0000"));
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() && !fid->is_boundary())
        {
            auto tets = _dsc->get_tets(fid.key());
            if(_dsc->get_label(tets[0]) == BOUND_LABEL
               || _dsc->get_label(tets[1]) == BOUND_LABEL)
            {
                auto nodes_on_face = _dsc->get_nodes(fid.key());
                auto norm = _dsc->get_normal(fid.key());
                std::bitset<4> direction = get_direction(norm);
                
                for (auto n : nodes_on_face)
                {
//                    bound_nodes += n; //
                    is_bound_vertex[(unsigned int)n] = 1;
                    direction_state[(unsigned int)n] = direction_state[(unsigned int)n] | direction;
                }
            }
        }
    }
    
    // For debuging
    d_direction_state = direction_state;
    d_is_image_boundary = is_bound_vertex;
    
    // 2. Align the boundary vertices to the boundary
    auto domain_dim = _img.dimension();
    auto max_displacement = 4;//_dsc->get_avg_edge_length()*1.5;
    
    double max_displacement_real = -INFINITY;
    
    // For debuging
    boundary_vertices_displacements = std::vector<vec3>(node_mem_size, vec3(0));
    
    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {
        
        if ( (nid->is_interface() or nid->is_crossing())
            && _dsc->exists(nid.key())
            and !nid->is_boundary())
        {
//            auto dis = _forces[nid.key()]*_dt;
//            auto dis = _internal_forces[nid.key()]*1;
            
            auto dis = (_internal_forces[nid.key()]*ALPHA + _forces[nid.key()] + _quality_control_forces[nid.key()]*QALPHA)*_dt;
            
//            auto dis = (_internal_forces[nid.key()]*ALPHA + _forces[nid.key()])*_dt;
            
            assert(!isnan(dis.length()));
            
            // limit it
            if (dis.length() > max_displacement)
            {
                dis = Util::normalize(dis)*max_displacement;
                cout << "Force is too large" << endl;
            }
            
            
            vec3 destination = nid->get_pos() + dis;
            auto threshold = _dsc->get_avg_edge_length();
            // Align boundary
            if (is_bound_vertex[nid.key()])
            {
                auto direct = direction_state[nid.key()];
                auto pos = nid->get_pos();
                if ((direct & X_direction).to_ulong() != 0)
                {
                    destination[0] = (std::abs(pos[0]) < threshold)? 0 : domain_dim[0] - 1;
                }
                if ((direct & Y_direction).to_ulong() != 0){
                    destination[1] = (std::abs(pos[1]) < threshold)? 0 : domain_dim[1] -1;
                }
                if ((direct & Z_direction).to_ulong() != 0){
                    destination[2] = (std::abs(pos[2]) < threshold)? 0 : domain_dim[2] -1;
                }
                
                boundary_vertices_displacements[nid.key()] = destination - nid->get_pos();
                
                dis = destination - nid->get_pos();
            }
            
            _dsc->set_destination(nid.key(), destination);
            
            if (max_displacement_real < dis.length())
            {
                max_displacement_real = dis.length();
            }
        }
    }
    
    cout << "Max displacement: " << max_displacement_real << endl;
}

void segment_function::compute_mesh_quality_control_force()
{
    // Adaptive force based on curvature
    std::vector<vec3> adaptive_force(_dsc->get_no_nodes_buffer(), vec3(0));
    
    for (auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (nit->is_interface() && !nit->is_crossing() && !nit->is_boundary())
        {
#ifdef DSC_CACHE
            auto nodes_around = *_dsc->get_nodes_cache(nit.key());
#else
            auto tids1 = get_tets(nid);
            auto fids1 = get_faces(tids1) - get_faces(nid);
            auto nodes_around = get_nodes(fids1);
#endif
            
            // It may be more complicated
            
            auto node_curvature = _internal_forces[nit.key()].length();
            int count = 0;
            vec3 ff(0.0);
            for (int idx = 0; idx < nodes_around.size(); idx++)
            {
                auto curvature =  _internal_forces[nodes_around[idx]].length();
                if(curvature > node_curvature)
                {
                    auto direct = _dsc->get_pos(nodes_around[idx]) - _dsc->get_pos(nit.key());
//                    double shortest_edge = _dsc->pars.MIN_EDGE_QUALITY*_dsc->AVG_LENGTH;
//                    if (direct.length() > shortest_edge)
//                    {
//                        direct = Util::normalize(direct)*(direct.length() - shortest_edge);
//                    }
                    ff += direct*(curvature - node_curvature);
                    count++;
                    
                    // What is its counter forces?
                }
            }
            
            if(count > 0)
            {
                ff = ff *0.1/ count;
                adaptive_force[nit.key()] = ff;
            }
        }
    }
    
    // Smooth force based on triange angle
    // Maximize minimal angle
    std::vector<vec3> angle_force(_dsc->get_no_nodes_buffer(), vec3(0));
    std::vector<int> count(_dsc->get_no_nodes_buffer(), 0);
    for (auto fit = _dsc->faces_begin(); fit != _dsc->faces_end(); fit++)
    {
        if (fit->is_interface() && !fit->is_boundary())
        {
#ifdef DSC_CACHE
            is_mesh::SimplexSet<is_mesh::NodeKey> nids = *_dsc->get_nodes_cache(fit.key());
#else
            is_mesh::SimplexSet<node_key> nids = _dsc->get_nodes(fit.key());
#endif
            static double cos_60 = cos(3.14159/3.0);
            auto nodes_pos = _dsc->get_pos(nids);
            
            for (int i = 0; i < 3; i++)
            {
                int i1 = (i+1)%3, i2 = (i+2)%3;
                auto cos_angle = Util::cos_angle<real>(nodes_pos[i], nodes_pos[i1], nodes_pos[i2]);
                
                vec3 f = ((nodes_pos[i1] + nodes_pos[i2])/2. - nodes_pos[i])*(cos_angle - cos_60);
                
                angle_force[nids[i]] += f;
                count[nids[i]]++;
            }
        }
    }
    // angle forces should be perpendicular with normal
    for (int i = 0; i < angle_force.size(); i++)
    {
        if(count[i] > 0)
        {
            auto n = _internal_forces[i];
            if (n.length() > 0.1)
            {
                n.normalize();
                auto f = angle_force[i];
                vec3 f1 = n*Util::dot(n, f);
                vec3 f2 = f - f1;
                
                angle_force[i] = f2/count[i];
            }
            
            angle_force[i] *= 0.1;

//            _quality_angle_forces[i] = angle_force[i];
        }
    }
    _quality_angle_forces = angle_force;
    _quality_control_forces = adaptive_force;
    
    for (int i = 0; i < _quality_control_forces.size(); i++)
    {
        _quality_control_forces[i] += _quality_angle_forces[i];
    }
}

// Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
void segment_function::compute_internal_force()
{
 
    std::vector<std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>>> interface_faces_around_node(_dsc->get_no_nodes_buffer(), std::vector<is_mesh::SimplexSet<is_mesh::FaceKey>>());
    
    
    // Find all connected region around nodes
    for(auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if(!(nit->is_interface() && !nit->is_boundary()))
            continue;
            
        
        auto tets = _dsc->get_tets(nit.key());
        
        while (tets.size() != 0)
        {
            std::queue<is_mesh::FaceKey> growing_queue;
            is_mesh::SimplexSet<is_mesh::TetrahedronKey> region;
            region.push_back(tets.back());
            tets -= tets.back();
            
            int label = _dsc->get_label(region.back());
            
            {
            auto fs = _dsc->get_faces(region.back());
            for(auto f : fs)
                growing_queue.push(f);
            }
            
            while (!growing_queue.empty())
            {
                auto f = growing_queue.front();
                growing_queue.pop();
                
                auto tets_f = _dsc->get_tets(f);
                //find
                for (auto t : tets_f)
                {
                    if (_dsc->get_label(t) != label)
                    {
                        continue;
                    }
                    
                    int idx = tets.index(t);
                    if (idx != -1)
                    {
                        region.push_back(t);
                        tets -= t;
                        
                        { // add to region
                            auto fs = _dsc->get_faces(t);
                            for(auto f : fs)
                                growing_queue.push(f);
                        }
                    }
                }
            }
            
            if (label != BOUND_LABEL)
            {
                is_mesh::SimplexSet<is_mesh::FaceKey> faceset = _dsc->get_faces(region) & _dsc->get_faces(nit.key()); // optimize later
                is_mesh::SimplexSet<is_mesh::FaceKey> interface_faceset;
                for (auto f : faceset)
                {
                    if (_dsc->get(f).is_interface())
                    {
                        interface_faceset.push_back(f);
                    }
                }
                
                interface_faces_around_node[nit.key()].push_back(interface_faceset);
            }
        }
    }
    
    //-----------------
    // Now compute the mean curvature
//    std::vector<double> node_cur(_dsc->get_no_nodes_buffer(), 0);
    
    std::vector<vec3> mean_curvature_norm_a(_dsc->get_no_nodes_buffer(), vec3(0));
    
    auto dsc = _dsc;
    
    for(auto n = dsc->nodes_begin(); n!= dsc->nodes_end(); n++)
    {
        if(n->is_interface() && !n->is_boundary())
        {
            // 1. Build nodes around
            
            if((long)n.key() == 1787)
            {
                
            }
            
            auto face_around = interface_faces_around_node[n.key()];
            
            for (auto fs : face_around)
            {
                auto node_around = dsc->extract_node_around(n.key(), fs);
                
                if (!node_around)
                {
                    continue;
                }
                
                // 2. Compute the curvature
                auto p = dsc->get_pos(n.key());
                // 2a. Mixed area
                //            vec3 normal(0.0);
                double area_mixed = 0;
                for (int i = 0; i < node_around->size(); i++)
                {
                    vec3 v0 = p;
                    vec3 v1 = dsc->get_pos((*node_around)[i]);
                    vec3 v2 = dsc->get_pos( (*node_around)[(i+1)% (node_around->size()) ] );
                    
                    double f_area = Util::area<real>(v0, v1, v2);
                    
                    double a0 = acos(dot(v1-v0, v2-v0)/(length(v1-v0)*length(v2-v0)));
                    double a1 = acos(dot(v2-v1, v0-v1)/(length(v2-v1)*length(v0-v1)));
                    double a2 = acos(dot(v0-v2, v1-v2)/(length(v0-v2)*length(v1-v2)));
                    
                    if(a0>(M_PI/2.0) && a1>(M_PI/2.0) && a2>(M_PI/2.0)) // f is non-obtuse
                    {
                        // Add Voronoi formula (see Section 3.3)
                        area_mixed += (1.0/8) *
                        ((1.0/tan(a1)) * sqr_length(v2-v0) +
                         (1.0/tan(a2)) * sqr_length(v1-v0));
                    }
                    else // Voronoi inappropriate
                    {
                        // Add either area(f)/4 or area(f)/2
                        area_mixed += f_area/3;
                    }
                }
                // 2b. unnormalize curvature normal
                auto curv_normal = vec3(0);
                for (int i = 0; i < node_around->size(); i++)
                {
                    auto nbr = dsc->get_pos((*node_around)[i]);
                    auto right = dsc->get_pos((*node_around)[(i+1)%node_around->size()]);
                    auto left = dsc->get_pos((*node_around)[( i-1 + node_around->size() ) %node_around->size()]);
                    
                    double d_left = Util::dot(Util::normalize(nbr-left), Util::normalize(p-left));
                    double d_right = Util::dot(Util::normalize(nbr-right), Util::normalize(p-right));
                    double cos_left = std::min(1.0, std::max(-1.0, d_left));
                    double cos_right = std::min(1.0, std::max(-1.0, d_right));
                    double sin_left = sqrt(1 - cos_left*cos_left);
                    double sin_right = sqrt(1 - cos_right*cos_right);
                    
                    double w = (sin_left*cos_right + sin_right*cos_left)/(1e-300 + sin_left*sin_right);
                    
                    curv_normal += w * (nbr-p);
                }
                
                // 2c. The curvature
                auto mean_curvature_norm = curv_normal / (4*area_mixed);
                
                mean_curvature_norm_a[n.key()] += mean_curvature_norm;
                
                
                delete node_around;
                
            }
            
         }
    }
    
    _internal_forces = mean_curvature_norm_a;
    
    // should we normalize it?
    double max_f = -INFINITY;
    for (auto f : _internal_forces)
    {
        if (max_f < f.length())
        {
            max_f = f.length();
        }
    }
    
    double scale = 5./max_f;
    for (auto & f : _internal_forces)
    {
        f *= scale;
    }
    cout << "Normalize curvature; scale = " << scale << endl;
    
//    // Now compute the real force
//    // Normalize the force
//    for(int i = 0; i < forces.size(); i++)
//    {
//        if (forces[i].length() > EPSILON )
//        {
//            forces[i] = Util::normalize(forces[i]) * node_cur[i];
//            
//            assert(!isnan(forces[i].length()));
//        }
//    }
//    
//    _internal_forces = forces;
}

void segment_function::compute_external_force()
{
    auto c = _mean_intensities;
    
    // Array to store temporary forces
    // Use fixxed array for better performance. Suppose that we have less than 10000 vertices
    std::vector<vec3> forces = std::vector<vec3>(_dsc->get_no_nodes_buffer(), vec3(0.0));
    
    // Loop on interface faces
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() && !fid->is_boundary())
        {
            auto tets = _dsc->get_tets(fid.key());
            if (_dsc->get_label(tets[0]) == BOUND_LABEL ||
                _dsc->get_label(tets[1]) == BOUND_LABEL )
            {
                // Ignore the faces on the boundary.
                //   These boundary vertices should only move along the boundary
                //   Forces faces on boundary will make them move perpendicullar to the boundary
                continue;
            }
            
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            double c0 = c[_dsc->get_label(tets[0])];
            double c1 = c[_dsc->get_label(tets[1])];
            
            // get normal
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            size_t tri_sample_index = std::ceil( sqrt(area) ) - 1;
            if (tri_sample_index >= tri_coord_size.size())
            {
                tri_sample_index = tri_coord_size.size() - 1;
            }
            if(tri_sample_index < 1)tri_sample_index = 1;
            
            auto a = tri_dis_coord[tri_sample_index - 1];
            
            if(a.size() > 1 && a.size() > _dsc->area(fid.key()))
            {
                std::cout << "--Error in triangle decomposition. Each sample point represents the area < 1 pixel^2 ----: "
                            << a.size() << " - " << _dsc->area(fid.key()) << std::endl;
            }
            
            for (auto coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);

                auto f = - Norm* ((2*g - c0 - c1) / (c1-c0) / area);
                
                // distribute
                forces[verts[0]] += f*coord[0];
                forces[verts[1]] += f*coord[1];
                forces[verts[2]] += f*coord[2];
            }
        }
    }

    _forces = forces;
}

void segment_function::face_split()
{
    // try to adapt the flat surface
    compute_internal_force();
    
    for(auto nit = _dsc->nodes_begin(); nit != _dsc->nodes_end(); nit++)
    {
        if (_dsc->exists(nit.key()) && nit->is_interface() && !nit->is_crossing() && !nit->is_boundary())
        {
            if (_internal_forces[nit.key()].length() < 0.1)
            {
                // collapsing edge
                auto edges = _dsc->get_edges(nit.key());
                int min_idx = -1;
                double shortest_length = INFINITY;
                for (int i = 0; i < edges.size(); i++)
                {
                    if (_dsc->get(edges[i]).is_interface() && !_dsc->get(edges[i]).is_boundary())
                    {
                        if (shortest_length > _dsc->length(edges[i]))
                        {
                            shortest_length = _dsc->length(edges[i]);
                            min_idx = i;
                        }
                    }
                }
                
                if(min_idx != -1)
                {
                    auto other_node = (_dsc->get_nodes(edges[min_idx]) - nit.key()).front();
                    _dsc->collapse(edges[min_idx]);
                }
            }
        }
    }
    
    return;
    
    auto c = _mean_intensities;
    
    std::vector<int> related_edges(_dsc->get_no_edges_buffer(), 0);
    
    int i = 0, j = 0;
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (_dsc->exists(fid.key()) && fid->is_interface())
        {
            auto tets = _dsc->get_tets(fid.key());
            if (_dsc->get_label(tets[0]) == BOUND_LABEL ||
                _dsc->get_label(tets[1]) == BOUND_LABEL )
            {
                // Ignore the faces on the boundary.
                // This ignorian should be considered later
                // We may also assess the boundary faces
                continue;
            }
            
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            double c0 = c[_dsc->get_label(tets[0])];
            double c1 = c[_dsc->get_label(tets[1])];
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            size_t tri_sample_index = std::ceil( sqrt(area) ) - 1;
            if (tri_sample_index >= tri_coord_size.size())
            {
                tri_sample_index = tri_coord_size.size() - 1;
            }
            if(tri_sample_index < 1)tri_sample_index = 1;
            
            auto a = tri_dis_coord[tri_sample_index - 1];
            
            if(a.size() > 1 && a.size() > _dsc->area(fid.key()))
            {
                std::cout << "--Error in triangle decomposition. Each sample point represents the area < 1 pixel^2 ----: "
                << a.size() << " - " << _dsc->area(fid.key()) << std::endl;
            }
            
            std::vector<double> forces_distribute(a.size(), 0.0);
            int idx = 0;
            
            for (auto coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);
                
                auto f = - ((2*g - c0 - c1) / (c1-c0) / area);
                
                forces_distribute[idx++] = f;
            }
            
            /*
             Analyze the force distribution
             */
            double magnitude_mean = 0, signed_mean = 0;
            for (auto & force : forces_distribute)
            {
                magnitude_mean += std::abs(force);
                signed_mean += force;
            }
            
            magnitude_mean /= forces_distribute.size();
            signed_mean /= forces_distribute.size();
            
            if(std::abs(signed_mean) < ratio_signed_and_mag_mean * magnitude_mean )
            {
                // This face cover an inhomogineous area
                // Consider split it
                _dsc->split_face(fid.key());
//                for(auto ee : _dsc->get_edges(fid.key()))
//                {
//                    related_edges[ee] = 1;
//                }
                i++;
                
            }
            j++;
        }
    }
    
    cout << "Adaptive split " << i << "/" << j << " triangles" << endl;
    
    _dsc->remove_interface_faces();

//    struct edge_sort{
//        double length;
//        is_mesh::EdgeKey ekey;
//    };
//    // Split the related edge, longest first
//    std::vector<edge_sort> edges_list;
//    for(int i = 0; i < related_edges.size(); i++)
//    {
//        if (related_edges[i] == 1)
//        {
//            edges_list.push_back({_dsc->length(is_mesh::EdgeKey(i)), is_mesh::EdgeKey(i)});
//        }
//    }
//    
//    std::sort(edges_list.begin(), edges_list.end(), [](edge_sort const& a, edge_sort const& b){return a.length > b.length;});
//    
//    for (auto ee : edges_list)
//    {
//        if (related_edges[ee.ekey] == 1 && _dsc->exists(ee.ekey))
//        {
//            auto edges_around = _dsc->get_edges(_dsc->get_faces(ee.ekey));
//            for(auto ea : edges_around)
//            {
//                related_edges[ea] = 0;
//            }
//            
//            _dsc->split(ee.ekey);
//        }
//    }
}


/**
 * Bounding box of point list
 */
void bounding_box(const std::vector<vec3> & pts, vec3 & ld, vec3 & ru)
{
    ld = vec3(INFINITY);
    ru = vec3(-INFINITY);
    for(auto & v : pts)
    {
        ld[0] = std::min(v[0], ld[0]);
        ld[1] = std::min(v[1], ld[1]);
        ld[2] = std::min(v[2], ld[2]);
        
        ru[0] = std::max(v[0], ru[0]);
        ru[1] = std::max(v[1], ru[1]);
        ru[2] = std::max(v[2], ru[2]);
    }
}

bool sort_intersect(intersect_pt p1, intersect_pt p2)
{
    {
        return p1.z < p2.z;
    }
}


void segment_function::update_average_intensity()
{
    // Using sum table. Much faster than normal loop
#ifdef LOG_DEBUG
    cout << "Computing average intensity with" << NB_PHASE << " phases " << endl;
#endif
    
    int nb_phase = NB_PHASE;
    
    // 1. Init the buffer for intersection
    auto dim = _img.dimension();

    std::vector<ray_z> init_rayz(dim[0] * dim[1]);
    for (int y = 0; y < dim[1]; y++)
    {
        for (int x = 0; x < dim[0]; x ++)
        {
            int idx = y*dim[0] + x;
            init_rayz[idx].x = x;
            init_rayz[idx].y = y;
        }
    }
    
    vector<std::vector<ray_z>> ray_intersect(nb_phase, init_rayz);

    // 2. Find intersection with interface
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface() and !fid->is_boundary())
        {
            auto tet = _dsc->get_tets(fid.key());
            auto phase0 = _dsc->get_label(tet[0]);
            auto phase1 = _dsc->get_label(tet[1]);
            
            // check all z-ray that intersect this triangle
            auto pts3 = _dsc->get_pos(_dsc->get_nodes(fid.key()));
            auto pts = pts3;
            pts[0][2] = 0; pts[1][2] = 0; pts[2][2] = 0;
            
            auto n = _dsc->get_normal(fid.key(), tet[0]);
            bool in_1 = Util::dot(n, vec3(0,0,1)) > 0;
            
            vec3 ld, ru;
            bounding_box(pts, ld, ru);
            for (int x = std::floor(ld[0]); x < std::round(ru[0]); x++)
            {
                for (int y = std::floor(ld[1]); y < std::round(ru[1]); y++)
                {
                    if(y < 0 || x < 0
                       || y >= dim[1] || x >= dim[0])
                        continue;
                    
                    try
                    {
                        bool bError;
                        auto bc = Util::barycentric_coords<double>(vec3(x+0.5, y+0.5, 0), pts[0], pts[1], pts[2], &bError);
                        
                        if (bError)
                        {
                            continue;
                        }
                        
                        if (bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                        { // inside
                            auto p = pts3[0]*bc[0] + pts3[1]*bc[1] + pts3[2]*bc[2];
                            
                            auto zz = std::nearbyint(p[2]);
                            
                            if(phase0 != BOUND_LABEL)
                            {
                                ray_intersect[phase0][y*dim[0] + x].intersects.push_back(intersect_pt(zz, !in_1));
                                
                            }
                            if(phase1 != BOUND_LABEL)
                            {
                                ray_intersect[phase1][y*dim[0] + x].intersects.push_back(intersect_pt(zz, in_1));
                                
                            }
                        }
                    }
                    catch (std::exception e)
                    {
                    
                    }
                }
            }
        }
    }

    // 3. Compute integral
    _d_rayz.clear();
    
    std::fill(_mean_intensities.begin(), _mean_intensities.end(), 0.0);
    vector<double> area(_mean_intensities.size(), 0.0);
    for (int i = 0; i < nb_phase; i++)
    {
        int count = 0;
        for (auto r : ray_intersect[i])
        {
            std::vector<intersect_pt> intersect_ps = r.intersects;
            if (intersect_ps.size() > 1)
            {
                std::sort(intersect_ps.begin(), intersect_ps.end(), sort_intersect);
                
                // remove identical intersections
                for (auto p = intersect_ps.begin()+1; p != intersect_ps.end(); p++)
                {
                    auto pre = p -1;
                    if (p->z == pre->z and
                        p->b_in == pre->b_in)
                    {
                        p = intersect_ps.erase(p);
                        
                        if(p == intersect_ps.end())
                        {
                            break;
                        }
                    }
                }
                // Now count
                if (intersect_ps.size() % 2 == 0) // Should be even
                {
                    count++;
                    auto newR = r;
                    newR.intersects = intersect_ps;
                    _d_rayz.push_back(newR);
                    
                    for (int j = 0; j < intersect_ps.size()/2; j++)
                    {
                        int z1 = intersect_ps[2*j].z;
                        int z2 = intersect_ps[2*j + 1].z;
                        
                        if(z1 < 0) z1 = 0;
                        if(z2 < z1)
                            z2 = z1;
                        
                        area[i] += z2 - z1;
                        _mean_intensities[i] += _img.sum_line_z(r.x, r.y, z1, z2);
                    }
                }
                else{}//???
            }
        }
    }
    
    _total_intensities = _mean_intensities;
    _phase_volume = area;
    
    for (int i  = 0; i < nb_phase; i++)
    {
        if(area[i]==0)
        {
            std::cout << "Zero; " << _mean_intensities[i] << std::endl;
        }
        _mean_intensities[i] /= area[i];
    }

    cout << "Mean intensity: ";
    for (int i = 0; i < nb_phase; i++)
    {
        cout << _mean_intensities[i] << " ; ";
    }
    cout << endl;

}

void segment_function::segment()
{
    int num_phases = NB_PHASE;
    
    static int iteration = 0;
    cout << "--------------- Iteration " << iteration++ << " ----------------" << endl;
    
    // 1. Compute average intensity
    profile t("Segment time");
    
    _mean_intensities.resize(num_phases);
    update_average_intensity();
    
    

    
    // 2. Compute external force
    compute_external_force();
    compute_internal_force();
    compute_mesh_quality_control_force(); // Must be after internal force
    
    // 3. Work around to align boundary vertices
    //  including set displacement for interface vertices
    work_around_on_boundary_vertices();
    
    _dsc->deform();
    
    /**
     4. RELABEL TETRAHEDRA
     */
    if (iteration % 5 == 0)
    {
//        face_split();
        relabel_tetrahedra();
    }
    
    t.done();
    profile::close();
}
