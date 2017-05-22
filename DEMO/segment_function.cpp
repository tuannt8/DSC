//
//  segment_function.cpp
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 12/9/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "segment_function.h"
#include "tet_dis_coord.hpp"

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

void segment_function::init()
{
    //_img.load("data/sphere_drill");
//    _img.load("../Large_data/hamster");
    cout << "Load 3D data" << endl;
#ifdef __APPLE__
    _img.load("../Large_data/fuel_cells_smaller");
#else
    _img.load("../../Large_data/fuel_cells_smaller");
#endif
    cout << "Done loading " << endl;
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
    Obstu_threshold(thres_T1, thres_T2, input, int((nb_phase-1)/2));
    
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
    
    // Optimize the labels of the tetrahedral
    int no_tets = _dsc->get_no_tets_buffer();
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
        
        auto tet_nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tit.key()));
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
    
    // 2. Remove tetrahedra that have high variations
    double portion_keep_for_opt = 0.4; // Depend how sparse the segmenting image. More sparse, more tetrahedra to be optimized
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
        
        const int mean_inten_tet = (int)(mean_inten_per_tet[tit.key()]*255);
        int label = 0;
        for (; label < thres_hold_array.size(); label++)
        {
            if (mean_inten_tet < thres_hold_array[label])
            {
                break;
            }
        }
        
        _dsc->set_label(tit.key(), label-1); // -1 because the thres_array start from 0
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

#define X_direction 0x001
#define Y_direction 0x010
#define Z_direction 0x100

inline std::uint8_t get_direction(vec3 a)
{
    std::uint8_t d = 0x0;
    if (std::abs(a[0]) > 0.8)
    {
        d = d & X_direction;
    }
    if (std::abs(a[1]) > 0.8)
    {
        d = d & Y_direction;
    }
    if (std::abs(a[2]) > 0.8)
    {
        d = d & Z_direction;
    }
    
    return d;
}

double segment_function::get_energy_tetrahedron(is_mesh::TetrahedronKey tkey, int assumed_label)
{
#ifdef DSC_CACHE
    auto nodes_pos = _dsc->get_pos(*_dsc->get_nodes_cache(tkey));
#else
    auto nodes_pos = _dsc->get_pos(tkey);
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
        auto nodes_pos = _dsc->get_pos(tid.key());
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
    auto node_mem_size = _dsc->get_no_nodes_buffer();
    std::vector<unsigned int> is_bound_vertex(node_mem_size,0); // Hard code. Assume that the number of vertex < 10000
    std::vector<std::uint8_t> direction_state(node_mem_size,0x0);
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
                std::uint8_t direction = get_direction(norm);
                
                for (auto n : nodes_on_face)
                {
//                    bound_nodes += n; //
                    is_bound_vertex[(unsigned int)n] = 1;
                    direction_state[(unsigned int)n] = direction_state[(unsigned int)n] & direction;
                }
            }
        }
    }
    
    // 2. Align the boundary vertices to the boundary
    auto domain_dim = _img.dimension();
    auto max_displacement = _dsc->get_avg_edge_length()*1.5;
    
    double max_displacement_real = -INFINITY;
    
    for (auto nid = _dsc->nodes_begin(); nid != _dsc->nodes_end(); nid++)
    {
        if ( (nid->is_interface() or nid->is_crossing())
            and _dsc->is_movable(nid.key())
            and !nid->is_boundary())
        {
            auto dis = _forces[nid.key()]*_dt;
            // limit it
            if (dis.length() > max_displacement)
            {
                dis = Util::normalize(dis)*max_displacement;
                cout << "Force is too large" << endl;
            }
            
            
            vec3 destination = nid->get_pos() + dis;
            auto threshold = _dsc->get_avg_edge_length();
            // Align boundary
            if (is_bound_vertex[nid.key()] == 1){
                auto direct = direction_state[nid.key()];
                if (direct | X_direction){
                    destination[0] = (std::abs(destination[0]) < threshold)? 0 : domain_dim[0];
                }
                if (direct | Y_direction){
                    destination[1] = (std::abs(destination[1]) < threshold)? 0 : domain_dim[1];
                }
                if (direct | Z_direction){
                    destination[2] = (std::abs(destination[2]) < threshold)? 0 : domain_dim[2];
                }
            }
            
            _dsc->set_destination(nid.key(), nid->get_pos() + dis);
            
            if (max_displacement_real < dis.length())
            {
                max_displacement_real = dis.length();
            }
        }
    }
    
    cout << "Max displacement: " << max_displacement_real << endl;
}
void segment_function::compute_external_force()
{
    auto c = _mean_intensities;
    
    // Array to store temporary forces
    // Use fixxed array for better performance. Suppose that we have less than 10000 vertices
    std::vector<vec3> forces = std::vector<vec3>(30000, vec3(0.0));
    
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
    auto c = _mean_intensities;
    for(auto fid = _dsc->faces_begin(); fid != _dsc->faces_end(); fid++)
    {
        if (fid->is_interface())
        {
            // check if it is stable first
            auto nodes = _dsc->get_nodes(fid.key());
            if (!(_vertex_stability_map[nodes[0]] == 1
                and _vertex_stability_map[nodes[1]] == 1
              and _vertex_stability_map[nodes[2]] == 1))
            {
                continue;
            }
            
            double var = 0.0;
            // Normal
            auto tets = _dsc->get_tets(fid.key());
            auto verts = _dsc->get_nodes(fid.key());
            auto pts = _dsc->get_pos(verts);
            
            double c0 = c[_dsc->get_label(tets[0])];
            double c1 = c[_dsc->get_label(tets[1])];
            
            vec3 Norm = _dsc->get_normal(fid.key());
            auto l01 = _dsc->barycenter(tets[1]) - _dsc->barycenter(tets[0]);
            Norm = Norm*dot(Norm, l01);// modify normal direction
            Norm.normalize();
            
            // Discretize the face
            double area = Util::area<double>(pts[0], pts[1], pts[2]);
            
            size_t n = std::ceil( sqrt(area) );
            if (n >= tri_coord_size.size())
            {
                n = tri_coord_size.size() - 1;
            }
            
            auto a = tri_dis_coord[n - 1];
            for (auto coord : a)
            {
                auto p = get_coord_tri(pts, coord);
                auto g = _img.get_value_f(p);
                
                //                auto f = - Norm* ((c1 - c0)*(2*g - c0 - c1) / area);
                auto f = - Norm* ((2*g - c0 - c1) / (c1-c0) / area);
                
                // add
                var += f.length();
            }
            
            var = var / a.size();
            
            if (var > 0.1) // The face want to move
            {
                cout << "Split face " << fid.key() << ", var = " << var << endl;
                _dsc->split(fid.key());
            }
        }
    }

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
#ifdef LOG_DEBUG
    cout << "Find intersection" << endl;
#endif
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
                        auto bc = Util::barycentric_coords<double>(vec3(x+0.5, y+0.5, 0), pts[0], pts[1], pts[2]);
                        
                        if (bc[0] > -EPSILON && bc[1] > -EPSILON && bc[2] > -EPSILON)
                        { // inside
                            auto p = pts3[0]*bc[0] + pts3[1]*bc[1] + pts3[2]*bc[2];
                            
                            
                            if(phase0 != BOUND_LABEL)
                                ray_intersect[phase0][y*dim[0] + x].intersects.push_back(intersect_pt(std::floor(p[2]), !in_1));
                            if(phase1 != BOUND_LABEL)
                                ray_intersect[phase1][y*dim[0] + x].intersects.push_back(intersect_pt(std::floor(p[2]), in_1));
                        }
                    }
                    catch (std::exception e)
                    {
                    
                    }
                }
            }
        }
    }
#ifdef LOG_DEBUG
    cout << "Count intersection" << endl;
#endif
    // 3. Compute integral
    _d_rayz.clear();
    
    std::fill(_mean_intensities.begin(), _mean_intensities.end(), 0.0);
    vector<double> area(_mean_intensities.size(), 0.0);
    for (int i = 0; i < nb_phase; i++) // we dont compute the background
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
                        
                        area[i] += z2 - z1;
                        _mean_intensities[i] += _img.sum_line_z(r.x, r.y, z1, z2);
                        
                    
                    }
                }
                else{}//???
            }
        }
#ifdef LOG_DEBUG
        cout << count << "Intersected rays" << endl;
#endif
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
        cout << _mean_intensities[i] << " -- ";
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
    
    // 3. Work around to align boundary vertices
    //  including set displacement for interface vertices
    work_around_on_boundary_vertices();
    
    _dsc->deform();
    
    /**
     4. RELABEL TETRAHEDRA
     */
    if (iteration % 20 == 0)
    {
//        initialization_discrete_opt();
        relabel_tetrahedra();
    }
    
    t.done();
    profile::close();
}
