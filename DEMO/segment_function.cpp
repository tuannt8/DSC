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
        auto mean_inten_tetra = _img.get_tetra_intensity(nodes_pos);
        
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
        }
    }
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
    cout << "Computing average intensity with" << NB_PHASE << " phases " << endl;
    
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

    cout << "Find intersection" << endl;
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
    
    cout << "Count intersection" << endl;
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
        cout << count << "Intersected rays" << endl;
    }
    
//    cout << "Update phase 0" << endl;
//    double s0 = _img.sum_area(dim[0]-1, dim[1]-1, dim[2]-1);
//    double v0 = dim[0]*dim[1]*dim[2];
//    for (int i = 2; i < nb_phase; i++) // we dont compute the background
//    {
//        s0 -= _mean_intensities[i];
//        v0 -= area[i];
//    }
//    _mean_intensities[1] = s0;
//    area[1] = v0;
    
    for (int i  = 0; i < nb_phase; i++)
    {
        if(area[i]==0)
        {
            std::cout << "Zero; " << _mean_intensities[i] << std::endl;
        }
        _mean_intensities[i] /= area[i];
    }
    
    cout << "Done mean intensity" << endl;
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
    cout << "--------------- Iteration " << ++iteration << "----------------" << endl;
    
    // 1. Compute average intensity
    profile t("averaging intensity");
    
    _mean_intensities.resize(num_phases);
    update_average_intensity();
    
    
    /**
     2. RELABEL TETRAHEDRA
     */
    t.change("Relabel tetrahedra");
    static int mesh_opt_counter = 0;
    if (mesh_opt_counter++ > 20)
    {

    }
    
    // 3. Compute external force
    t.change("Compute force");
    compute_external_force();
    
    // 3.1: Work around to align boundary vertices
    //  including set displacement for interface vertices
    work_around_on_boundary_vertices();
    
    _dsc->deform();
    
    // 4. Relabel tetrahedron
    // Asumme that the average intensitise do not change much
//    update_average_intensity();
//    switch (iteration) {
//        case 0:
//        case 2:
//        case 5:
//        case 10:
//            relabel_tetrahedra();
//            break;
//            
//        default:
//            break;
//    }
    if (iteration % 5 == 0)
    {
        relabel_tetrahedra();
    }


    
    iteration++;
}
