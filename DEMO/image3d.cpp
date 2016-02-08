//
//  image3d.cpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 12/7/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#include "image3d.h"
#include <stdio.h>
#include <iostream>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/filesystem.hpp>
#include "CImg.h"

using namespace std;


image3d::image3d()
{
    
}

image3d::~image3d()
{
    
}

void image3d::load(std::string path)
{
    if(!boost::filesystem::exists(path))
    {
        cout << path << " - Directory  does not exist" << endl;
        exit(0);
    }
    
    std::vector<cimg_byte> _images;
    
    vector<boost::filesystem::path> files;
    for ( boost::filesystem::directory_iterator it(path);
         it != boost::filesystem::directory_iterator(); ++it )
    {
        if (it->path().extension() == ".png")
        {
            files.push_back(it->path());
        }
    }
    
    int count = files.size();
    for (int i = 0; i < count; i++)
    {
        ostringstream name;
        name << path << "/im_" << i << ".png";
        files[i] = name.str();
    }
    

    cimg_byte im;
    im.load(files[0].c_str());
    
    _dim[Z] = (int)files.size();
    _dim[X] = im.width();
    _dim[Y] = im.height();
    _layer_size = _dim[X]*_dim[Y];
    
    _voxels.resize(_dim[X]*_dim[Y]*_dim[Z]);
    
    unsigned int idx = 0;
    for (int i = 0; i < files.size(); i++)
    {
        cimg_byte im;
        im.load(files[i].c_str());

        for (int j = 0; j < im.height(); j++)
            for(int i = 0; i < im.width(); i++)
            {
                _voxels[idx++] = (double)im(i,j) / 255.0;
            }
    }
}

double image3d::get_value_f(vec3 pt) const
{

    CGLA::Vec3i pti(floor(pt[0]), floor(pt[1]), floor(pt[2]));
    
    vec3 ptif(pti[0], pti[1], pti[2]);
    vec3 relative_coord = pt - ptif;
    
    // TODO: write interpolate function
    double c00 = get_value(pti[0], pti[1], pti[2]) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1], pti[2]) * relative_coord[0];
    double c01 = get_value(pti[0], pti[1], pti[2]+1) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1], pti[2] + 1) * relative_coord[0];
    double c10 = get_value(pti[0], pti[1]+1, pti[2]) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1] + 1, pti[2]) * relative_coord[0];
    double c11 = get_value(pti[0], pti[1] + 1, pti[2]  +1) * (1 - relative_coord[0])
                + get_value(pti[0] + 1, pti[1]+1, pti[2]+1) * relative_coord[0];
    
    
    double c0 = c00*(1-relative_coord[1]) + c10*relative_coord[1];
    double c1 = c01*(1-relative_coord[1]) + c11*relative_coord[1];
    
    double f =  c0*(1 - relative_coord[2]) + c1*relative_coord[2];
    assert(f < 1.01);
    return f;
}

double image3d::get_tetra_intensity(std::vector<vec3> tet_points, double * total_inten, double * volume)
{
    double v = Util::volume<double>(tet_points[0], tet_points[1], tet_points[2], tet_points[3]);
    
    int loop = 0;
    double v1 = v;
    while (v1 > 1.)
    {
        loop ++;
        v1 = v1/8.0;
    }
    
    *total_inten = 0;
    int deep = 0;
    get_integral_recur(tet_points, loop, total_inten, deep);
    
    *total_inten = *total_inten / pow(8, loop) * v;
    assert(*total_inten < 1000);
    
    if (volume)
    {
        *volume = v;
    }

    
    return *total_inten / v;
}

void image3d::get_integral_recur(std::vector<vec3> const & tet_points, int loops, double * total, int deep)
{
    if (deep >= loops)
    {
        auto midp = (tet_points[0] + tet_points[1] + tet_points[2] + tet_points[3]) / 4.0;
        *total += get_value_f(midp);
        assert(*total < 1000);
        return;
    }
    else{
        auto subdivisions = subdivide_tet(tet_points);
        for (auto e : subdivisions)
        {
            get_integral_recur(e, loops, total, deep + 1);
        }
    }
    
}

std::vector<std::vector<vec3>> image3d::subdivide_tet(std::vector<vec3> const & tet_points)
{
    std::vector<std::vector<vec3>> list;
    
    auto A0 = tet_points[0];
    auto A1 = tet_points[1];
    auto A2 = tet_points[2];
    auto A3 = tet_points[3];
    
    auto B0 = (A0 + A1)/2;
    auto B1 = (A0 + A3)/2;
    auto B2 = (A2 + A3)/2;
    auto B3 = (A1 + A3)/2;
    auto B4 = (A0 + A2)/2;
    auto B5 = (A1 + A2)/2;
    
    list.push_back({A0, B0, B4, B1});
    list.push_back({B1, B3, B2, A3});
    list.push_back({B4, B5, A2, B2});
    list.push_back({B0, A1, B5, B3});
    list.push_back({B4, B0, B2, B1});
    list.push_back({B0, B3, B2, B1});
    list.push_back({B0, B2, B5, B3});
    list.push_back({B0, B2, B4, B5});
    
    return list;
}

double * image3d::get_layer(const int idx)
{
    return &_voxels[idx*_dim[X]*_dim[Y]];
}

double image3d::get_value(const int x, const int y, const int z) const
{
    if(x >= 0 && y >= 0 && z >= 0
           && x < _dim[0] && y < _dim[1] && z < _dim[2])
    {
        double f = _voxels[index(x,y,z)];
        assert(f < 1.01 and f >= 0);
        return f;
    }

    return 0;
}

