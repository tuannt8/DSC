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
    
    float * cur = &_voxels[0];
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

float image3d::get_value_f(vec3 pt)
{
    if (pt[0] < 0 || pt[0] > _dim[0]
        || pt[1] < 0 || pt[1] > _dim[1]
        || pt[2] < 0 || pt[2] > _dim[2])
        return 0;
    
    CGLA::Vec3i pti(floor(pt[0]), floor(pt[1]), floor(pt[2]));
    
    vec3 ptif(pti[0], pti[1], pti[2]);
    vec3 relative_coord = pt - ptif;
    
    // TODO: write interpolate function
    float c00 = _voxels[index(pti[0], pti[1], pti[2])] * (1 - relative_coord[0])
                + _voxels[index(pti[0] + 1, pti[1], pti[2])] * relative_coord[0];
    float c01 = _voxels[index(pti[0], pti[1], pti[2]+1)] * (1 - relative_coord[0])
                + _voxels[index(pti[0] + 1, pti[1], pti[2] + 1)] * relative_coord[0];
    float c10 = _voxels[index(pti[0], pti[1]+1, pti[2])] * (1 - relative_coord[0])
                + _voxels[index(pti[0] + 1, pti[1] + 1, pti[2])] * relative_coord[0];
    float c11 = _voxels[index(pti[0], pti[1] + 1, pti[2]  +1)] * (1 - relative_coord[0])
                + _voxels[index(pti[0] + 1, pti[1]+1, pti[2]+1)] * relative_coord[0];
    
    
    float c0 = c00*(1-relative_coord[1]) + c10*relative_coord[1];
    float c1 = c01*(1-relative_coord[1]) + c11*relative_coord[1];
    
    return c0*(1 - relative_coord[2]) + c1*relative_coord[2];
}

float image3d::get_tetra_intensity(float * total_inten, float * area)
{
    return 0;
}

float * image3d::get_layer(const int idx)
{
    return &_voxels[idx*_dim[X]*_dim[Y]];
}

float image3d::get_value(const int & x, const int & y, const int & z) const
{
    if(x >= 0 and y >= 0 and z >= 0
           and x < _dim[0] and y < _dim[1] and x < _dim[2])
    {
        return _voxels[index(x,y,z)];
    }

    return 1;
}

