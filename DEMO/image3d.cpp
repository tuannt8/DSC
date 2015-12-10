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
    for (int i = 0; i < files.size(); i++)
    {
        cimg_byte im;
        im.load(files[i].c_str());
        
        for (int j = 0; j < im.height(); j++)
            for(int i = 0; i < im.width(); i++)
            {
                *(cur++) = (double)im(i,j) / 255.0;
            }
    }
    

}

float * image3d::get_layer(const int idx)
{
    return &_voxels[idx*_dim[X]*_dim[Y]];
}

float image3d::get_value(const int & x, const int & y, const int & z) const
{
    return _voxels[index(x,y,z)];
}

