//
//  image3d.hpp
//  DSC
//
//  Created by Tuan Nguyen Trung on 12/7/15.
//  Copyright © 2015 Asger Nyman Christiansen. All rights reserved.
//

#ifndef image3d_hpp
#define image3d_hpp

#include <stdio.h>
#include <string>
#define cimg_use_png
#define cimg_display 0
#include <CImg.h>
// #include <../Cellar/cimg/1.6.3/include/CImg.h>
#include <vector>
#include <cstdint>
#include "util.h"

typedef cimg_library::CImg<float> cimg_byte;

class image3d
{
    enum{X,Y,Z};
    
public:
    /** Constructor */
    image3d();
    ~image3d();
    
    /**
     Load images from directory.
     */
    void load(std::string path);
    
    /**
     Get - set voxel
     */
    float * get_layer(const int idx);
    float get_value(const int & x, const int & y, const int & z);
    inline int index(int x, int y, int z){
        return z*_dim[X]*_dim[Y] + y*_dim[X] + x;
    }
    
    int* dimension(){return _dim;}
    vec3 dimension_v(){return vec3(_dim[X], _dim[Y], _dim[Z]);}
private:
    // Currently hold all images.
    int _dim[3]; // x - y - z
    int _layer_size;
    std::vector<float> _voxels;
    
};

#endif /* image3d_hpp */

