//
//  define.h
//  DSC_compare
//
//  Created by Tuan Nguyen Trung on 27/11/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef define_h
#define define_h

#include <CGLA/Vec2d.h>
#include <CGLA/Vec3d.h>
#include <CGLA/Vec4d.h>
#include <CGLA/Mat3x3d.h>
#include <CGLA/Mat4x4d.h>
#include <CGLA/Vec3i.h>
//#include "util.h"


typedef CGLA::Vec3d vec3;
typedef CGLA::Vec3i vec3i;
typedef CGLA::Mat3x3d mat3x3d;

typedef CGLA::Mat4x4d mat4x4d;
typedef CGLA::Vec4d vec4;

#define index_cube(x,y,z) ((z)*NX*NY + (y)*NX + (x))

//#define for_nodes(dsc) for (auto nit = dsc->nodes_begin(); nit != dsc->nodes_end(); nit++)
//#define for_faces(dsc) for (auto nit = dsc->faces_begin(); nit != dsc->faces_end(); nit++)

#endif /* define_h */
