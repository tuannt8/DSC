//
//  eigen_wrapper.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 05/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "eigen_wrapper.hpp"

#include <Eigen/SVD>

using namespace Eigen;
void svd_solve(double * C, double *Q, double *L)
{
    Eigen::Matrix3d MC;
    MC << C[0] , C[1] , C[2] , C[3] , C[4] , C[5] , C[6] , C[7] , C[8];
    Eigen::JacobiSVD<Matrix3d> svd(MC, ComputeFullU|ComputeFullV);
    
    auto U = svd.matrixU();
    auto V = svd.matrixV();
    auto S = svd.singularValues();
    
    Map<Matrix3d>(Q, 3, 3) = V;
    L[0] = S.coeffRef(0);
    L[4] = S.coeffRef(1);
    L[8] = S.coeffRef(2);
}
