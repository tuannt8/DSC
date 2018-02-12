//
//  eigen_wrapper.cpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 05/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#include "eigen_wrapper.hpp"

#include "Eigen/SVD"
#include "Eigen/Dense"

using namespace Eigen;
void svd_solve(double * C, double *Q, double *L)
{
    Eigen::Matrix3d MC;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
          MC(i,j) = C[i*3 +j];
        }
    }
    
    Eigen::JacobiSVD<Matrix3d> svd(MC, ComputeFullU);
    
    auto U = svd.matrixU();
    auto S = svd.singularValues();
    
    // Copy to U to Q
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Q[i*3 + j] = U(i,j);
        }
    }
//    int i = 0;
//    for (int j = 0; j < 3; j++)
//    {
//        Q[i*3 + j] = U(2,j);
//    }
//    i = 2;
//    for (int j = 0; j < 3; j++)
//    {
//        Q[i*3 + j] = U(0,j);
//    }
    
    
    L[0] = S.coeffRef(0);
    L[4] = S.coeffRef(1);
    L[8] = S.coeffRef(2);
}

// A is colume major
void invert4x4(double * A, double *A_i)
{
    Matrix4d MA = Map<Matrix4d>(A);
    for (int i = 0; i < 4; i++)
    {
        for(int j = 0; j< 4; j++)
            MA(i,j) = A[i*4+j];
    }
    
    auto MA_i = MA.inverse();
    
    Map<MatrixXd>( A_i, 4, 4) =   MA_i;
}
