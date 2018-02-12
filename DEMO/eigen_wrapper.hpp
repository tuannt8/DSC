//
//  eigen_wrapper.hpp
//  DEMO
//
//  Created by Tuan Nguyen Trung on 05/12/2017.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef eigen_wrapper_hpp
#define eigen_wrapper_hpp

#include <stdio.h>
//#include "define.h"

// C = Q'LQ
void svd_solve(double * C, double *Q, double *L);

void invert4x4(double * A, double *A_i);

#endif /* eigen_wrapper_hpp */
