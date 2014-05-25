//
//  Matrices.h
//  RayFactor
//
//  Created by Trevor Walker on 26/08/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//
#include "Matrix4.h"

#ifndef RayFactor_Matrices_h
#define RayFactor_Matrices_h

void setIdentityMatrix(Matrix4 &x);
void preMultiply(Matrix4 &x, const Matrix4 &a);
void postMultiply(Matrix4 &x, const Matrix4 &a);


#endif
