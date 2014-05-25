//
//  Matrices.cpp
//  RayFactor
//
//  Created by Trevor Walker on 26/08/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//

#include "Matrices.h"


void setIdentityMatrix(Matrix4 &x)
{
    x.m[0] = x.m[5] = x.m[10] = x.m[15] = 1.0f;
    x.m[1] = x.m[2] = x.m[3] = x.m[4] = 0.0f;
    x.m[6] = x.m[7] = x.m[8] = x.m[9] = 0.0f;
    x.m[11] = x.m[12] = x.m[13] = x.m[14] = 0.0f;
}


// PreMultiply matrix x by matrix a
void preMultiply(Matrix4 &x, const Matrix4 &a)
{
    float sum;
    Matrix4 tmp = x;
    
    for(int c = 0; c < 4; c++)
        for(int r = 0; r < 4; r++)
        {
            sum = 0;
            for(int k = 0; k < 4; k++)
            {
                sum += a.m[4*k+r]*tmp.m[4*c+k];
            }
            x.m[4*c+r] = sum;
        }
}

// Postmultiple x by a
void postMultiply(Matrix4 &x, const Matrix4 &a)
{
    float sum;
    Matrix4 tmp = x;
    
    for(int c = 0; c < 4; c++)
        for(int r = 0; r < 4; r++)
        {
            sum = 0;
            for(int k = 0; k < 4; k++)
            {
                sum += tmp.m[4*k+r]*a.m[4*c+k];
            }
            x.m[4*c+r] = sum;
        }
}