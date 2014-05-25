//
//  Matrix4.h
//  RayFactor
//
//  Created by Trevor Walker on 26/08/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//

#ifndef RayFactor_Matrix4_h
#define RayFactor_Matrix4_h


class Matrix4 
{
public:
    float m[16] __attribute__ ((aligned( 16 )));						//The transformation matrix
    Matrix4() {
        m[0] = m[5] = m[10] = m[15] = 1.0f;
        m[1] = m[2] = m[3] = m[4] = 0.0f;
        m[6] = m[7] = m[8] = m[9] = 0.0f;
        m[11] = m[12] = m[13] = m[14] = 0.0f;
    }
};

#endif
