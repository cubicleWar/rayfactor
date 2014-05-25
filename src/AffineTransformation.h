#include "Matrices.h"
#include "Matrix4.h"
#include "Vector.h"

#ifndef RayFactor_AffineTransformation_h
#define RayFactor_AffineTransformation_h

namespace affineTransformation {

    static const float sEps = 0.00001f;
    static const float pi = 3.141592653589793f;
        
    void rotate(Matrix4 &x, float angle, Vector &u);
    void invRotate(Matrix4 &x, float angle, Vector &u);
    void scale(Matrix4 &x, float sx, float sy, float sz);
    void invScale(Matrix4 &x, float sx, float sy, float sz);
    void translate(Matrix4 &x, const Vector &d);
    void invTranslate(Matrix4 &x, const Vector &d);

}
#endif
