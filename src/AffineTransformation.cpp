#include "AffineTransformation.h"

namespace affineTransformation {
    
void rotate(Matrix4 &x, float angle, Vector &u)
{
    Matrix4 rm;
    u.normalise();
    float ang = (float)angle*pi/180.0f;					//Convert the angle to radians
    float c = cos(ang);
    float s = sin(ang);
    float mc = 1.0f - c;
    
    // This is clockwise rotation
    rm.m[0] = c + mc*u._x*u._x;
    rm.m[1] = mc*u._x*u._y + s*u._z;
    rm.m[2] = mc*u._x*u._z - s*u._y;
    rm.m[4] = mc*u._y*u._x - s*u._z;
    rm.m[5] = c + mc*u._y*u._y;
    rm.m[6] = mc*u._y*u._z + s*u._x;
    rm.m[8] = mc*u._z*u._x + s*u._y;
    rm.m[9] = mc*u._z*u._y - s*u._x;
    rm.m[10] = c + mc*u._z*u._z;
    
    postMultiply(x, rm);
}

void invRotate(Matrix4 &x, float angle, Vector &u)
{
    Matrix4 invRm;
    u.normalise();
    float ang = (float)angle*pi/180.0f;					//Convert the angle to radians
    float c = cos(ang);
    float s = sin(ang);
    float mc = 1.0f - c;
    
    // This is counter clockwise rotation
    invRm.m[0] = c + mc*u._x*u._x;
    invRm.m[1] = mc*u._x*u._y - s*u._z;
    invRm.m[2] = mc*u._x*u._z + s*u._y;
    invRm.m[4] = mc*u._y*u._x + s*u._z;
    invRm.m[5] = c + mc*u._y*u._y;
    invRm.m[6] = mc*u._y*u._z - s*u._x;
    invRm.m[8] = mc*u._z*u._x - s*u._y;
    invRm.m[9] = mc*u._z*u._y + s*u._x;
    invRm.m[10] = c + mc*u._z*u._z;
    
    preMultiply(x, invRm);
}

void scale(Matrix4 &x, float sx, float sy, float sz)
{
    
    Matrix4 scale;
    
    if(fabs(sx) < sEps || fabs(sy) < sEps || fabs(sz) < sEps)
    {
        cerr << "Degenerate scaling transformation!\n";
    }
    scale.m[0] = sx;
    scale.m[5] = sy;
    scale.m[10] =  sz;
	
    postMultiply(x, scale);
}

void invScale(Matrix4 &x, float sx, float sy, float sz)
{
    Matrix4 invScale;
    
    if(fabs(sx) < sEps || fabs(sy) < sEps || fabs(sz) < sEps)
    {
        cerr << "Degenerate scaling transformation!\n";
    }
    invScale.m[0] = 1.0f/sx;
    invScale.m[5] = 1.0f/sy;
    invScale.m[10] = 1.0f/sz;
    
    preMultiply(x, invScale);
}

void translate(Matrix4 &x, const Vector &d)
{
    Matrix4 tr;
    
    tr.m[3] = d._x;
    tr.m[7] = d._y;
    tr.m[11] = d._z;
    
    postMultiply(x, tr);
}

void invTranslate(Matrix4 &x, const Vector &d)
{
    Matrix4 invTr;
    
    invTr.m[3] = -d._x;
    invTr.m[7]= -d._y;
    invTr.m[11] = -d._z;
    
    preMultiply(x, invTr);
}
    
}