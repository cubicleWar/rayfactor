/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : Primitive
//
//	Description : ScnObject is the base class for any object in a given scene
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef __PRIMITIVE_H__
#define __PRIMITIVE_H__


#include <stdlib.h>
#include <cmath>
#include <omp.h>

#include "Rayx4.h"
#include "VFMatrix.h"
#include "RayfactorConstants.h"
#include "dSFMT.h"

#include "Matrix4.h"
#include "AffineTransformation.h"
#include "cycle.h"
#include "SSE_maths.h"

using namespace std;

class Primitive
{
private:
	int primitiveID;
	bool isBoundingElement;				// Whether the rays shoot outward (N) or inward (Y)
	bool willAnalysePrimitive;
protected:
	Vector scaleVector;
public:
    static int numThreads;
    __m128 iden;
	Primitive *next;
	int rayDensity;
	
	//AffineTransfrm *aS;

    Matrix4 affine;
    Matrix4 invAffine;
    
    
	Primitive( int primitiveID );
	virtual ~Primitive(){};
	int getID();
	
	void setNext( Primitive &n );
	
	void setRayDensity( int rD );
	void setIsBounding( bool bE );
	bool isBounding() { return isBoundingElement; }
	void setWillAnalyse( bool willAnalyse );
    bool willAnalyse() { return willAnalysePrimitive; };

	void rotate( float angle, Vector u );
	void scale( float sx, float sy, float sz );
	void translate( Vector d );
	void getFirstHit( Primitive* head, Rayx4& ray,__m128 &objectNo );
	
    
	void transfrmPoint(Vectorx4 &wldPoint, const Vectorx4 &genPoint);
	void xfrmRay( Rayx4& genRay, const Rayx4& r );
	void transfrmNormal(Vectorx4 &wldNormal, const Vectorx4 &genNormal);
    
    void alignRayDirection(Rayx4 &wsRay, const Rayx4 &osRay, const Vectorx4 &wsNormal);
    void getOSRayDirection(Rayx4 &osRay, const __m128& azimuthal, const __m128 &R, __m128 zenith);
    
	// Implemented by the subclasses
	//virtual float surfaceArea()=0;
	virtual void hit( const Rayx4& r, __m128& t, __m128& objectNo )=0;
	virtual void traceFactors( Primitive* head, VFMatrix* vfm )=0;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::xfrmRay
//
//	Comments : Converts a ray in world co-ordinates into a ray in generic co-ordinates using the inverse
//			   affine transformation matrix
//
//	Arguments: genRay is the ray in generic co-ordinates which is to be determined
//			   r is the ray in world co-ordinates to be converted to generic co-ordinates
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void FORCE_INLINE Primitive::xfrmRay( Rayx4& genRay, const Rayx4& r )
{
    /*
    float sx[4], sy[4], sz[4], cx[4], cy[4], cz[4], sxr[4], cxr[4];
    _mm_store_ps(sx, r.s._x);
    _mm_store_ps(sy, r.s._y);
    _mm_store_ps(sz, r.s._z);
    _mm_store_ps(cx, r.c._x);
    _mm_store_ps(cy, r.c._y);
    _mm_store_ps(cz, r.c._z);

    //cout << this->invAffine.m[4] << " " << this->invAffine.m[5] << " " << this->invAffine.m[6] << " " << this->invAffine.m[7] << endl;
   // cout << sy[0] << endl;
    for(int i=0; i <4; i++) {
        sxr[i] = invAffine.m[0]*sx[i] + invAffine.m[1]*sy[i] + invAffine.m[2]*sz[i] + invAffine.m[3];
        cxr[i] = invAffine.m[0]*cx[i] + invAffine.m[1]*cy[i] + invAffine.m[2]*cz[i];
    }
    */
    // genRay.s._x = invAffine.m[0]*r.s._x + invAffine.m[1]*r.s._y + invAffine.m[2]*r.s._z + invAffine.m[3];
    __m128 a0 = _mm_set1_ps(invAffine.m[0]);
    __m128 a1 = _mm_set1_ps(invAffine.m[1]);
    __m128 a2 = _mm_set1_ps(invAffine.m[2]);
    
    __m128 temp1 = _mm_mul_ps(r.s._x, a0);
    __m128 temp2 = _mm_mul_ps(r.s._y, a1);
    __m128 temp3 = _mm_mul_ps(r.s._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    __m128 temp4 = _mm_add_ps(temp3, _mm_set1_ps(invAffine.m[3]));
    genRay.s._x = _mm_add_ps(temp1, temp4);
    
    //genRay.c._x = invAffine.m[0]*r.c._x + invAffine.m[1]*r.c._y + invAffine.m[2]*r.c._z;
    
    temp1 = _mm_mul_ps(r.c._x, a0);
    temp2 = _mm_mul_ps(r.c._y, a1);
    temp3 = _mm_mul_ps(r.c._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    genRay.c._x = _mm_add_ps(temp1, temp3);
    
    //genRay.s._y = invAffine.m[4]*r.s._x + invAffine.m[5]*r.s._y + invAffine.m[6]*r.s._z + invAffine.m[7];
    a0 = _mm_set1_ps(invAffine.m[4]);
    a1 = _mm_set1_ps(invAffine.m[5]);
    a2 = _mm_set1_ps(invAffine.m[6]);
    //cout << invAffine.m[4] << " " << invAffine.m[5] << " " << invAffine.m[6] << " " << invAffine.m[7] << endl;
    temp1 = _mm_mul_ps(r.s._x, a0);
    temp2 = _mm_mul_ps(r.s._y, a1);
    temp3 = _mm_mul_ps(r.s._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    temp4 = _mm_add_ps(temp3, _mm_set1_ps(invAffine.m[7]));
    genRay.s._y = _mm_add_ps(temp1, temp4);
    //printf("r.s._y    = %vf\n", r.s._y);
    //printf("genRay.s._y    = %vf\n", genRay.s._y);
    
    //genRay.c._y = invAffine.m[4]*r.c._x + invAffine.m[5]*r.c._y + invAffine.m[6]*r.c._z;
    temp1 = _mm_mul_ps(r.c._x, a0);
    temp2 = _mm_mul_ps(r.c._y, a1);
    temp3 = _mm_mul_ps(r.c._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    genRay.c._y = _mm_add_ps(temp1, temp3);
    
    //genRay.s._z = invAffine.m[8]*r.s._x + invAffine.m[9]*r.s._y + invAffine.m[10]*r.s._z + invAffine.m[11];
    a0 = _mm_set1_ps(invAffine.m[8]);
    a1 = _mm_set1_ps(invAffine.m[9]);
    a2 = _mm_set1_ps(invAffine.m[10]);
    
    temp1 = _mm_mul_ps(r.s._x, a0);
    temp2 = _mm_mul_ps(r.s._y, a1);
    temp3 = _mm_mul_ps(r.s._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    temp4 = _mm_add_ps(temp3, _mm_set1_ps(invAffine.m[11]));
    genRay.s._z = _mm_add_ps(temp1, temp4);
    
    //genRay.c._z = invAffine.m[8]*r.c._x + invAffine.m[9]*r.c._y + invAffine.m[10]*r.c._z;
    temp1 = _mm_mul_ps(r.c._x, a0);
    temp2 = _mm_mul_ps(r.c._y, a1);
    temp3 = _mm_mul_ps(r.c._z, a2);
    
    temp1 = _mm_add_ps(temp1, temp2);
    genRay.c._z = _mm_add_ps(temp1, temp3);
    
    /*
     _mm_store_ps(sx, genRay.s._x);
     _mm_store_ps(cx, genRay.c._x);
    
    for(int k = 0; k < 4; k++) {
        if(fabs(sxr[k] - sx[k]) > 0.001) {
            cout << "s wrong " << sxr[k] << " " << sx[k] << endl;
        }
        if(fabs(cxr[k] - cx[k]) > 0.001) {
            cout << "c wrong" << endl;
        }
        
    }
    */
    
}


inline void Primitive::alignRayDirection(Rayx4 &wsRay, const Rayx4 &osRay, const Vectorx4 &wsNormal) {
    //const float d = wsNormal._z == -1.0 ? 1.0f/(1.0f-wsNormal._z) : 1.0f/(1.0f+wsNormal._z);
    /*
    __m128 t0 = _mm_add_ps(*(__m128*)_ps_1, wsNormal._z);           // t0 = 1.0 - wZ
    __m128 t1 = _mm_sub_ps(*(__m128*)_ps_1, wsNormal._z);           // t1 = 1.0 + wZ
    __m128 t2 = _mm_cmpneq_ps(t1, *(__m128*)_ps_0);                 // test for t1 != 0.0, i.e. wZ != -1.0
    
    __m128 d = _mm_blendv_ps(t0, t1, t2);
    d = _mm_div_ps(*(__m128*)_ps_1, d);
    */
    //replace above the line
    
    __m128 t0 = _mm_add_ps(*(__m128*)_ps_1, wsNormal._z);           // t1 = 1.0 + wZ
    __m128 d = _mm_div_ps(*(__m128*)_ps_1, t0);
    d = _mm_blendv_ps(d, *(__m128*)_ps_0p5, _mm_cmpeq_ps(d, *(__m128*)_ps_inf));
    
    
    
    //wsRay.c._x = wsNormal._z*osRay.c._x + wsNormal._y*d*(wsNormal._y*osRay.c._x - wsNormal._x*osRay.c._y) + wsNormal._x*osRay.c._z; 
    t0 = _mm_mul_ps(wsNormal._z, osRay.c._x);
    __m128 t1 = _mm_mul_ps(wsNormal._y,d);
    __m128 t2 = _mm_sub_ps(_mm_mul_ps(wsNormal._y, osRay.c._x), _mm_mul_ps(wsNormal._x, osRay.c._y));
    __m128 t3 = _mm_mul_ps(t1, t2);
    t0 = _mm_add_ps(t0, t3);
    t2 = _mm_mul_ps(wsNormal._x, osRay.c._z);
    wsRay.c._x = _mm_add_ps(t0, t2);
    
    //wsRay.c._y = wsNormal._x*d*(wsNormal._x*osRay.c._y - wsNormal._y*osRay.c._x) + wsNormal._z*osRay.c._y + wsNormal._y*osRay.c._z;
    t0 = _mm_mul_ps(wsNormal._x, d);
    t1 = _mm_sub_ps(_mm_mul_ps(wsNormal._x, osRay.c._y), _mm_mul_ps(wsNormal._y, osRay.c._x));
    t2 = _mm_mul_ps(t0, t1);
    t0 = _mm_mul_ps(wsNormal._z, osRay.c._y);
    t0 = _mm_add_ps(t0, t2);
    t3 = _mm_mul_ps(wsNormal._y, osRay.c._z);
    wsRay.c._y = _mm_add_ps(t0, t3);
    
    //wsRay.c._z = -wsNormal._x*osRay.c._x - wsNormal._y*osRay.c._y + wsNormal._z*osRay.c._z;
    t0 = _mm_mul_ps(wsNormal._x, osRay.c._x);
    t1 = _mm_mul_ps(wsNormal._y, osRay.c._y);
    t2 = _mm_mul_ps(wsNormal._z, osRay.c._z);
    t3 = _mm_sub_ps(t2, t1);
    wsRay.c._z = _mm_sub_ps(t3, t0);
    
    /*
    float nx[4], ny[4], nz[4], cx[4], cy[4], cz[4], wsR[4], dSSE[4];
    
    _mm_store_ps(nx, wsNormal._x);
    _mm_store_ps(ny, wsNormal._y);
    _mm_store_ps(nz, wsNormal._z);
    _mm_store_ps(cx, osRay.c._x);
    _mm_store_ps(cy, osRay.c._y);
    _mm_store_ps(cz, osRay.c._z);
    _mm_store_ps(wsR, wsRay.c._z);
    _mm_store_ps(dSSE, d);

    for(int i = 0; i < 4; i++) {
        float d = nz[i] == -1.0 ? 1.0/(1.0-nz[i]) : 1.0/(1.0+nz[i]);
        
        //float res = nz[i]*cx[i] + ny[i]*d*(ny[i]*cx[i] - nx[i]*cy[i]) +nx[i]*cz[i]; // x coordinate
        float res = -nx[i]*cx[i] - ny[i]*cy[i] + nz[i]*cz[i]; // z coordinate
        
        if(fabs(d - dSSE[i]) > 0.001) {
            cout << "Error d : " << d << " " << dSSE[i] << endl;
        }
        
        if(fabs(res - wsR[i]) > 0.001) {
            cout << "Error : " << res << " " << wsR[i] << endl;
        }
        
    }*/
    
}


// R is random number, azimuthal is the az
inline void Primitive::getOSRayDirection(Rayx4 &osRay, const __m128& azimuthal, const __m128 &R, __m128 zenith) {
    zenith = _mm_sqrtApprox_ps(R);
    _mm_sincos_ps(azimuthal, &osRay.c._y, &osRay.c._x);
    
    osRay.c._x = _mm_mul_ps(zenith,osRay.c._x); 
    osRay.c._y = _mm_mul_ps(zenith, osRay.c._y);
    osRay.c._z = _mm_sqrtApprox_ps(_mm_sub_ps(*(__m128*)_ps_1, R));
}

inline void Primitive::transfrmNormal(Vectorx4 &wldNormal, const Vectorx4 &genNormal)
{
    // wldNormal._x = invAffine.m[0]*genNormal._x + invAffine.m[4]*genNormal._y + invAffine.m[8]*genNormal._z;
    __m128 t1 = _mm_mul_ps(_mm_set1_ps(invAffine.m[0]), genNormal._x);
    __m128 t2 = _mm_mul_ps(_mm_set1_ps(invAffine.m[4]), genNormal._y);
    __m128 t3 = _mm_add_ps(t1, t2);
    __m128 t4 = _mm_mul_ps(_mm_set1_ps(invAffine.m[8]), genNormal._z);
    wldNormal._x = _mm_add_ps(t3, t4);
    
    
    //wldNormal._y = invAffine.m[1]*genNormal._x + invAffine.m[5]*genNormal._y + invAffine.m[9]*genNormal._z;
    t1 = _mm_mul_ps(_mm_set1_ps(invAffine.m[1]), genNormal._x);
    t2 = _mm_mul_ps(_mm_set1_ps(invAffine.m[5]), genNormal._y);
    t3 = _mm_add_ps(t1, t2);
    t4 = _mm_mul_ps(_mm_set1_ps(invAffine.m[9]), genNormal._z);
    wldNormal._y = _mm_add_ps(t3, t4);
    
    //wldNormal._z = invAffine.m[2]*genNormal._x + invAffine.m[6]*genNormal._y + invAffine.m[10]*genNormal._z;
    t1 = _mm_mul_ps(_mm_set1_ps(invAffine.m[2]), genNormal._x);
    t2 = _mm_mul_ps(_mm_set1_ps(invAffine.m[6]), genNormal._y);
    t3 = _mm_add_ps(t1, t2);
    t4 = _mm_mul_ps(_mm_set1_ps(invAffine.m[10]), genNormal._z);
    wldNormal._z = _mm_add_ps(t3, t4);
    
    wldNormal = wldNormal.normalise();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Primitive::transfrmPoint
//
//	Comments : Converts a Point from object space to world space using this objects affine
//			   transformation matrix
//
//	Arguments: wldPoint is a transformed point determined in this function
//			   genPoint is a point in object space
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	05/05/11	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void Primitive::transfrmPoint(Vectorx4 &wldPoint, const Vectorx4 &genPoint)
{
    //  wldPoint._x = affine.m[0]*genPoint._x + affine.m[1]*genPoint._y + affine.m[2]*genPoint._z + affine.m[3];
    __m128 t1 = _mm_mul_ps(_mm_set1_ps(affine.m[0]), genPoint._x);
    __m128 t2 = _mm_mul_ps(_mm_set1_ps(affine.m[1]), genPoint._y);
    __m128 t3 = _mm_add_ps(t1, t2);
    __m128 t4 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(affine.m[2]), genPoint._z), _mm_set1_ps(affine.m[3]));
    wldPoint._x = _mm_add_ps(t3, t4);
    
    //wldPoint._y = affine.m[4]*genPoint._x + affine.m[5]*genPoint._y + affine.m[6]*genPoint._z + affine.m[7];
    t1 = _mm_mul_ps(_mm_set1_ps(affine.m[4]), genPoint._x);
    t2 = _mm_mul_ps(_mm_set1_ps(affine.m[5]), genPoint._y);
    t3 = _mm_add_ps(t1, t2);
    t4 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(affine.m[6]), genPoint._z), _mm_set1_ps(affine.m[7]));
    wldPoint._y = _mm_add_ps(t3, t4);
    
    //wldPoint._z = affine.m[8]*genPoint._x + affine.m[9]*genPoint._y + affine.m[10]*genPoint._z + affine.m[11];
    t1 = _mm_mul_ps(_mm_set1_ps(affine.m[8]), genPoint._x);
    t2 = _mm_mul_ps(_mm_set1_ps(affine.m[9]), genPoint._y);
    t3 = _mm_add_ps(t1, t2);
    t4 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(affine.m[10]), genPoint._z), _mm_set1_ps(affine.m[11]));
    wldPoint._z = _mm_add_ps(t3, t4);
}

#endif
