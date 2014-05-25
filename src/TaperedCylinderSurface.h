/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : CylinderSurface
//
//	Description : CylinderSurface represents a open-ended cylinder of radius 1 whose base is centered 
//				  on the origin in the x-y plane. The cylinder surface extends for 1 unit along the z 
//				  axis.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __TAPEREDCYLINDERSURFACE_H__
#define __TAPEREDCYLINDERSURFACE_H__

#include "Primitive.h"

using namespace std;

class TaperedCylinderSurface : public Primitive
{
    private:
        __m128 smallRadius;
    public:
        TaperedCylinderSurface( int obID, float sR=1.0 );
		void setSmallRadius(float sR);
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);
};

inline float TaperedCylinderSurface::surfaceArea()
{
    // Note R should never be 1 here or the correction c += eps is required
	float c[4];
    _mm_store_ps(c, smallRadius);
    c[0] += epsilon; // Correction for calulating the number of rays to deal with the smallRadius = 0 case.
	
	return (sqrt((scaleVector._x*scaleVector._x+scaleVector._y*scaleVector._y)*c[0]*c[0])*
			(1+c[0])*scaleVector._z*PI)/(sqrt(2.0)*c[0]);
}

#endif
