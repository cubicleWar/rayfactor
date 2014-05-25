/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : CylinderSurface
//
//	Description : CylinderSurface represents a open-ended cylinder of radius 1 whose base is centered 
//				  on the origin in the x-y plane. The cylinder surface extends for 1 unit along the z 
//				  axis.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __CYLINDERSURFACE_H__
#define __CYLINDERSURFACE_H__

#include "Primitive.h"

using namespace std;

class CylinderSurface : public Primitive
{
    public:
        CylinderSurface( int obID );
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);
};

inline float CylinderSurface::surfaceArea()
{
	return 2.0*PI*pow(0.5*(pow(scaleVector._x,2)+pow(scaleVector._y, 2)),0.5);
}
#endif
