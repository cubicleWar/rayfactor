/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : Annulus
//
//	Description : Annulus represents a circular annulus of variable inner and outer radius, centered 
//				  on the origin in the x-y plane. The drection of the z-axis is the direction rays will 
//				  be fired by default
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __ANNULUS_H__
#define __ANNULUS_H__

#include "Primitive.h"

using namespace std;

class Annulus : public Primitive
{
    private:
        float rOuter;
        float rInner;
        __m128 rOuterSq;
        __m128 rInnerSq;
    public:
        Annulus(int obID, float rO=1.0f, float rI=0.5f);
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);

};

inline float Annulus::surfaceArea()
{
	return (PI*(rOuter*rOuter - rInner*rInner));
}

#endif
