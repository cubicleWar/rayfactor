
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : Disc
//
//	Description : Disc represents a circular disc of radius 1 and centered on the origin in the x-y 
//				  plane. The drection of the z-axis is the direction rays will be fired by default
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __DISC_H__
#define __DISC_H__

#include "Primitive.h"

using namespace std;

class Disc : public Primitive
{
    public:
        Disc(int obID): Primitive(obID) {};
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);

};

inline float Disc::surfaceArea()
{
	return (PI*scaleVector._x*scaleVector._y);
}

#endif
