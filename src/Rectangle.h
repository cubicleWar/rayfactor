
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : rectangle
//
//	Description : Square represents a square, plane like object. It lies in the x-y plane, has a 
//				  side length of 2 and is centered on the origin. Rays will be fired from this object in
//				  the direction of the z-ais by default.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __RECTANGLE_H__
#define __RECTANGLE_H__

#include "Primitive.h"

using namespace std;

class Rectangle : public Primitive
{
    public:
    Rectangle(int obID): Primitive(obID) {};
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);
};

inline float Rectangle::surfaceArea()
{
	return (2.0f*scaleVector._x*2.0f*scaleVector._y);
}
#endif
