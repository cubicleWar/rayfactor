
#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "Primitive.h"

using namespace std;

class Sphere : public Primitive
{
    public:
        Sphere(int obID);
		float surfaceArea();
        void hit(const Rayx4& r, __m128& tbest, __m128& objectNo);
		void traceFactors(Primitive *head, VFMatrix* vfm);

};

inline float Sphere::surfaceArea() 
{
	float p = 1.6075; // Knud Thomens's forumla (max error +-1.061%)
	return (float)(4.0f*PI*pow((pow(scaleVector._x,p)*pow(scaleVector._y,p)+pow(scaleVector._x,p)*pow(scaleVector._z,p)+pow(scaleVector._y,p)*pow(scaleVector._z,p))/3.0f,(1.0f/p)));
}

#endif
