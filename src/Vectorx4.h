/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : VECTOR
//
//	Description : VECTOR is an object which represents a three dimensional vector.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __VECTORX4_H__
#define __VECTORX4_H__

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "SSE_maths.h"

using namespace std;

class Vectorx4
{
	public:
		//constructors
		Vectorx4();
		Vectorx4(const Vectorx4& vec);
        Vectorx4(const __m128 x, const __m128 y, const __m128 z);
    
        void set(const __m128 x, const __m128 y, const __m128 z);
    
		//assignment operator
		Vectorx4& operator=(const Vectorx4& vec);
		
		//vector operations
		Vectorx4& operator+=(const Vectorx4& vec);
        Vectorx4 operator+(const Vectorx4& vec) const;
		Vectorx4& operator-=(const Vectorx4& vec);
        Vectorx4 operator-(const Vectorx4& vec) const;
		Vectorx4& operator*=(const Vectorx4& vec);
        Vectorx4 operator*(const Vectorx4& vec) const;
		
        __m128 dot(const Vectorx4& vec) const;
        Vectorx4& normalise();	
		
		//scalar operations
		Vectorx4& operator*=(const float scalar);
		Vectorx4 operator*(const float scalar) const;
		
		//member variables
		__m128 _x;
		__m128 _y;
		__m128 _z;
};


inline Vectorx4::Vectorx4() {
    _x = _mm_set1_ps(0.0f);
    _y = _mm_set1_ps(0.0f);
    _z = _mm_set1_ps(0.0f);
}
 
inline Vectorx4::Vectorx4(const Vectorx4& vec) { 
    _x = vec._x;
    _y = vec._y;
    _z = vec._z;
}

inline Vectorx4::Vectorx4(const __m128 x, const __m128 y, const __m128 z) { 
    _x = x;
    _y = y;
    _z = z;
}

inline void Vectorx4::set(const __m128 x, const __m128 y, const __m128 z) {
    _x = x;
    _y = y;
    _z = z;
}


inline Vectorx4& Vectorx4::operator=(const Vectorx4& vec)
{
	_x = vec._x;
	_y = vec._y;
	_z = vec._z;
	return *this;
}
	
inline Vectorx4& Vectorx4::operator+=(const Vectorx4& vec)
{
    _x = _mm_add_ps(_x, vec._x);
    _y = _mm_add_ps(_y, vec._y);
    _z = _mm_add_ps(_z, vec._z);
	
	return *this;
}

inline Vectorx4 Vectorx4::operator+(const Vectorx4& vec) const
{
	return Vectorx4(_mm_add_ps(_x, vec._x), _mm_add_ps(_y, vec._y), _mm_add_ps(_z, vec._z));
}

inline Vectorx4& Vectorx4::operator-=(const Vectorx4& vec)
{
	_x = _mm_sub_ps(_x, vec._x);
    _y = _mm_sub_ps(_y, vec._y);
    _z = _mm_sub_ps(_z, vec._z);
	
	return *this;
}

inline Vectorx4 Vectorx4::operator-(const Vectorx4& vec) const
{
	return Vectorx4(_mm_sub_ps(_x, vec._x), _mm_sub_ps(_y, vec._y), _mm_sub_ps(_z, vec._z));
}

//	Preforms scalar multiplication of this vector through the *= operator
inline Vectorx4& Vectorx4::operator*=(const float scalar)
{
	__m128 s = _mm_set1_ps(scalar);
    
    _x = _mm_mul_ps(_x, s);
    _y = _mm_mul_ps(_y, s);
    _z = _mm_mul_ps(_z, s);
    
	return *this;
}

inline Vectorx4 Vectorx4::operator*(const float scalar) const
{
    __m128 s = _mm_set1_ps(scalar);
    
	return Vectorx4(_mm_mul_ps(_x, s), _mm_mul_ps(_y, s), _mm_mul_ps(_z, s));
}



// Preform Cross product
inline Vectorx4& Vectorx4::operator*=(const Vectorx4& vec)
{
    _x = _mm_sub_ps(_mm_mul_ps(_y, vec._z), _mm_mul_ps(_z, vec._y));
    _y = _mm_sub_ps(_mm_mul_ps(vec._x, _z), _mm_mul_ps(_x, vec._z));
    _z = _mm_sub_ps(_mm_mul_ps(_x, vec._y), _mm_mul_ps(vec._x, _y));

	return *this;
}

inline Vectorx4 Vectorx4::operator*(const Vectorx4& vec) const
{
	return Vectorx4(_mm_sub_ps(_mm_mul_ps(_y, vec._z), _mm_mul_ps(_z, vec._y)),
                  _mm_sub_ps(_mm_mul_ps(vec._x, _z), _mm_mul_ps(_x, vec._z)), 
                  _mm_sub_ps(_mm_mul_ps(_x, vec._y), _mm_mul_ps(vec._x, _y)));
}

inline __m128 Vectorx4::dot(const Vectorx4& vec) const
{
    __m128 ret = _mm_add_ps(_mm_mul_ps(_x, vec._x), _mm_mul_ps(_y, vec._y));
    ret = _mm_add_ps(ret, _mm_mul_ps(_z, vec._z));
    
	return ret;
}

// Note a vector should never have a length of 0
inline Vectorx4& Vectorx4::normalise()
{
    // t = x^2 + y^2 + z^2
	__m128 t = _mm_add_ps(_mm_mul_ps(_x, _x), _mm_mul_ps(_y, _y));
    t = _mm_add_ps(t, _mm_mul_ps(_z, _z));
	
    // Could change this to just do the recipircol and see the difference
    __m128 tr = _mm_rSqrt1NR_ps(t);
        
    _x = _mm_mul_ps(_x, tr);
    _y = _mm_mul_ps(_y, tr);
    _z = _mm_mul_ps(_z, tr);

	return *this;
}

#endif
