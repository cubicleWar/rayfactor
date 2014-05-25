/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : VFMatrix
//
//	Description : VFMatrix is the view factor matrix and stores all the view factors for an n object
//				  system into an n x n array.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __VFMATRIX_H__
#define __VFMATRIX_H__

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

class VFMatrix 
{
	private:
		int noObjects;				// The number of Objects in the system
		float **vfm;				// The array of system viewfactors
	public:
		VFMatrix(int nS);
		~VFMatrix();
		int getNoObjects();
		float getViewFactor(int r, int c);
		void setViewFactor(float vF, int r, int c);
		void print();
};
#endif
