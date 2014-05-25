/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : Scene
//
//	Description : The Scene object stores all the data/objects for the system to be analysed
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __SCENE_H__
#define __SCENE_H__

#include <typeinfo>
#include <vector>
#include <sys/time.h>

#include "RayfactorConstants.h"

#include "Ray.h"
#include "Primitive.h"
#include "VFMatrix.h"
#include "tinyxml.h"

#include "CylinderSurface.h"
#include "TaperedCylinderSurface.h"
#include "Annulus.h"
#include "Rectangle.h"
#include "Disc.h"
#include "Sphere.h"
#include "Triangle.h"
#include "fastTriangle.h"

using namespace std;

class Scene {
    private:
		Primitive *tail;
		Primitive *head;
	
		float globalRayDensity;
		float numPrimitives;
		string sceneDescription;
    public:
        static int numThreads;
    
        Scene();
	
		// Methods for reading a scene saved in an xml document
		void readScene(const char* filename);
		bool loadSettings( TiXmlElement *settings );
		bool loadGeometry( TiXmlElement *geometry );
		bool initPrimitive( TiXmlElement *xmlElement, Primitive &primitive );
        bool addObject( Primitive &so );
		float getNumPrimitives();
	
        Primitive* getObject(int obID);
        void findViewFactors();		// Will change to radiative heat transfer -> look at shade
};

#endif
