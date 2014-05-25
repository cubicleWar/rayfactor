/*
 *  RayfactorConstants.h
 *  RayFactor 0.5
 *
 *  Created by Trevor Walker on 20/02/10.
 *  Copyright 2010 University Of Sydney. All rights reserved.
 *
 */

//#include <limits>
#include <climits>

#define PI 3.1415926535897932384626433832795f       // for float use 7 sig fig, for double 16
//#define dPI 6.283185f

//#define eps 1E-2        // Intersections before this distance are ignored
//#define inf 1E30    // Could also use std::numeric_limits<double>::infinity() and include <limits>

//when working do the same for maximum doubles etc
//const double inf = std::numeric_limits<float>::infinity();

#define kDefaultRayDensity 100000

// Parsing constants used when parsing the input xml file in Scene
#define kGlobalRayDensity	"globalRayDensity"
#define kDescription		"description"
#define kValue				"value"

