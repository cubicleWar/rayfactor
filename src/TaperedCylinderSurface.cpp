#include "TaperedCylinderSurface.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	TaperedCylinderSurface::TaperedCylinderSurface
//
//	Comments : Default constructor
//
//	Arguments : obID is the objects identification number
//			    sR is the radius of the top of the cylinder (default is 1.0)
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
TaperedCylinderSurface::TaperedCylinderSurface( int obID, float sR ): Primitive(obID)
{
    smallRadius = _mm_set1_ps(sR);
}

void TaperedCylinderSurface::setSmallRadius(float sR)
{
	smallRadius = _mm_set1_ps(sR);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	TaperedCylinderSurface::hit
//
//	Comments : Determines if the CylinderSurface object was hit and records the hit data in an 
//			   intersection object.
//
//	Arguments: r is the ray to test if it intersected with this object
//			   inter is an intersection object to record the hit data
//
//	Date		Developer		Ver		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	0.1		Created
//	05/03/10	Trevor Walker	0.5		Optimised and cut down to only record the single best intersection
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void TaperedCylinderSurface::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
{
    Rayx4 genRay;
    xfrmRay(genRay, r);
    
    const __m128 sm = _mm_sub_ps(smallRadius, *(__m128*)_ps_1);
    const __m128 fDir = _mm_mul_ps(sm, genRay.c._z);
    const __m128 fStart = _mm_add_ps(_mm_mul_ps(sm, genRay.s._z),*(__m128*)_ps_1);
    
    __m128 A = _mm_add_ps(_mm_mul_ps(genRay.c._x,genRay.c._x),_mm_mul_ps(genRay.c._y,genRay.c._y));
    A = _mm_sub_ps(A,_mm_mul_ps(fDir,fDir));
    
    __m128 B = _mm_add_ps(_mm_mul_ps(genRay.s._x,genRay.c._x),_mm_mul_ps(genRay.s._y,genRay.c._y));
    B = _mm_sub_ps(B,_mm_mul_ps(fStart,fDir));
    
    __m128 C = _mm_add_ps(_mm_mul_ps(genRay.s._x,genRay.s._x),_mm_mul_ps(genRay.s._y,genRay.s._y));
    C = _mm_sub_ps(C, _mm_mul_ps(fStart, fStart));
    
    __m128 discrim = _mm_sub_ps(_mm_mul_ps(B,B),_mm_mul_ps(A,C));
    
    __m128 mask = _mm_cmpgt_ps(discrim, *(__m128*)_ps_0);
    
    if(_mm_movemask_ps(mask) > 0) {
    
        discrim = _mm_sqrtApprox_ps(discrim);
        B = _mm_negate_ps(B);
        
        
        // float t1 = (-B - discrim)/A;
        __m128 t1 = _mm_div_ps(_mm_sub_ps(B, discrim),A);
        
        //float t2 = (-B + discrim)/A;
        __m128 t2 = _mm_div_ps(_mm_add_ps(B, discrim),A);
        
        //float zHit1 = genRay.s._z + genRay.c._z*t1;
        const __m128 zHit1 = _mm_add_ps(genRay.s._z,_mm_mul_ps(genRay.c._z, t1));
        
        //float zHit2 = genRay.s._z + genRay.c._z*t2;
        const __m128 zHit2 = _mm_add_ps(genRay.s._z,_mm_mul_ps(genRay.c._z, t2));
        
        
        // Check zHit1 >= 0.0 and zHit1 <= 1.0
        __m128 mask = _mm_and_ps(_mm_cmpge_ps(zHit1,*(__m128*)_ps_0), _mm_cmple_ps(zHit1,*(__m128*)_ps_1));
        mask = _mm_and_ps(mask, _mm_cmpgt_ps(t1, *(__m128*)_ps_EPS));
        
        t1 = _mm_blendv_ps(*(__m128*)_ps_inf,t1,mask);
        
        // Check zHit2 >= 0.0 and zHit1 <= 1.0 and t > 0
        mask =  _mm_and_ps(_mm_cmpge_ps(zHit2,*(__m128*)_ps_0), _mm_cmple_ps(zHit2,*(__m128*)_ps_1));
        mask = _mm_and_ps(mask, _mm_cmpgt_ps(t2, *(__m128*)_ps_EPS));
        
        t2 = _mm_blendv_ps(*(__m128*)_ps_inf,t2, mask);
        
        t1 = _mm_blendv_ps(t2,t1, _mm_cmplt_ps(t1, t2));
        
        // Another check for t1, t2 < tbest
        mask = _mm_cmplt_ps(t1, tbest);
        tbest = _mm_blendv_ps(tbest, t1,mask);
        objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	CylinderSurface::tracefactors
//
//	Comments : Iterates over the CylinderSurface surface using the monte carlo method to fire random 
//			   rays into the scene before developing the view factor matrix for this object based on the 
//		       intersection of those rays with the objects in the current scene.
//
//	Arguments: head is a pointer to the first object in the scene
//			   vfm is the view factor matrix for this scene
//
//	Date		Developer		Ver		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker			Created
//	21/09/08	Trevor Walker			Modified to accurately handle tapered cylinder using the trap method
//	07/11/08	Trevor Walker			Modified to accurately calculate the surface area of the scale eliptical
//										tapered cylinder.
//	16/01/10	Trevor Walker	0.4		Modified to use OpenMP parallel segments
//	05/03/10	Trevor Walker	0.5		Optimised to reduce computational time
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void TaperedCylinderSurface::traceFactors(Primitive *head, VFMatrix* vfm)
{
	const int noObjects = vfm->getNoObjects();
    
	unsigned long int noHits[noObjects];
	unsigned long int prv_noHits[noObjects];
	// Initiate the noHits array to prevent data corruption
	for(int z = 0; z < noObjects; z++) {
		noHits[z] = 0;
        prv_noHits[z] = 0;
	}	
	
	//Calculate the number of rays
	const unsigned long int noRaysCalc = ceil(rayDensity*surfaceArea()); // Maybe need error checking to make sure this is smaller than maximum
    
	const unsigned long int noRays = noRaysCalc > 0.0f ? noRaysCalc : ULONG_MAX;
	
	unsigned long int i = 0;
	
	const float isBounding = this->isBounding() ? -1.0f : 1.0f;
    
    if(Primitive::numThreads > 0) {
        omp_set_num_threads(Primitive::numThreads);
    }
    #pragma omp parallel firstprivate(i, prv_noHits)
	{
        Rayx4 wsRay, osRay;
        Vectorx4 osNormal, wsNormal;
        dsfmt_t dsfmtt;
        
        int o0, o1, o2, o3, k;
        __m128 objectNo, rand, azimuthal, temp, zenith, r;
        
		//dsfmt_init_gen_rand(&dsfmtt, time(NULL)^omp_get_thread_num());
        unsigned int seeds[5];
        time_t seconds = time(NULL);
        clock_t clocks = clock();
        seeds[0] = seconds>>   32;
        seeds[1] = (unsigned int)seconds;
        seeds[2] = clocks>>   32;
        seeds[3] = (unsigned int)clocks;
        seeds[4] = omp_get_thread_num();
        dsfmt_init_by_array(&dsfmtt, seeds, 5);
        
        #pragma omp for	
		for(i = 0; i < noRays; i+=4) {
            rand = _mm_set_ps(dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt));
            
            temp = _mm_add_ps(_mm_sub_ps(*(__m128*)_ps_1,rand),_mm_mul_ps(_mm_mul_ps(rand,smallRadius),smallRadius));
            temp = _mm_sub_ps(_mm_sqrtApprox_ps(temp),*(__m128*)_ps_1);
            osRay.s._z = _mm_div_ps(temp,_mm_sub_ps(smallRadius,*(__m128*)_ps_1));
            
            r = _mm_add_ps(*(__m128*)_ps_1,_mm_mul_ps(osRay.s._z,_mm_sub_ps(smallRadius,*(__m128*)_ps_1)));
            
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(0,3,2,1));
            temp = _mm_mul_ps(temp, *(__m128*)_ps_2PI);
            _mm_sincos_ps(temp, &osRay.s._y, &osRay.s._x);
            
            
            osRay.s._x = _mm_mul_ps(r, osRay.s._x);
            osRay.s._y = _mm_mul_ps(r, osRay.s._y);
			
			transfrmPoint(wsRay.s, osRay.s);
			
            // Find the ray direction as hemisphere lying in x-y plane
            azimuthal = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(1,0,3,2));
            azimuthal = _mm_mul_ps(*(__m128*)_ps_2PI, azimuthal);
            
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(2,1,0,3));
            
            getOSRayDirection(osRay, azimuthal, temp, zenith);
            
			// Calculate the normal at the point
			osNormal._x = osRay.s._x;
			osNormal._y = osRay.s._y;
            temp = _mm_add_ps(*(__m128*)_ps_1,_mm_mul_ps(_mm_sub_ps(smallRadius,*(__m128*)_ps_1),osRay.s._z));
            osNormal._z = _mm_mul_ps(_mm_sub_ps(*(__m128*)_ps_1, smallRadius),temp);
            
			
            osNormal *= isBounding;

            // Normal is transformed to world space and normalised.
            
			transfrmNormal(wsNormal, osNormal);
            alignRayDirection(wsRay, osRay, wsNormal);
            
			this->getFirstHit(head, wsRay, objectNo);
            
            o0 = _mm_extract_ps(objectNo, 0);
            o1 = _mm_extract_ps(objectNo, 1);
            o2 = _mm_extract_ps(objectNo, 2);
            o3 = _mm_extract_ps(objectNo, 3);
            
            
            // Note not 100% safe as could be accessing prv_noHits[-1]
            prv_noHits[o0] += o0 > -1 ? 1 : 0;
            prv_noHits[o1] += o1 > -1 ? 1 : 0;
            prv_noHits[o2] += o2 > -1 ? 1 : 0;
            prv_noHits[o3] += o3 > -1 ? 1 : 0;

		} // End of For loop
		#pragma omp critical
		{
			for(k = 0; k < noObjects; k++) {
				noHits[k] += prv_noHits[k];
			}
		} // End of critical section
	}
	
	// Populate the view factor matrix
	for(int j = 0; j < noObjects; j++)
	{
		float vf = 0;
		if(noHits[j] != 0)
		{
			vf = (float)noHits[j]/(float)noRays;
		}
		vfm->setViewFactor(vf, this->getID(), j);
	}
}


