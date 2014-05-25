#include "Disc.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Disc::hit
//
//	Comments : Determines if the disc object was hit and records the hit data in an intersection
//			   object.
//
//	Arguments: r is the ray to test if it intersected with this object
//			   inter is an intersection object to record the hit data
//
//	Date		Developer		Ver		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker			Created
//	05/03/10	Trevor Walker	0.5		Removed redundant hitInfo information
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Disc::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
{
    Rayx4 genRay;
    xfrmRay(genRay, r);
    
    const __m128 th = _mm_div_ps(_mm_negate_ps(genRay.s._z), genRay.c._z);
    
    
    //float expr = (genRay.s._x+genRay.c._x*th)*(genRay.s._x+genRay.c._x*th)+(genRay.s._y + genRay.c._y*th)*(genRay.s._y + genRay.c._y*th);
    __m128 expr1 = _mm_add_ps(genRay.s._x, _mm_mul_ps(genRay.c._x,th));
    __m128 expr2 = _mm_add_ps(genRay.s._y, _mm_mul_ps(genRay.c._y,th));
    
    expr1 = _mm_mul_ps(expr1, expr1);
    expr2 = _mm_mul_ps(expr2, expr2);
    expr1 = _mm_add_ps(expr1,expr2);
    
    const __m128 mask = _mm_and_ps(_mm_cmple_ps(expr1, *(__m128*)_ps_1), _mm_cmpgt_ps(th,*(__m128*)_ps_EPS));
    // Another compare to check its bigger than tbest
    tbest = _mm_blendv_ps(tbest, th, mask );
    objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
}    

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Disc::tracefactors
//
//	Comments : Iterates over the disc surface using the monte carlo method to fire random rays into 
//			   the scene before developing the view factor matrix for this object based on the 
//		       intersection of those rays with the objects in the current scene.
//
//	Arguments: head is a pointer to the first object in the scene
//			   vfm is the view factor matrix for this scene
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//	07/11/08	Trevor Walker	Modified to accurately calculate the surface area of the scaled disc.
//	25/11/11	Trevor Walker	Modified to reduce starting point run time.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Disc::traceFactors(Primitive *head, VFMatrix* vfm)
{
    const Vectorx4 osNormal(*(__m128*)_ps_0,*(__m128*)_ps_0,*(__m128*)_ps_1);
    Vectorx4 wsNormal;
    transfrmNormal(wsNormal, osNormal);

    const int noObjects = vfm->getNoObjects();
    
    unsigned long int i = 0;
    unsigned long int noHits[noObjects];
	unsigned long int prv_noHits[noObjects];
    
	// Initiate the noHits array to prevent data corruption
	for(i = 0; i < noObjects; i++) {
		noHits[i] = 0;
        prv_noHits[i] = 0;
	}
    
    const unsigned long int noRaysCalc = ceil(surfaceArea()*rayDensity);
	const unsigned long int noRays = noRaysCalc > 0 ? noRaysCalc : ULONG_MAX;
    
    //#pragma omp parallel num_threads(1)
    if(Primitive::numThreads > 0) {
        omp_set_num_threads(Primitive::numThreads);
    }
    #pragma omp parallel firstprivate( prv_noHits) private(i)
	{
        Rayx4 wsRay, osRay;
        __m128 objectNo, angle, azimuthal, temp, zenith, rand;
        int o0, o1, o2, o3, k;
        dsfmt_t dsfmtt;
        
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
        
        
		osRay.s._z = *(__m128*)_ps_0;
        
        #pragma omp for	
		for(i = 0; i < noRays; i+=4) {	
			//Choose the radius based on the area fraction
			
            rand = _mm_set_ps(dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt));
            
            angle = _mm_mul_ps(rand, *(__m128*)_ps_2PI);
            
            _mm_sincos_ps(angle, &osRay.s._x, &osRay.s._y);
			
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(0,3,2,1));
            temp = _mm_sqrtApprox_ps(rand);
            osRay.s._x = _mm_mul_ps(temp, osRay.s._x);
			osRay.s._y = _mm_mul_ps(temp, osRay.s._y);
            
			//Set the ray with this start point
			transfrmPoint(wsRay.s, osRay.s);
			
			//get an angle out from the generic plane
            azimuthal = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(1,0,3,2));
            
            azimuthal = _mm_mul_ps(azimuthal, *(__m128*)_ps_2PI);

            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(2,1,0,3));
            
            getOSRayDirection(osRay, azimuthal, temp, zenith);
            
            alignRayDirection(wsRay, osRay, wsNormal);
            
			getFirstHit(head, wsRay, objectNo);
			
            o0 = _mm_extract_ps(objectNo, 0);
            o1 = _mm_extract_ps(objectNo, 1);
            o2 = _mm_extract_ps(objectNo, 2);
            o3 = _mm_extract_ps(objectNo, 3);
            
            prv_noHits[o0] += o0 > -1 ? 1 : 0;
            prv_noHits[o1] += o1 > -1 ? 1 : 0;
            prv_noHits[o2] += o2 > -1 ? 1 : 0;
            prv_noHits[o3] += o3 > -1 ? 1 : 0;
            
		}	// end of for loop
        #pragma omp critical
		{
			for(k = 0; k < noObjects; k++) {
				noHits[k] += prv_noHits[k];
			}
		} // end of critical
	} // End of parallel
	
	//Populate the view factor matrix	
	for(int j = 0; j < noObjects; j++)
	{
		float vf = 0.0f;
		if(noHits[j] != 0)
		{
			vf = (float)noHits[j]/(float)noRays;
		}
		vfm->setViewFactor(vf, this->getID(), j);
	}

}


