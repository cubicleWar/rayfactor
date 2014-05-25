#include "Rectangle.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Rectangle::hit
//
//	Comments : Determines if the Rectangle object was hit and records the hit data in an intersection
//			   object.
//
//	Arguments: r is the ray to test if it intersected with this object
//			   inter is an intersection object to record the hit data
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Rectangle::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
{
    Rayx4 genRay;
    xfrmRay(genRay, r);

    __m128 th = _mm_div_ps(_mm_negate_ps(genRay.s._z), genRay.c._z);

    const __m128 hx = _mm_abs_ps(_mm_add_ps(genRay.s._x, _mm_mul_ps(genRay.c._x, th)));

    const __m128 hy = _mm_abs_ps(_mm_add_ps(genRay.s._y, _mm_mul_ps(genRay.c._y, th)));

    // th > EPS abs(hx) <= 1 abs(hy) <= 1
    __m128 mask = _mm_and_ps(_mm_cmpgt_ps(th, *(__m128*)_ps_EPS), _mm_cmple_ps(hx, *(__m128*)_ps_1));
    mask = _mm_and_ps(mask, _mm_cmple_ps(hy, *(__m128*)_ps_1));
    
    th = _mm_blendv_ps(*(__m128*)_ps_inf,th, mask);
    
    mask = _mm_cmplt_ps(th, tbest);
    tbest = _mm_blendv_ps(tbest, th,mask);
    objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Rectangle::tracefactors
//
//	Comments : Iterates over the Rectangle surface using the monte carlo method to fire random rays into 
//			   the scene before developing the view factor matrix for this object based on the 
//		       intersection of those rays with the objects in the current scene.
//
//	Arguments: head is a pointer to the first object in the scene
//			   vfm is the view factor matrix for this scene
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//	07/11/08	Trevor Walker	Modified to accurately calculate the surface area of the scaled Rectangle
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Rectangle::traceFactors(Primitive *head, VFMatrix* vfm)
{	
	const int noObjects = vfm->getNoObjects();
	unsigned long int noHits[noObjects];
	unsigned long int prv_noHits[noObjects];
	// Initiate the noHits array to prevent data corruption
	for(int z = 0; z < noObjects; z++) {
		noHits[z] = 0;
        prv_noHits[z] = 0;
	}
	
    // Calculate the Normal and the rotations from it
	// The normal is the scale regarless of the position on the disc
	const Vectorx4 osNormal(*(__m128*)_ps_0,*(__m128*)_ps_0,*(__m128*)_ps_1);
	Vectorx4 wsNormal;
	transfrmNormal(wsNormal, osNormal);
    
	const unsigned long int noRaysCalc = ceil(surfaceArea()*rayDensity);			//Ray density mulitplied by surface area
	const unsigned long int noRays = noRaysCalc > 0 ? noRaysCalc : ULONG_MAX;
    
    unsigned long int i;
    //#pragma omp parallel num_threads(1)
    if(Primitive::numThreads > 0) {
        omp_set_num_threads(Primitive::numThreads);
    }
    #pragma omp parallel firstprivate( prv_noHits, i)
    {
        Rayx4 wsRay, osRay;
        __m128 objectNo, azimuthal, temp, zenith, rand;
        int o0, o1, o2, o3, k;
        dsfmt_t dsfmtt;
	
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
		for(i = 0; i < noRays; i+=4)
		{
            rand = _mm_set_ps(dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt));
            
            temp = _mm_mul_ps(rand, *(__m128*)_ps_2);
            osRay.s._x = _mm_sub_ps(temp,*(__m128*)_ps_1);
            
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(0,3,2,1));
            temp = _mm_mul_ps(temp, *(__m128*)_ps_2);
            osRay.s._y = _mm_sub_ps(temp,*(__m128*)_ps_1);
 
            //Set the ray with this start point
            transfrmPoint(wsRay.s, osRay.s);
           			
            //Get the ray direction
            azimuthal = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(1,0,3,2));
            azimuthal = _mm_mul_ps(azimuthal, *(__m128*)_ps_2PI);

            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(2,1,0,3));
			
			getOSRayDirection(osRay, azimuthal, temp, zenith);
            // Rotate to align by the normal
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
		} // End for loop
		
		#pragma omp critical
		{
			for(k = 0; k < noObjects; k++) {
				noHits[k] += prv_noHits[k];
			}
		} // End of critical section						 
	} // End of parallel Section							 
								 
	//Populate the view factor matrix
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


