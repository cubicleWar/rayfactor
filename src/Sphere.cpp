#include "Sphere.h"

Sphere::Sphere(int obID): Primitive(obID)
{
}

void Sphere::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
{
    Rayx4 genRay;
    xfrmRay(genRay, r);

    const __m128 A = genRay.c.dot(genRay.c);
    const __m128 B = genRay.c.dot(genRay.s);
    const __m128 C = _mm_sub_ps(genRay.s.dot(genRay.s), *(__m128*)_ps_1);
    __m128 discrim = _mm_sub_ps(_mm_mul_ps(B,B), _mm_mul_ps(A,C));
    
    __m128 mask = _mm_cmpgt_ps(discrim, *(__m128*)_ps_0);
    
    if(_mm_movemask_ps(mask) > 0) {
        discrim = _mm_sqrtApprox_ps(discrim);
        // Could change these to use rcp with 
        const __m128 rAi = _mm_rcp_ps(A);
        const __m128 rA = _mm_sub_ps(_mm_add_ps(rAi,rAi), _mm_mul_ps(_mm_mul_ps(rAi, rAi), A));
        
        __m128 t1 = _mm_mul_ps(_mm_sub_ps(_mm_negate_ps(B),discrim),rA);
        __m128 t2 = _mm_mul_ps(_mm_sub_ps(discrim,B),rA);

        t1 = _mm_blendv_ps(*(__m128*)_ps_inf,t1,_mm_cmpgt_ps(t1, *(__m128*)_ps_EPS));
        t2 = _mm_blendv_ps(*(__m128*)_ps_inf,t2,_mm_cmpgt_ps(t2, *(__m128*)_ps_EPS));
        t1 = _mm_blendv_ps(t2, t1,_mm_cmplt_ps(t1, t2));
        
        // Another check to make sure its less than tbest
        mask = _mm_cmplt_ps(t1, tbest);
        tbest = _mm_blendv_ps(tbest, t1,mask);
        objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
    }
}  

void Sphere::traceFactors(Primitive *head, VFMatrix* vfm)
{
    const int noObjects = vfm->getNoObjects();
    
	unsigned long int noHits[noObjects];
	unsigned long int prv_noHits[noObjects];
	// Initiate the noHits array to prevent data corruption
	for(int z = 0; z < noObjects; z++) {
		noHits[z] = 0;
        prv_noHits[z] = 0;
	}	
    
    const unsigned long int noRaysCalc = ceil(rayDensity*surfaceArea()); 
    const unsigned long int noRays = noRaysCalc > 0 ? noRaysCalc : ULONG_MAX;
    
    unsigned long int i = 0;
	
	const float isBounding = this->isBounding() ? -1.0f : 1.0f;

    //#pragma omp parallel num_threads(1)
    if(Primitive::numThreads > 0) {
        omp_set_num_threads(Primitive::numThreads);
    }
    #pragma omp parallel firstprivate(i, prv_noHits)
	{
        Rayx4 wsRay, osRay;
        Vectorx4 osNormal, wsNormal;
        dsfmt_t dsfmtt;
        
        int o0, o1, o2, o3, k;
        __m128 objectNo, rand, azimuthal, temp, zenith;
        
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
            
            osRay.s._z = _mm_sub_ps(_mm_mul_ps(*(__m128*)_ps_2,rand),*(__m128*)_ps_1);
            //printf("z    = %vf\n", osRay.s._z);
            
            temp = _mm_mul_ps(osRay.s._z,osRay.s._z);
            temp = _mm_sub_ps(*(__m128*)_ps_1, temp);
            temp = _mm_sqrtApprox_ps(temp);
            
            
            azimuthal = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(0,3,2,1));
            azimuthal = _mm_mul_ps(*(__m128*)_ps_2PI, azimuthal);
            
            _mm_sincos_ps(azimuthal, &osRay.s._y, &osRay.s._x);
            
            osRay.s._x = _mm_mul_ps(temp,osRay.s._x);
            osRay.s._y = _mm_mul_ps(temp,osRay.s._y);
            
            transfrmPoint(wsRay.s, osRay.s);
            /*
            if(i < 2000) {
                float x[4], y[4], z[4];
                _mm_store_ps(x, wsRay.s._x);
                _mm_store_ps(y, wsRay.s._y);
                _mm_store_ps(z, wsRay.s._z);
            
                cout << x[0] << " " << y[0] << " " << z[0] << endl;
            }
            */
            // Find the ray direction as hemisphere lying in x-y plane
			azimuthal = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(1,0,3,2));
            azimuthal = _mm_mul_ps(*(__m128*)_ps_2PI, azimuthal);
            

            
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(2,1,0,3));
            
            getOSRayDirection(osRay, azimuthal, temp, zenith);
            
            
            osNormal._x = osRay.s._x;
			osNormal._y = osRay.s._y;
            osNormal._z = osRay.s._z;
            
            osNormal *= isBounding;
            
            
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
        #pragma omp critical    // Need to work from here
		{
			for(k = 0; k < noObjects; k++) {
				noHits[k] += prv_noHits[k];
			}
		} // End of critical section
	}
    
	// Populate the view factor matrix
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
    /*        
		#pragma omp for
		for(i = 0; i < noRays; i++) {
			// Find the start Point

            do {
                x1 = 2.0*dsfmt_genrand_open_open(&dsfmtt)-1.0;
                x2 = 2.0*dsfmt_genrand_open_open(&dsfmtt)-1.0;
                temp = x1*x1 + x2*x2;
            } while(temp > 1.0);
            
            
            osRay.s._x = 2.0*x1*sqrt(1.0-temp);
            osRay.s._y = 2.0*x2*sqrt(1.0-temp);
            osRay.s._z = 1.0-2.0*temp;

            //cout << sqrt(osRay.s._x*osRay.s._x + osRay.s._y*osRay.s._y +osRay.s._z*osRay.s._z) << endl;
            
			transfrmPoint(wsRay.s, osRay.s);
            
            // Find the ray direction as hemisphere lying in x-y plane
			azimuthal = dPI*dsfmt_genrand_close_open(&dsfmtt);
            temp = dsfmt_genrand_close_open(&dsfmtt);
            
            getOSRayDirection(osRay, azimuthal, temp, zenith);
            
            
            // Calculate the normal at the point
			osNormal._x = osRay.s._x;
			osNormal._y = osRay.s._y;
			osNormal._z = osRay.s._z;
			
            osNormal *= isBounding;
            
			// Normal is transformed and normalised.
			transfrmNormal(wsNormal, osNormal);
            alignRayDirection(wsRay, osRay, wsNormal);
			
			if(this->getFirstHit(head, wsRay, objectNo))
			{
				prv_noHits[objectNo]++;
			}
		} // End of For loop
		
		#pragma omp critical
		{
			for(k = 0; k < noObjects; k++) {
				noHits[k] += prv_noHits[k];
			}
		} // End of critical section						 
	} // End of parallel Section
	
	
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
*/

