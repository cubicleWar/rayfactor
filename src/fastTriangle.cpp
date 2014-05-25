#include "fastTriangle.h"
/*
 void Triangle::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
 {
 const __m128 nx4 = _mm_set1_ps(nx);
 const __m128 ny4 = _mm_set1_ps(ny);
 const __m128 nz4 = _mm_set1_ps(nz);
 
 
 const __m128 det = _mm_add_ps(_mm_add_ps(_mm_mul_ps(r.c._x, nx4), _mm_mul_ps(r.c._y, ny4)), _mm_mul_ps(r.c._z, nz4));
 
 //const __m128 dett = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_set1_ps(nd),_mm_mul_ps(r.s._x, nx4)), 
 //                                          _mm_mul_ps(r.s._y, ny4)), _mm_mul_ps(r.s._z, nz4));
 
 const __m128 dett = _mm_sub_ps(_mm_set1_ps(nd),_mm_add_ps(_mm_add_ps(_mm_mul_ps(r.s._x, nx4), _mm_mul_ps(r.s._y, ny4)),_mm_mul_ps(r.s._z, nz4)));
 
 __m128 mask = _mm_xor_ps(dett,_mm_sub_ps(_mm_mul_ps(tbest, det),dett));
 
 if(_mm_movemask_ps(mask) != 15) {
 __m128i maski = _mm_castps_si128(mask);                     // Convert it to an integer (No instruction)
 maski = _mm_srai_epi32(maski, 31);                          // Arithmetic shift the sign bit into all the bits
 maski = _mm_xor_si128(maski, *(__m128i *)_ps_inv);          // invert all the bits
 mask = _mm_castsi128_ps(maski);                             // Convert back to float (No instruction)
 
 const __m128 detpx = _mm_add_ps(_mm_mul_ps(det, r.s._x),_mm_mul_ps(dett, r.c._x));
 const __m128 detpy = _mm_add_ps(_mm_mul_ps(det, r.s._y),_mm_mul_ps(dett, r.c._y));
 const __m128 detpz = _mm_add_ps(_mm_mul_ps(det, r.s._z),_mm_mul_ps(dett, r.c._z));
 //detpw is just det
 
 
 const __m128 detu = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(detpx, _mm_set1_ps(ux)), _mm_mul_ps(detpy, _mm_set1_ps(uy))), 
 _mm_mul_ps(detpz, _mm_set1_ps(uz))), _mm_mul_ps(det, _mm_set1_ps(ud)));
 
 __m128 mask2 = _mm_xor_ps(detu,_mm_sub_ps(det,detu));
 
 if(_mm_movemask_ps(mask2) != 15) {
 maski = _mm_castps_si128(mask2);                     // Convert it to an integer (No instruction)
 maski = _mm_srai_epi32(maski, 31);                          // Arithmetic shift the sign bit into all the bits
 maski = _mm_xor_si128(maski, *(__m128i *)_ps_inv);          // invert all the bits
 mask = _mm_and_ps(mask, _mm_castsi128_ps(maski));
 
 
 const __m128 detv = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(detpx, _mm_set1_ps(vx)), _mm_mul_ps(detpy, _mm_set1_ps(vy))), 
 _mm_mul_ps(detpz, _mm_set1_ps(vz))), _mm_mul_ps(det, _mm_set1_ps(vd)));
 
 mask2 = _mm_xor_ps(detv,_mm_sub_ps(det,_mm_add_ps(detu, detv)));
 
 if(_mm_movemask_ps(mask2) != 15) {
 
 maski = _mm_castps_si128(mask2);                     // Convert it to an integer (No instruction)
 maski = _mm_srai_epi32(maski, 31);                          // Arithmetic shift the sign bit into all the bits
 maski = _mm_xor_si128(maski, *(__m128i *)_ps_inv);          // invert all the bits
 mask = _mm_and_ps(mask, _mm_castsi128_ps(maski));
 
 // Calculate 1/det with one networn iteration
 const __m128 rdeti = _mm_rcp_ps(det);
 const __m128 rdet = _mm_sub_ps(_mm_add_ps(rdeti,rdeti), _mm_mul_ps(_mm_mul_ps(rdeti, rdeti), det));
 
 //const __m128 u = _mm_mul_ps(rdet, detu);
 //const __m128 v = _mm_mul_ps(rdet, detv);
 const __m128 t = _mm_mul_ps(rdet, dett);
 
 mask = _mm_and_ps(mask, _mm_cmpgt_ps(t, *(__m128*)_ps_EPS));
 
 mask = _mm_and_ps(mask, _mm_cmplt_ps(t, tbest));
 
 tbest = _mm_blendv_ps(tbest,t, mask);
 objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
 }
 }
 }    
 }
 */

void fastTriangle::hit(const Rayx4& r, __m128& tbest, __m128& objectNo)
{
    
    const __m128 nx4 = _mm_set1_ps(nx);
    const __m128 ny4 = _mm_set1_ps(ny);
    const __m128 nz4 = _mm_set1_ps(nz);
    
    
    const __m128 det = _mm_add_ps(_mm_add_ps(_mm_mul_ps(r.c._x, nx4), _mm_mul_ps(r.c._y, ny4)), _mm_mul_ps(r.c._z, nz4));
    
    const __m128 dett = _mm_sub_ps(_mm_set1_ps(nd),_mm_add_ps(_mm_add_ps(_mm_mul_ps(r.s._x, nx4), _mm_mul_ps(r.s._y, ny4)),_mm_mul_ps(r.s._z, nz4)));
    
    const __m128 detpx = _mm_add_ps(_mm_mul_ps(det, r.s._x),_mm_mul_ps(dett, r.c._x));
    const __m128 detpy = _mm_add_ps(_mm_mul_ps(det, r.s._y),_mm_mul_ps(dett, r.c._y));
    const __m128 detpz = _mm_add_ps(_mm_mul_ps(det, r.s._z),_mm_mul_ps(dett, r.c._z));
    //detpw is just det
    
    
    const __m128 detu = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(detpx, _mm_set1_ps(ux)), _mm_mul_ps(detpy, _mm_set1_ps(uy))), 
                                              _mm_mul_ps(detpz, _mm_set1_ps(uz))), _mm_mul_ps(det, _mm_set1_ps(ud)));
    
    const __m128 detv = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(detpx, _mm_set1_ps(vx)), _mm_mul_ps(detpy, _mm_set1_ps(vy))), 
                                              _mm_mul_ps(detpz, _mm_set1_ps(vz))), _mm_mul_ps(det, _mm_set1_ps(vd)));
    
    
    // Calculate 1/det with one networn iteration
    const __m128 rdeti = _mm_rcp_ps(det);
    const __m128 rdet = _mm_sub_ps(_mm_add_ps(rdeti,rdeti), _mm_mul_ps(_mm_mul_ps(rdeti, rdeti), det));
    
    const __m128 u = _mm_mul_ps(rdet, detu);    // was commented
    const __m128 v = _mm_mul_ps(rdet, detv);    // was commented
    const __m128 t = _mm_mul_ps(rdet, dett);
    
    const __m128 uv = _mm_add_ps(u,v);
    
    __m128 mask = _mm_and_ps(_mm_cmplt_ps(t, tbest), _mm_cmpgt_ps(t, *(__m128*)_ps_EPS));
    
    mask = _mm_and_ps(mask, _mm_cmpge_ps(u, *(__m128*)_ps_0));
    
    mask = _mm_and_ps(mask, _mm_cmpge_ps(v, *(__m128*)_ps_0));
    
    mask = _mm_and_ps(mask, _mm_cmple_ps(uv, *(__m128*)_ps_1));
    
    tbest = _mm_blendv_ps(tbest,t, mask);
    objectNo = _mm_blendv_ps(objectNo, this->iden, mask);
    
}

void fastTriangle::traceFactors(Primitive *head, VFMatrix* vfm)
{
 	const int noObjects = vfm->getNoObjects();
    
	unsigned long int noHits[noObjects];
	unsigned long int prv_noHits[noObjects];
	// Initiate the noHits array to prevent data corruption
	for(int z = 0; z < noObjects; z++) {
		noHits[z] = 0;
        prv_noHits[z] = 0;
	}	
	const unsigned long int noRays = rayDensity;
	
	unsigned long int i = 0;
	
    const Vectorx4 wsNormal(_mm_set1_ps(nx*cw), _mm_set1_ps(ny*cw), _mm_set1_ps(nz*cw));
    
	//#pragma omp parallel num_threads(1)
    if(Primitive::numThreads > 0) {
        omp_set_num_threads(Primitive::numThreads);
    }
    #pragma omp parallel firstprivate(i, prv_noHits)
	{
        Rayx4 wsRay, osRay;
        dsfmt_t dsfmtt;
        
        int o0, o1, o2, o3, k;
        __m128 objectNo, rand, azimuthal, temp, zenith;
        
		// Seed RNG;
        unsigned int seeds[5];
        time_t seconds = time(NULL);
        clock_t clocks = clock();
        seeds[0] = seconds>>   32;
        seeds[1] = (unsigned int)seconds;
        seeds[2] = clocks>>   32;
        seeds[3] = (unsigned int)clocks;
        seeds[4] = omp_get_thread_num();
        dsfmt_init_by_array(&dsfmtt, seeds, 5);
        
        wsRay.s._x = _mm_set1_ps(cx);
        wsRay.s._y = _mm_set1_ps(cy);
        wsRay.s._z = _mm_set1_ps(cz);
        
        
        #pragma omp for	
		for(i = 0; i < noRays; i+=4) {
			// Find the start Point
			
			rand = _mm_set_ps(dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt),
                              dsfmt_genrand_close_open(&dsfmtt));
            
			// Find the ray direction as hemisphere lying in x-y plane
			//azimuthal = rand;
            azimuthal = _mm_mul_ps(*(__m128*)_ps_2PI, rand);
            
            temp = _mm_shuffle_ps(rand, rand, _MM_SHUFFLE(0,3,2,1));
            
            getOSRayDirection(osRay, azimuthal, temp, zenith);
            
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
