//
//  SSEMaths.h
//  SSEMaths
//
//  Created by Trevor Walker on 5/11/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//

#include <emmintrin.h>     // MMX instructions
#include <mmintrin.h>      // SSE
#include <xmmintrin.h>       // SSE 2
#include <pmmintrin.h>     // SSE 3
#include <tmmintrin.h>     // SSSE 3
#include <smmintrin.h>     // SSE 4.1
#include <nmmintrin.h>     // SSE 4.2


#ifndef SSEMaths_SSEMaths_h
#define SSEMaths_SSEMaths_h


#define SIMD_ALIGNMENT 16
#define epsilon 1E-4

#ifdef _MSC_VER /* visual c++ */
#define ALIGN16_BEG __declspec( align( SIMD_ALIGNMENT ) )
#define ALIGN16_END
#define FORCE_INLINE __forceinline
#else /* gcc or icc */
#define ALIGN16_BEG
#define ALIGN16_END __attribute__ ( (aligned ( SIMD_ALIGNMENT ) ) )
#define FORCE_INLINE __attribute__((__always_inline__))
#endif // _MSC_VER


#define _PS_CONST(Name, Val)                                            \
static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

#define _PS_CONST_TYPE(Name, Type, Val)                                 \
static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

#define _PI32_CONST(Name, Val)                                            \
static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

// Constant numbers
_PS_CONST(1  , 1.0f);
//_PS_CONST(m1  , -1.0f);   // RayFactor
_PS_CONST(2  , 2.0f);   // RayFactor
_PS_CONST(0p5, 0.5f);
_PS_CONST(0, 0.0f);
_PS_CONST(EPS, epsilon);    // RayFactor
_PS_CONST(PI, 3.1415926535897932384626433832795f);   // RayFactor
_PS_CONST(2PI, 2.0f*3.1415926535897932384626433832795f);  // RayFactor

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);



/* the smallest non denormalized float number */
_PS_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS_CONST_TYPE(sign_mask, int, 0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PS_CONST_TYPE(inf, int, 0x7f800000);   // inf could be replaced with mant_mask
_PS_CONST_TYPE(inv, int, 0xffffffff);   // To invert bits



inline __m128 FORCE_INLINE _mm_negate_ps( const __m128& x )
{
    //return _mm_sub_ps( _mm_setzero_ps(), x);
    //return _mm_xor_ps(x, _mm_set1_ps(-0.f));
    return _mm_xor_ps(x, *(__m128*)_ps_sign_mask);
}

inline __m128 FORCE_INLINE _mm_abs_ps( const __m128& x )
{
    //return _mm_sub_ps( _mm_setzero_ps(), x);
    //return _mm_xor_ps(x, _mm_set1_ps(-0.f));
    //return _mm_andnot_ps(x, *(__m128*)_ps_sign_mask);
    return _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
}

/////////////////////////////////////////////////////////////////
//
//      SSE Reciprocal square root functions
//
/////////////////////////////////////////////////////////////////


// SSE SIMD reciprocal square root with one Newton-Raphson iteration
inline __m128 FORCE_INLINE _mm_rSqrt1NR_ps( const __m128 &x ) 
{
    __m128 c0 = _mm_set1_ps(3.0f);
    __m128 c1 = _mm_set1_ps(-0.5f);
    
    __m128 r0 = x;
    __m128 r1 = _mm_rsqrt_ps( x );
    
    // The general Newton-Raphson reciprocal square root recurrence:
    // Y(n+1) = (3 - X * Yn * Yn) * (Yn / 2)
    
    r0 = _mm_mul_ps( r0, r1 );
    r0 = _mm_mul_ps( r0, r1 );
    r0 = _mm_sub_ps( r0, c0 );
    r1 = _mm_mul_ps( r1, c1 );
    r0 = _mm_mul_ps( r0, r1 );
    
    return r0;
}


/////////////////////////////////////////////////////////////////
//
//      SSE Square root functions
//
/////////////////////////////////////////////////////////////////


// SSE SIMD square root using rsqrt with 1 newton iteration
inline __m128 FORCE_INLINE _mm_sqrtApprox_ps( const __m128& x )
{
    return _mm_mul_ps(x, _mm_rSqrt1NR_ps(x));
}

_PS_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS_CONST(cephes_log_p0, 7.0376836292E-2);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS_CONST(cephes_log_p2, 1.1676998740E-1);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS_CONST(cephes_log_q1, -2.12194440e-4);
_PS_CONST(cephes_log_q2, 0.693359375);


inline __m128 _mm_log_ps(__m128 x) {
    
    __m128i emm0;
    
    __m128 one = *(__m128*)_ps_1;
    
    __m128 invalid_mask = _mm_cmple_ps(x, _mm_setzero_ps());
    
    x = _mm_max_ps(x, *(__m128*)_ps_min_norm_pos);  /* cut off denormalized stuff */
    
    emm0 = _mm_srli_epi32(_mm_castps_si128(x), 23);
    
    /* keep only the fractional part */
    x = _mm_and_ps(x, *(__m128*)_ps_inv_mant_mask);
    x = _mm_or_ps(x, *(__m128*)_ps_0p5);
    
    
    emm0 = _mm_sub_epi32(emm0, *(__m128i*)_pi32_0x7f);
    __m128 e = _mm_cvtepi32_ps(emm0);
    
    
    e = _mm_add_ps(e, one);
    
    /* part2: 
     if( x < SQRTHF ) {
     e -= 1;
     x = x + x - 1.0;
     } else { x = x - 1.0; }
     */
    __m128 mask = _mm_cmplt_ps(x, *(__m128*)_ps_cephes_SQRTHF);
    __m128 tmp = _mm_and_ps(x, mask);
    x = _mm_sub_ps(x, one);
    e = _mm_sub_ps(e, _mm_and_ps(one, mask));
    x = _mm_add_ps(x, tmp);
    
    
    __m128 z = _mm_mul_ps(x,x);
    
    __m128 y = *(__m128*)_ps_cephes_log_p0;
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p1);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p2);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p3);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p4);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p5);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p6);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p7);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_log_p8);
    y = _mm_mul_ps(y, x);
    
    y = _mm_mul_ps(y, z);
    
    
    tmp = _mm_mul_ps(e, *(__m128*)_ps_cephes_log_q1);
    y = _mm_add_ps(y, tmp);
    
    
    tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
    y = _mm_sub_ps(y, tmp);
    
    tmp = _mm_mul_ps(e, *(__m128*)_ps_cephes_log_q2);
    x = _mm_add_ps(x, y);
    x = _mm_add_ps(x, tmp);
    x = _mm_or_ps(x, invalid_mask); // negative arg will be NAN
    return x;
}

_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341);
_PS_CONST(cephes_exp_C1, 0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1);

inline __m128 _mm_exp_ps(__m128 x) {
    __m128 tmp = _mm_setzero_ps(), fx;
    
    __m128i emm0;
    
    __m128 one = *(__m128*)_ps_1;
    
    x = _mm_min_ps(x, *(__m128*)_ps_exp_hi);
    x = _mm_max_ps(x, *(__m128*)_ps_exp_lo);
    
    /* express exp(x) as exp(g + n*log(2)) */
    fx = _mm_mul_ps(x, *(__m128*)_ps_cephes_LOG2EF);
    fx = _mm_add_ps(fx, *(__m128*)_ps_0p5);
    
    /* how to perform a floorf with SSE: just below */
    
    emm0 = _mm_cvttps_epi32(fx);
    tmp  = _mm_cvtepi32_ps(emm0);
    
    /* if greater, substract 1 */
    __m128 mask = _mm_cmpgt_ps(tmp, fx);    
    mask = _mm_and_ps(mask, one);
    fx = _mm_sub_ps(tmp, mask);
    
    tmp = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C1);
    __m128 z = _mm_mul_ps(fx, *(__m128*)_ps_cephes_exp_C2);
    x = _mm_sub_ps(x, tmp);
    x = _mm_sub_ps(x, z);
    
    z = _mm_mul_ps(x,x);
    
    __m128 y = *(__m128*)_ps_cephes_exp_p0;
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p1);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p2);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p3);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p4);
    y = _mm_mul_ps(y, x);
    y = _mm_add_ps(y, *(__m128*)_ps_cephes_exp_p5);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, x);
    y = _mm_add_ps(y, one);
    
    /* build 2^n */
    emm0 = _mm_cvttps_epi32(fx);
    emm0 = _mm_add_epi32(emm0, *(__m128i*)_pi32_0x7f);
    emm0 = _mm_slli_epi32(emm0, 23);
    __m128 pow2n = _mm_castsi128_ps(emm0);
    
    y = _mm_mul_ps(y, pow2n);
    return y;
}


/////////////////////////////////////////////////////////////////
//
//      SSE Trignometric functions
//
/////////////////////////////////////////////////////////////////


_PS_CONST(minus_cephes_DP1, -0.78515625f);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4f);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8f);
_PS_CONST(sincof_p0, -1.9515295891E-4f);
_PS_CONST(sincof_p1,  8.3321608736E-3f);
_PS_CONST(sincof_p2, -1.6666654611E-1f);
_PS_CONST(coscof_p0,  2.443315711809948E-005f);
_PS_CONST(coscof_p1, -1.388731625493765E-003f);
_PS_CONST(coscof_p2,  4.166664568298827E-002f);
_PS_CONST(cephes_FOPI, 1.27323954473516f); // 4 / M


/* Do four sines at once
 
 The code is the exact rewriting of the cephes sinf function.
 Precision is excellent as long as x < 8192 (I did not bother to
 take into account the special handling they have for greater values
 -- it does not return garbage for arguments over 8192, though, but
 the extra precision is missing).
 
 Note that it is such that sinf((float)M_PI) = 8.74e-8, which is the
 surprising but correct result.
 
 Performance is also surprisingly good, 1.33 times faster than the
 macos vsinf SSE2 function, and 1.5 times faster than the
 __vrs4_sinf of amd's ACML (which is only available in 64 bits). Not
 too bad for an SSE1 function (with no special tuning) !
 However the latter libraries probably have a much better handling of NaN,
 Inf, denormalized and other special arguments..
 
 Since it is based on SSE intrinsics, it has to be compiled at -O2 to
 deliver full speed. 
 
 */
inline __m128 _mm_sin_ps(__m128 x) { // any x
    __m128 xmm1, xmm2 = _mm_setzero_ps(), xmm3, sign_bit, y;
    
    __m128i emm0, emm2;
    
    sign_bit = x;
    /* take the absolute value */
    x = _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
    /* extract the sign bit (upper one) */
    sign_bit = _mm_and_ps(sign_bit, *(__m128*)_ps_sign_mask);
    
    /* scale by 4/Pi */
    y = _mm_mul_ps(x, *(__m128*)_ps_cephes_FOPI);
    
    //printf("plop:"); print4(y); 
    /* store the integer part of y in mm0  (integer part is number of complete Pi/4 cycles in x) */
    emm2 = _mm_cvttps_epi32(y);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
    y = _mm_cvtepi32_ps(emm2);
    /* get the swap sign flag */
    emm0 = _mm_and_si128(emm2, *(__m128i*)_pi32_4);
    emm0 = _mm_slli_epi32(emm0, 29);
    /* get the polynom selection mask 
     there is one polynom for 0 <= x <= Pi/4
     and another one for Pi/4<x<=Pi/2
     
     Both branches will be computed.
     */
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2);
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
    
    __m128 swap_sign_bit = _mm_castsi128_ps(emm0);
    __m128 poly_mask = _mm_castsi128_ps(emm2);
    sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);
    
    /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = *(__m128*)_ps_minus_cephes_DP1;
    xmm2 = *(__m128*)_ps_minus_cephes_DP2;
    xmm3 = *(__m128*)_ps_minus_cephes_DP3;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);
    
    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    y = *(__m128*)_ps_coscof_p0;
    __m128 z = _mm_mul_ps(x,x);
    
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p1);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p2);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    __m128 tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, *(__m128*)_ps_1);
    
    /* Evaluate the second polynom  (Pi/4 <= x <= Pi/2) */
    
    __m128 y2 = *(__m128*)_ps_sincof_p0;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p1);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p2);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);
    
    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    y2 = _mm_and_ps(xmm3, y2); //, xmm3);
    y = _mm_andnot_ps(xmm3, y);
    y = _mm_add_ps(y,y2);
    /* update the sign */
    y = _mm_xor_ps(y, sign_bit);
    
    return y;
}

/* almost the same as sin_ps */
inline __m128 _mm_cos_ps(__m128 x) { // any x
    __m128 xmm1, xmm2 = _mm_setzero_ps(), xmm3, y;
    __m128i emm0, emm2;
    
    /* take the absolute value */
    x = _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
    
    /* scale by 4/Pi */
    y = _mm_mul_ps(x, *(__m128*)_ps_cephes_FOPI);
    
    
    /* store the integer part of y in mm0 */
    emm2 = _mm_cvttps_epi32(y);
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
    y = _mm_cvtepi32_ps(emm2);
    
    emm2 = _mm_sub_epi32(emm2, *(__m128i*)_pi32_2);
    
    /* get the swap sign flag */
    emm0 = _mm_andnot_si128(emm2, *(__m128i*)_pi32_4);
    emm0 = _mm_slli_epi32(emm0, 29);
    /* get the polynom selection mask */
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2);
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
    
    __m128 sign_bit = _mm_castsi128_ps(emm0);
    __m128 poly_mask = _mm_castsi128_ps(emm2);
    
    /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = *(__m128*)_ps_minus_cephes_DP1;
    xmm2 = *(__m128*)_ps_minus_cephes_DP2;
    xmm3 = *(__m128*)_ps_minus_cephes_DP3;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);
    
    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    y = *(__m128*)_ps_coscof_p0;
    __m128 z = _mm_mul_ps(x,x);
    
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p1);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p2);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    __m128 tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, *(__m128*)_ps_1);
    
    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    
    __m128 y2 = *(__m128*)_ps_sincof_p0;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p1);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p2);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);
    
    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    y2 = _mm_and_ps(xmm3, y2); //, xmm3);
    y = _mm_andnot_ps(xmm3, y);
    y = _mm_add_ps(y,y2);
    /* update the sign */
    y = _mm_xor_ps(y, sign_bit);
    
    return y;
}

/* since sin_ps and cos_ps are almost identical, sincos_ps could replace both of them..
 it is almost as fast, and gives you a free cosine with your sine */
inline void _mm_sincos_ps(__m128 x, __m128 *s, __m128 *c) {
    __m128 xmm1, xmm2, xmm3 = _mm_setzero_ps(), sign_bit_sin, y;
    
    __m128i emm0, emm2, emm4;
    
    sign_bit_sin = x;
    /* take the absolute value */
    x = _mm_and_ps(x, *(__m128*)_ps_inv_sign_mask);
    /* extract the sign bit (upper one) */
    sign_bit_sin = _mm_and_ps(sign_bit_sin, *(__m128*)_ps_sign_mask);
    
    /* scale by 4/Pi */
    y = _mm_mul_ps(x, *(__m128*)_ps_cephes_FOPI);   // Number of times y is divisable by Pi/4
    
    /* store the integer part of y in emm2 */
    emm2 = _mm_cvttps_epi32(y);                     // Whole number of times y is divisable by Pi/4
    
    /* j=(j+1) & (~1) (see the cephes sources) */
    emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1); // Add 1 to the number of times y is divisble by Pi/4 (j+1)
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);  // The bitwise AND of (j + 1) and the complement of 1 
    y = _mm_cvtepi32_ps(emm2);                      // Convert integer multiple to a floating point
    
    emm4 = emm2;
    
    /* get the swap sign flag for the sine */
    emm0 = _mm_and_si128(emm2, *(__m128i*)_pi32_4);
    emm0 = _mm_slli_epi32(emm0, 29);
    __m128 swap_sign_bit_sin = _mm_castsi128_ps(emm0);
    
    /* get the polynom selection mask for the sine*/
    emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2); // Bitwise AND with 2
    emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128()); // Create a mask on equality
    __m128 poly_mask = _mm_castsi128_ps(emm2);
    
    /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
    xmm1 = *(__m128*)_ps_minus_cephes_DP1;
    xmm2 = *(__m128*)_ps_minus_cephes_DP2;
    xmm3 = *(__m128*)_ps_minus_cephes_DP3;
    xmm1 = _mm_mul_ps(y, xmm1);
    xmm2 = _mm_mul_ps(y, xmm2);
    xmm3 = _mm_mul_ps(y, xmm3);
    x = _mm_add_ps(x, xmm1);
    x = _mm_add_ps(x, xmm2);
    x = _mm_add_ps(x, xmm3);
    
    emm4 = _mm_sub_epi32(emm4, *(__m128i*)_pi32_2);
    emm4 = _mm_andnot_si128(emm4, *(__m128i*)_pi32_4);
    emm4 = _mm_slli_epi32(emm4, 29);
    __m128 sign_bit_cos = _mm_castsi128_ps(emm4);
    
    sign_bit_sin = _mm_xor_ps(sign_bit_sin, swap_sign_bit_sin);
    
    
    /* Evaluate the first polynom  (0 <= x <= Pi/4) */
    __m128 z = _mm_mul_ps(x,x);
    y = *(__m128*)_ps_coscof_p0;
    
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p1);
    y = _mm_mul_ps(y, z);
    y = _mm_add_ps(y, *(__m128*)_ps_coscof_p2);
    y = _mm_mul_ps(y, z);
    y = _mm_mul_ps(y, z);
    __m128 tmp = _mm_mul_ps(z, *(__m128*)_ps_0p5);
    y = _mm_sub_ps(y, tmp);
    y = _mm_add_ps(y, *(__m128*)_ps_1);
    
    /* Evaluate the second polynom  (Pi/4 <= x <= 0) */
    
    __m128 y2 = *(__m128*)_ps_sincof_p0;
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p1);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_add_ps(y2, *(__m128*)_ps_sincof_p2);
    y2 = _mm_mul_ps(y2, z);
    y2 = _mm_mul_ps(y2, x);
    y2 = _mm_add_ps(y2, x);
    
    /*
     Essentually this function has one polynomial to map the cosine on 0 < Pi/4 and one for sine on 0 < x < Pi/4
     Then it breaks up the range 0 < x < Pi into Pi/4 increments and uses the symmetry to just rotate the sinhe and cosine approximations
     
     Range 0 : 0 < x < Pi/4     sine is approximated by the sine, cosine by the cos
     Range 1 : Pi/4 < x < Pi/2  sine is approximated by the cos, cosine by the sin
     Range 2 : Pi/2 < x < 3Pi/4 sine is approximated by the cos, cosine by the sin
     Range 3 : 3Pi/4 < x < Pi   sine is approximated by the sin, cosine by the cos
    
    
    */
    /* select the correct result from the two polynoms */  
    xmm3 = poly_mask;
    __m128 ysin2 = _mm_and_ps(xmm3, y2);
    __m128 ysin1 = _mm_andnot_ps(xmm3, y);
    y2 = _mm_sub_ps(y2,ysin2);
    y = _mm_sub_ps(y, ysin1);
    
    xmm1 = _mm_add_ps(ysin1,ysin2);
    xmm2 = _mm_add_ps(y,y2);
    
    /* update the sign */
    *s = _mm_xor_ps(xmm1, sign_bit_sin);
    *c = _mm_xor_ps(xmm2, sign_bit_cos);
}


#endif
