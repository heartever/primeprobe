#include <stdint.h>
#include <stdlib.h>
#include <sgx_trts.h>
#include <stdio.h>
#include <string.h>
#include "mpiheader.h"

unsigned short small_prime_numbers[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
    47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383,
    389, 397, 401, 409, 419, 421, 431, 433, 439, 443,
    449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
    587, 593, 599, 601, 607, 613, 617, 619, 631, 641,
    643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769,
    773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
    853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983,
    991, 997, 1009, 1013, 1019, 1021, 1031, 1033,
    1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
    1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
    1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213,
    1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
    1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307,
    1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399,
    1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
    1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493,
    1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559,
    1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609,
    1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667,
    1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
    1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789,
    1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871,
    1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931,
    1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997,
    1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
    2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111,
    2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161,
    2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243,
    2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297,
    2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
    2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411,
    2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473,
    2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551,
    2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633,
    2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
    2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729,
    2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791,
    2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851,
    2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917,
    2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
    3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061,
    3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137,
    3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209,
    3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271,
    3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
    3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391,
    3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467,
    3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533,
    3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583,
    3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
    3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709,
    3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779,
    3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851,
    3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917,
    3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
    4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049,
    4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111,
    4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177,
    4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243,
    4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
    4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391,
    4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457,
    4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519,
    4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597,
    4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
    4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729,
    4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799,
    4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889,
    4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951,
    4957, 4967, 4969, 4973, 4987, 4993, 4999,
    0
};

unsigned char __clz_tab[] =
{
  0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
  8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
  8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
  8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
};

int m_is_secure( const void *p )
{
	return 0;
}

mpi_limb_t
mpihelp_submul_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr,
		  mpi_size_t s1_size, mpi_limb_t s2_limb)
{
    mpi_limb_t cy_limb;
    mpi_size_t j;
    mpi_limb_t prod_high, prod_low;
    mpi_limb_t x;

    /* The loop counter and index J goes from -SIZE to -1.  This way
     * the loop becomes faster.  */
    j = -s1_size;
    res_ptr -= j;
    s1_ptr -= j;

    cy_limb = 0;
    do {
	umul_ppmm( prod_high, prod_low, s1_ptr[j], s2_limb);

	prod_low += cy_limb;
	cy_limb = (prod_low < cy_limb?1:0) + prod_high;

	x = res_ptr[j];
	prod_low = x - prod_low;
	cy_limb += prod_low > x?1:0;
	res_ptr[j] = prod_low;
    } while( ++j );

    return cy_limb;
}

G10_MPI_INLINE_DECL  mpi_limb_t
mpihelp_add_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr,
	       mpi_size_t s1_size, mpi_limb_t s2_limb)
{
    mpi_limb_t x;

    x = *s1_ptr++;
    s2_limb += x;
    *res_ptr++ = s2_limb;
    if( s2_limb < x ) { /* sum is less than the left operand: handle carry */
	while( --s1_size ) {
	    x = *s1_ptr++ + 1;	/* add carry */
	    *res_ptr++ = x;	/* and store */
	    if( x )		/* not 0 (no overflow): we can stop */
		goto leave;
	}
	return 1; /* return carry (size of s1 to small) */
    }

  leave:
    if( res_ptr != s1_ptr ) { /* not the same variable */
	mpi_size_t i;	       /* copy the rest */
	for( i=0; i < s1_size-1; i++ )
	    res_ptr[i] = s1_ptr[i];
    }
    return 0; /* no carry */
}

mpi_limb_t
mpihelp_add_n( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr,
	       mpi_ptr_t s2_ptr, mpi_size_t size)
{
    mpi_limb_t x, y, cy;
    mpi_size_t j;

    /* The loop counter and index J goes from -SIZE to -1.  This way
       the loop becomes faster.  */
    j = -size;

    /* Offset the base pointers to compensate for the negative indices. */
    s1_ptr -= j;
    s2_ptr -= j;
    res_ptr -= j;

    cy = 0;
    do {
	y = s2_ptr[j];
	x = s1_ptr[j];
	y += cy;		  /* add previous carry to one addend */
	cy = y < cy;		  /* get out carry from that addition */
	y += x; 		  /* add other addend */
	cy += y < x;		  /* get out carry from that add, combine */
	res_ptr[j] = y;
    } while( ++j );

    return cy;
}


G10_MPI_INLINE_DECL mpi_limb_t
mpihelp_add(mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size,
			       mpi_ptr_t s2_ptr, mpi_size_t s2_size)
{
    mpi_limb_t cy = 0;

    if( s2_size )
	cy = mpihelp_add_n( res_ptr, s1_ptr, s2_ptr, s2_size );

    if( s1_size - s2_size )
	cy = mpihelp_add_1( res_ptr + s2_size, s1_ptr + s2_size,
			    s1_size - s2_size, cy);
    return cy;
}

mpi_limb_t
mpihelp_sub_n( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr,
				  mpi_ptr_t s2_ptr, mpi_size_t size)
{
    mpi_limb_t x, y, cy;
    mpi_size_t j;

    /* The loop counter and index J goes from -SIZE to -1.  This way
       the loop becomes faster.  */
    j = -size;

    /* Offset the base pointers to compensate for the negative indices.  */
    s1_ptr -= j;
    s2_ptr -= j;
    res_ptr -= j;

    cy = 0;
    do {
	y = s2_ptr[j];
	x = s1_ptr[j];
	y += cy;		  /* add previous carry to subtrahend */
	cy = y < cy;		  /* get out carry from that addition */
	y = x - y;		  /* main subtract */
	cy += y > x;		  /* get out carry from the subtract, combine */
	res_ptr[j] = y;
    } while( ++j );

    return cy;
}

G10_MPI_INLINE_DECL mpi_limb_t
mpihelp_sub_1(mpi_ptr_t res_ptr,  mpi_ptr_t s1_ptr,
	      mpi_size_t s1_size, mpi_limb_t s2_limb )
{
    mpi_limb_t x;

    x = *s1_ptr++;
    s2_limb = x - s2_limb;
    *res_ptr++ = s2_limb;
    if( s2_limb > x ) {
	while( --s1_size ) {
	    x = *s1_ptr++;
	    *res_ptr++ = x - 1;
	    if( x )
		goto leave;
	}
	return 1;
    }

  leave:
    if( res_ptr != s1_ptr ) {
	mpi_size_t i;
	for( i=0; i < s1_size-1; i++ )
	    res_ptr[i] = s1_ptr[i];
    }
    return 0;
}  

G10_MPI_INLINE_DECL   mpi_limb_t
mpihelp_sub( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size,
				mpi_ptr_t s2_ptr, mpi_size_t s2_size)
{
    mpi_limb_t cy = 0;

    if( s2_size )
	cy = mpihelp_sub_n(res_ptr, s1_ptr, s2_ptr, s2_size);

    if( s1_size - s2_size )
	cy = mpihelp_sub_1(res_ptr + s2_size, s1_ptr + s2_size,
				      s1_size - s2_size, cy);
    return cy;
}

/* Multiply the natural numbers u (pointed to by UP) and v (pointed to by VP),
 * both with SIZE limbs, and store the result at PRODP.  2 * SIZE limbs are
 * always stored.  Return the most significant limb.
 *
 * Argument constraints:
 * 1. PRODP != UP and PRODP != VP, i.e. the destination
 *    must be distinct from the multiplier and the multiplicand.
 *
 *
 * Handle simple cases with traditional multiplication.
 *
 * This is the most critical code of multiplication.  All multiplies rely
 * on this, both small and huge.  Small ones arrive here immediately.  Huge
 * ones arrive here as this is the base case for Karatsuba's recursive
 * algorithm below.
 */
static mpi_limb_t
mul_n_basecase( mpi_ptr_t prodp, mpi_ptr_t up,
				 mpi_ptr_t vp, mpi_size_t size)
{
    mpi_size_t i;
    mpi_limb_t cy;
    mpi_limb_t v_limb;

    /* Multiply by the first limb in V separately, as the result can be
     * stored (not added) to PROD.  We also avoid a loop for zeroing.  */
    v_limb = vp[0];
    if( v_limb <= 1 ) {
	if( v_limb == 1 )
	    MPN_COPY( prodp, up, size );
	else
	    MPN_ZERO( prodp, size );
	cy = 0;
    }
    else
	cy = mpihelp_mul_1( prodp, up, size, v_limb );

    prodp[size] = cy;
    prodp++;

    /* For each iteration in the outer loop, multiply one limb from
     * U with one limb from V, and add it to PROD.  */
    for( i = 1; i < size; i++ ) {
	v_limb = vp[i];
	if( v_limb <= 1 ) {
	    cy = 0;
	    if( v_limb == 1 )
	       cy = mpihelp_add_n(prodp, prodp, up, size);
	}
	else
	    cy = mpihelp_addmul_1(prodp, up, size, v_limb);

	prodp[size] = cy;
	prodp++;
    }

    return cy;
}


G10_MPI_INLINE_DECL mpi_limb_t
mpihelp_add(mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size,
			       mpi_ptr_t s2_ptr, mpi_size_t s2_size);

MPI
mpi_alloc( unsigned nlimbs )
{
    MPI a;

    a = (MPI)malloc( sizeof *a );
    a->d = nlimbs? mpi_alloc_limb_space( nlimbs, 0 ) : NULL;

    a->alloced = nlimbs;
    a->nlimbs = 0;
    a->sign = 0;
    a->flags = 0;
    a->nbits = 0;
    return a;
}


/****************
 * Note: This copy function should not interpret the MPI
 *	 but copy it transparently.
 */
MPI
mpi_copy( MPI a )
{
    int i;
    MPI b;

    if( a && (a->flags & 4) ) {
	void *p = malloc( a->nbits );
	memcpy( p, a->d, a->nbits );
	b = mpi_set_opaque( NULL, p, a->nbits );
    }
    else if( a ) {

	b = mpi_alloc( a->nlimbs );
	b->nlimbs = a->nlimbs;
	b->sign = a->sign;
	b->flags  = a->flags;
	b->nbits = a->nbits;
	for(i=0; i < b->nlimbs; i++ )
	    b->d[i] = a->d[i];
    }
    else
	b = NULL;
    return b;
}

int mpi_is_secure(MPI a)
{
	return 0;
}

mpi_ptr_t
mpi_alloc_limb_space( unsigned nlimbs, int secure )
{
    size_t len = nlimbs * sizeof(mpi_limb_t);
    mpi_ptr_t p;

    //p = secure? malloc( len ):malloc( len );
    p = (mpi_ptr_t)malloc( len );

    return p;
}

void mpi_free_limb_space( mpi_ptr_t a )
{
    if( !a )
		return;
   
    free(a);
}

MPI
mpi_set_opaque( MPI a, void *p, unsigned int len )
{
    if( !a ) {
		a = mpi_alloc(0);
    }

    if( a->flags & 4 )
	free( a->d );
    else {
		mpi_free_limb_space(a->d);
    }

    a->d = (mpi_limb_t *)p;
    a->alloced = 0;
    a->nlimbs = 0;
    a->nbits = len;
    a->flags = 4;
    return a;
}

void
mpihelp_release_karatsuba_ctx( struct karatsuba_ctx *ctx )
{
    struct karatsuba_ctx *ctx2;

    if( ctx->tp )
	mpi_free_limb_space( ctx->tp );
    if( ctx->tspace )
	mpi_free_limb_space( ctx->tspace );
    for( ctx=ctx->next; ctx; ctx = ctx2 ) {
	ctx2 = ctx->next;
	if( ctx->tp )
	    mpi_free_limb_space( ctx->tp );
	if( ctx->tspace )
	    mpi_free_limb_space( ctx->tspace );
	free( ctx );
    }
}


/****************
 * Resize the array of A to NLIMBS. the additional space is cleared
 * (set to 0) [done by xrealloc()]
 */
void
mpi_resize( MPI a, unsigned nlimbs )
{
    if( nlimbs <= a->alloced )
		return; /* no need to do it */
    /* Note: a->secure is not used - instead the realloc functions
     * take care of it. Maybe we should drop a->secure completely
     * and rely on a mpi_is_secure function, which would be
     * a wrapper around m_is_secure
     */
    if( a->d )
		a->d = (mpi_limb_t *)realloc(a->d, nlimbs * sizeof(mpi_limb_t) );
    else
		a->d = (mpi_limb_t *)malloc( nlimbs * sizeof(mpi_limb_t) );

    a->alloced = nlimbs;
}

void
mpi_fdiv_r( MPI rem, MPI dividend, MPI divisor )
{
    int divisor_sign = divisor->sign;
    MPI temp_divisor = NULL;

    /* We need the original value of the divisor after the remainder has been
     * preliminary calculated.	We have to copy it to temporary space if it's
     * the same variable as REM.  */
    if( rem == divisor ) {
	temp_divisor = mpi_copy( divisor );
	divisor = temp_divisor;
    }

    mpi_tdiv_r( rem, dividend, divisor );

    if( ((divisor_sign?1:0) ^ (dividend->sign?1:0)) && rem->nlimbs )
	mpi_add( rem, rem, divisor);

    if( temp_divisor )
	mpi_free(temp_divisor);
}

unsigned long
mpi_fdiv_r_ui( MPI rem, MPI dividend, unsigned long divisor )
{
    mpi_limb_t rlimb;

    rlimb = mpihelp_mod_1( dividend->d, dividend->nlimbs, divisor );
    if( rlimb && dividend->sign )
	rlimb = divisor - rlimb;

    if( rem ) {
	rem->d[0] = rlimb;
	rem->nlimbs = rlimb? 1:0;
    }
    return rlimb;
}

void
mpi_add(MPI w, MPI u, MPI v)
{
    mpi_ptr_t wp, up, vp;
    mpi_size_t usize, vsize, wsize;
    int usign, vsign, wsign;

    if( u->nlimbs < v->nlimbs ) { /* Swap U and V. */
	usize = v->nlimbs;
	usign = v->sign;
	vsize = u->nlimbs;
	vsign = u->sign;
	wsize = usize + 1;
	RESIZE_IF_NEEDED(w, wsize);
	/* These must be after realloc (u or v may be the same as w).  */
	up    = v->d;
	vp    = u->d;
    }
    else {
	usize = u->nlimbs;
	usign = u->sign;
	vsize = v->nlimbs;
	vsign = v->sign;
	wsize = usize + 1;
	RESIZE_IF_NEEDED(w, wsize);
	/* These must be after realloc (u or v may be the same as w).  */
	up    = u->d;
	vp    = v->d;
    }
    wp = w->d;
    wsign = 0;

    if( !vsize ) {  /* simple */
	MPN_COPY(wp, up, usize );
	wsize = usize;
	wsign = usign;
    }
    else if( usign != vsign ) { /* different sign */
	/* This test is right since USIZE >= VSIZE */
	if( usize != vsize ) {
	    mpihelp_sub(wp, up, usize, vp, vsize);
	    wsize = usize;
	    MPN_NORMALIZE(wp, wsize);
	    wsign = usign;
	}
	else if( mpihelp_cmp(up, vp, usize) < 0 ) {
	    mpihelp_sub_n(wp, vp, up, usize);
	    wsize = usize;
	    MPN_NORMALIZE(wp, wsize);
	    if( !usign )
		wsign = 1;
	}
	else {
	    mpihelp_sub_n(wp, up, vp, usize);
	    wsize = usize;
	    MPN_NORMALIZE(wp, wsize);
	    if( usign )
		wsign = 1;
	}
    }
    else { /* U and V have same sign. Add them. */
	mpi_limb_t cy = mpihelp_add(wp, up, usize, vp, vsize);
	wp[usize] = cy;
	wsize = usize + cy;
	if( usign )
	    wsign = 1;
    }

    w->nlimbs = wsize;
    w->sign = wsign;
}

void
mpi_free( MPI a )
{
    if( !a )
		return;
    //MYDEBUG(88);
    //MYDEBUG(a->nlimbs);
    //MYDEBUG(a->flags);
    //MYDEBUG(33);
    if( a->flags & 4 )
    {
    	//MYDEBUG(55);
		free( a->d );
		//MYDEBUG(44);
	}
    else {
    	//MYDEBUG(77);
		mpi_free_limb_space(a->d);
		//MYDEBUG(66);
    }
	//MYDEBUG(99);
    free(a);
}

void
mpi_set( MPI w, MPI u)
{
    mpi_ptr_t wp, up;
    mpi_size_t usize = u->nlimbs;
    int usign = u->sign;

    RESIZE_IF_NEEDED(w, usize);
    wp = w->d;
    up = u->d;
    MPN_COPY( wp, up, usize );
    w->nlimbs = usize;
    w->nbits = u->nbits;
    w->flags = u->flags;
    w->sign = usign;
}

unsigned int
mpi_get_flags (MPI a)
{
  return a->flags;
}

#define mpi_is_opaque(a) ((a) && (mpi_get_flags (a)&4))
/****************
 * Sometimes we have MSL (most significant limbs) which are 0;
 * this is for some reasons not good, so this function removes them.
 */
void
mpi_normalize( MPI a )
{
    if( mpi_is_opaque (a) )
	return;

    for( ; a->nlimbs && !a->d[a->nlimbs-1]; a->nlimbs-- )
	;
}

/****************
 * Compare OP1_PTR/OP1_SIZE with OP2_PTR/OP2_SIZE.
 * There are no restrictions on the relative sizes of
 * the two arguments.
 * Return 1 if OP1 > OP2, 0 if they are equal, and -1 if OP1 < OP2.
 */
int
mpihelp_cmp( mpi_ptr_t op1_ptr, mpi_ptr_t op2_ptr, mpi_size_t size )
{
    mpi_size_t i;
    mpi_limb_t op1_word, op2_word;

    for( i = size - 1; i >= 0 ; i--) {
	op1_word = op1_ptr[i];
	op2_word = op2_ptr[i];
	if( op1_word != op2_word )
	    goto diff;
    }
    return 0;

  diff:
    /* This can *not* be simplified to
     *	 op2_word - op2_word
     * since that expression might give signed overflow.  */
    return (op1_word > op2_word) ? 1 : -1;
}


int
mpi_cmp( MPI u, MPI v )
{
    mpi_size_t usize, vsize;
    int cmp;

    mpi_normalize( u );
    mpi_normalize( v );
    usize = u->nlimbs;
    vsize = v->nlimbs;
    if( !u->sign && v->sign )
	return 1;
    if( u->sign && !v->sign )
	return -1;
    if( usize != vsize && !u->sign && !v->sign )
	return usize - vsize;
    if( usize != vsize && u->sign && v->sign )
	return vsize + usize;
    if( !usize )
	return 0;
    if( !(cmp=mpihelp_cmp( u->d, v->d, usize )) )
	return 0;
    if( (cmp < 0?1:0) == (u->sign?1:0))
	return 1;
    return -1;
}

int
mpi_cmp_ui( MPI u, unsigned long v )
{
    mpi_limb_t limb = v;

    mpi_normalize( u );
    if( !u->nlimbs && !limb )
	return 0;
    if( u->sign )
	return -1;
    if( u->nlimbs > 1 )
	return 1;

    if( u->d[0] == limb )
	return 0;
    else if( u->d[0] > limb )
	return 1;
    else
	return -1;
}



/* Return the number of limbs to store an MPI which is specified by
   the number of bytes to represent it. */
unsigned int
mpi_nlimb_hint_from_nbits (unsigned int nbits)
{
  return (nbits+BITS_PER_MPI_LIMB-1) / BITS_PER_MPI_LIMB;
}

void
mpihelp_mul_karatsuba_case( mpi_ptr_t prodp,
			    mpi_ptr_t up, mpi_size_t usize,
			    mpi_ptr_t vp, mpi_size_t vsize,
			    struct karatsuba_ctx *ctx )
{
    mpi_limb_t cy;

    if( !ctx->tspace || ctx->tspace_size < vsize ) {
	if( ctx->tspace )
	    mpi_free_limb_space( ctx->tspace );
	ctx->tspace = mpi_alloc_limb_space( 2 * vsize,
				       m_is_secure( up ) || m_is_secure( vp ) );
	ctx->tspace_size = vsize;
    }

    MPN_MUL_N_RECURSE( prodp, up, vp, vsize, ctx->tspace );

    prodp += vsize;
    up += vsize;
    usize -= vsize;
    if( usize >= vsize ) {
	if( !ctx->tp || ctx->tp_size < vsize ) {
	    if( ctx->tp )
		mpi_free_limb_space( ctx->tp );
	    ctx->tp = mpi_alloc_limb_space( 2 * vsize, m_is_secure( up )
						      || m_is_secure( vp ) );
	    ctx->tp_size = vsize;
	}

	do {
	    MPN_MUL_N_RECURSE( ctx->tp, up, vp, vsize, ctx->tspace );
	    cy = mpihelp_add_n( prodp, prodp, ctx->tp, vsize );
	    mpihelp_add_1( prodp + vsize, ctx->tp + vsize, vsize, cy );
	    prodp += vsize;
	    up += vsize;
	    usize -= vsize;
	} while( usize >= vsize );
    }

    if( usize ) {
	if( usize < KARATSUBA_THRESHOLD ) {
	    mpihelp_mul( ctx->tspace, vp, vsize, up, usize );
	}
	else {
	    if( !ctx->next ) {
		ctx->next = (karatsuba_ctx*)malloc( sizeof *ctx );
	    }
	    mpihelp_mul_karatsuba_case( ctx->tspace,
					vp, vsize,
					up, usize,
					ctx->next );
	}

	cy = mpihelp_add_n( prodp, prodp, ctx->tspace, vsize);
	mpihelp_add_1( prodp + vsize, ctx->tspace + vsize, usize, cy );
    }
}

static void
mul_n( mpi_ptr_t prodp, mpi_ptr_t up, mpi_ptr_t vp,
			mpi_size_t size, mpi_ptr_t tspace )
{
    if( size & 1 ) {
      /* The size is odd, and the code below doesn't handle that.
       * Multiply the least significant (size - 1) limbs with a recursive
       * call, and handle the most significant limb of S1 and S2
       * separately.
       * A slightly faster way to do this would be to make the Karatsuba
       * code below behave as if the size were even, and let it check for
       * odd size in the end.  I.e., in essence move this code to the end.
       * Doing so would save us a recursive call, and potentially make the
       * stack grow a lot less.
       */
      mpi_size_t esize = size - 1;	 /* even size */
      mpi_limb_t cy_limb;

      MPN_MUL_N_RECURSE( prodp, up, vp, esize, tspace );
      cy_limb = mpihelp_addmul_1( prodp + esize, up, esize, vp[esize] );
      prodp[esize + esize] = cy_limb;
      cy_limb = mpihelp_addmul_1( prodp + esize, vp, size, up[esize] );
      prodp[esize + size] = cy_limb;
    }
    else {
	/* Anatolij Alekseevich Karatsuba's divide-and-conquer algorithm.
	 *
	 * Split U in two pieces, U1 and U0, such that
	 * U = U0 + U1*(B**n),
	 * and V in V1 and V0, such that
	 * V = V0 + V1*(B**n).
	 *
	 * UV is then computed recursively using the identity
	 *
	 *	  2n   n	  n			n
	 * UV = (B  + B )U V  +  B (U -U )(V -V )  +  (B + 1)U V
	 *		  1 1	     1	0   0  1	      0 0
	 *
	 * Where B = 2**BITS_PER_MP_LIMB.
	 */
	mpi_size_t hsize = size >> 1;
	mpi_limb_t cy;
	int negflg;

	/* Product H.	   ________________  ________________
	 *		  |_____U1 x V1____||____U0 x V0_____|
	 * Put result in upper part of PROD and pass low part of TSPACE
	 * as new TSPACE.
	 */
	MPN_MUL_N_RECURSE(prodp + size, up + hsize, vp + hsize, hsize, tspace);

	/* Product M.	   ________________
	 *		  |_(U1-U0)(V0-V1)_|
	 */
	if( mpihelp_cmp(up + hsize, up, hsize) >= 0 ) {
	    mpihelp_sub_n(prodp, up + hsize, up, hsize);
	    negflg = 0;
	}
	else {
	    mpihelp_sub_n(prodp, up, up + hsize, hsize);
	    negflg = 1;
	}
	if( mpihelp_cmp(vp + hsize, vp, hsize) >= 0 ) {
	    mpihelp_sub_n(prodp + hsize, vp + hsize, vp, hsize);
	    negflg ^= 1;
	}
	else {
	    mpihelp_sub_n(prodp + hsize, vp, vp + hsize, hsize);
	    /* No change of NEGFLG.  */
	}
	/* Read temporary operands from low part of PROD.
	 * Put result in low part of TSPACE using upper part of TSPACE
	 * as new TSPACE.
	 */
	MPN_MUL_N_RECURSE(tspace, prodp, prodp + hsize, hsize, tspace + size);

	/* Add/copy product H. */
	MPN_COPY (prodp + hsize, prodp + size, hsize);
	cy = mpihelp_add_n( prodp + size, prodp + size,
			    prodp + size + hsize, hsize);

	/* Add product M (if NEGFLG M is a negative number) */
	if(negflg)
	    cy -= mpihelp_sub_n(prodp + hsize, prodp + hsize, tspace, size);
	else
	    cy += mpihelp_add_n(prodp + hsize, prodp + hsize, tspace, size);

	/* Product L.	   ________________  ________________
	 *		  |________________||____U0 x V0_____|
	 * Read temporary operands from low part of PROD.
	 * Put result in low part of TSPACE using upper part of TSPACE
	 * as new TSPACE.
	 */
	MPN_MUL_N_RECURSE(tspace, up, vp, hsize, tspace + size);

	/* Add/copy Product L (twice) */

	cy += mpihelp_add_n(prodp + hsize, prodp + hsize, tspace, size);
	if( cy )
	  mpihelp_add_1(prodp + hsize + size, prodp + hsize + size, hsize, cy);

	MPN_COPY(prodp, tspace, hsize);
	cy = mpihelp_add_n(prodp + hsize, prodp + hsize, tspace + hsize, hsize);
	if( cy )
	    mpihelp_add_1(prodp + size, prodp + size, size, 1);
    }
}



/////////////////////div
#define UDIV_NEEDS_NORMALIZATION 1

/* Define this unconditionally, so it can be used for debugging.  */
#define __udiv_qrnnd_c(q, r, n1, n0, d) \
  do {									\
    UWtype __d1, __d0, __q1, __q0, __r1, __r0, __m;			\
    __d1 = __ll_highpart (d);						\
    __d0 = __ll_lowpart (d);						\
									\
    __r1 = (n1) % __d1; 						\
    __q1 = (n1) / __d1; 						\
    __m = (UWtype) __q1 * __d0; 					\
    __r1 = __r1 * __ll_B | __ll_highpart (n0);				\
    if (__r1 < __m)							\
      { 								\
	__q1--, __r1 += (d);						\
	if (__r1 >= (d)) /* i.e. we didn't get carry when adding to __r1 */\
	  if (__r1 < __m)						\
	    __q1--, __r1 += (d);					\
      } 								\
    __r1 -= __m;							\
									\
    __r0 = __r1 % __d1; 						\
    __q0 = __r1 / __d1; 						\
    __m = (UWtype) __q0 * __d0; 					\
    __r0 = __r0 * __ll_B | __ll_lowpart (n0);				\
    if (__r0 < __m)							\
      { 								\
	__q0--, __r0 += (d);						\
	if (__r0 >= (d))						\
	  if (__r0 < __m)						\
	    __q0--, __r0 += (d);					\
      } 								\
    __r0 -= __m;							\
									\
    (q) = (UWtype) __q1 * __ll_B | __q0;				\
    (r) = __r0; 							\
  } while (0)
#define udiv_qrnnd __udiv_qrnnd_c
 
#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									\
    UWtype __x; 							\
    __x = (al) - (bl);							\
    (sh) = (ah) - (bh) - (__x > (al));					\
    (sl) = __x; 							\
  } while (0)

/* Divide the two-limb number in (NH,,NL) by D, with DI being the largest
 * limb not larger than (2**(2*BITS_PER_MP_LIMB))/D - (2**BITS_PER_MP_LIMB).
 * If this would yield overflow, DI should be the largest possible number
 * (i.e., only ones).  For correct operation, the most significant bit of D
 * has to be set.  Put the quotient in Q and the remainder in R.
 */
#define UDIV_QRNND_PREINV(q, r, nh, nl, d, di) \
    do {							    \
	mpi_limb_t _q, _ql, _r; 				    \
	mpi_limb_t _xh, _xl;					    \
	umul_ppmm (_q, _ql, (nh), (di));			    \
	_q += (nh);	/* DI is 2**BITS_PER_MPI_LIMB too small */  \
	umul_ppmm (_xh, _xl, _q, (d));				    \
	sub_ddmmss (_xh, _r, (nh), (nl), _xh, _xl);		    \
	if( _xh ) {						    \
	    sub_ddmmss (_xh, _r, _xh, _r, 0, (d));		    \
	    _q++;						    \
	    if( _xh) {						    \
		sub_ddmmss (_xh, _r, _xh, _r, 0, (d));		    \
		_q++;						    \
	    }							    \
	}							    \
	if( _r >= (d) ) {					    \
	    _r -= (d);						    \
	    _q++;						    \
	}							    \
	(r) = _r;						    \
	(q) = _q;						    \
    } while (0)

/****************
 * Divide (DIVIDEND_PTR,,DIVIDEND_SIZE) by DIVISOR_LIMB.
 * Write DIVIDEND_SIZE limbs of quotient at QUOT_PTR.
 * Return the single-limb remainder.
 * There are no constraints on the value of the divisor.
 *
 * QUOT_PTR and DIVIDEND_PTR might point to the same limb.
 */

mpi_limb_t
mpihelp_divmod_1( mpi_ptr_t quot_ptr,
		  mpi_ptr_t dividend_ptr, mpi_size_t dividend_size,
		  mpi_limb_t divisor_limb)
{
    mpi_size_t i;
    mpi_limb_t n1, n0, r;
    int dummy;

    if( !dividend_size )
	return 0;

    /* If multiplication is much faster than division, and the
     * dividend is large, pre-invert the divisor, and use
     * only multiplications in the inner loop.
     *
     * This test should be read:
     * Does it ever help to use udiv_qrnnd_preinv?
     * && Does what we save compensate for the inversion overhead?
     */
    if( UDIV_TIME > (2 * UMUL_TIME + 6)
	&& (UDIV_TIME - (2 * UMUL_TIME + 6)) * dividend_size > UDIV_TIME ) {
	int normalization_steps;

	count_leading_zeros( normalization_steps, divisor_limb );
	if( normalization_steps ) {
	    mpi_limb_t divisor_limb_inverted;

	    divisor_limb <<= normalization_steps;

	    /* Compute (2**2N - 2**N * DIVISOR_LIMB) / DIVISOR_LIMB.  The
	     * result is a (N+1)-bit approximation to 1/DIVISOR_LIMB, with the
	     * most significant bit (with weight 2**N) implicit.
	     */
	    /* Special case for DIVISOR_LIMB == 100...000.  */
	    if( !(divisor_limb << 1) )
		divisor_limb_inverted = ~(mpi_limb_t)0;
	    else
		udiv_qrnnd(divisor_limb_inverted, dummy,
			   -divisor_limb, 0, divisor_limb);

	    n1 = dividend_ptr[dividend_size - 1];
	    r = n1 >> (BITS_PER_MPI_LIMB - normalization_steps);

	    /* Possible optimization:
	     * if (r == 0
	     * && divisor_limb > ((n1 << normalization_steps)
	     *		       | (dividend_ptr[dividend_size - 2] >> ...)))
	     * ...one division less...
	     */
	    for( i = dividend_size - 2; i >= 0; i--) {
		n0 = dividend_ptr[i];
		UDIV_QRNND_PREINV( quot_ptr[i + 1], r, r,
				   ((n1 << normalization_steps)
			 | (n0 >> (BITS_PER_MPI_LIMB - normalization_steps))),
			      divisor_limb, divisor_limb_inverted);
		n1 = n0;
	    }
	    UDIV_QRNND_PREINV( quot_ptr[0], r, r,
			       n1 << normalization_steps,
			       divisor_limb, divisor_limb_inverted);
	    return r >> normalization_steps;
	}
	else {
	    mpi_limb_t divisor_limb_inverted;

	    /* Compute (2**2N - 2**N * DIVISOR_LIMB) / DIVISOR_LIMB.  The
	     * result is a (N+1)-bit approximation to 1/DIVISOR_LIMB, with the
	     * most significant bit (with weight 2**N) implicit.
	     */
	    /* Special case for DIVISOR_LIMB == 100...000.  */
	    if( !(divisor_limb << 1) )
		divisor_limb_inverted = ~(mpi_limb_t) 0;
	    else
		udiv_qrnnd(divisor_limb_inverted, dummy,
			   -divisor_limb, 0, divisor_limb);

	    i = dividend_size - 1;
	    r = dividend_ptr[i];

	    if( r >= divisor_limb )
		r = 0;
	    else
		quot_ptr[i--] = 0;

	    for( ; i >= 0; i-- ) {
		n0 = dividend_ptr[i];
		UDIV_QRNND_PREINV( quot_ptr[i], r, r,
				   n0, divisor_limb, divisor_limb_inverted);
	    }
	    return r;
	}
    }
    else {
	if(UDIV_NEEDS_NORMALIZATION) {
	    int normalization_steps;

	    count_leading_zeros (normalization_steps, divisor_limb);
	    if( normalization_steps ) {
		divisor_limb <<= normalization_steps;

		n1 = dividend_ptr[dividend_size - 1];
		r = n1 >> (BITS_PER_MPI_LIMB - normalization_steps);

		/* Possible optimization:
		 * if (r == 0
		 * && divisor_limb > ((n1 << normalization_steps)
		 *		   | (dividend_ptr[dividend_size - 2] >> ...)))
		 * ...one division less...
		 */
		for( i = dividend_size - 2; i >= 0; i--) {
		    n0 = dividend_ptr[i];
		    udiv_qrnnd (quot_ptr[i + 1], r, r,
			     ((n1 << normalization_steps)
			 | (n0 >> (BITS_PER_MPI_LIMB - normalization_steps))),
				divisor_limb);
		    n1 = n0;
		}
		udiv_qrnnd (quot_ptr[0], r, r,
			    n1 << normalization_steps,
			    divisor_limb);
		return r >> normalization_steps;
	    }
	}
	/* No normalization needed, either because udiv_qrnnd doesn't require
	 * it, or because DIVISOR_LIMB is already normalized.  */
	i = dividend_size - 1;
	r = dividend_ptr[i];

	if(r >= divisor_limb)
	    r = 0;
	else
	    quot_ptr[i--] = 0;

	for(; i >= 0; i--) {
	    n0 = dividend_ptr[i];
	    udiv_qrnnd( quot_ptr[i], r, r, n0, divisor_limb );
	}
	return r;
    }
}

void
mpi_tdiv_r( MPI rem, MPI num, MPI den)
{
    mpi_tdiv_qr(NULL, rem, num, den );
}

void
mpi_tdiv_qr( MPI quot, MPI rem, MPI num, MPI den)
{
    mpi_ptr_t np, dp;
    mpi_ptr_t qp, rp;
    mpi_size_t nsize = num->nlimbs;
    mpi_size_t dsize = den->nlimbs;
    mpi_size_t qsize, rsize;
    mpi_size_t sign_remainder = num->sign;
    mpi_size_t sign_quotient = num->sign ^ den->sign;
    unsigned normalization_steps;
    mpi_limb_t q_limb;
    mpi_ptr_t marker[5];
    int markidx=0;

    /* Ensure space is enough for quotient and remainder.
     * We need space for an extra limb in the remainder, because it's
     * up-shifted (normalized) below.  */
    rsize = nsize + 1;
    mpi_resize( rem, rsize);

    qsize = rsize - dsize;	  /* qsize cannot be bigger than this.	*/
    if( qsize <= 0 ) {
	if( num != rem ) {
	    rem->nlimbs = num->nlimbs;
	    rem->sign = num->sign;
	    MPN_COPY(rem->d, num->d, nsize);
	}
	if( quot ) {
	    /* This needs to follow the assignment to rem, in case the
	     * numerator and quotient are the same.  */
	    quot->nlimbs = 0;
	    quot->sign = 0;
	}
	return;
    }

    if( quot )
	mpi_resize( quot, qsize);

    /* Read pointers here, when reallocation is finished.  */
    np = num->d;
    dp = den->d;
    rp = rem->d;

    /* Optimize division by a single-limb divisor.  */
    if( dsize == 1 ) {
	mpi_limb_t rlimb;
	if( quot ) {
	    qp = quot->d;
	    rlimb = mpihelp_divmod_1( qp, np, nsize, dp[0] );
	    qsize -= qp[qsize - 1] == 0;
	    quot->nlimbs = qsize;
	    quot->sign = sign_quotient;
	}
	else
	    rlimb = mpihelp_mod_1( np, nsize, dp[0] );
	rp[0] = rlimb;
	rsize = rlimb != 0?1:0;
	rem->nlimbs = rsize;
	rem->sign = sign_remainder;
	return;
    }


    if( quot ) {
	qp = quot->d;
	/* Make sure QP and NP point to different objects.  Otherwise the
	 * numerator would be gradually overwritten by the quotient limbs.  */
	if(qp == np) { /* Copy NP object to temporary space.  */
	    np = marker[markidx++] = mpi_alloc_limb_space(nsize,
							  mpi_is_secure(quot));
	    MPN_COPY(np, qp, nsize);
	}
    }
    else /* Put quotient at top of remainder. */
	qp = rp + dsize;

    count_leading_zeros( normalization_steps, dp[dsize - 1] );

    /* Normalize the denominator, i.e. make its most significant bit set by
     * shifting it NORMALIZATION_STEPS bits to the left.  Also shift the
     * numerator the same number of steps (to keep the quotient the same!).
     */
    if( normalization_steps ) {
	mpi_ptr_t tp;
	mpi_limb_t nlimb;

	/* Shift up the denominator setting the most significant bit of
	 * the most significant word.  Use temporary storage not to clobber
	 * the original contents of the denominator.  */
	tp = marker[markidx++] = mpi_alloc_limb_space(dsize,mpi_is_secure(den));
	mpihelp_lshift( tp, dp, dsize, normalization_steps );
	dp = tp;

	/* Shift up the numerator, possibly introducing a new most
	 * significant word.  Move the shifted numerator in the remainder
	 * meanwhile.  */
	nlimb = mpihelp_lshift(rp, np, nsize, normalization_steps);
	if( nlimb ) {
	    rp[nsize] = nlimb;
	    rsize = nsize + 1;
	}
	else
	    rsize = nsize;
    }
    else {
	/* The denominator is already normalized, as required.	Copy it to
	 * temporary space if it overlaps with the quotient or remainder.  */
	if( dp == rp || (quot && (dp == qp))) {
	    mpi_ptr_t tp;

	    tp = marker[markidx++] = mpi_alloc_limb_space(dsize, mpi_is_secure(den));
	    MPN_COPY( tp, dp, dsize );
	    dp = tp;
	}

	/* Move the numerator to the remainder.  */
	if( rp != np )
	    MPN_COPY(rp, np, nsize);

	rsize = nsize;
    }

    q_limb = mpihelp_divrem( qp, 0, rp, rsize, dp, dsize );

    if( quot ) {
	qsize = rsize - dsize;
	if(q_limb) {
	    qp[qsize] = q_limb;
	    qsize += 1;
	}

	quot->nlimbs = qsize;
	quot->sign = sign_quotient;
    }

    rsize = dsize;
    MPN_NORMALIZE (rp, rsize);

    if( normalization_steps && rsize ) {
	mpihelp_rshift(rp, rp, rsize, normalization_steps);
	rsize -= rp[rsize - 1] == 0?1:0;
    }

    rem->nlimbs = rsize;
    rem->sign	= sign_remainder;
    while( markidx )
	mpi_free_limb_space(marker[--markidx]);
}

mpi_limb_t
mpihelp_mod_1(mpi_ptr_t dividend_ptr, mpi_size_t dividend_size,
				      mpi_limb_t divisor_limb)
{
    mpi_size_t i;
    mpi_limb_t n1, n0, r;
    int dummy;

    /* Botch: Should this be handled at all?  Rely on callers?	*/
    if( !dividend_size )
	return 0;

    /* If multiplication is much faster than division, and the
     * dividend is large, pre-invert the divisor, and use
     * only multiplications in the inner loop.
     *
     * This test should be read:
     *	 Does it ever help to use udiv_qrnnd_preinv?
     *	   && Does what we save compensate for the inversion overhead?
     */
    if( UDIV_TIME > (2 * UMUL_TIME + 6)
	&& (UDIV_TIME - (2 * UMUL_TIME + 6)) * dividend_size > UDIV_TIME ) {
	int normalization_steps;

	count_leading_zeros( normalization_steps, divisor_limb );
	if( normalization_steps ) {
	    mpi_limb_t divisor_limb_inverted;

	    divisor_limb <<= normalization_steps;

	    /* Compute (2**2N - 2**N * DIVISOR_LIMB) / DIVISOR_LIMB.  The
	     * result is a (N+1)-bit approximation to 1/DIVISOR_LIMB, with the
	     * most significant bit (with weight 2**N) implicit.
	     *
	     * Special case for DIVISOR_LIMB == 100...000.
	     */
	    if( !(divisor_limb << 1) )
		divisor_limb_inverted = ~(mpi_limb_t)0;
	    else
		udiv_qrnnd(divisor_limb_inverted, dummy,
			   -divisor_limb, 0, divisor_limb);

	    n1 = dividend_ptr[dividend_size - 1];
	    r = n1 >> (BITS_PER_MPI_LIMB - normalization_steps);

	    /* Possible optimization:
	     * if (r == 0
	     * && divisor_limb > ((n1 << normalization_steps)
	     *		       | (dividend_ptr[dividend_size - 2] >> ...)))
	     * ...one division less...
	     */
	    for( i = dividend_size - 2; i >= 0; i--) {
		n0 = dividend_ptr[i];
		UDIV_QRNND_PREINV(dummy, r, r,
				   ((n1 << normalization_steps)
			  | (n0 >> (BITS_PER_MPI_LIMB - normalization_steps))),
			  divisor_limb, divisor_limb_inverted);
		n1 = n0;
	    }
	    UDIV_QRNND_PREINV(dummy, r, r,
			      n1 << normalization_steps,
			      divisor_limb, divisor_limb_inverted);
	    return r >> normalization_steps;
	}
	else {
	    mpi_limb_t divisor_limb_inverted;

	    /* Compute (2**2N - 2**N * DIVISOR_LIMB) / DIVISOR_LIMB.  The
	     * result is a (N+1)-bit approximation to 1/DIVISOR_LIMB, with the
	     * most significant bit (with weight 2**N) implicit.
	     *
	     * Special case for DIVISOR_LIMB == 100...000.
	     */
	    if( !(divisor_limb << 1) )
		divisor_limb_inverted = ~(mpi_limb_t)0;
	    else
		udiv_qrnnd(divisor_limb_inverted, dummy,
			    -divisor_limb, 0, divisor_limb);

	    i = dividend_size - 1;
	    r = dividend_ptr[i];

	    if( r >= divisor_limb )
		r = 0;
	    else
		i--;

	    for( ; i >= 0; i--) {
		n0 = dividend_ptr[i];
		UDIV_QRNND_PREINV(dummy, r, r,
				  n0, divisor_limb, divisor_limb_inverted);
	    }
	    return r;
	}
    }
    else {
	if( UDIV_NEEDS_NORMALIZATION ) {
	    int normalization_steps;

	    count_leading_zeros(normalization_steps, divisor_limb);
	    if( normalization_steps ) {
		divisor_limb <<= normalization_steps;

		n1 = dividend_ptr[dividend_size - 1];
		r = n1 >> (BITS_PER_MPI_LIMB - normalization_steps);

		/* Possible optimization:
		 * if (r == 0
		 * && divisor_limb > ((n1 << normalization_steps)
		 *		   | (dividend_ptr[dividend_size - 2] >> ...)))
		 * ...one division less...
		 */
		for(i = dividend_size - 2; i >= 0; i--) {
		    n0 = dividend_ptr[i];
		    udiv_qrnnd (dummy, r, r,
				((n1 << normalization_steps)
			 | (n0 >> (BITS_PER_MPI_LIMB - normalization_steps))),
			 divisor_limb);
		    n1 = n0;
		}
		udiv_qrnnd (dummy, r, r,
			    n1 << normalization_steps,
			    divisor_limb);
		return r >> normalization_steps;
	    }
	}
	/* No normalization needed, either because udiv_qrnnd doesn't require
	 * it, or because DIVISOR_LIMB is already normalized.  */
	i = dividend_size - 1;
	r = dividend_ptr[i];

	if(r >= divisor_limb)
	    r = 0;
	else
	    i--;

	for(; i >= 0; i--) {
	    n0 = dividend_ptr[i];
	    udiv_qrnnd (dummy, r, r, n0, divisor_limb);
	}
	return r;
    }
}

/* Divide num (NP/NSIZE) by den (DP/DSIZE) and write
 * the NSIZE-DSIZE least significant quotient limbs at QP
 * and the DSIZE long remainder at NP.	If QEXTRA_LIMBS is
 * non-zero, generate that many fraction bits and append them after the
 * other quotient limbs.
 * Return the most significant limb of the quotient, this is always 0 or 1.
 *
 * Preconditions:
 * 0. NSIZE >= DSIZE.
 * 1. The most significant bit of the divisor must be set.
 * 2. QP must either not overlap with the input operands at all, or
 *    QP + DSIZE >= NP must hold true.	(This means that it's
 *    possible to put the quotient in the high part of NUM, right after the
 *    remainder in NUM.
 * 3. NSIZE >= DSIZE, even if QEXTRA_LIMBS is non-zero.
 */
mpi_limb_t
mpihelp_divrem( mpi_ptr_t qp, mpi_size_t qextra_limbs,
		mpi_ptr_t np, mpi_size_t nsize,
		mpi_ptr_t dp, mpi_size_t dsize)
{
    mpi_limb_t most_significant_q_limb = 0;

    switch(dsize) {
      case 0:
	/* We are asked to divide by zero, so go ahead and do it!  (To make
	   the compiler not remove this statement, return the value.)  */
	return 1 / dsize;

      case 1:
	{
	    mpi_size_t i;
	    mpi_limb_t n1;
	    mpi_limb_t d;

	    d = dp[0];
	    n1 = np[nsize - 1];

	    if( n1 >= d ) {
		n1 -= d;
		most_significant_q_limb = 1;
	    }

	    qp += qextra_limbs;
	    for( i = nsize - 2; i >= 0; i--)
		udiv_qrnnd( qp[i], n1, n1, np[i], d );
	    qp -= qextra_limbs;

	    for( i = qextra_limbs - 1; i >= 0; i-- )
		udiv_qrnnd (qp[i], n1, n1, 0, d);

	    np[0] = n1;
	}
	break;

      case 2:
	{
	    mpi_size_t i;
	    mpi_limb_t n1, n0, n2;
	    mpi_limb_t d1, d0;

	    np += nsize - 2;
	    d1 = dp[1];
	    d0 = dp[0];
	    n1 = np[1];
	    n0 = np[0];

	    if( n1 >= d1 && (n1 > d1 || n0 >= d0) ) {
		sub_ddmmss (n1, n0, n1, n0, d1, d0);
		most_significant_q_limb = 1;
	    }

	    for( i = qextra_limbs + nsize - 2 - 1; i >= 0; i-- ) {
		mpi_limb_t q;
		mpi_limb_t r;

		if( i >= qextra_limbs )
		    np--;
		else
		    np[0] = 0;

		if( n1 == d1 ) {
		    /* Q should be either 111..111 or 111..110.  Need special
		     * treatment of this rare case as normal division would
		     * give overflow.  */
		    q = ~(mpi_limb_t)0;

		    r = n0 + d1;
		    if( r < d1 ) {   /* Carry in the addition? */
			add_ssaaaa( n1, n0, r - d0, np[0], 0, d0 );
			qp[i] = q;
			continue;
		    }
		    n1 = d0 - (d0 != 0?1:0);
		    n0 = -d0;
		}
		else {
		    udiv_qrnnd (q, r, n1, n0, d1);
		    umul_ppmm (n1, n0, d0, q);
		}

		n2 = np[0];
	      q_test:
		if( n1 > r || (n1 == r && n0 > n2) ) {
		    /* The estimated Q was too large.  */
		    q--;
		    sub_ddmmss (n1, n0, n1, n0, 0, d0);
		    r += d1;
		    if( r >= d1 )    /* If not carry, test Q again.  */
			goto q_test;
		}

		qp[i] = q;
		sub_ddmmss (n1, n0, r, n2, n1, n0);
	    }
	    np[1] = n1;
	    np[0] = n0;
	}
	break;

      default:
	{
	    mpi_size_t i;
	    mpi_limb_t dX, d1, n0;

	    np += nsize - dsize;
	    dX = dp[dsize - 1];
	    d1 = dp[dsize - 2];
	    n0 = np[dsize - 1];

	    if( n0 >= dX ) {
		if(n0 > dX || mpihelp_cmp(np, dp, dsize - 1) >= 0 ) {
		    mpihelp_sub_n(np, np, dp, dsize);
		    n0 = np[dsize - 1];
		    most_significant_q_limb = 1;
		}
	    }

	    for( i = qextra_limbs + nsize - dsize - 1; i >= 0; i--) {
		mpi_limb_t q;
		mpi_limb_t n1, n2;
		mpi_limb_t cy_limb;

		if( i >= qextra_limbs ) {
		    np--;
		    n2 = np[dsize];
		}
		else {
		    n2 = np[dsize - 1];
		    MPN_COPY_DECR (np + 1, np, dsize - 1);
		    np[0] = 0;
		}

		if( n0 == dX ) {
		    /* This might over-estimate q, but it's probably not worth
		     * the extra code here to find out.  */
		    q = ~(mpi_limb_t)0;
		}
		else {
		    mpi_limb_t r;

		    udiv_qrnnd(q, r, n0, np[dsize - 1], dX);
		    umul_ppmm(n1, n0, d1, q);

		    while( n1 > r || (n1 == r && n0 > np[dsize - 2])) {
			q--;
			r += dX;
			if( r < dX ) /* I.e. "carry in previous addition?" */
			    break;
			n1 -= n0 < d1;
			n0 -= d1;
		    }
		}

		/* Possible optimization: We already have (q * n0) and (1 * n1)
		 * after the calculation of q.	Taking advantage of that, we
		 * could make this loop make two iterations less.  */
		cy_limb = mpihelp_submul_1(np, dp, dsize, q);

		if( n2 != cy_limb ) {
		    mpihelp_add_n(np, np, dp, dsize);
		    q--;
		}

		qp[i] = q;
		n0 = np[dsize - 1];
	    }
	}
    }

    return most_significant_q_limb;
}

/****************
 * Return the number of bits in A.
 */
unsigned
mpi_get_nbits( MPI a )
{
    unsigned n;

    mpi_normalize( a );
    if( a->nlimbs ) {
	mpi_limb_t alimb = a->d[a->nlimbs-1];
	if( alimb )
	    count_leading_zeros( n, alimb );
	else
	    n = BITS_PER_MPI_LIMB;
	n = BITS_PER_MPI_LIMB - n + (a->nlimbs-1) * BITS_PER_MPI_LIMB;
    }
    else
	n = 0;
    return n;
}

int
mpi_get_nlimbs (MPI a)
{
  return a->nlimbs;
}

int 
mpi_is_neg (MPI a)
{
  return a->sign;
}

/****************
 * Subtract the unsigned integer V from the mpi-integer U and store the
 * result in W.
 */
void
mpi_sub_ui(MPI w, MPI u, unsigned long v )
{
    mpi_ptr_t wp, up;
    mpi_size_t usize, wsize;
    int usign, wsign;

    usize = u->nlimbs;
    usign = u->sign;
    wsign = 0;

    /* If not space for W (and possible carry), increase space.  */
    wsize = usize + 1;
    if( w->alloced < wsize )
	mpi_resize(w, wsize);

    /* These must be after realloc (U may be the same as W).  */
    up = u->d;
    wp = w->d;

    if( !usize ) {  /* simple */
	wp[0] = v;
	wsize = v? 1:0;
	wsign = 1;
    }
    else if( usign ) {	/* mpi and v are negative */
	mpi_limb_t cy;
	cy = mpihelp_add_1(wp, up, usize, v);
	wp[usize] = cy;
	wsize = usize + cy;
    }
    else {  /* The signs are different.  Need exact comparison to determine
	     * which operand to subtract from which.  */
	if( usize == 1 && up[0] < v ) {
	    wp[0] = v - up[0];
	    wsize = 1;
	    wsign = 1;
	}
	else {
	    mpihelp_sub_1(wp, up, usize, v);
	    /* Size can decrease with at most one limb. */
	    wsize = usize - (wp[usize-1]==0);
	}
    }

    w->nlimbs = wsize;
    w->sign   = wsign;
}

void
mpi_sub(MPI w, MPI u, MPI v)
{
    if( w == v ) {
	MPI vv = mpi_copy(v);
	vv->sign = !vv->sign;
	mpi_add( w, u, vv );
	mpi_free(vv);
    }
    else {
	/* fixme: this is not thread-save (we temp. modify v) */
	v->sign = !v->sign;
	mpi_add( w, u, v );
	v->sign = !v->sign;
    }
}

/****************
 * Test whether bit N is set.
 */
int
mpi_test_bit( MPI a, unsigned n )
{
    unsigned limbno, bitno;
    mpi_limb_t limb;

    limbno = n / BITS_PER_MPI_LIMB;
    bitno  = n % BITS_PER_MPI_LIMB;

    if( limbno >= a->nlimbs )
	return 0; /* too far left: this is a 0 */
    limb = a->d[limbno];
    return (limb & (A_LIMB_1 << bitno))? 1: 0;
}

MPI
mpi_alloc_set_ui( unsigned long u)
{

    MPI w = mpi_alloc(1);

    w->d[0] = u;
    w->nlimbs = u? 1:0;
    w->sign = 0;
    return w;
}

/****************
 * Shift A by N bits to the right
 * FIXME: should use alloc_limb if X and A are same.
 */
void
mpi_rshift( MPI x, MPI a, unsigned n )
{
    mpi_ptr_t xp;
    mpi_size_t xsize;

    xsize = a->nlimbs;
    x->sign = a->sign;
    RESIZE_IF_NEEDED(x, xsize);
    xp = x->d;

    if( xsize ) {
	mpihelp_rshift( xp, a->d, xsize, n);
	MPN_NORMALIZE( xp, xsize);
    }
    x->nlimbs = xsize;
}

/****************
 * Calculate the multiplicative inverse X of A mod N
 * That is: Find the solution x for
 *		1 = (a*x) mod n
 */
void
mpi_invm( MPI x, MPI a, MPI n )
{
    /* Extended Euclid's algorithm (See TAOPC Vol II, 4.5.2, Alg X)
     * modified according to Michael Penk's solution for Exercice 35
     * with further enhancement */
    MPI u, v, u1, u2=NULL, u3, v1, v2=NULL, v3, t1, t2=NULL, t3;
    unsigned k;
    int sign;
    int odd ;

    u = mpi_copy(a);
    v = mpi_copy(n);

    for(k=0; !mpi_test_bit(u,0) && !mpi_test_bit(v,0); k++ ) {
	mpi_rshift(u, u, 1);
	mpi_rshift(v, v, 1);
    }
    odd = mpi_test_bit(v,0);

    u1 = mpi_alloc_set_ui(1);
    if( !odd )
	u2 = mpi_alloc_set_ui(0);
    u3 = mpi_copy(u);
    v1 = mpi_copy(v);
    if( !odd ) {
	v2 = mpi_alloc( mpi_get_nlimbs(u) );
	mpi_sub( v2, u1, u ); /* U is used as const 1 */
    }
    v3 = mpi_copy(v);
    if( mpi_test_bit(u, 0) ) { /* u is odd */
	t1 = mpi_alloc_set_ui(0);
	if( !odd ) {
	    t2 = mpi_alloc_set_ui(1); t2->sign = 1;
	}
	t3 = mpi_copy(v); t3->sign = !t3->sign;
	goto Y4;
    }
    else {
	t1 = mpi_alloc_set_ui(1);
	if( !odd )
	    t2 = mpi_alloc_set_ui(0);
	t3 = mpi_copy(u);
    }
    do {
	do {
	    if( !odd ) {
		if( mpi_test_bit(t1, 0) || mpi_test_bit(t2, 0) ) { /* one is odd */
		    mpi_add(t1, t1, v);
		    mpi_sub(t2, t2, u);
		}
		mpi_rshift(t1, t1, 1);
		mpi_rshift(t2, t2, 1);
		mpi_rshift(t3, t3, 1);
	    }
	    else {
		if( mpi_test_bit(t1, 0) )
		    mpi_add(t1, t1, v);
		mpi_rshift(t1, t1, 1);
		mpi_rshift(t3, t3, 1);
	    }
	  Y4:
	    ;
	} while( !mpi_test_bit( t3, 0 ) ); /* while t3 is even */

	if( !t3->sign ) {
	    mpi_set(u1, t1);
	    if( !odd )
		mpi_set(u2, t2);
	    mpi_set(u3, t3);
	}
	else {
	    mpi_sub(v1, v, t1);
	    sign = u->sign; u->sign = !u->sign;
	    if( !odd )
		mpi_sub(v2, u, t2);
	    u->sign = sign;
	    sign = t3->sign; t3->sign = !t3->sign;
	    mpi_set(v3, t3);
	    t3->sign = sign;
	}
	mpi_sub(t1, u1, v1);
	if( !odd )
	    mpi_sub(t2, u2, v2);
	mpi_sub(t3, u3, v3);
	if( t1->sign ) {
	    mpi_add(t1, t1, v);
	    if( !odd )
		mpi_sub(t2, t2, u);
	}
    } while( mpi_cmp_ui( t3, 0 ) ); /* while t3 != 0 */
    /* mpi_lshift( u3, k ); */
    mpi_set(x, u1);

    mpi_free(u1);
    mpi_free(v1);
    mpi_free(t1);
    if( !odd ) {
	mpi_free(u2);
	mpi_free(v2);
	mpi_free(t2);
    }
    mpi_free(u3);
    mpi_free(v3);
    mpi_free(t3);

    mpi_free(u);
    mpi_free(v);
}

/****************
 * Add the unsigned integer V to the mpi-integer U and store the
 * result in W. U and V may be the same.
 */
void
mpi_add_ui(MPI w, MPI u, unsigned long v )
{
    mpi_ptr_t wp, up;
    mpi_size_t usize, wsize;
    int usign, wsign;

    usize = u->nlimbs;
    usign = u->sign;
    wsign = 0;

    /* If not space for W (and possible carry), increase space.  */
    wsize = usize + 1;
    if( w->alloced < wsize )
	mpi_resize(w, wsize);

    /* These must be after realloc (U may be the same as W).  */
    up = u->d;
    wp = w->d;

    if( !usize ) {  /* simple */
	wp[0] = v;
	wsize = v? 1:0;
    }
    else if( !usign ) {  /* mpi is not negative */
	mpi_limb_t cy;
	cy = mpihelp_add_1(wp, up, usize, v);
	wp[usize] = cy;
	wsize = usize + cy;
    }
    else {  /* The signs are different.  Need exact comparison to determine
	     * which operand to subtract from which.  */
	if( usize == 1 && up[0] < v ) {
	    wp[0] = v - up[0];
	    wsize = 1;
	}
	else {
	    mpihelp_sub_1(wp, up, usize, v);
	    /* Size can decrease with at most one limb. */
	    wsize = usize - (wp[usize-1]==0);
	    wsign = 1;
	}
    }

    w->nlimbs = wsize;
    w->sign   = wsign;
}

void
mpi_set_ui( MPI w, unsigned long u)
{
    RESIZE_IF_NEEDED(w, 1);
    w->d[0] = u;
    w->nlimbs = u? 1:0;
    w->sign = 0;
    w->nbits = 0;
    w->flags = 0;
}

void
mpi_fdiv_qr( MPI quot, MPI rem, MPI dividend, MPI divisor )
{
    int divisor_sign = divisor->sign;
    MPI temp_divisor = NULL;

    if( quot == divisor || rem == divisor ) {
	temp_divisor = mpi_copy( divisor );
	divisor = temp_divisor;
    }

    mpi_tdiv_qr( quot, rem, dividend, divisor );

    if( (divisor_sign ^ dividend->sign) && rem->nlimbs ) {
	mpi_sub_ui( quot, quot, 1 );
	mpi_add( rem, rem, divisor);
    }

    if( temp_divisor )
	mpi_free(temp_divisor);
}

void
mpi_fdiv_q( MPI quot, MPI dividend, MPI divisor )
{
    MPI tmp = mpi_alloc( mpi_get_nlimbs(quot) );
    mpi_fdiv_qr( quot, tmp, dividend, divisor);
    mpi_free(tmp);
}

/****************
 * Find the greatest common divisor G of A and B.
 * Return: true if this 1, false in all other cases
 */
int
mpi_gcd( MPI g, MPI xa, MPI xb )
{
    MPI a, b;

    a = mpi_copy(xa);
    b = mpi_copy(xb);

    /* TAOCP Vol II, 4.5.2, Algorithm A */
    a->sign = 0;
    b->sign = 0;
    while( mpi_cmp_ui( b, 0 ) ) {
	mpi_fdiv_r( g, a, b ); /* g used as temorary variable */
	mpi_set(a,b);
	mpi_set(b,g);
    }
    mpi_set(g, a);

    mpi_free(a);
    mpi_free(b);
    return !mpi_cmp_ui( g, 1);
}

void
mpi_swap( MPI a, MPI b)
{
    struct gcry_mpi tmp;

    tmp = *a; *a = *b; *b = tmp;
}

/****************
 * Clear bit N of A.
 */
void
mpi_clear_bit( MPI a, unsigned n )
{
    unsigned limbno, bitno;

    limbno = n / BITS_PER_MPI_LIMB;
    bitno  = n % BITS_PER_MPI_LIMB;

    if( limbno >= a->nlimbs )
	return; /* don't need to clear this bit, it's to far to left */
    a->d[limbno] &= ~(A_LIMB_1 << bitno);
}

/****************
 * Set bit N of A. and clear all bits above
 */
void
mpi_set_highbit( MPI a, unsigned n )
{
    unsigned limbno, bitno;

    limbno = n / BITS_PER_MPI_LIMB;
    bitno  = n % BITS_PER_MPI_LIMB;

    if( limbno >= a->nlimbs ) { /* resize */
	if( a->alloced >= limbno )
	    mpi_resize(a, limbno+1 );
	a->nlimbs = limbno+1;
    }
    a->d[limbno] |= (A_LIMB_1<<bitno);
    for( bitno++; bitno < BITS_PER_MPI_LIMB; bitno++ )
	a->d[limbno] &= ~(A_LIMB_1 << bitno);
    a->nlimbs = limbno+1;
}

#define MPN_COPY_INCR( d, s, n) 	\
    do {				\
	mpi_size_t _i;			\
	for( _i = 0; _i < (n); _i++ )	\
	    (d)[_i] = (d)[_i];		\
    } while (0)
    
void
mpi_tdiv_q_2exp( MPI w, MPI u, unsigned count )
{
    mpi_size_t usize, wsize;
    mpi_size_t limb_cnt;

    usize = u->nlimbs;
    limb_cnt = count / BITS_PER_MPI_LIMB;
    wsize = usize - limb_cnt;
    if( limb_cnt >= usize )
	w->nlimbs = 0;
    else {
	mpi_ptr_t wp;
	mpi_ptr_t up;

	RESIZE_IF_NEEDED( w, wsize );
	wp = w->d;
	up = u->d;

	count %= BITS_PER_MPI_LIMB;
	if( count ) {
	    mpihelp_rshift( wp, up + limb_cnt, wsize, count );
	    wsize -= !wp[wsize - 1];
	}
	else {
	    MPN_COPY_INCR( wp, up + limb_cnt, wsize);
	}

	w->nlimbs = wsize;
    }
}

#define count_trailing_zeros(count, x) \
  do {									\
    UWtype __ctz_x = (x);						\
    UWtype __ctz_c;							\
    count_leading_zeros (__ctz_c, __ctz_x & -__ctz_x);			\
    (count) = W_TYPE_SIZE - 1 - __ctz_c;				\
  } while (0)

/****************
 * Count the number of zerobits at the low end of A
 */
unsigned int
mpi_trailing_zeros( MPI a )
{
    unsigned n, count = 0;

    for(n=0; n < a->nlimbs; n++ ) {
	if( a->d[n] ) {
	    unsigned nn;
	    mpi_limb_t alimb = a->d[n];

	    count_trailing_zeros( nn, alimb );
	    count += nn;
	    break;
	}
	count += BITS_PER_MPI_LIMB;
    }
    return count;

}

/****************
 * Return true if n is probably a prime
 */
static int
is_prime( MPI n, int steps, int *count )
{
    MPI x = mpi_alloc( mpi_get_nlimbs( n ) );
    MPI y = mpi_alloc( mpi_get_nlimbs( n ) );
    MPI z = mpi_alloc( mpi_get_nlimbs( n ) );
    MPI nminus1 = mpi_alloc( mpi_get_nlimbs( n ) );
    MPI a2 = mpi_alloc_set_ui( 2 );
    MPI q;
    unsigned i, j, k;
    int rc = 0;
    unsigned nbits = mpi_get_nbits( n );

    mpi_sub_ui( nminus1, n, 1 );

    /* find q and k, so that n = 1 + 2^k * q */
    q = mpi_copy( nminus1 );
    k = mpi_trailing_zeros( q );
    mpi_tdiv_q_2exp(q, q, k);

    for(i=0 ; i < steps; i++ ) {
	++*count;
	if( !i ) {
	    mpi_set_ui( x, 2 );
	}
	else {
            unsigned char *p;
            
            p = get_random_bits( nbits, 0, 0 );
            mpi_set_buffer( x, p, (nbits+7)/8, 0 );
            free(p);

	    /* Make sure that the number is smaller than the prime
	     * and keep the randomness of the high bit.  */
	    if( mpi_test_bit( x, nbits-2 ) ) {
		mpi_set_highbit( x, nbits-2 ); /* Clear all higher bits */
	    }
	    else {
		mpi_set_highbit( x, nbits-2 );
		mpi_clear_bit( x, nbits-2 );
	    }
	}
	mpi_powm( y, x, q, n);
	if( mpi_cmp_ui(y, 1) && mpi_cmp( y, nminus1 ) ) {
	    for( j=1; j < k && mpi_cmp( y, nminus1 ); j++ ) {
		mpi_powm(y, y, a2, n);
		if( !mpi_cmp_ui( y, 1 ) )
		    goto leave; /* not a prime */
	    }
	    if( mpi_cmp( y, nminus1 ) )
		goto leave; /* not a prime */
	}
    }
    rc = 1; /* may be a prime */

  leave:
    mpi_free( x );
    mpi_free( y );
    mpi_free( z );
    mpi_free( nminus1 );
    mpi_free( q );
    mpi_free (a2);

    return rc;
}

MPI
mpi_alloc_like( MPI a )
{
    MPI b;

    if( a && (a->flags & 4) ) {
	void *p = m_is_secure(a->d)? malloc( a->nbits )
				   : malloc( a->nbits );
	memcpy( p, a->d, a->nbits );
	b = mpi_set_opaque( NULL, p, a->nbits );
    }
    else if( a ) {

	b = mpi_alloc( a->nlimbs );

	b->nlimbs = 0;
	b->sign = 0;
	b->flags = a->flags;
	b->nbits = 0;
    }
    else
	b = NULL;
    return b;
}

/****************
 * Set bit N of A.
 */
void
mpi_set_bit( MPI a, unsigned n )
{
    unsigned limbno, bitno;

    limbno = n / BITS_PER_MPI_LIMB;
    bitno  = n % BITS_PER_MPI_LIMB;

    if( limbno >= a->nlimbs ) { /* resize */
	if( a->alloced >= limbno )
	    mpi_resize(a, limbno+1 );
	a->nlimbs = limbno+1;
    }
    a->d[limbno] |= (A_LIMB_1<<bitno);
}


static MPI
gen_prime( unsigned int nbits, int secret, int randomlevel )
{
    unsigned  nlimbs;
    MPI prime, ptest, pminus1, val_2, val_3, result;
    int i;
    unsigned x, step;
    int count1, count2;
    int *mods;

    if( !no_of_small_prime_numbers ) {
		for(i=0; small_prime_numbers[i]; i++ )
	    	no_of_small_prime_numbers++;
    }
    mods = (int *)malloc( no_of_small_prime_numbers * sizeof *mods );

    /* Make nbits fit into MPI implementation.  */
    nlimbs = mpi_nlimb_hint_from_nbits (nbits);
    val_2  = mpi_alloc_set_ui( 2 );
    val_3 = mpi_alloc_set_ui( 3);
    prime  = mpi_alloc( nlimbs );
    result = mpi_alloc_like( prime );
    pminus1= mpi_alloc_like( prime );
    ptest  = mpi_alloc_like( prime );
    count1 = count2 = 0; 

    for(;;) {  /* try forvever */
		int dotcount=0;

		/* generate a random number */
		{   
			unsigned char *p = get_random_bits( nbits, randomlevel, secret );
			mpi_set_buffer( prime, p, (nbits+7)/8, 0 );
			free(p);
		}

		/* Set high order bit to 1, set low order bit to 0.
		       If we are generating a secret prime we are most probably
		       doing that for RSA, to make sure that the modulus does have
		       the requested keysize we set the 2 high order bits */
		mpi_set_highbit( prime, nbits-1 );
		    if (secret)
		      mpi_set_bit (prime, nbits-2);
		mpi_set_bit( prime, 0 );

		/* calculate all remainders */
		for(i=0; (x = small_prime_numbers[i]); i++ )
			mods[i] = mpi_fdiv_r_ui(NULL, prime, x);
////MYDEBUG(11);
		/* now try some primes starting with prime */
		for(step=0; step < 20000; step += 2 ) {
			/* check against all the small primes we have in mods */
			count1++;
			for(i=0; (x = small_prime_numbers[i]); i++ ) {
				while( mods[i] + step >= x )
					mods[i] -= x;
				if( !(mods[i] + step) )
					break;
			}
			if( x )
			continue;   /* found a multiple of an already known prime */

			mpi_add_ui( ptest, prime, step );

			/* do a faster Fermat test */
			count2++;
			mpi_sub_ui( pminus1, ptest, 1);

			mpi_powm( result, val_2, pminus1, ptest );

			if( !mpi_cmp_ui( result, 1 ) ) { /* not composite */
				/* perform stronger tests */

				if( is_prime(ptest, 5, &count2 ) ) {
					if( !mpi_test_bit( ptest, nbits-1 ) ) {
				
					break; /* step loop, continue with a new prime */
					}

					mpi_free(val_2);
					mpi_free(val_3);
					mpi_free(result);
					mpi_free(pminus1);
					mpi_free(prime);
					free(mods);
					return ptest;
				}

			}
			if( ++dotcount == 10 ) {
		
			dotcount = 0;
			}
		}
		}
}
/****************
 * Generate a prime number (stored in secure memory)
 */
MPI
generate_secret_prime( unsigned  nbits )
{
    MPI prime;
	////MYDEBUG(3);
    prime = gen_prime( nbits, 1, 2 );
    ////MYDEBUG(3);
    return prime;
}


/****************
 * Use BUFFER to update MPI.
 */
void
mpi_set_buffer( MPI a, const byte *buffer, unsigned nbytes, int sign )
{
    const byte *p;
    mpi_limb_t alimb;
    int nlimbs;
    int i;

    nlimbs = (nbytes + BYTES_PER_MPI_LIMB - 1) / BYTES_PER_MPI_LIMB;
    RESIZE_IF_NEEDED(a, nlimbs);
    a->sign = sign;

    for(i=0, p = buffer+nbytes-1; p >= buffer+BYTES_PER_MPI_LIMB; ) {
    alimb = 0;
#if BYTES_PER_MPI_LIMB == 4
	alimb |= (mpi_limb_t)*p-- ;
	alimb |= (mpi_limb_t)*p-- <<  8 ;
	alimb |= (mpi_limb_t)*p-- << 16 ;
	alimb |= (mpi_limb_t)*p-- << 24 ;
#elif BYTES_PER_MPI_LIMB == 8
	alimb |= (mpi_limb_t)*p--	;     //TRACE(13, alimb); TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- <<  8 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 16 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 24 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 32 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 40 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 48 ; //TRACE(13, alimb);TRACE(13, *p);
	alimb |= (mpi_limb_t)*p-- << 56 ; 
	//TRACE(14, alimb);
#else
#error please implement for this limb size.
#endif
	a->d[i++] = alimb;
    }
    if( p >= buffer ) {
    alimb = 0;
#if BYTES_PER_MPI_LIMB == 4
	alimb  |= *p--	    ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- <<  8 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 16 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 24 ;
#elif BYTES_PER_MPI_LIMB == 8
	alimb  |= (mpi_limb_t)*p-- ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- <<	8 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 16 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 24 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 32 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 40 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 48 ;
	if( p >= buffer ) alimb |= (mpi_limb_t)*p-- << 56 ;
#else
#error please implement for this limb size.
#endif
	a->d[i++] = alimb;
    }
    a->nlimbs = i;
}


/****************
 * Return a pointer to a randomized buffer of level 0 and LENGTH bits
 * caller must free the buffer.
 * Note: The returned value is rounded up to bytes.
 */
byte *
get_random_bits( size_t nbits, int level, int secure )
{
	size_t nbytes = (nbits+7)/8;
    byte *buf = new byte[nbytes];
    
	sgx_read_rand((unsigned char*)buf, nbytes);
	
	return buf;
}

#define mpi_is_secure(a) ((a) && (mpi_get_flags (a)&1))

mpi_limb_t
mpihelp_mul_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size,
						    mpi_limb_t s2_limb)
{
    mpi_limb_t cy_limb;
    mpi_size_t j;
    mpi_limb_t prod_high, prod_low;

    /* The loop counter and index J goes from -S1_SIZE to -1.  This way
     * the loop becomes faster.  */
    j = -s1_size;

    /* Offset the base pointers to compensate for the negative indices.  */
    s1_ptr -= j;
    res_ptr -= j;

    cy_limb = 0;
    do {
	umul_ppmm( prod_high, prod_low, s1_ptr[j], s2_limb );
	prod_low += cy_limb;
	cy_limb = (prod_low < cy_limb?1:0) + prod_high;
	res_ptr[j] = prod_low;
    } while( ++j );

    return cy_limb;
}

mpi_limb_t
mpihelp_addmul_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr,
		  mpi_size_t s1_size, mpi_limb_t s2_limb)
{
    mpi_limb_t cy_limb;
    mpi_size_t j;
    mpi_limb_t prod_high, prod_low;
    mpi_limb_t x;

    /* The loop counter and index J goes from -SIZE to -1.  This way
     * the loop becomes faster.  */
    j = -s1_size;
    res_ptr -= j;
    s1_ptr -= j;

    cy_limb = 0;
    do {
	umul_ppmm( prod_high, prod_low, s1_ptr[j], s2_limb );

	prod_low += cy_limb;
	cy_limb = (prod_low < cy_limb?1:0) + prod_high;

	x = res_ptr[j];
	prod_low = x + prod_low;
	cy_limb += prod_low < x?1:0;
	res_ptr[j] = prod_low;
    } while ( ++j );
    return cy_limb;
}

//int count = 0;
//bool flagprint = false;
void
mpih_sqr_n_basecase( mpi_ptr_t prodp, mpi_ptr_t up, mpi_size_t size )
{
//	if(flagprint) MYDEBUG(count++);
    mpi_size_t i;
    mpi_limb_t cy_limb;
    mpi_limb_t v_limb;
	
    /* Multiply by the first limb in V separately, as the result can be
     * stored (not added) to PROD.  We also avoid a loop for zeroing.  */
    v_limb = up[0];
    if( v_limb <= 1 ) {
		if( v_limb == 1 )
			MPN_COPY( prodp, up, size );  //size == 3 here
		else
			MPN_ZERO(prodp, size);
		cy_limb = 0;
    }
    else
	cy_limb = mpihelp_mul_1( prodp, up, size, v_limb );

    prodp[size] = cy_limb;
    prodp++;

    /* For each iteration in the outer loop, multiply one limb from
     * U with one limb from V, and add it to PROD.  */
    for( i=1; i < size; i++) {
		v_limb = up[i];
		if( v_limb <= 1 ) {
			cy_limb = 0;
			if( v_limb == 1 )
			cy_limb = mpihelp_add_n(prodp, prodp, up, size);
		}
		else
			cy_limb = mpihelp_addmul_1(prodp, up, size, v_limb);

		prodp[size] = cy_limb;
		prodp++;
    }
}

void
mpi_mulm( MPI w, MPI u, MPI v, MPI m)
{
    mpi_mul(w, u, v);
    mpi_fdiv_r( w, w, m );
}


void
mpi_mul( MPI w, MPI u, MPI v)
{
    mpi_size_t usize, vsize, wsize;
    mpi_ptr_t up, vp, wp;
    mpi_limb_t cy;
    int usign, vsign, usecure, vsecure, sign_product;
    int assign_wp=0;
    mpi_ptr_t tmp_limb=NULL;


    if( u->nlimbs < v->nlimbs ) { /* Swap U and V. */
	usize = v->nlimbs;
	usign = v->sign;
	usecure = mpi_is_secure(v);
	up    = v->d;
	vsize = u->nlimbs;
	vsign = u->sign;
	vsecure = mpi_is_secure(u);
	vp    = u->d;
    }
    else {
	usize = u->nlimbs;
	usign = u->sign;
	usecure = mpi_is_secure(u);
	up    = u->d;
	vsize = v->nlimbs;
	vsign = v->sign;
	vsecure = mpi_is_secure(v);
	vp    = v->d;
    }
    sign_product = usign ^ vsign;
    wp = w->d;

    /* Ensure W has space enough to store the result.  */
    wsize = usize + vsize;
    if ( !mpi_is_secure (w) && (mpi_is_secure (u) || mpi_is_secure (v)) ) {
        /* w is not allocated in secure space but u or v is.  To make sure
         * that no temporray results are stored in w, we temporary use 
         * a newly allocated limb space for w */
        wp = mpi_alloc_limb_space( wsize, 1 );
        assign_wp = 2; /* mark it as 2 so that we can later copy it back to
                        * mormal memory */
    }
    else if( w->alloced < wsize ) {
	if( wp == up || wp == vp ) {
	    wp = mpi_alloc_limb_space( wsize, mpi_is_secure(w) );
	    assign_wp = 1;
	}
	else {
	    mpi_resize(w, wsize );
	    wp = w->d;
	}
    }
    else { /* Make U and V not overlap with W.	*/
	if( wp == up ) {
	    /* W and U are identical.  Allocate temporary space for U.	*/
	    up = tmp_limb = mpi_alloc_limb_space( usize, usecure  );
	    /* Is V identical too?  Keep it identical with U.  */
	    if( wp == vp )
		vp = up;
	    /* Copy to the temporary space.  */
	    MPN_COPY( up, wp, usize );
	}
	else if( wp == vp ) {
	    /* W and V are identical.  Allocate temporary space for V.	*/
	    vp = tmp_limb = mpi_alloc_limb_space( vsize, vsecure );
	    /* Copy to the temporary space.  */
	    MPN_COPY( vp, wp, vsize );
	}
    }

    if( !vsize )
	wsize = 0;
    else {
	cy = mpihelp_mul( wp, up, usize, vp, vsize );
	wsize -= cy? 0:1;
    }

    if( assign_wp ) {
        if (assign_wp == 2) {
            /* copy the temp wp from secure memory back to normal memory */
	    mpi_ptr_t tmp_wp = mpi_alloc_limb_space (wsize, 0);
	    MPN_COPY (tmp_wp, wp, wsize);
            mpi_free_limb_space (wp);
            wp = tmp_wp;
        }
	mpi_assign_limb_space( w, wp, wsize );
    }
    w->nlimbs = wsize;
    w->sign = sign_product;
    if( tmp_limb )
	mpi_free_limb_space( tmp_limb );
}
  
/* Multiply the natural numbers u (pointed to by UP, with USIZE limbs)
 * and v (pointed to by VP, with VSIZE limbs), and store the result at
 * PRODP.  USIZE + VSIZE limbs are always stored, but if the input
 * operands are normalized.  Return the most significant limb of the
 * result.
 *
 * NOTE: The space pointed to by PRODP is overwritten before finished
 * with U and V, so overlap is an error.
 *
 * Argument constraints:
 * 1. USIZE >= VSIZE.
 * 2. PRODP != UP and PRODP != VP, i.e. the destination
 *    must be distinct from the multiplier and the multiplicand.
 */

mpi_limb_t
mpihelp_mul( mpi_ptr_t prodp, mpi_ptr_t up, mpi_size_t usize,
			      mpi_ptr_t vp, mpi_size_t vsize)
{
    mpi_ptr_t prod_endp = prodp + usize + vsize - 1;
    mpi_limb_t cy;
    struct karatsuba_ctx ctx;

    if( vsize < KARATSUBA_THRESHOLD ) {
	mpi_size_t i;
	mpi_limb_t v_limb;

	if( !vsize )
	    return 0;

	/* Multiply by the first limb in V separately, as the result can be
	 * stored (not added) to PROD.	We also avoid a loop for zeroing.  */
	v_limb = vp[0];
	if( v_limb <= 1 ) {
	    if( v_limb == 1 )
		MPN_COPY( prodp, up, usize );
	    else
		MPN_ZERO( prodp, usize );
	    cy = 0;
	}
	else
	    cy = mpihelp_mul_1( prodp, up, usize, v_limb );

	prodp[usize] = cy;
	prodp++;

	/* For each iteration in the outer loop, multiply one limb from
	 * U with one limb from V, and add it to PROD.	*/
	for( i = 1; i < vsize; i++ ) {
	    v_limb = vp[i];
	    if( v_limb <= 1 ) {
		cy = 0;
		if( v_limb == 1 )
		   cy = mpihelp_add_n(prodp, prodp, up, usize);
	    }
	    else
		cy = mpihelp_addmul_1(prodp, up, usize, v_limb);

	    prodp[usize] = cy;
	    prodp++;
	}

	return cy;
    }

    memset( &ctx, 0, sizeof ctx );
    mpihelp_mul_karatsuba_case( prodp, up, usize, vp, vsize, &ctx );
    mpihelp_release_karatsuba_ctx( &ctx );
    return *prod_endp;
}

void
mpih_sqr_n( mpi_ptr_t prodp, mpi_ptr_t up, mpi_size_t size, mpi_ptr_t tspace)
{
    if( size & 1 ) {
	/* The size is odd, and the code below doesn't handle that.
	 * Multiply the least significant (size - 1) limbs with a recursive
	 * call, and handle the most significant limb of S1 and S2
	 * separately.
	 * A slightly faster way to do this would be to make the Karatsuba
	 * code below behave as if the size were even, and let it check for
	 * odd size in the end.  I.e., in essence move this code to the end.
	 * Doing so would save us a recursive call, and potentially make the
	 * stack grow a lot less.
	 */
	mpi_size_t esize = size - 1;	   /* even size */
	mpi_limb_t cy_limb;

	MPN_SQR_N_RECURSE( prodp, up, esize, tspace );
	cy_limb = mpihelp_addmul_1( prodp + esize, up, esize, up[esize] );
	prodp[esize + esize] = cy_limb;
	cy_limb = mpihelp_addmul_1( prodp + esize, up, size, up[esize] );

	prodp[esize + size] = cy_limb;
    }
    else {
	mpi_size_t hsize = size >> 1;
	mpi_limb_t cy;

	/* Product H.	   ________________  ________________
	 *		  |_____U1 x U1____||____U0 x U0_____|
	 * Put result in upper part of PROD and pass low part of TSPACE
	 * as new TSPACE.
	 */
	MPN_SQR_N_RECURSE(prodp + size, up + hsize, hsize, tspace);

	/* Product M.	   ________________
	 *		  |_(U1-U0)(U0-U1)_|
	 */
	if( mpihelp_cmp( up + hsize, up, hsize) >= 0 )
	    mpihelp_sub_n( prodp, up + hsize, up, hsize);
	else
	    mpihelp_sub_n (prodp, up, up + hsize, hsize);

	/* Read temporary operands from low part of PROD.
	 * Put result in low part of TSPACE using upper part of TSPACE
	 * as new TSPACE.  */
	MPN_SQR_N_RECURSE(tspace, prodp, hsize, tspace + size);

	/* Add/copy product H  */
	MPN_COPY(prodp + hsize, prodp + size, hsize);
	cy = mpihelp_add_n(prodp + size, prodp + size,
			   prodp + size + hsize, hsize);

	/* Add product M (if NEGFLG M is a negative number).  */
	cy -= mpihelp_sub_n (prodp + hsize, prodp + hsize, tspace, size);

	/* Product L.	   ________________  ________________
	 *		  |________________||____U0 x U0_____|
	 * Read temporary operands from low part of PROD.
	 * Put result in low part of TSPACE using upper part of TSPACE
	 * as new TSPACE.  */
	MPN_SQR_N_RECURSE (tspace, up, hsize, tspace + size);

	/* Add/copy Product L (twice).	*/
	cy += mpihelp_add_n (prodp + hsize, prodp + hsize, tspace, size);
	if( cy )
	    mpihelp_add_1(prodp + hsize + size, prodp + hsize + size,
							    hsize, cy);

	MPN_COPY(prodp, tspace, hsize);
	cy = mpihelp_add_n (prodp + hsize, prodp + hsize, tspace + hsize, hsize);
	if( cy )
	    mpihelp_add_1 (prodp + size, prodp + size, size, 1);
    }
}


void
mpi_assign_limb_space( MPI a, mpi_ptr_t ap, unsigned nlimbs )
{
    mpi_free_limb_space(a->d);
    a->d = ap;
    a->alloced = nlimbs;
}

/* Shift U (pointed to by UP and USIZE digits long) CNT bits to the left
 * and store the USIZE least significant digits of the result at WP.
 * Return the bits shifted out from the most significant digit.
 *
 * Argument constraints:
 * 1. 0 < CNT < BITS_PER_MP_LIMB
 * 2. If the result is to be written over the input, WP must be >= UP.
 */
mpi_limb_t
mpihelp_lshift( mpi_ptr_t wp, mpi_ptr_t up, mpi_size_t usize,
					    unsigned int cnt)
{
    mpi_limb_t high_limb, low_limb;
    unsigned sh_1, sh_2;
    mpi_size_t i;
    mpi_limb_t retval;

    sh_1 = cnt;
    wp += 1;
    sh_2 = BITS_PER_MPI_LIMB - sh_1;
    i = usize - 1;
    low_limb = up[i];
    retval = low_limb >> sh_2;
    high_limb = low_limb;
    while( --i >= 0 ) {
	low_limb = up[i];
	wp[i] = (high_limb << sh_1) | (low_limb >> sh_2);
	high_limb = low_limb;
    }
    wp[i] = high_limb << sh_1;

    return retval;
}

/* Shift U (pointed to by UP and USIZE limbs long) CNT bits to the right
 * and store the USIZE least significant limbs of the result at WP.
 * The bits shifted out to the right are returned.
 *
 * Argument constraints:
 * 1. 0 < CNT < BITS_PER_MP_LIMB
 * 2. If the result is to be written over the input, WP must be <= UP.
 */
mpi_limb_t
mpihelp_rshift( mpi_ptr_t wp, mpi_ptr_t up, mpi_size_t usize, unsigned cnt)
{
    mpi_limb_t high_limb, low_limb;
    unsigned sh_1, sh_2;
    mpi_size_t i;
    mpi_limb_t retval;

    sh_1 = cnt;
    wp -= 1;
    sh_2 = BITS_PER_MPI_LIMB - sh_1;
    high_limb = up[0];
    retval = high_limb << sh_2;
    low_limb = high_limb;
    for( i=1; i < usize; i++) {
	high_limb = up[i];
	wp[i] = (low_limb >> sh_1) | (high_limb << sh_2);
	low_limb = high_limb;
    }
    wp[i] = low_limb >> sh_1;

    return retval;
}

/****************
 * RES = BASE ^ EXP mod MOD
 */
void
mpi_powm( MPI res, MPI base, MPI exponent, MPI mod)
{
    mpi_ptr_t  rp, ep, mp, bp;
    mpi_size_t esize, msize, bsize, rsize;
    int               msign, bsign, rsign;
    int        esec,  msec,  bsec,  rsec;
    mpi_size_t size;
    int mod_shift_cnt;
    int negative_result;
    mpi_ptr_t mp_marker=NULL, bp_marker=NULL, ep_marker=NULL;
    mpi_ptr_t xp_marker=NULL;
    int assign_rp=0;
    mpi_ptr_t tspace = NULL;
    mpi_size_t tsize=0;   /* to avoid compiler warning */
			  /* fixme: we should check that the warning is void*/

    esize = exponent->nlimbs;
    msize = mod->nlimbs;
    size = 2 * msize;
    msign = mod->sign;

    esec = mpi_is_secure(exponent);
    msec = mpi_is_secure(mod);
    bsec = mpi_is_secure(base);
    rsec = mpi_is_secure(res);

    rp = res->d;
    ep = exponent->d;

    if( !msize )
		msize = 1 / msize;	    /* provoke a signal */

    if( !esize ) {
		/* Exponent is zero, result is 1 mod MOD, i.e., 1 or 0
		 * depending on if MOD equals 1.  */
		rp[0] = 1;
		res->nlimbs = (msize == 1 && mod->d[0] == 1) ? 0 : 1;
		res->sign = 0;
		goto leave;
    }
	
    /* Normalize MOD (i.e. make its most significant bit set) as required by
     * mpn_divrem.  This will make the intermediate values in the calculation
     * slightly larger, but the correct result is obtained after a final
     * reduction using the original MOD value.	*/
    mp = mp_marker = mpi_alloc_limb_space(msize, msec); 
    count_leading_zeros( mod_shift_cnt, mod->d[msize-1] );
    if( mod_shift_cnt )
		mpihelp_lshift( mp, mod->d, msize, mod_shift_cnt );
    else
		MPN_COPY( mp, mod->d, msize );

    bsize = base->nlimbs;
    bsign = base->sign;
    if( bsize > msize ) { /* The base is larger than the module. Reduce it. */
	/* Allocate (BSIZE + 1) with space for remainder and quotient.
	 * (The quotient is (bsize - msize + 1) limbs.)  */
	bp = bp_marker = mpi_alloc_limb_space( bsize + 1, bsec );
	MPN_COPY( bp, base->d, bsize );
	/* We don't care about the quotient, store it above the remainder,
	 * at BP + MSIZE.  */
	mpihelp_divrem( bp + msize, 0, bp, bsize, mp, msize );
	bsize = msize;
	/* Canonicalize the base, since we are going to multiply with it
	 * quite a few times.  */
	MPN_NORMALIZE( bp, bsize );
    }
    else
	bp = base->d;

    if( !bsize ) {
	res->nlimbs = 0;
	res->sign = 0;
	goto leave;
    }

    if( res->alloced < size ) {
	/* We have to allocate more space for RES.  If any of the input
	 * parameters are identical to RES, defer deallocation of the old
	 * space.  */
	if( rp == ep || rp == mp || rp == bp ) {
	    rp = mpi_alloc_limb_space( size, rsec );
	    assign_rp = 1;
	}
	else {
	    mpi_resize( res, size );
	    rp = res->d;
	}
    }
    else { /* Make BASE, EXPONENT and MOD not overlap with RES.  */
	if( rp == bp ) {
	    /* RES and BASE are identical.  Allocate temp. space for BASE.  */
	    bp = bp_marker = mpi_alloc_limb_space( bsize, bsec );
	    MPN_COPY(bp, rp, bsize);
	}
	if( rp == ep ) {
	    /* RES and EXPONENT are identical.
               Allocate temp. space for EXPONENT.  */
	    ep = ep_marker = mpi_alloc_limb_space( esize, esec );
	    MPN_COPY(ep, rp, esize);
	}
	if( rp == mp ) {
	    /* RES and MOD are identical.  Allocate temporary space for MOD.*/
	    mp = mp_marker = mpi_alloc_limb_space( msize, msec );
	    MPN_COPY(mp, rp, msize);
	}
    }

    MPN_COPY( rp, bp, bsize );
    rsize = bsize;
    rsign = bsign;

    {
	mpi_size_t i;
	mpi_ptr_t xp = xp_marker = mpi_alloc_limb_space( 2 * (msize + 1), msec );
	int c;
	mpi_limb_t e;
	mpi_limb_t carry_limb;
	struct karatsuba_ctx karactx;

	memset( &karactx, 0, sizeof karactx );
	negative_result = (ep[0] & 1) && base->sign;

	i = esize - 1;
	e = ep[i];
	count_leading_zeros (c, e);
	e = (e << c) << 1;     /* shift the exp bits to the left, lose msb */
	c = BITS_PER_MPI_LIMB - 1 - c;

	/* Main loop.
	 *
	 * Make the result be pointed to alternately by XP and RP.  This
	 * helps us avoid block copying, which would otherwise be necessary
	 * with the overlap restrictions of mpihelp_divmod. With 50% probability
	 * the result after this loop will be in the area originally pointed
	 * by RP (==RES->d), and with 50% probability in the area originally
	 * pointed to by XP.
	 */

	for(;;) {
	    while( c ) {
		mpi_ptr_t tp;
		mpi_size_t xsize;

		/*mpihelp_mul_n(xp, rp, rp, rsize);*/
		if( rsize < KARATSUBA_THRESHOLD )
		    mpih_sqr_n_basecase( xp, rp, rsize );
		else {
		    if( !tspace ) {
			tsize = 2 * rsize;
			tspace = mpi_alloc_limb_space( tsize, 0 );
		    }
		    else if( tsize < (2*rsize) ) {
			mpi_free_limb_space( tspace );
			tsize = 2 * rsize;
			tspace = mpi_alloc_limb_space( tsize, 0 );
		    }
		    mpih_sqr_n( xp, rp, rsize, tspace );
		}

		xsize = 2 * rsize;
		if( xsize > msize ) {
		    mpihelp_divrem(xp + msize, 0, xp, xsize, mp, msize);
		    xsize = msize;
		}

		tp = rp; rp = xp; xp = tp;
		rsize = xsize;

		if( (mpi_limb_signed_t)e < 0 ) {
		    /*mpihelp_mul( xp, rp, rsize, bp, bsize );*/
		    if( bsize < KARATSUBA_THRESHOLD ) {
			mpihelp_mul( xp, rp, rsize, bp, bsize );
		    }
		    else {
			mpihelp_mul_karatsuba_case(
				     xp, rp, rsize, bp, bsize, &karactx );
		    }

		    xsize = rsize + bsize;
		    if( xsize > msize ) {
			mpihelp_divrem(xp + msize, 0, xp, xsize, mp, msize);
			xsize = msize;
		    }

		    tp = rp; rp = xp; xp = tp;
		    rsize = xsize;
		}
		e <<= 1;
		c--;
	    }

	    i--;
	    if( i < 0 )
		break;
	    e = ep[i];
	    c = BITS_PER_MPI_LIMB;
	}

	/* We shifted MOD, the modulo reduction argument, left MOD_SHIFT_CNT
	 * steps.  Adjust the result by reducing it with the original MOD.
	 *
	 * Also make sure the result is put in RES->d (where it already
	 * might be, see above).
	 */
	if( mod_shift_cnt ) {
	    carry_limb = mpihelp_lshift( res->d, rp, rsize, mod_shift_cnt);
	    rp = res->d;
	    if( carry_limb ) {
		rp[rsize] = carry_limb;
		rsize++;
	    }
	}
	else {
	    MPN_COPY( res->d, rp, rsize);
	    rp = res->d;
	}

	if( rsize >= msize ) {
	    mpihelp_divrem(rp + msize, 0, rp, rsize, mp, msize);
	    rsize = msize;
	}

	/* Remove any leading zero words from the result.  */
	if( mod_shift_cnt )
	    mpihelp_rshift( rp, rp, rsize, mod_shift_cnt);
	MPN_NORMALIZE (rp, rsize);

	mpihelp_release_karatsuba_ctx( &karactx );
    }

    if( negative_result && rsize ) {
	if( mod_shift_cnt )
	    mpihelp_rshift( mp, mp, msize, mod_shift_cnt);
	mpihelp_sub( rp, mp, msize, rp, rsize);
	rsize = msize;
	rsign = msign;
	MPN_NORMALIZE(rp, rsize);
    }
    res->nlimbs = rsize;
    res->sign = rsign;

  leave:
    if( assign_rp ) mpi_assign_limb_space( res, rp, size );
    if( mp_marker ) mpi_free_limb_space( mp_marker );
    if( bp_marker ) mpi_free_limb_space( bp_marker );
    if( ep_marker ) mpi_free_limb_space( ep_marker );
    if( xp_marker ) mpi_free_limb_space( xp_marker );
    if( tspace )    mpi_free_limb_space( tspace );
}

/****************
 * clear bit N of A and all bits above
 */
void
mpi_clear_highbit( MPI a, unsigned n )
{
    unsigned limbno, bitno;

    limbno = n / BITS_PER_MPI_LIMB;
    bitno  = n % BITS_PER_MPI_LIMB;

    if( limbno >= a->nlimbs )
	return; /* not allocated, so need to clear bits :-) */

    for( ; bitno < BITS_PER_MPI_LIMB; bitno++ )
	a->d[limbno] &= ~(A_LIMB_1 << bitno);
    a->nlimbs = limbno+1;
}

static void
m_out_of_n( unsigned char *array, int m, int n )
{
    int i=0, i1=0, j=0, jp=0,  j1=0, k1=0, k2=0;

    if( !m || m >= n )
	return;

    if( m == 1 ) { /* special case */
	for(i=0; i < n; i++ )
	    if( array[i] ) {
		array[i++] = 0;
		if( i >= n )
		    i = 0;
		array[i] = 1;
		return;
	    }
    }

    for(j=1; j < n; j++ ) {
	if( array[n-1] == array[n-j-1] )
	    continue;
	j1 = j;
	break;
    }

    if( m & 1 ) { /* m is odd */
	if( array[n-1] ) {
	    if( j1 & 1 ) {
		k1 = n - j1;
		k2 = k1+2;
		if( k2 > n )
		    k2 = n;
		goto leave;
	    }
	    goto scan;
	}
	k2 = n - j1 - 1;
	if( k2 == 0 ) {
	    k1 = i;
	    k2 = n - j1;
	}
	else if( array[k2] && array[k2-1] )
	    k1 = n;
	else
	    k1 = k2 + 1;
    }
    else { /* m is even */
	if( !array[n-1] ) {
	    k1 = n - j1;
	    k2 = k1 + 1;
	    goto leave;
	}

	if( !(j1 & 1) ) {
	    k1 = n - j1;
	    k2 = k1+2;
	    if( k2 > n )
		k2 = n;
	    goto leave;
	}
scan:
	jp = n - j1 - 1;
	for(i=1; i <= jp; i++ ) {
	    i1 = jp + 2 - i;
	    if( array[i1-1]  ) {
		if( array[i1-2] ) {
		    k1 = i1 - 1;
		    k2 = n - j1;
		}
		else {
		    k1 = i1 - 1;
		    k2 = n + 1 - j1;
		}
		goto leave;
	    }
	}
	k1 = 1;
	k2 = n + 1 - m;
    }
  leave:
    array[k1-1] = !array[k1-1];
    array[k2-1] = !array[k2-1];
}

void
mpi_mul_ui( MPI prod, MPI mult, unsigned long small_mult )
{
    mpi_size_t size, prod_size;
    mpi_ptr_t  prod_ptr;
    mpi_limb_t cy;
    int sign;

    size = mult->nlimbs;
    sign = mult->sign;

    if( !size || !small_mult ) {
	prod->nlimbs = 0;
	prod->sign = 0;
	return;
    }

    prod_size = size + 1;
    if( prod->alloced < prod_size )
	mpi_resize( prod, prod_size );
    prod_ptr = prod->d;

    cy = mpihelp_mul_1( prod_ptr, mult->d, size, (mpi_limb_t)small_mult );
    if( cy )
	prod_ptr[size++] = cy;
    prod->nlimbs = size;
    prod->sign = sign;
}

/****************
 * Check whether dividend is divisible by divisor
 * (note: divisor must fit into a limb)
 */
int
mpi_divisible_ui(MPI dividend, unsigned long divisor )
{
    return !mpihelp_mod_1( dividend->d, dividend->nlimbs, divisor );
}

/****************
 * Returns: true if this may be a prime
 */
static int
check_prime( MPI prime, MPI val_2 )
{
    int i;
    unsigned x;
    int count=0;

    /* check against small primes */
    for(i=0; (x = small_prime_numbers[i]); i++ ) {
	if( mpi_divisible_ui( prime, x ) )
	    return 0;
    }

    /* a quick fermat test */
    {
	MPI result = mpi_alloc_like( prime );
	MPI pminus1 = mpi_alloc_like( prime );
	mpi_sub_ui( pminus1, prime, 1);
	mpi_powm( result, val_2, pminus1, prime );
	mpi_free( pminus1 );
	if( mpi_cmp_ui( result, 1 ) ) { /* if composite */
	    mpi_free( result );
	    return 0;
	}
	mpi_free( result );
    }

    /* perform stronger tests */
    if( is_prime(prime, 5, &count ) )
		return 1; /* is probably a prime */
    return 0;
}

/****************
 * We do not need to use the strongest RNG because we gain no extra
 * security from it - The prime number is public and we could also
 * offer the factors for those who are willing to check that it is
 * indeed a strong prime.
 *
 * mode 0: Standard
 *	1: Make sure that at least one factor is of size qbits.
 */
MPI
generate_elg_prime( int mode, unsigned pbits, unsigned qbits,
		    MPI g, MPI **ret_factors )
{
    int n;  /* number of factors */
    int m;  /* number of primes in pool */
    unsigned fbits; /* length of prime factors */
    MPI *factors; /* current factors */
    MPI *pool;	/* pool of primes */
    MPI q;	/* first prime factor (variable)*/
    MPI prime;	/* prime test value */
    MPI q_factor; /* used for mode 1 */
    byte *perms = NULL;
    int i, j;
    int count1, count2;
    unsigned nprime;
    unsigned req_qbits = qbits; /* the requested q bits size */
    MPI val_2  = mpi_alloc_set_ui( 2 );

    /* find number of needed prime factors */
    for(n=1; (pbits - qbits - 1) / n  >= qbits; n++ )
	;
    n--;
   
    if( mode == 1 ) {
		n--;
		fbits = (pbits - 2*req_qbits -1) / n;
		qbits =  pbits - req_qbits - n*fbits;
    }
    else {
		fbits = (pbits - req_qbits -1) / n;
		qbits = pbits - n*fbits;
    }

    prime = mpi_alloc ( mpi_nlimb_hint_from_nbits (pbits) );
    q = gen_prime( qbits, 0, 0 );
    q_factor = mode==1? gen_prime( req_qbits, 0, 0 ) : NULL;

    /* allocate an array to hold the factors + 2 for later usage */
    factors = (MPI *)malloc( (n+2) * sizeof *factors );

    /* make a pool of 3n+5 primes (this is an arbitrary value) */
    m = n*3+5;
    if( mode == 1 )
		m += 5; /* need some more for DSA */
    if( m < 25 )
		m = 25;
    pool = (MPI *)malloc( m * sizeof *pool );

    /* permutate over the pool of primes */
    count1=count2=0;
    do {
      next_try:
	if( !perms ) {
	    /* allocate new primes */
	    for(i=0; i < m; i++ ) {
//			mpi_free(pool[i]);
			pool[i] = NULL;
	    }
	    /* init m_out_of_n() */
	    perms = (byte *)malloc( m );
	    for(i=0; i < n; i++ ) {
			perms[i] = 1;
			pool[i] = gen_prime( fbits, 0, 0 );
			factors[i] = pool[i];
	    }
	}
	else {
	    m_out_of_n( perms, n, m );
	    for(i=j=0; i < m && j < n ; i++ )
		if( perms[i] ) {
		    if( !pool[i] )
		    {
				pool[i] = gen_prime( fbits, 0, 0 );
			}
		    factors[j++] = pool[i];
		}
	    if( i == n ) {
		free(perms); perms = NULL;
		
		goto next_try;	/* allocate new primes */
	    }
	}
	mpi_set( prime, q );
	mpi_mul_ui( prime, prime, 2 );
	if( mode == 1 )
	    mpi_mul( prime, prime, q_factor );
	for(i=0; i < n; i++ )
	    mpi_mul( prime, prime, factors[i] );
	mpi_add_ui( prime, prime, 1 );
	nprime = mpi_get_nbits(prime);
	if( nprime < pbits ) {
	    if( ++count1 > 20 ) {
		count1 = 0;
		qbits++;
		
        mpi_free (q);
		q = gen_prime( qbits, 0, 0 );
		goto next_try;
	    }
	}
	else
	    count1 = 0;
	if( nprime > pbits ) {
	    if( ++count2 > 20 ) {
		count2 = 0;
		qbits--;
		mpi_free (q);
		q = gen_prime( qbits, 0, 0 );
		goto next_try;
	    }
	}
	else
	    count2 = 0;
    } while( !(nprime == pbits && check_prime( prime, val_2 )) );

    
    if( ret_factors ) { /* caller wants the factors */
	*ret_factors = (MPI *)malloc( (n+2) * sizeof **ret_factors);
    i = 0;
	if( mode == 1 ) {
	    (*ret_factors)[i++] = mpi_copy( q_factor );
	    for(; i <= n; i++ )
		(*ret_factors)[i] = mpi_copy( factors[i-1] );
	}
	else {
	    for(; i < n; i++ )
		(*ret_factors)[i] = mpi_copy( factors[i] );
	}
    }

    if( g ) { /* create a generator (start with 3)*/
	MPI tmp   = mpi_alloc( mpi_get_nlimbs(prime) );
	MPI b	  = mpi_alloc( mpi_get_nlimbs(prime) );
	MPI pmin1 = mpi_alloc( mpi_get_nlimbs(prime) );
	
	factors[n] = q;
	factors[n+1] = mpi_alloc_set_ui(2);
	mpi_sub_ui( pmin1, prime, 1 );
	mpi_set_ui(g,2);
	do {
	    mpi_add_ui(g, g, 1);

	    for(i=0; i < n+2; i++ ) {
		/*fputc('~', stderr);*/
		mpi_fdiv_q(tmp, pmin1, factors[i] );
		/* (no mpi_pow(), but it is okay to use this with mod prime) */
		mpi_powm(b, g, tmp, prime );
		if( !mpi_cmp_ui(b, 1) )
		    break;
	    }
	   
	} while( i < n+2 );
	mpi_free(factors[n+1]);
	mpi_free(tmp);
	mpi_free(b);
	mpi_free(pmin1);
    }
    
    free( factors );	/* (factors are shallow copies) */
    for(i=0; i < m; i++ )
		mpi_free( pool[i] );
    free( pool );
    free(perms);
    mpi_free(val_2);
    mpi_free(q);
    return prime;
}

MPI mpi_alloc_secure( unsigned nlimbs )
{
    MPI a;

    a = (MPI)malloc( sizeof *a );
    a->d = nlimbs? mpi_alloc_limb_space( nlimbs, 1 ) : NULL;

    a->alloced = nlimbs;
    a->flags = 1;
    a->nlimbs = 0;
    a->sign = 0;
    a->nbits = 0;
    return a;
}
