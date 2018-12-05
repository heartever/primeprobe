#ifndef MPIHEADER_H
#define MPIHEADER_H
#include <stdint.h>
#include <sgx_trts.h>
#include <stdio.h>
#include <string.h>
#include "libgcryptEnclave_t.h"
#define BITS_PER_LONG 64
#define BYTES_PER_MPI_LIMB      (BITS_PER_LONG / 8)

typedef unsigned long long mpi_limb_t;
typedef mpi_limb_t *mpi_ptr_t;	/* pointer to a limb */
typedef int mpi_size_t;        /* (must be a signed type) */
typedef   signed long int mpi_limb_signed_t;
typedef unsigned char byte;

#define BITS_PER_MPI_LIMB   64
#define W_TYPE_SIZE BITS_PER_MPI_LIMB
  typedef mpi_limb_t   UWtype;
  typedef unsigned int UHWtype;
  typedef uint64_t     UDWtype;
  
#  define G10_MPI_INLINE_DECL  extern __inline__
#define A_LIMB_1 ((mpi_limb_t)1)

struct gcry_mpi {
    int alloced;    /* array size (# of allocated limbs) */
    int nlimbs;     /* number of valid limbs */
    unsigned int nbits; /* the real number of valid bits (info only) */
    int sign;	    /* indicates a negative number */
    unsigned flags; /* bit 0: array must be allocated in secure memory space */
		    /* bit 1: not used */
		    /* bit 2: the limb is a pointer to some xmalloced data */
    mpi_limb_t *d;  /* array with the limbs */
};
typedef struct gcry_mpi *MPI;


#define __BITS4 (W_TYPE_SIZE / 4)
#define __ll_B ((UWtype) 1 << (W_TYPE_SIZE / 2))
#define __ll_lowpart(t) ((UWtype) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((UWtype) (t) >> (W_TYPE_SIZE / 2))

#define UMUL_TIME 1
#define UDIV_TIME UMUL_TIME

#define KARATSUBA_THRESHOLD 16
struct karatsuba_ctx {
    struct karatsuba_ctx *next;
    mpi_ptr_t tspace;
    mpi_size_t tspace_size;
    mpi_ptr_t tp;
    mpi_size_t tp_size;
};

static int no_of_small_prime_numbers = 668;


void mpi_powm( MPI res, MPI base, MPI exponent, MPI mod);

void mpi_mulm( MPI w, MPI u, MPI v, MPI m);
void mpi_mul( MPI w, MPI u, MPI v);

////////////////////// declaration here  
mpi_ptr_t mpi_alloc_limb_space( unsigned nlimbs, int secure );
MPI mpi_set_opaque( MPI a, void *p, unsigned int len );
void mpi_free( MPI a );
void mpi_add(MPI w, MPI u, MPI v);
int mpihelp_cmp( mpi_ptr_t op1_ptr, mpi_ptr_t op2_ptr, mpi_size_t size );
extern mpi_limb_t mpihelp_mul( mpi_ptr_t prodp, mpi_ptr_t up, mpi_size_t usize, mpi_ptr_t vp, mpi_size_t vsize);
static void mul_n( mpi_ptr_t prodp, mpi_ptr_t up, mpi_ptr_t vp,	mpi_size_t size, mpi_ptr_t tspace )	;		
extern mpi_limb_t mpihelp_addmul_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size, mpi_limb_t s2_limb);
extern mpi_limb_t mpihelp_mul_1( mpi_ptr_t res_ptr, mpi_ptr_t s1_ptr, mpi_size_t s1_size, mpi_limb_t s2_limb);  
void mpi_tdiv_qr( MPI quot, MPI rem, MPI num, MPI den);
void mpi_tdiv_r( MPI rem, MPI num, MPI den);
extern mpi_limb_t mpihelp_lshift( mpi_ptr_t wp, mpi_ptr_t up, mpi_size_t usize, unsigned int cnt);
extern mpi_limb_t mpihelp_rshift( mpi_ptr_t wp, mpi_ptr_t up, mpi_size_t usize, unsigned cnt);
mpi_limb_t mpihelp_mod_1(mpi_ptr_t dividend_ptr, mpi_size_t dividend_size, mpi_limb_t divisor_limb);
mpi_limb_t mpihelp_divrem( mpi_ptr_t qp, mpi_size_t qextra_limbs, mpi_ptr_t np, mpi_size_t nsize, mpi_ptr_t dp, mpi_size_t dsize)	;
byte *get_random_bits( size_t nbits, int level, int secure );
void mpi_set_buffer( MPI a, const byte *buffer, unsigned nbytes, int sign );
void mpih_sqr_n_basecase( mpi_ptr_t prodp, mpi_ptr_t up, mpi_size_t size );
MPI mpi_copy( MPI a );
int mpi_cmp_ui( MPI u, unsigned long v );
MPI generate_elg_prime( int mode, unsigned pbits, unsigned qbits, MPI g, MPI **ret_factors );
void mpi_clear_highbit( MPI a, unsigned n );

unsigned mpi_get_nbits( MPI a );
int mpi_cmp( MPI u, MPI v );
MPI mpi_alloc( unsigned nlimbs );
unsigned int mpi_nlimb_hint_from_nbits (unsigned int nbits);
int mpi_get_nlimbs (MPI a);
int mpi_is_neg (MPI a);
void mpi_sub(MPI w, MPI u, MPI v);
void mpi_fdiv_r( MPI rem, MPI dividend, MPI divisor );
void mpi_sub_ui(MPI w, MPI u, unsigned long v );
void mpi_set( MPI w, MPI u);
void mpi_invm( MPI x, MPI a, MPI n );
void mpi_add_ui(MPI w, MPI u, unsigned long v );
void mpi_fdiv_q( MPI quot, MPI dividend, MPI divisor );
void mpi_set_ui( MPI w, unsigned long u);
int mpi_gcd( MPI g, MPI xa, MPI xb );
void mpi_swap( MPI a, MPI b);
MPI generate_secret_prime( unsigned  nbits );
void mpi_assign_limb_space( MPI a, mpi_ptr_t ap, unsigned nlimbs );

int rsa_generate( int algo, unsigned nbits, MPI *skey, MPI **retfactors );
int rsa_check_secret_key( int algo, MPI *skey );
int rsa_encrypt( int algo, MPI *resarr, MPI data, MPI *pkey );
int rsa_decrypt( int algo, MPI *result, MPI *data, MPI *skey );
int rsa_sign( int algo, MPI *resarr, MPI data, MPI *skey );
int rsa_verify( int algo, MPI hash, MPI *data, MPI *pkey );
unsigned rsa_get_nbits( int algo, MPI *pkey );
const char *rsa_get_info( int algo, int *npkey, int *nskey,
				    int *nenc, int *nsig, int *use );
				    
int elg_generate( int algo, unsigned nbits, MPI *skey, MPI **retfactors );
int elg_check_secret_key( int algo, MPI *skey );
int elg_encrypt( int algo, MPI *resarr, MPI data, MPI *pkey );
int elg_decrypt( int algo, MPI *result, MPI *data, MPI *skey );
int elg_sign( int algo, MPI *resarr, MPI data, MPI *skey );
int elg_verify( int algo, MPI hash, MPI *data, MPI *pkey );
unsigned elg_get_nbits( int algo, MPI *pkey );
const char *elg_get_info( int algo, int *npkey, int *nskey,
				    int *nenc, int *nsig, int *use );
MPI mpi_alloc_secure( unsigned nlimbs );
#define RESIZE_IF_NEEDED(a,b) \
    do {			   \
	if( (a)->alloced < (b) )   \
	    mpi_resize((a), (b));  \
    } while(0)

#define MPN_MUL_N_RECURSE(prodp, up, vp, size, tspace) \
    do {						\
	if( (size) < KARATSUBA_THRESHOLD )		\
	    mul_n_basecase (prodp, up, vp, size);	\
	else						\
	    mul_n (prodp, up, vp, size, tspace);	\
    } while (0);

#define MPN_SQR_N_RECURSE(prodp, up, size, tspace) \
    do {					    \
	if ((size) < KARATSUBA_THRESHOLD)	    \
	    mpih_sqr_n_basecase (prodp, up, size);	 \
	else					    \
	    mpih_sqr_n (prodp, up, size, tspace);	 \
    } while (0);

#define MPN_NORMALIZE(d, n)  \
    do {		       \
	while( (n) > 0 ) {     \
	    if( (d)[(n)-1] ) \
		break;	       \
	    (n)--;	       \
	}		       \
    } while(0)

/* Zero N limbs at D */
#define MPN_ZERO(d, n) \
    do {				  \
	int  _i;			  \
	for( _i = 0; _i < (n); _i++ )  \
	    (d)[_i] = 0;		    \
    } while (0)
    
/* Copy N limbs from S to D.  */
#define MPN_COPY( d, s, n) \
    do {				\
	mpi_size_t _i;			\
	for( _i = 0; _i < (n); _i++ )	\
	    (d)[_i] = (s)[_i];		\
    } while(0)
   
   
#define count_leading_zeros(count, x) \
  do {									\
    UWtype __xr = (x);							\
    UWtype __a; 							\
							\
    if (W_TYPE_SIZE <= 32)						\
      { 								\
		__a = __xr < ((UWtype) 1 << 2*__BITS4)				\
		  ? (__xr < ((UWtype) 1 << __BITS4) ? 0 : __BITS4)		\
		  : (__xr < ((UWtype) 1 << 3*__BITS4) ?  2*__BITS4 : 3*__BITS4);\
      } 								\
		else								\
		  { 								\
		for (__a = W_TYPE_SIZE - 8; __a > 0; __a -= 8)			\
		  if (((__xr >> __a) & 0xff) != 0)				\
			break;							\
      } 								\
								\
    (count) = W_TYPE_SIZE - (__clz_tab[__xr >> __a] + __a);		\
  } while (0)

/////////////////////////////////
#define UDIV_NEEDS_NORMALIZATION 1
#define udiv_qrnnd __udiv_qrnnd_c
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

#define MPN_COPY_DECR( d, s, n ) \
    do {				\
	mpi_size_t _i;			\
	for( _i = (n)-1; _i >= 0; _i--) \
	   (d)[_i] = (s)[_i];		\
    } while(0)
    
#define sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									\
    UWtype __x; 							\
    __x = (al) - (bl);							\
    (sh) = (ah) - (bh) - (__x > (al));					\
    (sl) = __x; 							\
  } while (0)

#define add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									\
    UWtype __x; 							\
    __x = (al) + (bl);							\
    (sh) = (ah) + (bh) + (__x < (al));					\
    (sl) = __x; 							\
  } while (0)
 	
 	
#define umul_ppmm(w1, w0, u, v) 					\
  do {									\
    UWtype __x0, __x1, __x2, __x3;					\
    UHWtype __ul, __vl, __uh, __vh;					\
    UWtype __u = (u), __v = (v);					\
									\
    __ul = __ll_lowpart (__u);						\
    __uh = __ll_highpart (__u); 					\
    __vl = __ll_lowpart (__v);						\
    __vh = __ll_highpart (__v); 					\
									\
    __x0 = (UWtype) __ul * __vl;					\
    __x1 = (UWtype) __ul * __vh;					\
    __x2 = (UWtype) __uh * __vl;					\
    __x3 = (UWtype) __uh * __vh;					\
									\
    __x1 += __ll_highpart (__x0);/* this can't give carry */            \
    __x1 += __x2;		/* but this indeed can */		\
    if (__x1 < __x2)		/* did we get it? */			\
      __x3 += __ll_B;		/* yes, add it in the proper pos. */	\
									\
    (w1) = __x3 + __ll_highpart (__x1); 				\
    (w0) = (__ll_lowpart (__x1) << W_TYPE_SIZE/2) + __ll_lowpart (__x0);\
  } while (0)
    
#endif
