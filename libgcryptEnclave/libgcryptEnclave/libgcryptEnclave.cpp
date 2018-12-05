#include "libgcryptEnclave_t.h"
#include "sgx_trts.h"

#include "mpiheader.h"

using namespace std;

int nbits = 1024;
MPI skey[4]; MPI *retfactors;

MPI test = mpi_alloc( 0 );
MPI out1[2];

MPI pk[3];
extern bool flagprint;

int FM(uint32_t input)
{
	elg_generate(0, nbits, skey, &retfactors);

	for(int i = 0; i < skey[3]->nlimbs; i++)   // x
	{
		TRACE(14, skey[3]->d[i]);
	}
	MYDEBUG(skey[3]->nlimbs);                // 4
	MYDEBUG(elg_check_secret_key(1, skey));  // 0
	TRACE(14, (uint64_t)&mpih_sqr_n_basecase);
	TRACE(14, (uint64_t)&elg_decrypt);
	
    pk[0] = skey[0];
    pk[1] = skey[1];
    pk[2] = skey[2];
	{	
    	unsigned char *p = get_random_bits( nbits, 0, 0 );
		mpi_set_buffer( test, p, (nbits+7)/8, 0 );
		free(p);
    }	
	
    elg_encrypt(1, out1, test, pk);
    return 1;
}

/*int testvar[1000] = {10};

void maccess2(void *p)
{
//	asm volatile ("clflush 0(%0)\n" :: "r" (p) : "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
	asm volatile ("movq (%0), %%rax" :: "r" (p) : "memory", "rax");
	asm volatile ("mfence");
}*/

void decrypt()
{
//	while(1)
	{  
		MPI out2;
		int testvar2 = 0;
//		flagprint = true;
	//for(int i = 0; i < 1; i++)
	//{
	//	TRACE(14, (uint64_t)&testvar[500]);
	//	maccess2(&testvar[500]);
	//for(int i = 0; i < 20; i++)
		elg_decrypt(1, &out2, out1, skey);
	//}
//	MYDEBUG(1000);
//		for(int i = 0; i < out2->nlimbs; i++)   // secret exponent
	//	{
//			TRACE(14, out2->d[i]);
//		}
//		mpi_free(out2);
//		enclsleep(1000000);
//		MYDEBUG(out2->nlimbs);
	}
    
	
}
