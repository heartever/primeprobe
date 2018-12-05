// ReadMapping.cpp : Defines the entry point for the console application.
//
#include "sgx_urts.h"
#include "libgcryptEnclave_u.h"
#include <stdio.h>
#include <unistd.h>
#include <sched.h>
#include <time.h>
#define ONE_OVER_BILLION 1E-9

#if defined(_MSC_VER)
# define ENCLAVE_FILENAME "libgcryptEnclave.signed.dll"
#elif defined(__GNUC__)
# define ENCLAVE_FILENAME "libgcryptEnclave.signed.so"
#endif

void MYDEBUG(int out)
{
	printf("Running here: %d.\n", out);
}

void TRACE(uint32_t type, uint64_t addr)
{
	switch(type)
	{
	case 0:
		printf("mpih_sqr_n_basecase address: %016llx.\n", addr); break;
	case 1:
		printf("suffixes[0] address: %016llx.\n", addr); break;
	case 2:
		printf("suffixes[sampled->rank1(j)-1] address: %016llx.\n", addr); break;
	case 3:
		printf("loadFromBuffer() function address: %016llx.\n", addr); break;
	case 4:
		printf("C[c] address: %016llx.\n", addr); break;
	case 5:
		printf("remap[] address: %016llx.\n", addr); break;
	case 6:
		printf("nw() address: %016llx.\n", addr); break;
	case 7:
		printf("dpm_init() address: %016llx.\n", addr); break;
	case 8:
		printf("nw_align() address: %016llx.\n", addr); break;
	case 9:
		printf("F[][]: %016llx.\n", addr); break;
	case 10:
		printf("T_bwt: %016llx.\n", addr); break;
	case 11:
		printf("T_bwt->bitstring[..]: %016llx.\n", addr); break;
	case 12:
		printf("T_bwt->occ: %016llx.\n", addr); break;
	case 13:
		printf("%016llx.\n", addr); break;
	case 14:
		printf("-----------------  %016llx.\n", addr); break;
	}
}

void enclsleep(int usecond)
{
	usleep(usecond);
}

int main(int argc, char* argv[])
{
	/* enclave setup. */
	sgx_enclave_id_t eid;
	sgx_status_t ret = SGX_ERROR_UNEXPECTED;
	sgx_launch_token_t token = {0};
	int updated = 0;

	ret = sgx_create_enclave(ENCLAVE_FILENAME, 1, &token, &updated, &eid, NULL);  // 1 for debug

	if(ret != SGX_SUCCESS)
	{
		printf("Create Enclave Failed. Error code: 0x%02x.\n", ret);
		return -1;
	}
	/* enclave setup finished. */
	
	/* pin the process to core 2 */
	cpu_set_t mask;
  	CPU_ZERO(&mask);
  	CPU_SET(1, &mask);
  
	int result; int pn = 0;
	FM(eid, &result, pn);
	
	uint32_t count = 0;
	struct timespec begin, end;
	
//	decrypt(eid);
	
	clock_gettime(CLOCK_REALTIME, &begin);
	while(1)
	{
		//char c = getchar();
		//if(c == '1')
		decrypt(eid);
		
//		printf("%d.\n", count++);
//		usleep(1000000);
	}
	clock_gettime(CLOCK_REALTIME, &end);
	// Calculate time it took
	double accum = ( end.tv_sec - begin.tv_sec ) + ( end.tv_nsec - begin.tv_nsec ) * ONE_OVER_BILLION;
	printf("Time: %lf\n", accum);
	return 0;
}

