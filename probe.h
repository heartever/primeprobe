#ifndef __PROBE_H__
#define __PROBE_H__


#define PAGE_BITS	12
#define PAGE_SIZE	(1<<PAGE_BITS)
#define PAGE_MASK	(PAGE_SIZE - 1)

#define SHMKEY 1000

typedef struct listEntry *listEntry_t;


int probe(pageset_t pages, int cand, int num);
uint64_t probe_npages();
void probe_init(void * buf, uint64_t size);
int probe2(pageset_t pages, int cand, int num);

/*
 * Creates a listEntry on cache lines of pages that may conflict with adrs.
 */
listEntry_t probe_makelist(pageset_t pages, void *adrs);

int probe_timelist(listEntry_t le);
listEntry_t doPrime(listEntry_t le);
void prime(listEntry_t le);
int probe_time(void *p);


void randline();


#endif //__PROBE_H__
