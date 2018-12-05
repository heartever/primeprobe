#include <string.h>

#include "analyse.h"



int match1(char *buf, char *pattern) {
  int w = 0;
  while (*pattern) {
    w += (*pattern == *buf);
    pattern++;
    buf++;
  }
  return w;
}



int match(char *buf, int len, char *pattern) {
  int patternlen = strlen(pattern);

  int rv = 0;
  for (int i = 0; i < len; i++)
    buf[i] = (buf[i]==0) ? '0' : '1';

    /*finding a pattern in a flow, return number of matching possibilities*/
  for (int i = 0; i < len - patternlen; i++) {
    int w = match1(buf + i, pattern);
    if (w >= patternlen-3) 
      rv++;
  }
  return rv;
}
