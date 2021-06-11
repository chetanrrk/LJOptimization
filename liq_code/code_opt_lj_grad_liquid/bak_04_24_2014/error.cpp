#include <stdio.h>
#include <stdlib.h>

void NAMD_bug(const char *errmsg)
{
  fprintf(stderr, "%s\n", errmsg);
  exit(1);
}

void NAMD_die(const char *errmsg)
{
  fprintf(stderr, "%s\n", errmsg);
  exit(1);
}
