/* this is a C-callable wrapper for the C++ pme library */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pmetest.h"
#include "ComputePme.h"


/* default values */
#define DEFAULT_PME_TOLERANCE  1e-6
#define DEFAULT_PME_GRIDSIZEX  64
#define DEFAULT_PME_GRIDSIZEY  64
#define DEFAULT_PME_GRIDSIZEZ  64
#define DEFAULT_PME_INTERPORDER 4

extern double erfc(double x);
extern double pmeDirect(int natoms, Vector f[], Vector pos[], const double q[], const PmetestParams *prm);

/* define Pme internal data structure */
struct Pmetest_Tag {
  PmetestParams savePmeParams;
  ComputePme *computePme;
  Vector *position;
  Vector *direct_force;
  Vector *recip_force;
  int natoms;
};


Pmetest *pmetest_create(PmetestParams *params)
{
  double cutoff, tolerance, ewaldcof, ewaldcof_lo, ewaldcof_hi;
  int i;
  Pmetest *pme;

  /* set default values */
  if (params->tolerance == 0) {
    params->tolerance = DEFAULT_PME_TOLERANCE;
  }
  if (params->nxspacings == 0) {
    params->nxspacings = DEFAULT_PME_GRIDSIZEX;
  }
  if (params->nyspacings == 0) {
    params->nyspacings = DEFAULT_PME_GRIDSIZEY;
  }
  if (params->nzspacings == 0) {
    params->nzspacings = DEFAULT_PME_GRIDSIZEZ;
  }
  if (params->interporder == 0) {
    params->interporder = DEFAULT_PME_INTERPORDER;
  }

  /* initialize ewaldcof based on PMETolerance and cutoff */
  cutoff = params->cutoff;
  tolerance = params->tolerance;
  ewaldcof = 1.0;
  while ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) ewaldcof *= 2.0;
  ewaldcof_lo = 0.;
  ewaldcof_hi = ewaldcof;
  for (i = 0; i < 100; ++i ) {
    ewaldcof = 0.5 * ( ewaldcof_lo + ewaldcof_hi );
    if ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) {
      ewaldcof_lo = ewaldcof;
    }
    else {
      ewaldcof_hi = ewaldcof;
    }
  }
  params->ewaldcof = ewaldcof;

  if ((pme = (Pmetest *) calloc(1, sizeof(Pmetest))) == NULL) {
    return NULL;
  }
  pme->savePmeParams = *params;
  pme->computePme = new ComputePme(pme->savePmeParams);
  return pme;
}


int pmetest_compute(Pmetest *pme, PmetestSystem *sys)
{
  int i;

//  TEXT("initialization");

  /* see if we need resize or initialize our arrays */
  if (pme->natoms != pme->savePmeParams.natoms) {
    int natoms = pme->savePmeParams.natoms;
    void *tmp;
    if ((tmp = realloc(pme->position, natoms * sizeof(Vector))) == NULL) {
      return -1;
    }
    pme->position = (Vector *) tmp;
    if ((tmp = realloc(pme->direct_force, natoms * sizeof(Vector))) == NULL) {
      return -1;
    }
    pme->direct_force = (Vector *) tmp;
    if ((tmp = realloc(pme->recip_force, natoms * sizeof(Vector))) == NULL) {
      return -1;
    }
    pme->recip_force = (Vector *) tmp;
    pme->natoms = natoms;
  }

  /* init position vector */
  for (i = 0;  i < pme->natoms;  i++) {
    pme->position[i] = sys->pos[i];
  }

  /* zero force arrays, virial array, and potentials */
  memset(pme->direct_force, 0, pme->natoms * sizeof(Vector));
  memset(pme->recip_force, 0, pme->natoms * sizeof(Vector));
  memset(sys->virial_recip, 0, 9 * sizeof(double));
  sys->u_direct = 0;
  sys->u_recip = 0;

#if 0
  /* compute the real space part */
//  TEXT("compute real space part");
  sys->u_direct = pmeDirect(pme->natoms, pme->direct_force,
      pme->position, sys->charge, &(pme->savePmeParams));

  memset(pme->direct_force, 0, pme->natoms * sizeof(Vector));
#endif

  /* compute the reciprocal space part */
//  TEXT("compute reciprocal space part");
  pme->computePme->doWork(pme->natoms, pme->position, sys->charge,
      pme->recip_force, &(sys->u_recip), sys->virial_recip);

  /* add up the results */
//  TEXT("accumulate results");
  if (sys->f_direct != NULL) {
    for (i = 0;  i < pme->natoms;  i++) {
      sys->f_direct[i].x = pme->direct_force[i].x;
      sys->f_direct[i].y = pme->direct_force[i].y;
      sys->f_direct[i].z = pme->direct_force[i].z;
    }
  }
  if (sys->f_recip != NULL) {
    for (i = 0;  i < pme->natoms;  i++) {
      sys->f_recip[i].x = pme->recip_force[i].x;
      sys->f_recip[i].y = pme->recip_force[i].y;
      sys->f_recip[i].z = pme->recip_force[i].z;
    }
  }
  for (i = 0;  i < pme->natoms;  i++) {
    sys->f_elec[i].x = pme->direct_force[i].x + pme->recip_force[i].x;
    sys->f_elec[i].y = pme->direct_force[i].y + pme->recip_force[i].y;
    sys->f_elec[i].z = pme->direct_force[i].z + pme->recip_force[i].z;
  }
  sys->u_elec = sys->u_direct + sys->u_recip;

//  TEXT("all done");
  return 0;
}


void pmetest_destroy(Pmetest *pme)
{
  delete pme->computePme;
  free(pme->position);
  free(pme->direct_force);
  free(pme->recip_force);
  free(pme);
}
