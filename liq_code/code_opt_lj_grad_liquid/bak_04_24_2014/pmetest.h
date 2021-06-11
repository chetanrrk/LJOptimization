#ifndef PMETEST_H
#define PMETEST_H

#include "mdtypes.h"

#ifdef __cplusplus
extern "C" {
#endif


  /*
   * just init cell* vectors and cutoff
   * the rest have good defaults - set to 0 to get the default values
   *
   * for now i'm leaving this stuff the way it is
   * i'll just try to get the code to work and clean it up later
   */
  typedef struct PmetestParams_Tag {
    MD_Dvec center;    /* center of cell */
    MD_Dvec cellvec1;  /* provide 3 orientation vectors for the cell */
    MD_Dvec cellvec2;  /*   they define the dimensions of the cell */
    MD_Dvec cellvec3;  /*   they don't have to be orthogonal */
    double cutoff;     /* split distance between direct and recip parts */
    double tolerance;  /* defaults to 1e-6 - good value */
    double ewaldcof;   /* this is init'ed based on PMETolerance and cutoff */
    int32 nxspacings;  /* these each default to 64 - good values */
    int32 nyspacings;  /*   I think these need to be a power of 2 */
    int32 nzspacings;  /*   is it best if they are the same value? */
    int32 interporder; /* defaults to 4 - a good value */
    int32 natoms;      /* number of atoms (expected array length below) */
  } PmetestParams;


  /*
   * user is responsible for all arrays here
   * make sure to init position and charge before calling compute
   */
  typedef struct PmetestSystem_Tag {
    /* output, user supplies arrays */
    double u_elec;
    double u_direct;
    double u_recip;
    double virial_recip[9];
    MD_Dvec *f_elec;
    MD_Dvec *f_direct;
    MD_Dvec *f_recip;
    /* input */
    MD_Dvec *pos;
    double *charge;
  } PmetestSystem;


  /*
   * forward declaration of internal state type
   */
  struct Pmetest_Tag;
  typedef struct Pmetest_Tag Pmetest;


  Pmetest *pmetest_create(PmetestParams *);
  int pmetest_compute(Pmetest *, PmetestSystem *);
  void pmetest_destroy(Pmetest *);

#ifdef __cplusplus
}
#endif

#endif /* PMETEST_H */
