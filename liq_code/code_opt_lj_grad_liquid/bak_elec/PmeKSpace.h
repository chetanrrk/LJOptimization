/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PME_K_SPACE_H__
#define PME_K_SPACE_H__

#include "PmeBase.h"
#include "lattice.h"

class PmeKSpace {

public:
  PmeKSpace(PmeGrid grid, int K2_start, int K2_end);
  ~PmeKSpace();

  double compute_energy(double q_arr[], Lattice lattice, double ewald,
                        double virial[]);
  
private:
  // b-spline moduli
  double *bm1, *bm2, *bm3; 
  double *exp1, *exp2, *exp3;
  double i_pi_volume, piob;

  const PmeGrid myGrid;
  const int k2_start, k2_end;

  void init_exp(double *xp, int K, double recip);
};

#endif
