/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include "PmeRealSpace.h"

PmeRealSpace::PmeRealSpace(PmeGrid grid,int natoms)
  : myGrid(grid), N(natoms) {
  int order = myGrid.order;
  M = new double[3*N*order];
  dM = new double[3*N*order];
}

PmeRealSpace::~PmeRealSpace() {
  delete [] M;
  delete [] dM;
}


void PmeRealSpace::fill_b_spline(PmeParticle p[]) {
  double fr[3]; 
  double *Mi, *dMi;
  int i, stride;
  int order = myGrid.order;

  stride = 3*order;
  Mi = M; dMi = dM;
  for (i=0; i<N; i++) {
    fr[0] = p[i].x; fr[0] -= (int)(fr[0]);
    fr[1] = p[i].y; fr[1] -= (int)(fr[1]);
    fr[2] = p[i].z; fr[2] -= (int)(fr[2]);
    compute_b_spline(fr, Mi, dMi, order);
    Mi += stride;
    dMi += stride;
  }
}

void PmeRealSpace::fill_charges(double **q_arr, char *f_arr, char *fz_arr, PmeParticle p[]) {
  
  int i, j, k, l;
  int stride;
  int K1, K2, K3, dim2, dim3, order;
  double *Mi;
  Mi = M;
  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2; dim3=myGrid.dim3;
  order = myGrid.order;
  stride = 3*order;

  fill_b_spline(p);

  for (i=0; i<N; i++) {
    double q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    for (j=0; j<order; j++) {
      double m1;
      int ind1;
      m1 = Mi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2;
      u2 = u2i;
      for (k=0; k<order; k++) {
        double m1m2;
	int ind2;
        m1m2 = m1*Mi[order+k];
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	double *qline = q_arr[ind2];
	if ( ! qline ) {
	  q_arr[ind2] = qline = new double[dim3];
	  memset( (void*) qline, 0, dim3 * sizeof(double) );
	}
	f_arr[ind2] = 1;
        for (l=0; l<order; l++) {
	  double m3;
	  int ind;
	  m3 = Mi[2*order + l];
	  int u3 = u3i + l;
          ind = u3 + (u3 < 0 ? K3 : 0);
          qline[ind] += m1m2*m3; 
        }
      }
    }
    Mi += stride;
    for (l=0; l<order; l++) {
      int u3 = u3i + l;
      int ind = u3 + (u3 < 0 ? K3 : 0);
      fz_arr[ind] = 1;
    }
  }
}

void PmeRealSpace::compute_forces(const double * const *q_arr,
				const PmeParticle p[], Vector f[]) {
  
  int i, j, k, l, stride;
  double f1, f2, f3;
  double *Mi, *dMi;
  int K1, K2, K3, dim2, order;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3; dim2=myGrid.dim2;
  order = myGrid.order;
  stride=3*order;
  Mi = M; dMi = dM;
 
  for (i=0; i<N; i++) {
    double q;
    int u1, u2, u2i, u3i;
    q = p[i].cg;
    f1=f2=f3=0.0;
    u1 = (int)(p[i].x);
    u2i = (int)(p[i].y);
    u3i = (int)(p[i].z);
    u1 -= order;
    u2i -= order;
    u3i -= order;
    u3i += 1;
    for (j=0; j<order; j++) {
      double m1, d1;
      int ind1;
      m1=Mi[j]*q;
      d1=K1*dMi[j]*q;
      u1++;
      ind1 = (u1 + (u1 < 0 ? K1 : 0))*dim2; 
      u2 = u2i;
      for (k=0; k<order; k++) {
        double m2, d2, m1m2, m1d2, d1m2;
	int ind2;
        m2=Mi[order+k];
	d2=K2*dMi[order+k];
	m1m2=m1*m2;
	m1d2=m1*d2;
	d1m2=d1*m2;
	u2++;
	ind2 = ind1 + (u2 + (u2 < 0 ? K2 : 0));
	const double *qline = q_arr[ind2];
        for (l=0; l<order; l++) {
	  double term, m3, d3;
	  int ind;
	  m3=Mi[2*order+l];
	  d3=K3*dMi[2*order+l];
	  int u3 = u3i + l;
	  ind = u3 + (u3 < 0 ? K3 : 0);
	  term = qline[ind];
	  f1 -= d1m2 * m3 * term;
	  f2 -= m1d2 * m3 * term;
	  f3 -= m1m2 * d3 * term;
        }
      }
    }
    Mi += stride;
    dMi += stride;
    f[i].x = f1;
    f[i].y = f2;
    f[i].z = f3;
  }
}
   
