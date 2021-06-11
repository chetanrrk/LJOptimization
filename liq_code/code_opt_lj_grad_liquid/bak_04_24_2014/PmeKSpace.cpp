/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h>
#include <math.h>
#include "PmeKSpace.h"

static void dftmod(double *bsp_mod, double *bsp_arr, int nfft) {
  int j, k;
  double twopi, arg, sum1, sum2;
  double infft = 1.0/nfft;
/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
  twopi =  2.0 * M_PI;

  for (k = 0; k <nfft; ++k) {
    sum1 = 0.;
    sum2 = 0.;
    for (j = 0; j < nfft; ++j) {
      arg = twopi * k * j * infft;
      sum1 += bsp_arr[j] * cos(arg);
      sum2 += bsp_arr[j] * sin(arg);
    }
    bsp_mod[k] = sum1*sum1 + sum2*sum2;
  }
}

static void compute_b_moduli(double *bm, int K, int order) {
  int i;
  double fr[3];

  double *M = new double[3*order];
  double *dM = new double[3*order];
  double *scratch = new double[K];

  fr[0]=fr[1]=fr[2]=0.0;
  compute_b_spline(fr,M,dM,order);  
  for (i=0; i<order; i++) bm[i] = M[i];
  for (i=order; i<K; i++) bm[i] = 0.0;
  dftmod(scratch, bm, K);
  for (i=0; i<K; i++) bm[i] = 1.0/scratch[i];

  delete [] scratch;
  delete [] dM;
  delete [] M;
}

PmeKSpace::PmeKSpace(PmeGrid grid, int K2_start, int K2_end) 
  : myGrid(grid), k2_start(K2_start), k2_end(K2_end) {
  int K1, K2, K3, order;
  K1=myGrid.K1; K2=myGrid.K2, K3=myGrid.K3; order=myGrid.order;

  bm1 = new double[K1];
  bm2 = new double[K2];
  bm3 = new double[K3];

  exp1 = new double[K1/2 + 1];
  exp2 = new double[K2/2 + 1];
  exp3 = new double[K3/2 + 1];

  compute_b_moduli(bm1, K1, order);
  compute_b_moduli(bm2, K2, order);
  compute_b_moduli(bm3, K3, order);
}

PmeKSpace::~PmeKSpace() {
  delete [] bm1;
  delete [] bm2;
  delete [] bm3;
  
  delete [] exp1;
  delete [] exp2;
  delete [] exp3;
}

double PmeKSpace::compute_energy(double *q_arr, Lattice lattice, double ewald, double *virial) {
  double energy = 0.0;

  int pad2, pad3, n;
  int k1, k2, k3, ind;
  int K1, K2, K3;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3;
  pad2 = 0;
  pad3 = myGrid.dim3-K3-(K3 & 1 ? 1 : 2);

  i_pi_volume = 1.0/(M_PI * lattice_volume(&lattice));
  piob = M_PI/ewald;
  piob *= piob;

  for (n=0; n<6; virial[n++] = 0.0);

  if ( lattice_is_orthogonal(&lattice) ) {
    double recipx = lattice.b1.x;
    double recipy = lattice.b2.y;
    double recipz = lattice.b3.z;
    init_exp(exp1, K1, recipx);
    init_exp(exp2, K2, recipy);
    init_exp(exp3, K3, recipz);

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double m1, m11, b1, xp1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recipx;
      m11 = m1*m1;
      xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double m2, m22, b1b2, xp2;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2 = k2_s*recipy;
        m22 = m2*m2;
        xp2 = exp2[abs(k2_s)]*xp1;
        if ( k1==0 && k2==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          double m3, m33, xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          theta3 = bm3[k3] *b1b2;
          m3 = k3*recipz;
          m33 = m3*m3;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          msq = m11 + m22 + m33;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          virial[0] += fac*(1.0+vir*m11);
          virial[1] += fac*vir*m1*m2;
          virial[2] += fac*vir*m1*m3;
          virial[3] += fac*(1.0+vir*m22);
          virial[4] += fac*vir*m2*m3;
          virial[5] += fac*(1.0+vir*m33);
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  } else if ( cross(lattice.a1,lattice.a2).unit() == lattice.a3.unit() ) {

    Vector recip1 = lattice.b1;
    Vector recip2 = lattice.b2;
    Vector recip3 = lattice.b3;
    double recip3_x = recip3.x;
    double recip3_y = recip3.y;
    double recip3_z = recip3.z;
    init_exp(exp3, K3, recip3.length());

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double b1; Vector m1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recip1;
      // xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double xp2, b1b2, m2_x, m2_y, m2_z;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2_x = m1.x + k2_s*recip2.x;
        m2_y = m1.y + k2_s*recip2.y;
        m2_z = m1.z + k2_s*recip2.z;
        // xp2 = exp2[abs(k2_s)]*xp1;
        xp2 = i_pi_volume*exp(-piob*(m2_x*m2_x+m2_y*m2_y+m2_z*m2_z));
        if ( k1==0 && k2==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  } else {

    Vector recip1 = lattice.b1;
    Vector recip2 = lattice.b2;
    Vector recip3 = lattice.b3;
    double recip3_x = recip3.x;
    double recip3_y = recip3.y;
    double recip3_z = recip3.z;

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double b1; Vector m1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recip1;
      // xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double b1b2, m2_x, m2_y, m2_z;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2_x = m1.x + k2_s*recip2.x;
        m2_y = m1.y + k2_s*recip2.y;
        m2_z = m1.z + k2_s*recip2.z;
        // xp2 = exp2[abs(k2_s)]*xp1;
        if ( k1==0 && k2==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          // xp3 = exp3[k3];
          xp3 = i_pi_volume*exp(-piob*msq);
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  }

  for (n=0; n<6; ++n) virial[n] *= 0.5;
  return 0.5*energy;
}


void PmeKSpace::init_exp(double *xp, int K, double recip) {
  int i;
  double fac;
  fac = -piob*recip*recip;
  for (i=0; i<= K/2; i++)
    xp[i] = exp(fac*i*i);
} 
