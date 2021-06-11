/* modified from NAMD */

#ifndef COMPUTEPME_H
#define COMPUTEPME_H

#include "pmetest.h"
#include "PmeBase.h"
#include "Vector.h"

class PmeRealSpace;
class ComputePmeMgr;


class ComputePme { //: public ComputeHomePatches {
public:
  ComputePme(const PmetestParams &); //ComputeID c);
  virtual ~ComputePme();
  void doWork(int natoms, const Vector *position, const double *charge,
              Vector *force, double *energy, double *virial);
            // energy is length 1, virial is length 9
  void ungridForces();

  PmetestParams *simParams;
 private:
  ComputePmeMgr *computePmeMgr;
  PmeGrid myGrid;
  int qsize, fsize, bsize;
public:
  double **q_arr;
private:
  char *f_arr;
  char *fz_arr;
public:
  double energy;
  double virial[6];
private:
  int resultsRemaining;
  PmeRealSpace *myRealSpace;
  int numLocalAtoms;
  PmeParticle *localData;

  Vector *pmeForce;
  double *pmeVirial;
  PmetestParams storePmeParams;

  int useAvgPositions;
};

#endif

