/* modified from NAMD */

/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include <stdio.h>

#include "fftw.h"
#include "rfftw.h"

#include "ComputePme.h"
#include "PmeRealSpace.h"
#include "PmeKSpace.h"

#ifndef SQRT_PI
#define SQRT_PI 1.7724538509055160273 /* mathematica 15 digits*/
#endif


int CkNumPes(void) { return 1; }  // only have 1 processor
int CkMyPe(void) { return 0; }    // always have processor id 0



struct LocalPmeInfo {
	int nx, x_start;
	int ny_after_transpose, y_start_after_transpose;
};


class ComputePmeMgr { //: public BOCclass {
public:
	ComputePmeMgr();
	~ComputePmeMgr();
	
	void initialize(); //CkQdMsg*);
	void gridCalc1(void);
	void gridCalc2(void);
	void gridCalc3(void);
	void ungridCalc(void);
	
	void setCompute(ComputePme *c) { pmeCompute = c; }
	
private:
	ComputePme *pmeCompute;
	PmeGrid myGrid;
	Lattice lattice;
	PmeKSpace *myKSpace;
	double *qgrid;
	
	fftw_plan forward_plan_x, backward_plan_x;
	rfftwnd_plan forward_plan_yz, backward_plan_yz;
	fftw_complex *work;
	
	LocalPmeInfo *localInfo;
	int qgrid_start;
	int qgrid_len;
	int fgrid_start;
	int fgrid_len;
	
	int numSources;
	int numGridPes;
	int numTransPes;
	int numRecipPes;
	int numDestRecipPes;
	int firstDestRecipPe;
	int myRecipPe;
	int *recipPeMap;
	int *recipPeDest;
	int grid_count;
	int trans_count;
	int untrans_count;
	int ungrid_count;
	int trans_buf_len;
	int untrans_buf_len;
	double recipEnergy;
	double recip_vir[6];
	double recipEnergy2;
	double recip_vir2[6];
};

ComputePmeMgr::ComputePmeMgr() : /*pmeProxy(thisgroup),*/ pmeCompute(0) {
	myKSpace = 0;
	localInfo = new LocalPmeInfo[CkNumPes()];
	recipPeMap = new int[CkNumPes()];
	recipPeDest = new int[CkNumPes()];
	qgrid = 0;
	work = 0;
	grid_count = 0;
	trans_count = 0;
	untrans_count = 0;
	ungrid_count = 0;
	trans_buf_len = 0;
	untrans_buf_len = 0;
}

void ComputePmeMgr::initialize() { //CkQdMsg *msg) {
	PmetestParams *simParams = pmeCompute->simParams;
	
	{  // decide how many pes to use for reciprocal sum
		int nrp = 1;
		
		// rules based on work available
		int minslices = 1;
		int dimx = simParams->nxspacings;
		int nrpx = ( dimx + minslices - 1 ) / minslices;
		if ( nrpx > nrp ) nrp = nrpx;
		int dimy = simParams->nyspacings;
		int nrpy = ( dimy + minslices - 1 ) / minslices;
		if ( nrpy > nrp ) nrp = nrpy;
		
		// rules based on processors available
		int nrpp = CkNumPes();
		// if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
		if ( nrpp < nrp ) nrp = nrpp;
		
		// user override
		//    int nrps = simParams->PMEProcessors;
		int nrps = CkNumPes();
		if ( nrps > CkNumPes() ) nrps = CkNumPes();
		if ( nrps > 0 ) nrp = nrps;
		
		// make sure there aren't any totally empty processors
		int bx = ( dimx + nrp - 1 ) / nrp;
		int nrpbx = ( dimx + bx - 1 ) / bx;
		int by = ( dimy + nrp - 1 ) / nrp;
		int nrpby = ( dimy + by - 1 ) / by;
		nrp = ( nrpby > nrpbx ? nrpby : nrpbx );
		if ( bx != ( dimx + nrp - 1 ) / nrp )
			NAMD_bug("Error in selecting number of PME processors.");
		if ( by != ( dimy + nrp - 1 ) / nrp )
			NAMD_bug("Error in selecting number of PME processors.");
		
		numGridPes = nrpbx;
		numTransPes = nrpby;
		numRecipPes = nrp;
	}
	
	myRecipPe = 0;
	
	myGrid.K1 = simParams->nxspacings;
	myGrid.K2 = simParams->nyspacings;
	myGrid.K3 = simParams->nzspacings;
	myGrid.order = simParams->interporder;
	myGrid.dim2 = myGrid.K2;
	myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
	myGrid.block1 = ( myGrid.K1 + numRecipPes - 1 ) / numRecipPes;
	myGrid.block2 = ( myGrid.K2 + numRecipPes - 1 ) / numRecipPes;
	
	int nx = 0;
	int ny = 0;
	for ( int pe = 0; pe < numRecipPes; ++pe ) {
		localInfo[pe].x_start = nx;
		nx += myGrid.block1;
		if ( nx > myGrid.K1 ) nx = myGrid.K1;
		localInfo[pe].nx = nx - localInfo[pe].x_start;
		localInfo[pe].y_start_after_transpose = ny;
		ny += myGrid.block2;
		if ( ny > myGrid.K2 ) ny = myGrid.K2;
		localInfo[pe].ny_after_transpose =
			ny - localInfo[pe].y_start_after_transpose;
	}
	
	lattice_set(&lattice,
		simParams->cellvec1, simParams->cellvec2,
		simParams->cellvec3, simParams->center);
	
	
	ungrid_count = numDestRecipPes;
	
	if ( myRecipPe < 0 ) return;
	// the following only for nodes doing reciprocal sum
	
	int k2_start = localInfo[myRecipPe].y_start_after_transpose;
	int k2_end = k2_start + localInfo[myRecipPe].ny_after_transpose;
	myKSpace = new PmeKSpace(myGrid, k2_start, k2_end);
	
	int local_size = myGrid.block1 * myGrid.K2 * myGrid.dim3;
	int local_size_2 = myGrid.block2 * myGrid.K1 * myGrid.dim3;
	if ( local_size < local_size_2 ) local_size = local_size_2;
	qgrid = new double[local_size];
	
	qgrid_start = localInfo[myRecipPe].x_start * myGrid.K2 * myGrid.dim3;
	qgrid_len = localInfo[myRecipPe].nx * myGrid.K2 * myGrid.dim3;
	fgrid_start = localInfo[myRecipPe].x_start * myGrid.K2;
	fgrid_len = localInfo[myRecipPe].nx * myGrid.K2;
	
	int n[3]; n[0] = myGrid.K1; n[1] = myGrid.K2; n[2] = myGrid.K3;
	
	work = new fftw_complex[n[0]];
	
	/*
	*  see if using FFTW_ESTIMATE makes start up faster
	*/
	
	//  if ( ! CkMyPe() ) iout << iINFO << "Optimizing 4 FFT steps.  1..." << endi;
//	printf("#!! Optimizing 4 FFT steps.  1...");  fflush(stdout);
	forward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_REAL_TO_COMPLEX,
		FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
	//  if ( ! CkMyPe() ) iout << " 2..." << endi;
//	printf(" 2...");  fflush(stdout);
	forward_plan_x = fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
		FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) qgrid,
		localInfo[myRecipPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
	//  if ( ! CkMyPe() ) iout << " 3..." << endi;
//	printf(" 3...");  fflush(stdout);
	backward_plan_x = fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
		FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) qgrid,
		localInfo[myRecipPe].ny_after_transpose * myGrid.dim3 / 2, work, 1);
	//  if ( ! CkMyPe() ) iout << " 4..." << endi;
//	printf(" 4...");  fflush(stdout);
	backward_plan_yz = rfftwnd_create_plan_specific(2, n+1, FFTW_COMPLEX_TO_REAL,
		FFTW_ESTIMATE | FFTW_IN_PLACE | FFTW_USE_WISDOM, qgrid, 1, 0, 0);
	//  if ( ! CkMyPe() ) iout << "   Done.\n" << endi;
//	printf("   Done.\n");  fflush(stdout);
	
	memset( (void*) qgrid, 0, qgrid_len * sizeof(double) );
	if ( myRecipPe >= numGridPes ) trans_count = numGridPes;
}

ComputePmeMgr::~ComputePmeMgr() {
	delete myKSpace;
	delete [] localInfo;
	delete [] recipPeMap;
	delete [] recipPeDest;
	delete [] qgrid;
	delete [] work;
}


void ComputePmeMgr::gridCalc1(void) {
	int i, j, k=0;
	for (i = 0;  i < myGrid.K1 * myGrid.K2;  i++) {
		for (j = 0;  j < myGrid.dim3;  j++) {
			qgrid[k++] = pmeCompute->q_arr[i][j];
		}
	}
	
	rfftwnd_real_to_complex(forward_plan_yz, localInfo[myRecipPe].nx,
		qgrid, 1, myGrid.dim2 * myGrid.dim3, 0, 0, 0);
	
	gridCalc2();
}


void ComputePmeMgr::gridCalc2(void) {
	int zdim = myGrid.dim3;
	int ny = localInfo[myRecipPe].ny_after_transpose;
	
	// finish forward FFT (x dimension)
	fftw(forward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
		ny * zdim / 2, 1, work, 1, 0);
	
	// reciprocal space portion of PME
	BigReal ewaldcof = pmeCompute->simParams->ewaldcof;
	recipEnergy = myKSpace->compute_energy(qgrid, lattice, ewaldcof, recip_vir);
	// CkPrintf("Ewald reciprocal energy = %f\n", recipEnergy2);
	
	// start backward FFT (x dimension)
	fftw(backward_plan_x, ny * zdim / 2, (fftw_complex *) qgrid,
		ny * zdim / 2, 1, work, 1, 0);
	
	gridCalc3();
}

void ComputePmeMgr::gridCalc3(void) {
	// finish backward FFT
	rfftwnd_complex_to_real(backward_plan_yz, localInfo[myRecipPe].nx,
		(fftw_complex *) qgrid, 1, myGrid.dim2 * myGrid.dim3 / 2, 0, 0, 0);
	
	
	int i, j, k=0;
	for (i = 0;  i < myGrid.K1 * myGrid.K2;  i++) {
		for (j = 0;  j < myGrid.dim3;  j++) {
			pmeCompute->q_arr[i][j] = qgrid[k++];
		}
	}
	pmeCompute->energy += recipEnergy;
	for (i = 0;  i < 6;  i++) {
		pmeCompute->virial[i] += recip_vir[i];
	}
	
	ungridCalc();
}

void ComputePmeMgr::ungridCalc(void) {
	pmeCompute->ungridForces();
	
	ungrid_count = numDestRecipPes;
}


static void scale_coordinates(PmeParticle p[], int N, Lattice lattice, PmeGrid grid) {
	Vector origin = lattice.o;
	Vector recip1 = lattice.b1;
	Vector recip2 = lattice.b2;
	Vector recip3 = lattice.b3;
	double ox = origin.x;
	double oy = origin.y;
	double oz = origin.z;
	double r1x = recip1.x;
	double r1y = recip1.y;
	double r1z = recip1.z;
	double r2x = recip2.x;
	double r2y = recip2.y;
	double r2z = recip2.z;
	double r3x = recip3.x;
	double r3y = recip3.y;
	double r3z = recip3.z;
	int K1 = grid.K1;
	int K2 = grid.K2;
	int K3 = grid.K3;
	
	for (int i=0; i<N; i++) {
		double px = p[i].x - ox;
		double py = p[i].y - oy;
		double pz = p[i].z - oz;
		double sx = px*r1x + py*r1y + pz*r1z;
		double sy = px*r2x + py*r2y + pz*r2z;
		double sz = px*r3x + py*r3y + pz*r3z;
		p[i].x = K1 * ( sx - floor(sx) );
		p[i].y = K2 * ( sy - floor(sy) );
		p[i].z = K3 * ( sz - floor(sz) );
		//  Check for rare rounding condition where K * ( 1 - epsilon ) == K
		//  which was observed with g++ on Intel x86 architecture.
		if ( p[i].x == K1 ) p[i].x = 0;
		if ( p[i].y == K2 ) p[i].y = 0;
		if ( p[i].z == K3 ) p[i].z = 0;
	}
}

static void scale_forces(Vector f[], int N, Lattice lattice) {
	Vector recip1 = lattice.b1;
	Vector recip2 = lattice.b2;
	Vector recip3 = lattice.b3;
	double r1x = recip1.x;
	double r1y = recip1.y;
	double r1z = recip1.z;
	double r2x = recip2.x;
	double r2y = recip2.y;
	double r2z = recip2.z;
	double r3x = recip3.x;
	double r3y = recip3.y;
	double r3z = recip3.z;
	
	for (int i=0; i<N; i++) {
		double f1 = f[i].x;
		double f2 = f[i].y;
		double f3 = f[i].z;
		f[i].x = f1*r1x + f2*r2x + f3*r3x;
		f[i].y = f1*r1y + f2*r2y + f3*r3y;
		f[i].z = f1*r1z + f2*r2z + f3*r3z;
	}
}


ComputePme::ComputePme(const PmetestParams &pmeParams) //ComputeID c) :
{
	computePmeMgr = new ComputePmeMgr;
	computePmeMgr->setCompute(this);
	
	useAvgPositions = 1;
	
	storePmeParams = pmeParams;
	simParams = &storePmeParams;
	myGrid.K1 = simParams->nxspacings;
	myGrid.K2 = simParams->nyspacings;
	myGrid.K3 = simParams->nzspacings;
	myGrid.order = simParams->interporder;
	myGrid.dim2 = myGrid.K2;
	myGrid.dim3 = 2 * (myGrid.K3/2 + 1);
	qsize = myGrid.K1 * myGrid.dim2 * myGrid.dim3;
	fsize = myGrid.K1 * myGrid.dim2;
	q_arr = new double*[fsize];
	memset( (void*) q_arr, 0, fsize * sizeof(double*) );
	// kludge so we won't segfault
	for (int i = 0;  i < fsize;  i++) {
		q_arr[i] = new double[myGrid.dim3];
		memset( (void*) q_arr[i], 0, myGrid.dim3 * sizeof(double) );
	}
	// end kludge
	f_arr = new char[fsize];
	fz_arr = new char[myGrid.K3];
	
	computePmeMgr->initialize();
}

ComputePme::~ComputePme()
{
	for (int i=0; i<fsize; ++i) {
		if ( q_arr[i] ) {
			delete [] q_arr[i];
		}
	}
	delete [] q_arr;
	delete [] f_arr;
	delete [] fz_arr;
	
	delete computePmeMgr;
}

void ComputePme::doWork(int natoms, const Vector *position,
                        const double *charge, Vector *resultForce,
                        double *resultEnergy, double *resultVirial)
{
	Lattice lattice;
	lattice_set(&lattice,
		simParams->cellvec1, simParams->cellvec2,
		simParams->cellvec3, simParams->center);
	numLocalAtoms = natoms;
	
	localData = new PmeParticle[numLocalAtoms];
	
	// get positions and charges
	// TAGIT
	// cout << "scaling= " << storePmeParams.scaling << endl;
	// cout << "dielectric_1= " << storePmeParams.dielectric_1 << endl;
	// cout << "ewaldcof= " << storePmeParams.ewaldcof << endl;
	
	PmeParticle * data_ptr = localData;
	
    int numAtoms = natoms, i;
	
    for(i=0; i<numAtoms; ++i)
    {
		data_ptr->x = position[i].x;
		data_ptr->y = position[i].y;
		data_ptr->z = position[i].z;
		data_ptr->cg = charge[i];
		++data_ptr;
    }
	
	for ( i=0; i<6; ++i ) { virial[i] = 0; }
	energy = 0;
	
	// calculate self energy
	BigReal ewaldcof = simParams->ewaldcof;
	BigReal selfEnergy = 0;
	data_ptr = localData;
	for(i=0; i<numLocalAtoms; ++i)
	{
		selfEnergy += data_ptr->cg * data_ptr->cg;
		++data_ptr;
	}
	selfEnergy *= -1. * ewaldcof / SQRT_PI;
	energy += selfEnergy;
	
	for (i=0; i<fsize; ++i) {
		if ( q_arr[i] ) {
			memset( (void*) (q_arr[i]), 0, myGrid.dim3 * sizeof(double) );
		}
	}
	memset( (void*) f_arr, 0, fsize * sizeof(char) );
	memset( (void*) fz_arr, 0, myGrid.K3 * sizeof(char) );
	myRealSpace = new PmeRealSpace(myGrid,numLocalAtoms);
	scale_coordinates(localData, numLocalAtoms, lattice, myGrid);
	myRealSpace->fill_charges(q_arr, f_arr, fz_arr, localData);
	
	pmeForce = resultForce;
	pmeVirial = resultVirial;
	computePmeMgr->gridCalc1();
	*resultEnergy = energy;
}

void ComputePme::ungridForces() {
	
    Vector *localResults = new Vector[numLocalAtoms];
    myRealSpace->compute_forces(q_arr, localData, localResults);
    delete [] localData;
    delete myRealSpace;
    Lattice lattice;
    lattice_set(&lattice,
		simParams->cellvec1, simParams->cellvec2,
		simParams->cellvec3, simParams->center);
    scale_forces(localResults, numLocalAtoms, lattice);
	
    Vector *results_ptr = localResults;
	
    // add in forces
	Vector *f = pmeForce; //r->f[Results::slow];
	int numAtoms = numLocalAtoms; //(*ap).p->getNumAtoms();
	
	for(int i=0; i<numAtoms; ++i)
	{
        f[i].x += results_ptr->x;
        f[i].y += results_ptr->y;
        f[i].z += results_ptr->z;
        ++results_ptr;
	}
	
    delete [] localResults;
	
	pmeVirial[0] += virial[0];
	pmeVirial[1] += virial[1];
	pmeVirial[2] += virial[2];
	pmeVirial[3] += virial[1];
	pmeVirial[4] += virial[3];
	pmeVirial[5] += virial[4];
	pmeVirial[6] += virial[2];
	pmeVirial[7] += virial[4];
	pmeVirial[8] += virial[5];
	
}
