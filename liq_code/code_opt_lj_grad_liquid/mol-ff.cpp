#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ff.h"

#ifndef _WIN32
#define MPI_ON 1
#endif

#define MAX_PROC	(72)

#ifdef MPI_ON
#include <mpi.h>
#include <unistd.h>	//for linux
#else
#include <direct.h>	//for win
#endif

#define MAX_N	(12288)
//#define MAX_N	(8192)
#define MAX_LJ_TYPE	(256)

enum {R_SELF, R_BOND, R_ANGLE, R_DIHEDRAL, R_CLOSE, R_SWITCH, R_FAR};

//start	DCD data
#define MOLFILE_SUCCESS           0   /**< succeeded in reading file      */
#define MOLFILE_EOF              -1   /**< end of file                    */
#define MOLFILE_ERROR            -1   /**< error reading/opening a file   */

#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2

#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */


typedef struct {
  FILE *fd;
  int natoms;
  int nsets;
  int setsread;
  int istart;
  int nsavc;
  double delta;
  int nfixed;
  float *x, *y, *z;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;  
  int first;
  int with_unitcell;
} dcdhandle;
//end	DCD data

// PBC box info
double L_Box=27.2, L_Half;
float A, B, C, alpha_loc, beta, my_gamma; 
dcdhandle *dcd;

void Quit_With_Error_Msg(char szMsg[]);
void Read_xpsf(char szName[]);
void Read_Dcd(char szName[]);
void Assign_LJ_Parameters(void);

static int read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART, 
                   int *NSAVC, double *DELTA, int *NAMNF, int **FREEINDEXES, float **fixedcoords, int *reverseEndian, int *charmm);
static int read_charmm_extrablock(FILE *fd, int charmm, int reverseEndian, float *unitcell);
static int read_next_timestep(void *v, int natoms);
void Free_Memory(void);
double Cal_E_VDW(double r_switch, double r_cut);
double Cal_dE_d_VDW_Param(double r_switch, double r_cut);
double Cal_E_VDW_LRC(void);
double Cal_dE_LRC_d_VDW_Param(void);
double Distance_1D_PBC(double x1, double x2);
void Accumulate_Gradients_VDW_Param(int iFrame);
void Accumulate_Gradients_V_VDW_Param(void);
void Output_Gradients_VDW_Param(int nFrames);
void Determine_Mol_Size(void);

FILE *fFile_Run_Log;	// will be shared by other source code

CMol Mol;
CForceField ForceField;

int nAtom, nTotal, nMol, MolSize;
char szName_FF[256], szName_Psf[256], szName_Dcd[256];

float x[MAX_N], y[MAX_N], z[MAX_N];
unsigned char Dist_Mat[MAX_N][MAX_N];
char szChemName[MAX_N][12], szAtomName[MAX_N][12], szMolName[MAX_N][12];
int  ResID[MAX_N];

//start	data related with VDW LRC
double r_Switch=10.0, r_Cutoff=12.0, E_LRC;
int LJtypecount, numAtomsByLJType[MAX_LJ_TYPE], List_LJ_Type[MAX_N];
double ATable[MAX_LJ_TYPE][MAX_LJ_TYPE], BTable[MAX_LJ_TYPE][MAX_LJ_TYPE];
double Sigma_12_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE], Sigma_6_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE], Sigma_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE];
double Sigma_11_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE], Sigma_5_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE];
double Epsilon_Sqrt[MAX_LJ_TYPE], Epsilon_14_Sqrt[MAX_LJ_TYPE], Epsilon_IJ_Mat[MAX_LJ_TYPE][MAX_LJ_TYPE];
double Epsilon_List[MAX_LJ_TYPE], Epsilon_14_List[MAX_LJ_TYPE];
//end	data related with VDW LRC

double dE_d_LJEmin[MAX_LJ_TYPE], dE_d_LJRmin[MAX_LJ_TYPE];
double dE_d_LJEmin_Acc[MAX_LJ_TYPE], dE_d_LJRmin_Acc[MAX_LJ_TYPE];
double E_dE_d_LJEmin_Acc[MAX_LJ_TYPE], E_dE_d_LJRmin_Acc[MAX_LJ_TYPE];

double V_dE_d_LJEmin_Acc[MAX_LJ_TYPE], V_dE_d_LJRmin_Acc[MAX_LJ_TYPE];
double T_Sim;	// simulation temperature
int Active_LJ[N_REC_MAX], nActiveLJ;
double Volume, Volume_Mol, Volume_Mol_Acc=0.0;
double *E_Liquid, E_Liquid_Avg;
int nFrame_Skip;
int printLRCFlag = 0; // chetan: print flag for vdw A/B coeffs and lrc energies to compare with NAMD

int nProc=1, ProgID=0;

int main(int argc, char *argv[])
{
#ifdef MPI_ON
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
#endif

	if(argc != 9)	{
		printf("Usage: gas_sim ff.str mol.xpsf traj.dcd r_Switch r_Cutoff nFrame_Skip T_Sim L_Box\n");
#ifdef MPI_ON
		MPI_Finalize();
#endif
		exit(1);
	}


	strcpy(szName_FF, argv[1]);
	strcpy(szName_Psf, argv[2]);
	strcpy(szName_Dcd, argv[3]);
	r_Switch = atof(argv[4]);
	r_Cutoff = atof(argv[5]);
	nFrame_Skip = atoi(argv[6]);
	T_Sim = atoi(argv[7]);
	L_Box = atof(argv[8]); // the length of box changes in NPT simulation and is based on DCD

	fFile_Run_Log = stdout;

	ForceField.ReadForceField(szName_FF);
	Read_xpsf(szName_Psf);

	Assign_LJ_Parameters();

	memset(dE_d_LJEmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);
	memset(dE_d_LJRmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);
	memset(E_dE_d_LJEmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);
	memset(E_dE_d_LJRmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);
	memset(V_dE_d_LJEmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);
	memset(V_dE_d_LJRmin_Acc, 0, sizeof(double)*MAX_LJ_TYPE);

	Read_Dcd(szName_Dcd);

	Free_Memory();

#ifdef MPI_ON
	MPI_Finalize();
#endif

	return 0;
}

void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);
	exit(1);
}


void Read_xpsf(char szName[])
{
	FILE *fIn;
	char szErrorMsg[256], szLine[256], *ReadLine, szBuff[256], szMol[256], szBuff2[256];
	int i, j, k, l, ReadItem, iTmp, nBond, nBondRead, IdxList[16];
	int AtomBond[MAX_N], Bond_Array[MAX_N][8], Atom_i, Atom_j, Atom_k, Atom_l, Atom_Host;
	double fTmp=0.0, mass[MAX_N];

	memset(AtomBond, 0, sizeof(int)*MAX_N);
	memset(Bond_Array, 0, sizeof(int)*MAX_N*8);

	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		sprintf(szErrorMsg, "Fail to open file: %s\nQuit\n", szName);
		Quit_With_Error_Msg(szErrorMsg);
	}

	nAtom = -1;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine==NULL)	{
			break;
		}
		ReadItem = sscanf(szLine, "%d%s", &nAtom, szBuff);
		if( (ReadItem==2) && (strcmp(szBuff, "!NATOM")==0) )	{
			break;
		}
	}

	if(nAtom < 0)	{
		fclose(fIn);
		Quit_With_Error_Msg("nAtom < 0!\nQuit\n");
	}

	if(nAtom > MAX_N)	{
		Quit_With_Error_Msg("Fatal error!\nnAtom > MAX_N\n");
	}

	for(i=0; i<nAtom; i++)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine==NULL)	{
			fclose(fIn);
			Quit_With_Error_Msg("Error in reading xpsf file.\n");
		}
		ReadItem = sscanf(szLine, "%d%s%d%s%s%s%lf%lf%d%lf%lf", &iTmp, szMolName[i], &(ResID[i]), szMol, szAtomName[i], szChemName[i], &fTmp, &(mass[i]), &iTmp, &fTmp, &fTmp);
		if(ReadItem != 11)	{
			Quit_With_Error_Msg("Error in reading xpsf file.\n");
		}
	}

	for(i=0; i<nAtom; i++)	{
		for(j=i; j<nAtom; j++)	{
			Dist_Mat[i][j] = Dist_Mat[j][i] = 100;	// long distance
		}
		Dist_Mat[i][i] = 0;
	}

	//start	reading entries of bond
	while(1)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine==NULL)	{
			fclose(fIn);
			Quit_With_Error_Msg("Fail to find the entry for bond info in xpsf file.\n");
		}
		ReadItem = sscanf(szLine, "%d%s%s", &nBond, szBuff, szBuff2);
		if(strncmp(szBuff, "!NBOND", 5)==0)	{
			break;
		}
	}

	nBondRead = 0;
	while(nBondRead < nBond)	{
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine==NULL)	{
			fclose(fIn);
			Quit_With_Error_Msg("Error in reading the entry for bond info in xpsf file.\n");
		}
		ReadItem = sscanf(szLine, "%d%d%d%d%d%d%d%d", &(IdxList[0]), &(IdxList[1]), &(IdxList[2]), &(IdxList[3]), &(IdxList[4]), &(IdxList[5]), &(IdxList[6]), &(IdxList[7]));

		if( (ReadItem > 0) && (ReadItem%2==0) )	{
			for(i=0; i<ReadItem; i++)	{
				IdxList[i] -= 1;
			}
			for(i=0; i<ReadItem; i+=2)	{
				Dist_Mat[IdxList[i]][IdxList[i+1]] = 1;
				Dist_Mat[IdxList[i+1]][IdxList[i]] = 1;

				Bond_Array[IdxList[i]][AtomBond[IdxList[i]]] = IdxList[i+1];
				Bond_Array[IdxList[i+1]][AtomBond[IdxList[i+1]]] = IdxList[i];
				AtomBond[IdxList[i]]++;
				AtomBond[IdxList[i+1]]++;
			}
			nBondRead += (ReadItem/2);
		}
		else	{
			fclose(fIn);
			Quit_With_Error_Msg("Error in reading the entry for bond info in xpsf file.\n");
		}
	}
	//end	reading entries of bond

	fclose(fIn);

	//start	to constrcut distance matrix
	for(i=0; i<nAtom; i++)	{	//start	enumeration
		Atom_i = i;

		for(j=0; j<AtomBond[Atom_i]; j++)	{	//i->j
			Atom_j = Bond_Array[Atom_i][j];

			if(Dist_Mat[Atom_i][Atom_j] > 1)	{
				Dist_Mat[Atom_i][Atom_j] = 1;
				Dist_Mat[Atom_j][Atom_i] = 1;
			}

			for(k=0; k<AtomBond[Atom_j]; k++)	{	//i->j->k
				Atom_k = Bond_Array[Atom_j][k];

				if(Dist_Mat[Atom_i][Atom_k] > 2)	{
					Dist_Mat[Atom_i][Atom_k] = 2;
					Dist_Mat[Atom_k][Atom_i] = 2;
				}
				
				for(l=0; l<AtomBond[Atom_k]; l++)	{	//i->j->k->l
					Atom_l = Bond_Array[Atom_k][l];
					
					if(Dist_Mat[Atom_i][Atom_l] > 3)	{
						Dist_Mat[Atom_i][Atom_l] = 3;
						Dist_Mat[Atom_l][Atom_i] = 3;
					}
				}
			}
		}
	}
	//end	to constrcut distance matrix

	//start	to substitute the exclusion list of host atom to drudes' list
	for(i=0; i<nAtom; i++)	{
		if( (mass[i] > 1.0E-10) && (mass[i] < 0.8) )	{	// Drude particle
			Atom_i = i;
			Atom_Host = Bond_Array[Atom_i][0];

			for(j=0; j<nAtom; j++)	{
				if(Dist_Mat[Atom_Host][j] < Dist_Mat[Atom_i][j])	{
					Dist_Mat[Atom_i][j] = Dist_Mat[Atom_Host][j];
					Dist_Mat[j][Atom_i] = Dist_Mat[Atom_i][j];
				}
			}
		}
	}
	//end	to substitute the exclusion list of host atom to drudes' list

	//start	to substitute the exclusion list of host atom to LPs' list
	for(i=0; i<nAtom; i++)	{
		if(mass[i] < 1.0E-10)	{	// zero mass, it should be a lone pair.
			Atom_i = i;
			Atom_Host = Bond_Array[Atom_i][0];

			for(j=0; j<nAtom; j++)	{
				if(Dist_Mat[Atom_Host][j] < Dist_Mat[Atom_i][j])	{
					Dist_Mat[Atom_i][j] = Dist_Mat[Atom_Host][j];
					Dist_Mat[j][Atom_i] = Dist_Mat[Atom_Host][j];
				}
			}
		}
	}
	//end	to substitute the exclusion list of host atom to LPs' list

	Determine_Mol_Size();
}


void Read_Dcd(char szName[])
{
	int rc, iFrame, ReadItem;
	FILE *fd, *fIn;
	char szName_E_NAMD[256]="E_namd_dcd.txt"; //Make sure output steps of energy print and dcd matches!!
	double E_VDW, E_LRC;
	
	fd = fopen(szName, "rb");
	
	if(fd == NULL)	{
		Quit_With_Error_Msg("Fail to open dcd file.\nQuit.\n");
	}
	
	dcd = new dcdhandle;
	memset(dcd, 0, sizeof(dcdhandle));
	dcd->fd = fd;
	
	if ((rc = read_dcdheader(dcd->fd, &dcd->natoms, &dcd->nsets, &dcd->istart, 
		&dcd->nsavc, &dcd->delta, &dcd->nfixed, &dcd->freeind, 
		&dcd->fixedcoords, &dcd->reverse, &dcd->charmm))) {
		fclose(dcd->fd);
		free(dcd);

		Quit_With_Error_Msg("Error in reading dcd head file.\n");
	}

	dcd->first = 1;
	nTotal = dcd->natoms;

	if( nTotal != nAtom )	{
		Quit_With_Error_Msg("nTotal != nAtom\nQuit!\n");
	}

	E_Liquid = new double[dcd->nsets];	// Please make sure NAMD energy output and DCD are synchronized.
	E_Liquid_Avg = 0.0;

	fIn = fopen(szName_E_NAMD, "r");
	if(fIn == NULL)	{
		printf("Fail to open file %s\nQuit.\n", szName_E_NAMD);
		exit(1);
	}
	for(iFrame=0; iFrame<dcd->nsets; iFrame++)	{
		ReadItem = fscanf(fIn, "%lf", &(E_Liquid[iFrame]));
		if(ReadItem != 1)	{
			printf("Error in reading file: %s at %d step.\nQuit.\n", szName_E_NAMD, iFrame+1);
			exit(1);
		}
		if(iFrame >= nFrame_Skip)	{
			E_Liquid_Avg += (E_Liquid[iFrame]);
		}
	}
	fclose(fIn);
	E_Liquid_Avg /= (dcd->nsets - nFrame_Skip);

	
	dcd->x = new float[dcd->natoms];
	dcd->y = new float[dcd->natoms];
	dcd->z = new float[dcd->natoms];

	if (!dcd->x || !dcd->y || !dcd->z) {
		fclose(dcd->fd);
		Free_Memory();
		Quit_With_Error_Msg("Unable to allocate space for reading dcd file.\n");
	}

	for(iFrame=0; iFrame<nFrame_Skip; iFrame++)	{
		read_next_timestep(dcd, dcd->natoms);
	}

	for(; iFrame<dcd->nsets; iFrame++)	{
		read_next_timestep(dcd, dcd->natoms);

		memcpy(x, dcd->x, sizeof(float)*nAtom);
		memcpy(y, dcd->y, sizeof(float)*nAtom);
		memcpy(z, dcd->z, sizeof(float)*nAtom);

		Volume = 1.0 * A*B*C * 1.0;	// make sure double precision is used !!
		Volume_Mol = Volume / nMol;
		L_Box = A;	// update the box size info

		E_VDW = Cal_dE_d_VDW_Param(r_Switch, r_Cutoff);
		E_LRC = Cal_dE_LRC_d_VDW_Param();
		printf("E_VDW = %f E_LRC = %f Total_VDW = %f\n",E_VDW,E_LRC,E_VDW+E_LRC);
		Accumulate_Gradients_VDW_Param(iFrame);
		Accumulate_Gradients_V_VDW_Param();
	}	
	fclose(fd);

	delete []E_Liquid;

	if(ProgID == 0)	Output_Gradients_VDW_Param(dcd->nsets - nFrame_Skip);
}


static int read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART, 
                   int *NSAVC, double *DELTA, int *NAMNF, 
                   int **FREEINDEXES, float **fixedcoords, int *reverseEndian, 
                   int *charmm)
{
  int input_integer[2];  /* buffer space */
  int i, ret_val, rec_scale;
  char hdrbuf[84];    /* char buffer used to store header */
  int NTITLE;
  int dcdcordmagic;
  char *corp = (char *) &dcdcordmagic;

  /* coordinate dcd file magic string 'CORD' */
  corp[0] = 'C';
  corp[1] = 'O';
  corp[2] = 'R';
  corp[3] = 'D';

  /* First thing in the file should be an 84.
   * some 64-bit compiles have a 64-bit record length indicator,
   * so we have to read two ints and check in a more complicated 
   * way. :-( */

  ret_val = fread(input_integer, sizeof(int), 2, fd);

  /* Check magic number in file header and determine byte order*/
  if ((input_integer[0]+input_integer[1]) == 84) {
    *reverseEndian=0;
    rec_scale=RECSCALE64BIT;
    printf("dcdplugin) detected CHARMM -i8 64-bit DCD file of native endianness\n");
  } else if (input_integer[0] == 84 && input_integer[1] == dcdcordmagic) {
    *reverseEndian=0;
    rec_scale=RECSCALE32BIT;
//    printf("dcdplugin) detected standard 32-bit DCD file of native endianness\n");
  } else {
    /* now try reverse endian */
//    swap4_aligned(input_integer, 2); /* will have to unswap magic if 32-bit */
    if ((input_integer[0]+input_integer[1]) == 84) {
      *reverseEndian=1;
      rec_scale=RECSCALE64BIT;
      printf("dcdplugin) detected CHARMM -i8 64-bit DCD file of opposite endianness\n");
    } else {
 //     swap4_aligned(&input_integer[1], 1); /* unswap magic (see above) */
      if (input_integer[0] == 84 && input_integer[1] == dcdcordmagic) {
        *reverseEndian=1;
        rec_scale=RECSCALE32BIT;
        printf("dcdplugin) detected standard 32-bit DCD file of opposite endianness\n");
      } else {
        /* not simply reversed endianism or -i8, something rather more evil */
        printf("dcdplugin) unrecognized DCD header:\n");
        printf("dcdplugin)   [0]: %10d  [1]: %10d\n", input_integer[0], input_integer[1]);
        printf("dcdplugin)   [0]: 0x%08x  [1]: 0x%08x\n", input_integer[0], input_integer[1]);
        return DCD_BADFORMAT;

      }
    }
  }

  /* check for magic string, in case of long record markers */
  if (rec_scale == RECSCALE64BIT) { 
    ret_val = fread(input_integer, sizeof(int), 1, fd);
    if (input_integer[0] != dcdcordmagic) {
      printf("dcdplugin) failed to find CORD magic in CHARMM -i8 64-bit DCD file\n");
      return DCD_BADFORMAT;
    }
  }

  /* Buffer the entire header for random access */
  ret_val = fread(hdrbuf, sizeof(char), 80, fd);

  /* CHARMm-genereate DCD files set the last integer in the     */
  /* header, which is unused by X-PLOR, to its version number.  */
  /* Checking if this is nonzero tells us this is a CHARMm file */
  /* and to look for other CHARMm flags.                        */
  if (*((int *) (hdrbuf + 76)) != 0) {
    (*charmm) = DCD_IS_CHARMM;
    if (*((int *) (hdrbuf + 40)) != 0)
      (*charmm) |= DCD_HAS_EXTRA_BLOCK;

    if (*((int *) (hdrbuf + 44)) == 1)
      (*charmm) |= DCD_HAS_4DIMS;

    if (rec_scale == RECSCALE64BIT)
      (*charmm) |= DCD_HAS_64BIT_REC;
  
  } else {
    (*charmm) = DCD_IS_XPLOR; /* must be an X-PLOR format DCD file */
  }

  if (*charmm & DCD_IS_CHARMM) {
    /* CHARMM and NAMD versions 2.1b1 and later */
//    printf("dcdplugin) CHARMM format DCD file (also NAMD 2.1 and later)\n");
  } else {
    /* CHARMM and NAMD versions prior to 2.1b1  */
    printf("dcdplugin) X-PLOR format DCD file (also NAMD 2.0 and earlier)\n");
  }


  /* Store the number of sets of coordinates (NSET) */
  (*NSET) = *((int *) (hdrbuf));

  /* Store ISTART, the starting timestep */
  (*ISTART) = *((int *) (hdrbuf + 4));

  /* Store NSAVC, the number of timesteps between dcd saves */
  (*NSAVC) = *((int *) (hdrbuf + 8));

  /* Store NAMNF, the number of fixed atoms */
  (*NAMNF) = *((int *) (hdrbuf + 32));



  /* Read in the timestep, DELTA */
  /* Note: DELTA is stored as a double with X-PLOR but as a float with CHARMm */
  if ((*charmm) & DCD_IS_CHARMM) {
    float ftmp;
    ftmp = *((float *)(hdrbuf+36)); /* is this safe on Alpha? */

    *DELTA = (double)ftmp;
  } else {
    (*DELTA) = *((double *)(hdrbuf + 36));
  }

  /* Get the end size of the first block */
  ret_val = fread(input_integer, sizeof(int), rec_scale, fd);

  if (input_integer[0] != 84) {
      return DCD_BADFORMAT;
  }
  
  /* Read in the size of the next block */
  input_integer[1] = 0;
  ret_val = fread(input_integer, sizeof(int), rec_scale, fd);

  if ((((input_integer[0]+input_integer[1])-4) % 80) == 0) {
    /* Read NTITLE, the number of 80 character title strings there are */
    ret_val = fread(&NTITLE, sizeof(int), 1, fd);

    if (NTITLE < 0) {
      printf("dcdplugin) WARNING: Bogus NTITLE value: %d (hex: %08x)\n", 
             NTITLE, NTITLE);
      return DCD_BADFORMAT;
    }

    if (NTITLE > 1000) {
      printf("dcdplugin) WARNING: Bogus NTITLE value: %d (hex: %08x)\n", 
             NTITLE, NTITLE);
      if (NTITLE == 1095062083) {
        printf("dcdplugin) WARNING: Broken Vega ZZ 2.4.0 DCD file detected\n");
        printf("dcdplugin) Assuming 2 title lines, good luck...\n");
        NTITLE = 2;
      } else {
        printf("dcdplugin) Assuming zero title lines, good luck...\n");
        NTITLE = 0;
      }
    }

    for (i=0; i<NTITLE; i++) {
      fseek(fd, 80, SEEK_CUR);
    }

    /* Get the ending size for this block */
    ret_val = fread(input_integer, sizeof(int), rec_scale, fd);
  } else {
    return DCD_BADFORMAT;
  }

  /* Read in an integer '4' */
  input_integer[1] = 0;
  ret_val = fread(input_integer, sizeof(int), rec_scale, fd);
  
  if ((input_integer[0]+input_integer[1]) != 4) {
    return DCD_BADFORMAT;
  }

  /* Read in the number of atoms */
  ret_val = fread(N, sizeof(int), 1, fd);

  /* Read in an integer '4' */
  input_integer[1] = 0;
  ret_val = fread(input_integer, sizeof(int), rec_scale, fd);

  if ((input_integer[0]+input_integer[1]) != 4) {
    return DCD_BADFORMAT;
  }

  *FREEINDEXES = NULL;
  *fixedcoords = NULL;
  if (*NAMNF != 0) {
    (*FREEINDEXES) = (int *) calloc(((*N)-(*NAMNF)), sizeof(int));
    if (*FREEINDEXES == NULL)
      return DCD_BADMALLOC;

    *fixedcoords = (float *) calloc((*N)*4 - (*NAMNF), sizeof(float));
    if (*fixedcoords == NULL)
      return DCD_BADMALLOC;

    /* Read in index array size */
    input_integer[1]=0;
    ret_val = fread(input_integer, sizeof(int), rec_scale, fd);

    if ((input_integer[0]+input_integer[1]) != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }

    ret_val = fread((*FREEINDEXES), sizeof(int), ((*N)-(*NAMNF)), fd);

    input_integer[1]=0;
    ret_val = fread(input_integer, sizeof(int), rec_scale, fd);

    if ((input_integer[0]+input_integer[1]) != ((*N)-(*NAMNF))*4) {
      return DCD_BADFORMAT;
    }
  }


  return DCD_SUCCESS;
}

static int read_dcdstep(FILE *fd, int N, float *X, float *Y, float *Z, 
                        float *unitcell, int num_fixed,
                        int first, int *indexes, float *fixedcoords, 
                        int reverseEndian, int charmm) {
	int ret_val, rec_scale;   /* Return value from read */
	int ReadItem;
	
	if (charmm & DCD_HAS_64BIT_REC) {
		rec_scale=RECSCALE64BIT;
	} else {
		rec_scale=RECSCALE32BIT;
	}
	
	if ((num_fixed==0) || first) {
		/* temp storage for reading formatting info */
		/* note: has to be max size we'll ever use  */
		int tmpbuf[6*RECSCALEMAX]; 
		int i;
		
		/* if there are no fixed atoms or this is the first timestep read */
		/* then we read all coordinates normally.                         */
		
		/* read the charmm periodic cell information */
		/* XXX this too should be read together with the other items in a */
		/*     single fio_readv() call in order to prevent lots of extra  */
		/*     kernel/user context switches.                              */
		ret_val = read_charmm_extrablock(fd, charmm, reverseEndian, unitcell);
		if (ret_val) return ret_val;
		
		ReadItem = fread(&tmpbuf[0], sizeof(int), rec_scale, fd);
		if(ReadItem != rec_scale)	return DCD_BADREAD;
		ReadItem = fread(X, sizeof(float), N, fd);
		if(ReadItem != N)	return DCD_BADREAD;
		ReadItem = fread(&tmpbuf[1*rec_scale], sizeof(int), rec_scale*2, fd);
		if(ReadItem != rec_scale*2)	return DCD_BADREAD;
		ReadItem = fread(Y, sizeof(float), N, fd);
		if(ReadItem != N)	return DCD_BADREAD;
		ReadItem = fread(&tmpbuf[3*rec_scale], sizeof(int), rec_scale*2, fd);
		if(ReadItem != rec_scale*2)	return DCD_BADREAD;
		ReadItem = fread(Z, sizeof(float), N, fd);
		if(ReadItem != N)	return DCD_BADREAD;
		ReadItem = fread(&tmpbuf[5*rec_scale], sizeof(int), rec_scale, fd);
		if(ReadItem != rec_scale)	return DCD_BADREAD;
		
		/* convert endianism if necessary */
		if (reverseEndian) {
			//     swap4_aligned(&tmpbuf[0], rec_scale*6);
			//     swap4_aligned(X, N);
			//     swap4_aligned(Y, N);
			//     swap4_aligned(Z, N);
		}
		
		/* double-check the fortran format size values for safety */
		if(rec_scale == 1) {
			for (i=0; i<6; i++) {
				if (tmpbuf[i] != (signed int)sizeof(float)*N) return DCD_BADFORMAT;
			}
		} else {
			for (i=0; i<6; i++) {
				if ((tmpbuf[2*i]+tmpbuf[2*i+1]) != (signed int)sizeof(float)*N) return DCD_BADFORMAT;
			}
		}
		
		/* copy fixed atom coordinates into fixedcoords array if this was the */
		/* first timestep, to be used from now on.  We just copy all atoms.   */
		if (num_fixed && first) {
			memcpy(fixedcoords, X, N*sizeof(float));
			memcpy(fixedcoords+N, Y, N*sizeof(float));
			memcpy(fixedcoords+2*N, Z, N*sizeof(float));
		}
		
		/* read in the optional charmm 4th array */
		/* XXX this too should be read together with the other items in a */
		/*     single fio_readv() call in order to prevent lots of extra  */
		/*     kernel/user context switches.                              */
//		ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
//		if (ret_val) return ret_val;
	} else {
		/* if there are fixed atoms, and this isn't the first frame, then we */
		/* only read in the non-fixed atoms for all subsequent timesteps.    */
//		ret_val = read_charmm_extrablock(fd, charmm, reverseEndian, unitcell);
//		if (ret_val) return ret_val;
//		ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//			fixedcoords, fixedcoords+3*N, X, charmm);
//		if (ret_val) return ret_val;
//		ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//			fixedcoords+N, fixedcoords+3*N, Y, charmm);
//		if (ret_val) return ret_val;
//		ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//			fixedcoords+2*N, fixedcoords+3*N, Z, charmm);
//		if (ret_val) return ret_val;
//		ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
//		if (ret_val) return ret_val;
	}
	
	return DCD_SUCCESS;
}


static int read_charmm_extrablock(FILE *fd, int charmm, int reverseEndian, float *unitcell)
{
  int i, input_integer[2], rec_scale;

  if (charmm & DCD_HAS_64BIT_REC) {
    rec_scale = RECSCALE64BIT;
  } else {
    rec_scale = RECSCALE32BIT;
  }

  if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
    /* Leading integer must be 48 */
    input_integer[1] = 0;
    if (fread(input_integer, sizeof(int), rec_scale, fd) != (unsigned int)rec_scale)
      return DCD_BADREAD; 
//    if (reverseEndian) swap4_aligned(input_integer, rec_scale);
    if ((input_integer[0]+input_integer[1]) == 48) {
      double tmp[6];
      if (fread(tmp, 48, 1, fd) != 1) return DCD_BADREAD;
//      if (reverseEndian) 
//        swap8_aligned(tmp, 6);
      for (i=0; i<6; i++) unitcell[i] = (float)tmp[i];
    } else {
      /* unrecognized block, just skip it */
      if (fseek(fd, (input_integer[0]+input_integer[1]), SEEK_CUR)) return DCD_BADREAD;
    }
    if (fread(input_integer, sizeof(int), rec_scale, fd) != (unsigned int)rec_scale) return DCD_BADREAD; 
  } 

  return DCD_SUCCESS;
}

#define PI_HALF (1.57079632679489661922)
static int read_next_timestep(void *v, int natoms) {
  dcdhandle *dcd;
  int rc;
  float unitcell[6];
  unitcell[0] = unitcell[2] = unitcell[5] = 1.0f;
  unitcell[1] = unitcell[3] = unitcell[4] = 90.0f;
  dcd = (dcdhandle *)v;

  /* Check for EOF here; that way all EOF's encountered later must be errors */
  if (dcd->setsread == dcd->nsets) return MOLFILE_EOF;
  dcd->setsread++;
  rc = read_dcdstep(dcd->fd, dcd->natoms, dcd->x, dcd->y, dcd->z, unitcell,
             dcd->nfixed, dcd->first, dcd->freeind, dcd->fixedcoords, 
             dcd->reverse, dcd->charmm);
  dcd->first = 0;
  if (rc < 0) {  
    Quit_With_Error_Msg("read_dcdstep");
  }


  A = unitcell[0];
  B = unitcell[2];
  C = unitcell[5];

  if (unitcell[1] >= -1.0 && unitcell[1] <= 1.0 &&
      unitcell[3] >= -1.0 && unitcell[3] <= 1.0 &&
      unitcell[4] >= -1.0 && unitcell[4] <= 1.0) {
    /* This file was generated by CHARMM, or by NAMD > 2.5, with the angle */
    /* cosines of the periodic cell angles written to the DCD file.        */ 
    /* This formulation improves rounding behavior for orthogonal cells    */
    /* so that the angles end up at precisely 90 degrees, unlike acos().   */
    alpha_loc = (float)(90.0 - asin(unitcell[4]) * 90.0 / PI_HALF); /* cosBC */
    beta  = (float)(90.0 - asin(unitcell[3]) * 90.0 / PI_HALF); /* cosAC */
    my_gamma = (float)(90.0 - asin(unitcell[1]) * 90.0 / PI_HALF); /* cosAB */
  } else {
    /* This file was likely generated by NAMD 2.5 and the periodic cell    */
    /* angles are specified in degrees rather than angle cosines.          */
    alpha_loc = unitcell[4]; /* angle between B and C */
    beta  = unitcell[3]; /* angle between A and C */
    my_gamma = unitcell[1]; /* angle between A and B */
  }
 
  return 1;
}
#undef PI_HALF


void Free_Memory(void)
{
	if(dcd)	{
		if (dcd->x)	delete [](dcd->x);
		if (dcd->y)	delete [](dcd->y);
		if (dcd->z)	delete [](dcd->z);
		delete dcd;
	}
}


double Cal_E_VDW(double r_switch, double r_cut)
{
	int i, j, LJ_ID_I, LJ_ID_J;
	double E_VDW=0.0, r_switch_SQ, r_cut_SQ, dx, dy, dz, r_SQ;
	double A, B, Inv_R2, Inv_R_6, Inv_R_12, dE, switch_coef, switch_coef0, dr_switch_cut, switch_coef_2;

	L_Half = 0.5 * L_Box;

	r_switch_SQ = r_switch*r_switch;
	r_cut_SQ = r_cut*r_cut;
	dr_switch_cut = r_cut_SQ - r_switch_SQ;
	switch_coef0 = 1.0/( dr_switch_cut * dr_switch_cut * dr_switch_cut );

	for(i=0; i<nAtom; i++)	{
		LJ_ID_I = List_LJ_Type[i];
		if( (LJ_ID_I < 0) || (Epsilon_Sqrt[LJ_ID_I] < 1.0E-50) )	{
			continue;
		}
		for(j=i+1; j<nAtom; j++)	{
			LJ_ID_J = List_LJ_Type[j];
			if( (LJ_ID_J < 0) || (Epsilon_Sqrt[LJ_ID_J] < 1.0E-50) )	{
				continue;
			}
			if(Dist_Mat[i][j] <= R_ANGLE)	{
				continue;
			}

			dx = Distance_1D_PBC(x[i], x[j]);
			dy = Distance_1D_PBC(y[i], y[j]);
			dz = Distance_1D_PBC(z[i], z[j]);
			r_SQ = dx*dx + dy*dy + dz*dz;

			if(r_SQ <= r_cut_SQ)	{
				A = ATable[LJ_ID_I][LJ_ID_J];
				B = BTable[LJ_ID_I][LJ_ID_J];
				
				Inv_R2 = 1.0 / r_SQ;
				Inv_R_6 = Inv_R2 * Inv_R2 * Inv_R2;
				Inv_R_12 = Inv_R_6*Inv_R_6;
				dE = A * Inv_R_12 - B * Inv_R_6;
				
				if(r_SQ > r_switch_SQ)	{
					switch_coef = switch_coef0 * (r_cut_SQ-r_SQ) * (r_cut_SQ-r_SQ)
						* (r_cut_SQ + 2.0*r_SQ - 3.0*r_switch_SQ);
					switch_coef_2 = 12.0 * (r_cut_SQ-r_SQ) * (r_SQ - r_switch_SQ) * switch_coef0;
					dE *= switch_coef;
				}
				else	{
				}

				if(Dist_Mat[i][j] == R_DIHEDRAL)	{	// 1-4 pair. Assuming Emin_14 = 0.5*Emin or Emin_14 = Emin
					dE *= ( (Epsilon_14_Sqrt[LJ_ID_I]*Epsilon_14_Sqrt[LJ_ID_J])/(Epsilon_Sqrt[LJ_ID_I]*Epsilon_Sqrt[LJ_ID_J]) );
				}

				E_VDW += dE;
			}
			
		}
	}

	return E_VDW;
}

double Cal_dE_d_VDW_Param(double r_switch, double r_cut)
{
	int i, j, LJ_ID_I, LJ_ID_J;
	double E_VDW=0.0, r_switch_SQ, r_cut_SQ, dx, dy, dz, r_SQ;
	double A, B, Inv_R2, Inv_R_6, Inv_R_12, dE, switch_coef, switch_coef0, dr_switch_cut;
	double dE_d_Emin, dE_d_Rmin, Scale_Emin_14;
	double dE_d_LJEmin_Sum[MAX_LJ_TYPE], dE_d_LJRmin_Sum[MAX_LJ_TYPE];
	double E_dE_d_LJEmin_Sum[MAX_LJ_TYPE], E_dE_d_LJRmin_Sum[MAX_LJ_TYPE];

	memset(dE_d_LJEmin, 0, sizeof(double)*ForceField.n_Rec_LJ);
	memset(dE_d_LJRmin, 0, sizeof(double)*ForceField.n_Rec_LJ);

	printf("%s %f\n","L_BOX_vdw",L_Box);

	L_Half = 0.5 * L_Box;

	r_switch_SQ = r_switch*r_switch;
	r_cut_SQ = r_cut*r_cut;
	dr_switch_cut = r_cut_SQ - r_switch_SQ;
	switch_coef0 = 1.0/( dr_switch_cut * dr_switch_cut * dr_switch_cut );

	for(i=0; i<nAtom; i++)	{
		if(i%nProc != ProgID)	{	// parallel
			continue;
		}

		LJ_ID_I = List_LJ_Type[i];
		if( (LJ_ID_I < 0) || (Epsilon_Sqrt[LJ_ID_I] < 1.0E-50) )	{
			continue;
		}
		for(j=i+1; j<nAtom; j++)	{
			LJ_ID_J = List_LJ_Type[j];
			if( (LJ_ID_J < 0) || (Epsilon_Sqrt[LJ_ID_J] < 1.0E-50) )	{
				continue;
			}
			if(Dist_Mat[i][j] <= R_ANGLE)	{
				continue;
			}

			dx = Distance_1D_PBC(x[i], x[j]);
			dy = Distance_1D_PBC(y[i], y[j]);
			dz = Distance_1D_PBC(z[i], z[j]);
			r_SQ = dx*dx + dy*dy + dz*dz;

			if(r_SQ <= r_cut_SQ)	{
				A = ATable[LJ_ID_I][LJ_ID_J];
				B = BTable[LJ_ID_I][LJ_ID_J];
				
				Inv_R2 = 1.0 / r_SQ;
				Inv_R_6 = Inv_R2 * Inv_R2 * Inv_R2;
				Inv_R_12 = Inv_R_6*Inv_R_6;
				dE = A * Inv_R_12 - B * Inv_R_6;
				
				dE_d_Emin = ( ( Sigma_12_Mat[LJ_ID_I][LJ_ID_J]*Inv_R_12 - 2.0*(Sigma_6_Mat[LJ_ID_I][LJ_ID_J]*Inv_R_6) ) * 0.5);
				dE_d_Rmin = ( Epsilon_IJ_Mat[LJ_ID_I][LJ_ID_J] * ( Sigma_12_Mat[LJ_ID_I][LJ_ID_J]*Inv_R_12 - (Sigma_6_Mat[LJ_ID_I][LJ_ID_J]*Inv_R_6) ) * 12.0 / (Sigma_Mat[LJ_ID_I][LJ_ID_J]));

				if(r_SQ > r_switch_SQ)	{
					switch_coef = switch_coef0 * (r_cut_SQ-r_SQ) * (r_cut_SQ-r_SQ)
						* (r_cut_SQ + 2.0*r_SQ - 3.0*r_switch_SQ);
//					switch_coef_2 = 12.0 * (r_cut_SQ-r_SQ) * (r_SQ - r_switch_SQ) * switch_coef0;	// used for force only
					dE *= switch_coef;
					dE_d_Emin *= switch_coef;
					dE_d_Rmin *= switch_coef;
					//printf("dE running: %f\n",dE);
				}
				else	{
				}

				if(Dist_Mat[i][j] == R_DIHEDRAL)	{	// 1-4 pair. Assuming Emin = 0.5*Emin
					Scale_Emin_14 = ( (Epsilon_14_Sqrt[LJ_ID_I]*Epsilon_14_Sqrt[LJ_ID_J])/(Epsilon_Sqrt[LJ_ID_I]*Epsilon_Sqrt[LJ_ID_J]) );
					dE *= Scale_Emin_14;
//					dE_d_Emin *= Scale_Emin_14;	// ??
					dE_d_Rmin *= Scale_Emin_14;

					dE_d_LJEmin[LJ_ID_I] += ( dE_d_Emin * Epsilon_14_Sqrt[LJ_ID_J] * Epsilon_14_List[LJ_ID_I] / (Epsilon_14_Sqrt[LJ_ID_I] * Epsilon_List[LJ_ID_I]) );
					dE_d_LJEmin[LJ_ID_J] += ( dE_d_Emin * Epsilon_14_Sqrt[LJ_ID_I] * Epsilon_14_List[LJ_ID_J] / (Epsilon_14_Sqrt[LJ_ID_J] * Epsilon_List[LJ_ID_J]) );
				}
				else	{
					dE_d_LJEmin[LJ_ID_I] += (dE_d_Emin * Epsilon_Sqrt[LJ_ID_J] / Epsilon_Sqrt[LJ_ID_I]);
					dE_d_LJEmin[LJ_ID_J] += (dE_d_Emin * Epsilon_Sqrt[LJ_ID_I] / Epsilon_Sqrt[LJ_ID_J]);
				}


				dE_d_LJRmin[LJ_ID_I] += dE_d_Rmin;
				dE_d_LJRmin[LJ_ID_J] += dE_d_Rmin;

				E_VDW += dE;
				//printf("E_VDW_Counting: %f\n",E_VDW);
			}
		}
	}


#ifdef MPI_ON
	MPI_Reduce(dE_d_LJEmin, dE_d_LJEmin_Sum, ForceField.n_Rec_LJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(dE_d_LJRmin, dE_d_LJRmin_Sum, ForceField.n_Rec_LJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
	if(ProgID == 0)	{
		for(i=0; i<ForceField.n_Rec_LJ; i++)	{
			dE_d_LJEmin[i] = dE_d_LJEmin_Sum[i];
			dE_d_LJRmin[i] = dE_d_LJRmin_Sum[i];
		}
	}
	MPI_Bcast(dE_d_LJEmin, ForceField.n_Rec_LJ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	MPI_Bcast(dE_d_LJRmin, ForceField.n_Rec_LJ, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
#endif

	return E_VDW;
}

inline double Distance_1D_PBC(double x1, double x2)
{
	double dx;
	dx = x1 - x2;
	if(dx < -L_Half)	{
		dx += L_Box;
	}
	if(dx > L_Half)	{
		dx -= L_Box;
	}
	return dx;
}


void Assign_LJ_Parameters(void)
{
	char szChemName_Local[8][N_LEN_CHEM_NAME];
	double Para_List[MAX_DIH_ITEM*3];
	int i, j, Idx;
	double Sigma_I, Sigma_J, Sigma, Sigma_2, Sigma_5, Sigma_6, Sigma_11, Sigma_12, Epsilon, Epsilon_I, Epsilon_J;
	double Epsilon_14;

	LJtypecount = ForceField.n_Rec_LJ;
	for(i=0; i<LJtypecount; i++)	{
		numAtomsByLJType[i] = 0;
	}

	for(i=0; i<LJtypecount; i++)	{
		Epsilon_I = -ForceField.LJ_Rec[i].para[1];
		Sigma_I = ForceField.LJ_Rec[i].para[2];

		for(j=i; j<LJtypecount; j++)	{
			Epsilon_J = -ForceField.LJ_Rec[j].para[1];
			Sigma_J = ForceField.LJ_Rec[j].para[2];

			Sigma = Sigma_I + Sigma_J;
			Epsilon = sqrt(Epsilon_I * Epsilon_J);
	
			Sigma_2 = Sigma * Sigma;
			Sigma_5 = Sigma_2 * Sigma_2 * Sigma;
			Sigma_6 = Sigma_2 * Sigma_2 * Sigma_2;
			Sigma_11 = Sigma_6 * Sigma_2 * Sigma_2 * Sigma;
			Sigma_12 = Sigma_6 * Sigma_6;

			ATable[i][j] = ATable[j][i] = Epsilon * Sigma_12;
			BTable[i][j] = BTable[j][i] = Epsilon * Sigma_6 * 2.0;

			Sigma_11_Mat[i][j] = Sigma_11_Mat[j][i] = Sigma_11;
			Sigma_12_Mat[i][j] = Sigma_12_Mat[j][i] = Sigma_12;
			Sigma_5_Mat[i][j] = Sigma_5_Mat[j][i] = Sigma_5;
			Sigma_6_Mat[i][j] = Sigma_6_Mat[j][i] = Sigma_6;
			Sigma_Mat[i][j] = Sigma_Mat[j][i] = Sigma;
			Epsilon_IJ_Mat[i][j] = Epsilon_IJ_Mat[j][i] = Epsilon;
		}
	}

	for(i=0; i<LJtypecount; i++)	{
		Epsilon_I = -ForceField.LJ_Rec[i].para[1];
		Epsilon_14 = -ForceField.LJ_Rec[i].para[4];

		// Approximation to decrease the number of variables !!!
		// assume Epsilon_14 == Epsilon .or.  Epsilon_14 == 0.5*Epsilon
		if(Epsilon_14 != Epsilon_I)	{
			Epsilon_14 = 0.5 * Epsilon_I;
		}
		Epsilon_List[i] = Epsilon_I;
		Epsilon_14_List[i] = Epsilon_14;

		Epsilon_Sqrt[i] = sqrt(Epsilon_List[i]);
		Epsilon_14_Sqrt[i] = sqrt(Epsilon_14_List[i]);
	}

	//start	to assign parameters for LJ
	memset(Active_LJ, 0, sizeof(int)*N_REC_MAX);
	for(i=0; i<nAtom; i++)	{
		strcpy(szChemName_Local[0], szChemName[i]);

		Idx = ForceField.GetPara_LJ(szChemName_Local, Para_List);
		List_LJ_Type[i] = Idx;

		if(Idx >= 0)	{
			numAtomsByLJType[Idx]++;
			Active_LJ[Idx] = 1;
		}
		else	{
		}
	}

	nActiveLJ = 0;
	for(i=0; i<LJtypecount; i++)	{
		if(Active_LJ[i])	{
			printf("VDW type: %s\n", ForceField.LJ_Rec[i].Chem);
			nActiveLJ++;
		}
	}
	printf("Total %d VDW type.\n", nActiveLJ);

	//Cal_E_VDW_LRC(); //Chetan: Useless and utter waste of time to call this function here
	//Cal_dE_LRC_d_VDW_Param(); //Chetan: Useless and utter wate of time to call this function here
}

double Cal_E_VDW_LRC(void)
{
	int numLJsites;
        double npairs;
	int i, j, count;
	double sumOfAs, sumOfBs;
	double LJAvgA, LJAvgB, tail_corr_ener, tail_corr_virial;

	numLJsites = 0;
        npairs = 0.0;
	count = 0;
	sumOfAs = sumOfBs = 0.0;

	/*	
	for(i=0; i<LJtypecount; i++)	{
		if(numAtomsByLJType[i])	{
			for (j=0; j<LJtypecount; j++) {
				if(i == j)	{
					sumOfAs += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j] * ATable[i][j];
					sumOfBs += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j] * BTable[i][j];
					count += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j];
				}
				else	{
					sumOfAs += (numAtomsByLJType[i]) * numAtomsByLJType[j] * ATable[i][j];
					sumOfBs += (numAtomsByLJType[i]) * numAtomsByLJType[j] * BTable[i][j];
					count += (numAtomsByLJType[i]) * numAtomsByLJType[j];
				}
			}
		}
	}
	*/

	//Chetan: Previously nAtoms was being used instead of numLJsites to prescale A and B. This is fixed here.
	//Chetan: The L_Box is hardset to 27.2 A (if not called by read_dcd routune)... updating L_Box is safer.
	for (int i=0; i < LJtypecount; i++) {
                for (int j=0; j < LJtypecount; j++) {
                     A = ATable[i][j];
                     B = BTable[i][j];
                     if (!A && !B) continue; // as in NAMD, don't count if both A and B are "0"     
                     if(i==j) {
                        npairs = (numAtomsByLJType[i] - 1)*(numAtomsByLJType[j]);
                        sumOfAs += npairs*A;
                        sumOfBs += npairs*B;
                        count += npairs;
                        numLJsites += numAtomsByLJType[i];
                     }
                     else {
                        npairs = (numAtomsByLJType[i])*(numAtomsByLJType[j]);
                        sumOfAs += npairs*A;
                        sumOfBs += npairs*B;
                        count += npairs;
                     }
                }
        }

	
	LJAvgA = sumOfAs / count;
	LJAvgB = sumOfBs / count;

	double rcut = r_Cutoff;
	double rcut2= rcut * rcut;
	double rcut3= rcut2 * rcut;
	double rcut4= rcut2 * rcut2;
	double rcut5= rcut2 * rcut3;
	double rswitch = r_Switch;
	double rswitch2 = rswitch * rswitch;
	double rswitch3 = rswitch2 * rswitch;
	double rswitch4 = rswitch2 * rswitch2;
	double rswitch5 = rswitch2 * rswitch3;

	//tail_corr_ener = tail_corr_virial = (16*nAtom*nAtom*PI*(-105*LJAvgB*rcut5*rswitch5 + LJAvgA*(3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
	
	///****Chetan: scaling avgA and avgB bu nAtom*nAtom doesn't match NAMD for drudes. Since A and B for drudes are 0, they shouldn't be used for scaling and numLJsites*numLJsites should be used instead.***// 
	tail_corr_ener = tail_corr_virial = (16*numLJsites*numLJsites*PI*(-105*LJAvgB*rcut5*rswitch5 + LJAvgA*(3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
	
	E_LRC = tail_corr_ener/(L_Box*L_Box*L_Box);

	return E_LRC;
}

double Cal_dE_LRC_d_VDW_Param(void)
{
	int numLJsites;
	double npairs; 
	int i, j, count;
	double sumOfAs, sumOfBs;
	double LJAvgA, LJAvgB, tail_corr_ener, tail_corr_virial;
	
	numLJsites = 0;
	npairs = 0.0;
	count = 0;
	sumOfAs = sumOfBs = 0.0;
	
	/*
	for(i=0; i<LJtypecount; i++)	{
		if(numAtomsByLJType[i])	{
			for (j=0; j<LJtypecount; j++) {
				if(i == j)	{
					sumOfAs += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j] * ATable[i][j];
					sumOfBs += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j] * BTable[i][j];
					count += (numAtomsByLJType[i] - 1) * numAtomsByLJType[j];
				}
				else	{
					sumOfAs += (numAtomsByLJType[i]) * numAtomsByLJType[j] * ATable[i][j];
					sumOfBs += (numAtomsByLJType[i]) * numAtomsByLJType[j] * BTable[i][j];
					count += (numAtomsByLJType[i]) * numAtomsByLJType[j];
				}
			}
		}
	}
	*/
	
	//Chetan: Previously nAtoms was being used instead of numLJsites to prescale A and B. This is fixed here.
	//Chetan: If not called from read_dcd routine, L_Box is now set to user fed value.
	for (int i=0; i < LJtypecount; i++) {
    		for (int j=0; j < LJtypecount; j++) {
		     A = ATable[i][j];
                     B = BTable[i][j];	
		     if (!A && !B) continue; // as in NAMD, don't count if both A and B are "0"	
		     if(i==j) {
      			npairs = (numAtomsByLJType[i] - 1)*(numAtomsByLJType[j]);
      			sumOfAs += npairs*A;
      			sumOfBs += npairs*B;
      			count += npairs;
      			numLJsites += numAtomsByLJType[i];
		     }
		     else {
			npairs = (numAtomsByLJType[i])*(numAtomsByLJType[j]);
			sumOfAs += npairs*A;
			sumOfBs += npairs*B;
			count += npairs;		
		     }		
    		}
  	}

	LJAvgA = sumOfAs / count; //matches NAMD
	LJAvgB = sumOfBs / count; //matches NAMD

	if (printLRCFlag==0) {
		printf("LONG-RANGE LJ: AVERAGE A AND B COEFFICIENTS %f %f\n",LJAvgA,LJAvgB);
	}

	double rcut = r_Cutoff;
	double rcut2= rcut * rcut;
	double rcut3= rcut2 * rcut;
	double rcut4= rcut2 * rcut2;
	double rcut5= rcut2 * rcut3;
	double rswitch = r_Switch;
	double rswitch2 = rswitch * rswitch;
	double rswitch3 = rswitch2 * rswitch;
	double rswitch4 = rswitch2 * rswitch2;
	double rswitch5 = rswitch2 * rswitch3;


	//tail_corr_ener = tail_corr_virial = (16*nAtom*nAtom*PI*(-105*LJAvgB*rcut5*rswitch5 + LJAvgA*(3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
	
	///****Chetan: scaling avgA and avgB bu nAtom*nAtom doesn't match NAMD for drudes. Since A and B for drudes are 0, they shouldn't be used for scaling and numLJsites*numLJsites should be used instead.***// 
	tail_corr_ener = tail_corr_virial = (16*numLJsites*numLJsites*PI*(-105*LJAvgB*rcut5*rswitch5 + LJAvgA*(3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));


	if (printLRCFlag==0) {
		printf("TAIL_CORRECTION_ENERGY: %f\n",tail_corr_ener);
	}


	E_LRC = tail_corr_ener/(L_Box*L_Box*L_Box);


	double d_LJAvgA_d_Emin, d_LJAvgA_d_Rmin, d_LJAvgB_d_Emin, d_LJAvgB_d_Rmin;
	double nij_Pair, Coeff_dE_d_LJAvgA, Coeff_dE_d_LJAvgB;


	//Coeff_dE_d_LJAvgA = (16*nAtom*nAtom*PI*( (3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
	//Coeff_dE_d_LJAvgB = (16*nAtom*nAtom*PI*(-105*rcut5*rswitch5))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));

	//Chetan: Applying numLJsites fix here as well to be consistent with NAMD
	Coeff_dE_d_LJAvgA = (16*numLJsites*numLJsites*PI*( (3*rcut4 + 9*rcut3*rswitch + 11*rcut2*rswitch2 + 9*rcut*rswitch3 + 3*rswitch4)))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));
	Coeff_dE_d_LJAvgB = (16*numLJsites*numLJsites*PI*(-105*rcut5*rswitch5))/(315*rcut5*rswitch5*((rcut + rswitch)*(rcut + rswitch)*(rcut + rswitch)));

	Coeff_dE_d_LJAvgA /= (L_Box*L_Box*L_Box);
	Coeff_dE_d_LJAvgB /= (L_Box*L_Box*L_Box);

//	memset(dE_d_LJEmin, 0, sizeof(double)*ForceField.n_Rec_LJ);
//	memset(dE_d_LJRmin, 0, sizeof(double)*ForceField.n_Rec_LJ);

	for(i=0; i<LJtypecount; i++)	{
		if(Active_LJ[i])	{
			d_LJAvgA_d_Emin = d_LJAvgA_d_Rmin = 0.0;
			d_LJAvgB_d_Emin = d_LJAvgB_d_Rmin = 0.0;

			for(j=0; j<LJtypecount; j++)	{
				if(Active_LJ[j])	{
					if(i == j)	{
						nij_Pair = 1.0 * (numAtomsByLJType[i] - 1) * numAtomsByLJType[j] * 2.0;	// accounting for two components in derivative of [ii]
					}
					else	{
						nij_Pair = 1.0 * numAtomsByLJType[i] * numAtomsByLJType[j] * 2.0;	// accounting for ij and ji pair
					}

					if(Epsilon_Sqrt[i] > 1.0E-50)	{
						d_LJAvgA_d_Emin += ( nij_Pair * Sigma_12_Mat[i][j] * 0.5 *  Epsilon_Sqrt[j] / Epsilon_Sqrt[i]);
						d_LJAvgB_d_Emin += ( nij_Pair * Sigma_6_Mat[i][j] * 1.0 *  Epsilon_Sqrt[j] / Epsilon_Sqrt[i]);
					}

					d_LJAvgA_d_Rmin += ( nij_Pair * Epsilon_IJ_Mat[i][j] * 12.0 * Sigma_11_Mat[i][j]);
					d_LJAvgB_d_Rmin += ( nij_Pair * Epsilon_IJ_Mat[i][j] * 12.0 * Sigma_5_Mat[i][j]);
				}
			}

			d_LJAvgA_d_Emin /= count;
			d_LJAvgB_d_Emin /= count;
			d_LJAvgA_d_Rmin /= count;
			d_LJAvgB_d_Rmin /= count;

			dE_d_LJEmin[i] += (Coeff_dE_d_LJAvgA*d_LJAvgA_d_Emin + Coeff_dE_d_LJAvgB*d_LJAvgB_d_Emin);
			dE_d_LJRmin[i] += (Coeff_dE_d_LJAvgA*d_LJAvgA_d_Rmin + Coeff_dE_d_LJAvgB*d_LJAvgB_d_Rmin);
		}
	}
	printLRCFlag = 1; //Chetan: will print A/B coeffs and tail_corr_ener only once  
	return E_LRC;
}

void Accumulate_Gradients_VDW_Param(int iFrame)
{
	int i;

	for(i=0; i<LJtypecount; i++)	{
		dE_d_LJEmin_Acc[i] += dE_d_LJEmin[i];
		dE_d_LJRmin_Acc[i] += dE_d_LJRmin[i];

		E_dE_d_LJEmin_Acc[i] += (E_Liquid[iFrame] * dE_d_LJEmin[i]);
		E_dE_d_LJRmin_Acc[i] += (E_Liquid[iFrame] * dE_d_LJRmin[i]);
	}

//	printf("\n");
//	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
//		if(Active_LJ[i])	{
//			printf("%20.10lf %20.10lf\n", dE_d_LJEmin[i]/nMol,dE_d_LJRmin[i]/nMol);
//		}
//	}
}

void Accumulate_Gradients_V_VDW_Param(void)
{
	int i;

	Volume_Mol_Acc += Volume_Mol;

	for(i=0; i<LJtypecount; i++)	{
		V_dE_d_LJEmin_Acc[i] += (Volume_Mol * dE_d_LJEmin[i]);
		V_dE_d_LJRmin_Acc[i] += (Volume_Mol * dE_d_LJRmin[i]);
	}
}

void Output_Gradients_VDW_Param(int nFrames)
{
	FILE *fOut;
	int i;
	double Volume_Mol_Avg;

	//printf("nFrames %d\n",nFrames);
	//printf("E_Liquid_Avg %f\n",E_Liquid_Avg); 

	fOut = fopen("dE_dvdw_liquid.txt", "w");
	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
		if(Active_LJ[i])	{
//			fprintf(fOut, "%20.10lf %20.10lf\n", dE_d_LJEmin_Acc[i]/(nFrames*nMol),dE_d_LJRmin_Acc[i]/(nFrames*nMol));
			fprintf(fOut, "%-8s %20.10lf %20.10lf\n", ForceField.LJ_Rec[i].RealChem, 
				dE_d_LJEmin_Acc[i]/(nFrames*nMol) - ( E_dE_d_LJEmin_Acc[i]/nFrames - E_Liquid_Avg*dE_d_LJEmin_Acc[i]/nFrames ) / (T_Sim*KBOLTZ*nMol),
				dE_d_LJRmin_Acc[i]/(nFrames*nMol) - ( E_dE_d_LJRmin_Acc[i]/nFrames - E_Liquid_Avg*dE_d_LJRmin_Acc[i]/nFrames ) / (T_Sim*KBOLTZ*nMol) );
		}
//		else	{
//			fprintf(fOut, "%20.10lf %20.10lf\n", 0.0, 0.0);
//		}
	}
	fclose(fOut);

	Volume_Mol_Avg = Volume_Mol_Acc/nFrames;
	fOut = fopen("dV_dvdw_liquid.txt", "w");
	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
		if(Active_LJ[i])	{
			fprintf(fOut, "%-8s %20.10lf %20.10lf\n", ForceField.LJ_Rec[i].RealChem, 
				Volume_Mol_Avg*dE_d_LJEmin_Acc[i]/nFrames - V_dE_d_LJEmin_Acc[i]/nFrames, 
				Volume_Mol_Avg*dE_d_LJRmin_Acc[i]/nFrames - V_dE_d_LJRmin_Acc[i]/nFrames);
		}
//		else	{
//			fprintf(fOut, "%20.10lf %20.10lf\n", 0.0, 0.0);
//		}
	}
	fclose(fOut);

}

void Determine_Mol_Size(void)	// assume neat liquid only has one component
{
	int ResFirst;

	ResFirst = ResID[0];

	MolSize = 0;
	while(ResID[MolSize] == ResFirst)	{
		MolSize++;
	}

	nMol = nAtom / MolSize;
}

