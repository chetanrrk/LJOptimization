#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ff-gas.h"
#include <dirent.h>

//#define MPI_ON  0

//#ifdef MPI_ON
//#include <mpi.h>
//#endif

//start DCD data
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
//end   DCD data

float A, B, C, alpha_loc, beta, my_gamma;
dcdhandle *dcd;


void findDCDFiles(); //new
void Read_Dcd(char szName[]); //new
static int read_dcdheader(FILE *fd, int *N, int *NSET, int *ISTART,
                   int *NSAVC, double *DELTA, int *NAMNF, int **FREEINDEXES, float **fixedcoords, int *reverseEndian, int *charmm);
static int read_next_timestep(void *v, int natoms); //new
void Free_Memory(void); //new
static int read_charmm_extrablock(FILE *fd, int charmm, int reverseEndian, float *unitcell);



void Quit_With_Error_Msg(char szMsg[]);
void ReadSnapshot(void);
void ReadNamdEne(char szName[]); // reads pot ene. of namd

FILE *fFile_Run_Log;	// will be shared by other source code

CMol Mol;
CForceField ForceField;

int ProgID=0, nProc=1;

int nFrame_Skip;
int nMol, nAtom, nTotal;
char szName_FF[256], szName_Psf[256], szName_Snapshot[256];
char szName_Data[256];
double *px=NULL, *py=NULL, *pz=NULL;
double L_Box=100.0;

#define MAX_N_Mols     (500) // new
char dcds[MAX_N_Mols][256]; // new; stores names of dcd files
double *potEne = new double [20001]; //stores all the NAMD pot energies

int main(int argc, char *argv[])
{
	double T_Sim, E_Mean, Time_Step;
	int i, nSteps, IdxMol, Idx;
	FILE *fOut;
/*
#ifdef MPI_ON
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProgID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
#endif
*/
	if(argc != 8)	{
		printf("Usage: gas_sim ff.str mol.xpsf T_Sim MD_Step Time_Step(ps) L_Box nFrame_Skip \n");
/*
#ifdef MPI_ON
		MPI_Finalize();
#endif
*/
		exit(1);
	}

	strcpy(szName_FF, argv[1]);
	strcpy(szName_Psf, argv[2]);
	T_Sim = atof(argv[3]);
	nSteps = atoi(argv[4]);
	Time_Step = atof(argv[5]);
	L_Box = atof(argv[6]);
	nFrame_Skip = atoi(argv[7]);

	fFile_Run_Log = stdout;

	ForceField.ReadForceField(szName_FF);

	Mol.ReadPSF(szName_Psf, 0);
	Mol.Transfer_Drude_Mass();
	nAtom = Mol.nAtom;

	Mol.AssignForceFieldParameters(&ForceField);

	memset(Mol.dE_d_LJEmin_Acc, 0, sizeof(double)*N_REC_MAX);
	memset(Mol.dE_d_LJRmin_Acc, 0, sizeof(double)*N_REC_MAX);
	memset(Mol.E_dE_d_LJEmin_Acc, 0, sizeof(double)*N_REC_MAX);
	memset(Mol.E_dE_d_LJRmin_Acc, 0, sizeof(double)*N_REC_MAX);
	//printf("******here0******\n");

	/****************************************
	//Chetan
	//MD is run separately with NAMD 
	//MD is run by randomly selecting "N" molecules and checking if enthalpy converge...increase "N" until convergence 
	//For each mol in the box, dcd is loaded and  Mol obj created using  mol in every "n" MD steps
	//For each mol in the box, every "n" steps of MD, the vdw grads wrt sigma and epsilon are computed
	//Gradients are accumulated and printed to a File 
	******************************************/

	//initializing dcds arrays
        for(int i=0; i<sizeof(dcds)/sizeof(dcds[0]); i++){
                strcpy(dcds[i],"null");
        }

	//printf("******here1******\n");

	findDCDFiles(); // finds all dcd files available from NAMD run
	//printf("******here2******\n");
	
	//Call this for each dcd to get snapshot energies
	//ReadNamdEne(); // reads the stored pot. energies from  NAMD run

	// This for loop has replaced the old for loop of Lei's MD
        for (int i=ProgID; i<sizeof(dcds)/sizeof(dcds[0]); i+=nProc){
                if(strcmp("null",dcds[i]) != 0){
			ReadNamdEne(dcds[i]); // reads the stored pot. energies from  NAMD run for each dcd

			char myDcd[256];
			strcpy(myDcd,dcds[i]);
			strcat(myDcd,".dcd");

			//printf("******here3******\n");

		        Read_Dcd(myDcd); //read coords from each snapshot of each dcd
			//printf("******here4******\n");
                        Free_Memory();

			// reinitializing the pot energies array for each dcd
			for (int i=0; i< 20001; i++){
			    potEne[i]=0.0;
			}
		}
	}

	//printf("all done!");	
	
	//exit(0); // looks good until here
/*
#ifdef MPI_ON
	//Chetan: The Pot. energy is now read from NAMD logs and averaged
	printf("running mpi");
	if(ProgID == 0)	{
		FILE *fIn;
		char buff[255];
		fIn = fopen("data_E_gas.txt", "r");	
		fscanf(fIn,"%s",buff);
		double E_Mean = atof(buff);
		//printf("Energy: %f\n", E_Mean);
		fclose(fIn);
	}
	
	
	//exit(0);
	

	double *pData_dE_dEmin, *pData_dE_dRmin, *pData_dE_dEmin_Sum, *pData_dE_dRmin_Sum;
	double *pData_E_dE_dEmin, *pData_E_dE_dRmin, *pData_E_dE_dEmin_Sum, *pData_E_dE_dRmin_Sum;

	pData_dE_dEmin = new double[Mol.nActiveLJ];
	pData_dE_dRmin = new double[Mol.nActiveLJ];
	pData_dE_dEmin_Sum = new double[Mol.nActiveLJ];
	pData_dE_dRmin_Sum = new double[Mol.nActiveLJ];

	pData_E_dE_dEmin = new double[Mol.nActiveLJ];
	pData_E_dE_dRmin = new double[Mol.nActiveLJ];
	pData_E_dE_dEmin_Sum = new double[Mol.nActiveLJ];
	pData_E_dE_dRmin_Sum = new double[Mol.nActiveLJ];

	Idx = 0;
	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
		if(Mol.Active_LJ[i])	{
			// dE_d_LJEmin_Acc are summed over all molecule copies and snapsots
			// nAcc_dE_d_LJ is mol copies * MDsteps (i.e. nMols * MDSteps/call_made_for_derivative_calc)
			// pData_dE_dEmin etc are average over dE_d_LJEmin_Acc scaled by nAcc_dE_d_LJ

			//printf("%s\n",ForceField.LJ_Rec[i].RealChem);
			//printf("%s%i %f %f \n","Mol.dE_d_LJEmin_Acc",i,Mol.dE_d_LJEmin_Acc[i],Mol.nAcc_dE_d_LJ);
			//printf("%s%i %f %f \n","Mol.dE_d_LJRmin_Acc",i,Mol.dE_d_LJRmin_Acc[i],Mol.nAcc_dE_d_LJ);
			//printf("%s%i %f %f \n","Mol.E_dE_d_LJEmin_Acc",i,Mol.E_dE_d_LJEmin_Acc[i],Mol.nAcc_dE_d_LJ);
			//printf("%s%i %f %f \n","Mol.E_dE_d_LJRmin_Acc",i,Mol.E_dE_d_LJRmin_Acc[i],Mol.nAcc_dE_d_LJ);
		
			pData_dE_dEmin[Idx] = Mol.dE_d_LJEmin_Acc[i]/Mol.nAcc_dE_d_LJ;
			pData_dE_dRmin[Idx] = Mol.dE_d_LJRmin_Acc[i]/Mol.nAcc_dE_d_LJ;

			pData_E_dE_dEmin[Idx] = Mol.E_dE_d_LJEmin_Acc[i]/Mol.nAcc_dE_d_LJ;
			pData_E_dE_dRmin[Idx] = Mol.E_dE_d_LJRmin_Acc[i]/Mol.nAcc_dE_d_LJ;

			Idx++;
		}
	}

	//pData_dE_dEmin_Sum is sum of pData_dE_dEmin etc across number of processors
	MPI_Reduce(pData_dE_dEmin, pData_dE_dEmin_Sum, Mol.nActiveLJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(pData_dE_dRmin, pData_dE_dRmin_Sum, Mol.nActiveLJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(pData_E_dE_dEmin, pData_E_dE_dEmin_Sum, Mol.nActiveLJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(pData_E_dE_dRmin, pData_E_dE_dRmin_Sum, Mol.nActiveLJ, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

	//Master computes the average of accumulated derivatives scaled by the number of processors used 
	
	//exit(0);

	if(ProgID == 0)	{
		fOut = fopen("dE_dvdw_gas.txt", "w");
		Idx = 0;
		for(i=0; i<ForceField.n_Rec_LJ; i++)	{
			if(Mol.Active_LJ[i])	{
                        	//printf("%s\n",ForceField.LJ_Rec[i].RealChem);
                        	//printf("pData_dE_dEmin_Sum %f\n",pData_dE_dEmin_Sum[Idx]);
                        	//printf("pData_E_dE_dEmin_Sum %f\n",pData_E_dE_dEmin_Sum[Idx]);
                        	//printf("pData_dE_dRmin_Sum %f\n",pData_dE_dRmin_Sum[Idx]);
                        	//printf("pData_E_dE_dRmin_Sum %f\n",pData_E_dE_dRmin_Sum[Idx]);
				//printf("E_Mean %f\n",E_Mean);				

				fprintf(fOut, "%-8s %20.10lf %20.10lf\n", ForceField.LJ_Rec[i].RealChem, 
					pData_dE_dEmin_Sum[Idx]/nProc - (pData_E_dE_dEmin_Sum[Idx]/nProc - E_Mean*pData_dE_dEmin_Sum[Idx]/nProc)/(T_Sim*KBOLTZ), 
					pData_dE_dRmin_Sum[Idx]/nProc - (pData_E_dE_dRmin_Sum[Idx]/nProc - E_Mean*pData_dE_dRmin_Sum[Idx]/nProc)/(T_Sim*KBOLTZ));

				Idx++;
			}
		}
		fclose(fOut);
	}


	
	//exit(0);

	delete []pData_dE_dEmin;
	delete []pData_dE_dRmin;
	delete []pData_dE_dEmin_Sum;
	delete []pData_dE_dRmin_Sum;

	delete []pData_E_dE_dEmin;
	delete []pData_E_dE_dRmin;
	delete []pData_E_dE_dEmin_Sum;
	delete []pData_E_dE_dRmin_Sum;

#else
*/

	//printf("running on a single proc");

	FILE *fIn;
        char buff[255];
        fIn = fopen("../../data_E_gas.txt", "r");     
        fscanf(fIn,"%s",buff);
        //double E_Mean = atof(buff);
        E_Mean = atof(buff);
        //printf("Energy: %f\n", E_Mean);
        fclose(fIn);

	//exit(0);


	fOut = fopen("dE_dvdw_gas.txt", "w");
	for(i=0; i<ForceField.n_Rec_LJ; i++)	{
		if(Mol.Active_LJ[i])	{
			//printf("%s\n",ForceField.LJ_Rec[i].RealChem);
			//printf("dE_d_LJEmin_Acc %f\n",Mol.dE_d_LJEmin_Acc[i]);
			//printf("E_dE_d_LJEmin_Acc %f\n",Mol.E_dE_d_LJEmin_Acc[i]);
			//printf("dE_d_LJRmin_Acc %f\n",Mol.dE_d_LJRmin_Acc[i]);
			//printf("E_dE_d_LJRmin_Acc %f\n",Mol.E_dE_d_LJRmin_Acc[i]);

                        //printf("%s %s %s %s %f %f %d\n",ForceField.LJ_Rec[i].RealChem,"Mol.dE_d_LJEmin_Acc[i]","Mol.dE_d_LJRmin_Acc[i]","Mol.nAcc_dE_d_LJ",Mol.dE_d_LJEmin_Acc[i],Mol.dE_d_LJRmin_Acc[i],Mol.nAcc_dE_d_LJ);
                        //printf("%s %s %f\n",ForceField.LJ_Rec[i].RealChem,"Mol.E_dE_d_LJEmin_Acc[i]",Mol.E_dE_d_LJEmin_Acc[i]);
                        //printf("%s %s %f\n",ForceField.LJ_Rec[i].RealChem,"Mol.E_dE_d_LJRmin_Acc[i]",Mol.E_dE_d_LJRmin_Acc[i]);
	
			fprintf(fOut, "%-8s %20.10lf %20.10lf\n", ForceField.LJ_Rec[i].RealChem, 
				Mol.dE_d_LJEmin_Acc[i]/Mol.nAcc_dE_d_LJ - (Mol.E_dE_d_LJEmin_Acc[i]/Mol.nAcc_dE_d_LJ - E_Mean*Mol.dE_d_LJEmin_Acc[i]/Mol.nAcc_dE_d_LJ)/(T_Sim*KBOLTZ), 
				Mol.dE_d_LJRmin_Acc[i]/Mol.nAcc_dE_d_LJ - (Mol.E_dE_d_LJRmin_Acc[i]/Mol.nAcc_dE_d_LJ - E_Mean*Mol.dE_d_LJRmin_Acc[i]/Mol.nAcc_dE_d_LJ)/(T_Sim*KBOLTZ));
		}
	}
	fclose(fOut);
//#endif

	if(px)	delete []px;
	if(py)	delete []py;
	if(pz)	delete []pz;
	delete Mol.myrand;
	delete []potEne;	
/*
#ifdef MPI_ON
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
*/
	return 0;
}

void ReadNamdEne(char szName[]){
     FILE *fIn;
     char buff[255];
     char *prefix;
     prefix = strtok(szName,".");
     char fName[80];
     strcpy(fName,"../logs/E_namd_dcd_");    
     strcat(fName,prefix);
     strcat(fName,".txt");
     fIn = fopen(fName,"r");
     if(fIn==NULL){
        printf("Fail to open file: %s\nQuit\n",fName);
        exit(1);
     }
     int lineCount=-1; // we want to skip Energy at step 0s since, in namd12, energy at step 1 correspnds to dcd at idx 1
     while(fgets(buff,sizeof(buff),fIn)){
        char* err;
        if(feof(fIn)!=0){
                printf("Invalid energy value: %s\nQuit\n","NULL Found");
                exit(1);
        }
        else {
                double pot ; //= strtod(buff, &err);
		pot = atof(buff);
		if (pot==0){
		   break;	
		}
                else {
			if (lineCount >= 0) {
				potEne[lineCount] = pot; 
			}
		}
                //if(potEne[lineCount]==0){
                //if(pot==0){
                //        printf("Invalid energy value: %s\nQuit\n","Junk Value");
                //        exit(1);
                //}
        }
        lineCount++;
     }
     fclose(fIn);	

     /*for(int i=0; i<20001;i++){
	printf("pot %i stored %f\n",i, potEne[i]);
     }*/	

     /*int i=0;
     while(true){
         if(potEne[i]==0){
                printf("num of vals: %d",i);
                break;
         }
         else{
                printf("%f,",potEne[i]);
         }

         i++;
     }
     printf("\n");*/

}


//Chetan: Reads the dcds automatically and grabs the coords
void findDCDFiles(){
   DIR *dp;
   struct dirent *ep;
   const char* extension = ".dcd";
   int i = 0; //dcd count       
   dp = opendir ("./");
   if (dp != NULL)
    { 
      while (ep = readdir (dp))
        for (char*p = strtok(ep->d_name,".");p != NULL; p = strtok(NULL, ".") )
        {
           if (strcmp(p,"dcd") == 0){
                char file_name[256];
                strcpy(file_name,ep->d_name);
                strcat(file_name,".dcd");
                strcpy(dcds[i],file_name);
                i+=1;     
           }    
        }  
      (void) closedir (dp);
    } 
   else
    perror ("Couldn't open the directory");
    
}


//Chetan: Reads each dcd and grabs the coords
void Read_Dcd(char szName[])
{
        int rc, iFrame, ReadItem;
        FILE *fd, *fIn;
        fd = fopen(szName, "rb");

        //printf("=====reading dcd: %s======\n",szName);
        if(fd == NULL)  {
                Quit_With_Error_Msg("Fail to open dcd file.\nQuit.\n");
        }

        dcd = new dcdhandle;
        memset(dcd, 0, sizeof(dcdhandle));
        dcd->fd = fd;

	
	//printf("============TEST1==========\n");


        if ((rc = read_dcdheader(dcd->fd, &dcd->natoms, &dcd->nsets, &dcd->istart,
                &dcd->nsavc, &dcd->delta, &dcd->nfixed, &dcd->freeind,
                &dcd->fixedcoords, &dcd->reverse, &dcd->charmm))) {
                fclose(dcd->fd);
                free(dcd);

                Quit_With_Error_Msg("Error in reading dcd head file.\n");
        }

        dcd->first = 1;
        nTotal = dcd->natoms;

        if( nTotal != nAtom )   {
                Quit_With_Error_Msg("nTotal != nAtom\nQuit!\n");
        }


	//printf("============TEST2==========\n");

        dcd->x = new float[dcd->natoms];
        dcd->y = new float[dcd->natoms];
        dcd->z = new float[dcd->natoms];

        if (!dcd->x || !dcd->y || !dcd->z) {
                fclose(dcd->fd);
                Free_Memory();
                Quit_With_Error_Msg("Unable to allocate space for reading dcd file.\n");
        }


	//printf("============TEST3==========\n");


        //for(iFrame=0; iFrame<nFrame_Skip; iFrame++)     {
        //        read_next_timestep(dcd, dcd->natoms);
        //}

	//printf("dcd->nsets %i\n",dcd->nsets);
	
	
	//printf("============TEST4==========\n");
	iFrame=0;
	//printf("iFrame %i\n",iFrame);

        for(; iFrame<dcd->nsets; iFrame++)      {
                read_next_timestep(dcd, dcd->natoms);
		
		for(int ic=0; ic<dcd->natoms;ic++){
			Mol.x[ic] = dcd->x[ic];
			Mol.y[ic] = dcd->y[ic];
			Mol.z[ic] = dcd->z[ic];
			//printf("%s %f ","x ",dcd->x[ic]);
			//printf("%s %f ","x ",Mol.x[ic]);
			//printf("%s %f ","y ",dcd->y[ic]);
			//printf("%s %f ","y ",Mol.y[ic]);
			//printf("%s %f \n","z ",dcd->z[ic]);
			//printf("%s %f \n","z ",Mol.z[ic]);
		}
		

		//printf("iFrame %i\n",iFrame);
		Mol.E_Total = potEne[iFrame]; // read stored pot. ene from NAMD md
		//printf("Mol.E_Total: %f\n",Mol.E_Total);

		//TODO: Check this derivative here between lei and my calculation for additive ff example		
		//printf("============TEST==========\n");
		Mol.Cal_dE_d_VDW_Param(); //derivatives computed here

        }
	//exit(0);
        fclose(fd);

}

//Chetan: Copied as implemented in Charmm (borrowed from Lei's code) 
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
//  printf("dcdplugin) detected standard 32-bit DCD file of native endianness\n");
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



//Chetan: Copied as implemented in Charmm (borrowed from Lei's code) 
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
                if(ReadItem != rec_scale)       return DCD_BADREAD;
                ReadItem = fread(X, sizeof(float), N, fd);
                if(ReadItem != N)       return DCD_BADREAD;
                ReadItem = fread(&tmpbuf[1*rec_scale], sizeof(int), rec_scale*2, fd);
                if(ReadItem != rec_scale*2)     return DCD_BADREAD;
                ReadItem = fread(Y, sizeof(float), N, fd);
                if(ReadItem != N)       return DCD_BADREAD;
                ReadItem = fread(&tmpbuf[3*rec_scale], sizeof(int), rec_scale*2, fd);
                if(ReadItem != rec_scale*2)     return DCD_BADREAD;
                ReadItem = fread(Z, sizeof(float), N, fd);
                if(ReadItem != N)       return DCD_BADREAD;
                ReadItem = fread(&tmpbuf[5*rec_scale], sizeof(int), rec_scale, fd);
                if(ReadItem != rec_scale)       return DCD_BADREAD;

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
//              ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
//              if (ret_val) return ret_val;
        } else {
                /* if there are fixed atoms, and this isn't the first frame, then we */
                /* only read in the non-fixed atoms for all subsequent timesteps.    */
//              ret_val = read_charmm_extrablock(fd, charmm, reverseEndian, unitcell);
//              if (ret_val) return ret_val;
//              ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//              fixedcoords, fixedcoords+3*N, X, charmm);
//              if (ret_val) return ret_val;
//              ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//              	fixedcoords+N, fixedcoords+3*N, Y, charmm);
//              if (ret_val) return ret_val;
//              ret_val = read_fixed_atoms(fd, N, N-num_fixed, indexes, reverseEndian,
//                      fixedcoords+2*N, fixedcoords+3*N, Z, charmm);
//              if (ret_val) return ret_val;
//              ret_val = read_charmm_4dim(fd, charmm, reverseEndian);
//              if (ret_val) return ret_val;
        }

        return DCD_SUCCESS;
}


//Chetan: Copied as implemented in Charmm (borrowed from Lei's code) 
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
        if(dcd) {
                if (dcd->x)     delete [](dcd->x);
                if (dcd->y)     delete [](dcd->y);
                if (dcd->z)     delete [](dcd->z);
                delete dcd;
        }
}


//Chetan: Copied as implemented in Charmm (borrowed from Lei's code) 
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


//Chetan: Kept intact
void Quit_With_Error_Msg(char szMsg[])
{
	fprintf(fFile_Run_Log, "%s", szMsg);
	fflush(fFile_Run_Log);
	exit(1);
}

//Chetan: Kept intact...can be useful if we ever have to read pdbs
void ReadSnapshot(void)
{
	FILE *fIn;
	int nLine=0, ReadItem;
	char szLine[256], *ReadLine;
	double fTmp=0.0;

	fIn = fopen(szName_Snapshot, "r");

	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szName_Snapshot);
		exit(1);
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if(strncmp(szLine, "ATOM", 4)==0)	{
			nLine++;
		}
	}

	if(nLine%nAtom)	{
		printf("There are %d atoms in %s.\nThere are %.2lf molecules.\n", nLine, szName_Snapshot, 1.0*nLine/nAtom);
		fclose(fIn);
		exit(1);
	}

	px = new double[nLine];
	py = new double[nLine];
	pz = new double[nLine];

	fseek(fIn, 0, SEEK_SET);

	nLine = 0;
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine == NULL)	{
			break;
		}
		if(strncmp(szLine, "ATOM", 4)==0)	{
			ReadItem = sscanf(szLine+30, "%lf%lf%lf", &(px[nLine]), &(py[nLine]), &(pz[nLine]));
			if(ReadItem != 3)	{
				printf("Error in reading file: %s\n%s\n", szName_Snapshot, szLine);
				fclose(fIn);
				exit(1);
			}
			nLine++;
		}
	}

	nMol = nLine/nAtom;
	if(ProgID == 0)	printf("%d molecules in %s.\n", nMol, szName_Snapshot);

	fclose(fIn);
}

//Chetan: kept intact if we ever have to check our structures
void Iterative_Check_BondLength(int Idx, int pVisited[])
{
	int i, Idx_j;
	double *px, *py, *pz, L_Half, dx, dy, dz;

	px = Mol.x;
	py = Mol.y;
	pz = Mol.z;

	L_Half = 0.5 * L_Box;

	for(i=0; i<Mol.AtomBond[Idx]; i++)	{
		Idx_j = Mol.Bond_Array[Idx][i];
		if(pVisited[Idx_j] == 0)	{
			pVisited[Idx_j] = 1;

			dx = px[Idx_j] - px[Idx];
			dy = py[Idx_j] - py[Idx];
			dz = pz[Idx_j] - pz[Idx];

			if(dx > L_Half)	{
				px[Idx_j] -= L_Box;
			}
			else if(dx < -L_Half)	{
				px[Idx_j] += L_Box;
			}

			if(dy > L_Half)	{
				py[Idx_j] -= L_Box;
			}
			else if(dy < -L_Half)	{
				py[Idx_j] += L_Box;
			}

			if(dz > L_Half)	{
				pz[Idx_j] -= L_Box;
			}
			else if(dz < -L_Half)	{
				pz[Idx_j] += L_Box;
			}

			Iterative_Check_BondLength(Idx_j, pVisited);
		}
	}
}

