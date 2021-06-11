#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ff-gas.h"
            
#define PI	(3.14159265358979323846)
#define PI2	(6.28318530717958647692)
#define radian	(57.29577951308232088)
#define radianInv	(0.017453292519943295)

#define MAX_LEN_LINE	(196)
#define ITEM_IN_LINE_8	(8)
#define ITEM_IN_LINE_9	(9)

#ifndef min(a,b)
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_dvector(double *v, int nl, int nh);



//start	to data and subroutines for CMap
#define CMAP_DIM	(24)

double mctp[4][CMAP_DIM][CMAP_DIM];	//24X24 data from parameter / toppar file
void setcmap(int num, int xm, double gmap[][CMAP_DIM][CMAP_DIM], double dx);
void gcstup2(double ggrd[][CMAP_DIM][CMAP_DIM], double ty[4], double ty1[4], double ty2[4], double ty12[4], int ip1, int ip2);
void gcscf(double ty[4], double ty1[4], double ty2[4], double ty12[4], double gres1, double gres2, double tc[][4]);
double E_Phi_Psi_QM[CMAP_DIM*CMAP_DIM], E_Phi_Psi_MM[CMAP_DIM*CMAP_DIM];
//end	to data and subroutines for CMap
extern FILE *fFile_Run_Log;	// will be shared by other source code
extern void Quit_With_Error_Msg(char szMsg[]);

#define NONBOND_CUTOFF	(150.0)
#define FAR_DIST		(NONBOND_CUTOFF+250.0)
double Dist_Cutoff_SQ=NONBOND_CUTOFF*NONBOND_CUTOFF;

void ReadIntDataFromBuff(char szBuff[], int Data[], int n);
void ReadDoubleDataFromBuff(char szBuff[], double Data[], int n);

int FindString(char szBuff[], char szTag[]);
int FindString(char szBuff[], char szTag[])
{
	int Len_Buff, Len_Tag, iPos, iPos_Max;

	Len_Tag = strlen(szTag);
	Len_Buff = strlen(szBuff);

	iPos_Max = Len_Buff-Len_Tag;

	for(iPos=0; iPos<iPos_Max; iPos++)	{
		if(strncmp(szBuff+iPos, szTag, Len_Tag) == 0)	{
			return iPos;
		}
	}
	return (-1);
}

void CMol::ReadPSF(char szNamePsf[], int MultiSegment)
{
        FILE *fIn;
        char szLine[256], szTmp[256], szMolNameLocal[8], DrudeFlag, szErrorMsg[256];
        int i, j, n_Max, iTmp, ReadItem, n_Item, n_Line, n_Left, Count;
        int n_LP_Host[MAX_ATOM], Pointer_LP[MAX_ATOM], LP_Host_List[MAX_ATOM*MAX_LP_POS_HOST], SegID_Last;
        int n_LPH_Host[MAX_ATOM], Pointer_LPH[MAX_ATOM];
        double dTmp=0.0;
        char *token; // stores xpsf token with LPHs host info 
        int LPHIdx;

        nAtom = nBond = nAngle = nDihedral = nImpro = nLP = nLPH = nDrude = nAniso = 0;

        myForceField = NULL;
        memset(IsDrudeHost, 0, sizeof(int)*MAX_ATOM);

        StartLast = EndLast = -1;       // invalid

        fIn = fopen(szNamePsf, "r");
        if(fIn == NULL) {
//              printf("CMol::ReadPSF> Fail to open file %s\nQuit\n", szNamePsf);
                fprintf(fFile_Run_Log, "CMol::ReadPSF> Fail to open file %s\nQuit\n", szNamePsf);
                fflush(fFile_Run_Log);
                exit(1);
        }

        while(1)        {
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for atom information in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NATOM") >= 0)   {
                        break;
                }
        }
        sscanf(szLine, "%d", &nAtom);   //atom number

        if(nAtom > MAX_ATOM)    {
                fprintf(fFile_Run_Log, "nAtom = %d\nnAtom> MAX_ATOM\nQuit\n", nAtom);
                fflush(fFile_Run_Log);
                exit(1);
        }

        nDrude = 0;
        for(i=0; i<nAtom; i++)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadItem = sscanf(szLine, "%d %s %d %s %s %s %lf %lf %d %lf %lf",
                        &iTmp, szMolNameLocal, &(SegID[i]), ResName[i], AtomName[i], ChemName[i], &(CG[i]),
                        &(mass[i]), &(IsFixed[i]), &(alpha[i]), &(thole[i]));

                if(ReadItem == 9)       {
                        fprintf(fFile_Run_Log, "The format of xpsf is not supported. It is supposed to have 11 columns for atom information, with parameters of drude model included.\n");
                        fflush(fFile_Run_Log);
                        exit(1);
                }

                if(ReadItem != 11)      {
                        fprintf(fFile_Run_Log, "Error in reading psf file: %s\n%s\nQuit\n", szNamePsf, szLine);
                        fflush(fFile_Run_Log);
                        exit(1);
                }

                if( (mass[i] < 1.0) && (AtomName[i][0]=='D') )  {
                        IsDrude[i] = 1;
                        IsDrudeHost[i-1] = 1;
                        DrudeList[nDrude] = i;
                        nDrude++;
                }
                else    {
                        IsDrude[i] = 0;
                }

                if(IsFixed[i] == -1)    {       //the flag for lone pair
                        if(strcmp(ChemName[i],"LPH")==0)        {
                                IsLPH[i] = 1;
                                LPHList[nLPH] = i;
                                nLPH++;
                        }
                        else    {
                                IsLonePair[i] = 1;
                                LPList[nLP] = i;
                                nLP++;
                        }
                }
                else    {
                        IsLonePair[i] = 0;
                }
        }
//      strcpy(MolName, szMolNameLocal);        //the name of this molecule

        /*if(To_Delete_Last_Atom) {
                nAtom--;        // to delete the last atom. Designed for delete the test charge
        }*/

        //start to identify the range of the last segment
        if(MultiSegment)        {
                EndLast = nAtom - 1;
                SegID_Last = SegID[EndLast];
                for(i=EndLast; i>=0; i--)       {
                        if(SegID[i] != SegID_Last)      {
                                break;
                        }
                }
                if(i<0) {       // There is only one segment
                        sprintf(szErrorMsg, "Error in reading xpsf file: %s\nIt is expected that there are more than one segment for molecule-molecule interactions.\nQuit\n", szNamePsf);
                        Quit_With_Error_Msg(szErrorMsg);
                }
                else    {
                        StartLast = i+1;
                }
        }
        //end   to identify the range of the last segment

        while(1)        {       //to find the entry for bond list
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for bond list in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NBOND") >= 0)   {
                        break;
                }
        }
        sscanf(szLine, "%d", &nBond);   //number of bonds

        if( (nDrude > 0) && (nBond < nDrude) )  {
                fprintf(fFile_Run_Log, "Error in xpsf file.\n(nDrude > 0) && (nBond < nDrude)\nQuit\n");
                fflush(fFile_Run_Log);
                exit(1);
        }

        n_Item = 2*nBond;
        n_Line = n_Item/ITEM_IN_LINE_8;
        n_Left = n_Item - n_Line*ITEM_IN_LINE_8;
        Count = 0;
        for(i=0; i<n_Line; i++) {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, BondList+Count, ITEM_IN_LINE_8);
                Count += ITEM_IN_LINE_8;
        }
        if(n_Left > 0)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, BondList+Count, n_Left);
                Count += n_Left;
        }

        while(1)        {       //to find the entry for angle list
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for angle list in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NTHETA") >= 0)  {
                        break;
                }
        }
        sscanf(szLine, "%d", &nAngle);  //number of angles

        n_Item = 3*nAngle;
        n_Line = n_Item/ITEM_IN_LINE_9;
        n_Left = n_Item - n_Line*ITEM_IN_LINE_9;
        Count = 0;
        for(i=0; i<n_Line; i++) {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, AngleList+Count, ITEM_IN_LINE_9);
                Count += ITEM_IN_LINE_9;
        }
        if(n_Left > 0)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, AngleList+Count, n_Left);
                Count += n_Left;
        }

        while(1)        {
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for dihedral list in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NPHI:") >= 0)   {
                        break;
                }
        }
        ReadItem = sscanf(szLine, "%d", &nDihedral);    //number of dihedrals

        if(nDihedral >= MAX_DIHEDRAL)   {
                fprintf(fFile_Run_Log, "nDihedral = %d\nnDihedral >= MAX_DIHEDRAL\nQuit.\n", nDihedral);
                fflush(fFile_Run_Log);
                exit(1);
        }

        n_Item = 4*nDihedral;
        n_Line = n_Item/ITEM_IN_LINE_8;
        n_Left = n_Item - n_Line*ITEM_IN_LINE_8;
        Count = 0;
        for(i=0; i<n_Line; i++) {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, DihedralList+Count, ITEM_IN_LINE_8);
                Count += ITEM_IN_LINE_8;
        }
        if(n_Left > 0)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, DihedralList+Count, n_Left);
                Count += n_Left;
        }

        while(1)        {
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for improper dihedral list in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NIMPHI") >= 0)  {
                        break;
                }
        }
        ReadItem = sscanf(szLine, "%d", &nImpro);       //number of improper dihedrals

        n_Item = 4*nImpro;
        n_Line = n_Item/ITEM_IN_LINE_8;
        n_Left = n_Item - n_Line*ITEM_IN_LINE_8;
        Count = 0;
        for(i=0; i<n_Line; i++) {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, ImprDihedralList+Count, ITEM_IN_LINE_8);
                Count += ITEM_IN_LINE_8;
        }
        if(n_Left > 0)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, ImprDihedralList+Count, n_Left);
                Count += n_Left;
        }

        while(1)        {       //to find the entry for lone pair 
                if(feof(fIn))   {
                        Quit_With_Error_Msg("Fail to find the entry for lone pair list in xpsf file.\nQuit\n");
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(FindString(szLine, "!NUMLP") >= 0)   {       // find the entry for lone pair list
                        break;
                }
        }

        //start to read the host atoms information of lone pairs
        int totalLP = nLP+nLPH;
        int lpHost[totalLP];
        int pointerLP[totalLP];  // stores pointer and assing to LP or LPH
        double numLinesToRead=0; // to read LPs/LPHs indices and host info
        int maxIdx = totalLP*4; // assuming all LPs are relative or bisector
        int tmpIndices[maxIdx]; // stores indices records from PSF      
        int idxCount=0; // keeps count of LP and host indices
        int nextLPIdx; // next LP index in pointerLP list
        int lphIdx=0;  // colinear lone pairs counter
        int lpIdx=0;   // relative/bisector lone pairs counter 

        double dist,angle,dihed; // temporarily stores info

        for(i=0; i<totalLP; i++)        {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadItem = sscanf(szLine, "%d %d %s %lf %lf %lf",
                        &lpHost[i], &pointerLP[i], szTmp,&dist,&angle,&dihed);

                if(ReadItem != 6)       {
//                      printf("Error in reading lone pair information.\n%s\nQuit.\n", szLine);
                        fprintf(fFile_Run_Log, "Error in reading lone pair information.\n%s\nQuit.\n", szLine);
                        fflush(fFile_Run_Log);
                        exit(1);
                }

                if(lpHost[i] == 3)      {       //LP Relative or Bisector. 
                        n_LP_Host[lpIdx] = lpHost[i];
                        Pointer_LP[lpIdx] = pointerLP[i];
                        LP_Dist[lpIdx] = dist;
                        LP_Theta[lpIdx] = angle;
                        LP_Phi[lpIdx] = dihed;
                        LP_Sin_Theta[lpIdx] = sin(LP_Theta[i]*radianInv);
                        LP_Cos_Theta[lpIdx] = cos(LP_Theta[i]*radianInv);
                        LP_Sin_Phi[lpIdx] = sin(LP_Phi[i]*radianInv);
                        LP_Cos_Phi[lpIdx] = cos(LP_Phi[i]*radianInv);
                        lpIdx++;
                }
                else if(lpHost[i] == 2)   {       //LP Colinear 
                        n_LPH_Host[lphIdx] = lpHost[i];
                        Pointer_LPH[lphIdx] = pointerLP[i];
                        LPH_Dist[lphIdx] = dist;
                        LPH_Scale[lphIdx] = angle;
                        LPH_Phi[lphIdx] = dihed;
                        lphIdx++;
                }
                else {
                        fprintf(fFile_Run_Log,"%s\nThe number of host atoms of the LPH is not 2 or 3!\nQuit\n", szLine);
                        fflush(fFile_Run_Log);
                        exit(1);
                }

        }
        //end   to read the host atoms information of lone pairs

        //printf("finding lines to read \n");

        //start to read LP and LPHs indices and host indices
        for(int i=0; i<totalLP;i++)    {
                if (lpHost[i]==3) {
                        numLinesToRead+=4;
                }
                else if(lpHost[i]==2)     {
                        numLinesToRead+=3;
                }
        }
        numLinesToRead = ceil(numLinesToRead/8);
        //printf("reading lp and host indices\n");

        for(int j=0; j<(int)numLinesToRead; j++)        {
                fgets(szLine, MAX_LEN_LINE, fIn);
                token = strtok(szLine," ");
                while(token !=NULL)     {
                        tmpIndices[idxCount]=atoi(token);
                        token = strtok(NULL, " ");
                        idxCount++;
                }
        }


        //for(int i=0; i<idxCount; i++) {
        //      printf("tmpidx %d=%d\n",i,tmpIndices[i]);
        //}

        //printf("error here1\n");

        j =0; // temporary index counter
        int currentHostIdx=0;
        lphIdx=0;
        lpIdx=0;
        int hostIdx; // 1st or 2nd or 3rd host for LP or LPH
        for(int i=0; i<totalLP-1;i++)  {
                hostIdx=0; // reset for each lone pair
                nextLPIdx = pointerLP[i+1]-1;
                //printf("currenthostidx =%d, nextLPIdx =%d\n",currentHostIdx,nextLPIdx);
                if (nextLPIdx-currentHostIdx==4) {
                        LP_Host[lpIdx] = tmpIndices[j]-1; // indices start from 0 not 1
                        j++;
                        //printf("lpidx %d\n",LP_Host[lpIdx]);
                        while(j < nextLPIdx)    {
                                LP_PosHost[lpIdx][hostIdx] = tmpIndices[j]-1; // indices start from 0 not 1
                                j++;
                                //printf("hostlp %d ",LP_PosHost[lpIdx][hostIdx]);
                                hostIdx++;

                        }
                        lpIdx++;
                }
                else if (nextLPIdx-currentHostIdx==3) {
                        LPH_Host[lphIdx] = tmpIndices[j]-1;
                        j++;
                        //printf("lphidx %d\n",LPHList[lphIdx]);
                        while(j < nextLPIdx)    {
                                LPH_PosHost[lphIdx][hostIdx] = tmpIndices[j]-1;
                                j++;
                                //printf("hostlph %d ",LPH_PosHost[lphIdx][hostIdx]);
                                hostIdx++;
                        }
                        lphIdx++;
                        printf("\n");
                }
                else if(nLPH==0){
                        printf(" ");
                }
                else {
                        fprintf(fFile_Run_Log,"%s\nThe index of Lone pair and is not either 3 or 4! Please check PSF.\nQuit\n");
                        fflush(fFile_Run_Log);
                        exit(1);

                }
                currentHostIdx = nextLPIdx;
        } // end of for

        //printf("reading last remianing host\n");

        //Gets the last remaining lonepair+Host Record  
        hostIdx=0;
        if (idxCount-j==4)      {
                LP_Host[lpIdx] = tmpIndices[j]-1; // indices start from 0 not 1
                LP_PosHost[lpIdx][hostIdx] = tmpIndices[j+1]-1; // indices start from 0 not 1
                LP_PosHost[lpIdx][hostIdx+1] = tmpIndices[j+2]-1; // indices start from 0 not 1
                LP_PosHost[lpIdx][hostIdx+2] = tmpIndices[j+3]-1; // indices start from 0 not 1
        }
        else if(idxCount-j==3)  {
                LPH_Host[lphIdx] = tmpIndices[j]-1;
                LPH_PosHost[lphIdx][hostIdx] = tmpIndices[j+1]-1;
                LPH_PosHost[lphIdx][hostIdx+1] = tmpIndices[j+2]-1;
        }
        else if(nLPH==0){
                printf(" ");
        }
        else {
                fprintf(fFile_Run_Log,"%s\nThe index of Lone pair and is not either 3 or 4! Please check PSF.\nQuit\n");
                fflush(fFile_Run_Log);
                exit(1);

        }


/*
        //TESTING

        printf("Testing Lone Pairs\n");
        for(i=0; i<nLP; i++)    {
                printf("lone pair %d = %d\n",i,LP_Host[i]);
                printf("Hosts are:\n");
                for(j=0; j<3; j++)      {
                        printf("%d ",LP_PosHost[i][j]);
                }
                printf("\n");   
        }
        
        printf("Testing Lone Pairs Halos\n");
        for(i=0; i<nLPH; i++)   {
                printf("lph %d =%d\n",i,LPHList[i]);
                printf("Hosts are:\n");
                for(j=0; j<2; j++)  {
                        printf("%d ",LPH_PosHost[i][j]);                
                }
                printf("\n");
        }

*/

        //end to read LP and LPHs indices and host indices



        //end   to read the list of lone pairs and their hosts


        ///////////////////////////////////////////////////////////////////////////////// Check !!!!!!
        while(1)        {       //to find the entry for anisotropy 
                if(feof(fIn))   {
                        fprintf(fFile_Run_Log, "Fail to find the entry for anisotropy parameters in xpsf file. Please use correct version of CHARMM, e.g., c36a5\nQuit\n");
                        fflush(fFile_Run_Log);
                        exit(1);
                }
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(strncmp(szLine+11, "!NUMANISO", 9) == 0)     {       //anisotropy entry
                        break;
                }
        }
        sscanf(szLine, "%d", &nAniso);
        for(i=0; i<nAniso; i++) {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadItem = sscanf(szLine, "%lf %lf %lf", &(Para_Aniso[i][0]) , &(Para_Aniso[i][1]), &(Para_Aniso[i][2]));
                if(ReadItem !=3)        {
//                      printf("Error in reading parameters for anisotropy.\n%s\nQuit.\n", szLine);
                        fprintf(fFile_Run_Log, "Error in reading parameters for anisotropy.\n%s\nQuit.\n", szLine);
                        fflush(fFile_Run_Log);
                        exit(1);
                }

//              Para_Aniso[i][0] *= 2.0;
//              Para_Aniso[i][1] *= 2.0;
//              Para_Aniso[i][2] *= 2.0;
        }
        n_Item = 4*nAniso;
        n_Line = n_Item/ITEM_IN_LINE_8;
        n_Left = n_Item - n_Line*ITEM_IN_LINE_8;
        Count=0;
        for(i=0; i<n_Line; i++, Count+=2)       {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, AnisoList[Count], ITEM_IN_LINE_8);  //automatically read two entries
        }
        if(n_Left > 0)  {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, AnisoList[Count], 4);       //read one entry
        }
        for(i=0; i<nAniso; i++) {
                AnisoList[i][0]--;
                AnisoList[i][1]--;
                AnisoList[i][2]--;
                AnisoList[i][3]--;
        }


        nCMapTerm = 0;
        memset(CMapList, 0, sizeof(int)*MAX_ATOM*8);
        while(1)        {       //to find the entry for CMap 
                fgets(szLine, MAX_LEN_LINE, fIn);
                if(strncmp(szLine+11, "!NCRTERM: cross-terms", 21) == 0)        {       //lone pair list
                        sscanf(szLine, "%d", &nCMapTerm);
                        break;
                }
        }
        for(i=0; i<nCMapTerm; i++)      {
                fgets(szLine, MAX_LEN_LINE, fIn);
                ReadIntDataFromBuff(szLine, CMapList[i], 8);
                for(j=0; j<8; j++)      {
                        CMapList[i][j]--;       //the atom index starting from 0 !
                }
        }

        fclose(fIn);

        //start to make all index of atom starting from 0 instead of 1
        n_Max = nBond * 2;
        for(i=0; i<n_Max; i++)  {
                BondList[i]--;
        }
        n_Max = nAngle * 3;
        for(i=0; i<n_Max; i++)  {
                AngleList[i]--;
        }
        n_Max = nDihedral * 4;
        for(i=0; i<n_Max; i++)  {
                DihedralList[i]--;
        }
        n_Max = nImpro * 4;
        for(i=0; i<n_Max; i++)  {
                ImprDihedralList[i]--;
        }
        //end   to make all index of atom starting from 0 instead of 1

        nRealAtom = 0;
        for(i=0; i<nAtom; i++)  {
                if(mass[i] > 0.8)       {       // a real atom
                        RealAtom_List[nRealAtom] = i;
                        nRealAtom++;
                }
        }


        BuildDistanceMatrix();

        Geo_Opt_Drude_Only = 0;
        Is_Phi_Psi_Constrained = 0;
        E_CMap_On = 1;
        nBond_Fixed = nAngle_Fixed = nPhi_Fixed = 0;
        BuildActiveEnergyTermList();

        return;
}


void CMol::AssignForceFieldParameters(CForceField* pForceField)
{
	int i, j, iPos, Atom_1, Atom_2, Atom_3, Atom_4;
	char szChemName[8][N_LEN_CHEM_NAME];
	double Para_List[MAX_DIH_ITEM*3];

	myForceField = pForceField;

	//start	to set parameters zero
	memset(Para_k_b, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);
	memset(Para_b0, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);

	memset(Para_k_a, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);
	memset(Para_theta0, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);

	memset(Para_k_Dih, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM*MAX_DIH_ITEM);
	memset(Para_phi, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM*MAX_DIH_ITEM);

	memset(Para_k_ImpDih, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);
	memset(Para_Imp_phi, 0, sizeof(double)*MAX_FACTOR*MAX_ATOM);

	memset(Para_LJ_Epsilon, 0, sizeof(double)*MAX_ATOM);
	memset(Para_LJ_Sigma, 0, sizeof(double)*MAX_ATOM);
	//end	to set parameters zero


	//start	to assign parameters for bond
	iPos = 0;
	for(i=0; i<nBond; i++)	{
		Atom_1 = BondList[iPos  ];
		Atom_2 = BondList[iPos+1];

		if(IsDrude[Atom_2] == 1)	{	//parameters for drude
			Para_k_b[i] = K_DRUDE;	//the spring between host and drudes
			Para_b0[i] = 0.0;
			
			iPos+=2;
			continue;
		}

		strcpy(szChemName[0], ChemName[Atom_1]);
		strcpy(szChemName[1], ChemName[Atom_2]);

		myForceField->GetPara_Bond(szChemName, Para_List);
		Para_k_b[i] = Para_List[0];
		Para_b0[i] = Para_List[1];

		iPos+=2;
	}
	//end	to assign parameters for bond

	//start	to assign parameters for angle
	iPos = 0;
	for(i=0; i<nAngle; i++)	{
		Atom_1 = AngleList[iPos  ];
		Atom_2 = AngleList[iPos+1];
		Atom_3 = AngleList[iPos+2];

		strcpy(szChemName[0], ChemName[Atom_1]);
		strcpy(szChemName[1], ChemName[Atom_2]);
		strcpy(szChemName[2], ChemName[Atom_3]);

		myForceField->GetPara_Angle(szChemName, Para_List);
		Para_k_a[i] = Para_List[0];
		Para_theta0[i] = Para_List[1] * radianInv;

		if(fabs(Para_List[2]) > 1.0E-100)	{	//not zero
			Para_k_Urey[i] = Para_List[2];	//K_UB
			Para_b0_Urey[i] = Para_List[3];	//S_0
		}

		iPos+=3;	//angle
	}
	//end	to assign parameters for angle


	//start	to assign parameters for dihedral
	iPos = 0;
	for(i=0; i<nDihedral; i++)	{
		Atom_1 = DihedralList[iPos  ];
		Atom_2 = DihedralList[iPos+1];
		Atom_3 = DihedralList[iPos+2];
		Atom_4 = DihedralList[iPos+3];

		strcpy(szChemName[0], ChemName[Atom_1]);
		strcpy(szChemName[1], ChemName[Atom_2]);
		strcpy(szChemName[2], ChemName[Atom_3]);
		strcpy(szChemName[3], ChemName[Atom_4]);

		memset(Para_List, 0, sizeof(double)*MAX_DIH_ITEM*3);
		myForceField->GetPara_Dihedral(szChemName, Para_List);

		for(j=1; j<MAX_DIH_ITEM; j++)	{
			if(fabs(Para_List[j*3]) > 1.0E-50)	{	//valid parameter
				Para_k_Dih[i][j] = Para_List[j*3];
				Para_phi[i][j] = Para_List[j*3+2] * radianInv;
			}
		}

		iPos+=4;	//dihedral
	}
	//end	to assign parameters for dihedral


	//start	to assign parameters for improper dihedral
	iPos = 0;
	for(i=0; i<nImpro; i++)	{
		Atom_1 = ImprDihedralList[iPos  ];
		Atom_2 = ImprDihedralList[iPos+1];
		Atom_3 = ImprDihedralList[iPos+2];
		Atom_4 = ImprDihedralList[iPos+3];

		strcpy(szChemName[0], ChemName[Atom_1]);
		strcpy(szChemName[1], ChemName[Atom_2]);
		strcpy(szChemName[2], ChemName[Atom_3]);
		strcpy(szChemName[3], ChemName[Atom_4]);

		myForceField->GetPara_ImproDIhedral(szChemName, Para_List);
		Para_k_ImpDih[i] = Para_List[0];
		Para_Imp_phi[i] = Para_List[1] * radianInv;
		Para_Type_ImpDih[i] = Para_List[2];

		iPos+=4;	//dihedral
	}
	//end	to assign parameters for improper dihedral

	//start	to assign parameters for LJ
	for(i=0; i<nAtom; i++)	{
		strcpy(szChemName[0], ChemName[i]);

		LJ_Para_Rec[i] = myForceField->GetPara_LJ(szChemName, Para_List);	// if LJ_Para_Rec[i] == 0, be not fitted !

//		if(LJ_Para_Rec[i] >= 0)	{
//		if(fabs(myForceField->LJ_Rec[LJ_Para_Rec[i]].para[1]) > 1.0E-5)	{	// non-zero
//			fprintf(fFile_Run_Log, "Found LJ type %8s for atom %8s\n", ChemName[i], AtomName[i]);
//		}
	}
	fflush(fFile_Run_Log);
	//end	to assign parameters for LJ


	Setup_NBFix();
	Setup_NonBondParameters();

	Identify_Bonds_With_H();
}

void CMol::ReadPDB(char szNamePdb[])
{
	FILE *fIn;
	int LineRead, ReadItem;
	char szLine[256];

	fIn = fopen(szNamePdb, "r");
	if(fIn == NULL)	{
//		printf("Fail to open file %s\nQuit\n", szNamePdb);
		fprintf(fFile_Run_Log, "Fail to open file %s\nQuit\n", szNamePdb);
		fflush(fFile_Run_Log);
		exit(1);
	}

	LineRead = 0;
	while(1)	{
		fgets(szLine, MAX_LEN_LINE, fIn);

		if(feof(fIn))	{
			if(LineRead < nAtom)	{
//				printf("From CMol::ReadPDB() > LineRead < nAtom\nQuit\n");
				fprintf(fFile_Run_Log, "From CMol::ReadPDB() > LineRead < nAtom\nQuit\n");
				fflush(fFile_Run_Log);
				fclose(fIn);
				exit(1);
			}
			break;
		}
		if(strncmp(szLine, "ATOM", 4) == 0)	{
			ReadItem = sscanf(szLine+30, "%lf %lf %lf", &(x[LineRead]), &(y[LineRead]), &(z[LineRead]));
			if(ReadItem == 3)	{
				LineRead++;
			}
			else	{
//				printf("Error in reading file: %s\n%s\nQuit\n", szNamePdb, szLine);
				fprintf(fFile_Run_Log, "Error in reading file: %s\n%s\nQuit\n", szNamePdb, szLine);
				fflush(fFile_Run_Log);
				fclose(fIn);
				exit(1);
			}
		}
	}
	fclose(fIn);
	Position_LonePair();
}

void CMol::ReadCRD(char szNameCrd[])
{
	FILE *fIn;
	int i, ReadItem;
	char szLine[256], szTmp[256];

	fIn = fopen(szNameCrd, "r");
	if(fIn == NULL)	{
//		printf("Fail to open file %s\nQuit\n", szNameCrd);
		fprintf(fFile_Run_Log, "Fail to open file %s\nQuit\n", szNameCrd);
		fflush(fFile_Run_Log);
		exit(1);
	}
	while(1)	{
		fgets(szLine, MAX_LEN_LINE, fIn);

		if(strncmp(szLine+10, "  EXT", 5) == 0)	{	//find the flag for coordinates
			break;
		}
	}

	for(i=0; i<nAtom; i++)	{
		fgets(szLine, MAX_LEN_LINE, fIn);
//		ReadItem = sscanf(szLine+40, "%lf %lf %lf", &(x[i]), &(y[i]), &(z[i]));
		ReadItem = sscanf(szLine, "%s%s%s%s%lf %lf %lf %s%s%s", 
			szTmp, szTmp, szTmp, szTmp, &(x[i]), &(y[i]), &(z[i]), szTmp, szTmp, szTmp);
			if(ReadItem != 10)	{
//				printf("Error in reading file: %s\n%s\nQuit\n", szNameCrd, szLine);
				fprintf(fFile_Run_Log, "Error in reading file: %s\n%s\nQuit\n", szNameCrd, szLine);
				fflush(fFile_Run_Log);
				fclose(fIn);
				exit(1);
			}
	}
	fclose(fIn);
	Position_LonePair();
}

//#define DeltaX	(1.0E-8)
#define DeltaX	(5.0E-8)
//#define DeltaX	(1.0E-6)

void CMol::TestFirstDerivative(void)
{
	int i, AtmIdx;
	double gx_Analytical[MAX_ATOM], gy_Analytical[MAX_ATOM], gz_Analytical[MAX_ATOM];
	double gx_Numerical[MAX_ATOM], gy_Numerical[MAX_ATOM], gz_Numerical[MAX_ATOM];
	double x_Save, y_Right, y_Left;
	
	Cal_E(0);
	
	for(i=0; i<nAtom; i++)	{	//calculate the analytical first derivative
		gx_Analytical[i] = grad_x[i];
		gy_Analytical[i] = grad_y[i];
		gz_Analytical[i] = grad_z[i];
	}
	
	for(AtmIdx=0; AtmIdx<nAtom; AtmIdx++)	{
		x_Save = x[AtmIdx];
		
		x[AtmIdx] = x_Save + DeltaX;
		Cal_E(0);
		y_Right = E_Total;
		
		x[AtmIdx] = x_Save - DeltaX;
		Cal_E(0);
		y_Left = E_Total;
		
		gx_Numerical[AtmIdx] = (y_Right-y_Left)/(2.0*DeltaX);
		
		
		x_Save = y[AtmIdx];
		
		y[AtmIdx] = x_Save + DeltaX;
		Cal_E(0);
		y_Right = E_Total;
		
		y[AtmIdx] = x_Save - DeltaX;
		Cal_E(0);
		y_Left = E_Total;
		
		gy_Numerical[AtmIdx] = (y_Right-y_Left)/(2.0*DeltaX);
		
		
		x_Save = z[AtmIdx];
		
		z[AtmIdx] = x_Save + DeltaX;
		Cal_E(0);
		y_Right = E_Total;
		
		z[AtmIdx] = x_Save - DeltaX;
		Cal_E(0);
		y_Left = E_Total;
		
		gz_Numerical[AtmIdx] = (y_Right-y_Left)/(2.0*DeltaX);
		
//		printf("%d   %18.14lf  %18.14lf %18.14lf\n", 
		fprintf(fFile_Run_Log, "%d   %18.14lf  %18.14lf %18.14lf\n", 
			AtmIdx, 
			fabs(gx_Numerical[AtmIdx]-gx_Analytical[AtmIdx]), 
			fabs(gx_Numerical[AtmIdx]-gx_Analytical[AtmIdx]), 
			fabs(gx_Numerical[AtmIdx]-gx_Analytical[AtmIdx]));
		fflush(fFile_Run_Log);
	}
}

#undef DeltaX


#define	MUpdate		(5)
#define DIAG_LBFGS	(1.0)
#define	TINY		(1.0E-200)

#define	RMSTol		(2.0E-8)

double CMol::FullGeometryOptimization_LBFGS(int Only_Opt_Drude)
{
	//start	data for LBFGS
	double YS, YY, SQ, YR, BETA, GNORM, GNORM_Min=1.0E100, x_Best[MAX_ATOM], y_Best[MAX_ATOM], z_Best[MAX_ATOM];
	double Gtmp[3*MAX_ATOM], Diag[3*MAX_ATOM], Stp[3*MAX_ATOM], Step, Dot_G_Gtmp;
	double Rho1[MUpdate], Alpha[MUpdate];
	double SearchStep[MUpdate][3*MAX_ATOM], GDif[MUpdate][3*MAX_ATOM], G_1D[3*MAX_ATOM];
	int Point, PointLast=0, Bound, CP, NIterDone, i, j, Idx, nDim;
	int NIterMax;
	double MaxStep;
	//end	data for LBFGS

	if(Only_Opt_Drude)	{
		NIterMax = 200;
		MaxStep = 0.006;
	}
	else	{
		NIterMax = 8000;
		MaxStep = 0.25;
	}

	nDim = 3*nAtom;	//the real dimensions of variables
	Cal_E(0);

	//start	to pack gradient into 1D array
	Idx = 0;
	for(i=0; i<nAtom; i++, Idx+=3)	{
		G_1D[Idx  ] = grad_x[i];
		G_1D[Idx+1] = grad_y[i];
		G_1D[Idx+2] = grad_z[i];

		x_Best[i] = x[i];	//init the best coordinates
		y_Best[i] = y[i];
		z_Best[i] = z[i];
	}
	//end	to pack gradient into 1D array

	for(NIterDone=1; NIterDone<=NIterMax; NIterDone++)	{
		if(NIterDone==1)	{
			Point=0;

			GNORM=0.0;
			for(i=0; i<nDim; i++)	{
				Diag[i]=DIAG_LBFGS;
				SearchStep[0][i]=-G_1D[i]*Diag[i];
				Gtmp[i]=SearchStep[0][i];
				GNORM+=(G_1D[i]*G_1D[i]);
			}
			GNORM=sqrt(GNORM);
//			printf("GNORM = %lf\n", GNORM);
			Step=min(1.0/GNORM, GNORM);
			for(i=0; i<nDim; i++)	{
				Stp[i]=Step;
			}
		}
		else	{
			Bound=NIterDone-1;
			if(NIterDone > MUpdate)	{
				Bound=MUpdate;
			}
			YS=YY=0.0;
			for(i=0; i<nDim; i++)	{
				YS+=(GDif[PointLast][i]*SearchStep[PointLast][i]);
				YY+=(GDif[PointLast][i]*GDif[PointLast][i]);
			}
			if(YS < TINY)	{
				YS=1.0;
			}
			if(YY < TINY)	{
				YY=1.0;
			}
			for(i=0; i<nDim; i++)	{
				Diag[i]=YS/YY;
			}

			Rho1[PointLast]=1.0/YS;
			for(i=0; i<nDim; i++)	{
				Gtmp[i]=-G_1D[i];
			}
			CP=Point;
			for(j=1; j<=Bound; j++)	{
				CP--;
				if(CP <0)	{
					CP=MUpdate-1;
				}
				SQ=0.0;
				for(i=0; i<nDim; i++)	{
					SQ+=(SearchStep[CP][i]*Gtmp[i]);
				}
				Alpha[CP]=Rho1[CP]*SQ;
				for(i=0; i<nDim; i++)	{
					Gtmp[i]-=(Alpha[CP]*GDif[CP][i]);
				}
			}

			for(i=0; i<nDim; i++)	{
				Gtmp[i]*=Diag[i];
			}

			for(j=1; j<=Bound; j++)	{
				YR=0.0;
				for(i=0; i<nDim; i++)	{
					YR+=(GDif[CP][i]*Gtmp[i]);
				}
				BETA=Rho1[CP]*YR;
				BETA=Alpha[CP]-BETA;
				for(i=0; i<nDim; i++)	{
					Gtmp[i]+=(BETA*SearchStep[CP][i]);
				}
				CP++;
				if(CP==MUpdate)	{
					CP=0;
				}
			}

			for(i=0; i<nDim; i++)	{
				Stp[i]=1.0;
			}
		}

		if(NIterDone > 1)	{
			for(i=0; i<nDim; i++)	{
				SearchStep[Point][i]=Gtmp[i];
			}
		}


		Dot_G_Gtmp=0.0;
		for(i=0; i<nDim; i++)	{
			Dot_G_Gtmp+=(G_1D[i]*Gtmp[i]);
		}
		if(Dot_G_Gtmp > 0.0)	{	//Reversing step
			for(i=0; i<nDim; i++)	{
				Gtmp[i]=-Gtmp[i];
				SearchStep[Point][i]=Gtmp[i];
			}
		}
		for(i=0; i<nDim; i++)	{
			Gtmp[i]=G_1D[i];
		}

		//start	to check the step size to make sure it is smaller than MaxStep
		Step=0.0;
		for(i=0; i<nDim; i++)	{
			Step+=(SearchStep[Point][i]*SearchStep[Point][i]);
		}
		Step=sqrt(Step);
		if(Step > MaxStep)	{
			for(i=0; i<nDim; i++)	{
				Stp[i]=MaxStep/Step;
			}
		}
		//end	to check the step size to make sure it is smaller than MaxStep

		//start	to update the coordinates
		Idx = 0;
		for(i=0; i<nAtom; i++)	{
			SearchStep[Point][Idx]*=Stp[Idx];
			x[i]+=(SearchStep[Point][Idx]);
			Idx++;

			SearchStep[Point][Idx]*=Stp[Idx];
			y[i]+=(SearchStep[Point][Idx]);
			Idx++;

			SearchStep[Point][Idx]*=Stp[Idx];
			z[i]+=(SearchStep[Point][Idx]);
			Idx++;
		}
		//end	to update the coordinates

		
		Cal_E(0);
		
		
		//start	to pack gradient into 1D array
		Idx = 0;
		for(i=0; i<nAtom; i++, Idx+=3)	{
			G_1D[Idx  ] = grad_x[i];
			G_1D[Idx+1] = grad_y[i];
			G_1D[Idx+2] = grad_z[i];
		}
		//end	to pack gradient into 1D array
		

		for(i=0; i<nDim; i++)	{
			GDif[Point][i]=G_1D[i]-Gtmp[i];
		}
		PointLast=Point;
		Point=(++Point)%MUpdate;	//Increase first, then use the value

		GNORM=0.0;
		for(i=0; i<nDim; i++)	{
			GNORM+=(G_1D[i]*G_1D[i]);
		}
		GNORM=sqrt(GNORM);
//		printf("GNORM = %lf\n", GNORM);

		if(GNORM < GNORM_Min)	{	//better coordinates found
			for(i=0; i<nAtom; i++)	{
				x_Best[i] = x[i];	//save the best coordinates
				y_Best[i] = y[i];
				z_Best[i] = z[i];
			}
			GNORM_Min = GNORM;
		}

		if(GNORM < RMSTol)	{
//			printf("Successfully converge in %d steps.\nQuit.\n", NIterDone);
			break;
		}
	}
	
	for(i=0; i<nAtom; i++)	{
		x[i] = x_Best[i];	//make sure the best coordinates are used
		y[i] = y_Best[i];
		z[i] = z_Best[i];
	}

//	printf("NIterDone = %d\n", NIterDone);
	if(NIterDone > NIterMax)	{
		fprintf(fFile_Run_Log, "NIterDone > NIterMax in optimization. You may want to tune the parameters in L-BFGS.\n");
		fflush(fFile_Run_Log);
	}

	return (Cal_E(0));
}
//#undef NIterMax
#undef MUpdate
#undef DIAG_LBFGS
//#undef MaxStep
#undef TINY
#undef RMSTol


#define BOHRR (0.529177249)
#define DEBYEC (2.541766 / BOHRR)

double CMol::Cal_Dipole(double& Dipole_x, double& Dipole_y, double& Dipole_z, double& Dipole)
{
	int i;
	double TotalCharge;

	TotalCharge = 0.0;
	Dipole = Dipole_x = Dipole_y = Dipole_z = 0.0;

	for(i=0; i<nAtom; i++)	{
		Dipole_x += (CG[i]*x[i]);
		Dipole_y += (CG[i]*y[i]);
		Dipole_z += (CG[i]*z[i]);
		TotalCharge += CG[i];
	}
	Dipole_x *= DEBYEC;
	Dipole_y *= DEBYEC;
	Dipole_z *= DEBYEC;
	Dipole = sqrt(Dipole_x*Dipole_x + Dipole_y*Dipole_y + Dipole_z*Dipole_z);

	if(fabs(TotalCharge) > 1.0E-4)	{
//		printf("Warning! Total charge is not zero!\nCG = %lf\n", TotalCharge);
		fprintf(fFile_Run_Log, "Warning! Total charge is not zero!\nCG = %lf\n", TotalCharge);
		fflush(fFile_Run_Log);
	}
//	printf("%7.5lf %7.5lf %7.5lf %7.5lf\n", Dipole_x, Dipole_y, Dipole_z, Dipole);

	return Dipole;
}

void CMol::Cal_Force_Short(void)
{
	E_Bond = E_UreyB = E_Angle = E_Dihedral = E_ImproDihedral = E_CMap = E_Aniso = E_Thole = E_DrudeHyper = 0.0;
	E_Constrain_Phi = E_Constrain_Psi = 0.0;

	memset(grad_x, 0, sizeof(double)*nAtom);	//zero the gradient
	memset(grad_y, 0, sizeof(double)*nAtom);
	memset(grad_z, 0, sizeof(double)*nAtom);

	Cal_E_Bond();
	Cal_E_UreyB();
	Cal_E_Angle();
	Cal_E_Dihedral(0);
	Cal_E_ImproperDihedral();
//	Cal_E_VDW_ELEC();

//	Cal_E_User();

	if(E_CMap_On)	{
		Cal_E_CMap();
	}

	Cal_E_Anisotropy();
	Cal_E_Thole();
	Cal_E_DrudeHyper();

//	if(Is_Phi_Psi_Constrained == 1)	{
//		Cal_E_Constrain_Phi();
//		Cal_E_Constrain_Psi();
//	}

	for(int i=0; i<nAtom; i++)	{
		Bonded_fx[i] = -grad_x[i];
		Bonded_fy[i] = -grad_y[i];
		Bonded_fz[i] = -grad_z[i];
	}


	return;
}

void CMol::Cal_Force_Long(void)
{
	E_VDW = E_Elec = E_User = 0.0;

	memset(grad_x, 0, sizeof(double)*nAtom);	//zero the gradient
	memset(grad_y, 0, sizeof(double)*nAtom);
	memset(grad_z, 0, sizeof(double)*nAtom);

	Position_LonePair();

	Cal_E_VDW_ELEC();
	Cal_E_User();

	E_Bond += E_Aniso;	// as CHARMM does
	E_Elec += E_Thole;	// as CHARMM does

	E_Total = E_Bond + E_UreyB + E_Angle + E_Dihedral + E_ImproDihedral + E_VDW + E_Elec + E_CMap + E_DrudeHyper + E_User;
	E_Total += (E_Constrain_Phi+E_Constrain_Psi);

	LonePairForceReDistribute();

	for(int i=0; i<nAtom; i++)	{
		Nonbonded_fx[i] = -grad_x[i];
		Nonbonded_fy[i] = -grad_y[i];
		Nonbonded_fz[i] = -grad_z[i];
	}


//	printf("EBond   = %26.18lf\nE_Angle =%26.18lf\nE_UreyB = %26.18lf\nE_Dihed = %26.18lf\nE_Impro = %26.18lf\nE_VDW   = %26.18lf\nE_Elec  = %26.18lf\nE_CMap  = %26.18lf\nE_User  = %26.18lf\n\nE_Total = %26.18lf\n\n", 
//		E_Bond, E_Angle, E_UreyB, E_Dihedral, E_ImproDihedral, E_VDW, E_Elec, E_CMap, E_User, E_Total);

	return;
}

double CMol::Cal_E(int PrintE)
{
	int i;

	E_Total = E_Bond = E_UreyB = E_Angle = E_Dihedral = E_ImproDihedral = E_VDW = E_Elec = E_CMap = E_Aniso = E_Thole = E_DrudeHyper = E_User = 0.0;
	E_Constrain_Phi = E_Constrain_Psi = 0.0;

	memset(grad_x, 0, sizeof(double)*nAtom);	//zero the gradient
	memset(grad_y, 0, sizeof(double)*nAtom);
	memset(grad_z, 0, sizeof(double)*nAtom);

	Position_LonePair();

	Cal_E_Bond();
	Cal_E_UreyB();
	Cal_E_Angle();
	Cal_E_Dihedral(0);
	Cal_E_ImproperDihedral();
	Cal_E_VDW_ELEC();

//	Cal_E_User();

	if(E_CMap_On)	{
		Cal_E_CMap();
	}

	Cal_E_Anisotropy();
	Cal_E_Thole();
	Cal_E_DrudeHyper();

//	if(Is_Phi_Psi_Constrained == 1)	{
//		Cal_E_Constrain_Phi();
//		Cal_E_Constrain_Psi();
//	}


	E_Bond += E_Aniso;	// as CHARMM does
	E_Elec += E_Thole;	// as CHARMM does
	
	E_Total = E_Bond + E_UreyB + E_Angle + E_Dihedral + E_ImproDihedral + E_VDW + E_Elec + E_CMap + E_DrudeHyper + E_User;

	E_Total += (E_Constrain_Phi+E_Constrain_Psi);

	LonePairForceReDistribute();

	if(Geo_Opt_Drude_Only)	{
		for(i=0; i<nAtom; i++)	{
			if(IsFixed[i] != 0)	{	//set all constrained atoms with zero forces
				grad_x[i] = 0.0;
				grad_y[i] = 0.0;
				grad_z[i] = 0.0;
			}
		}
	}


	if(PrintE)	{
		printf("EBond   = %26.18lf\nE_Angle =%26.18lf\nE_UreyB = %26.18lf\nE_Dihed = %26.18lf\nE_Impro = %26.18lf\nE_VDW   = %26.18lf\nE_Elec  = %26.18lf\nE_CMap  = %26.18lf\nE_User  = %26.18lf\n\nE_Total = %26.18lf\n\n", 
			E_Bond, E_Angle, E_UreyB, E_Dihedral, E_ImproDihedral, E_VDW, E_Elec, E_CMap, E_User, E_Total);
	}

//	for(int i=0; i<nAtom; i++)	{
//		printf("%21.18e %21.18e %21.18e\n", grad_x[i], grad_y[i], grad_z[i]);
//		printf("%25.18lf %25.18lf %25.18lf\n", grad_x[i], grad_y[i], grad_z[i]);
//		printf("%25.18lf %25.18lf %25.18lf\n", x[i], y[i], z[i]);
//	}

//	Cal_Dipole();

	return E_Total;
}

extern int Restr_Atom_1, Restr_Atom_2;

void CMol::Cal_E_User(void)
{
	int ia, ib;
	double Rij, dx, dy, dz, dr, dedr, dedrx, dedry, dedrz;
	double r0_12=3.0, r0_34=7.0, kspring=10.0;

	E_User = 0.0;

//	return;

	ia = 1;
	ib = 3;
	
	dx=x[ia]-x[ib];
	dy=y[ia]-y[ib];
	dz=z[ia]-z[ib];
	Rij=sqrt(dx*dx+dy*dy+dz*dz);
	
	d_Restrain = Rij;
	
	dr=Rij-r0_12;

	E_User+=(kspring*dr*dr);
	
	dedr=(2.0*kspring*dr)/Rij;
	dedrx=dedr*dx;
	dedry=dedr*dy;
	dedrz=dedr*dz;
	
	grad_x[ia] += dedrx;	//gradient
	grad_y[ia] += dedry;
	grad_z[ia] += dedrz;
	grad_x[ib] -= dedrx;
	grad_y[ib] -= dedry;
	grad_z[ib] -= dedrz;


	ia = 0;
	ib = 4;
	
	dx=x[ia]-x[ib];
	dy=y[ia]-y[ib];
	dz=z[ia]-z[ib];
	Rij=sqrt(dx*dx+dy*dy+dz*dz);
	
	d_Restrain = Rij;
	
	dr=Rij-r0_34;

	E_User+=(kspring*dr*dr);
	
	dedr=(2.0*kspring*dr)/Rij;
	dedrx=dedr*dx;
	dedry=dedr*dy;
	dedrz=dedr*dz;
	
	grad_x[ia] += dedrx;	//gradient
	grad_y[ia] += dedry;
	grad_z[ia] += dedrz;
	grad_x[ib] -= dedrx;
	grad_y[ib] -= dedry;
	grad_z[ib] -= dedrz;

}

void CMol::Cal_E_Bond(void)
{
	int iActive, i, iPos, ia, ib;
	double Rij, dx, dy, dz, dr, dedr, dedrx, dedry, dedrz;

	E_Bond = 0.0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_Bond; iActive++)	{
			i = Active_List_Bond[iActive];	//the index of bond
			iPos = 2*i;
			
			ia = BondList[iPos  ];
			ib = BondList[iPos+1];
			
			dx=x[ia]-x[ib];
			dy=y[ia]-y[ib];
			dz=z[ia]-z[ib];
			Rij=sqrt(dx*dx+dy*dy+dz*dz);
			
			
			if(Rij < 1.0E-20)	{	//zero
				continue;
			}
			
			dr=Rij-Para_b0[i];
			
			E_Bond+=(Para_k_b[i]*dr*dr);
			
			dedr=(2.0*Para_k_b[i]*dr)/Rij;
			dedrx=dedr*dx;
			dedry=dedr*dy;
			dedrz=dedr*dz;
			
			grad_x[ia] += dedrx;	//gradient
			grad_y[ia] += dedry;
			grad_z[ia] += dedrz;
			grad_x[ib] -= dedrx;
			grad_y[ib] -= dedry;
			grad_z[ib] -= dedrz;
		}
	}
	else	{
		for(iActive=0; iActive<nBond; iActive++)	{
			i = iActive;	//the index of bond
			iPos = 2*i;
			
			ia = BondList[iPos  ];
			ib = BondList[iPos+1];
			
			dx=x[ia]-x[ib];
			dy=y[ia]-y[ib];
			dz=z[ia]-z[ib];
			Rij=sqrt(dx*dx+dy*dy+dz*dz);
			
			
			if(Rij < 1.0E-20)	{	//zero
				continue;
			}
			
			dr=Rij-Para_b0[i];
			
			E_Bond+=(Para_k_b[i]*dr*dr);
			
			dedr=(2.0*Para_k_b[i]*dr)/Rij;
			dedrx=dedr*dx;
			dedry=dedr*dy;
			dedrz=dedr*dz;
			
			grad_x[ia] += dedrx;	//gradient
			grad_y[ia] += dedry;
			grad_z[ia] += dedrz;
			grad_x[ib] -= dedrx;
			grad_y[ib] -= dedry;
			grad_z[ib] -= dedrz;
		}
	}




	return;
}

void CMol::Cal_E_UreyB(void)
{
	int iActive, i, iPos, ia, ib;
	double Rij, dx, dy, dz, dr, dedr, dedrx, dedry, dedrz;

	E_UreyB = 0.0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_Angle; iActive++)	{
			i = Active_List_Angle[iActive];	//the index of angle
			iPos = 3*i;
			
			ia = AngleList[iPos  ];
			ib = AngleList[iPos+2];
			
			dx=x[ia]-x[ib];
			dy=y[ia]-y[ib];
			dz=z[ia]-z[ib];
			Rij=sqrt(dx*dx+dy*dy+dz*dz);
			
			
			if(Rij < 1.0E-20)	{	//zero
				continue;
			}
			
			dr=Rij-Para_b0_Urey[i];
			
			E_UreyB+=(Para_k_Urey[i]*dr*dr);
			
			dedr=(2.0*Para_k_Urey[i]*dr)/Rij;
			dedrx=dedr*dx;
			dedry=dedr*dy;
			dedrz=dedr*dz;
			
			grad_x[ia] += dedrx;	//gradient
			grad_y[ia] += dedry;
			grad_z[ia] += dedrz;
			grad_x[ib] -= dedrx;
			grad_y[ib] -= dedry;
			grad_z[ib] -= dedrz;
		}
	}
	else	{
		for(iActive=0; iActive<nAngle; iActive++)	{
			i = iActive;	//the index of angle
			iPos = 3*i;
			
			ia = AngleList[iPos  ];
			ib = AngleList[iPos+2];
			
			dx=x[ia]-x[ib];
			dy=y[ia]-y[ib];
			dz=z[ia]-z[ib];
			Rij=sqrt(dx*dx+dy*dy+dz*dz);
			
			
			if(Rij < 1.0E-20)	{	//zero
				continue;
			}
			
			dr=Rij-Para_b0_Urey[i];
			
			E_UreyB+=(Para_k_Urey[i]*dr*dr);
			
			dedr=(2.0*Para_k_Urey[i]*dr)/Rij;
			dedrx=dedr*dx;
			dedry=dedr*dy;
			dedrz=dedr*dz;
			
			grad_x[ia] += dedrx;	//gradient
			grad_y[ia] += dedry;
			grad_z[ia] += dedrz;
			grad_x[ib] -= dedrx;
			grad_y[ib] -= dedry;
			grad_z[ib] -= dedrz;
		}
	}


	return;
}

void CMol::Cal_E_Angle(void)
{
	int iActive, i, iPos, ia, ib, ic;
	double xia, yia, zia, xib, yib, zib, xic, yic, zic;
	double xab, yab, zab, xcb, ycb, zcb, rab2, rcb2;
	double xp, yp, zp, rp, dot, cosine, angle;
	double dt, deddt;
	double terma, termc, dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic;
	
	E_Angle=0.0;

	if(Geo_Opt_Drude_Only)	{	// optimize drude positions only
		for(iActive=0; iActive<nActive_Angle; iActive++)	{
			i = Active_List_Angle[iActive];	//the index of angle
			iPos = 3*i;
			
			ia = AngleList[iPos  ];
			ib = AngleList[iPos+1];
			ic = AngleList[iPos+2];
			
			xia = x[ia];	yia = y[ia];	zia = z[ia];
			xib = x[ib];	yib = y[ib];	zib = z[ib];
			xic = x[ic];	yic = y[ic];	zic = z[ic];
			
			xab = xia - xib;	yab = yia - yib;	zab = zia - zib;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			
			rab2 = xab*xab + yab*yab + zab*zab;
			rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
			
			rab2=max(rab2, 0.000000001);
			rcb2=max(rcb2, 0.000000001);
			
			xp = ycb*zab - zcb*yab;
			yp = zcb*xab - xcb*zab;
			zp = xcb*yab - ycb*xab;
			rp = sqrt(xp*xp + yp*yp + zp*zp);
			rp = max(rp,0.000000001);
			dot = xab*xcb + yab*ycb + zab*zcb;
			cosine = dot / sqrt(rab2*rcb2);
			cosine = min(1.0,max(-1.0,cosine));
			angle = acos(cosine);
			dt = angle - Para_theta0[i];
			
			E_Angle+=(Para_k_a[i]*dt*dt);
			
			deddt=2.0*Para_k_a[i]*dt;
			
			terma = -deddt / (rab2*rp);
			termc = deddt / (rcb2*rp);
			dedxia = terma * (yab*zp-zab*yp);
			dedyia = terma * (zab*xp-xab*zp);
			dedzia = terma * (xab*yp-yab*xp);
			dedxic = termc * (ycb*zp-zcb*yp);
			dedyic = termc * (zcb*xp-xcb*zp);
			dedzic = termc * (xcb*yp-ycb*xp);
			dedxib = -dedxia - dedxic;
			dedyib = -dedyia - dedyic;
			dedzib = -dedzia - dedzic;
			
			grad_x[ia] += dedxia;
			grad_y[ia] += dedyia;
			grad_z[ia] += dedzia;
			
			grad_x[ib] += dedxib;
			grad_y[ib] += dedyib;
			grad_z[ib] += dedzib;
			
			grad_x[ic] += dedxic;
			grad_y[ic] += dedyic;
			grad_z[ic] += dedzic;
		}
	}
	else	{	// full optimizations
		for(iActive=0; iActive<nAngle; iActive++)	{
			i = iActive;
			iPos = 3*iActive;
			
			ia = AngleList[iPos  ];
			ib = AngleList[iPos+1];
			ic = AngleList[iPos+2];
			
			xia = x[ia];	yia = y[ia];	zia = z[ia];
			xib = x[ib];	yib = y[ib];	zib = z[ib];
			xic = x[ic];	yic = y[ic];	zic = z[ic];
			
			xab = xia - xib;	yab = yia - yib;	zab = zia - zib;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			
			rab2 = xab*xab + yab*yab + zab*zab;
			rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
			
			rab2=max(rab2, 0.000000001);
			rcb2=max(rcb2, 0.000000001);
			
			xp = ycb*zab - zcb*yab;
			yp = zcb*xab - xcb*zab;
			zp = xcb*yab - ycb*xab;
			rp = sqrt(xp*xp + yp*yp + zp*zp);
			rp = max(rp,0.000000001);
			dot = xab*xcb + yab*ycb + zab*zcb;
			cosine = dot / sqrt(rab2*rcb2);
			cosine = min(1.0,max(-1.0,cosine));
			angle = acos(cosine);
			dt = angle - Para_theta0[i];
			
			E_Angle+=(Para_k_a[i]*dt*dt);
			
			deddt=2.0*Para_k_a[i]*dt;
			
			terma = -deddt / (rab2*rp);
			termc = deddt / (rcb2*rp);
			dedxia = terma * (yab*zp-zab*yp);
			dedyia = terma * (zab*xp-xab*zp);
			dedzia = terma * (xab*yp-yab*xp);
			dedxic = termc * (ycb*zp-zcb*yp);
			dedyic = termc * (zcb*xp-xcb*zp);
			dedzic = termc * (xcb*yp-ycb*xp);
			dedxib = -dedxia - dedxic;
			dedyib = -dedyia - dedyic;
			dedzib = -dedzia - dedzic;
			
			grad_x[ia] += dedxia;
			grad_y[ia] += dedyia;
			grad_z[ia] += dedzia;
			
			grad_x[ib] += dedxib;
			grad_y[ib] += dedyib;
			grad_z[ib] += dedzib;
			
			grad_x[ic] += dedxic;
			grad_y[ic] += dedyic;
			grad_z[ic] += dedzic;
		}
	}

	
	return;
}

/*
void CMol::Gen_Phi_Psi_Map(void)
{
	char szName[256];
	FILE *fOut;
	int i;

	fOut = fopen("coor-cmap_cal-576.dat", "w");
	for(Phi0_Constrain=-180.0; Phi0_Constrain<180.0; Phi0_Constrain+=15.0)	{
		for(Psi0_Constrain=-180.0; Psi0_Constrain<180.0; Psi0_Constrain+=15.0)	{
			sprintf(szName, "F:/Keti-Chicago/Kinase/code/Parameterization/Auto-para/CMap/from_Benoit/coord-phi-psi/alad-dr_phi%d_psi%d_start.cor", 
				(int)floor(Phi0_Constrain+0.1), (int)floor(Psi0_Constrain+0.1));
			ReadCRD(szName);
			FullGeometryOptimization_LBFGS();
			printf("E(%5.0lf,%5.0lf) = %12.7lf   Phi = %7.2lf  Psi = %7.2lf\n", 
				Phi0_Constrain, Psi0_Constrain, E_Total-E_Constrain_Phi-E_Constrain_Psi, Phi_Real, Psi_Real);

			for(i=0; i<nAtom; i++)	{
				fprintf(fOut, "%20.15lf %20.15lf %20.15lf\n", x[i], y[i], z[i]);
			}
		}
	}

	fclose(fOut);

}
*/



void CMol::Gen_Phi_Psi_Map(void)
{
	FILE *fIn, *fOut;
	int i, Count=0;
	double E_Min=1.0E100, fTmp;
	double E_CMap[CMAP_DIM*CMAP_DIM];

	fIn = fopen("coor-cmap_cal-576.dat", "r");
	for(Phi0_Constrain=-180.0; Phi0_Constrain<180.0; Phi0_Constrain+=15.0)	{
		for(Psi0_Constrain=-180.0; Psi0_Constrain<180.0; Psi0_Constrain+=15.0)	{
			for(i=0; i<nAtom; i++)	{	//read in the optimized parameters
				fscanf(fIn, "%lf %lf %lf\n", &(x[i]), &(y[i]), &(z[i]));
			}

			FullGeometryOptimization_LBFGS(0);
//			printf("E(%5.0lf,%5.0lf) = %12.7lf   Phi = %7.2lf  Psi = %7.2lf\n", 
			fprintf(fFile_Run_Log, "E(%5.0lf,%5.0lf) = %12.7lf   Phi = %7.2lf  Psi = %7.2lf\n", 
				Phi0_Constrain, Psi0_Constrain, E_Total-E_Constrain_Phi-E_Constrain_Psi, Phi_Real, Psi_Real);
			fflush(fFile_Run_Log);

			E_Phi_Psi_MM[Count] = E_Total-E_Constrain_Phi-E_Constrain_Psi;
			if(E_Phi_Psi_MM[Count] < E_Min)	{
				E_Min = E_Phi_Psi_MM[Count];
			}
			Count++;
		}
	}
	fclose(fIn);

	fIn = fopen("E_Phi_Psi_QM-Pedro-576.dat", "r");
	if(fIn == NULL)	{
//		printf("Fail to open file: E_Phi_Psi_QM-Pedro-576.dat\nQuit.\n");
		fprintf(fFile_Run_Log, "Fail to open file: E_Phi_Psi_QM-Pedro-576.dat\nQuit.\n");
		fflush(fFile_Run_Log);
		exit(1);
	}
	fOut = fopen("my-new-cmap.dat", "w");

	Count = 0;
	for(Phi0_Constrain=-180.0; Phi0_Constrain<180.0; Phi0_Constrain+=15.0)	{	//shift the MM energies
		for(Psi0_Constrain=-180.0; Psi0_Constrain<180.0; Psi0_Constrain+=15.0)	{
			fscanf(fIn, "%lf %lf %lf", &fTmp, &fTmp, &(E_Phi_Psi_QM[Count]));
			E_Phi_Psi_MM[Count] -= E_Min;
			E_CMap[Count] = E_Phi_Psi_QM[Count]-E_Phi_Psi_MM[Count];
			fprintf(fOut, "%lf\n", E_CMap[Count]);
			Count++;
		}
	}
	fclose(fOut);
	fclose(fIn);


	fOut = fopen("new-cmap-charmm.txt", "w");
	Count = 0;
	for(Phi0_Constrain=-180.0; Phi0_Constrain<180.0; Phi0_Constrain+=15.0)	{	//shift the MM energies
		fprintf(fOut, "!%d\n", (int)Phi0_Constrain);
		fprintf(fOut, "%12.5lf%12.5lf%12.5lf%12.5lf%12.5lf\n", 
			E_CMap[Count  ], E_CMap[Count+1], E_CMap[Count+2], E_CMap[Count+3], E_CMap[Count+4]);
		Count+=5;
		fprintf(fOut, "%12.5lf%12.5lf%12.5lf%12.5lf%12.5lf\n", 
			E_CMap[Count  ], E_CMap[Count+1], E_CMap[Count+2], E_CMap[Count+3], E_CMap[Count+4]);
		Count+=5;
		fprintf(fOut, "%12.5lf%12.5lf%12.5lf%12.5lf%12.5lf\n", 
			E_CMap[Count  ], E_CMap[Count+1], E_CMap[Count+2], E_CMap[Count+3], E_CMap[Count+4]);
		Count+=5;
		fprintf(fOut, "%12.5lf%12.5lf%12.5lf%12.5lf%12.5lf\n", 
			E_CMap[Count  ], E_CMap[Count+1], E_CMap[Count+2], E_CMap[Count+3], E_CMap[Count+4]);
		Count+=5;
		fprintf(fOut, "%12.5lf%12.5lf%12.5lf%12.5lf\n\n", 
			E_CMap[Count  ], E_CMap[Count+1], E_CMap[Count+2], E_CMap[Count+3]);
		Count+=4;
	}
	fclose(fOut);

}


#define Soft_k_Constrain_Torsion	(0.6)
#define FLAT_PHI					(15.0)
void CMol::Cal_E_Dihedral(int ToSave_Phi)
{
	int iActive, i, iPos, ia, ib, ic, id, Idx_n;
	double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
	double xba, yba, zba, xca, yca, zca, xcb, ycb, zcb, xdc, ydc, zdc, xdb, ydb, zdb, rcb;
	double xt, yt, zt, xu, yu, zu, xtu, ytu, ztu, rt2, ru2, rtru;
	double cosine, sine, fai, d_fai_In, d_phi;
	double dedphi;
	double dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
	double dedxt, dedyt, dedzt, dedxu, dedyu, dedzu;
	
	E_Dihedral=0.0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_Dihedral; iActive++)	{
			i = Active_List_Dihedral[iActive];	//the index of dihedral
			iPos = 4*i;
			
			ia = DihedralList[iPos  ];
			ib = DihedralList[iPos+1];
			ic = DihedralList[iPos+2];
			id = DihedralList[iPos+3];
			
			xia = x[ia];		yia = y[ia];		zia = z[ia];
			xib = x[ib];		yib = y[ib];		zib = z[ib];
			xic = x[ic];		yic = y[ic];		zic = z[ic];
			xid = x[id];		yid = y[id];		zid = z[id];
			xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
			
			xt = yba*zcb - ycb*zba;
			yt = zba*xcb - zcb*xba;
			zt = xba*ycb - xcb*yba;
			xu = ycb*zdc - ydc*zcb;
			yu = zcb*xdc - zdc*xcb;
			zu = xcb*ydc - xdc*ycb;
			xtu = yt*zu - yu*zt;
			ytu = zt*xu - zu*xt;
			ztu = xt*yu - xu*yt;
			rt2 = xt*xt + yt*yt + zt*zt;
			ru2 = xu*xu + yu*yu + zu*zu;
			rt2=max(rt2, 0.000000001);
			ru2=max(ru2, 0.000000001);
			rtru = sqrt(rt2 * ru2);
			
			if(rtru!=0.0)	{
				rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
				cosine = (xt*xu + yt*yu + zt*zu) / rtru;
				sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
				
				
				fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
				if(cosine < 0.0)	{	//[90,270]
					if(sine < 0.0)	{	//[-180,-90]
						fai = -180.0 - fai;
					}
					else	{	//[90,180]
						fai = 180.0 - fai;
					}
				}
				dedphi = 0.0;


				if(ToSave_Phi)	{	// to save the initial phi into All_Phi0 used for constrain
					All_Phi0[i] = fai;
				}
				
				if(Soft_Restrain_Dihedrals && Is_Phi_Constrained[i] )	{	// to apply soft restraints on all torsions. 
					d_phi = (fai - All_Phi0[i]);
					if(d_phi > 180.0)	{	//period is 360
						d_phi -= 360.0;
					}
					else if(d_phi < -180.0)	{
						d_phi += 360.0;
					}
					if(d_phi < -FLAT_PHI)	{
						d_phi += FLAT_PHI;
					}
					else if(d_phi > FLAT_PHI)	{
						d_phi -= FLAT_PHI;
					}
					else	{
						d_phi = 0.0;
					}
					d_phi *= radianInv;
					dedphi += (Soft_k_Constrain_Torsion*d_phi);
					E_Dihedral += (0.5*Soft_k_Constrain_Torsion*d_phi*d_phi);
				}
				
				fai *= radianInv;
				for(Idx_n=1; Idx_n<MAX_DIH_ITEM; Idx_n++)	{
					if(fabs(Para_k_Dih[i][Idx_n]) < 1.0E-50)	{	//zero
						continue;
					}
					
					d_fai_In = Idx_n*fai - Para_phi[i][Idx_n];	// n*phi - phi_0
					E_Dihedral += (Para_k_Dih[i][Idx_n] * (1.0 + cos(d_fai_In)));
					
					dedphi -= ( Idx_n * Para_k_Dih[i][Idx_n] * sin(d_fai_In) );	//accumulate all terms
				}
				
				xca = xic - xia;
				yca = yic - yia;
				zca = zic - zia;
				xdb = xid - xib;
				ydb = yid - yib;
				zdb = zid - zib;
				
				dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
				dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
				dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
				dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
				dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
				dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
				
				dedxia = zcb*dedyt - ycb*dedzt;
				dedyia = xcb*dedzt - zcb*dedxt;
				dedzia = ycb*dedxt - xcb*dedyt;
				dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
				dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
				dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
				dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
				dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
				dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
				dedxid = zcb*dedyu - ycb*dedzu;
				dedyid = xcb*dedzu - zcb*dedxu;
				dedzid = ycb*dedxu - xcb*dedyu;
				
				grad_x[ia] += dedxia;
				grad_y[ia] += dedyia;
				grad_z[ia] += dedzia;
				grad_x[ib] += dedxib;
				grad_y[ib] += dedyib;
				grad_z[ib] += dedzib;
				grad_x[ic] += dedxic;
				grad_y[ic] += dedyic;
				grad_z[ic] += dedzic;
				grad_x[id] += dedxid;
				grad_y[id] += dedyid;
				grad_z[id] += dedzid;
			}
			else	{
				//			printf("Error in calculation of dihedral angle.\n");
				fprintf(fFile_Run_Log, "Error in calculation of dihedral angle.\n");
				fflush(fFile_Run_Log);
			}
		}
	}
	else	{
		for(iActive=0; iActive<nDihedral; iActive++)	{
			i = iActive;	//the index of dihedral
			iPos = 4*i;
			
			ia = DihedralList[iPos  ];
			ib = DihedralList[iPos+1];
			ic = DihedralList[iPos+2];
			id = DihedralList[iPos+3];
			
			xia = x[ia];		yia = y[ia];		zia = z[ia];
			xib = x[ib];		yib = y[ib];		zib = z[ib];
			xic = x[ic];		yic = y[ic];		zic = z[ic];
			xid = x[id];		yid = y[id];		zid = z[id];
			xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
			
			xt = yba*zcb - ycb*zba;
			yt = zba*xcb - zcb*xba;
			zt = xba*ycb - xcb*yba;
			xu = ycb*zdc - ydc*zcb;
			yu = zcb*xdc - zdc*xcb;
			zu = xcb*ydc - xdc*ycb;
			xtu = yt*zu - yu*zt;
			ytu = zt*xu - zu*xt;
			ztu = xt*yu - xu*yt;
			rt2 = xt*xt + yt*yt + zt*zt;
			ru2 = xu*xu + yu*yu + zu*zu;
			rt2=max(rt2, 0.000000001);
			ru2=max(ru2, 0.000000001);
			rtru = sqrt(rt2 * ru2);
			
			if(rtru!=0.0)	{
				rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
				cosine = (xt*xu + yt*yu + zt*zu) / rtru;
				sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
				
				
				fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
				if(cosine < 0.0)	{	//[90,270]
					if(sine < 0.0)	{	//[-180,-90]
						fai = -180.0 - fai;
					}
					else	{	//[90,180]
						fai = 180.0 - fai;
					}
				}
				
				fai *= radianInv;
				dedphi = 0.0;
				for(Idx_n=1; Idx_n<MAX_DIH_ITEM; Idx_n++)	{
					if(fabs(Para_k_Dih[i][Idx_n]) < 1.0E-50)	{	//zero
						continue;
					}
					
					d_fai_In = Idx_n*fai - Para_phi[i][Idx_n];	// n*phi - phi_0
					E_Dihedral += (Para_k_Dih[i][Idx_n] * (1.0 + cos(d_fai_In)));
					
					dedphi -= ( Idx_n * Para_k_Dih[i][Idx_n] * sin(d_fai_In) );	//accumulate all terms
				}
				
				xca = xic - xia;
				yca = yic - yia;
				zca = zic - zia;
				xdb = xid - xib;
				ydb = yid - yib;
				zdb = zid - zib;
				
				dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
				dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
				dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
				dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
				dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
				dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
				
				dedxia = zcb*dedyt - ycb*dedzt;
				dedyia = xcb*dedzt - zcb*dedxt;
				dedzia = ycb*dedxt - xcb*dedyt;
				dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
				dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
				dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
				dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
				dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
				dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
				dedxid = zcb*dedyu - ycb*dedzu;
				dedyid = xcb*dedzu - zcb*dedxu;
				dedzid = ycb*dedxu - xcb*dedyu;
				
				grad_x[ia] += dedxia;
				grad_y[ia] += dedyia;
				grad_z[ia] += dedzia;
				grad_x[ib] += dedxib;
				grad_y[ib] += dedyib;
				grad_z[ib] += dedzib;
				grad_x[ic] += dedxic;
				grad_y[ic] += dedyic;
				grad_z[ic] += dedzic;
				grad_x[id] += dedxid;
				grad_y[id] += dedyid;
				grad_z[id] += dedzid;
			}
			else	{
				//			printf("Error in calculation of dihedral angle.\n");
				fprintf(fFile_Run_Log, "Error in calculation of dihedral angle.\n");
				fflush(fFile_Run_Log);
			}
		}
	}

	
	return;
}


#define k_Constrain_Phi	(10000.0)
void CMol::Cal_E_Constrain_Phi(void)
{
	int iActive;
	int ia, ib, ic, id;
	double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
	double xba, yba, zba, xca, yca, zca, xcb, ycb, zcb, xdc, ydc, zdc, xdb, ydb, zdb, rcb;
	double xt, yt, zt, xu, yu, zu, xtu, ytu, ztu, rt2, ru2, rtru;
	double cosine, sine, fai;
	double dedphi, dPhi;
	double dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
	double dedxt, dedyt, dedzt, dedxu, dedyu, dedzu;
	
	E_Constrain_Phi=0.0;
	
	for(iActive=0; iActive<nActive_CMap; iActive++)	{
		ia = CMapList[iActive][0];
		ib = CMapList[iActive][1];
		ic = CMapList[iActive][2];
		id = CMapList[iActive][3];
		
		xia = x[ia];		yia = y[ia];		zia = z[ia];
		xib = x[ib];		yib = y[ib];		zib = z[ib];
		xic = x[ic];		yic = y[ic];		zic = z[ic];
		xid = x[id];		yid = y[id];		zid = z[id];
		xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
		xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
		xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
		
		xt = yba*zcb - ycb*zba;
		yt = zba*xcb - zcb*xba;
		zt = xba*ycb - xcb*yba;
		xu = ycb*zdc - ydc*zcb;
		yu = zcb*xdc - zdc*xcb;
		zu = xcb*ydc - xdc*ycb;
		xtu = yt*zu - yu*zt;
		ytu = zt*xu - zu*xt;
		ztu = xt*yu - xu*yt;
		rt2 = xt*xt + yt*yt + zt*zt;
		ru2 = xu*xu + yu*yu + zu*zu;
		rt2=max(rt2, 0.000000001);
		ru2=max(ru2, 0.000000001);
		rtru = sqrt(rt2 * ru2);
		
		if(rtru!=0.0)	{
			rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
			cosine = (xt*xu + yt*yu + zt*zu) / rtru;
			sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
			
			fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
			if(cosine < 0.0)	{	//[90,270]
				if(sine < 0.0)	{	//[-180,-90]
					fai = -180.0 - fai;
				}
				else	{	//[90,180]
					fai = 180.0 - fai;
				}
			}
			
			Phi_Real = fai;
			
			dPhi = (fai - Phi0_Constrain);
			if(dPhi > 180.0)	{	//period is 360
				dPhi -= 360.0;
			}
			else if(dPhi < -180.0)	{
				dPhi += 360.0;
			}
			dPhi *= radianInv;
			E_Constrain_Phi = 0.5*k_Constrain_Phi*dPhi*dPhi;
			dedphi = k_Constrain_Phi*dPhi;
			
			xca = xic - xia;
			yca = yic - yia;
			zca = zic - zia;
			xdb = xid - xib;
			ydb = yid - yib;
			zdb = zid - zib;
			
			dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
			dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
			dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
			dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
			dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
			dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
			
			dedxia = zcb*dedyt - ycb*dedzt;
			dedyia = xcb*dedzt - zcb*dedxt;
			dedzia = ycb*dedxt - xcb*dedyt;
			dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
			dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
			dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
			dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
			dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
			dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
			dedxid = zcb*dedyu - ycb*dedzu;
			dedyid = xcb*dedzu - zcb*dedxu;
			dedzid = ycb*dedxu - xcb*dedyu;
			
			grad_x[ia] += dedxia;
			grad_y[ia] += dedyia;
			grad_z[ia] += dedzia;
			grad_x[ib] += dedxib;
			grad_y[ib] += dedyib;
			grad_z[ib] += dedzib;
			grad_x[ic] += dedxic;
			grad_y[ic] += dedyic;
			grad_z[ic] += dedzic;
			grad_x[id] += dedxid;
			grad_y[id] += dedyid;
			grad_z[id] += dedzid;
		}
		else	{
//			printf("Error in Cal_E_Constrain_Phi().\n");
			fprintf(fFile_Run_Log, "Error in Cal_E_Constrain_Phi().\n");
			fflush(fFile_Run_Log);
		}
	}	
	
	
	return;
}
#undef k_Constrain_Phi


#define k_Constrain_Psi	(10000.0)
void CMol::Cal_E_Constrain_Psi(void)	//design for alad only
{
	int iActive;
	int ia, ib, ic, id;
	double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
	double xba, yba, zba, xca, yca, zca, xcb, ycb, zcb, xdc, ydc, zdc, xdb, ydb, zdb, rcb;
	double xt, yt, zt, xu, yu, zu, xtu, ytu, ztu, rt2, ru2, rtru;
	double cosine, sine, fai;
	double dedphi, dPhi;
	double dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
	double dedxt, dedyt, dedzt, dedxu, dedyu, dedzu;
	
	E_Constrain_Psi=0.0;
	
	for(iActive=0; iActive<nActive_CMap; iActive++)	{
		ia = CMapList[iActive][4];
		ib = CMapList[iActive][5];
		ic = CMapList[iActive][6];
		id = CMapList[iActive][7];
		
		xia = x[ia];		yia = y[ia];		zia = z[ia];
		xib = x[ib];		yib = y[ib];		zib = z[ib];
		xic = x[ic];		yic = y[ic];		zic = z[ic];
		xid = x[id];		yid = y[id];		zid = z[id];
		xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
		xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
		xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
		
		xt = yba*zcb - ycb*zba;
		yt = zba*xcb - zcb*xba;
		zt = xba*ycb - xcb*yba;
		xu = ycb*zdc - ydc*zcb;
		yu = zcb*xdc - zdc*xcb;
		zu = xcb*ydc - xdc*ycb;
		xtu = yt*zu - yu*zt;
		ytu = zt*xu - zu*xt;
		ztu = xt*yu - xu*yt;
		rt2 = xt*xt + yt*yt + zt*zt;
		ru2 = xu*xu + yu*yu + zu*zu;
		rt2=max(rt2, 0.000000001);
		ru2=max(ru2, 0.000000001);
		rtru = sqrt(rt2 * ru2);
		
		if(rtru!=0.0)	{
			rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
			cosine = (xt*xu + yt*yu + zt*zu) / rtru;
			sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
			
			fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
			if(cosine < 0.0)	{	//[90,270]
				if(sine < 0.0)	{	//[-180,-90]
					fai = -180.0 - fai;
				}
				else	{	//[90,180]
					fai = 180.0 - fai;
				}
			}
			Psi_Real = fai;
			
			dPhi = (fai - Psi0_Constrain);
			if(dPhi > 180.0)	{	//period is 360
				dPhi -= 360.0;
			}
			else if(dPhi < -180.0)	{
				dPhi += 360.0;
			}
			dPhi *= radianInv;
			
			E_Constrain_Psi = 0.5*k_Constrain_Psi*dPhi*dPhi;
			dedphi = k_Constrain_Psi*dPhi;
			
			xca = xic - xia;
			yca = yic - yia;
			zca = zic - zia;
			xdb = xid - xib;
			ydb = yid - yib;
			zdb = zid - zib;
			
			dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
			dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
			dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
			dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
			dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
			dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
			
			dedxia = zcb*dedyt - ycb*dedzt;
			dedyia = xcb*dedzt - zcb*dedxt;
			dedzia = ycb*dedxt - xcb*dedyt;
			dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
			dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
			dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
			dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
			dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
			dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
			dedxid = zcb*dedyu - ycb*dedzu;
			dedyid = xcb*dedzu - zcb*dedxu;
			dedzid = ycb*dedxu - xcb*dedyu;
			
			grad_x[ia] += dedxia;
			grad_y[ia] += dedyia;
			grad_z[ia] += dedzia;
			grad_x[ib] += dedxib;
			grad_y[ib] += dedyib;
			grad_z[ib] += dedzib;
			grad_x[ic] += dedxic;
			grad_y[ic] += dedyic;
			grad_z[ic] += dedzic;
			grad_x[id] += dedxid;
			grad_y[id] += dedyid;
			grad_z[id] += dedzid;
		}
		else	{
//			printf("Error in Cal_E_Constrain_Psi().\n");
			fprintf(fFile_Run_Log, "Error in Cal_E_Constrain_Psi().\n");
			fflush(fFile_Run_Log);
		}
	}
	
	return;
}
#undef k_Constrain_Psi


void CMol::Cal_E_ImproperDihedral(void)
{
	int iActive, i, iPos, ia, ib, ic, id;
	double xia, yia, zia, xib, yib, zib, xic, yic, zic, xid, yid, zid;
	double xba, yba, zba, xca, yca, zca, xcb, ycb, zcb, xdc, ydc, zdc, xdb, ydb, zdb, rcb;
	double xt, yt, zt, xu, yu, zu, xtu, ytu, ztu, rt2, ru2, rtru;
	double cosine, sine, fai, dFai;
	double dedphi;
	double dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic, dedyic, dedzic, dedxid, dedyid, dedzid;
	double dedxt, dedyt, dedzt, dedxu, dedyu, dedzu;
	
	E_ImproDihedral = 0.0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_ImpDih; iActive++)	{
			i = Active_List_ImpDih[iActive];	//the index of improper dihedral
			iPos = 4*i;
			
			ia = ImprDihedralList[iPos  ];
			ib = ImprDihedralList[iPos+1];
			ic = ImprDihedralList[iPos+2];
			id = ImprDihedralList[iPos+3];
			
			xia = x[ia];		yia = y[ia];		zia = z[ia];
			xib = x[ib];		yib = y[ib];		zib = z[ib];
			xic = x[ic];		yic = y[ic];		zic = z[ic];
			xid = x[id];		yid = y[id];		zid = z[id];
			xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
			
			xt = yba*zcb - ycb*zba;
			yt = zba*xcb - zcb*xba;
			zt = xba*ycb - xcb*yba;
			xu = ycb*zdc - ydc*zcb;
			yu = zcb*xdc - zdc*xcb;
			zu = xcb*ydc - xdc*ycb;
			xtu = yt*zu - yu*zt;
			ytu = zt*xu - zu*xt;
			ztu = xt*yu - xu*yt;
			rt2 = xt*xt + yt*yt + zt*zt;
			ru2 = xu*xu + yu*yu + zu*zu;
			rt2=max(rt2, 0.000000001);
			ru2=max(ru2, 0.000000001);
			rtru = sqrt(rt2 * ru2);
			
			if(rtru!=0.0)	{
				rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
				cosine = (xt*xu + yt*yu + zt*zu) / rtru;
				sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
				
				
				fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
				if(cosine < 0.0)	{	//[90,270]
					if(sine < 0.0)	{	//[-180,-90]
						fai = -180.0 - fai;
					}
					else	{	//[90,180]
						fai = 180.0 - fai;
					}
				}
				
				fai *= radianInv;
				
				if((int)(Para_Type_ImpDih[i]+1.0E-20) == 0)	{	// 0 - CHARMM potential
					dFai = fai - Para_Imp_phi[i];
					
					//start	to make sure the absolute value of dFai is smaller than PI
					if(dFai > PI)	{	// period 2*PI
						dFai -= PI2;
					}
					else if(dFai < -PI)	{
						dFai += PI2;
					}
					//end	to make sure the absolute value of dFai is smaller than PI
					
					E_ImproDihedral += (Para_k_ImpDih[i] * dFai * dFai);	// E = k * (phi - phi_0)^2
					dedphi = 2.0 * Para_k_ImpDih[i] * dFai;	//accumulate all terms
				}
				else if((int)(Para_Type_ImpDih[i]+1.0E-20) == 2)	{	// 2 - AMBER's potential. The same form as dihedral potential
					dFai = Para_Type_ImpDih[i]*fai - Para_Imp_phi[i];	// n*phi - phi_0
					E_ImproDihedral += (Para_k_ImpDih[i] * (1.0 + cos(dFai)));
					dedphi = -( Para_Type_ImpDih[i] * Para_k_ImpDih[i] * sin(dFai) );	//accumulate all terms
				}
				else	{
					fprintf(fFile_Run_Log, "Unknown improper angle parameters. Only 0 (CHARMM) and 2 (AMBER) types are supported.\nQuit\n");
					fflush(fFile_Run_Log);
					exit(1);
				}
				
				
				xca = xic - xia;
				yca = yic - yia;
				zca = zic - zia;
				xdb = xid - xib;
				ydb = yid - yib;
				zdb = zid - zib;
				
				dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
				dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
				dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
				dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
				dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
				dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
				
				dedxia = zcb*dedyt - ycb*dedzt;
				dedyia = xcb*dedzt - zcb*dedxt;
				dedzia = ycb*dedxt - xcb*dedyt;
				dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
				dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
				dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
				dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
				dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
				dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
				dedxid = zcb*dedyu - ycb*dedzu;
				dedyid = xcb*dedzu - zcb*dedxu;
				dedzid = ycb*dedxu - xcb*dedyu;
				
				grad_x[ia] += dedxia;
				grad_y[ia] += dedyia;
				grad_z[ia] += dedzia;
				grad_x[ib] += dedxib;
				grad_y[ib] += dedyib;
				grad_z[ib] += dedzib;
				grad_x[ic] += dedxic;
				grad_y[ic] += dedyic;
				grad_z[ic] += dedzic;
				grad_x[id] += dedxid;
				grad_y[id] += dedyid;
				grad_z[id] += dedzid;
			}
			else	{
				//			printf("Error in calculation of improper dihedral.\n");
				fprintf(fFile_Run_Log, "Error in calculation of improper dihedral.\n");
				fflush(fFile_Run_Log);
			}
		}
	}
	else	{
		for(iActive=0; iActive<nImpro; iActive++)	{
			i = iActive;	//the index of improper dihedral
			iPos = 4*i;
			
			ia = ImprDihedralList[iPos  ];
			ib = ImprDihedralList[iPos+1];
			ic = ImprDihedralList[iPos+2];
			id = ImprDihedralList[iPos+3];
			
			xia = x[ia];		yia = y[ia];		zia = z[ia];
			xib = x[ib];		yib = y[ib];		zib = z[ib];
			xic = x[ic];		yic = y[ic];		zic = z[ic];
			xid = x[id];		yid = y[id];		zid = z[id];
			xba = xib - xia;	yba = yib - yia;	zba = zib - zia;
			xcb = xic - xib;	ycb = yic - yib;	zcb = zic - zib;
			xdc = xid - xic;	ydc = yid - yic;	zdc = zid - zic;
			
			xt = yba*zcb - ycb*zba;
			yt = zba*xcb - zcb*xba;
			zt = xba*ycb - xcb*yba;
			xu = ycb*zdc - ydc*zcb;
			yu = zcb*xdc - zdc*xcb;
			zu = xcb*ydc - xdc*ycb;
			xtu = yt*zu - yu*zt;
			ytu = zt*xu - zu*xt;
			ztu = xt*yu - xu*yt;
			rt2 = xt*xt + yt*yt + zt*zt;
			ru2 = xu*xu + yu*yu + zu*zu;
			rt2=max(rt2, 0.000000001);
			ru2=max(ru2, 0.000000001);
			rtru = sqrt(rt2 * ru2);
			
			if(rtru!=0.0)	{
				rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
				cosine = (xt*xu + yt*yu + zt*zu) / rtru;
				sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);		//as in Tinker
				
				
				fai = asin(min(max(sine,-1.0),1.0))*radian;	// [-90,90]
				if(cosine < 0.0)	{	//[90,270]
					if(sine < 0.0)	{	//[-180,-90]
						fai = -180.0 - fai;
					}
					else	{	//[90,180]
						fai = 180.0 - fai;
					}
				}
				
				fai *= radianInv;
				
				if((int)(Para_Type_ImpDih[i]+1.0E-20) == 0)	{	// 0 - CHARMM potential
					dFai = fai - Para_Imp_phi[i];
					
					//start	to make sure the absolute value of dFai is smaller than PI
					if(dFai > PI)	{	// period 2*PI
						dFai -= PI2;
					}
					else if(dFai < -PI)	{
						dFai += PI2;
					}
					//end	to make sure the absolute value of dFai is smaller than PI
					
					E_ImproDihedral += (Para_k_ImpDih[i] * dFai * dFai);	// E = k * (phi - phi_0)^2
					dedphi = 2.0 * Para_k_ImpDih[i] * dFai;	//accumulate all terms
				}
				else if((int)(Para_Type_ImpDih[i]+1.0E-20) == 2)	{	// 2 - AMBER's potential. The same form as dihedral potential
					dFai = Para_Type_ImpDih[i]*fai - Para_Imp_phi[i];	// n*phi - phi_0
					E_ImproDihedral += (Para_k_ImpDih[i] * (1.0 + cos(dFai)));
					dedphi = -( Para_Type_ImpDih[i] * Para_k_ImpDih[i] * sin(dFai) );	//accumulate all terms
				}
				else	{
					fprintf(fFile_Run_Log, "Unknown improper angle parameters. Only 0 (CHARMM) and 2 (AMBER) types are supported.\nQuit\n");
					fflush(fFile_Run_Log);
					exit(1);
				}
				
				
				xca = xic - xia;
				yca = yic - yia;
				zca = zic - zia;
				xdb = xid - xib;
				ydb = yid - yib;
				zdb = zid - zib;
				
				dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
				dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
				dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
				dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
				dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
				dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
				
				dedxia = zcb*dedyt - ycb*dedzt;
				dedyia = xcb*dedzt - zcb*dedxt;
				dedzia = ycb*dedxt - xcb*dedyt;
				dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
				dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
				dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
				dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
				dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
				dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
				dedxid = zcb*dedyu - ycb*dedzu;
				dedyid = xcb*dedzu - zcb*dedxu;
				dedzid = ycb*dedxu - xcb*dedyu;
				
				grad_x[ia] += dedxia;
				grad_y[ia] += dedyia;
				grad_z[ia] += dedzia;
				grad_x[ib] += dedxib;
				grad_y[ib] += dedyib;
				grad_z[ib] += dedzib;
				grad_x[ic] += dedxic;
				grad_y[ic] += dedyic;
				grad_z[ic] += dedzic;
				grad_x[id] += dedxid;
				grad_y[id] += dedyid;
				grad_z[id] += dedzid;
			}
			else	{
				//			printf("Error in calculation of improper dihedral.\n");
				fprintf(fFile_Run_Log, "Error in calculation of improper dihedral.\n");
				fflush(fFile_Run_Log);
			}
		}
	}
	

	return;
}

void CMol::Cal_E_VDW_ELEC(void)	//non-bonded interactions
{
	int iActive, iPair, i, j;
	double Rij_SQ, Rij, dxij, dyij, dzij, R6, R12;
	double drdEdr, dedxij, dedyij, dedzij, dE;
	
	E_VDW = E_Elec = 0.0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_NonBond; iActive++)	{
			iPair = Active_List_NonBond[iActive];	//the index of non-bonded pair
			
			i = NB_List_i[iPair];
			j = NB_List_j[iPair];

			if(j >= nAtom)	{
				continue;
			}
			
			dxij=x[j]-x[i];
			dyij=y[j]-y[i];
			dzij=z[j]-z[i];
			
			Rij_SQ=dxij*dxij+dyij*dyij+dzij*dzij;
			
			if(Rij_SQ > Dist_Cutoff_SQ)	{	// two atoms are too far, just skip. Designed for molecule molecule interactions. 
				continue;
			}
		
				
			//		if(Para_LJ_Epsilon_IJ[iPair] > 1.0E-50)	{	//not zero
			//start	to calculate VDW energy and gradient
			R6=Rij_SQ*Rij_SQ*Rij_SQ;
			R12=R6*R6;
			
			
			dE = Para_LJ_Epsilon_IJ[iPair] * ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
			
			// R_12 is turned off here !!!
			//			dE = Para_LJ_Epsilon_IJ[iPair] * ( - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
			E_VDW += (dE);
			
			//			printf("E_VDW:  %3d %3d  %20.15lf\n", i+1, j+1, dE);
			
			drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( Para_LJ_Sigma_IJ_pow_12[iPair]/(R12*Rij_SQ) - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
			
			// R_12 is turned off here !!!
			//			drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
			
			dedxij=drdEdr*dxij;
			dedyij=drdEdr*dyij;
			dedzij=drdEdr*dzij;
			
			
			grad_x[i] += dedxij;
			grad_y[i] += dedyij;
			grad_z[i] += dedzij;
			
			grad_x[j] -= dedxij;
			grad_y[j] -= dedyij;
			grad_z[j] -= dedzij;
			//end	to calculate VDW energy and gradient
			//		}
			
			
			//start	to calculate elec energy and gradient
			Rij = sqrt(Rij_SQ);
			dE = Para_Elec_Pair[iPair] / Rij;
			E_Elec += dE;
			
			//		printf("E_Elec:  %3d  %3d  %25.18lf\n", i+1, j+1, dE);
			
			drdEdr = Para_Elec_Pair[iPair] / (Rij*Rij_SQ);
			
			dedxij=drdEdr*dxij;
			dedyij=drdEdr*dyij;
			dedzij=drdEdr*dzij;
			
			
			grad_x[i] += dedxij;
			grad_y[i] += dedyij;
			grad_z[i] += dedzij;
			
			grad_x[j] -= dedxij;
			grad_y[j] -= dedyij;
			grad_z[j] -= dedzij;
			//end	to calculate elec energy and gradient
		}
	}
	else	{
		for(iActive=0; iActive<n_NB_Pair; iActive++)	{
			iPair = iActive;	//the index of non-bonded pair
			
			i = NB_List_i[iPair];
			j = NB_List_j[iPair];

			if(j >= nAtom)	{
				continue;
			}
						
			dxij=x[j]-x[i];
			dyij=y[j]-y[i];
			dzij=z[j]-z[i];
			
			Rij_SQ=dxij*dxij+dyij*dyij+dzij*dzij;
			
			if(Rij_SQ > Dist_Cutoff_SQ)	{	// two atoms are too far, just skip. Designed for molecule molecule interactions. 
				continue;
			}
			
			//		if(Para_LJ_Epsilon_IJ[iPair] > 1.0E-50)	{	//not zero
			//start	to calculate VDW energy and gradient
			R6=Rij_SQ*Rij_SQ*Rij_SQ;
			R12=R6*R6;
			
			
			dE = Para_LJ_Epsilon_IJ[iPair] * ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
			
			// R_12 is turned off here !!!
			//			dE = Para_LJ_Epsilon_IJ[iPair] * ( - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
			E_VDW += (dE);
			
			//			printf("E_VDW:  %3d %3d  %20.15lf\n", i+1, j+1, dE);
			
			drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( Para_LJ_Sigma_IJ_pow_12[iPair]/(R12*Rij_SQ) - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
			
			// R_12 is turned off here !!!
			//			drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
			
			dedxij=drdEdr*dxij;
			dedyij=drdEdr*dyij;
			dedzij=drdEdr*dzij;
			
			
			grad_x[i] += dedxij;
			grad_y[i] += dedyij;
			grad_z[i] += dedzij;
			
			grad_x[j] -= dedxij;
			grad_y[j] -= dedyij;
			grad_z[j] -= dedzij;
			//end	to calculate VDW energy and gradient
			//		}
			
			
			//start	to calculate elec energy and gradient
			Rij = sqrt(Rij_SQ);
			dE = Para_Elec_Pair[iPair] / Rij;
			E_Elec += dE;
			
			//		printf("E_Elec:  %3d  %3d  %25.18lf\n", i+1, j+1, dE);
			
			drdEdr = Para_Elec_Pair[iPair] / (Rij*Rij_SQ);
			
			dedxij=drdEdr*dxij;
			dedyij=drdEdr*dyij;
			dedzij=drdEdr*dzij;
			
			
			grad_x[i] += dedxij;
			grad_y[i] += dedyij;
			grad_z[i] += dedzij;
			
			grad_x[j] -= dedxij;
			grad_y[j] -= dedyij;
			grad_z[j] -= dedzij;
			//end	to calculate elec energy and gradient
		}
	}

	
	return;
}

void CMol::Cal_dE_d_VDW_Param(void)
{
	int iActive, iPair, i, j, LJ_Type_i, LJ_Type_j;
	double Rij_SQ, dxij, dyij, dzij, R6, R12;
	double dE_d_Emin, dE_d_Rmin;
	double E_min_i, E_min_j, R_min_i, R_min_j;

	memset(dE_d_LJEmin, 0, sizeof(double)*myForceField->n_Rec_LJ);
	memset(dE_d_LJRmin, 0, sizeof(double)*myForceField->n_Rec_LJ);

	//int count = 0;

	// this is computing LJ derivatives of all types i and j by summing over all NB pairs i.e 1-4 and beyond only
	//printf("%s %d\n","n_NB_Pair",n_NB_Pair);
	for(iActive=0; iActive<n_NB_Pair; iActive++)	{
		iPair = iActive;	//the index of non-bonded pair
		
		i = NB_List_i[iPair];
		j = NB_List_j[iPair];

		//printf("%s %s\n","n_NB_Pair_i",ChemName[i]);
		//printf("%s %s\n","n_NB_Pair_j",ChemName[j]);

		if(j >= nAtom)	{
			continue;
		}


		//printf("%s %i \n","count",count);

		E_min_i = Para_LJ_Epsilon[i];
		E_min_j = Para_LJ_Epsilon[j];
		R_min_i = Para_LJ_Sigma[i];
		R_min_j = Para_LJ_Sigma[j];


		//printf("%s %f \n","E_min_i",E_min_i);
		//printf("%s %f \n","E_min_j",E_min_j);
		//printf("%s %f \n","R_min_i",R_min_i);
		//printf("%s %f \n","R_min_j",R_min_j);

		
		if( (fabs(E_min_i) < 1.0E-10) || (fabs(E_min_j) < 1.0E-10) )	{	// any one is zero, then skip
			continue;
		}
	
		if(Para_LJ_Sqrt_Epsilon[i] < 1.0E-20)	{
			continue;
		}
		LJ_Type_i = LJ_Para_Rec[i];
		if(Para_LJ_Sqrt_Epsilon[j] < 1.0E-20)	{
			continue;
		}
		LJ_Type_j = LJ_Para_Rec[j];

			
		//printf("%s %s %f %f \n","xi","xj",x[i],x[j]);
		//printf("%s %s %f %f \n","yi","yj",y[i],y[j]);
		//printf("%s %s %f %f \n","zi","zj",z[i],z[j]);
		
		dxij=x[j]-x[i];
		dyij=y[j]-y[i];
		dzij=z[j]-z[i];
		
		Rij_SQ=dxij*dxij+dyij*dyij+dzij*dzij;
		
		if(Rij_SQ > Dist_Cutoff_SQ)	{	// two atoms are too far, just skip. Designed for molecule molecule interactions. 
			continue;
		}
	
		//Chetan: This was needed for drude lj opt runs 
		//if(Rij_SQ <=0.0) {
		//	continue;
		//}
	
		//start	to calculate VDW energy and gradient
		R6=Rij_SQ*Rij_SQ*Rij_SQ;
		R12=R6*R6;

		//printf("%s %f\n","Rij_SQ",Rij_SQ);
		//printf("%s %f %s %f \n","R6",R6,"R12",R12);
		//printf("%s %f \n","Para_LJ_Sigma_IJ_pow_12[iPair]",Para_LJ_Sigma_IJ_pow_12[iPair]);
		//printf("%s %f \n","Para_LJ_Sigma_IJ_pow_6[iPair]",Para_LJ_Sigma_IJ_pow_6[iPair]);
	
		if(DistMatrix[i][j] == 3)	{	// 1-4 pairs
			dE_d_Emin = ( ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) ) * 0.5);
			dE_d_LJEmin[LJ_Type_i] += (dE_d_Emin * Para_LJ_Sqrt_Epsilon_14[j] * Para_LJ_Epsilon_14[i] / (Para_LJ_Sqrt_Epsilon_14[i] * Para_LJ_Epsilon[i]) );
			dE_d_LJEmin[LJ_Type_j] += (dE_d_Emin * Para_LJ_Sqrt_Epsilon_14[i] * Para_LJ_Epsilon_14[j] / (Para_LJ_Sqrt_Epsilon_14[j] * Para_LJ_Epsilon[j]));
		}
		else	{
			dE_d_Emin = ( ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) ) * 0.5);
			dE_d_LJEmin[LJ_Type_i] += (dE_d_Emin * Para_LJ_Sqrt_Epsilon[j] / Para_LJ_Sqrt_Epsilon[i]);
			dE_d_LJEmin[LJ_Type_j] += (dE_d_Emin * Para_LJ_Sqrt_Epsilon[i] / Para_LJ_Sqrt_Epsilon[j]);
		}
		

		//printf("%s %i \n","count",count);

		//printf("%s %f \n","dE_d_Emin",dE_d_Emin);
		//printf("%s %s %f \n","dE_d_LJEmin_i",ChemName[i],dE_d_LJEmin[LJ_Type_i]);
		//printf("%s %s %f \n","dE_d_LJEmin_j",ChemName[j],dE_d_LJEmin[LJ_Type_j]);
		


		dE_d_Rmin = ( Para_LJ_Epsilon_IJ[iPair] * ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - (Para_LJ_Sigma_IJ_pow_6[iPair]/R6) ) * 12.0 / (Para_LJ_Sigma_IJ[iPair]));
		//dE_d_Rmin = ( Para_LJ_Epsilon_IJ[iPair] * ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - (Para_LJ_Sigma_IJ_pow_6[iPair]/(2*R6)) ) * 12.0 / (Para_LJ_Sigma_IJ[iPair]));
		dE_d_LJRmin[LJ_Type_i] += dE_d_Rmin;
		dE_d_LJRmin[LJ_Type_j] += dE_d_Rmin;
		
		//printf("%s %f \n","dE_d_Rmin",dE_d_Rmin);	
		//printf("%s %s %f \n","dE_d_LJRmin_i",ChemName[i],dE_d_LJRmin[LJ_Type_i]);	
		//printf("%s %s %f \n","dE_d_LJRmin_j",ChemName[j],dE_d_LJRmin[LJ_Type_j]);	
		//exit(0);

		//count = count+1;

	}

	//exit(0);

	//accumulated over all the number_of_copies_of_a_molecule*(MDSteps/num_md_steps_taken_to_call_derivative_calc) for each LJ Type
	//notice that the sum is over each "i" that is done separately for each LJ type

	for(i=0; i<myForceField->n_Rec_LJ; i++)	{
		//printf("%s%i %f \n","dE_d_LJEmin",i,dE_d_LJEmin[i]);
		//printf("%s%i %f \n","dE_d_LJRmin",i,dE_d_LJRmin[i]);
		//printf("%s%i %f \n","dE_d_LJEmin",i,dE_d_LJEmin_Acc[i]);
		//printf("%s%i %f \n","dE_d_LJRmin",i,dE_d_LJRmin_Acc[i]);
		//printf("%s%i %f %f \n","E_dE_d_LJEmin_Acc",i,E_dE_d_LJEmin_Acc[i],E_Total);
		//printf("%s%i %f %f \n","E_dE_d_LJRmin_Acc",i,E_dE_d_LJRmin_Acc[i],E_Total);


		dE_d_LJEmin_Acc[i] += dE_d_LJEmin[i];
		dE_d_LJRmin_Acc[i] += dE_d_LJRmin[i];
                //printf("%s%i %s %f \n","dE_d_LJEmin_Acc",i,myForceField->LJ_Rec[i].RealChem,dE_d_LJEmin_Acc[i]);
                //printf("%s%i %s %f \n","dE_d_LJRmin_Acc",i,myForceField->LJ_Rec[i].RealChem,dE_d_LJRmin_Acc[i]);

		
		E_dE_d_LJEmin_Acc[i] += (dE_d_LJEmin[i] * E_Total);
		E_dE_d_LJRmin_Acc[i] += (dE_d_LJRmin[i] * E_Total);

		//printf("%s%i %f %f \n","E_dE_d_LJEmin_Acc",i,E_dE_d_LJEmin_Acc[i],E_Total);
		//printf("%s%i %f %f \n","E_dE_d_LJRmin_Acc",i,E_dE_d_LJRmin_Acc[i],E_Total);
	}

	//exit(0);

	nAcc_dE_d_LJ++; //number of times this function was called
        //printf("nAcc_dE_d_LJ %d\n",nAcc_dE_d_LJ);

	return;
}

void CMol::Cal_E_VDW_Only(void)
{
	int iActive, iPair, i, j;
	double Rij_SQ, Rij, dxij, dyij, dzij, R6, R12;
	double drdEdr, dedxij, dedyij, dedzij, dE;
	
	E_VDW = 0.0;
	for(iActive=0; iActive<n_NB_Pair; iActive++)	{
		iPair = iActive;	//the index of non-bonded pair
		
		i = NB_List_i[iPair];
		j = NB_List_j[iPair];
		
		if(j >= nAtom)	{
			continue;
		}
		
		dxij=x[j]-x[i];
		dyij=y[j]-y[i];
		dzij=z[j]-z[i];
		
		Rij_SQ=dxij*dxij+dyij*dyij+dzij*dzij;
		
		if(Rij_SQ > Dist_Cutoff_SQ)	{	// two atoms are too far, just skip. Designed for molecule molecule interactions. 
			continue;
		}
		
		//		if(Para_LJ_Epsilon_IJ[iPair] > 1.0E-50)	{	//not zero
		//start	to calculate VDW energy and gradient
		R6=Rij_SQ*Rij_SQ*Rij_SQ;
		R12=R6*R6;
		
		
		dE = Para_LJ_Epsilon_IJ[iPair] * ( Para_LJ_Sigma_IJ_pow_12[iPair]/R12 - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
		
		// R_12 is turned off here !!!
		//			dE = Para_LJ_Epsilon_IJ[iPair] * ( - 2.0*(Para_LJ_Sigma_IJ_pow_6[iPair]/R6) );
		E_VDW += (dE);
		
		//			printf("E_VDW:  %3d %3d  %20.15lf\n", i+1, j+1, dE);
		
		drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( Para_LJ_Sigma_IJ_pow_12[iPair]/(R12*Rij_SQ) - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
		
		// R_12 is turned off here !!!
		//			drdEdr = Para_LJ_Epsilon_IJ[iPair] * 12.0 * ( - Para_LJ_Sigma_IJ_pow_6[iPair]/(R6*Rij_SQ) );
		
		dedxij=drdEdr*dxij;
		dedyij=drdEdr*dyij;
		dedzij=drdEdr*dzij;
		
		
		grad_x[i] += dedxij;
		grad_y[i] += dedyij;
		grad_z[i] += dedzij;
		
		grad_x[j] -= dedxij;
		grad_y[j] -= dedyij;
		grad_z[j] -= dedzij;
		//end	to calculate VDW energy and gradient
		//		}
		
		
		//start	to calculate elec energy and gradient
		Rij = sqrt(Rij_SQ);
		dE = Para_Elec_Pair[iPair] / Rij;
		E_Elec += dE;
		
		//		printf("E_Elec:  %3d  %3d  %25.18lf\n", i+1, j+1, dE);
		
		drdEdr = Para_Elec_Pair[iPair] / (Rij*Rij_SQ);
		
		dedxij=drdEdr*dxij;
		dedyij=drdEdr*dyij;
		dedzij=drdEdr*dzij;
		
		
		grad_x[i] += dedxij;
		grad_y[i] += dedyij;
		grad_z[i] += dedzij;
		
		grad_x[j] -= dedxij;
		grad_y[j] -= dedyij;
		grad_z[j] -= dedzij;
		//end	to calculate elec energy and gradient
	}
}

void CMol::Cal_E_CMap(void)
{
	int iActive;
	double mres,mres1;
	double fx1,fy1,fz1,gx1,gy1,gz1,hx1,hy1,hz1;
	double fx2,fy2,fz2,gx2,gy2,gz2,hx2,hy2,hz2;
	double ax1,ay1,az1,bx1,by1,bz1;
	double ax2,ay2,az2,bx2,by2,bz2;
	double ra21,rb21,ra2r1,rb2r1,rg21,rg1,rgr1;
	double ra22,rb22,ra2r2,rb2r2,rg22,rg2,rgr2;
	double rabr1,cp1,sp1,rabr2,cp2,sp2;
	double phi1,phi2,xphi1,xphi2;
	double ty[4],ty1[4],ty2[4],ty12[4];
	double tc[4][4],tt,tu;
	double e,df1,df2;
	
	int i1,j1,k1,l1,i2,j2,k2,l2,ict;
	int nwarn,nwarnx;
	int iphi1,iphi2;
	double gaa1,gbb1,fg1,hg1,fga1,hgb1;
	double dfx1,dfy1,dfz1,dhx1,dhy1,dhz1,dgx1,dgy1,dgz1;
	double dtfx1,dtfy1,dtfz1;
	double dthx1,dthy1,dthz1,dtgx1,dtgy1,dtgz1;
	double gaa2,gbb2,fg2,hg2,fga2,hgb2;
	double dfx2,dfy2,dfz2,dhx2,dhy2,dhz2,dgx2,dgy2,dgz2;
	double dtfx2,dtfy2,dtfz2;
	double dthx2,dthy2,dthz2,dtgx2,dtgy2,dtgz2;
	int ii;
	double ddf1,ddf2,ddf12,fac;
	const double rxmin=0.005, rxmin2=0.000025;
	
	E_CMap = 0.0;
	
	if(nActive_CMap == 0) return;
	
	nwarn=0;
	nwarnx=0;

	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_CMap; iActive++)	{
			ict = Active_List_CMap[iActive];	//the index of CMap
			
			i1=CMapList[ict][0];
			j1=CMapList[ict][1];
			k1=CMapList[ict][2];
			l1=CMapList[ict][3];
			
			i2=CMapList[ict][4];
			j2=CMapList[ict][5];
			k2=CMapList[ict][6];
			l2=CMapList[ict][7];
			
			fx1=x[i1]-x[j1];
			fy1=y[i1]-y[j1];
			fz1=z[i1]-z[j1];
			gx1=x[j1]-x[k1];
			gy1=y[j1]-y[k1];
			gz1=z[j1]-z[k1];
			hx1=x[l1]-x[k1];
			hy1=y[l1]-y[k1];
			hz1=z[l1]-z[k1];
			
			ax1=fy1*gz1-fz1*gy1;
			ay1=fz1*gx1-fx1*gz1;
			az1=fx1*gy1-fy1*gx1;
			bx1=hy1*gz1-hz1*gy1;
			by1=hz1*gx1-hx1*gz1;
			bz1=hx1*gy1-hy1*gx1;
			
			ra21=ax1*ax1+ay1*ay1+az1*az1;
			rb21=bx1*bx1+by1*by1+bz1*bz1;
			rg21=gx1*gx1+gy1*gy1+gz1*gz1;
			rg1=sqrt(rg21);
			
			if((ra21 <= rxmin2) || (rb21 <= rxmin2) || (rg1 <=rxmin))	{
				nwarn=nwarn+1;
				if( nwarn <= 5 ) {
					//				printf("ECMap: warning.  dihedral %5d is almost linear. Derivatives may be affected for atoms: %5d %5d %5d %5d\n", 
					fprintf(fFile_Run_Log, "ECMap: warning.  dihedral %5d is almost linear. Derivatives may be affected for atoms: %5d %5d %5d %5d\n", 
						ict,i1,j1,k1,l1);
					fflush(fFile_Run_Log);
				}
				continue;
			}
			
			rgr1=1.0/rg1;
			ra2r1=1.0/ra21;
			rb2r1=1.0/rb21;
			rabr1=sqrt(ra2r1*rb2r1);
			cp1=(ax1*bx1+ay1*by1+az1*bz1)*rabr1;
			sp1=rg1*rabr1*(ax1*hx1+ay1*hy1+az1*hz1);
			
			
			if( (cp1 < -0.5) || (cp1 > 0.5) ) {
				phi1=asin(sp1)*radian;
				
				if(cp1 < 0.0) {
					if(phi1 > 0.0) {
						phi1=180.0-phi1;
					}
					else	{
						phi1=-180.0-phi1;
					}
				}
			}
			else	{
				phi1=acos(cp1)*radian;
				if (sp1 < 0.0) {
					phi1=-phi1;
				}
			}
			
			fx2=x[i2]-x[j2];
			fy2=y[i2]-y[j2];
			fz2=z[i2]-z[j2];
			gx2=x[j2]-x[k2];
			gy2=y[j2]-y[k2];
			gz2=z[j2]-z[k2];
			hx2=x[l2]-x[k2];
			hy2=y[l2]-y[k2];
			hz2=z[l2]-z[k2];
			
			ax2=fy2*gz2-fz2*gy2;
			ay2=fz2*gx2-fx2*gz2;
			az2=fx2*gy2-fy2*gx2;
			bx2=hy2*gz2-hz2*gy2;
			by2=hz2*gx2-hx2*gz2;
			bz2=hx2*gy2-hy2*gx2;
			
			ra22=ax2*ax2+ay2*ay2+az2*az2;
			rb22=bx2*bx2+by2*by2+bz2*bz2;
			rg22=gx2*gx2+gy2*gy2+gz2*gz2;
			rg2=sqrt(rg22);
			
			
			if( (ra22<=rxmin2) || (rb22<=rxmin2) || (rg2<=rxmin) )	{
				nwarn=nwarn+1;
				//			printf("Warning in CMap: %d %d %d %d %d\n", ict,i2,j2,k2,l2);
				fprintf(fFile_Run_Log, "Warning in CMap: %d %d %d %d %d\n", ict,i2,j2,k2,l2);
				fflush(fFile_Run_Log);
				continue;
			}
			
			rgr2=1.0/rg2;
			ra2r2=1.0/ra22;
			rb2r2=1.0/rb22;
			rabr2=sqrt(ra2r2*rb2r2);
			cp2=(ax2*bx2+ay2*by2+az2*bz2)*rabr2;
			sp2=rg2*rabr2*(ax2*hx2+ay2*hy2+az2*hz2);
			
			if( (cp2<-0.5) || (cp2>0.5) )	{
				phi2=asin(sp2)*radian;
				if(cp2<0.0)	{
					if(phi2>0.0) {
						phi2=180.0-phi2;
					}
					else	{
						phi2=-180.0-phi2;
					}
				}
			}
			else	{
				phi2=acos(cp2)*radian;
				if(sp2<0.0) {
					phi2=-phi2;
				}
			}
			
			xphi1=phi1+180.0;
			xphi2=phi2+180.0;
			
			if (xphi1 < 0.0) {
				xphi1=xphi1+360.0;
			}
			else if (xphi1>=360.0) {
				xphi1=xphi1-360.0;
			}
			
			if (xphi2<0.0) {
				xphi2=xphi2+360.0;
			}
			else if (xphi2>=360.0) {
				xphi2=xphi2-360.0;
			}
			
			mres=360.0/CMAP_DIM;
			
			//     iphi1=int(xphi1/mres)+1
			//     iphi2=int(xphi2/mres)+1
			
			
			iphi1=int(xphi1/mres);
			iphi2=int(xphi2/mres);
			
			gcstup2(mctp, ty,ty1,ty2,ty12,iphi1,iphi2);
			
			mres1=mres;
			gcscf(ty,ty1,ty2,ty12,mres1,mres,tc);
			
			tt=(xphi1-iphi1*mres)/mres;
			tu=(xphi2-iphi2*mres)/mres;
			
			e=0.0;
			df1=0.0;
			df2=0.0;
			ddf1=0.0;
			ddf2=0.0;
			ddf12=0.0;
			
			for(ii=3; ii>=0; ii--)	{
				e=tt*e+((tc[3][ii]*tu+tc[2][ii])*tu+tc[1][ii])*tu+tc[0][ii];
				df1=tu*df1+(3.0*tc[ii][3]*tt+2.0*tc[ii][2])*tt+tc[ii][1];
				df2=tt*df2+(3.0*tc[3][ii]*tu+2.0*tc[2][ii])*tu+tc[1][ii];
				ddf1=tu*ddf1+6.0*tc[ii][3]*tt+2*tc[ii][2];
				ddf2=tt*ddf2+6.0*tc[3][ii]*tu+2.0*tc[2][ii];
			}
			
			
			ddf12=tc[1][1]+2.0*tc[1][2]*tt+3.0*tc[1][3]*tt*tt+2.0*tu*(tc[2][1]+2.0*tc[2][2]*tt+3.0*tc[2][3]*tt*tt)+3.0*tu*tu*(tc[3][1]+2.0*tc[3][2]*tt+3.0*tc[3][3]*tt*tt);
			fac=radian/mres;
			
			df1=df1*fac;
			df2=df2*fac;
			
			ddf1=ddf1*fac*fac;
			ddf2=ddf2*fac*fac;
			ddf12=ddf12*fac*fac;
			
			E_CMap += e;
			
			fg1=fx1*gx1+fy1*gy1+fz1*gz1;
			hg1=hx1*gx1+hy1*gy1+hz1*gz1;
			fga1=fg1*ra2r1*rgr1;
			hgb1=hg1*rb2r1*rgr1;
			gaa1=-ra2r1*rg1;
			gbb1=rb2r1*rg1;
			
			dtfx1=gaa1*ax1;
			dtfy1=gaa1*ay1;
			dtfz1=gaa1*az1;
			dtgx1=fga1*ax1-hgb1*bx1;
			dtgy1=fga1*ay1-hgb1*by1;
			dtgz1=fga1*az1-hgb1*bz1;
			dthx1=gbb1*bx1;
			dthy1=gbb1*by1;
			dthz1=gbb1*bz1;
			
			dfx1=df1*dtfx1;
			dfy1=df1*dtfy1;
			dfz1=df1*dtfz1;
			dgx1=df1*dtgx1;
			dgy1=df1*dtgy1;
			dgz1=df1*dtgz1;
			dhx1=df1*dthx1;
			dhy1=df1*dthy1;
			dhz1=df1*dthz1;
			
			grad_x[i1] += dfx1;
			grad_y[i1] += dfy1;
			grad_z[i1] += dfz1;
			grad_x[j1] += (-dfx1+dgx1);
			grad_y[j1] += (-dfy1+dgy1);
			grad_z[j1] += (-dfz1+dgz1);
			grad_x[k1] += (-dhx1-dgx1);
			grad_y[k1] += (-dhy1-dgy1);
			grad_z[k1] += (-dhz1-dgz1);
			grad_x[l1] += dhx1;
			grad_y[l1] += dhy1;
			grad_z[l1] += dhz1;
			
			fg2=fx2*gx2+fy2*gy2+fz2*gz2;
			hg2=hx2*gx2+hy2*gy2+hz2*gz2;
			fga2=fg2*ra2r2*rgr2;
			hgb2=hg2*rb2r2*rgr2;
			gaa2=-ra2r2*rg2;
			gbb2=rb2r2*rg2;
			
			dtfx2=gaa2*ax2;
			dtfy2=gaa2*ay2;
			dtfz2=gaa2*az2;
			dtgx2=fga2*ax2-hgb2*bx2;
			dtgy2=fga2*ay2-hgb2*by2;
			dtgz2=fga2*az2-hgb2*bz2;
			dthx2=gbb2*bx2;
			dthy2=gbb2*by2;
			dthz2=gbb2*bz2;
			
			dfx2=df2*dtfx2;
			dfy2=df2*dtfy2;
			dfz2=df2*dtfz2;
			dgx2=df2*dtgx2;
			dgy2=df2*dtgy2;
			dgz2=df2*dtgz2;
			dhx2=df2*dthx2;
			dhy2=df2*dthy2;
			dhz2=df2*dthz2;
			
			grad_x[i2] += dfx2;
			grad_y[i2] += dfy2;
			grad_z[i2] += dfz2;
			grad_x[j2] += (-dfx2+dgx2);
			grad_y[j2] += (-dfy2+dgy2);
			grad_z[j2] += (-dfz2+dgz2);
			grad_x[k2] += (-dhx2-dgx2);
			grad_y[k2] += (-dhy2-dgy2);
			grad_z[k2] += (-dhz2-dgz2);
			grad_x[l2] += dhx2;
			grad_y[l2] += dhy2;
			grad_z[l2] += dhz2;
		}
	}
	else	{
		for(iActive=0; iActive<nCMapTerm; iActive++)	{
			ict = iActive;	//the index of CMap
			
			i1=CMapList[ict][0];
			j1=CMapList[ict][1];
			k1=CMapList[ict][2];
			l1=CMapList[ict][3];
			
			i2=CMapList[ict][4];
			j2=CMapList[ict][5];
			k2=CMapList[ict][6];
			l2=CMapList[ict][7];
			
			fx1=x[i1]-x[j1];
			fy1=y[i1]-y[j1];
			fz1=z[i1]-z[j1];
			gx1=x[j1]-x[k1];
			gy1=y[j1]-y[k1];
			gz1=z[j1]-z[k1];
			hx1=x[l1]-x[k1];
			hy1=y[l1]-y[k1];
			hz1=z[l1]-z[k1];
			
			ax1=fy1*gz1-fz1*gy1;
			ay1=fz1*gx1-fx1*gz1;
			az1=fx1*gy1-fy1*gx1;
			bx1=hy1*gz1-hz1*gy1;
			by1=hz1*gx1-hx1*gz1;
			bz1=hx1*gy1-hy1*gx1;
			
			ra21=ax1*ax1+ay1*ay1+az1*az1;
			rb21=bx1*bx1+by1*by1+bz1*bz1;
			rg21=gx1*gx1+gy1*gy1+gz1*gz1;
			rg1=sqrt(rg21);
			
			if((ra21 <= rxmin2) || (rb21 <= rxmin2) || (rg1 <=rxmin))	{
				nwarn=nwarn+1;
				if( nwarn <= 5 ) {
					//				printf("ECMap: warning.  dihedral %5d is almost linear. Derivatives may be affected for atoms: %5d %5d %5d %5d\n", 
					fprintf(fFile_Run_Log, "ECMap: warning.  dihedral %5d is almost linear. Derivatives may be affected for atoms: %5d %5d %5d %5d\n", 
						ict,i1,j1,k1,l1);
					fflush(fFile_Run_Log);
				}
				continue;
			}
			
			rgr1=1.0/rg1;
			ra2r1=1.0/ra21;
			rb2r1=1.0/rb21;
			rabr1=sqrt(ra2r1*rb2r1);
			cp1=(ax1*bx1+ay1*by1+az1*bz1)*rabr1;
			sp1=rg1*rabr1*(ax1*hx1+ay1*hy1+az1*hz1);
			
			
			if( (cp1 < -0.5) || (cp1 > 0.5) ) {
				phi1=asin(sp1)*radian;
				
				if(cp1 < 0.0) {
					if(phi1 > 0.0) {
						phi1=180.0-phi1;
					}
					else	{
						phi1=-180.0-phi1;
					}
				}
			}
			else	{
				phi1=acos(cp1)*radian;
				if (sp1 < 0.0) {
					phi1=-phi1;
				}
			}
			
			fx2=x[i2]-x[j2];
			fy2=y[i2]-y[j2];
			fz2=z[i2]-z[j2];
			gx2=x[j2]-x[k2];
			gy2=y[j2]-y[k2];
			gz2=z[j2]-z[k2];
			hx2=x[l2]-x[k2];
			hy2=y[l2]-y[k2];
			hz2=z[l2]-z[k2];
			
			ax2=fy2*gz2-fz2*gy2;
			ay2=fz2*gx2-fx2*gz2;
			az2=fx2*gy2-fy2*gx2;
			bx2=hy2*gz2-hz2*gy2;
			by2=hz2*gx2-hx2*gz2;
			bz2=hx2*gy2-hy2*gx2;
			
			ra22=ax2*ax2+ay2*ay2+az2*az2;
			rb22=bx2*bx2+by2*by2+bz2*bz2;
			rg22=gx2*gx2+gy2*gy2+gz2*gz2;
			rg2=sqrt(rg22);
			
			
			if( (ra22<=rxmin2) || (rb22<=rxmin2) || (rg2<=rxmin) )	{
				nwarn=nwarn+1;
				//			printf("Warning in CMap: %d %d %d %d %d\n", ict,i2,j2,k2,l2);
				fprintf(fFile_Run_Log, "Warning in CMap: %d %d %d %d %d\n", ict,i2,j2,k2,l2);
				fflush(fFile_Run_Log);
				continue;
			}
			
			rgr2=1.0/rg2;
			ra2r2=1.0/ra22;
			rb2r2=1.0/rb22;
			rabr2=sqrt(ra2r2*rb2r2);
			cp2=(ax2*bx2+ay2*by2+az2*bz2)*rabr2;
			sp2=rg2*rabr2*(ax2*hx2+ay2*hy2+az2*hz2);
			
			if( (cp2<-0.5) || (cp2>0.5) )	{
				phi2=asin(sp2)*radian;
				if(cp2<0.0)	{
					if(phi2>0.0) {
						phi2=180.0-phi2;
					}
					else	{
						phi2=-180.0-phi2;
					}
				}
			}
			else	{
				phi2=acos(cp2)*radian;
				if(sp2<0.0) {
					phi2=-phi2;
				}
			}
			
			xphi1=phi1+180.0;
			xphi2=phi2+180.0;
			
			if (xphi1 < 0.0) {
				xphi1=xphi1+360.0;
			}
			else if (xphi1>=360.0) {
				xphi1=xphi1-360.0;
			}
			
			if (xphi2<0.0) {
				xphi2=xphi2+360.0;
			}
			else if (xphi2>=360.0) {
				xphi2=xphi2-360.0;
			}
			
			mres=360.0/CMAP_DIM;
			
			//     iphi1=int(xphi1/mres)+1
			//     iphi2=int(xphi2/mres)+1
			
			
			iphi1=int(xphi1/mres);
			iphi2=int(xphi2/mres);
			
			gcstup2(mctp, ty,ty1,ty2,ty12,iphi1,iphi2);
			
			mres1=mres;
			gcscf(ty,ty1,ty2,ty12,mres1,mres,tc);
			
			tt=(xphi1-iphi1*mres)/mres;
			tu=(xphi2-iphi2*mres)/mres;
			
			e=0.0;
			df1=0.0;
			df2=0.0;
			ddf1=0.0;
			ddf2=0.0;
			ddf12=0.0;
			
			for(ii=3; ii>=0; ii--)	{
				e=tt*e+((tc[3][ii]*tu+tc[2][ii])*tu+tc[1][ii])*tu+tc[0][ii];
				df1=tu*df1+(3.0*tc[ii][3]*tt+2.0*tc[ii][2])*tt+tc[ii][1];
				df2=tt*df2+(3.0*tc[3][ii]*tu+2.0*tc[2][ii])*tu+tc[1][ii];
				ddf1=tu*ddf1+6.0*tc[ii][3]*tt+2*tc[ii][2];
				ddf2=tt*ddf2+6.0*tc[3][ii]*tu+2.0*tc[2][ii];
			}
			
			
			ddf12=tc[1][1]+2.0*tc[1][2]*tt+3.0*tc[1][3]*tt*tt+2.0*tu*(tc[2][1]+2.0*tc[2][2]*tt+3.0*tc[2][3]*tt*tt)+3.0*tu*tu*(tc[3][1]+2.0*tc[3][2]*tt+3.0*tc[3][3]*tt*tt);
			fac=radian/mres;
			
			df1=df1*fac;
			df2=df2*fac;
			
			ddf1=ddf1*fac*fac;
			ddf2=ddf2*fac*fac;
			ddf12=ddf12*fac*fac;
			
			E_CMap += e;
			
			fg1=fx1*gx1+fy1*gy1+fz1*gz1;
			hg1=hx1*gx1+hy1*gy1+hz1*gz1;
			fga1=fg1*ra2r1*rgr1;
			hgb1=hg1*rb2r1*rgr1;
			gaa1=-ra2r1*rg1;
			gbb1=rb2r1*rg1;
			
			dtfx1=gaa1*ax1;
			dtfy1=gaa1*ay1;
			dtfz1=gaa1*az1;
			dtgx1=fga1*ax1-hgb1*bx1;
			dtgy1=fga1*ay1-hgb1*by1;
			dtgz1=fga1*az1-hgb1*bz1;
			dthx1=gbb1*bx1;
			dthy1=gbb1*by1;
			dthz1=gbb1*bz1;
			
			dfx1=df1*dtfx1;
			dfy1=df1*dtfy1;
			dfz1=df1*dtfz1;
			dgx1=df1*dtgx1;
			dgy1=df1*dtgy1;
			dgz1=df1*dtgz1;
			dhx1=df1*dthx1;
			dhy1=df1*dthy1;
			dhz1=df1*dthz1;
			
			grad_x[i1] += dfx1;
			grad_y[i1] += dfy1;
			grad_z[i1] += dfz1;
			grad_x[j1] += (-dfx1+dgx1);
			grad_y[j1] += (-dfy1+dgy1);
			grad_z[j1] += (-dfz1+dgz1);
			grad_x[k1] += (-dhx1-dgx1);
			grad_y[k1] += (-dhy1-dgy1);
			grad_z[k1] += (-dhz1-dgz1);
			grad_x[l1] += dhx1;
			grad_y[l1] += dhy1;
			grad_z[l1] += dhz1;
			
			fg2=fx2*gx2+fy2*gy2+fz2*gz2;
			hg2=hx2*gx2+hy2*gy2+hz2*gz2;
			fga2=fg2*ra2r2*rgr2;
			hgb2=hg2*rb2r2*rgr2;
			gaa2=-ra2r2*rg2;
			gbb2=rb2r2*rg2;
			
			dtfx2=gaa2*ax2;
			dtfy2=gaa2*ay2;
			dtfz2=gaa2*az2;
			dtgx2=fga2*ax2-hgb2*bx2;
			dtgy2=fga2*ay2-hgb2*by2;
			dtgz2=fga2*az2-hgb2*bz2;
			dthx2=gbb2*bx2;
			dthy2=gbb2*by2;
			dthz2=gbb2*bz2;
			
			dfx2=df2*dtfx2;
			dfy2=df2*dtfy2;
			dfz2=df2*dtfz2;
			dgx2=df2*dtgx2;
			dgy2=df2*dtgy2;
			dgz2=df2*dtgz2;
			dhx2=df2*dthx2;
			dhy2=df2*dthy2;
			dhz2=df2*dthz2;
			
			grad_x[i2] += dfx2;
			grad_y[i2] += dfy2;
			grad_z[i2] += dfz2;
			grad_x[j2] += (-dfx2+dgx2);
			grad_y[j2] += (-dfy2+dgy2);
			grad_z[j2] += (-dfz2+dgz2);
			grad_x[k2] += (-dhx2-dgx2);
			grad_y[k2] += (-dhy2-dgy2);
			grad_z[k2] += (-dhz2-dgz2);
			grad_x[l2] += dhx2;
			grad_y[l2] += dhy2;
			grad_z[l2] += dhz2;
		}
	}
	
	
	nwarn=nwarn+nwarnx;
	
	if(nwarn > 5)	{
//		printf("Total of %5d warnings from ecmap.\n", nwarn);
		fprintf(fFile_Run_Log, "Total of %5d warnings from ecmap.\n", nwarn);
		fflush(fFile_Run_Log);
	}
	
	return;
}

#define KHYPER	(40000.0)
#define HORDER	(4.0)
#define RHYPER	(0.2)
void CMol::Cal_E_DrudeHyper(void)
{
	double dx, dy, dz, rij, rij_SQ, dudr, RHYPER_SQ;
	int i, j;

	E_DrudeHyper = 0.0;

	RHYPER_SQ = RHYPER*RHYPER;
	for(i=1; i<nAtom; i++)	{
		if(IsDrude[i])	{
			j = i - 1;
			dx=x[i]-x[j];
			dy=y[i]-y[j];
			dz=z[i]-z[j];
			rij_SQ = dx*dx + dy*dy + dz*dz;
			
			if(rij_SQ > RHYPER_SQ)	{
				rij = sqrt(rij_SQ);
				dudr = pow( KHYPER*HORDER*(rij-RHYPER),  HORDER-2.0);
				E_DrudeHyper += ( KHYPER*pow(rij-RHYPER, HORDER) );
				
				grad_x[i] += dudr*dx;
				grad_y[i] += dudr*dy;
				grad_z[i] += dudr*dz;
				
				grad_x[j] -= dudr*dx;
				grad_y[j] -= dudr*dy;
				grad_z[j] -= dudr*dz;
			}
		}
	}
}

void CMol::Cal_E_Anisotropy(void)
{
	int iActive, Idx;
	int i, j, l, m, n;
	double u1x, u1y, u1z, u2x, u2y, u2z, drx, dry, drz, dpar, dperp;
	double u1mag, u2mag, kpar0, kperp0, kiso0;
	
	E_Aniso = 0.0;
	
	if(Geo_Opt_Drude_Only)	{
		for(iActive=0; iActive<nActive_Aniso; iActive++)	{
			Idx = Active_List_Aniso[iActive];	//the index of anisotropy
			
			i = AnisoList[Idx][0];
			j = i + 1;
			l = AnisoList[Idx][1];
			m = AnisoList[Idx][2];
			n = AnisoList[Idx][3];
			
			u1x = x[i] - x[l];
			u1y = y[i] - y[l];
			u1z = z[i] - z[l];
			u2x = x[m] - x[n];
			u2y = y[m] - y[n];
			u2z = z[m] - z[n];
			
			u1mag = 1.0/sqrt( u1x*u1x + u1y*u1y + u1z*u1z );
			u2mag = 1.0/sqrt( u2x*u2x + u2y*u2y + u2z*u2z );
			
			u1x *= u1mag;
			u1y *= u1mag;
			u1z *= u1mag;
			u2x *= u2mag;
			u2y *= u2mag;
			u2z *= u2mag;
			
			drx = x[j] - x[i];
			dry = y[j] - y[i];
			drz = z[j] - z[i];
			
			kpar0  = 2.0*Para_Aniso[Idx][0];
			kperp0 = 2.0*Para_Aniso[Idx][1];
			kiso0  = 2.0*Para_Aniso[Idx][2];
			
			dpar = drx*u1x + dry*u1y + drz*u1z;
			dperp = drx*u2x + dry*u2y + drz*u2z;
			
			E_Aniso += ( 0.5*kpar0*dpar*dpar + 0.5*kperp0*dperp*dperp + 0.5*kiso0*(drx*drx+dry*dry+drz*drz) );
			
			grad_x[i] -= kiso0*drx;
			grad_y[i] -= kiso0*dry;
			grad_z[i] -= kiso0*drz;
			
			grad_x[j] += kiso0*drx;
			grad_y[j] += kiso0*dry;
			grad_z[j] += kiso0*drz;
			
			grad_x[i]=grad_x[i] + kpar0*dpar*(-u1x + (drx - u1x*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2x;
			grad_y[i]=grad_y[i] + kpar0*dpar*(-u1y + (dry - u1y*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2y;
			grad_z[i]=grad_z[i] + kpar0*dpar*(-u1z + (drz - u1z*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2z;
			
			grad_x[j] =  grad_x[j] + kpar0*dpar*u1x + kperp0*dperp*u2x;
			grad_y[j] =  grad_y[j] + kpar0*dpar*u1y + kperp0*dperp*u2y;
			grad_z[j] =  grad_z[j] + kpar0*dpar*u1z + kperp0*dperp*u2z;
			
			grad_x[l] =  grad_x[l] + kpar0*dpar*(-drx + u1x*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			grad_y[l] =  grad_y[l] + kpar0*dpar*(-dry + u1y*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			grad_z[l] =  grad_z[l] + kpar0*dpar*(-drz + u1z*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			
			grad_x[m] = grad_x[m]+kperp0*dperp*(drx - u2x*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_y[m] = grad_y[m]+kperp0*dperp*(dry - u2y*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_z[m] = grad_z[m]+kperp0*dperp*(drz - u2z*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			
			grad_x[n] = grad_x[n]+kperp0*dperp*(-drx + u2x*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_y[n] = grad_y[n]+kperp0*dperp*(-dry + u2y*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_z[n] = grad_z[n]+kperp0*dperp*(-drz + u2z*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
		}
	}
	else	{
		for(iActive=0; iActive<nAniso; iActive++)	{
			Idx = iActive;	//the index of anisotropy
			
			i = AnisoList[Idx][0];
			j = i + 1;
			l = AnisoList[Idx][1];
			m = AnisoList[Idx][2];
			n = AnisoList[Idx][3];
			
			u1x = x[i] - x[l];
			u1y = y[i] - y[l];
			u1z = z[i] - z[l];
			u2x = x[m] - x[n];
			u2y = y[m] - y[n];
			u2z = z[m] - z[n];
			
			u1mag = 1.0/sqrt( u1x*u1x + u1y*u1y + u1z*u1z );
			u2mag = 1.0/sqrt( u2x*u2x + u2y*u2y + u2z*u2z );
			
			u1x *= u1mag;
			u1y *= u1mag;
			u1z *= u1mag;
			u2x *= u2mag;
			u2y *= u2mag;
			u2z *= u2mag;
			
			drx = x[j] - x[i];
			dry = y[j] - y[i];
			drz = z[j] - z[i];
			
			kpar0  = 2.0*Para_Aniso[Idx][0];
			kperp0 = 2.0*Para_Aniso[Idx][1];
			kiso0  = 2.0*Para_Aniso[Idx][2];
			
			dpar = drx*u1x + dry*u1y + drz*u1z;
			dperp = drx*u2x + dry*u2y + drz*u2z;
			
			E_Aniso += ( 0.5*kpar0*dpar*dpar + 0.5*kperp0*dperp*dperp + 0.5*kiso0*(drx*drx+dry*dry+drz*drz) );
			
			grad_x[i] -= kiso0*drx;
			grad_y[i] -= kiso0*dry;
			grad_z[i] -= kiso0*drz;
			
			grad_x[j] += kiso0*drx;
			grad_y[j] += kiso0*dry;
			grad_z[j] += kiso0*drz;
			
			grad_x[i]=grad_x[i] + kpar0*dpar*(-u1x + (drx - u1x*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2x;
			grad_y[i]=grad_y[i] + kpar0*dpar*(-u1y + (dry - u1y*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2y;
			grad_z[i]=grad_z[i] + kpar0*dpar*(-u1z + (drz - u1z*(drx*u1x+dry*u1y+drz*u1z))*u1mag) - kperp0*dperp*u2z;
			
			grad_x[j] =  grad_x[j] + kpar0*dpar*u1x + kperp0*dperp*u2x;
			grad_y[j] =  grad_y[j] + kpar0*dpar*u1y + kperp0*dperp*u2y;
			grad_z[j] =  grad_z[j] + kpar0*dpar*u1z + kperp0*dperp*u2z;
			
			grad_x[l] =  grad_x[l] + kpar0*dpar*(-drx + u1x*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			grad_y[l] =  grad_y[l] + kpar0*dpar*(-dry + u1y*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			grad_z[l] =  grad_z[l] + kpar0*dpar*(-drz + u1z*(drx*u1x+dry*u1y+drz*u1z))*u1mag;
			
			grad_x[m] = grad_x[m]+kperp0*dperp*(drx - u2x*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_y[m] = grad_y[m]+kperp0*dperp*(dry - u2y*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_z[m] = grad_z[m]+kperp0*dperp*(drz - u2z*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			
			grad_x[n] = grad_x[n]+kperp0*dperp*(-drx + u2x*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_y[n] = grad_y[n]+kperp0*dperp*(-dry + u2y*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
			grad_z[n] = grad_z[n]+kperp0*dperp*(-drz + u2z*(drx*u2x+dry*u2y+drz*u2z))*u2mag;
		}
	}


	return;	
}

double Fun_E_Thole(double qq, double aa, double rij, double& df)
{
	double norm, au, expau, polyau, e_thole_1;
	
	if(aa < 0.0)	{	// only used for test
//		printf("Error in Fun_E_Thole(), aa < 0.0 !\nQuit.\n");
		fprintf(fFile_Run_Log, "Error in Fun_E_Thole(), aa < 0.0 !\nQuit.\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	norm = pow(aa, -1.0/6.0);
	au = norm * rij;
	expau = exp(-au);
	polyau = 1.0 + au*0.5;
	e_thole_1 = (qq/rij) * (1.0-polyau*expau);
	polyau = 1.0 + au*(1.0+0.5*au);
	df = qq/(rij*rij*rij) * (polyau*expau-1.0);

	return e_thole_1;
}

void CMol::Cal_E_Thole(void)
{
	int Idx, i, j;
	double dfaa, dfda, dfad, dfdd;
	double aa, qq;
	double rx,ry,rz,rij;
	double dxi, dyi, dzi;

	E_Thole = 0.0;

	for(Idx=0; Idx<n_Thole_Pair; Idx++)	{
		i = Thole_List_i[Idx];
		j = Thole_List_j[Idx];
		
//		printf(" %3d %3d\n", i+1, j+1);
		
		aa = Para_Thole_aa[Idx];
		qq = Para_Thole_qq[Idx];
		
		
		rx = x[i]-x[j];		//atom-atom
		ry = y[i]-y[j];
		rz = z[i]-z[j];
		rij = sqrt(rx*rx + ry*ry + rz*rz);
		E_Thole += Fun_E_Thole(qq,aa,rij, dfaa);
		
		dxi = rx*dfaa;
		dyi = ry*dfaa;
		dzi = rz*dfaa;
		grad_x[i] += dxi;
		grad_y[i] += dyi;
		grad_z[i] += dzi;
		grad_x[j] -= dxi;
		grad_y[j] -= dyi;
		grad_z[j] -= dzi;
		
		
		rx = x[i]-x[j+1];	//atom-drude
		ry = y[i]-y[j+1];
		rz = z[i]-z[j+1];
		rij = sqrt(rx*rx + ry*ry + rz*rz);
		E_Thole += Fun_E_Thole(-qq,aa,rij, dfad);
		
		dxi = rx*dfad;
		dyi = ry*dfad;
		dzi = rz*dfad;
		
		grad_x[i  ] += dxi;
		grad_y[i  ] += dyi;
		grad_z[i  ] += dzi;
		grad_x[j+1] -= dxi;
		grad_y[j+1] -= dyi;
		grad_z[j+1] -= dzi;
		
		
		rx = x[i+1]-x[j];	//drude-atom
		ry = y[i+1]-y[j];
		rz = z[i+1]-z[j];
		rij = sqrt(rx*rx + ry*ry + rz*rz);
		E_Thole += Fun_E_Thole(-qq,aa,rij, dfda);
		
		dxi = rx*dfda;
		dyi = ry*dfda;
		dzi = rz*dfda;
		
		grad_x[i+1] += dxi;
		grad_y[i+1] += dyi;
		grad_z[i+1] += dzi;
		grad_x[j  ] -= dxi;
		grad_y[j  ] -= dyi;
		grad_z[j  ] -= dzi;
		
		
		rx = x[i+1]-x[j+1];	//drude-drude
		ry = y[i+1]-y[j+1];
		rz = z[i+1]-z[j+1];
		rij = sqrt(rx*rx + ry*ry + rz*rz);
		E_Thole += Fun_E_Thole(qq,aa,rij, dfdd);
		
		dxi = rx*dfdd;
		dyi = ry*dfdd;
		dzi = rz*dfdd;
		
		grad_x[i+1] += dxi;
		grad_y[i+1] += dyi;
		grad_z[i+1] += dzi;
		grad_x[j+1] -= dxi;
		grad_y[j+1] -= dyi;
		grad_z[j+1] -= dzi;
		
	}

}



void CMol::Position_LonePair(void)
{
	int Idx, i, j, k, l;
	double Dist;

	for(Idx=0; Idx<nLP; Idx++)	{
		i = LPList[Idx];	//the atom of lone pair

		j = LP_PosHost[Idx][0];
		k = LP_PosHost[Idx][1];
		l = LP_PosHost[Idx][2];

		Dist = LP_Dist[Idx];

		if(Dist == 0.0)	{	// zero. center of mass/geometry 
			x[i] = (x[j] + x[k] + x[l]) / 3.0;
			y[i] = (y[j] + y[k] + y[l]) / 3.0;
			z[i] = (z[j] + z[k] + z[l]) / 3.0;
		}
		else if(Dist < 0.0)	{
			Dist = -Dist;
			x[i] = (x[k] + x[l]) * 0.5;
			y[i] = (y[k] + y[l]) * 0.5;
			z[i] = (z[k] + z[l]) * 0.5;

			CartCV_Charmm(l,i,j,i,Idx, Dist);
		}
		else	{
			CartCV_Charmm(l,k,j,i,Idx, Dist);
		}

	}

	return;
}

void CMol::CartCV_Charmm(int i, int j, int k, int l, int Idx_LP, double dist)	//translated from charmm
{
	double cst, snt, csp, snp, xa, ya, za, xb, yb, zb, xc, yc, zc, ra, rc, wa, wb, wc;
	
//	dist = LP_Dist[Idx_LP];
	cst = LP_Cos_Theta[Idx_LP];
	snt = LP_Sin_Theta[Idx_LP];
	csp = LP_Cos_Phi[Idx_LP];
	snp = LP_Sin_Phi[Idx_LP];

    xa=x[j]-x[k];
    ya=y[j]-y[k];
    za=z[j]-z[k];
    xb=x[i]-x[j];
    yb=y[i]-y[j];
    zb=z[i]-z[j];

	ra = 1.0/sqrt(xa*xa + ya*ya + za*za);
    xa*=ra;
    ya*=ra;
    za*=ra;

    xc=yb*za-zb*ya;
    yc=zb*xa-xb*za;
    zc=xb*ya-yb*xa;

	rc = 1.0/sqrt(xc*xc + yc*yc + zc*zc);

	if(rc < 1.0E-9)	{
//		printf(" Note from CARTCV: I-J-K is linear.\nQuit\n");
		fprintf(fFile_Run_Log, " Note from CARTCV: I-J-K is linear.\nQuit\n");
		fflush(fFile_Run_Log);
		exit(1);
	}

	xc *= rc;
	yc *= rc;
	zc *= rc;

    xb=ya*zc-za*yc;
    yb=za*xc-xa*zc;
    zb=xa*yc-ya*xc;

    wa=dist*cst;
    wb=dist*snt*csp;
    wc=dist*snt*snp;

    x[l]=x[k]+wa*xa+wb*xb+wc*xc;
    y[l]=y[k]+wa*ya+wb*yb+wc*yc;
    z[l]=z[k]+wa*za+wb*zb+wc*zc;

	return;
}

void CMol::LonePairForceReDistribute(void)
{
	int Idx, QBis=0;
	int i,j,k,l;
	double sdot_SQ;
    double vx,vy,vz,qx,qy,qz,px,py,pz;
    double rrdot,srdot,stdot,ttdot,fdot,pdot,qdot,rdot,sdot,tdot,fpdot;
    double fpx,fpy,fpz,fqx,fqy,fqz,frx,fry,frz;
    double ftx,fty,ftz,rrx,rry,rrz,ttx,tty,ttz;
    double rx,ry,rz,sx,sy,sz,tx,ty,tz,fact;
	
	for(Idx=0; Idx<nLP; Idx++)	{
		i = LPList[Idx];	//the atom of lone pair
		
		j = LP_PosHost[Idx][0];
		k = LP_PosHost[Idx][1];
		l = LP_PosHost[Idx][2];
		
		QBis = 0;	// default: relative
		if(LP_Dist[Idx] < 0)	{	// bisector
			QBis = 1;
		}
		//		Dist = LP_Dist[Idx];
		//		Theta = LP_Theta[Idx];
		//		Phi = LP_Phi[Idx];
		
		rx=x[i]-x[j];
		ry=y[i]-y[j];
		rz=z[i]-z[j];
		rdot=rx*rx+ry*ry+rz*rz;
		fdot=(grad_x[i]*rx+grad_y[i]*ry+grad_z[i]*rz)/rdot;
		rdot=sqrt(rdot);
		frx=rx*fdot;
		fry=ry*fdot;
		frz=rz*fdot;
		
		grad_x[j] += frx;
		grad_y[j] += fry;
		grad_z[j] += frz;
		
		if(QBis)	{
			sx=x[j]-(x[k]+x[l])*0.5;
			sy=y[j]-(y[k]+y[l])*0.5;
			sz=z[j]-(z[k]+z[l])*0.5;
			tx=(x[k]-x[l])*0.5;
			ty=(y[k]-y[l])*0.5;
			tz=(z[k]-z[l])*0.5;
		}
		else	{
			sx=x[j]-x[k];
			sy=y[j]-y[k];
			sz=z[j]-z[k];
			tx=x[k]-x[l];
			ty=y[k]-y[l];
			tz=z[k]-z[l];
		}
		
		sdot=sqrt(sx*sx+sy*sy+sz*sz);
		sdot_SQ = sdot*sdot;
		tdot=sqrt(tx*tx+ty*ty+tz*tz);
		px=ry*sz-rz*sy;
		py=rz*sx-rx*sz;
		pz=rx*sy-ry*sx;
		pdot=sqrt(px*px+py*py+pz*pz);
		qx=sy*tz-sz*ty;
		qy=sz*tx-sx*tz;
		qz=sx*ty-sy*tx;
		qdot=sqrt(qx*qx+qy*qy+qz*qz);
		
		if(pdot < 1.0E-10)	{	//very small
			ftx=grad_x[i]-frx;
			fty=grad_y[i]-fry;
			ftz=grad_z[i]-frz;
			fact=(sdot-rdot)/sdot;
			grad_x[j] += (ftx*fact);
			grad_y[j] += (fty*fact);
			grad_z[j] += (ftz*fact);
			fact=rdot/sdot;
			
			if(QBis)	{
				fact *= 0.5;
				grad_x[k] += (ftx*fact);
				grad_y[k] += (fty*fact);
				grad_z[k] += (ftz*fact);
				grad_x[l] += (ftx*fact);
				grad_y[l] += (fty*fact);
				grad_z[l] += (ftz*fact);
			}
			else	{
				grad_x[k] += (ftx*fact);
				grad_y[k] += (fty*fact);
				grad_z[k] += (ftz*fact);
			}
		}
		else	{
			fpdot=(grad_x[i]*px+grad_y[i]*py+grad_z[i]*pz)/pdot;
			fpx=px*fpdot;
			fpy=py*fpdot;
			fpz=pz*fpdot;
			ftx=grad_x[i]-frx-fpx;
			fty=grad_y[i]-fry-fpy;
			ftz=grad_z[i]-frz-fpz;
			
			grad_x[j] += ftx;
			grad_y[j] += fty;
			grad_z[j] += ftz;
			
			vx=ry*ftz-rz*fty;
			vy=rz*ftx-rx*ftz;
			vz=rx*fty-ry*ftx;
			fact=1.0/sdot_SQ;
			ftx=(sy*vz-sz*vy)*fact;
			fty=(sz*vx-sx*vz)*fact;
			ftz=(sx*vy-sy*vx)*fact;
			
			grad_x[j] -= ftx;
			grad_y[j] -= fty;
			grad_z[j] -= ftz;
			
			if(QBis)	{
				grad_x[k] += (ftx*0.5);
				grad_y[k] += (fty*0.5);
				grad_z[k] += (ftz*0.5);
				grad_x[l] += (ftx*0.5);
				grad_y[l] += (fty*0.5);
				grad_z[l] += (ftz*0.5);
			}
			else	{
				grad_x[k] += ftx;
				grad_y[k] += fty;
				grad_z[k] += ftz;
			}
			
			srdot=(sx*rx+sy*ry+sz*rz)/sdot_SQ;
			rrx=rx-sx*srdot;
			rry=ry-sy*srdot;
			rrz=rz-sz*srdot;
			rrdot=sqrt(rrx*rrx+rry*rry+rrz*rrz);
			stdot=(sx*tx+sy*ty+sz*tz)/sdot_SQ;
			ttx=tx-sx*stdot;
			tty=ty-sy*stdot;
			ttz=tz-sz*stdot;
			ttdot=sqrt(ttx*ttx+tty*tty+ttz*ttz);
			fact=(rrdot*fpdot)/(ttdot*qdot);
			fqx=qx*fact;
			fqy=qy*fact;
			fqz=qz*fact;
			
			grad_x[l] += fqx;
			grad_y[l] += fqy;
			grad_z[l] += fqz;
			
			grad_x[j] += (fpx*(1.0+srdot)+fqx*stdot);
			grad_y[j] += (fpy*(1.0+srdot)+fqy*stdot);
			grad_z[j] += (fpz*(1.0+srdot)+fqz*stdot);
			
			ftx=fqx*(1.0+stdot)+fpx*srdot;
			fty=fqy*(1.0+stdot)+fpy*srdot;
			ftz=fqz*(1.0+stdot)+fpz*srdot;
			
			if(QBis)	{
				grad_x[k] -= (ftx*0.5);
				grad_y[k] -= (fty*0.5);
				grad_z[k] -= (ftz*0.5);
				grad_x[l] -= (ftx*0.5);
				grad_y[l] -= (fty*0.5);
				grad_z[l] -= (ftz*0.5);
			}
			else	{
				grad_x[k] -= ftx;
				grad_y[k] -= fty;
				grad_z[k] -= ftz;
			}
		}
		
		
		grad_x[i]=0.0;
		grad_y[i]=0.0;
		grad_z[i]=0.0;
	}
	
	return;	
}

/*
void CMol::LonePairForceReDistribute(void)
{
	int Idx;
	int i,j,k,l;
	double sdot_SQ;
    double vx,vy,vz,qx,qy,qz,px,py,pz;
    double rrdot,srdot,stdot,ttdot,fdot,pdot,qdot,rdot,sdot,tdot,fpdot;
    double fpx,fpy,fpz,fqx,fqy,fqz,frx,fry,frz;
    double ftx,fty,ftz,rrx,rry,rrz,ttx,tty,ttz;
    double rx,ry,rz,sx,sy,sz,tx,ty,tz,fact;

	for(Idx=0; Idx<nLP; Idx++)	{
		i = LPList[Idx];	//the atom of lone pair

		j = LP_PosHost[Idx][0];
		k = LP_PosHost[Idx][1];
		l = LP_PosHost[Idx][2];

//		Dist = LP_Dist[Idx];
//		Theta = LP_Theta[Idx];
//		Phi = LP_Phi[Idx];
		
		rx=x[i]-x[j];
		ry=y[i]-y[j];
		rz=z[i]-z[j];
		rdot=rx*rx+ry*ry+rz*rz;
		fdot=(grad_x[i]*rx+grad_y[i]*ry+grad_z[i]*rz)/rdot;
		rdot=sqrt(rdot);
		frx=rx*fdot;
		fry=ry*fdot;
		frz=rz*fdot;
		
		grad_x[j] += frx;
		grad_y[j] += fry;
		grad_z[j] += frz;
		
		sx=x[j]-x[k];
		sy=y[j]-y[k];
		sz=z[j]-z[k];
		tx=x[k]-x[l];
		ty=y[k]-y[l];
		tz=z[k]-z[l];
		
		sdot=sqrt(sx*sx+sy*sy+sz*sz);
		sdot_SQ = sdot*sdot;
		tdot=sqrt(tx*tx+ty*ty+tz*tz);
		px=ry*sz-rz*sy;
		py=rz*sx-rx*sz;
		pz=rx*sy-ry*sx;
		pdot=sqrt(px*px+py*py+pz*pz);
		qx=sy*tz-sz*ty;
		qy=sz*tx-sx*tz;
		qz=sx*ty-sy*tx;
		qdot=sqrt(qx*qx+qy*qy+qz*qz);
		
		if(pdot < 1.0E-10)	{	//very small
			ftx=grad_x[i]-frx;
			fty=grad_y[i]-fry;
			ftz=grad_z[i]-frz;
			fact=(sdot-rdot)/sdot;
			grad_x[j] += (ftx*fact);
			grad_y[j] += (fty*fact);
			grad_z[j] += (ftz*fact);
			fact=rdot/sdot;
			grad_x[k] += (ftx*fact);
			grad_y[k] += (fty*fact);
			grad_z[k] += (ftz*fact);
		}
		else	{
                fpdot=(grad_x[i]*px+grad_y[i]*py+grad_z[i]*pz)/pdot;
                fpx=px*fpdot;
                fpy=py*fpdot;
                fpz=pz*fpdot;
                ftx=grad_x[i]-frx-fpx;
                fty=grad_y[i]-fry-fpy;
                ftz=grad_z[i]-frz-fpz;

                grad_x[j] += ftx;
                grad_y[j] += fty;
                grad_z[j] += ftz;

                vx=ry*ftz-rz*fty;
                vy=rz*ftx-rx*ftz;
                vz=rx*fty-ry*ftx;
                fact=1.0/sdot_SQ;
                ftx=(sy*vz-sz*vy)*fact;
                fty=(sz*vx-sx*vz)*fact;
                ftz=(sx*vy-sy*vx)*fact;

                grad_x[j] -= ftx;
                grad_y[j] -= fty;
                grad_z[j] -= ftz;
                grad_x[k] += ftx;
                grad_y[k] += fty;
                grad_z[k] += ftz;

                srdot=(sx*rx+sy*ry+sz*rz)/sdot_SQ;
                rrx=rx-sx*srdot;
                rry=ry-sy*srdot;
                rrz=rz-sz*srdot;
                rrdot=sqrt(rrx*rrx+rry*rry+rrz*rrz);
                stdot=(sx*tx+sy*ty+sz*tz)/sdot_SQ;
                ttx=tx-sx*stdot;
                tty=ty-sy*stdot;
                ttz=tz-sz*stdot;
                ttdot=sqrt(ttx*ttx+tty*tty+ttz*ttz);
                fact=(rrdot*fpdot)/(ttdot*qdot);
                fqx=qx*fact;
                fqy=qy*fact;
                fqz=qz*fact;

                grad_x[l] += fqx;
                grad_y[l] += fqy;
                grad_z[l] += fqz;

                grad_x[j] += (fpx*(1.0+srdot)+fqx*stdot);
                grad_y[j] += (fpy*(1.0+srdot)+fqy*stdot);
                grad_z[j] += (fpz*(1.0+srdot)+fqz*stdot);
		
                ftx=fqx*(1.0+stdot)+fpx*srdot;
                fty=fqy*(1.0+stdot)+fpy*srdot;
                ftz=fqz*(1.0+stdot)+fpz*srdot;

                   grad_x[k] -= ftx;
                   grad_y[k] -= fty;
                   grad_z[k] -= ftz;
		}
		
	
       grad_x[i]=0.0;
       grad_y[i]=0.0;
       grad_z[i]=0.0;
	}

	return;	
}
*/

void CMol::BuildDistanceMatrix(void)
{
	int i, j, k, l, iPos, Atom_i, Atom_j, Atom_k, Atom_l, Atom_Host;

	//start	to build the bonded atom arrays for each atom
	memset(AtomBond, 0, sizeof(int)*nAtom);
	memset(Bond_Array, 0, sizeof(int)*MAX_ATOM*MAX_ATOM_BOND);

	iPos = 0;
	for(i=0; i<nBond; i++, iPos+=2)	{
		Atom_i = BondList[iPos  ];
		Atom_j = BondList[iPos+1];

		Bond_Array[Atom_i][AtomBond[Atom_i]] = Atom_j;	//add atom_j into the list of atom_i
		AtomBond[Atom_i]++;

		Bond_Array[Atom_j][AtomBond[Atom_j]] = Atom_i;	//add atom_i into the list of atom_j
		AtomBond[Atom_j]++;
	}
	//end	to build the bonded atom arrays for each atom


	//start	to constrcut distance matrix
	memset(DistMatrix, 99, sizeof(char)*MAX_ATOM*MAX_ATOM);	//99, an arbitrary large numer
	for(i=0; i<nAtom; i++)	{	// i-i, itself. |d| = zero
		DistMatrix[i][i] = 0;
	}

	for(i=0; i<nAtom; i++)	{	//start	enumeration
		Atom_i = i;

		for(j=0; j<AtomBond[Atom_i]; j++)	{	//i->j
			Atom_j = Bond_Array[Atom_i][j];

			if(DistMatrix[Atom_i][Atom_j] > 1)	{
				DistMatrix[Atom_i][Atom_j] = 1;
				DistMatrix[Atom_j][Atom_i] = 1;
			}

			for(k=0; k<AtomBond[Atom_j]; k++)	{	//i->j->k
				Atom_k = Bond_Array[Atom_j][k];

				if(DistMatrix[Atom_i][Atom_k] > 2)	{
					DistMatrix[Atom_i][Atom_k] = 2;
					DistMatrix[Atom_k][Atom_i] = 2;
				}
				
				for(l=0; l<AtomBond[Atom_k]; l++)	{	//i->j->k->l
					Atom_l = Bond_Array[Atom_k][l];
					
					if(DistMatrix[Atom_i][Atom_l] > 3)	{
						DistMatrix[Atom_i][Atom_l] = 3;
						DistMatrix[Atom_l][Atom_i] = 3;
					}
				}
			}
		}
	}
	//end	to constrcut distance matrix

	//start	to substitute the exclusion list of host atom to drudes' list
	for(i=0; i<nDrude; i++)	{
		Atom_i = DrudeList[i];
		Atom_Host = Atom_i - 1;	//the previous one is supposed to be the host atom

		for(j=0; j<nAtom; j++)	{
			if(DistMatrix[Atom_Host][j] < DistMatrix[Atom_i][j])	{
				DistMatrix[Atom_i][j] = DistMatrix[Atom_Host][j];
				DistMatrix[j][Atom_i] = DistMatrix[Atom_i][j];
			}
		}
	}
	//end	to substitute the exclusion list of host atom to drudes' list

	//start	to substitute the exclusion list of host atom to LPs' list
	for(i=0; i<nLP; i++)	{
		Atom_i = LPList[i];
		Atom_Host = LP_Host[i];	//the previous one is supposed to be the host atom

		for(j=0; j<nAtom; j++)	{
			if(DistMatrix[Atom_Host][j] < DistMatrix[Atom_i][j])	{
				DistMatrix[Atom_i][j] = DistMatrix[Atom_Host][j];
				DistMatrix[j][Atom_i] = DistMatrix[Atom_Host][j];
			}
		}
	}
	//end	to substitute the exclusion list of host atom to LPs' list


	n_NB_Pair = 0;
	n_Thole_Pair = 0;
	for(i=0; i<nAtom; i++)	{
		for(j=i+1; j<nAtom; j++)	{
			if(DistMatrix[i][j] >= 3)	{	//only 1-2, 1-3 are excluded 
				NB_List_i[n_NB_Pair] = i;
				NB_List_j[n_NB_Pair] = j;
				n_NB_Pair++;
			}
			else if(DistMatrix[i][j] > 0)	{	//not themself
				if( ((i+1) < nAtom) && ((j+1) < nAtom) )	{
					if( IsDrude[i+1] &&  IsDrude[j+1] )	{
						Thole_List_i[n_Thole_Pair] = i;	//save the index of host atom 1
						Thole_List_j[n_Thole_Pair] = j;	//save the index of host atom 2

						Para_Thole_aa[n_Thole_Pair] = alpha[i]*alpha[j]/pow(thole[i]+thole[j], 6.0);
						Para_Thole_qq[n_Thole_Pair] = MD_COULOMB*CG[i+1]*CG[j+1];

//						printf("%d %d %lf  %lf\n", i+1, j+1, Para_Thole_aa[n_Thole_Pair], Para_Thole_qq[n_Thole_Pair]);

						n_Thole_Pair++;
					}
				}
			}
		}

//		printf("%3d  %6d\n", i+1, n_NB_Pair);
	}

}

void CMol::FixAllRealAtoms(void)
{
	int i;

	for(i=0; i<nAtom; i++)	{
		if(IsDrude[i])	{	//allow drude to relax
			IsFixed[i] = 0;
		}
		else if(IsLonePair[i])	{	//Lone pairs are always allowed to relax.
			IsFixed[i] = 0;
		}
		else	{
			IsFixed[i] = 1;	//fix all other atoms
		}
	}
}


void CMol::BuildActiveEnergyTermList()
{
	int i, iPos;

	//start	to setup IsFixed[]
	if(Geo_Opt_Drude_Only)	{	//optimize drude only, fix all other atoms
		FixAllRealAtoms();
	}
	else	{	//all atoms are allowed to freely relax (In fact, lone pairs will be constrained by Position_LonePair().)
		for(i=0; i<nAtom; i++)	{
			IsFixed[i] = 0;
		}
	}
	//end	to setup IsFixed[]


	nActive_Bond = nActive_Angle = nActive_Dihedral = nActive_ImpDih = nActive_NonBond = nActive_CMap = nActive_Aniso = 0;

	//start	to build the list of active bond
	iPos = 0;
	for(i=0; i<nBond; i++, iPos+=2)	{
		if( (IsFixed[BondList[iPos]]) && (IsFixed[BondList[iPos+1]]))	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_Bond[nActive_Bond] = i;	//add to the active list
			nActive_Bond++;
		}
	}
	//end	to build the list of active bond


	//start	to build the list of active angle
	iPos = 0;
	for(i=0; i<nAngle; i++, iPos+=3)	{
		if( (IsFixed[AngleList[iPos]]) && (IsFixed[AngleList[iPos+1]]) && (IsFixed[AngleList[iPos+2]]))	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_Angle[nActive_Angle] = i;	//add to the active list
			nActive_Angle++;
		}
	}
	//end	to build the list of active angle


	//start	to build the list of active dihedral
	iPos = 0;
	for(i=0; i<nDihedral; i++, iPos+=4)	{
		if( (IsFixed[DihedralList[iPos]]) && (IsFixed[DihedralList[iPos+1]]) && (IsFixed[DihedralList[iPos+2]]) && (IsFixed[DihedralList[iPos+3]]))	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_Dihedral[nActive_Dihedral] = i;	//add to the active list
			nActive_Dihedral++;
		}
	}
	//end	to build the list of active dihedral


	//start	to build the list of active improper dihedral
	iPos = 0;
	for(i=0; i<nImpro; i++, iPos+=4)	{
		if( (IsFixed[ImprDihedralList[iPos]]) && (IsFixed[ImprDihedralList[iPos+1]]) && (IsFixed[ImprDihedralList[iPos+2]]) && (IsFixed[ImprDihedralList[iPos+3]]))	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_ImpDih[nActive_ImpDih] = i;	//add to the active list
			nActive_ImpDih++;
		}
	}
	//end	to build the list of active improper dihedral


	//start	to build the list of active non-bonded pairs
	for(i=0; i<n_NB_Pair; i++)	{
		if( (IsFixed[NB_List_i[i]]) && (IsFixed[NB_List_j[i]]) )	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_NonBond[nActive_NonBond] = i;	//add to the active list
			nActive_NonBond++;
		}
	}
	//end	to build the list of active non-bonded pairs

	//start	to build the list of active CMap
	for(i=0; i<nCMapTerm; i++)	{
		if( (IsFixed[CMapList[i][0]]) && (IsFixed[CMapList[i][1]]) && (IsFixed[CMapList[i][2]]) && (IsFixed[CMapList[i][3]]) && 
			 (IsFixed[CMapList[i][4]]) && (IsFixed[CMapList[i][5]]) &&  (IsFixed[CMapList[i][6]]) && (IsFixed[CMapList[i][7]]) )	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_CMap[nActive_CMap] = i;	//add to the active list
			nActive_CMap++;
		}
	}
	//end	to build the list of active CMap

	//start	to build the list of anisotropy
	for(i=0; i<nAniso; i++)	{
		if( (IsFixed[AnisoList[i][0]]) && (IsFixed[AnisoList[i][0]+1]) && (IsFixed[AnisoList[i][1]]) && (IsFixed[AnisoList[i][2]]) && 
			 (IsFixed[AnisoList[i][3]]) )	{	//both atoms are fixed, do not need calculation
		}
		else	{
			Active_List_Aniso[nActive_Aniso] = i;	//add to the active list
			nActive_Aniso++;
		}
	}
	//end	to build the list of anisotropy

	return;
}

void CMol::Setup_NonBondParameters(void)
{
	int i, Atom_i, Atom_j, NBFix_Idx, LJ_Idx;
	LJ_REC *pLJ_Rec;

	memset(Active_LJ, 0, sizeof(int)*N_REC_MAX);

	//start	to assign parameters for LJ
	pLJ_Rec = myForceField->LJ_Rec;
	for(i=0; i<nAtom; i++)	{
		LJ_Idx = LJ_Para_Rec[i];
		if(LJ_Idx >= 0)	{
			Para_LJ_Epsilon[i] = pLJ_Rec[LJ_Idx].para[1];		// E_Min
			Para_LJ_Sigma[i] = pLJ_Rec[LJ_Idx].para[2];			// r_Min
			Para_LJ_Epsilon_14[i] = pLJ_Rec[LJ_Idx].para[4];	// E_Min_14
			Para_LJ_Sigma_14[i] = pLJ_Rec[LJ_Idx].para[5];		// r_Min_14

			// Approximation to decrease the number of variables !!!
			// assume Epsilon_14 == Epsilon .or.  Epsilon_14 == 0.5*Epsilon
			if(Para_LJ_Epsilon_14[i] != Para_LJ_Epsilon[i])	{
				Para_LJ_Epsilon_14[i] = 0.5 * Para_LJ_Epsilon[i];
			}

			Active_LJ[LJ_Idx] = 1;
			Para_LJ_Sqrt_Epsilon[i] = sqrt(fabs(Para_LJ_Epsilon[i]));
			Para_LJ_Sqrt_Epsilon_14[i] = sqrt(fabs(Para_LJ_Epsilon_14[i]));
		}
		else	{
			Para_LJ_Epsilon[i] = Para_LJ_Sigma[i] = Para_LJ_Epsilon_14[i] = Para_LJ_Sigma_14[i] = 0.0;		// all zero
		}
	}
	//end	to assign parameters for LJ

	//start	to check active LJ record
	nActiveLJ = 0;
	for(i=0; i<myForceField->n_Rec_LJ; i++)	{
		if(Active_LJ[i])	{
			printf("Active LJ: %s\n", pLJ_Rec[i].Chem);
			nActiveLJ++;
		}
	}
	printf("There are %3d active LJ types.\n", nActiveLJ);
	//end	to check active LJ record


	for(i=0; i<n_NB_Pair; i++)	{
		Atom_i = NB_List_i[i];
		Atom_j = NB_List_j[i];

		if(DistMatrix[Atom_i][Atom_j] == 3)	{	//1-4 interactions
			Para_LJ_Epsilon_IJ[i] = sqrt( Para_LJ_Epsilon_14[Atom_i] * Para_LJ_Epsilon_14[Atom_j] );
			Para_LJ_Sigma_IJ[i] = ( Para_LJ_Sigma_14[Atom_i] + Para_LJ_Sigma_14[Atom_j] );
		}
		else	{
			Para_LJ_Epsilon_IJ[i] = sqrt( Para_LJ_Epsilon[Atom_i] * Para_LJ_Epsilon[Atom_j] );
			Para_LJ_Sigma_IJ[i] = ( Para_LJ_Sigma[Atom_i] + Para_LJ_Sigma[Atom_j] );
		}

		NBFix_Idx = NBFix_Rec[Atom_i][Atom_j];
		if(NBFix_Idx >= 0)	{	// setup the NBFix parameters from the force field
			Para_LJ_Epsilon_IJ[i] = fabs(myForceField->NBFix_Rec[NBFix_Idx].para[0]);
			Para_LJ_Sigma_IJ[i] = myForceField->NBFix_Rec[NBFix_Idx].para[1];
		}

		//printf("%s %f \n","Para_LJ_Sigma_IJ[i]",Para_LJ_Sigma_IJ[i]);
		Para_LJ_Sigma_IJ_pow_6[i] = Para_LJ_Sigma_IJ[i]*Para_LJ_Sigma_IJ[i]*Para_LJ_Sigma_IJ[i];
		Para_LJ_Sigma_IJ_pow_6[i] *= Para_LJ_Sigma_IJ_pow_6[i];	// ^6
		Para_LJ_Sigma_IJ_pow_12[i] = Para_LJ_Sigma_IJ_pow_6[i]*Para_LJ_Sigma_IJ_pow_6[i];	// ^12
	
		Para_Elec_Pair[i] = MD_COULOMB * CG[Atom_i] * CG[Atom_j];
	}
	return;
}

void CMol::Setup_TholePairParameters(void)
{
	int iPair, Atom_1, Atom_2;

	for(iPair=0; iPair<n_Thole_Pair; iPair++)	{
		Atom_1 = Thole_List_i[iPair];
		Atom_2 = Thole_List_j[iPair];
		
		Para_Thole_aa[iPair] = alpha[Atom_1]*alpha[Atom_2]/pow(thole[Atom_1]+thole[Atom_2], 6.0);
		Para_Thole_qq[iPair] = MD_COULOMB*CG[Atom_1+1]*CG[Atom_2+1];

//		printf("%d %d %lf  %lf\n", Atom_1+1, Atom_2+1, Para_Thole_aa[iPair], Para_Thole_qq[iPair]);
	}
}

void CMol::TranslateLastSegment(void)
{
	int i;

	if(StartLast< 0)	{
		Quit_With_Error_Msg("Fatal error!\nTrying to translate the last segment, but there is only one segment.\nQuit.\n");
	}
	for(i=StartLast; i<=EndLast; i++)	{
		x[i] += FAR_DIST;
//		y[i] += FAR_DIST;
//		z[i] += FAR_DIST;
	}
}

void CMol::TranslateLastSegmentBack(void)
{
	int i;

	if(StartLast< 0)	{
		Quit_With_Error_Msg("Fatal error!\nTrying to translate the last segment back, but there is only one segment.\nQuit.\n");
	}
	for(i=StartLast; i<=EndLast; i++)	{
		x[i] = x_save[i];
		y[i] = y_save[i];
		z[i] = z_save[i];
	}
}

void CMol::BackupCoordinates(void)
{
	int i;

	for(i=0; i<nAtom; i++)	{
		x_save[i] = x[i];
		y_save[i] = y[i];
		z_save[i] = z[i];
	}
}

void CMol::RestoreCoordinates(void)
{
	int i;

	for(i=0; i<nAtom; i++)	{
		x[i] = x_save[i];
		y[i] = y_save[i];
		z[i] = z_save[i];
	}
}

void CMol::WriteCRDFile(char szName[])
{
	FILE *fOut;
	int i;

	fOut = fopen(szName, "w");

	fprintf(fOut, "* CRD file generated by Lei Huang's code.\n");
	fprintf(fOut, "%10d  EXT\n", nAtom);

	for(i=0; i<nAtom; i++)	{
		fprintf(fOut, "%10d%10d  %-7s   %-8s%20.10lf%20.10lf%20.10lf  %-7s   1               0.0000000000\n", 
			i+1, SegID[i], ResName[i], AtomName[i], x[i], y[i], z[i], ResName[i]);
	}

	fclose(fOut);
}


void CMol::Setup_NBFix(void)
{
	int n_NBFix_Rec;
	int i, j, k, IsExistInNBFixRec[MAX_ATOM];
	NBFix_REC *pRec;

	n_NBFix_Rec = myForceField->n_Rec_NBFix;
	pRec = myForceField->NBFix_Rec;

	for(i=0; i<nAtom; i++)	{
		IsExistInNBFixRec[i] = 0;
		for(j=0; j<n_NBFix_Rec; j++)	{
			if( (strcmp(pRec[j].Chem_1, ChemName[i]) == 0) || (strcmp(pRec[j].Chem_2, ChemName[i]) == 0) )	{
				IsExistInNBFixRec[i]++;
			}
		}
	}

	for(i=0; i<nAtom; i++)	{
		for(j=0; j<nAtom; j++)	{
			NBFix_Rec[i][j] = -1;	// no NBFix parameters defined
		}
	}
	for(i=0; i<nAtom; i++)	{
		for(j=i; j<nAtom; j++)	{
			if( IsExistInNBFixRec[i] && IsExistInNBFixRec[j] )	{
				for(k=0; k<n_NBFix_Rec; k++)	{
					if( ( (strcmp(pRec[k].Chem_1, ChemName[i]) == 0) && (strcmp(pRec[k].Chem_2, ChemName[j]) == 0) ) || ( (strcmp(pRec[k].Chem_2, ChemName[i]) == 0) && (strcmp(pRec[k].Chem_1, ChemName[j]) == 0) ) )	{
						NBFix_Rec[i][j] = NBFix_Rec[j][i] = k;
						fprintf(fFile_Run_Log, "Found NBFix parameters for %5s and %5s\n", pRec[k].Chem_1, pRec[k].Chem_2);
					}
				}
			}
		}
	}
}


void ReadDoubleDataFromBuff(char szBuff[], double Data[], int n)
{
	int n_Item;
	int i, iPos, nLen;
	char szData[MAX_LEN_LINE];

	iPos = 0;
	for(i=0; i<n; i++)	{
		while(szBuff[iPos] == ' ')	{	//to find the first non-blank character
			iPos++;
		}

		n_Item = sscanf(szBuff+iPos, "%s", szData);
		if(n_Item != 1)	{
//			printf("Error in reading string:\n%s\nQuit.\n", szBuff);
			fprintf(fFile_Run_Log, "Error in reading string:\n%s\nQuit.\n", szBuff);
			fflush(fFile_Run_Log);
			exit(1);
		}
		nLen = strlen(szData);
		Data[i] = atof(szData);
		iPos += nLen;
	}

	return;
}

void ReadIntDataFromBuff(char szBuff[], int Data[], int n)
{
	int n_Item;
	int i, iPos, nLen;
	char szData[MAX_LEN_LINE];

	iPos = 0;
	for(i=0; i<n; i++)	{
		while(szBuff[iPos] == ' ')	{	//to find the first non-blank character
			iPos++;
		}

		n_Item = sscanf(szBuff+iPos, "%s", szData);
		if(n_Item != 1)	{
//			printf("Error in reading string:\n%s\nQuit.\n", szBuff);
			fprintf(fFile_Run_Log, "Error in reading string:\n%s\nQuit.\n", szBuff);
			fflush(fFile_Run_Log);
			exit(1);
		}
		nLen = strlen(szData);
		Data[i] = atoi(szData);
		iPos += nLen;
	}

	return;
}

CForceField::CForceField()
{
	n_Rec_Bond = n_Rec_Angle = n_Rec_Dihedral = n_Rec_ImproDihedral = n_Rec_LJ = n_Rec_NBFix = 0;
}

CForceField::~CForceField()
{
}

void CForceField::ReadForceField(char szNamePara[])
{
	FILE *fIn;
	char szLine[256], ErrorMsg[256], szComment[256];
	int n_Len, ReadItem, i, j, iTmp, iCount, Is_CMap_On;

	fIn = fopen(szNamePara, "r");
	if(fIn == NULL)	{
//		printf("CForceField::ReadForceField> Fail to open file %s\nQuit\n", szNamePara);
		fprintf(fFile_Run_Log, "CForceField::ReadForceField> Fail to open file %s\nQuit\n", szNamePara);
		fflush(fFile_Run_Log);
		exit(1);
	}

	while(1)	{	//bond parameters
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Error in reading file %s. Fail to find the entry for BONDS.\nQuit\n", szNamePara);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}
		fgets(szLine, MAX_LEN_LINE, fIn);
		if( strncmp(szLine, "BONDS", 5) == 0 )	{
			break;
		}
	}

	while(1)	{
		if(feof(fIn))	{
			sprintf(ErrorMsg, "Error in reading file %s. Fail to find the entry for ANGLES.\nQuit\n", szNamePara);
			fclose(fIn);
			Quit_With_Error_Msg(ErrorMsg);
		}

		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "ANGLES", 6) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %lf %lf", 
			Bond_Rec[n_Rec_Bond].Chem[0], Bond_Rec[n_Rec_Bond].Chem[1], 
			 &(Bond_Rec[n_Rec_Bond].para[0]), &(Bond_Rec[n_Rec_Bond].para[1]));
		if( ReadItem == 4 )	{
			n_Rec_Bond++;
		}
	}
	if(n_Rec_Bond > N_REC_MAX)	{
//		printf("n_Rec_Bond > N_REC_MAX. n_Rec_Bond = %d\nQuit\n", n_Rec_Bond);
		fprintf(fFile_Run_Log, "n_Rec_Bond > N_REC_MAX. n_Rec_Bond = %d\nQuit\n", n_Rec_Bond);
		fflush(fFile_Run_Log);
		exit(1);
	}

	while(1)	{
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "DIHEDRALS", 9) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %lf %lf %lf %lf", 
			Angle_Rec[n_Rec_Angle].Chem[0], Angle_Rec[n_Rec_Angle].Chem[1], Angle_Rec[n_Rec_Angle].Chem[2], 
			 &(Angle_Rec[n_Rec_Angle].para[0]), &(Angle_Rec[n_Rec_Angle].para[1]), 
			 &(Angle_Rec[n_Rec_Angle].para[2]), &(Angle_Rec[n_Rec_Angle].para[3]));
		if( ReadItem == 5 )	{
			Angle_Rec[n_Rec_Angle].para[2] = 0.0;
			Angle_Rec[n_Rec_Angle].para[3] = 0.0;

			n_Rec_Angle++;
		}
		else if( ReadItem == 7 )	{
			n_Rec_Angle++;
		}
	}
	if(n_Rec_Angle > N_REC_MAX)	{
//		printf("n_Rec_Angle > N_REC_MAX. n_Rec_Angle = %d\nQuit\n", n_Rec_Angle);
		fprintf(fFile_Run_Log, "n_Rec_Angle > N_REC_MAX. n_Rec_Angle = %d\nQuit\n", n_Rec_Angle);
		fflush(fFile_Run_Log);
		exit(1);
	}


	while(1)	{
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "IMPROPER", 8) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %s %lf %lf %lf", 
			Dihedral_Rec[n_Rec_Dihedral].Chem[0], Dihedral_Rec[n_Rec_Dihedral].Chem[1], Dihedral_Rec[n_Rec_Dihedral].Chem[2], Dihedral_Rec[n_Rec_Dihedral].Chem[3], 
			 &(Dihedral_Rec[n_Rec_Dihedral].para[0]), &(Dihedral_Rec[n_Rec_Dihedral].para[1]), &(Dihedral_Rec[n_Rec_Dihedral].para[2]));
		if(ReadItem != 7)	{
			continue;
		}
		n_Rec_Dihedral++;
	}
	if(n_Rec_Dihedral > N_REC_MAX)	{
//		printf("n_Rec_Dihedral > N_REC_MAX. n_Rec_Dihedral = %d\nQuit\n", n_Rec_Dihedral);
		fprintf(fFile_Run_Log, "n_Rec_Dihedral > N_REC_MAX. n_Rec_Dihedral = %d\nQuit\n", n_Rec_Dihedral);
		fflush(fFile_Run_Log);
		exit(1);
	}


	while(1)	{
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		if(strncmp(szLine, "CMAP", 4) == 0 )	{
			Is_CMap_On = 1;
			break;
		}
		if(strncmp(szLine, "NONBONDED", 9) == 0 )	{
			Is_CMap_On = 0;
			break;
		}

		ReadItem = sscanf(szLine, "%s %s %s %s %lf %lf %lf", 
			ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[0], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[1], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[2], ImproDihedral_Rec[n_Rec_ImproDihedral].Chem[3], 
			 &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[0]), &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[1]), &(ImproDihedral_Rec[n_Rec_ImproDihedral].para[2]));
		if(ReadItem == 7)	{
			n_Rec_ImproDihedral++;
		}
	}
	if(n_Rec_ImproDihedral > N_REC_MAX)	{
		fprintf(fFile_Run_Log, "n_Rec_ImproDihedral > N_REC_MAX. n_Rec_ImproDihedral = %d\nQuit\n", n_Rec_ImproDihedral);
		fflush(fFile_Run_Log);
		exit(1);
	}

	if(Is_CMap_On)	{
		// only read the first entry for CMAP
		//start	to read CMap parameters 24X24 from topar file
		double CMap_data[5];	//five items per line
		for(i=0; i<CMAP_DIM; i++)	{
			while(1)	{
				fgets(szLine, MAX_LEN_LINE, fIn);
				ReadItem = sscanf(szLine+1, "%d", &iTmp);
				if( (ReadItem == 1) && (szLine[0] == '!') )	{	//the entry of CMap data
					break;
				}
			}
			
			iCount = 0;
			while(1)	{
				fgets(szLine, MAX_LEN_LINE, fIn);
				ReadItem = sscanf(szLine, "%lf %lf %lf %lf %lf", 
					&(CMap_data[0]), &(CMap_data[1]), &(CMap_data[2]), &(CMap_data[3]), &(CMap_data[4]));
				if( ReadItem == 5 )	{
					for(j=0; j<5; j++)	{
						mctp[0][i][iCount+j] = CMap_data[j];
					}
					iCount += 5;
				}
				else if( (iCount == 20) && (ReadItem == 4) )	{
					for(j=0; j<5; j++)	{
						mctp[0][i][iCount+j] = CMap_data[j];
					}
					iCount += 4;
					break;
				}
				else	{
					fprintf(fFile_Run_Log, "Please check your parameter file around the entry of CMap.\n");
					fflush(fFile_Run_Log);
					exit(1);
				}
			}
		}
		//	setcmap(24, 12, mctp, 360.0/24.0);
		setcmap(24, 12, mctp, 15.0);
		//end	to read CMap parameters 24X24 from topar file
		
		while(1)	{
			if(feof(fIn))	{
				fprintf(fFile_Run_Log, "Fail to find the entry for NONBONDED parameters.\nQuit\n");
				fflush(fFile_Run_Log);
				exit(1);
			}
			fgets(szLine, MAX_LEN_LINE, fIn);
			
			if(strncmp(szLine, "NONBONDED", 9) == 0 )	{
				break;
			}
		}
	}


	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}
		if(strncmp(szLine, "NBFIX", 5) == 0 )	{
			break;
		}

		ReadItem = sscanf(szLine, "%s %lf %lf %lf %lf %lf %lf", 
			LJ_Rec[n_Rec_LJ].Chem, 
			&(LJ_Rec[n_Rec_LJ].para[0]), &(LJ_Rec[n_Rec_LJ].para[1]), &(LJ_Rec[n_Rec_LJ].para[2]), 	//LJ parameters
			&(LJ_Rec[n_Rec_LJ].para[3]), &(LJ_Rec[n_Rec_LJ].para[4]), &(LJ_Rec[n_Rec_LJ].para[5]));	//LJ 1-4 parameters
		if(ReadItem == 4)	{
			LJ_Rec[n_Rec_LJ].para[3] = LJ_Rec[n_Rec_LJ].para[0];	//assign 1-4 parameters same as non-1-4 parameters defaultly
			LJ_Rec[n_Rec_LJ].para[4] = LJ_Rec[n_Rec_LJ].para[1];
			LJ_Rec[n_Rec_LJ].para[5] = LJ_Rec[n_Rec_LJ].para[2];

			if(Extract_Real_VDW_Type(szLine, szComment, LJ_Rec[n_Rec_LJ].RealChem) == 0)	{
				strcpy(LJ_Rec[n_Rec_LJ].RealChem, LJ_Rec[n_Rec_LJ].Chem);
			}

			n_Rec_LJ++;
		}
		else if(ReadItem == 7)	{	//with full parameters
			if(Extract_Real_VDW_Type(szLine, szComment, LJ_Rec[n_Rec_LJ].RealChem) == 0)	{
				strcpy(LJ_Rec[n_Rec_LJ].RealChem, LJ_Rec[n_Rec_LJ].Chem);
			}

			n_Rec_LJ++;
		}
		if(n_Rec_LJ > N_REC_MAX)	{
			fprintf(fFile_Run_Log, "n_Rec_LJ > N_REC_MAX. n_Rec_LJ = %d\nQuit\n", n_Rec_LJ);
			fflush(fFile_Run_Log);
			exit(1);
		}
	}

	while(1)	{
		if(feof(fIn))	{
			break;
		}
		fgets(szLine, MAX_LEN_LINE, fIn);
		if(szLine[0] == '!')	{
			continue;
		}
		n_Len = strlen(szLine);
		if(n_Len < 2)	{ //empty line
			continue;
		}

		ReadItem = sscanf(szLine, "%s %s %lf %lf", 
			NBFix_Rec[n_Rec_NBFix].Chem_1, NBFix_Rec[n_Rec_NBFix].Chem_2, 
			&(NBFix_Rec[n_Rec_NBFix].para[0]), &(NBFix_Rec[n_Rec_NBFix].para[1]));
		if(ReadItem == 4)	{
			n_Rec_NBFix++;
			if(n_Rec_NBFix > N_REC_MAX)	{
				fprintf(fFile_Run_Log, "n_Rec_NBFix > N_REC_MAX. n_Rec_NBFix = %d\nQuit\n", n_Rec_NBFix);
				fflush(fFile_Run_Log);
				exit(1);
			}
		}
	}


	fclose(fIn);
}

int CForceField::Extract_Real_VDW_Type(char szBuff[], char szComment[], char szRealName[])
{
	int i, nLen, ReadItem;
	char szTmp[256], szTmp_2[256];

	nLen = strlen(szBuff);
	for(i=nLen-1; i>=0; i--)	{
		if(szBuff[i] == '!')	{
			break;
		}
	}

	szComment[0] = 0;
	szRealName[0] = 0;

	if(i<=3)	{
		return 0;
	}

	strcpy(szComment, szBuff+i);

	ReadItem = sscanf(szBuff+i+1, "%s%s", szTmp, szTmp_2);
	if( (ReadItem == 2) && (strncmp(szTmp, "ORG_VDW", 7)==0) )	{
		strcpy(szRealName, szTmp_2);
		return 1;
	}

	return 0;
}

int Compare_Chem_Type(char Type_1[], char Type_2[])
{
	if( (strcmp(Type_1, "X")==0) || (strcmp(Type_2, "X")==0) )	{	// same type
		return 0;
	}
	return (strcmp(Type_1, Type_2));
}

void CForceField::GetPara_Bond(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	if( (strcmp(szChemName[0], "LP")==0) || (strcmp(szChemName[1], "LP")==0) )	{
		Para[0] = Para[1] = 0.0;
		return;
	}

	for(i=0; i<n_Rec_Bond; i++)	{
		if( (strcmp(szChemName[0], Bond_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Bond_Rec[i].Chem[1])==0) )	{
			Para[0] = Bond_Rec[i].para[0];
			Para[1] = Bond_Rec[i].para[1];
			return;
		}
	}

	for(i=0; i<n_Rec_Bond; i++)	{
		if( (strcmp(szChemName[1], Bond_Rec[i].Chem[0])==0) && (strcmp(szChemName[0], Bond_Rec[i].Chem[1])==0) )	{
			Para[0] = Bond_Rec[i].para[0];
			Para[1] = Bond_Rec[i].para[1];
			return;
		}
	}

//	printf("\n\nFatal error!\nCan't find the parameters for bond %s - %s\nQuit.\n", szChemName[0], szChemName[1]);
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for bond %s - %s\nQuit.\n", szChemName[0], szChemName[1]);
	fflush(fFile_Run_Log);
	exit(1);

	return;
}

void CForceField::ReadUpdatedCMap(void)
{
	int i, j;
	
	//start	to read newly generated CMap
	FILE *fIn;
	fIn = fopen("my-new-cmap.dat", "r");
	if(fIn != NULL)	{	// this updated cmap exists
		for(i=0; i<CMAP_DIM; i++)	{
			for(j=0; j<CMAP_DIM; j++)	{
				fscanf(fIn, "%lf", &(mctp[0][i][j]));
			}
		}
	}
	fclose(fIn);
	//end	to read newly generated CMap

	setcmap(24, 12, mctp, 15.0);

}

void CForceField::GetPara_Angle(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (strcmp(szChemName[0], Angle_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Angle_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (strcmp(szChemName[2], Angle_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Angle_Rec[i].Chem[1])==0) && (strcmp(szChemName[0], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (Compare_Chem_Type(szChemName[0], Angle_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Angle_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];
			return;
		}
	}

	for(i=0; i<n_Rec_Angle; i++)	{
		if( (Compare_Chem_Type(szChemName[2], Angle_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Angle_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[0], Angle_Rec[i].Chem[2])==0) )	{
			Para[0] = Angle_Rec[i].para[0];
			Para[1] = Angle_Rec[i].para[1];
			Para[2] = Angle_Rec[i].para[2];
			Para[3] = Angle_Rec[i].para[3];
			return;
		}
	}


//	printf("\n\nFatal error!\nCan't find the parameters for angle %s, %s, %s\nQuit.\n", 
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for angle %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2]);
	fflush(fFile_Run_Log);
	exit(1);

	return;
}

void CForceField::GetPara_Dihedral(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i, Idx, Count, iPos;

	Count = 0;
	iPos = 0;

	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (strcmp(szChemName[0], Dihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], Dihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], Dihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[3], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}
	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (strcmp(szChemName[3], Dihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[2], Dihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[1], Dihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[0], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}

	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[0], Dihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], Dihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], Dihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[3], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}
	for(i=0; i<n_Rec_Dihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[3], Dihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[2], Dihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[1], Dihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[0], Dihedral_Rec[i].Chem[3])==0))	{
			Idx = (int)(Dihedral_Rec[i].para[1] + 1.0E-6);
			iPos = Idx*3;
			Para[iPos  ] = Dihedral_Rec[i].para[0];
			Para[iPos+1] = Dihedral_Rec[i].para[1];
			Para[iPos+2] = Dihedral_Rec[i].para[2];

			Count++;
		}
	}
	if(Count > 0)	{
		return;
	}


//	printf("\n\nFatal error!\nCan't find the parameters for dihedral %s, %s, %s, %s\nQuit.\n", 
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for dihedral %s, %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2], szChemName[3]);
	fflush(fFile_Run_Log);
	exit(1);

	return;
}

void CForceField::GetPara_ImproDIhedral(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i;

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (strcmp(szChemName[0], ImproDihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[1], ImproDihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[2], ImproDihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[3], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (strcmp(szChemName[3], ImproDihedral_Rec[i].Chem[0])==0) && (strcmp(szChemName[2], ImproDihedral_Rec[i].Chem[1])==0) && (strcmp(szChemName[1], ImproDihedral_Rec[i].Chem[2])==0) && (strcmp(szChemName[0], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[0], ImproDihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[1], ImproDihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[2], ImproDihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[3], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type
			return;
		}
	}

	for(i=0; i<n_Rec_ImproDihedral; i++)	{
		if( (Compare_Chem_Type(szChemName[3], ImproDihedral_Rec[i].Chem[0])==0) && (Compare_Chem_Type(szChemName[2], ImproDihedral_Rec[i].Chem[1])==0) && (Compare_Chem_Type(szChemName[1], ImproDihedral_Rec[i].Chem[2])==0) && (Compare_Chem_Type(szChemName[0], ImproDihedral_Rec[i].Chem[3])==0))	{
			Para[0] = ImproDihedral_Rec[i].para[0];
			Para[1] = ImproDihedral_Rec[i].para[2];
			Para[2] = ImproDihedral_Rec[i].para[1];	// potential type
			return;
		}
	}

//	printf("\n\nFatal error!\nCan't find the parameters for improper dihedral %s, %s, %s, %s\nQuit.\n", 
	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for improper dihedral %s, %s, %s, %s\nQuit.\n", 
		szChemName[0], szChemName[1], szChemName[2], szChemName[3]);
	fflush(fFile_Run_Log);
	exit(1);

	return;
}

int CForceField::GetPara_LJ(char szChemName[][N_LEN_CHEM_NAME], double Para[])
{
	int i, j;

	if(strcmp(szChemName[0], "-----") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	if(strcmp(szChemName[0], "DOH2") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}
	if(strcmp(szChemName[0], "DRUD") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	if(strcmp(szChemName[0], "LP") == 0)	{	//drude particle
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return (-1);
	}

	for(i=0; i<n_Rec_LJ; i++)	{
		if( strcmp(szChemName[0], LJ_Rec[i].Chem)==0 )	{
			Para[0] = LJ_Rec[i].para[0];
//			Para[1] = LJ_Rec[i].para[1];
			Para[1] =-LJ_Rec[i].para[1];
			Para[2] = LJ_Rec[i].para[2];
			Para[3] = LJ_Rec[i].para[3];
//			Para[4] = LJ_Rec[i].para[4];
			Para[4] =-LJ_Rec[i].para[4];
			Para[5] = LJ_Rec[i].para[5];

			for(j=0; j<n_Rec_LJ; j++)	{
				if(strcmp(LJ_Rec[j].RealChem, LJ_Rec[i].RealChem)==0)	{	// find the first LJ type with desired real vdw name
					return j;
				}
			}

			return i;
		}
	}




	if(strcmp(szChemName[0], "CALD")==0)	{	// it might be the test charge
		fprintf(fFile_Run_Log, "Warning: Can't find the parameters for LJ, CALD. The parameters will be set as zero.\n\n");
		Para[0] = Para[1] = Para[2] = 0.0;
		Para[3] = Para[4] = Para[5] = 0.0;
		return -1;
	}

	fprintf(fFile_Run_Log, "\n\nFatal error!\nCan't find the parameters for LJ, %s\nQuit.\n", 
		szChemName[0]);
	fflush(fFile_Run_Log);
//	exit(1);

	return (-1);
}

int  CForceField::QueryNBFix(char szChem_1[], char szChem_2[])
{
	int i;

	for(i=0; i<n_Rec_NBFix; i++)	{
		if( ((strcmp(szChem_1, NBFix_Rec[i].Chem_1)==0) && (strcmp(szChem_2, NBFix_Rec[i].Chem_2)==0)) || ((strcmp(szChem_2, NBFix_Rec[i].Chem_1)==0) && (strcmp(szChem_1, NBFix_Rec[i].Chem_2)==0)) )	{
			return i;
		}
	}
	return -1;
}

void cmapspl(double dx, double y[], int n, double u[], double y2[])
{
	int i, iMax;
    double pinv;
    double dxinv;
	
    y2[0]=0.0;
    u[0]=0.0;
	
    dxinv=1.0/dx;
	iMax = n-2;
	for(i=1; i<=iMax; i++)	{
		pinv=1.0/(y2[i-1]+4.0);
		y2[i] = -pinv;
		u[i]=((6.0*y[i+1]-12.0*y[i]+6.0*y[i-1])*dxinv*dxinv-u[i-1])*pinv;
	}
	
    y2[n-1]=0.0;
	
	for(i=n-2; i>=0; i--)	{
		y2[i]=y2[i]*y2[i+1]+u[i];
	}
	
	return;
}

void cmapspi(double xmin, double dx, double ya[], double y2a[], double x, double& y, double& y1)
{
	int inx;
	double a,b;
	
	//    inx=nint((x-xmin)/dx-half)+1		//Please check!!!
    inx=(int)floor((x-xmin)/dx);
	
    a=(xmin+(inx+1)*dx-x)/dx;
	//    b=(x-xmin-(inx-1)*dx)/dx;
    b=1.0-a;
	
    y=a*ya[inx]+b*ya[inx+1]+ ((a*a*a-a)*y2a[inx]+(b*b*b-b)*y2a[inx+1])*(dx*dx)/6.0;
    y1=(ya[inx+1]-ya[inx])/dx - (3.0*a*a-1.0)/6.0*dx*y2a[inx]+ (3.0*b*b-1.0)/6.0*dx*y2a[inx+1];
	
    return;
}


void setcmap(int num, int xm, double gmap[][CMAP_DIM][CMAP_DIM], double dx)
{
	double **tgmap;
	double **t2;
	double *u,*u2;
	double *yytmp,*y1tmp;
    double phi,psi;
    double xmin;
    double v,v1,v2,v12;
    int i,iMax,j,jMax,k,kMax,ii,jj;

	//start	to for test only
//	FILE *fOut;
//	fOut = fopen("cmap-tmp.dat", "w");
//	for(i=0; i<CMAP_DIM; i++)	{
//		for(j=0; j<CMAP_DIM; j++)	{
//			fprintf(fOut, "%lf\n", gmap[0][i][j]);
//		}
//	}
//	fclose(fOut);
	//end	to for test only

	
	tgmap = dmatrix(1,num+xm+xm,1,num+xm+xm);	//like Fortran, starting from 1 !!!
	t2 = dmatrix(1,num+xm+xm,1,num+xm+xm);	//like Fortran, starting from 1 !!!
	u = dvector(1,num+xm+xm);
	u2 = dvector(1,num+xm+xm);
	yytmp = dvector(1,num+xm+xm);
	y1tmp = dvector(1,num+xm+xm);
	
	memset(u+1, 0, sizeof(double)*(num+xm+xm));
	memset(u2+1, 0, sizeof(double)*(num+xm+xm));
    xmin=-180.0-xm*dx;
	iMax = num+xm+xm;
	jMax = num+xm+xm;
	for(i=1; i<=iMax; i++)	{
		ii=(i+num-xm-1)%num;
		for(j=1; j<=jMax; j++)	{
			jj=(j+num-xm-1)%num;
			tgmap[i][j] = gmap[0][ii][jj];
		}
	}
	
	for(i=1; i<=iMax; i++)	{
		cmapspl(dx,tgmap[i]+1,num+xm+xm,u+1,t2[i]+1);	//please check
	}
	
	iMax = num+xm;
	jMax = num+xm;
	kMax = num+xm+xm;
	for(i=1+xm; i<=iMax; i++)	{
		phi=(i-xm-1)*dx-180.0;
		for(j=1+xm; j<=jMax; j++)	{
			psi=(j-xm-1)*dx-180.0;
			for(k=1; k<=kMax; k++)	{
				cmapspi(xmin,dx,tgmap[k]+1, t2[k]+1,psi,yytmp[k],y1tmp[k]);	//please check
			}
			
			cmapspl(dx,yytmp+1,num+xm+xm,u+1,u2+1);	//please check
			cmapspi(xmin,dx,yytmp+1,u2+1,phi,v,v1);
			cmapspl(dx,y1tmp+1,num+xm+xm,u+1,u2+1);
			cmapspi(xmin,dx,y1tmp+1,u2+1,phi,v2,v12);
			gmap[1][i-xm-1][j-xm-1] = v1;
			gmap[2][i-xm-1][j-xm-1] = v2;
			gmap[3][i-xm-1][j-xm-1] = v12;
		}
	}
	
	free_dmatrix(tgmap,1,num+xm+xm,1,num+xm+xm);
	free_dmatrix(t2,1,num+xm+xm,1,num+xm+xm);
	free_dvector(u,1,num+xm+xm);
	free_dvector(u2,1,num+xm+xm);
	free_dvector(yytmp,1,num+xm+xm);
	free_dvector(y1tmp,1,num+xm+xm);
	
	return;
}


void gcstup2(double ggrd[][CMAP_DIM][CMAP_DIM], double ty[4], double ty1[4], double ty2[4], double ty12[4], int ip1, int ip2)
{
	int ip1p1,ip2p1;
	
    ip1p1 = (ip1+1)%CMAP_DIM; //????????????????
    ip2p1 = (ip2+1)%CMAP_DIM; //????????????????
	
    ty[0]=ggrd[0][ip1][ip2];
    ty[1]=ggrd[0][ip1p1][ip2];
    ty[2]=ggrd[0][ip1p1][ip2p1];
    ty[3]=ggrd[0][ip1][ip2p1];
	
    ty1[0]=ggrd[1][ip1][ip2];
    ty1[1]=ggrd[1][ip1p1][ip2];
    ty1[2]=ggrd[1][ip1p1][ip2p1];
    ty1[3]=ggrd[1][ip1][ip2p1];
	
    ty2[0]=ggrd[2][ip1][ip2];
    ty2[1]=ggrd[2][ip1p1][ip2];
    ty2[2]=ggrd[2][ip1p1][ip2p1];
    ty2[3]=ggrd[2][ip1][ip2p1];
	
    ty12[0]=ggrd[3][ip1][ip2];
    ty12[1]=ggrd[3][ip1p1][ip2];
    ty12[2]=ggrd[3][ip1p1][ip2p1];
    ty12[3]=ggrd[3][ip1][ip2p1];
	
}


void gcscf(double ty[4], double ty1[4], double ty2[4], double ty12[4], double gres1, double gres2, double tc[][4])
{
	int i,j,k,i_n;
	double xx,tx[16];
	int wt[16][16]={ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
		-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0, 
		2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 
		0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
		0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 
		0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 
		-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0, 
		9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2, 
		-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2, 
		2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 
		-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1, 
		4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};	//please check!!! Traspose needed!
	
	for(i=0; i<4; i++)	{
		tx[i  ]=ty[i];
		tx[i+4]=ty1[i]*gres1;
		tx[i+8]=ty2[i]*gres2;
		tx[i+12]=ty12[i]*gres1*gres2;
	}
	
	i_n = 0;
	for(i=0; i<4; i++)	{
		for(j=0; j<4; j++)	{
			xx = 0.0;
			
			for(k=0; k<16; k++)	{
				xx += wt[i_n][k]*tx[k];
			}
			i_n ++;
			tc[j][i] = xx;
		}
	}
	
	return;
}


void CMol::SavePdb(char szName[])
{
	int n;
	FILE *fOut;
	char szAtomName[16], szResName[16];

	fOut = fopen(szName, "w");
	fprintf(fOut, "REMARK PDB file generated by Lei's code.\n");
	for(n=0; n<nAtom; n++)	{
		strcpy(szAtomName, AtomName[n]);
		szAtomName[3] = 0;
		strcpy(szResName, ResName[n]);
		szResName[3] = 0;
		fprintf(fOut, "ATOM%7d  %-3s %-3s A%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n", 
			n+1, szAtomName, szResName, SegID[n], x[n], y[n], z[n]);
	}
	fprintf(fOut, "END\n");
	fclose(fOut);
}



int CMol::Query_Dihedral_Index(int ia, int ib, int ic, int id)
{
	int iPos, i;

	for(i=0; i<nDihedral; i++)	{
		iPos = 4*i;

		if( (ib==DihedralList[iPos+1]) && (ic==DihedralList[iPos+2]) )	{
			if((ia==DihedralList[iPos]) && (id==DihedralList[iPos+3]))	{
				return i;
			}
		}
		else if((ib==DihedralList[iPos+2]) && (ic==DihedralList[iPos+1]))	{
			if((ia==DihedralList[iPos+3]) && (id==DihedralList[iPos]))	{
				return i;
			}
		}
	}
	return (-1);
}


void CMol::Restrain_All_Torsions(int Flag)
{
	if(Flag == 0)	{	// turn off the restrain on forces
		Soft_Restrain_Dihedrals = 0;	// no restrain
	}
	else	{
		Cal_E_Dihedral(1);	// setup phi0 as the initial value of phi
		Soft_Restrain_Dihedrals = 1;
	}
}

void CMol::Identify_Bonds_With_H(void)
{
	int i, j, ia, ib, iPos, MaxBond=6;
	int *pBondCount, *pBondList, *pBondIdxList, iOffset_a, iOffset_b, H_Count;

	nRigidBond = 0;

	pBondCount = new int[nAtom];
	if(pBondCount == NULL)	{
		printf("pBondCount == NULL\n");
		exit(1);
	}
	pBondList = new int[nAtom*MaxBond];
	if(pBondList == NULL)	{
		printf("pBondList == NULL\n");
		exit(1);
	}
	pBondIdxList = new int[nAtom*MaxBond];
	if(pBondIdxList == NULL)	{
		printf("pBondIdxList == NULL\n");
		exit(1);
	}

	for(i=0; i<nAtom; i++)	{
		pBondCount[i] = 0;
	}

	for(i=0; i<nBond; i++)	{
		iPos = i * 2;
		ia = BondList[iPos  ];
		ib = BondList[iPos+1];

		iOffset_a = ia * MaxBond;
		iOffset_b = ib * MaxBond;

		pBondList[iOffset_a+pBondCount[ia]] = ib;
		pBondIdxList[iOffset_a+pBondCount[ia]] = i;
		pBondCount[ia]++;

		pBondList[iOffset_b+pBondCount[ib]] = ia;
		pBondIdxList[iOffset_b+pBondCount[ib]] = i;
		pBondCount[ib]++;
	}

	nHGrp = 0;
	for(i=0; i<nAtom; i++)	{
		H_Count = 0;
		if( mass[i] > 1.5 )	{	// a heavy atom, serving as a center of H group
			iOffset_a = i * MaxBond;
			for(j=0; j<pBondCount[i]; j++)	{
				ia = pBondList[iOffset_a+j];
				if( (mass[ia] < 1.1) && (mass[ia] > 0.9) )	{	// H atom
					if(H_Count == 0)	{	// a new H group
						HGrpFirst[nHGrp] = nRigidBond;
						HGrpSize[nHGrp] = 0;
						nHGrp++;
					}

					// record this rigid bond
					RigidBondList[nRigidBond][0] = i;
					RigidBondList[nRigidBond][1] = ia;
					Rigid_b0[nRigidBond] = Para_b0[pBondIdxList[iOffset_a+j]];
					Rigid_b0_SQ[nRigidBond] = Rigid_b0[nRigidBond]*Rigid_b0[nRigidBond];

					HGrpSize[nHGrp-1]++;
					nRigidBond++;
					H_Count++;
				}
			}
		}
	}

	// !!!!!!!!  H-H group is ignored !!


	delete []pBondCount;
	delete []pBondList;
	delete []pBondIdxList;
}


#define TOL_RATTLE	(1.0E-4)
#define TOL_RATTLE_SQ	(2.0E-8)
#define MAX_H_GRP	(8)
void CMol::Rattle(double dt)
{
	int Idx, iOffset, ia, ib, i, icnt, IdxBond;
	int done, consFailure, iter, maxiter=100;
	double pab_x, pab_y, pab_z, pabsq, diffsq, invdt;
	double refab_x[MAX_H_GRP], refab_y[MAX_H_GRP], refab_z[MAX_H_GRP], rpab;
	double rma, rmb, gab;
	double dp_x, dp_y, dp_z;
	double df_x, df_y, df_z;
	double x_Save[MAX_ATOM], y_Save[MAX_ATOM], z_Save[MAX_ATOM];
	double netdp_x[MAX_ATOM], netdp_y[MAX_ATOM], netdp_z[MAX_ATOM];
	int HGrp_List[MAX_H_GRP], HGrpSize_Local;

	if(Shake_On == 0)	{
		return;
	}

	if(dt != 0.0)	{
		invdt = 1.0/dt;
	}
	else	{
		invdt = 0.0;
	}

	memcpy(x_Save, x, sizeof(double)*nAtom);
	memcpy(y_Save, y, sizeof(double)*nAtom);
	memcpy(z_Save, z, sizeof(double)*nAtom);


	for(Idx=0; Idx<nHGrp; Idx++)	{
		icnt = HGrpSize[Idx];
		HGrpSize_Local = icnt + 1;
		iOffset = HGrpFirst[Idx];

		HGrp_List[0] = RigidBondList[iOffset][0];	// the heavy atom

		for(i=0; i<icnt; i++)	{
			HGrp_List[i+1] = RigidBondList[iOffset+i][1];	// H atom
		}

		for(i=0; i<HGrpSize_Local; i++)	{
			ia = HGrp_List[i];

			netdp_x[ia] = 0.0;
			netdp_y[ia] = 0.0;
			netdp_z[ia] = 0.0;

			x[ia] += vv[0][ia]*dt;
			y[ia] += vv[1][ia]*dt;
			z[ia] += vv[2][ia]*dt;
		}

		for(i=0; i<icnt; i++)	{
			ia = RigidBondList[iOffset+i][0];
			ib = RigidBondList[iOffset+i][1];

			refab_x[i] = x_Save[ia] - x_Save[ib];
			refab_y[i] = y_Save[ia] - y_Save[ib];
			refab_z[i] = z_Save[ia] - z_Save[ib];
		}

		for(iter=0; iter<maxiter; ++iter)	{
			done = 1;
			consFailure = 0;

			for(i=0; i<icnt; i++)	{
				IdxBond = iOffset+i;
				ia = RigidBondList[IdxBond][0];
				ib = RigidBondList[IdxBond][1];
				
				pab_x = x[ia] - x[ib];
				pab_y = y[ia] - y[ib];
				pab_z = z[ia] - z[ib];
				pabsq = pab_x*pab_x + pab_y*pab_y  + pab_z*pab_z;
				diffsq = Rigid_b0_SQ[IdxBond] - pabsq;
				
				if( fabs(diffsq) > (Rigid_b0_SQ[IdxBond] * TOL_RATTLE_SQ) )	{
					rpab = refab_x[i]*pab_x + refab_y[i]*pab_y + refab_z[i]*pab_z;
					if ( rpab < ( Rigid_b0_SQ[IdxBond] * 1.0e-6 ) ) {
						done = 0;
						consFailure = 1;
						printf("Constraint failure in rattle algorithm for atom %d and %d", ia, ib);
						continue;
					}
					rma = Inv_mass[ia];
					rmb = Inv_mass[ib];
					gab = diffsq / ( 2.0 * ( rma + rmb ) * rpab );
					dp_x = refab_x[i] * gab;
					dp_y = refab_y[i] * gab;
					dp_z = refab_z[i] * gab;
					
					x[ia] += rma * dp_x;
					y[ia] += rma * dp_y;
					z[ia] += rma * dp_z;
					
					x[ib] -= rmb * dp_x;
					y[ib] -= rmb * dp_y;
					z[ib] -= rmb * dp_z;
					
					if ( invdt != 0. ) {
						dp_x *= invdt;
						dp_y *= invdt;
						dp_z *= invdt;
						
						netdp_x[ia] += dp_x;
						netdp_y[ia] += dp_y;
						netdp_z[ia] += dp_z;
						
						netdp_x[ib] -= dp_x;
						netdp_y[ib] -= dp_y;
						netdp_z[ib] -= dp_z;
					}
					done = 0;
				}
				
				if ( done ) break;
			}
		}

		if(invdt != 0.)	{
			for(i=0; i<HGrpSize_Local; i++)	{
				ia = HGrp_List[i];
				
				vv[0][ia] += netdp_x[ia]*Inv_mass[ia];
				vv[1][ia] += netdp_y[ia]*Inv_mass[ia];
				vv[2][ia] += netdp_z[ia]*Inv_mass[ia];
				
				df_x = -netdp_x[ia] * invdt;
				df_y = -netdp_y[ia] * invdt;
				df_z = -netdp_z[ia] * invdt;
				
				grad_x[ia] += df_x;
				grad_y[ia] += df_y;
				grad_z[ia] += df_z;
				
				x[ia] = x_Save[ia];
				y[ia] = y_Save[ia];
				z[ia] = z_Save[ia];
			}
		}
	}
}

/*
void CMol::Init_LangevinDynamics(double Temperature, double Time_Step, int Shake)
{
	int i;
	double fact, alpha;

	if(Shake)	{
		Shake_On = 1;
	}
	else	{
		Shake_On = 0;
	}

	dT_ps = Time_Step;
	Kelvin = Temperature;
	DeltaT=dT_ps/TIMFAC;	// 0.5 fs
//	DeltaT=0.01/TIMFAC;	// 0.5 fs
	DELTA2=DeltaT*DeltaT;
	DELTAS = 0.5*DELTA2;
//	iseed = time(NULL);
//	iseed = -67;
	kbt = Kelvin*KBOLTZ;

	NATOM2=nAtom+nAtom;
	NATOM3=NATOM2+nAtom;


	
	for(i=0; i<nAtom; i++)	{
		xcomp[i] = x[i];
		ycomp[i] = y[i];
		zcomp[i] = z[i];
	}
	
	myrand = new Random(iseed);
	myrand->split(1, 2);
	
	Init_Velocity();
	Cal_E(0);
//	Lngfill();
//	DlnGV();
	
	for(i=0; i<nAtom; i++)	{
		xold[i] = 0.0;
		yold[i] = 0.0;
		zold[i] = 0.0;
		xnew[i] = 0.0;
		ynew[i] = 0.0;
		znew[i] = 0.0;
	}
	
	for(i=0; i<nAtom; i++)	{
		fact = DELTAS/mass[i];
		alpha = 2.0*GAMMA[i+NATOM2]*GAMMA[i+NATOM3]*DELTA2;
		xold[i] = vv[0][i]*alpha-grad_x[i]*fact;
		yold[i] = vv[1][i]*alpha-grad_y[i]*fact;
		zold[i] = vv[2][i]*alpha-grad_z[i]*fact;
		alpha = alpha/GAMMA[i+NATOM2];
		xnew[i] = vv[0][i]*alpha+grad_x[i]*fact;
		ynew[i] = vv[1][i]*alpha+grad_y[i]*fact;
		znew[i] = vv[2][i]*alpha+grad_z[i]*fact;
	}
	
//	v_c1 = exp(-1.0*dt);
	v_c1 = exp(-FBETA*dT_ps);
	v_c2 = sqrt(1.0-v_c1*v_c1);

	for(i=0; i<nAtom; i++)	{
		if(mass[i] != 0.0)	{
//			vmean[i] = sqrt(kbt/mass[i]);
			vmean[i] = sqrt(2.0 * dT_ps * FBETA * kbt / mass[i]);
		}
		else	{
			vmean[i] = 0.0;
		}
	}


	return;
}
*/

void CMol::Transfer_Drude_Mass(void)
{
	int i;

	if(nDrude == 0)	{
		return;
	}

	for(i=0; i<nAtom; i++)	{
		if(IsDrude[i])	{
			mass[i-1] = mass[i-1] + mass[i];
			mass[i] = 0.0;

			Inv_mass[i-1] = 1.0 / mass[i-1];
			Inv_mass[i  ] = 0.0;	// not used anyway
		}
	}
}

void CMol::Optimize_Drude_Position(void)
{
	int i;

	Geo_Opt_Drude_Only = 1;
	FullGeometryOptimization_LBFGS(1);
	Geo_Opt_Drude_Only = 0;

	Cal_E(0);

	for(i=0; i<nAtom; i++)	{	// assign the same velocity for Drude and nucleus
		if(IsDrude[i])	{
			vv[0][i] = vv[0][i-1];
			vv[1][i] = vv[1][i-1];
			vv[2][i] = vv[2][i-1];
		}
	}
}

void CMol::Pre_NVT_Velocity_Verlet(void)
{
	int i;
	double DeltaT_SQ, dt_gamma, scale_v, scale2_v;

	DeltaT_SQ = DeltaT * DeltaT;
	dt_gamma = DeltaT * FBETA;
	scale_v = 1.0 - 0.5 * dt_gamma;
	scale2_v = 1.0 + 0.5 * dt_gamma;

	Rattle(0.0);

	Cal_E(0);

	for(i=0; i<nAtom; i++)	{	// update velocity with -0.5 * dT
		if(IsDrude[i])	{
			vv[0][i] += (  0.5*DeltaT*grad_x[i-1]*Inv_mass[i-1] );
			vv[1][i] += (  0.5*DeltaT*grad_y[i-1]*Inv_mass[i-1] );
			vv[2][i] += (  0.5*DeltaT*grad_z[i-1]*Inv_mass[i-1] );
		}
		else	{
			vv[0][i] += (  0.5*DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += (  0.5*DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += (  0.5*DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}

	Rattle(-DeltaT);

	for(i=0; i<nAtom; i++)	{	// update forces with 1.0 * dT
		if(IsDrude[i])	{
			vv[0][i] += (  -DeltaT*grad_x[i-1]*Inv_mass[i-1] );
			vv[1][i] += (  -DeltaT*grad_y[i-1]*Inv_mass[i-1] );
			vv[2][i] += (  -DeltaT*grad_z[i-1]*Inv_mass[i-1] );
		}
		else	{
			vv[0][i] += (  -DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += (  -DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += (  -DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}

	Rattle(DeltaT);

	for(i=0; i<nAtom; i++)	{	// update velocity with -0.5 * dT
		if(IsDrude[i])	{
			vv[0][i] += (  0.5*DeltaT*grad_x[i-1]*Inv_mass[i-1] );
			vv[1][i] += (  0.5*DeltaT*grad_y[i-1]*Inv_mass[i-1] );
			vv[2][i] += (  0.5*DeltaT*grad_z[i-1]*Inv_mass[i-1] );
		}
		else	{
			vv[0][i] += (  0.5*DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += (  0.5*DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += (  0.5*DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}
}

/*
void CMol::NVT_Velocity_Verlet(void)
{
	int i, j;
	double DeltaT_SQ, vmean_i, dt_gamma, scale_v, scale2_v;

	DeltaT_SQ = DeltaT * DeltaT;
	dt_gamma = dT_ps * FBETA;
	scale_v = 1.0 - 0.5 * dt_gamma;
	scale2_v = 1.0 + 0.5 * dt_gamma;

	for(i=0; i<nAtom; i++)	{	// update velocity with 0.5 * dT
		if(IsDrude[i])	{
//			vv[0][i] += ( - 0.5*DeltaT*grad_x[i-1]*Inv_mass[i-1] );
//			vv[1][i] += ( - 0.5*DeltaT*grad_y[i-1]*Inv_mass[i-1] );
//			vv[2][i] += ( - 0.5*DeltaT*grad_z[i-1]*Inv_mass[i-1] );
			vv[0][i] = vv[0][i-1];
			vv[1][i] = vv[1][i-1];
			vv[2][i] = vv[2][i-1];
		}
		else	{
			vv[0][i] += ( - 0.5*DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += ( - 0.5*DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += ( - 0.5*DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}

	for(i=0; i<nAtom; i++)	{
		x[i] += ( vv[0][i]*DeltaT);
		y[i] += ( vv[1][i]*DeltaT);
		z[i] += ( vv[2][i]*DeltaT);
	}

	Cal_E(0);	// update forces
	
	for(i=0; i<nAtom; i++)	{	// langevinVelocitiesBBK1()
        vv[0][i] *= scale_v;
        vv[1][i] *= scale_v;
        vv[2][i] *= scale_v;
	}

	for(i=0; i<nAtom; i++)	{	// update forces with 1.0 * dT
		if(IsDrude[i])	{
//			vv[0][i] += (  -DeltaT*grad_x[i-1]*Inv_mass[i-1] );
//			vv[1][i] += (  -DeltaT*grad_y[i-1]*Inv_mass[i-1] );
//			vv[2][i] += (  -DeltaT*grad_z[i-1]*Inv_mass[i-1] );
			vv[0][i] = vv[0][i-1];
			vv[1][i] = vv[1][i-1];
			vv[2][i] = vv[2][i-1];
		}
		else	{
			vv[0][i] += (  -DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += (  -DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += (  -DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}
	
	Rattle(DeltaT);
	for(i=0; i<nAtom; i++)	{	// langevinVelocitiesBBK2()
		if(IsDrude[i])	{
			continue;
		}
		vmean_i = vmean[i];
		vv[2][i] += ( myrand->gaussian()*vmean_i );
		vv[1][i] += ( myrand->gaussian()*vmean_i );
		vv[0][i] += ( myrand->gaussian()*vmean_i );
		
        vv[0][i] /= scale2_v;
        vv[1][i] /= scale2_v;
        vv[2][i] /= scale2_v;
	}

//	for(i=0; i<nAtom; i++)	{
//		vmean_i = vmean[i];
//		vv[0][i] = v_c1*vv[0][i] + v_c2*grand_charmm(vmean_i);
//		vv[1][i] = v_c1*vv[1][i] + v_c2*grand_charmm(vmean_i);
//		vv[2][i] = v_c1*vv[2][i] + v_c2*grand_charmm(vmean_i);
//	}
	
	Rattle(DeltaT);
	
	if(nDrude > 0)	Optimize_Drude_Position();
	
	for(i=0; i<nAtom; i++)	{	// update forces with -0.5 * dT
		if(IsDrude[i])	{
			vv[0][i] += (   0.5*DeltaT*grad_x[i-1]*Inv_mass[i-1] );
			vv[1][i] += (   0.5*DeltaT*grad_y[i-1]*Inv_mass[i-1] );
			vv[2][i] += (   0.5*DeltaT*grad_z[i-1]*Inv_mass[i-1] );
		}
		else	{
			vv[0][i] += (   0.5*DeltaT*grad_x[i]*Inv_mass[i] );
			vv[1][i] += (   0.5*DeltaT*grad_y[i]*Inv_mass[i] );
			vv[2][i] += (   0.5*DeltaT*grad_z[i]*Inv_mass[i] );
		}
	}
	

	E_Kinetic = 0.0;
	for(i=0; i<nAtom; i++)	{
		for(j=0; j<=2; j++)	{
			E_Kinetic += (0.5*mass[i]*vv[j][i]*vv[j][i]);
		}
	}
}
*/

void CMol::NVE_Velocity_Verlet(void)
{
	int i, j;
	double DeltaT_SQ;

	DeltaT_SQ = DeltaT * DeltaT;

	for(i=0; i<nAtom; i++)	{
		x[i] += ( vv[0][i]*DeltaT - 0.5*DeltaT_SQ*grad_x[i]*Inv_mass[i] );
		y[i] += ( vv[1][i]*DeltaT - 0.5*DeltaT_SQ*grad_y[i]*Inv_mass[i] );
		z[i] += ( vv[2][i]*DeltaT - 0.5*DeltaT_SQ*grad_z[i]*Inv_mass[i] );
	}

	for(i=0; i<nAtom; i++)	{
		vv[0][i] += ( - 0.5*DeltaT*grad_x[i]*Inv_mass[i] );
		vv[1][i] += ( - 0.5*DeltaT*grad_y[i]*Inv_mass[i] );
		vv[2][i] += ( - 0.5*DeltaT*grad_z[i]*Inv_mass[i] );
	}

//	Rattle(DeltaT);

	Cal_E(0);

	E_Kinetic = 0.0;
	for(i=0; i<nAtom; i++)	{
		vv[0][i] += (  -0.5*DeltaT*grad_x[i]*Inv_mass[i] );
		vv[1][i] += (  -0.5*DeltaT*grad_y[i]*Inv_mass[i] );
		vv[2][i] += (  -0.5*DeltaT*grad_z[i]*Inv_mass[i] );

		for(j=0; j<=2; j++)	{
			E_Kinetic += (0.5*mass[i]*vv[j][i]*vv[j][i]);
		}

	}

}

void CMol::NVE_Velocity_Verlet_MTS(int LongStep)
{
	int i, j, Step;
	double DeltaT_SQ, LongDeltaT, LongDeltaT_SQ;

	DeltaT_SQ = DeltaT * DeltaT;
	LongDeltaT = DeltaT * LongStep;
	LongDeltaT_SQ = LongDeltaT * LongDeltaT;

	for(i=0; i<nAtom; i++)	{	// long
		vv[0][i] += ( 0.5*LongDeltaT*Nonbonded_fx[i]*Inv_mass[i] );
		vv[1][i] += ( 0.5*LongDeltaT*Nonbonded_fy[i]*Inv_mass[i] );
		vv[2][i] += ( 0.5*LongDeltaT*Nonbonded_fz[i]*Inv_mass[i] );
	}
//	Cal_Force_Short();	// swtich to force short

	for(Step=1; Step<=LongStep; Step++)	{	// force short
		for(i=0; i<nAtom; i++)	{
			vv[0][i] += ( 0.5*DeltaT*Bonded_fx[i]*Inv_mass[i] );
			vv[1][i] += ( 0.5*DeltaT*Bonded_fy[i]*Inv_mass[i] );
			vv[2][i] += ( 0.5*DeltaT*Bonded_fz[i]*Inv_mass[i] );
		}

		for(i=0; i<nAtom; i++)	{
			x[i] += ( vv[0][i]*DeltaT);
			y[i] += ( vv[1][i]*DeltaT);
			z[i] += ( vv[2][i]*DeltaT);
		}

		Rattle(DeltaT);
		Cal_Force_Short();

		for(i=0; i<nAtom; i++)	{
			vv[0][i] += ( 0.5*DeltaT*Bonded_fx[i]*Inv_mass[i] );
			vv[1][i] += ( 0.5*DeltaT*Bonded_fy[i]*Inv_mass[i] );
			vv[2][i] += ( 0.5*DeltaT*Bonded_fz[i]*Inv_mass[i] );
		}
	}

	Cal_Force_Long();

	E_Kinetic = 0.0;
	for(i=0; i<nAtom; i++)	{
		vv[0][i] += ( 0.5*LongDeltaT*Nonbonded_fx[i]*Inv_mass[i] );
		vv[1][i] += ( 0.5*LongDeltaT*Nonbonded_fy[i]*Inv_mass[i] );
		vv[2][i] += ( 0.5*LongDeltaT*Nonbonded_fz[i]*Inv_mass[i] );

		for(j=0; j<=2; j++)	{
			E_Kinetic += (0.5*mass[i]*vv[j][i]*vv[j][i]);
		}
	}
}


void CMol::LangevinDynamics(int nSteps)
{
	int i;
	double fact, alpha;
	
	for(i=0; i<nAtom; i++)	{
		x[i] = xcomp[i] + xold[i];
		y[i] = ycomp[i] + yold[i];
		z[i] = zcomp[i] + zold[i];
	}

	Cal_E(0);

	for(i=0; i<nAtom; i++)	{
		xcomp[i] = x[i];
		ycomp[i] = y[i];
		zcomp[i] = z[i];
	}
	
	DlnGV();
	
	
	for(i=0; i<nAtom; i++)	{
		fact=GAMMA[i+nAtom];
		alpha=GAMMA[i+NATOM2];
		xnew[i]=alpha*xold[i]-grad_x[i]*fact;
		ynew[i]=alpha*yold[i]-grad_y[i]*fact;
		znew[i]=alpha*zold[i]-grad_z[i]*fact;
	}
	
	for(i=0; i<nAtom; i++)	{
		fact=GAMMA[i+NATOM3];
		vv[0][i] = (xnew[i] + xold[i])*fact;
		vv[1][i] = (ynew[i] + yold[i])*fact;
		vv[2][i] = (znew[i] + zold[i])*fact;
	}	
	
	for(i=0; i<nAtom; i++)	{
		fact=xold[i];
		xold[i]=xnew[i];
		xnew[i]=fact;
		fact=yold[i];
		yold[i]=ynew[i];
		ynew[i]=fact;
		fact=zold[i];
		zold[i]=znew[i];
		znew[i]=fact;
	}
		
	return;
}


void CMol::DlnGV(void)
{
	int i, k;
	double A, B, RDX, RDY, RDZ;
	
	k = 0;
	
	for(i=0; i<nAtom; i++)	{
		if(k==0)	{
			A=GAMMA[i]*sqrt(-2.0*log(rand_charmm(iseed)));
			B=PI2*rand_charmm(iseed);
			k=1;
			RDX=A*cos(B);
			RDY=A*sin(B);
			grad_x[i] += RDX;
			grad_y[i] += RDY;
			A=sqrt(-2.0*log(rand_charmm(iseed)));
			B=PI2*rand_charmm(iseed);
			RDZ=GAMMA[i]*A*cos(B);
			grad_z[i] += RDZ;
		}
		else	{
			k=0;
			RDX=GAMMA[i]*A*sin(B);
			grad_x[i] += RDX;
			A=GAMMA[i]*sqrt(-2.0*log(rand_charmm(iseed)));
			B=PI2*rand_charmm(iseed);
			RDY=A*cos(B);
			RDZ=A*sin(B);
			grad_y[i] += RDY;
			grad_z[i] += RDZ;
		}
	}
}

#define DIVIS	(2147483647.0)
#define DENOM	(2147483711.0)
#define MULTIP	(16807.0)

double CMol::rand_charmm(int &iseed)
{
	double dseed, r;
	
	if(iseed <= 1)	{
		iseed = 314159;
	}
	dseed=MULTIP*iseed;
	
	dseed=fmod(dseed,DIVIS);
	r = dseed/DENOM;
	iseed=(int)dseed;
	
	return r;
}

double CMol::grand_charmm(double Variance)	//Variance indicate T
{
	double A, B, r;
	
	A=sqrt(-2.0*log(rand_charmm(iseed)));
	B=PI2*rand_charmm(iseed);
	r=Variance*A*cos(B);
	
	return r;
}


void CMol::Lngfill(void)
{
	int i;
	double GAM;

	for(i=0; i<nAtom; i++)	{
		GAM=TIMFAC*FBETA*DeltaT;
		GAMMA[i]=sqrt(2.0*mass[i]*GAM*kbt)/DeltaT;
		GAMMA[i+nAtom]=DeltaT*DeltaT/((1.0+GAM*0.5)*mass[i]);
		GAMMA[i+NATOM2]=(1.0-GAM*0.5)/(1.0+GAM*0.5);
		GAMMA[i+NATOM3]=0.5*sqrt(1.0+GAM*0.5)/DeltaT;
	}
	
	return;
}

/*
void CMol::Init_Velocity(void)
{
	int i, j;
	double SD_Vel, randnum, vx_Sum=0.0, vy_Sum=0.0, vz_Sum=0.0, mass_Sum=0.0;
	Random random(iseed);

//	
//	E_Kinetic = 0.0;
//	for(j=0; j<nAtom; j++)	{
//		SD_Vel = sqrt(KBOLTZ * Kelvin / mass[j]);
//		for(i=0; i<=2; i++)	{
//			vv[i][j]=grand_charmm(SD_Vel);
//			E_Kinetic += (0.5*mass[j]*vv[i][j]*vv[i][j]);
//		}
//	}
//

  for(i=0; i<nAtom; i++)
  {
	SD_Vel = sqrt(KBOLTZ * Kelvin / mass[i]);

    //  The following comment was stolen from X-PLOR where
    //  the following section of code was adapted from.
    
    //  This section generates a Gaussian random
    //  deviate of 0.0 mean and standard deviation RFD for
    //  each of the three spatial dimensions.
    //  The algorithm is a "sum of uniform deviates algorithm"
    //  which may be found in Abramowitz and Stegun,
    //  "Handbook of Mathematical Functions", pg 952.
    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += random.uniform();
    }

    randnum -= 6.0;

    vv[0][i] = randnum*SD_Vel;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += random.uniform();
    }

    randnum -= 6.0;

    vv[1][i] = randnum*SD_Vel;

    for (randnum=0.0, j=0; j<12; j++)
    {
      randnum += random.uniform();
    }

    randnum -= 6.0;
    
    vv[2][i] = randnum*SD_Vel;
  }


  for(i=0; i<nAtom; i++)	{
	  mass_Sum += mass[i];
	  vx_Sum += (vv[0][i] * mass[i]);
	  vy_Sum += (vv[1][i] * mass[i]);
	  vz_Sum += (vv[2][i] * mass[i]);
  }
  vx_Sum /= mass_Sum;
  vy_Sum /= mass_Sum;
  vz_Sum /= mass_Sum;

  for(i=0; i<nAtom; i++)	{
	  vv[0][i] -= vx_Sum;
	  vv[1][i] -= vy_Sum;
	  vv[2][i] -= vz_Sum;

	  if(IsDrude[i])	{	// alway same as the host atom
		  vv[0][i] = vv[0][i-1];
		  vv[1][i] = vv[1][i-1];
		  vv[2][i] = vv[2][i-1];
	  }
  }


}
*/

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) Quit_With_Error_Msg("allocation failure in dvector()");
	return v-nl+1;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
	if (!m) Quit_With_Error_Msg("allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
	if (!m[nrl]) Quit_With_Error_Msg("allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	return m;
}

void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl-1));
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	free((char*) (m[nrl]+ncl-1));
	free((char*) (m+nrl-1));
}

