#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nlopt.h"

#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>	//for linux
#endif

#define MAX_VDW		(512)
#define MAX_VAR		(MAX_VDW*2)
#define K_BOTZMAN	(0.001987191)
#define MAX_N_MOL	(512)

#ifndef MINAB
#define MINAB 
#define	min(a,b)	((a<b)?(a):(b))
#define	max(a,b)	((a>b)?(a):(b))
#endif

double func_Obj_Grad(unsigned n, const double *x, double *grad, void *my_func_data);
double Cal_Objective_Function_Grad(void);
void Prepare_New_Force_Field(int ID);
void Submit_All_Jobs(double *grad);

double Get_Largest_Grad(double *grad);

int count = 0, n_Para=0, nVDW_Type=0, JobIdx;
nlopt_opt opt;
double x[MAX_VAR];
double f_Cur;

double lb[MAX_VAR];
double ub[MAX_VAR];
char szLog[]="lj-log.txt";
char szName_JobScript[256]="lj-sim-C3R-HC3R.py";

//double w_qm=2.0, w_v=1000.0, w_hvap=100.0, w_sfe=1.0;	// the weights for three terms
//double Corr, Volume, Hvap, dG_SFE;
//double Volume_0[MAX_N_MOL];
//double Hvap_0[MAX_N_MOL];	// target value for ethane
double Grad_List[MAX_VAR];
char ChemName[MAX_VAR][16];

void Read_VDW_Type(void);

void Read_VDW_Type(void)
{
	FILE *fIn;
	char szName[]="vdw-to-fit.txt", szLine[256], *ReadLine, szType[256];
	double para[8];
	int ReadItem;

	fIn = fopen(szName, "r");

	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit.\n", szName);
		exit(1);
	}

	nVDW_Type = n_Para = 0;
	fgets(szLine, 256, fIn);
	
	while(1)	{
		if(feof(fIn))	{
			break;
		}
		ReadLine = fgets(szLine, 256, fIn);
		if(ReadLine)	{
			ReadItem = sscanf(szLine, "%s%lf%lf%lf%lf%lf%lf", szType, &(para[0]), &(para[1]), &(para[2]), &(para[3]), &(para[4]), &(para[5]));
			if( (ReadItem == 7) && (szType[0] != '!'))	{
				strcpy(ChemName[nVDW_Type], szType);
				x[n_Para  ] = para[0];	// Emin
				x[n_Para+1] = para[3];	// Rmin
				lb[n_Para  ] = para[1];
				ub[n_Para  ] = para[2];
				lb[n_Para+1] = para[4];
				ub[n_Para+1] = para[5];
				nVDW_Type++;
				n_Para+=2;
			}
		}
	}

	fclose(fIn);
}

int main(int argc, char *argv[])
{
	FILE *fOut;
	//double minf=1.0E100; /* the minimum objective value, upon return */
        double minf;
	if(argc == 2)	{
		strcpy(szName_JobScript, argv[1]);
	}

//	Cal_Objective_Function_Grad();

	Read_VDW_Type();

	opt = nlopt_create(NLOPT_LD_LBFGS, n_Para); /* algorithm and dimensionality */
        int i;
        for (i=0;i<n_Para;i++) {
          printf("Param: %f ; lb: %f; ub: %f\n",x[i],lb[i],ub[i]);
        }

	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, ub);

	nlopt_set_min_objective(opt, func_Obj_Grad, NULL);
	nlopt_set_maxeval(opt,30);
	nlopt_set_xtol_rel(opt, 1e-3); //new tolerence
	//nlopt_set_ftol_rel(opt,1e-3); //tolerance on objective function

	
	fOut = fopen(szLog, "w");	// to create an empty log file
	fprintf(fOut, "There are %3d vdw type (%3d parameters) to be fitted.\n", nVDW_Type, n_Para);
	fclose(fOut);
  
        int check;
        check = nlopt_optimize(opt, x, &minf);
	if (check < 0) {
		printf("nlopt failed! %d\n",check);
	}
	else {
		printf("found minimum after %d evaluations\n", count);
		printf("Objective function = %0.10g\n", minf);
		printf("success return code = %d\n",check);
	}

	nlopt_destroy(opt);

	return 0;
}


double func_Obj_Grad(unsigned n, const double *x, double *grad, void *my_func_data)
{
	FILE *fOut;

	count++;

	Prepare_New_Force_Field(count);
	Submit_All_Jobs(grad);	// submit jobs and wait until all jobs are done
	f_Cur = Cal_Objective_Function_Grad();

	memcpy(grad, Grad_List, sizeof(double)*n);

	fOut = fopen(szLog, "a+");
	fseek(fOut, SEEK_END, 0);
	if(grad)	{
		fprintf(fOut, "%5d %.5lf  %.5lf\n", count, f_Cur, Get_Largest_Grad(grad));
	}
	else	{
		fprintf(fOut, "%5d %.5lf\n", count, f_Cur);
	}
	fclose(fOut);

    return f_Cur;	// objective function
}

double Get_Largest_Grad(double *grad)
{
	double f_Max=0.0, f_abs;

	if(grad)	{
		for(int i=0; i<n_Para; i++)	{
			f_abs = fabs(grad[i]);

			if(f_abs > f_Max)	{
				f_Max = f_abs;
			}
		}

		return f_Max;
	}
	else	{
		return 100.0;
	}
}

void Prepare_New_Force_Field(int ID)
{
	FILE *fOut;
	char szDirName[256], szCmd[256], szName[256];
	int i;

	sprintf(szDirName, "run_%d", count);
	sprintf(szCmd, "mkdir -p %s", szDirName);
	system(szCmd);

	sprintf(szName, "%s/vdw-param.txt", szDirName);
	fOut = fopen(szName, "w");

	for(i=0; i<n_Para; i+=2)	{
//		fprintf(fOut, "%-8s 0.00   %14.9lf  %14.9lf    0.00   %14.9lf  %14.9lf\n", 
//			ChemName[i/2], -x[i], x[i+1], -x[i]*0.5, x[i+1]);	// Emin, Rmin, Emin_14, Rmin_14
		fprintf(fOut, "%-8s  %14.9lf  %14.9lf\n", 
			ChemName[i/2], -x[i], x[i+1]);	// Emin, Rmin, Emin_14, Rmin_14
	}

	fclose(fOut);

	return;	
}

void Submit_All_Jobs(double *grad)
{
	char szDirName[256], szCurDir[256], szCmd[256];
	FILE *fIn;

	getcwd(szCurDir, 256);
	
	sprintf(szDirName, "%s/run_%d", szCurDir, count);
	chdir(szDirName);
	
	fIn = fopen("done.txt", "r");
	if(fIn == NULL)	{
		sprintf(szCmd, "python ../%s > log-job.txt 2>&1 &", szName_JobScript);
		system(szCmd);
		
		//start	to wait for all jobs finish
		while(1)	{
			fIn = fopen("done.txt", "r");
			if(fIn == NULL)	{
				system("sleep 15");
			}
			else	{
				fclose(fIn);
				break;
			}
		}
		//end	to wait for all jobs finish
	}
	else	{
		fclose(fIn);
	}
	

	chdir(szCurDir);
}


double Cal_Objective_Function_Grad(void)	// the objective function and gradient will be calculated in python script
{
	FILE *fIn;
	int i, nCount, iPos;
	int ReadItem;
	double f, grad[2];
	char szDirName[256], szCurDir[256], szName[]="obj_func.txt", szName_Grad[]="grad_list.txt", szBuff[256];

	getcwd(szCurDir, 256);

	sprintf(szDirName, "%s/run_%d", szCurDir, count);
	chdir(szDirName);

	//start	reading objective function
	fIn = fopen(szName, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szName);
		exit(1);
	}

	ReadItem = fscanf(fIn, "%lf", &f);

	if(ReadItem != 1)	{
		printf("Error in reading %s\nQuit\n", szName);
		fclose(fIn);
		exit(1);
	}

	fclose(fIn);
	//end	reading objective function

	//start	reading objective function
	fIn = fopen(szName_Grad, "r");
	if(fIn == NULL)	{
		printf("Fail to open file: %s\nQuit\n", szName_Grad);
		exit(1);
	}

	nCount = 0;
	iPos = 0;

	for(i=0; i<nVDW_Type; i++, iPos+=2)	{
		ReadItem = fscanf(fIn, "%s%lf%lf", szBuff, &(grad[0]), &(grad[1]));
		if(ReadItem != 3)	{
			printf("Error in reading file: %s\nQuit\n", szName_Grad);
			fclose(fIn);
			exit(1);
		}
		if(strcmp(ChemName[i], szBuff)!=0)	{
			printf("Error in reading gradient from file: %s\nQuit\n", szName_Grad);
			fclose(fIn);
			exit(1);
		}
		Grad_List[iPos  ] = grad[0];
		Grad_List[iPos+1] = grad[1];
	}

	fclose(fIn);
	//end	reading objective function


	chdir(szCurDir);

	
	return f;
}
