



#include <stdlib.h> 
#include <omp.h>
#include <stdio.h>
#include <mpi.h>

#include "Resolution.h"

extern int X_np;

void Initial_condition
(
// ======================================================== //
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*U1q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5q)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ======================================================== //
)
{
	
#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"


	double (*U1out)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U2out)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U3out)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U4out)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U5out)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];

	double (*U1out_old)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U2out_old)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U3out_old)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U4out_old)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];
	double (*U5out_old)[Y_out][Z_out] = new double[X_out][Y_out][Z_out];

	char data[100];

	FILE *fptr;

	sprintf(data,"initial.bin");
	fptr = fopen(data,"rb");

	fread(U1out,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U2out,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U3out,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U4out,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U5out,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	
	fread(U1out_old,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U2out_old,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U3out_old,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U4out_old,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);
	fread(U5out_old,sizeof(double),(nx-1)*(ny-1)*(nz-1),fptr);

	fclose(fptr);

	
	double rho, U, V, W, VV, P;


	int ii = 2;

	istart = gstart[myid];        
	iend = gend0[myid]; 

	for (i = istart; i <= iend; i++) {
		
				ii = ii+1;

#pragma omp parallel for private(k)

		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {  

				
				U1[ii][j][k] = U1out[i][j-2][k-2]/J[ii][j][k];
				U2[ii][j][k] = U2out[i][j-2][k-2]/J[ii][j][k];
				U3[ii][j][k] = U3out[i][j-2][k-2]/J[ii][j][k];
				U4[ii][j][k] = U4out[i][j-2][k-2]/J[ii][j][k];
				U5[ii][j][k] = U5out[i][j-2][k-2]/J[ii][j][k];

				U1q[ii][j][k] = U1out_old[i][j-2][k-2]/J[ii][j][k];
				U2q[ii][j][k] = U2out_old[i][j-2][k-2]/J[ii][j][k];
				U3q[ii][j][k] = U3out_old[i][j-2][k-2]/J[ii][j][k];
				U4q[ii][j][k] = U4out_old[i][j-2][k-2]/J[ii][j][k];
				U5q[ii][j][k] = U5out_old[i][j-2][k-2]/J[ii][j][k];
				
			}
		}
		
		
	}


	delete [] U1out;
	delete [] U2out;
	delete [] U3out;
	delete [] U4out;
	delete [] U5out;

	delete [] U1out_old;
	delete [] U2out_old;
	delete [] U3out_old;
	delete [] U4out_old;
	delete [] U5out_old;

}