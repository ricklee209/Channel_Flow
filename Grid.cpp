






#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define PI 3.141592653589793

#include "Resolution.h"

extern int X_np;

void Grid
(
// ======================================================== //
int myid,

double (*X_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*Y_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*Z_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ======================================================== //
)

{

#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

	double temp;

double gamma = 2.8;

// ======================= //
	istart = 0;            //
	iend = gend[myid]+3;   //
// ======================= //
	

//#pragma omp parallel for private(j,k,temp)
		for (i = istart; i <= iend; i++) {
			for (j = 2; j < nyy; j++) { 
				for (k = 2; k < nzz; k++) {

				Y_point[i][j][k] = (0.5*high)*(1-1./tanh(gamma)*tanh(gamma*(1-2*(j-1.5)*deltaET)));
				
			}
		}
	}  


if (myid == 0) {
	for (j = 2; j <= ny; j++) {
		for (k = 2; k <= nz; k++) {   

			istart = 2;

			Y_point[istart][j][k]= Y_point[istart+1][j][k];
			
		}
	}
}

if (myid == np-1) {
	for (j = 2; j <= ny; j++) {
		for (k = 2; k <= nz; k++) {   

			iend = gend[myid];

			Y_point[iend+1][j][k]= Y_point[iend][j][k];

		}
	}
}


//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;    			  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {   

			Y_point[i][1][k]= -Y_point[i][2][k];
			Y_point[i][nyy][k]= Y_point[i][ny][k]+2*Y_point[i][2][k];


		}
	}


	
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;     			  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (j =1 ; j <= nyy; j++) {  

			Y_point[i][j][1]= Y_point[i][j][2];
			Y_point[i][j][nzz]= Y_point[i][j][nz];

		}
	}



//// ========================= ////
		 istart = 3;		   ////	
//// ========================= ////
		iend = gend[myid];	   ////
//// ========================= ////

#pragma omp parallel for private(j,k)

		for (i = istart; i <= iend; i++) {
			for (j = 2; j < nyy; j++) { 
				for (k = 2; k < nzz; k++) {

				X_point[i][j][k]=deltaXI*(i-1.5);
				
			}
		}
	}  

/**** output XYZ location ****/
// ============================================================================ //
	/*char data[100];
	FILE *fptr;
	sprintf(data,"PointXYZ.m");
	fptr = fopen(data,"w");

	if (fptr != NULL)  
	{
		for (k = 2; k <= nz; k++) {
			fprintf(fptr,"X_point(:,:,%d)=[\n",k-1);
			for (i = 2; i <= nx; i++) {
				for (j = 2; j <= ny; j++) {     
					fprintf(fptr,"%.16f\t",X_point[i][j][k]);
				}
				fprintf(fptr,";");
				fprintf(fptr,"\n");
			}
			fprintf(fptr,"];");
			fprintf(fptr,"\n");
		}



		for (k = 2; k <= nz; k++) {
			fprintf(fptr,"Y_point(:,:,%d)=[\n",k-1);
			for (i = 2; i <= nx; i++) {
				for (j = 2; j <= ny; j++) {     
					fprintf(fptr,"%.16f\t",Y_point[i][j][k]);
				}
				fprintf(fptr,";");
				fprintf(fptr,"\n");
			}
			fprintf(fptr,"];");
			fprintf(fptr,"\n");
		}

		fclose(fptr);
	}

	else
		printf("File opening Failure\n");*/

// ============================================================================ //
/**** output XYZ location end ****/


}