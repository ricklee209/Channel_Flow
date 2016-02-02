




#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <omp.h>
#include "Resolution.h"

extern int X_np;

void Y_boundary_condition
(
// ============================================================================ //
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"
	
	double rho,U,V,W,VV,P,C,T,h,H;
	double temp;

//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;            			  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;				  ////
//// ============================================ ////

#pragma omp parallel for private(k,rho,U,V,W,VV,P,temp,T)
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
			
			/*
			U1_[i][1][k] = U1_[i][2][k]*J[i][2][k]/J[i][1][k];
			U2_[i][1][k] = -U2_[i][2][k]*J[i][2][k]/J[i][1][k];
			U3_[i][1][k] = -U3_[i][2][k]*J[i][2][k]/J[i][1][k];
			U4_[i][1][k] = -U4_[i][2][k]*J[i][2][k]/J[i][1][k];
			U5_[i][1][k] = U5_[i][2][k]*J[i][2][k]/J[i][1][k];

			U1_[i][0][k] = U1_[i][3][k]*J[i][3][k]/J[i][0][k];
			U2_[i][0][k] = -U2_[i][3][k]*J[i][3][k]/J[i][0][k];
			U3_[i][0][k] = -U3_[i][3][k]*J[i][3][k]/J[i][0][k];
			U4_[i][0][k] = -U4_[i][3][k]*J[i][3][k]/J[i][0][k];
			U5_[i][0][k] = U5_[i][3][k]*J[i][3][k]/J[i][0][k];


			U1_[i][nyy][k] = U1_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U2_[i][nyy][k] = -U2_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U3_[i][nyy][k] = -U3_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U4_[i][nyy][k] = -U4_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];
			U5_[i][nyy][k] = U5_[i][ny][k]*J[i][ny][k]/J[i][nyy][k];

			U1_[i][nyyy][k] = U1_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U2_[i][nyyy][k] = -U2_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U3_[i][nyyy][k] = -U3_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U4_[i][nyyy][k] = -U4_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			U5_[i][nyyy][k] = U5_[i][ny-1][k]*J[i][ny-1][k]/J[i][nyyy][k];
			*/

			rho = U1_[i][2][k]*J[i][2][k];
			U = U2_[i][2][k]/U1_[i][2][k];
			V = U3_[i][2][k]/U1_[i][2][k];
			W = U4_[i][2][k]/U1_[i][2][k];     
			VV = U*U+V*V+W*W;
			P = (U5_[i][2][k]*J[i][2][k]-0.5*rho*VV)*(K-1);
			temp = P/rho/R;

			T = 2*298.0592-temp;

			rho = P/R/T;

			U1_[i][1][k] = P/R/T/J[i][1][k];
			U2_[i][1][k] = -P/R/T*U/J[i][1][k];
			U3_[i][1][k] = -P/R/T*V/J[i][1][k];
			U4_[i][1][k] = -P/R/T*W/J[i][1][k];
			U5_[i][1][k] = (rho*R*T/(K-1)+0.5*P/R/T*VV)/J[i][1][k];




			rho = U1_[i][ny][k]*J[i][ny][k];
			U = U2_[i][ny][k]/U1_[i][ny][k];
			V = U3_[i][ny][k]/U1_[i][ny][k];
			W = U4_[i][ny][k]/U1_[i][ny][k];     
			VV = U*U+V*V+W*W;
			P = (U5_[i][ny][k]*J[i][ny][k]-0.5*rho*VV)*(K-1);
			temp = P/rho/R;

			T = 2*298.0592-temp;

			rho = P/R/T;

			U1_[i][nyy][k] = P/R/T/J[i][nyy][k];
			U2_[i][nyy][k] = -P/R/T*U/J[i][nyy][k];
			U3_[i][nyy][k] = -P/R/T*V/J[i][nyy][k];
			U4_[i][nyy][k] = -P/R/T*W/J[i][nyy][k];
			U5_[i][nyy][k] = (rho*R*T/(K-1)+0.5*P/R/T*VV)/J[i][nyy][k];
			
		}
	}
#pragma omp barrier

}
