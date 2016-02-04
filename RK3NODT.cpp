



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int X_np;

void RK3NODT
(
// ============================================================================ //
int myid,

int RK,

double deltaT,

double *er1,
double *er2,
double *er3,
double *er4,
double *er5,

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

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFx1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFz1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*vF2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*vF5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 


/**** MR = RHS ****/
double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
/**** MR = RHS-end ****/


// ============================================================================ //
)

{

	
#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

	double tmp;

	double Rp1,Rp2,Rp3,Rp4,Rp5,
		   Rf1,Rf2,Rf3,Rf4,Rf5,
		   Rk1,Rk2,Rk3,Rk4,Rk5;

	double c1[3] = { 1.0, 0.75, 1.0/3.0};
	double c2[3] = { 0.0, 0.25, 2.0/3.0};
	double c3[3] = { 1.0, 0.25, 2.0/3.0};



	

//// ============================================ ////
			istart = 3; 	            		  ////	
//// ============================================ ////
			iend = gend[myid];		    		  ////
//// ============================================ ////


	for (i = istart ; i <= iend; i++) {
		
#pragma omp parallel for private(\
k,_k,tmp,\
Rf1,Rf2,Rf3,Rf4,Rf5,Rk1,Rk2,Rk3,Rk4,Rk5\
)
		for (j = 2; j < nyy; j++) {
			for (k = 2, _k = 1; k < nzz; k++, _k++) {

				tmp = U2_[i][j][k]/U1_[i][j][k];

				Rf1 = -((inFx1[i][j-1][_k]-inFx1[i-1][j-1][_k])/deltaXI+
					(inFy1[i-1][j][_k]-inFy1[i-1][j-1][_k])/deltaET+
					(inFz1[i-1][j-1][k]-inFz1[i-1][j-1][_k])/deltaZT);
					
				Rf2 = -((inFx2[i][j-1][_k]-inFx2[i-1][j-1][_k])/deltaXI+
					(inFy2[i-1][j][_k]-inFy2[i-1][j-1][_k])/deltaET+
					(inFz2[i-1][j-1][k]-inFz2[i-1][j-1][_k])/deltaZT);
					
					
				Rf3 = -((inFx3[i][j-1][_k]-inFx3[i-1][j-1][_k])/deltaXI+
					(inFy3[i-1][j][_k]-inFy3[i-1][j-1][_k])/deltaET+
					(inFz3[i-1][j-1][k]-inFz3[i-1][j-1][_k])/deltaZT);
					
				Rf4 = -((inFx4[i][j-1][_k]-inFx4[i-1][j-1][_k])/deltaXI+
					(inFy4[i-1][j][_k]-inFy4[i-1][j-1][_k])/deltaET+
					(inFz4[i-1][j-1][k]-inFz4[i-1][j-1][_k])/deltaZT);
					
				Rf5 = -((inFx5[i][j-1][_k]-inFx5[i-1][j-1][_k])/deltaXI+
					(inFy5[i-1][j][_k]-inFy5[i-1][j-1][_k])/deltaET+
					(inFz5[i-1][j-1][k]-inFz5[i-1][j-1][_k])/deltaZT);

					
				Rk1 = Rf1;
				Rk2 = Rf2+vF2[i][j][k]+f/J[i][j][k];
				Rk3 = Rf3+vF3[i][j][k];
				Rk4 = Rf4+vF4[i][j][k];
				Rk5 = Rf5+vF5[i][j][k]+f*tmp/J[i][j][k];

				MR1[i][j][k] = deltaT*Rk1;
				MR2[i][j][k] = deltaT*Rk2;
				MR3[i][j][k] = deltaT*Rk3;
				MR4[i][j][k] = deltaT*Rk4;
				MR5[i][j][k] = deltaT*Rk5;
				

				U1_[i][j][k] = c1[RK-1]*U1[i][j][k] + c2[RK-1]*U1_[i][j][k] + c3[RK-1]*MR1[i][j][k];
				U2_[i][j][k] = c1[RK-1]*U2[i][j][k] + c2[RK-1]*U2_[i][j][k] + c3[RK-1]*MR2[i][j][k];
				U3_[i][j][k] = c1[RK-1]*U3[i][j][k] + c2[RK-1]*U3_[i][j][k] + c3[RK-1]*MR3[i][j][k];
				U4_[i][j][k] = c1[RK-1]*U4[i][j][k] + c2[RK-1]*U4_[i][j][k] + c3[RK-1]*MR4[i][j][k];
				U5_[i][j][k] = c1[RK-1]*U5[i][j][k] + c2[RK-1]*U5_[i][j][k] + c3[RK-1]*MR5[i][j][k];

			}
		}
	}
#pragma omp barrier



}