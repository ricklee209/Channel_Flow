#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "Resolution.h"

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

extern int X_np;

void Filter
(
// ============================================================================ //
int myid,

double Roe_criterion,
double E_high,
double E_low,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*Y_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpX)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*EpZ)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]

// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

int fX_np = gcount[myid]+6;

double (*G1p)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G2p)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G3p)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G4p)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G5p)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*G1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*G5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*GG1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GG2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GG3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GG4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GG5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*GGG1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGG2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGG3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGG4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGG5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*GGGG1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGG2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGG3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGG4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGG5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*GGGGG1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGGG2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGGG3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGGG4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*GGGGG5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*Gx1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gx2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gx3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gx4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gx5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*Gy1)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gy2)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gy3)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gy4)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*Gy5)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*fU1_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU2_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU3_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU4_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU5_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*fU1_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU2_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU3_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU4_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*fU5_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];

double (*E_2h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double (*E_h)[Y_m][Z_m] = new double[fX_np][Y_m][Z_m];
double Emax = 0;
double Ue, Ve, We;



double a,b,c;

// ======================= //
	istart = 0;            //
	iend = gend[myid]+3;   //
// ======================= //
	
for (i = istart; i <= iend; i++) {
		
#pragma omp parallel for private(k)

	for (j = 1; j <= nyy; j++) {
		for (k = 1; k <= nzz; k++) {

			G1[i][j][k] = GG1[i][j][k] = GGG1[i][j][k] = GGGG1[i][j][k] = GGGGG1[i][j][k] = Gx1[i][j][k] = Gy1[i][j][k] = G1p[i][j][k] =  fU1_h[i][j][k] = fU1_2h[i][j][k] = U1_[i][j][k];
			G2[i][j][k] = GG2[i][j][k] = GGG2[i][j][k] = GGGG2[i][j][k] = GGGGG2[i][j][k] = Gx2[i][j][k] = Gy2[i][j][k] = G2p[i][j][k] =  fU2_h[i][j][k] = fU2_2h[i][j][k] = U2_[i][j][k];
			G3[i][j][k] = GG3[i][j][k] = GGG3[i][j][k] = GGGG3[i][j][k] = GGGGG3[i][j][k] = Gx3[i][j][k] = Gy3[i][j][k] = G3p[i][j][k] =  fU3_h[i][j][k] = fU3_2h[i][j][k] = U3_[i][j][k];
			G4[i][j][k] = GG4[i][j][k] = GGG4[i][j][k] = GGGG4[i][j][k] = GGGGG4[i][j][k] = Gx4[i][j][k] = Gy4[i][j][k] = G4p[i][j][k] =  fU4_h[i][j][k] = fU4_2h[i][j][k] = U4_[i][j][k];
			G5[i][j][k] = GG5[i][j][k] = GGG5[i][j][k] = GGGG5[i][j][k] = GGGGG5[i][j][k] = Gx5[i][j][k] = Gy5[i][j][k] = G5p[i][j][k] =  fU5_h[i][j][k] = fU5_2h[i][j][k] = U5_[i][j][k];

		}
	}
}

#pragma omp barrier




/******* The primary filter *******/




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(U1_[i-1][j][k]+6*U1_[i][j][k]+U1_[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(U2_[i-1][j][k]+6*U2_[i][j][k]+U2_[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(U3_[i-1][j][k]+6*U3_[i][j][k]+U3_[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(U4_[i-1][j][k]+6*U4_[i][j][k]+U4_[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(U5_[i-1][j][k]+6*U5_[i][j][k]+U5_[i+1][j][k]);

			}
		}
	}

#pragma omp barrier




	
//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {


				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier





//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)

		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				G1p[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				G2p[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				G3p[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				G4p[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				G5p[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);

			}
		}
	}

#pragma omp barrier


/******* The primary filter end *******/




/******* First filter *******/


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	

	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(G1p[i-1][j][k]+6*G1p[i][j][k]+G1p[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(G2p[i-1][j][k]+6*G2p[i][j][k]+G2p[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(G3p[i-1][j][k]+6*G3p[i][j][k]+G3p[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(G4p[i-1][j][k]+6*G4p[i][j][k]+G4p[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(G5p[i-1][j][k]+6*G5p[i][j][k]+G5p[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier



//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				G1[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				G2[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				G3[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				G4[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				G5[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);

			}
		}
	}

#pragma omp barrier

/******* First filter end *******/




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
/******* Second filter *******/

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(G1[i-1][j][k]+6*G1[i][j][k]+G1[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(G2[i-1][j][k]+6*G2[i][j][k]+G2[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(G3[i-1][j][k]+6*G3[i][j][k]+G3[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(G4[i-1][j][k]+6*G4[i][j][k]+G4[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(G5[i-1][j][k]+6*G5[i][j][k]+G5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GG1[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GG2[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GG3[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GG4[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GG5[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

//#pragma omp barrier

/******* Second filter end *******/


MPI_Comm comm;
comm=MPI_COMM_WORLD;
MPI_Status istat[8];

istart=3;
iend = gend[myid];

l_nbr = myid - 1;
r_nbr = myid + 1;
if(myid == 0) l_nbr=MPI_PROC_NULL;
if(myid == nproc-1) r_nbr=MPI_PROC_NULL;

icount = 3*Y_m*Z_m;

itag=110;
MPI_Sendrecv((void *)&GG1[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG1[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=120;
MPI_Sendrecv((void *)&GG2[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG2[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=130;
MPI_Sendrecv((void *)&GG3[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG3[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=140;
MPI_Sendrecv((void *)&GG4[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG4[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=150;
MPI_Sendrecv((void *)&GG5[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG5[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);

itag=160;
MPI_Sendrecv((void *)&GG1[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG1[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=170;
MPI_Sendrecv((void *)&GG2[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG2[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=180;
MPI_Sendrecv((void *)&GG3[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG3[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=190;
MPI_Sendrecv((void *)&GG4[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG4[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=200;
MPI_Sendrecv((void *)&GG5[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG5[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);

MPI_Barrier(MPI_COMM_WORLD);




/******* Third filter *******/


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(GG1[i-1][j][k]+6*GG1[i][j][k]+GG1[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(GG2[i-1][j][k]+6*GG2[i][j][k]+GG2[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(GG3[i-1][j][k]+6*GG3[i][j][k]+GG3[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(GG4[i-1][j][k]+6*GG4[i][j][k]+GG4[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(GG5[i-1][j][k]+6*GG5[i][j][k]+GG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGG1[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGG2[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGG3[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGG4[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGG5[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Third filter end *******/




/******* Forth filter *******/


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(GGG1[i-1][j][k]+6*GGG1[i][j][k]+GGG1[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(GGG2[i-1][j][k]+6*GGG2[i][j][k]+GGG2[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(GGG3[i-1][j][k]+6*GGG3[i][j][k]+GGG3[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(GGG4[i-1][j][k]+6*GGG4[i][j][k]+GGG4[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(GGG5[i-1][j][k]+6*GGG5[i][j][k]+GGG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGGG1[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGGG2[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGGG3[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGGG4[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGGG5[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Forth filter end *******/



/******* Fifth filter *******/


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.125*(GGGG1[i-1][j][k]+6*GGGG1[i][j][k]+GGGG1[i+1][j][k]);
				Gx2[i][j][k] = 0.125*(GGGG2[i-1][j][k]+6*GGGG2[i][j][k]+GGGG2[i+1][j][k]);
				Gx3[i][j][k] = 0.125*(GGGG3[i-1][j][k]+6*GGGG3[i][j][k]+GGGG3[i+1][j][k]);
				Gx4[i][j][k] = 0.125*(GGGG4[i-1][j][k]+6*GGGG4[i][j][k]+GGGG4[i+1][j][k]);
				Gx5[i][j][k] = 0.125*(GGGG5[i-1][j][k]+6*GGGG5[i][j][k]+GGGG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.125*(2*b/a*Gx1[i][j-1][k]*J[i][j-1][k]+6*Gx1[i][j][k]*J[i][j][k]+2*c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.125*(2*b/a*Gx2[i][j-1][k]*J[i][j-1][k]+6*Gx2[i][j][k]*J[i][j][k]+2*c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.125*(2*b/a*Gx3[i][j-1][k]*J[i][j-1][k]+6*Gx3[i][j][k]*J[i][j][k]+2*c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.125*(2*b/a*Gx4[i][j-1][k]*J[i][j-1][k]+6*Gx4[i][j][k]*J[i][j][k]+2*c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.125*(2*b/a*Gx5[i][j-1][k]*J[i][j-1][k]+6*Gx5[i][j][k]*J[i][j][k]+2*c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGGGG1[i][j][k] = 0.125*(Gy1[i][j][k-1]+6*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGGGG2[i][j][k] = 0.125*(Gy2[i][j][k-1]+6*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGGGG3[i][j][k] = 0.125*(Gy3[i][j][k-1]+6*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGGGG4[i][j][k] = 0.125*(Gy4[i][j][k-1]+6*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGGGG5[i][j][k] = 0.125*(Gy5[i][j][k-1]+6*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Fifth filter end *******/

// 6-15 G+20 G^2-15 G^3+6 G^4-G^5 

// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				fU1_h[i][j][k] = 6*G1p[i][j][k]-15*G1[i][j][k]+20*GG1[i][j][k]-15*GGG1[i][j][k]+6*GGGG1[i][j][k]-GGGGG1[i][j][k];
				fU2_h[i][j][k] = 6*G2p[i][j][k]-15*G2[i][j][k]+20*GG2[i][j][k]-15*GGG2[i][j][k]+6*GGGG2[i][j][k]-GGGGG2[i][j][k];
				fU3_h[i][j][k] = 6*G3p[i][j][k]-15*G3[i][j][k]+20*GG3[i][j][k]-15*GGG3[i][j][k]+6*GGGG3[i][j][k]-GGGGG3[i][j][k];
				fU4_h[i][j][k] = 6*G4p[i][j][k]-15*G4[i][j][k]+20*GG4[i][j][k]-15*GGG4[i][j][k]+6*GGGG4[i][j][k]-GGGGG4[i][j][k];
				fU5_h[i][j][k] = 6*G5p[i][j][k]-15*G5[i][j][k]+20*GG5[i][j][k]-15*GGG5[i][j][k]+6*GGGG5[i][j][k]-GGGGG5[i][j][k];

			}
		}
	}
	
#pragma omp barrier


	


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k,Ue,Ve,We)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Ue = U2_[i][j][k]/U1_[i][j][k]-fU2_h[i][j][k]/fU1_h[i][j][k];
				Ve = U3_[i][j][k]/U1_[i][j][k]-fU3_h[i][j][k]/fU1_h[i][j][k];
				We = U4_[i][j][k]/U1_[i][j][k]-fU4_h[i][j][k]/fU1_h[i][j][k];

				E_h[i][j][k] = Ue*Ue+Ve*Ve+We*We;

			}
		}
	}
	
#pragma omp barrier


///////////////////////////////// h ending /////////////////////////////////











///////////////////////////////// 2h starting /////////////////////////////////



/******* The primary filter *******/
	
//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(U1_[i-1][j][k]+2*U1_[i][j][k]+U1_[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(U2_[i-1][j][k]+2*U2_[i][j][k]+U2_[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(U3_[i-1][j][k]+2*U3_[i][j][k]+U3_[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(U4_[i-1][j][k]+2*U4_[i][j][k]+U4_[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(U5_[i-1][j][k]+2*U5_[i][j][k]+U5_[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				
			}
		}
	}
	
#pragma omp barrier


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				G1p[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				G2p[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				G3p[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				G4p[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				G5p[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);

			}
		}
	}

#pragma omp barrier


/******* The primary filter end *******/




/******* First filter *******/

//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(G1p[i-1][j][k]+2*G1p[i][j][k]+G1p[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(G2p[i-1][j][k]+2*G2p[i][j][k]+G2p[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(G3p[i-1][j][k]+2*G3p[i][j][k]+G3p[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(G4p[i-1][j][k]+2*G4p[i][j][k]+G4p[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(G5p[i-1][j][k]+2*G5p[i][j][k]+G5p[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
	for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier

	


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				G1[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				G2[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				G3[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				G4[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				G5[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);

			}
		}
	}

#pragma omp barrier

/******* First filter end *******/




/******* Second filter *******/


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(G1[i-1][j][k]+2*G1[i][j][k]+G1[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(G2[i-1][j][k]+2*G2[i][j][k]+G2[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(G3[i-1][j][k]+2*G3[i][j][k]+G3[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(G4[i-1][j][k]+2*G4[i][j][k]+G4[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(G5[i-1][j][k]+2*G5[i][j][k]+G5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier

	



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GG1[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GG2[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GG3[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GG4[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GG5[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Second filter end *******/





istart=3;
iend = gend[myid];

l_nbr = myid - 1;
r_nbr = myid + 1;
if(myid == 0) l_nbr=MPI_PROC_NULL;
if(myid == nproc-1) r_nbr=MPI_PROC_NULL;

icount = 3*Y_m*Z_m;

itag=110;
MPI_Sendrecv((void *)&GG1[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG1[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=120;
MPI_Sendrecv((void *)&GG2[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG2[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=130;
MPI_Sendrecv((void *)&GG3[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG3[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=140;
MPI_Sendrecv((void *)&GG4[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG4[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
itag=150;
MPI_Sendrecv((void *)&GG5[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&GG5[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);

itag=160;
MPI_Sendrecv((void *)&GG1[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG1[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=170;
MPI_Sendrecv((void *)&GG2[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG2[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=180;
MPI_Sendrecv((void *)&GG3[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG3[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=190;
MPI_Sendrecv((void *)&GG4[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG4[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
itag=200;
MPI_Sendrecv((void *)&GG5[iend-2][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&GG5[istart-3][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);

MPI_Barrier(MPI_COMM_WORLD);






/******* Third filter *******/


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(GG1[i-1][j][k]+2*GG1[i][j][k]+GG1[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(GG2[i-1][j][k]+2*GG2[i][j][k]+GG2[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(GG3[i-1][j][k]+2*GG3[i][j][k]+GG3[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(GG4[i-1][j][k]+2*GG4[i][j][k]+GG4[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(GG5[i-1][j][k]+2*GG5[i][j][k]+GG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier

	


//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 1;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+2;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGG1[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGG2[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGG3[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGG4[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGG5[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Third filter end *******/




/******* Forth filter *******/

//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(GGG1[i-1][j][k]+2*GGG1[i][j][k]+GGG1[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(GGG2[i-1][j][k]+2*GGG2[i][j][k]+GGG2[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(GGG3[i-1][j][k]+2*GGG3[i][j][k]+GGG3[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(GGG4[i-1][j][k]+2*GGG4[i][j][k]+GGG4[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(GGG5[i-1][j][k]+2*GGG5[i][j][k]+GGG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier




//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier

	

//// ============================================ ////
		if (myid ==0) istart = 3;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+1;    			  ////
//// ============================================ ////
	

for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGGG1[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGGG2[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGGG3[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGGG4[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGGG5[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Forth filter end *******/



/******* Fifth filter *******/

// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Gx1[i][j][k] = 0.25*(GGGG1[i-1][j][k]+2*GGGG1[i][j][k]+GGGG1[i+1][j][k]);
				Gx2[i][j][k] = 0.25*(GGGG2[i-1][j][k]+2*GGGG2[i][j][k]+GGGG2[i+1][j][k]);
				Gx3[i][j][k] = 0.25*(GGGG3[i-1][j][k]+2*GGGG3[i][j][k]+GGGG3[i+1][j][k]);
				Gx4[i][j][k] = 0.25*(GGGG4[i-1][j][k]+2*GGGG4[i][j][k]+GGGG4[i+1][j][k]);
				Gx5[i][j][k] = 0.25*(GGGG5[i-1][j][k]+2*GGGG5[i][j][k]+GGGG5[i+1][j][k]);

			}
		}
	}

#pragma omp barrier




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(a,b,c,k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				a = Y_point[i][j+1][k]-Y_point[i][j-1][k]; // H_L+H_R

				b = Y_point[i][j+1][k]-Y_point[i][j][k]; // H_R

				c = Y_point[i][j][k]-Y_point[i][j-1][k]; // H_L

				Gy1[i][j][k] = 0.5*(b/a*Gx1[i][j-1][k]*J[i][j-1][k]+Gx1[i][j][k]*J[i][j][k]+c/a*Gx1[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy2[i][j][k] = 0.5*(b/a*Gx2[i][j-1][k]*J[i][j-1][k]+Gx2[i][j][k]*J[i][j][k]+c/a*Gx2[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy3[i][j][k] = 0.5*(b/a*Gx3[i][j-1][k]*J[i][j-1][k]+Gx3[i][j][k]*J[i][j][k]+c/a*Gx3[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy4[i][j][k] = 0.5*(b/a*Gx4[i][j-1][k]*J[i][j-1][k]+Gx4[i][j][k]*J[i][j][k]+c/a*Gx4[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				Gy5[i][j][k] = 0.5*(b/a*Gx5[i][j-1][k]*J[i][j-1][k]+Gx5[i][j][k]*J[i][j][k]+c/a*Gx5[i][j+1][k]*J[i][j+1][k])/J[i][j][k];
				

			}
		}
	}

#pragma omp barrier

	


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				GGGGG1[i][j][k] = 0.25*(Gy1[i][j][k-1]+2*Gy1[i][j][k]+Gy1[i][j][k+1]);
				GGGGG2[i][j][k] = 0.25*(Gy2[i][j][k-1]+2*Gy2[i][j][k]+Gy2[i][j][k+1]);
				GGGGG3[i][j][k] = 0.25*(Gy3[i][j][k-1]+2*Gy3[i][j][k]+Gy3[i][j][k+1]);
				GGGGG4[i][j][k] = 0.25*(Gy4[i][j][k-1]+2*Gy4[i][j][k]+Gy4[i][j][k+1]);
				GGGGG5[i][j][k] = 0.25*(Gy5[i][j][k-1]+2*Gy5[i][j][k]+Gy5[i][j][k+1]);
				
			}
		}
	}

#pragma omp barrier

/******* Fifth filter end *******/

// 6-15 G+20 G^2-15 G^3+6 G^4-G^5 

// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				fU1_2h[i][j][k] = 6*G1p[i][j][k]-15*G1[i][j][k]+20*GG1[i][j][k]-15*GGG1[i][j][k]+6*GGGG1[i][j][k]-GGGGG1[i][j][k];
				fU2_2h[i][j][k] = 6*G2p[i][j][k]-15*G2[i][j][k]+20*GG2[i][j][k]-15*GGG2[i][j][k]+6*GGGG2[i][j][k]-GGGGG2[i][j][k];
				fU3_2h[i][j][k] = 6*G3p[i][j][k]-15*G3[i][j][k]+20*GG3[i][j][k]-15*GGG3[i][j][k]+6*GGGG3[i][j][k]-GGGGG3[i][j][k];
				fU4_2h[i][j][k] = 6*G4p[i][j][k]-15*G4[i][j][k]+20*GG4[i][j][k]-15*GGG4[i][j][k]+6*GGGG4[i][j][k]-GGGGG4[i][j][k];
				fU5_2h[i][j][k] = 6*G5p[i][j][k]-15*G5[i][j][k]+20*GG5[i][j][k]-15*GGG5[i][j][k]+6*GGGG5[i][j][k]-GGGGG5[i][j][k];

			}
		}
	}
	
#pragma omp barrier




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k, Ue,Ve,We)
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				Ue = U2_[i][j][k]/U1_[i][j][k]-fU2_2h[i][j][k]/fU1_2h[i][j][k];
				Ve = U3_[i][j][k]/U1_[i][j][k]-fU3_2h[i][j][k]/fU1_2h[i][j][k];
				We = U4_[i][j][k]/U1_[i][j][k]-fU4_2h[i][j][k]/fU1_2h[i][j][k];

				E_2h[i][j][k] = Ue*Ue+Ve*Ve+We*We;

			}
		}
	}
	
#pragma omp barrier





if (myid == 0) {

	for (j = 2; j <= ny; j++) {
		for (k = 2; k <= nz; k++) {

			istart = 2;

			E_h[istart][j][k] = E_h[istart+1][j][k];
			E_2h[istart][j][k] = E_2h[istart+1][j][k];

		}
	}

}



if (myid == nproc-1) {

	for (j = 2; j <= ny; j++) {
		for (k = 2; k <= nz; k++) {

			iend = gend[myid];

			E_h[iend+1][j][k] = E_h[iend][j][k];
			E_2h[iend+1][j][k] = E_2h[iend][j][k];

		}
	}
}



// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

for (i = istart; i <= iend; i++) {
	for (k = 2; k <= nz; k++) {

		E_h[i][1][k] = E_h[i][2][k];
		E_h[i][nyy][k] = E_h[i][ny][k];

		E_2h[i][1][k] = E_2h[i][2][k];
		E_2h[i][nyy][k] = E_2h[i][ny][k];

	}
}




// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

for (i = istart; i <= iend; i++) {
	for (j = 2; j <= ny; j++) {

		E_h[i][j][1] = E_h[i][j][2];
		E_h[i][j][nzz] = E_h[i][j][nz];

		E_2h[i][j][1] = E_2h[i][j][2];
		E_2h[i][j][nzz] = E_2h[i][j][nz];

	}
}




	


// ===================== //
	istart = 2;          //
	iend = gend[myid];   //
// ===================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k,Emax)

		for (j = 1; j < ny; j++) {
			for (k = 1; k < nz; k++) {

				Emax = (E_h[i][j+1][k+1]+E_h[i+1][j+1][k+1])/(E_2h[i][j+1][k+1]+E_2h[i+1][j+1][k+1]);


				if ( Emax > E_high ) EpX[i][j][k] = min(EpX[i][j][k]+Roe_criterion, 1.0);

				if ( Emax < E_low  ) EpX[i][j][k] = max(EpX[i][j][k]-Roe_criterion, 0.0); 
				

			}
		}
	}


	


// ===================== //
	istart = 2;          //
	iend = gend[myid]-1; //
// ===================== //


for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k,Emax)

		for (j = 1; j < nyy; j++) {
			for (k = 1; k < nz; k++) {

				Emax = (E_h[i+1][j][k+1]+E_h[i+1][j+1][k+1])/(E_2h[i+1][j][k+1]+E_2h[i+1][j+1][k+1]);


				
				if ( Emax > E_high ) EpY[i][j][k] = min(EpY[i][j][k]+Roe_criterion, 1.0);

				if ( Emax < E_low  ) EpY[i][j][k] = max(EpY[i][j][k]-Roe_criterion, 0.0); 
				

			}
		}
	}



	

// ===================== //
	istart = 2;          //
	iend = gend[myid]-1; //
// ===================== //
	
for (i = istart; i <= iend; i++) {

#pragma omp parallel for private(k,Emax)
	
		for (j = 1; j < ny; j++) {
			for (k = 1; k < nzz; k++) {

				Emax = (E_h[i+1][j+1][k]+E_h[i+1][j+1][k+1])/(E_2h[i+1][j+1][k]+E_2h[i+1][j+1][k+1]);

				
				if ( Emax > E_high ) EpZ[i][j][k] = min(EpZ[i][j][k]+Roe_criterion, 1.0);

				if ( Emax < E_low  ) EpZ[i][j][k] = max(EpZ[i][j][k]-Roe_criterion, 0.0); 
				
					

			}
		}
	}











	delete [] G1p;
	delete [] G2p;
	delete [] G3p;
	delete [] G4p;
	delete [] G5p;

	delete [] G1;
	delete [] G2;
	delete [] G3;
	delete [] G4;
	delete [] G5;

	delete [] GG1;
	delete [] GG2;
	delete [] GG3;
	delete [] GG4;
	delete [] GG5;

	delete [] GGG1;
	delete [] GGG2;
	delete [] GGG3;
	delete [] GGG4;
	delete [] GGG5;

	delete [] GGGG1;
	delete [] GGGG2;
	delete [] GGGG3;
	delete [] GGGG4;
	delete [] GGGG5;

	delete [] GGGGG1;
	delete [] GGGGG2;
	delete [] GGGGG3;
	delete [] GGGGG4;
	delete [] GGGGG5;

	delete [] Gx1;
	delete [] Gx2;
	delete [] Gx3;
	delete [] Gx4;
	delete [] Gx5;

	delete [] Gy1;
	delete [] Gy2;
	delete [] Gy3;
	delete [] Gy4;
	delete [] Gy5;

	delete [] fU1_h;
	delete [] fU2_h;
	delete [] fU3_h;
	delete [] fU4_h;
	delete [] fU5_h;

	delete [] fU1_2h;
	delete [] fU2_2h;
	delete [] fU3_2h;
	delete [] fU4_2h;
	delete [] fU5_2h;

	delete [] E_2h;
	delete [] E_h;


}