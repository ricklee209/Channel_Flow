



#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

#define PI 3.141592653589793

void Jacobian_v	
(
// ============================================================================ //
int myid,

double (*J_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// =========================================================================== //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double temp,temp1,temp2,temp3,temp4;

	double xdxi_v,xdet_v,xdzt_v,ydxi_v,ydet_v,ydzt_v,zdxi_v,zdet_v,zdzt_v;
	
	
	double gamma = 1.5;


	
//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////


#pragma omp parallel for private(j,k,temp,temp1,temp2)
	for (i = istart; i <=  iend; i++) {
		for (j = 1; j <= ny; j++) {
			for (k = 1; k < nz; k++) {   

				xidx_v[i][j][k]=1;
				xidy_v[i][j][k]=0;
				xidz_v[i][j][k]=0;
				xdxi_v=1;
				xdet_v=0;
				xdzt_v=0;

				etdx_v[i][j][k]=0;

				temp = (0.5*high)*(1-1./tanh(gamma)*tanh(gamma*(1-2*(j-1.0)*deltaET)));
				temp1 = (1-temp/(0.5*high))*tanh(gamma);
				
				etdy_v[i][j][k]= tanh(gamma)/(2*(0.5*high)*gamma*(1-temp1*temp1));

				etdz_v[i][j][k]=0;
				ydxi_v=0;

				temp2 = 1/cosh(gamma*(1-2*(j-1.0)*deltaET));

				ydet_v = high*gamma*1./tanh(gamma)*temp2*temp2;

				ydzt_v=0;

				ztdx_v[i][j][k]=0;
				ztdy_v[i][j][k]=0;
				ztdz_v[i][j][k]=1;
				zdxi_v=0;
				zdet_v=0;
				zdzt_v=1;

				J_v[i][j][k]=1./(xdxi_v*ydet_v*zdzt_v);
				
			}
		}
	}
	

//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

	for (i = istart; i <=  iend; i++) {
		for (k = 1; k < nz; k++) {   


			etdx_v[i][1][k]=etdx[i][1][k];
			etdy_v[i][1][k]=etdy[i][1][k];
			etdz_v[i][1][k]=etdz[i][1][k];
			
			etdx_v[i][ny][k]=etdx[i][nyy][k];
			etdy_v[i][ny][k]=etdy[i][nyy][k];
			etdz_v[i][ny][k]=etdz[i][nyy][k];
			
			J_v[i][ny][k]=J[i][nyy][k];

			J_v[i][1][k]=J[i][1][k];

		}
	}


}

