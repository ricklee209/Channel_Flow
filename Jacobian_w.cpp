


#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

#define PI 3.141592653589793

void Jacobian_w	
(
// ============================================================================ //
int myid,

double (*J_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// =========================================================================== //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"


	double temp,temp1,temp2,temp3,temp4;

	double xdxi_w,xdet_w,xdzt_w,ydxi_w,ydet_w,ydzt_w,zdxi_w,zdet_w,zdzt_w;
	double gamma = 1.5;

//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

#pragma omp parallel for private(j,k,temp,temp1,temp2)

	for (i = istart; i <=  iend; i++) {
		for (j = 1; j < ny; j++) {
			for (k = 1; k <= nz; k++) {   

				xidx_w[i][j][k]=1;
				xidy_w[i][j][k]=0;
				xidz_w[i][j][k]=0;
				xdxi_w=1;
				xdet_w=0;
				xdzt_w=0;

				etdx_w[i][j][k]=0;

				
				temp = (0.5*high)*(1-1./tanh(gamma)*tanh(gamma*(1-2*(j-0.5)*deltaET)));
				temp1 = (1-temp/(0.5*high))*tanh(gamma);
				

				etdy_w[i][j][k]= tanh(gamma)/(2*(0.5*high)*gamma*(1-temp1*temp1));

				etdz_w[i][j][k]=0;
				ydxi_w=0;

				temp2 = 1/cosh(gamma*(1-2*(j-0.5)*deltaET));

				ydet_w = high*gamma*1./tanh(gamma)*temp2*temp2;

				ydzt_w = 0;


				ztdx_w[i][j][k]=0;
				ztdy_w[i][j][k]=0;
				ztdz_w[i][j][k]=1;
				zdxi_w=0;
				zdet_w=0;
				zdzt_w=1;

				J_w[i][j][k]=1./(xdxi_w*ydet_w*zdzt_w);
				
			}
		}
	}


}

