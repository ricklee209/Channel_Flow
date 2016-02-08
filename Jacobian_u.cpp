



#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

#define PI 3.141592653589793

void Jacobian_u	
// ============================================================================ //
(
int myid,

double (*J_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
// =========================================================================== //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double temp,temp1,temp2,temp3,temp4;

	double xdxi_u,xdet_u,xdzt_u,ydxi_u,ydet_u,ydzt_u,zdxi_u,zdet_u,zdzt_u;
	
	double gamma = 1.5;

//// ============================================ ////
		istart = 2;							      ////
//// ============================================ ////
		iend = gend[myid];			     		  ////
//// ============================================ ////

//#pragma omp parallel for private(j,k,temp,temp1,temp2,temp4)
	for (i = istart; i <=  iend; i++) {
		for (j = 1; j < ny; j++) {
			for (k = 1; k < nz; k++) {
				
				xidx_u[i][j][k]=1;
				xidy_u[i][j][k]=0;
				xidz_u[i][j][k]=0;
				xdxi_u=1;
				xdet_u=0;
				xdzt_u=0;

				etdx_u[i][j][k]=0;

				
				temp = (0.5*high)*(1-1./tanh(gamma)*tanh(gamma*(1-2*(j-0.5)*deltaET)));
				temp1 = (1-temp/(0.5*high))*tanh(gamma);
				

				etdy_u[i][j][k]= tanh(gamma)/(2*(0.5*high)*gamma*(1-temp1*temp1));

				etdz_u[i][j][k]=0;
				ydxi_u=0;

				temp2 = 1/cosh(gamma*(1-2*(j-0.5)*deltaET));

				ydet_u = high*gamma*1./tanh(gamma)*temp2*temp2;

				ydzt_u=0;

				ztdx_u[i][j][k]=0;
				ztdy_u[i][j][k]=0;
				ztdz_u[i][j][k]=1;
				zdxi_u=0;
				zdet_u=0;
				zdzt_u=1;

				J_u[i][j][k]=1./(xdxi_u*ydet_u*zdzt_u);
				
			}
		}
	}
#pragma omp barrier



}

