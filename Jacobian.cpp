



#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

#define PI 3.141592653589793

void Jacobian	
(
// ===================================================== //
int myid,

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
// ===================================================== //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double temp,temp1,temp2,temp3,temp4;

	double xdxi,xdet,xdzt,ydxi,ydet,ydzt,zdxi,zdet,zdzt;
	
	
	double gamma = 1.5;
	
//// ============================================ ////
		if (myid ==0) istart = 3;		          ////
		else istart = 0;			              ////
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid];    ////
		else iend = gend[myid]+3;	     		  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {   

				
				xidx[i][j][k]=1;
				xidy[i][j][k]=0;
				xidz[i][j][k]=0;
				
				
				xdxi=1;
				xdet=0;
				xdzt=0;

				etdx[i][j][k]=0;

				temp = (0.5*high)*(1-1./tanh(gamma)*tanh(gamma*(1-2*(j-1.5)*deltaET)));
				temp1 = (1-temp/(0.5*high))*tanh(gamma);
				
				etdy[i][j][k] = tanh(gamma)/(2*(0.5*high)*gamma*(1-temp1*temp1));

				etdz[i][j][k]=0;
				ydxi=0;

				temp2 = 1/cosh(gamma*(1-2*(j-1.5)*deltaET));

				ydet = high*gamma*1./tanh(gamma)*temp2*temp2;

				ydzt=0;

				ztdx[i][j][k]=0;
				ztdy[i][j][k]=0;
				ztdz[i][j][k]=1;
				
				zdxi=0;
				zdet=0;
				zdzt=1;

				/*J(j,i,k)=1/(xdxi(j,i,k)*ydET(j,i,k)*zdzt(j,i,k)+xdzt(j,i,k)*ydxi(j,i,k)*zdET(j,i,k)+xdET(j,i,k)*ydzt(j,i,k)*zdxi(j,i,k)-xdxi(j,i,k)*ydzt(j,i,k)*zdET(j,i,k)-xdET(j,i,k)*ydxi(j,i,k)*zdzt(j,i,k)-xdzt(j,i,k)*ydET(j,i,k)*zdxi(j,i,k));*/
				J[i][j][k]=1./(xdxi*ydet*zdzt);
				
				
				
				//if (myid == 1 && i == 2 && k == 10)printf("%f\t%f\n",J[i][j][k],etdy[i][j][k]);
				

			}
		}
	}

	
	if (myid == 0) {

		for (j = 1; j <= nyy; j++) {
			for (k = 1; k <= nzz; k++) {   

				xidx[2][j][k]=xidx[3][j][k];
				xidy[2][j][k]=xidy[3][j][k];
				xidz[2][j][k]=xidz[3][j][k];

				etdx[2][j][k]=etdx[3][j][k];
				etdy[2][j][k]=etdy[3][j][k];
				etdz[2][j][k]=etdz[3][j][k];

				ztdx[2][j][k] = ztdx[3][j][k];
				ztdx[2][j][k] = ztdx[3][j][k];
				ztdx[2][j][k] = ztdx[3][j][k];

				J[2][j][k]=J[3][j][k];

			}
		}

	}

	if (myid == nproc-1) {

		iend = gend[myid]+1;

		for (j = 1; j <= nyy; j++) {
			for (k = 1; k <= nzz; k++) {   

				xidx[iend][j][k]=xidx[iend-1][j][k];
				xidy[iend][j][k]=xidy[iend-1][j][k];
				xidz[iend][j][k]=xidz[iend-1][j][k];

				etdx[iend][j][k]=etdx[iend-1][j][k];
				etdy[iend][j][k]=etdy[iend-1][j][k];
				etdz[iend][j][k]=etdz[iend-1][j][k];

				ztdx[iend][j][k] = ztdx[iend-1][j][k];
				ztdx[iend][j][k] = ztdx[iend-1][j][k];
				ztdx[iend][j][k] = ztdx[iend-1][j][k];

				J[iend][j][k]=J[iend-1][j][k];

			}
		}

	}


//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;			              ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;	     		  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {   

			xidx[i][1][k]=xidx[i][2][k];
			xidy[i][1][k]=xidy[i][2][k];
			xidz[i][1][k]=xidz[i][2][k];

			etdx[i][1][k]=etdx[i][2][k];
			etdy[i][1][k]=etdy[i][2][k];
			etdz[i][1][k]=etdz[i][2][k];

			ztdx[i][1][k]=ztdx[i][2][k];
			ztdy[i][1][k]=ztdx[i][2][k];
			ztdz[i][1][k]=ztdx[i][2][k];
			
			xidx[i][nyy][k]=xidx[i][ny][k];
			xidy[i][nyy][k]=xidy[i][ny][k];
			xidz[i][nyy][k]=xidz[i][ny][k];

			etdx[i][nyy][k]=etdx[i][ny][k];
			etdy[i][nyy][k]=etdy[i][ny][k];
			etdz[i][nyy][k]=etdz[i][ny][k];

			ztdx[i][nyy][k]=ztdx[i][ny][k];
			ztdy[i][nyy][k]=ztdx[i][ny][k];
			ztdz[i][nyy][k]=ztdx[i][ny][k];
			
			J[i][1][k]=J[i][2][k];
			J[i][nyy][k]=J[i][ny][k];
			
		}
	}



//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;			              ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;	     		  ////
//// ============================================ ////

	for (i = istart; i <= iend; i++) {
		for (j = 1; j <= nyy; j++) {   

			xidx[i][j][1]=xidx[i][j][2];
			xidy[i][j][1]=xidy[i][j][2];
			xidz[i][j][1]=xidz[i][j][2];

			etdx[i][j][1]=etdx[i][j][2];
			etdy[i][j][1]=etdy[i][j][2];
			etdz[i][j][1]=etdz[i][j][2];

			ztdx[i][j][1]=ztdx[i][j][2];
			ztdy[i][j][1]=ztdx[i][j][2];
			ztdz[i][j][1]=ztdx[i][j][2];
			
			xidx[i][j][nzz]=xidx[i][j][nz];
			xidy[i][j][nzz]=xidy[i][j][nz];
			xidz[i][j][nzz]=xidz[i][j][nz];

			etdx[i][j][nzz]=etdx[i][j][nz];
			etdy[i][j][nzz]=etdy[i][j][nz];
			etdz[i][j][nzz]=etdz[i][j][nz];

			ztdx[i][j][nzz]=ztdx[i][j][nz];
			ztdy[i][j][nzz]=ztdx[i][j][nz];
			ztdz[i][j][nzz]=ztdx[i][j][nz];
			
			J[i][j][1]=J[i][j][2];
			J[i][j][nzz]=J[i][j][nz];
			
		}
	}
	

}

