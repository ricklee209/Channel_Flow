




#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <omp.h>

#include "Resolution.h"

extern int X_np;

void Z_boundary_condition
// ============================================================================ //
(
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
				
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 0;						  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+3;				  ////
//// ============================================ ////
		
//#pragma omp parallel for private(j)
	for (i = istart; i <= iend; i++) {
		for (j = 0; j <= nyyy; j++) {
			/* back side */
			U1_[i][j][1] = U1_[i][j][nz];
			U2_[i][j][1] = U2_[i][j][nz];
			U3_[i][j][1] = U3_[i][j][nz];
			U4_[i][j][1] = U4_[i][j][nz];
			U5_[i][j][1] = U5_[i][j][nz];
			

			/* front side */
			U1_[i][j][nzz] = U1_[i][j][2];
			U2_[i][j][nzz] = U2_[i][j][2];
			U3_[i][j][nzz] = U3_[i][j][2];
			U4_[i][j][nzz] = U4_[i][j][2];
			U5_[i][j][nzz] = U5_[i][j][2];
			
		}
	}
#pragma omp barrier

}
