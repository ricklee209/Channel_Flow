



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int X_np;

void Hyperplane 
(
// =========================================================== //
int myid,

int (*Nijk)[2] = new int[(X_np-3)+ny+nz+7][2],

int (*ijk)[3] = new int[(X_np-6)*Y_out*Z_out*3][3]

// =========================================================== //
)

{

#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"


	int m, mi;

	int ii = 0;

	int Ntemp;

//// ======================= ////
		istart = 3;		     ////	
//// ======================= ////
		iend = gend[myid];   ////
//// ======================= ////

	for (m = 7 ; m <= iend+ny+nz; m++) { 

		mi = 0;
		for (i = istart ; i <= iend; i++) { 
			for (j = 2; j < nyy; j++) {
				for (k = 2; k < nzz; k++) {

					if (m == i+j+k) {

						mi = mi+1;

						ii = ii+1;

						ijk[ii][0] = i;
						ijk[ii][1] = j;
						ijk[ii][2] = k;

					}    // ---- if (m == i+j+k) ---- //

				}
			}
		}

		Nijk[m][0] = ii;

		Nijk[m][1] = mi;

		
	}    // ---- for (m = 7 ; m <= (Xnp-3)+ny+nz; i++) ---- //


}