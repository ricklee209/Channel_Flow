



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>

#include "Resolution.h"

extern int X_np;

void Driven_force
// ============================================================================ //
(
int myid,

double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
)
// ============================================================================ //
{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double ER_Q;
	double Q = 0;
	double Qall = 0;
	double temp = 0;
	double Jsum = 0;

	for (j = 2; j <= ny; j++) {

		Jsum = Jsum + 1./J[3][j][3];

	}


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //

	for (i = istart; i <= iend; i++) {
		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {     
				
				temp = (U2_[i][j][k]*J[i][j][k])*1./J[i][j][k]/Jsum;

				Q = Q+temp;
				
			}
		}
	}

	Q = Q/X_out/Z_out;

	MPI_Comm comm;
	comm=MPI_COMM_WORLD;

	idest = 0;

	MPI_Reduce ((void*)&Q, (void*)&Qall, 1, MPI_DOUBLE ,MPI_SUM, idest, comm);

	ER_Q = fabs(Q0-Qall)/Q0;

	if (myid == 0) printf("Q0=%f\tQall=%f\tER_Q=%f\n",Q0,Qall,ER_Q);

}