



#include <omp.h>
#include "Resolution.h"

void X_boundary_condition
// ============================================================================ //
(
double (*U1_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U2_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U3_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U4_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 
double (*U5_)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory], 

double (*J)[Y_memory][Z_memory] = new double[X_memory][Y_memory][Z_memory]
)
// ============================================================================ //
{

#include "ijk.h"

#pragma omp parallel for private(k)
	for (j = 2; j < nyy; j++) {
		for (k = 2; k < nzz; k++) {

			U1_[1][j][k] = U1_[nx][j][k];
			U2_[1][j][k] = U2_[nx][j][k];
			U3_[1][j][k] = U3_[nx][j][k];
			U4_[1][j][k] = U4_[nx][j][k];
			U5_[1][j][k] = U5_[nx][j][k];

			U1_[nxx][j][k] = U1_[2][j][k];
			U2_[nxx][j][k] = U2_[2][j][k];
			U3_[nxx][j][k] = U3_[2][j][k];
			U4_[nxx][j][k] = U4_[2][j][k];
			U5_[nxx][j][k] = U5_[2][j][k];

		}
	}
#pragma omp barrier

}
