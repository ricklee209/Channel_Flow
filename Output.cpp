



#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#include "Resolution.h"

extern int X_np;

void Output
(
// ============================================================================ //
int step,

int myid,

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

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "MPI_prm.h"
#include "Mpi_division.h"

	double (*U1out)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U2out)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U3out)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U4out)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U5out)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];

	double (*U1out_old)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U2out_old)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U3out_old)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U4out_old)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	double (*U5out_old)[Y_out][Z_out] = new double[X_out+1][Y_out][Z_out];
	
	        
	int ii = 2;

//// ========================= ////
	  istart = gstart[myid];   ////	
//// ========================= ////
	  iend = gend0[myid];      ////
//// ========================= ////

		for (i = istart; i <= iend; i++) {

			ii = ii+1;

#pragma omp parallel for private(k)
			for (j = 0; j < ny-1; j++) { 
				for (k = 0; k < nz-1; k++) {


					U1out[i][j][k] = U1[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U2out[i][j][k] = U2[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U3out[i][j][k] = U3[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U4out[i][j][k] = U4[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U5out[i][j][k] = U5[ii][j+2][k+2]*J[ii][j+2][k+2];

					U1out_old[i][j][k] = U1q[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U2out_old[i][j][k] = U2q[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U3out_old[i][j][k] = U3q[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U4out_old[i][j][k] = U4q[ii][j+2][k+2]*J[ii][j+2][k+2]; 
					U5out_old[i][j][k] = U5q[ii][j+2][k+2]*J[ii][j+2][k+2];

					
				}
			}
		}


		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];

		if (myid > 0) {

			istart=gstart[myid];
			icount=gcount[myid]*Y_out*Z_out;
			idest=0;

			itag = 210;
			MPI_Send((void *)&U1out[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 220;
			MPI_Send((void *)&U2out[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 230;
			MPI_Send((void *)&U3out[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 240;
			MPI_Send((void *)&U4out[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 250;
			MPI_Send((void *)&U5out[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);

			itag = 260;
			MPI_Send((void *)&U1out_old[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 270;
			MPI_Send((void *)&U2out_old[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 280;
			MPI_Send((void *)&U3out_old[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 290;
			MPI_Send((void *)&U4out_old[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);
			itag = 300;
			MPI_Send((void *)&U5out_old[istart][0][0], icount, MPI_DOUBLE, idest, itag, comm);

		}

		else {

			for ( isrc=1; isrc < nproc; isrc++ ) {

				istart=gstart[isrc];
				icount=gcount[isrc]*Y_out*Z_out;

				itag = 210;
				MPI_Recv((void *)&U1out[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 220;
				MPI_Recv((void *)&U2out[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 230;
				MPI_Recv((void *)&U3out[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 240;
				MPI_Recv((void *)&U4out[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 250;
				MPI_Recv((void *)&U5out[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);

				itag = 260;
				MPI_Recv((void *)&U1out_old[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 270;
				MPI_Recv((void *)&U2out_old[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 280;
				MPI_Recv((void *)&U3out_old[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 290;
				MPI_Recv((void *)&U4out_old[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);
				itag = 300;
				MPI_Recv((void *)&U5out_old[istart][0][0], icount, MPI_DOUBLE, isrc, itag, comm, istat);

			}

		}



		if (myid == 0) {

			char LESdata[100];
			FILE *fptr;
			sprintf(LESdata,"LES""%0.5d"".bin",step);
			fptr = fopen(LESdata,"wb");

			if (fptr != NULL)
			{

				fwrite(U1out,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U2out,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U3out,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U4out,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U5out,sizeof(double),X_out*Y_out*Z_out,fptr);

				
				fwrite(U1out_old,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U2out_old,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U3out_old,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U4out_old,sizeof(double),X_out*Y_out*Z_out,fptr);
				fwrite(U5out_old,sizeof(double),X_out*Y_out*Z_out,fptr);
				

				fclose(fptr);

			}

			else printf("File opening Failure\n");

		}


	delete [] U1out;
	delete [] U2out;
	delete [] U3out;
	delete [] U4out;
	delete [] U5out;

	delete [] U1out_old;
	delete [] U2out_old;
	delete [] U3out_old;
	delete [] U4out_old;
	delete [] U5out_old;


}