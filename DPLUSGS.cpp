



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 

#include "Resolution.h"

extern int X_np;

void DPLUSGS
(
// ============================================================================ //
int myid,

double deltaT,

double *er1,
double *er2,
double *er3,
double *er4,
double *er5,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFx1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*inFz1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*vF2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*vF5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 


/**** MR = RHS ****/
double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** MR = RHS-end ****/


/**** NR = dU*p2 ****/
double (*NR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** NR = dU*p2-end ****/


/**** NL = dU*p1 ****/
double (*NL1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
/**** NL = dU*p1-end ****/

double (*Residual1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*Residual2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*Residual3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*Residual4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*Residual5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]

// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"
#include "MPI_prm.h"
#include "Mpi_division.h"



	double rho,U,V,W,VV,P,C,T,H;
	double u,v,w,beta;
	double temp,temp1,temp2;

	double lambda,lambdaX,lambdaY,lambdaZ;

	double LL1,LL2,LL3,LL4,LL5,
	       UU1,UU2,UU3,UU4,UU5,

		   Rp1,Rp2,Rp3,Rp4,Rp5,
		   Rf1,Rf2,Rf3,Rf4,Rf5,
		   Rk1,Rk2,Rk3,Rk4,Rk5,
		   Sx,Sy,Sz;

	double d11,d12,d13,d14,d15,
	       d21,d22,d23,d24,d25,
	       d31,d32,d33,d34,d35,
	       d41,d42,d43,d44,d45,
	       d51,d52,d53,d54,d55;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

	double Aplus11,Aplus12,Aplus13,Aplus14,
           Aplus21,Aplus22,
	       Aplus31,        Aplus33,
		   Aplus41,          	   Aplus44,
	       Aplus51,Aplus52,Aplus53,Aplus54,Aplus55;

	double Bplus11,Bplus12,Bplus13,Bplus14,
           Bplus21,Bplus22,
	       Bplus31,        Bplus33,
		   Bplus41,          	   Bplus44,
	       Bplus51,Bplus52,Bplus53,Bplus54,Bplus55;

	double Cplus11,Cplus12,Cplus13,Cplus14,
           Cplus21,Cplus22,
	       Cplus31,        Cplus33,
		   Cplus41,          	   Cplus44,
	       Cplus51,Cplus52,Cplus53,Cplus54,Cplus55;

	double Amius11,Amius22,Amius33,Amius44,Amius55,
		   Bmius11,Bmius22,Bmius33,Bmius44,Bmius55,
		   Cmius11,Cmius22,Cmius33,Cmius44,Cmius55;

	double rhoold,Uold,Vold,Wold,VVold,Pold,Told;

	double e1 = 0;
	double e2 = 0;
	double e3 = 0;
	double e4 = 0;
	double e5 = 0;

	double et1,et2,et3,et4,et5;
	
	
	double DN = 1./(X_out*Y_out*Z_out);


	if (myid == 0) {
		//#pragma omp parallel for private(k)
		for (j = 1; j <= nyy; j++) {
			for (k = 1; k <= nzz; k++) {

				NL1[2][j][k] = 0;
				NL2[2][j][k] = 0;
				NL3[2][j][k] = 0;
				NL4[2][j][k] = 0;
				NL5[2][j][k] = 0;

				NR1[2][j][k] = 0;
				NR2[2][j][k] = 0;
				NR3[2][j][k] = 0;
				NR4[2][j][k] = 0;
				NR5[2][j][k] = 0;

			}
		}
	}

	if (myid == nproc-1) {

		for (j = 1; j <= nyy; j++) {
			for (k = 1; k <= nzz; k++) {

				iend = gend[myid];

				NL1[iend+1][j][k] = 0;
				NL2[iend+1][j][k] = 0;
				NL3[iend+1][j][k] = 0;
				NL4[iend+1][j][k] = 0;
				NL5[iend+1][j][k] = 0;

				NR1[iend+1][j][k] = 0;
				NR2[iend+1][j][k] = 0;
				NR3[iend+1][j][k] = 0;
				NR4[iend+1][j][k] = 0;
				NR5[iend+1][j][k] = 0;

			}
		}
	}

//// ============================================ ////
		istart = 2;		                          ////	
//// ============================================ ////
		iend = gend[myid]+1;                      ////
//// ============================================ ////

//#pragma omp parallel for private(k)
	for (i = istart; i <= iend; i++) {
		for (k = 1; k <= nzz; k++) {


			NL1[i][1][k] = 0;
			NL2[i][1][k] = 0;
			NL3[i][1][k] = 0;
			NL4[i][1][k] = 0;
			NL5[i][1][k] = 0;

			NR1[i][1][k] = 0;
			NR2[i][1][k] = 0;
			NR3[i][1][k] = 0;
			NR4[i][1][k] = 0;
			NR5[i][1][k] = 0;


			NL1[i][nyy][k] = 0;
			NL2[i][nyy][k] = 0;
			NL3[i][nyy][k] = 0;
			NL4[i][nyy][k] = 0;
			NL5[i][nyy][k] = 0;

			NR1[i][nyy][k] = 0;
			NR2[i][nyy][k] = 0;
			NR3[i][nyy][k] = 0;
			NR4[i][nyy][k] = 0;
			NR5[i][nyy][k] = 0;

		}
	}
#pragma omp barrier

	
//// ============================================ ////
		istart = 2;		                          ////	
//// ============================================ ////
		iend = gend[myid]+1;                      ////
//// ============================================ ////

//#pragma omp parallel for private(j)
	for (i = istart; i <= iend; i++) {
		for (j = 1; j <= nyy; j++) {

			NL1[i][j][1] = 0;
			NL2[i][j][1] = 0;
			NL3[i][j][1] = 0;
			NL4[i][j][1] = 0;
			NL5[i][j][1] = 0;

			NR1[i][j][1] = 0;
			NR2[i][j][1] = 0;
			NR3[i][j][1] = 0;
			NR4[i][j][1] = 0;
			NR5[i][j][1] = 0;

			NL1[i][j][nzz] = 0;
			NL2[i][j][nzz] = 0;
			NL3[i][j][nzz] = 0;
			NL4[i][j][nzz] = 0;
			NL5[i][j][nzz] = 0;

			NR1[i][j][nzz] = 0;
			NR2[i][j][nzz] = 0;
			NR3[i][j][nzz] = 0;
			NR4[i][j][nzz] = 0;
			NR5[i][j][nzz] = 0;

		}
	}
	//#pragma omp barrier


//// ============================================ ////
			istart =  3;	            		  ////	
//// ============================================ ////
			iend = gend[myid];		    		  ////
//// ============================================ ////

#pragma omp parallel for private(\
j,k,_j,_k,\
rho,u,v,w,U,V,W,VV,P,C,T,H,temp,temp1,temp2,\
d11,d12,d13,d14,d15,\
d21,d22,d23,d24,d25,\
d31,d32,d33,d34,d35,\
d41,d42,d43,d44,d45,\
d51,d52,d53,d54,d55,\
Rp1,Rp2,Rp3,Rp4,Rp5,Rf1,Rf2,Rf3,Rf4,Rf5,Rk1,Rk2,Rk3,Rk4,Rk5,\
xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
Sx,Sy,Sz,beta,\
lambda,lambdaX,lambdaY,lambdaZ\
)

	for (i = istart ; i <= iend; i++) {
		for (j = 2, _j = 1; j < nyy; j++, _j++) {
			for (k = 2, _k = 1; k < nzz; k++, _k++) {

				/* flux parameter */
				rho = U1_[i][j][k]*J[i][j][k];
				u = U2_[i][j][k]/U1_[i][j][k];
				v = U3_[i][j][k]/U1_[i][j][k];
				w = U4_[i][j][k]/U1_[i][j][k];     
				VV = u*u+v*v+w*w;
				P = (U5_[i][j][k]*J[i][j][k]-0.5*rho*VV)*(K-1);
				C = K*P/rho; /**** C = sqrt(K*P/rho); ****/
				T = P/rho;
				H = 0.5*VV+C/(K-1);
				

				/* preconditioning */
				/*
				if (VV/C <= e)
					beta = e;
				else
					beta = VV/C;
				*/

				beta = max(VV/C,e);

				/* Inv(Gamma) */
				temp = rho*K;

				d11 = -(H*(K-1)+VV-K*(T+VV))*beta;
				d12 = -(K-1)*u*beta;
				d13 = -(K-1)*v*beta;
				d14 = -(K-1)*w*beta;
				d15 = (K-1)*beta;

				d21 = -u/rho;                         
				d22 = 1/rho;    
				d23 = 0;
				d24 = 0;
				d25 = 0;

				d31 = -v/rho;                      
				d32 = 0;
				d33 = 1/rho;   
				d34 = 0;
				d35 = 0;

				d41 = -w/rho;                       
				d42 = 0;
				d43 = 0;
				d44 = 1/rho;   
				d45 = 0;

				d51 = (K-1)*(VV+K*T*beta+(K-1)*VV*beta+H*(-1+beta-K*beta))/temp;
				d52 = -(K-1)*u*(1+(K-1)*beta)/temp;
				d53 = -(K-1)*v*(1+(K-1)*beta)/temp;
				d54 = -(K-1)*w*(1+(K-1)*beta)/temp;
				d55 = (K-1)*(1+(K-1)*beta)/temp;

				Rp1 = -(3*U1_[i][j][k]-4*U1[i][j][k]+U1q[i][j][k])/(2*deltaT);
				Rp2 = -(3*U2_[i][j][k]-4*U2[i][j][k]+U2q[i][j][k])/(2*deltaT);
				Rp3 = -(3*U3_[i][j][k]-4*U3[i][j][k]+U3q[i][j][k])/(2*deltaT);
				Rp4 = -(3*U4_[i][j][k]-4*U4[i][j][k]+U4q[i][j][k])/(2*deltaT);
				Rp5 = -(3*U5_[i][j][k]-4*U5[i][j][k]+U5q[i][j][k])/(2*deltaT);
				
				Rf1 = -((inFx1[i][_j][_k]-inFx1[i-1][_j][_k])/deltaXI+
					(inFy1[i-1][j][_k]-inFy1[i-1][_j][_k])/deltaET+
					(inFz1[i-1][_j][k]-inFz1[i-1][_j][_k])/deltaZT);
				Rf2 = -((inFx2[i][_j][_k]-inFx2[i-1][_j][_k])/deltaXI+
					(inFy2[i-1][j][_k]-inFy2[i-1][_j][_k])/deltaET+
					(inFz2[i-1][_j][k]-inFz2[i-1][_j][_k])/deltaZT);
				Rf3 = -((inFx3[i][_j][_k]-inFx3[i-1][_j][_k])/deltaXI+
					(inFy3[i-1][j][_k]-inFy3[i-1][_j][_k])/deltaET+
					(inFz3[i-1][_j][k]-inFz3[i-1][_j][_k])/deltaZT);
				Rf4 = -((inFx4[i][_j][_k]-inFx4[i-1][_j][_k])/deltaXI+
					(inFy4[i-1][j][_k]-inFy4[i-1][_j][_k])/deltaET+
					(inFz4[i-1][_j][k]-inFz4[i-1][_j][_k])/deltaZT);
				Rf5 = -((inFx5[i][_j][_k]-inFx5[i-1][_j][_k])/deltaXI+
					(inFy5[i-1][j][_k]-inFy5[i-1][_j][_k])/deltaET+
					(inFz5[i-1][_j][k]-inFz5[i-1][_j][_k])/deltaZT);

				Rk1 = Rp1+Rf1;
				Rk2 = Rp2+Rf2+vF2[i][j][k]+f/J[i][j][k];
				Rk3 = Rp3+Rf3+vF3[i][j][k];
				Rk4 = Rp4+Rf4+vF4[i][j][k];
				Rk5 = Rp5+Rf5+vF5[i][j][k]+f*u/J[i][j][k];

				MR1[i][j][k] = (d11*Rk1+d12*Rk2+d13*Rk3+d14*Rk4+d15*Rk5);
				MR2[i][j][k] = (d21*Rk1+d22*Rk2+d23*Rk3+d24*Rk4+d25*Rk5);
				MR3[i][j][k] = (d31*Rk1+d32*Rk2+d33*Rk3+d34*Rk4+d35*Rk5);
				MR4[i][j][k] = (d41*Rk1+d42*Rk2+d43*Rk3+d44*Rk4+d45*Rk5);
				MR5[i][j][k] = (d51*Rk1+d52*Rk2+d53*Rk3+d54*Rk4+d55*Rk5);
				
				/**** Hybrid LUSGS ****/

				xix = xidx[i][j][k];
				xiy = xidy[i][j][k];
				xiz = xidz[i][j][k];

				etx = etdx[i][j][k];
				ety = etdy[i][j][k];
				etz = etdz[i][j][k];

				ztx = ztdx[i][j][k];
				zty = ztdy[i][j][k];
				ztz = ztdz[i][j][k];


				U = xix*u+xiy*v+xiz*w;
				V = etx*u+ety*v+etz*w;
				W = ztx*u+zty*v+ztz*w;

				Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*(xix*xix+xiy*xiy+xiz*xiz));
				Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*(etx*etx+ety*ety+etz*etz));
				Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C*(ztx*ztx+zty*zty+ztz*ztz));

				lambdaX = ((beta+1)*fabs(U)+Sx)/2;
				lambdaY = ((beta+1)*fabs(V)+Sy)/2;
				lambdaZ = ((beta+1)*fabs(W)+Sz)/2;


				/*
				lambdaX = 0.5*((beta+1)*fabs(U)+Sx)+2*K*mu_L/Pr_L/rho/(xidx[i][j][k]/deltaXI+xidy[i][j][k]/deltaET+xidz[i][j][k]/deltaZT);
				lambdaY = 0.5*((beta+1)*fabs(V)+Sy)+2*K*mu_L/Pr_L/rho/(etdx[i][j][k]/deltaXI+etdy[i][j][k]/deltaET+etdz[i][j][k]/deltaZT);
				lambdaZ = 0.5*((beta+1)*fabs(W)+Sz)+2*K*mu_L/Pr_L/rho/(ztdx[i][j][k]/deltaXI+ztdy[i][j][k]/deltaET+ztdz[i][j][k]/deltaZT);
				*/

				lambda = lambdaX/deltaXI+lambdaY/deltaET+lambdaZ/deltaZT;

				d11 = 2*deltaT/(3*beta+2*lambda*deltaT);

				d22 = d33 = d44 = d55 = 2*deltaT/(3+2*lambda*deltaT);

				temp1 = 6*(K-1)*(beta-1)*deltaT;
				temp2 = K*(3+2*lambda*deltaT)*(3*beta+2*lambda*deltaT)*rho;
				temp = -temp1/temp2;

				d51 = temp;

				NL1[i][j][k] = d11*MR1[i][j][k];
				NL2[i][j][k] = d22*MR2[i][j][k];
				NL3[i][j][k] = d33*MR3[i][j][k];
				NL4[i][j][k] = d44*MR4[i][j][k];
				NL5[i][j][k] = d51*MR1[i][j][k]+d55*MR5[i][j][k];

			}
		}
	}
#pragma omp barrier





// ============================================ //
	for (int sweep = 1; sweep <= 4; sweep++) {  //
// ============================================ //
		
		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		MPI_Status istat[8];

		/**** implicit data communation ****/
		istart=3;
		iend = gend[myid];

		l_nbr = myid - 1;
		r_nbr = myid + 1;
		if(myid == 0) l_nbr=MPI_PROC_NULL;
		if(myid == nproc-1) r_nbr=MPI_PROC_NULL;


		icount = Y_m*Z_m;
		
		
		itag=460;
		MPI_Sendrecv((void *)&NL1[iend][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&NL1[istart-1][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
		itag=470;
		MPI_Sendrecv((void *)&NL2[iend][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&NL2[istart-1][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
		itag=480;
		MPI_Sendrecv((void *)&NL3[iend][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&NL3[istart-1][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
		itag=490;
		MPI_Sendrecv((void *)&NL4[iend][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&NL4[istart-1][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
		itag=500;
		MPI_Sendrecv((void *)&NL5[iend][0][0], icount, MPI_DOUBLE, r_nbr, itag, (void *)&NL5[istart-1][0][0], icount, MPI_DOUBLE, l_nbr, itag, comm, istat);
		
		
		itag=510;
		MPI_Sendrecv((void *)&NL1[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&NL1[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
		itag=520;
		MPI_Sendrecv((void *)&NL2[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&NL2[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
		itag=530;
		MPI_Sendrecv((void *)&NL3[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&NL3[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
		itag=540;
		MPI_Sendrecv((void *)&NL4[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&NL4[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
		itag=550;
		MPI_Sendrecv((void *)&NL5[istart][0][0], icount, MPI_DOUBLE, l_nbr, itag, (void *)&NL5[iend+1][0][0], icount, MPI_DOUBLE, r_nbr, itag, comm, istat);
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		/**** implicit data communation-end ****/
		
//// ========================= ////
		 istart = 3;		   ////	
//// ========================= ////
		iend = gend[myid];	   ////
//// ========================= ////

#pragma omp parallel for private(\
	rho,u,v,w,U,V,W,VV,P,C,T,H,\
	LL1,LL2,LL3,LL4,LL5,\
	UU1,UU2,UU3,UU4,UU5,\
	Sx,Sy,Sz,beta,\
	lambda,lambdaX,lambdaY,lambdaZ,\
	d11,d22,d33,d44,d51,d55,\
	temp,temp1,temp2,\
	Aplus11,Aplus12,Aplus13,Aplus14,\
	Aplus21,Aplus22,\
	Aplus31,        Aplus33,\
	Aplus41,          	    Aplus44,\
	Aplus51,Aplus52,Aplus53,Aplus54,Aplus55,\
	Bplus11,Bplus12,Bplus13,Bplus14,\
	Bplus21,Bplus22,\
	Bplus31,        Bplus33,\
	Bplus41,          	    Bplus44,\
	Bplus51,Bplus52,Bplus53,Bplus54,Bplus55,\
	Cplus11,Cplus12,Cplus13,Cplus14,\
	Cplus21,Cplus22,\
	Cplus31,        Cplus33,\
	Cplus41,          	    Cplus44,\
	Cplus51,Cplus52,Cplus53,Cplus54,Cplus55,\
	Amius11,Amius22,Amius33,Amius44,Amius55,\
	Bmius11,Bmius22,Bmius33,Bmius44,Bmius55,\
	Cmius11,Cmius22,Cmius33,Cmius44,Cmius55,\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
	j,k,j_,k_,_j,_k\
	)

	
		for (i = istart ; i <= iend; i++) { 
			for (j = 2, _j = 1; j < nyy; j++, _j++) {
				for (k = 2, _k = 1; k < nzz; k++, _k++) {



//// ========================================================================= ////
/**** lower part ****/														   ////
//// ========================================================================= ////




					rho = U1_[i-1][j][k]*J[i-1][j][k];
					u = U2_[i-1][j][k]/U1_[i-1][j][k];
					v = U3_[i-1][j][k]/U1_[i-1][j][k];
					w = U4_[i-1][j][k]/U1_[i-1][j][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i-1][j][k]*J[i-1][j][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho; /**** C = sqrt(K*P/rho); ****/
					T = P/rho;
					H = 0.5*VV+C/(K-1);
					/*
					if (VV/C <= e)
						beta = e;
					else
						beta = VV/C;
					*/

					beta = max(VV/C,e);


					xix = xidx[i-1][j][k];
					xiy = xidy[i-1][j][k];
					xiz = xidz[i-1][j][k];


					U = xix*u+xiy*v+xiz*w;

					Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*(xix*xix+xiy*xiy+xiz*xiz));

					lambdaX = ((beta+1)*fabs(U)+Sx)/2;

					d11 = d22 = d33 = d44 = d55 = fabs(lambdaX);
					
					
					Aplus11 = 0.5*(beta*U+d11);
					Aplus12 = 0.5*P*K*beta*xix;
					Aplus13 = 0.5*P*K*beta*xiy;
					Aplus14 = 0.5*P*K*beta*xiz;

					Aplus21 = 0.5*xix/rho;
					Aplus22 = 0.5*(U+d22);

					Aplus31 = 0.5*xiy/rho;
					Aplus33 = 0.5*(U+d33);

					Aplus41 = 0.5*xiz/rho;
					Aplus44 = 0.5*(U+d44);

					Aplus51 = 0.5*U*(K-1)*(beta-1)/rho/K;
					Aplus52 = 0.5*xix*beta*T*(K-1);
					Aplus53 = 0.5*xiy*beta*T*(K-1);
					Aplus54 = 0.5*xiz*beta*T*(K-1);
					Aplus55 = 0.5*(U+d55);
					
					
					

					rho = U1_[i][j-1][k]*J[i][j-1][k];
					u = U2_[i][j-1][k]/U1_[i][j-1][k];
					v = U3_[i][j-1][k]/U1_[i][j-1][k];
					w = U4_[i][j-1][k]/U1_[i][j-1][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j-1][k]*J[i][j-1][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);
					
					beta = max(VV/C,e);


					etx = etdx[i][j-1][k];
					ety = etdy[i][j-1][k];
					etz = etdz[i][j-1][k];

					V = etx*u+ety*v+etz*w;

					Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*(etx*etx+ety*ety+etz*etz));

					lambdaY = ((beta+1)*fabs(V)+Sy)/2;

					d11 = d22 = d33 = d44 = d55 = fabs(lambdaY);
					
					/* flux splitting Bplus*/
					Bplus11 = 0.5*(beta*V+d11);
					Bplus12 = 0.5*P*K*beta*etx;
					Bplus13 = 0.5*P*K*beta*ety;
					Bplus14 = 0.5*P*K*beta*etz;

					Bplus21 = 0.5*etx/rho;
					Bplus22 = 0.5*(V+d22);

					Bplus31 = 0.5*ety/rho;
					Bplus33 = 0.5*(V+d33);

					Bplus41 = 0.5*etz/rho;
					Bplus44 = 0.5*(V+d44);

					Bplus51 = 0.5*V*(K-1)*(beta-1)/rho/K;
					Bplus52 = 0.5*beta*T*(K-1)*etx;
					Bplus53 = 0.5*beta*T*(K-1)*ety;
					Bplus54 = 0.5*beta*T*(K-1)*etz;
					Bplus55 = 0.5*(V+d55);



					rho = U1_[i][j][k-1]*J[i][j][k-1];
					u = U2_[i][j][k-1]/U1_[i][j][k-1];
					v = U3_[i][j][k-1]/U1_[i][j][k-1];
					w = U4_[i][j][k-1]/U1_[i][j][k-1];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j][k-1]*J[i][j][k-1]-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);


					beta = max(VV/C,e);


					ztx = ztdx[i][j][k-1]; 
					zty = ztdy[i][j][k-1]; 
					ztz = ztdz[i][j][k-1]; 

					W = ztx*u+zty*v+ztz*w;

					Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C*(ztx*ztx+zty*zty+ztz*ztz));

					lambdaZ = ((beta+1)*fabs(W)+Sz)/2;

					d11 = d22 = d33 = d44 = d55 = fabs(lambdaZ);
					

					Cplus11 = 0.5*(beta*W+d11);
					Cplus12 = 0.5*P*K*beta*ztx;
					Cplus13 = 0.5*P*K*beta*zty;
					Cplus14 = 0.5*P*K*beta*ztz;

					Cplus21 = 0.5*ztx/rho;
					Cplus22 = 0.5*(W+d22);

					Cplus31 = 0.5*zty/rho;
					Cplus33 = 0.5*(W+d33);

					Cplus41 = 0.5*ztz/rho;
					Cplus44 = 0.5*(W+d44);

					Cplus51 = 0.5*W*(K-1)*(beta-1)/rho/K;
					Cplus52 = 0.5*beta*T*(K-1)*ztx;
					Cplus53 = 0.5*beta*T*(K-1)*zty;
					Cplus54 = 0.5*beta*T*(K-1)*ztz;
					Cplus55 = 0.5*(W+d55);
					
					LL1 = (-Aplus11*NL1[i-1][j][k]
					       -Aplus12*NL2[i-1][j][k]
						   -Aplus13*NL3[i-1][j][k]
					       -Aplus14*NL4[i-1][j][k])/deltaXI+

						 (-Bplus11*NL1[i][_j][k]
						  -Bplus12*NL2[i][_j][k]
						  -Bplus13*NL3[i][_j][k]
						  -Bplus14*NL4[i][_j][k])/deltaET+

						 (-Cplus11*NL1[i][j][_k]
						  -Cplus12*NL2[i][j][_k]
						  -Cplus13*NL3[i][j][_k]
						  -Cplus14*NL4[i][j][_k])/deltaZT;
						  
				
						  
					LL2 = (-Aplus21*NL1[i-1][j][k]
						   -Aplus22*NL2[i-1][j][k])/deltaXI+

						  (-Bplus21*NL1[i][_j][k]
						   -Bplus22*NL2[i][_j][k])/deltaET+

						  (-Cplus21*NL1[i][j][_k]
					       -Cplus22*NL2[i][j][_k])/deltaZT;

					LL3 = (-Aplus31*NL1[i-1][j][k]
						   -Aplus33*NL3[i-1][j][k])/deltaXI+

						 (-Bplus31*NL1[i][_j][k]
					      -Bplus33*NL3[i][_j][k])/deltaET+

						 (-Cplus31*NL1[i][j][_k]
				      	  -Cplus33*NL3[i][j][_k])/deltaZT;

					LL4 = (-Aplus41*NL1[i-1][j][k]
					       -Aplus44*NL4[i-1][j][k])/deltaXI+

						  (-Bplus41*NL1[i][_j][k]
					       -Bplus44*NL4[i][_j][k])/deltaET+

						  (-Cplus41*NL1[i][j][_k]
					       -Cplus44*NL4[i][j][_k])/deltaZT;

					LL5 = (-Aplus51*NL1[i-1][j][k]
						   -Aplus52*NL2[i-1][j][k]
						   -Aplus53*NL3[i-1][j][k]
						   -Aplus54*NL4[i-1][j][k]
						   -Aplus55*NL5[i-1][j][k])/deltaXI+

						  (-Bplus51*NL1[i][_j][k]
						   -Bplus52*NL2[i][_j][k]
						   -Bplus53*NL3[i][_j][k]
						   -Bplus54*NL4[i][_j][k]
						   -Bplus55*NL5[i][_j][k])/deltaET+

						  (-Cplus51*NL1[i][j][_k]
						   -Cplus52*NL2[i][j][_k]
						   -Cplus53*NL3[i][j][_k]
						   -Cplus54*NL4[i][j][_k]
						   -Cplus55*NL5[i][j][_k])/deltaZT;

		   
//// ========================================================================= ////
/**** lower part-end/     													   ////
//// ========================================================================= ////




//// ========================================================================= ////
/**** upper part ****/														   ////
//// ========================================================================= ////


						   rho = U1_[i+1][j][k]*J[i+1][j][k];
					u = U2_[i+1][j][k]/U1_[i+1][j][k];
					v = U3_[i+1][j][k]/U1_[i+1][j][k];
					w = U4_[i+1][j][k]/U1_[i+1][j][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i+1][j][k]*J[i+1][j][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho; /**** C = sqrt(K*P/rho); ****/
					T = P/rho;
					H = 0.5*VV+C/(K-1);
					/*
					if (VV/C <= e)
						beta = e;
					else
						beta = VV/C;
					*/

					beta = max(VV/C,e);


					xix = xidx[i+1][j][k];
					xiy = xidy[i+1][j][k];
					xiz = xidz[i+1][j][k];


					U = xix*u+xiy*v+xiz*w;

					Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*(xix*xix+xiy*xiy+xiz*xiz));

					lambdaX = ((beta+1)*fabs(U)+Sx)/2;

					d11 = d22 = d33 = d44 = d55 =  fabs(lambdaX);
					
					
					Aplus11 = 0.5*(beta*U+d11);
					Aplus12 = 0.5*P*K*beta*xix;
					Aplus13 = 0.5*P*K*beta*xiy;
					Aplus14 = 0.5*P*K*beta*xiz;

					Aplus21 = 0.5*xix/rho;
					Aplus22 = 0.5*(U+d22);

					Aplus31 = 0.5*xiy/rho;
					Aplus33 = 0.5*(U+d33);

					Aplus41 = 0.5*xiz/rho;
					Aplus44 = 0.5*(U+d44);

					Aplus51 = 0.5*U*(K-1)*(beta-1)/rho/K;
					Aplus52 = 0.5*xix*beta*T*(K-1);
					Aplus53 = 0.5*xiy*beta*T*(K-1);
					Aplus54 = 0.5*xiz*beta*T*(K-1);
					Aplus55 = 0.5*(U+d55);

					/* flux splitting mius*/
					Amius11 = 0.5*(beta*U-d11);
					Amius22 = 0.5*(U-d22);
					Amius33 = 0.5*(U-d33);
					Amius44 = 0.5*(U-d44);
					Amius55 = 0.5*(U-d55);


					
					rho = U1_[i][j+1][k]*J[i][j+1][k];
					u = U2_[i][j+1][k]/U1_[i][j+1][k];
					v = U3_[i][j+1][k]/U1_[i][j+1][k];
					w = U4_[i][j+1][k]/U1_[i][j+1][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j+1][k]*J[i][j+1][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);
					
					beta = max(VV/C,e);


					etx = etdx[i][j+1][k];
					ety = etdy[i][j+1][k];
					etz = etdz[i][j+1][k];

					V = etx*u+ety*v+etz*w;

					Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*(etx*etx+ety*ety+etz*etz));

					lambdaY = ((beta+1)*fabs(V)+Sy)/2;

					d11 = d22 = d33 = d44 = d55 =  fabs(lambdaY);
					

					/* flux splitting Bplus*/
					Bplus11 = 0.5*(beta*V+d11);
					Bplus12 = 0.5*P*K*beta*etx;
					Bplus13 = 0.5*P*K*beta*ety;
					Bplus14 = 0.5*P*K*beta*etz;

					Bplus21 = 0.5*etx/rho;
					Bplus22 = 0.5*(V+d22);

					Bplus31 = 0.5*ety/rho;
					Bplus33 = 0.5*(V+d33);

					Bplus41 = 0.5*etz/rho;
					Bplus44 = 0.5*(V+d44);

					Bplus51 = 0.5*V*(K-1)*(beta-1)/rho/K;
					Bplus52 = 0.5*beta*T*(K-1)*etx;
					Bplus53 = 0.5*beta*T*(K-1)*ety;
					Bplus54 = 0.5*beta*T*(K-1)*etz;
					Bplus55 = 0.5*(V+d55);

					
					/* flux splitting mius*/
					Bmius11 = 0.5*(beta*V-d11);
					Bmius22 = 0.5*(V-d22);
					Bmius33 = 0.5*(V-d33);
					Bmius44 = 0.5*(V-d44);
					Bmius55 = 0.5*(V-d55);



					rho = U1_[i][j][k+1]*J[i][j][k+1];
					u = U2_[i][j][k+1]/U1_[i][j][k+1];
					v = U3_[i][j][k+1]/U1_[i][j][k+1];
					w = U4_[i][j][k+1]/U1_[i][j][k+1];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j][k+1]*J[i][j][k+1]-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);


					beta = max(VV/C,e);

					
					ztx = ztdx[i][j][k+1]; 
					zty = ztdy[i][j][k+1]; 
					ztz = ztdz[i][j][k+1]; 

					W = ztx*u+zty*v+ztz*w;

					Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C*(ztx*ztx+zty*zty+ztz*ztz));

					lambdaZ = ((beta+1)*fabs(W)+Sz)/2;

					d11 = d22 = d33 = d44 = d55 =  fabs(lambdaZ);
					

					Cplus11 = 0.5*(beta*W+d11);
					Cplus12 = 0.5*P*K*beta*ztx;
					Cplus13 = 0.5*P*K*beta*zty;
					Cplus14 = 0.5*P*K*beta*ztz;

					Cplus21 = 0.5*ztx/rho;
					Cplus22 = 0.5*(W+d22);

					Cplus31 = 0.5*zty/rho;
					Cplus33 = 0.5*(W+d33);

					Cplus41 = 0.5*ztz/rho;
					Cplus44 = 0.5*(W+d44);

					Cplus51 = 0.5*W*(K-1)*(beta-1)/rho/K;
					Cplus52 = 0.5*beta*T*(K-1)*ztx;
					Cplus53 = 0.5*beta*T*(K-1)*zty;
					Cplus54 = 0.5*beta*T*(K-1)*ztz;
					Cplus55 = 0.5*(W+d55);


					
					/* flux splitting mius*/
					Cmius11 = 0.5*(beta*W-d11);
					Cmius22 = 0.5*(W-d22);
					Cmius33 = 0.5*(W-d33);
					Cmius44 = 0.5*(W-d44);
					Cmius55 = 0.5*(W-d55);




					UU1 = (+Amius11*NL1[i+1][j][k]
						   +Aplus12*NL2[i+1][j][k]
						   +Aplus13*NL3[i+1][j][k]
					       +Aplus14*NL4[i+1][j][k])/deltaXI+

						  (+Bmius11*NL1[i][j_][k]
					       +Bplus12*NL2[i][j_][k]
					       +Bplus13*NL3[i][j_][k]
					       +Bplus14*NL4[i][j_][k])/deltaET+

						  (+Cmius11*NL1[i][j][k_]
					       +Cplus12*NL2[i][j][k_]
					       +Cplus13*NL3[i][j][k_]
					       +Cplus14*NL4[i][j][k_])/deltaZT;
						   

					UU2 = (+Aplus21*NL1[i+1][j][k]
						   +Amius22*NL2[i+1][j][k])/deltaXI+

						  (+Bplus21*NL1[i][j_][k]
					       +Bmius22*NL2[i][j_][k])/deltaET+

						  (+Cplus21*NL1[i][j][k_]
					       +Cmius22*NL2[i][j][k_])/deltaZT;


					UU3 = (+Aplus31*NL1[i+1][j][k]
					       +Amius33*NL3[i+1][j][k])/deltaXI+

						  (+Bplus31*NL1[i][j_][k]
					       +Bmius33*NL3[i][j_][k])/deltaET+

						  (+Cplus31*NL1[i][j][k_]
					       +Cmius33*NL3[i][j][k_])/deltaZT;


					UU4 = (+Aplus41*NL1[i+1][j][k]
					       +Amius44*NL4[i+1][j][k])/deltaXI+

						  (+Bplus41*NL1[i][j_][k]
				      	   +Bmius44*NL4[i][j_][k])/deltaET+

						  (+Cplus41*NL1[i][j][k_]
					       +Cmius44*NL4[i][j][k_])/deltaZT;

					UU5 = (+Aplus51*NL1[i+1][j][k]
					       +Aplus52*NL2[i+1][j][k]
					       +Aplus53*NL3[i+1][j][k]
					       +Aplus54*NL4[i+1][j][k]
					       +Amius55*NL5[i+1][j][k])/deltaXI+

						  (+Bplus51*NL1[i][j_][k]
					       +Bplus52*NL2[i][j_][k]
					       +Bplus53*NL3[i][j_][k]
					       +Bplus54*NL4[i][j_][k]
					       +Bmius55*NL5[i][j_][k])/deltaET+

						  (+Cplus51*NL1[i][j][k_]
					       +Cplus52*NL2[i][j][k_]
					       +Cplus53*NL3[i][j][k_]
					       +Cplus54*NL4[i][j][k_]
					       +Cmius55*NL5[i][j][k_])/deltaZT;


						   
//// ========================================================================= ////
/**** upper part-end ****/   												   ////
//// ========================================================================= ////


					rho = U1_[i][j][k]*J[i][j][k];
					u = U2_[i][j][k]/U1_[i][j][k];
					v = U3_[i][j][k]/U1_[i][j][k];
					w = U4_[i][j][k]/U1_[i][j][k];     
					VV = u*u+v*v+w*w;
					P = (U5_[i][j][k]*J[i][j][k]-0.5*rho*VV)*(K-1);
					C = K*P/rho;
					T = P/rho;
					H = 0.5*VV+C/(K-1);
					
					beta = max(VV/C,e);

					xix = xidx[i][j][k];
					xiy = xidy[i][j][k];
					xiz = xidz[i][j][k];

					etx = etdx[i][j][k];
					ety = etdy[i][j][k];
					etz = etdz[i][j][k];

					ztx = ztdx[i][j][k];
					zty = ztdy[i][j][k];
					ztz = ztdz[i][j][k];

					U = xix*u+xiy*v+xiz*w;
					V = etx*u+ety*v+etz*w;
					W = ztx*u+zty*v+ztz*w;

					Sx = sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*(xix*xix+xiy*xiy+xiz*xiz));
					Sy = sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*(etx*etx+ety*ety+etz*etz));
					Sz = sqrt(W*W*(beta-1)*(beta-1)+4*beta*C*(ztx*ztx+zty*zty+ztz*ztz));

					lambdaX = ((beta+1)*fabs(U)+Sx)/2;
					lambdaY = ((beta+1)*fabs(V)+Sy)/2;
					lambdaZ = ((beta+1)*fabs(W)+Sz)/2;


					/*
					lambdaX = 0.5*((beta+1)*fabs(U)+Sx)+2*K*mu_L/Pr_L/rho/(xidx[i][j][k]/deltaXI+xidy[i][j][k]/deltaET+xidz[i][j][k]/deltaZT);
					lambdaY = 0.5*((beta+1)*fabs(V)+Sy)+2*K*mu_L/Pr_L/rho/(etdx[i][j][k]/deltaXI+etdy[i][j][k]/deltaET+etdz[i][j][k]/deltaZT);
					lambdaZ = 0.5*((beta+1)*fabs(W)+Sz)+2*K*mu_L/Pr_L/rho/(ztdx[i][j][k]/deltaXI+ztdy[i][j][k]/deltaET+ztdz[i][j][k]/deltaZT);
					*/

					lambda = lambdaX/deltaXI+lambdaY/deltaET+lambdaZ/deltaZT;

					d11 = 2*deltaT/(3*beta+2*lambda*deltaT);

					d22 = d33 =  d44 = d55 = 2*deltaT/(3+2*lambda*deltaT);

					temp1 = 6*(K-1)*(beta-1)*deltaT;
					temp2 = K*(3+2*lambda*deltaT)*(3*beta+2*lambda*deltaT)*rho;
					temp = -temp1/temp2;

					d51 = temp;
					
					 NR1[i][j][k] = d11*(MR1[i][j][k]-LL1-UU1);
					 NR2[i][j][k] = d22*(MR2[i][j][k]-LL2-UU2);
					 NR3[i][j][k] = d33*(MR3[i][j][k]-LL3-UU3);
					 NR4[i][j][k] = d44*(MR4[i][j][k]-LL4-UU4);
					 NR5[i][j][k] = d51*(MR1[i][j][k]-LL1-UU1)+d55*(MR5[i][j][k]-LL5-UU5);
					
				}
			}
		}

		
//// ============================================ ////
		 istart = 3;		    		  ////	
//// ============================================ ////
			iend = gend[myid];		    		  ////
//// ============================================ ////

#pragma omp parallel for private(j,k)
		for (i = istart; i <= iend; i++) {
			for (j = 2; j <= ny; j++) { 
				for (k = 2; k <= nz; k++) { 

					NL1[i][j][k] = NR1[i][j][k];
					NL2[i][j][k] = NR2[i][j][k];
					NL3[i][j][k] = NR3[i][j][k];
					NL4[i][j][k] = NR4[i][j][k];
					NL5[i][j][k] = NR5[i][j][k];
					
				}
			}
		}

// ============================================ //
		}  /**** sweep-end ****/				//
// ============================================ //
		

		
//// ============================== ////
		 istart = 3;		        ////	
//// ============================== ////
			iend = gend[myid];	    ////
//// ============================== ////

#pragma omp parallel for private(j, k, rhoold, P, U, V, W, T, VV, rho, Uold, Vold, Wold, VVold, Pold, Told)

	for (i = istart; i <= iend; i++) {
		for (j = 2; j < nyy; j++) {
			for (k = 2; k < nzz; k++) {     
				
				/* flux parameter */
				rhoold = U1_[i][j][k]*J[i][j][k];
				Uold = U2_[i][j][k]/U1_[i][j][k];
				Vold = U3_[i][j][k]/U1_[i][j][k];
				Wold = U4_[i][j][k]/U1_[i][j][k];
				VVold = Uold*Uold+Vold*Vold+Wold*Wold;
				Pold = (U5_[i][j][k]*J[i][j][k]-0.5*rhoold*VVold)*(K-1);
				Told = Pold/rhoold;
				
				
				P = Pold+NR1[i][j][k]*J[i][j][k];
				U = Uold+NR2[i][j][k]*J[i][j][k];
				V = Vold+NR3[i][j][k]*J[i][j][k];
				W = Wold+NR4[i][j][k]*J[i][j][k];
				T = Told+NR5[i][j][k]*J[i][j][k];
				

				rho = P/T;
				VV = U*U+V*V+W*W;
				U1_[i][j][k] = rho/J[i][j][k];
				U2_[i][j][k] = rho*U/J[i][j][k];
				U3_[i][j][k] = rho*V/J[i][j][k];
				U4_[i][j][k] = rho*W/J[i][j][k];
				U5_[i][j][k] = (P/(K-1)+0.5*rho*(U*U+V*V+W*W))/J[i][j][k];
				

				Residual1[i][j][k] = (P-Pold)*(P-Pold);
				Residual2[i][j][k] = (U-Uold)*(U-Uold);
				Residual3[i][j][k] = (V-Vold)*(V-Vold);
				Residual4[i][j][k] = (W-Wold)*(W-Wold);
				Residual5[i][j][k] = (T-Told)*(T-Told);
				
				
			}
		}
	}

//// ======================== ////
		 istart = 3;		  ////	
//// ======================== ////
		iend = gend[myid];	  ////
//// ======================== ////

	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
				
				e1 = e1+Residual1[i][j][k]*Residual1[i][j][k];
				e2 = e2+Residual2[i][j][k]*Residual2[i][j][k];
				e3 = e3+Residual3[i][j][k]*Residual3[i][j][k];
				e4 = e4+Residual4[i][j][k]*Residual4[i][j][k];
				e5 = e5+Residual5[i][j][k]*Residual5[i][j][k];
					
			}
		}
	}
	
	e1 = sqrt(e1)*DN;
	e2 = sqrt(e2)*DN;
	e3 = sqrt(e3)*DN;
	e4 = sqrt(e4)*DN;
	e5 = sqrt(e5)*DN;
		
	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	
	MPI_Allreduce ((void*)&e1, (void*)&*er1, 1, MPI_DOUBLE, MPI_SUM, comm );
	MPI_Allreduce ((void*)&e2, (void*)&*er2, 1, MPI_DOUBLE, MPI_SUM, comm );
	MPI_Allreduce ((void*)&e3, (void*)&*er3, 1, MPI_DOUBLE, MPI_SUM, comm );
	MPI_Allreduce ((void*)&e4, (void*)&*er4, 1, MPI_DOUBLE, MPI_SUM, comm );
	MPI_Allreduce ((void*)&e5, (void*)&*er5, 1, MPI_DOUBLE, MPI_SUM, comm );
		
}
// ============================================================================ //
