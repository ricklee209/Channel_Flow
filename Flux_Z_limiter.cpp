



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Flux_Z_limiter
(
// ============================================================================ //
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*J_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFz1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFz5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpZ)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

#include "ijk.h"
#include "prm.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	double rho,U,V,W,VV,P,C,T,h,H;
	double u,v,w;
	double temp,temp1,temp2,temp3,temp4,temp5;
	double beta,S,_S_;

	double m11,m12,m13,m14,m15,m22,m23,m24,m25,m33,m34,m35,m44,m45,m55;
	double thedac1,thedac2,thedac3,thedad1,thedad2,thedad3;

	
	
	double u1ii,u2ii,u3ii,u4ii,u5ii;
	double u1i,u2i,u3i,u4i,u5i;
	double u1,u2,u3,u4,u5;
	double iu1,iu2,iu3,iu4,iu5;
	double iiu1,iiu2,iiu3,iiu4,iiu5;

	double irL,rL,rLi,bL;
	double irR,rR,rRi,bR;
	double minL3,minR3;
	double jL,jR;




	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2;
	double d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35;
	double d41,d42,d43,d44,d45,d51,d52,d53,d54,d55, Fav1,Fav2,Fav3,Fav4,Fav5;

/**** MUSCL 5th-order ****/
// ============================================================================ //
	



//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;						  ////	
//// ============================================ ////
		iend = gend[myid];						  ////
//// ============================================ ////


	for (i = istart; i <= iend; i++) {
#pragma omp parallel for private(k,\
		u1ii,u2ii,u3ii,u4ii,u5ii,u1i,u2i,u3i,u4i,u5i,\
		u1,u2,u3,u4,u5,\
		iu1,iu2,iu3,iu4,iu5,iiu1,iiu2,iiu3,iiu4,iiu5,\
		irL,rL,rLi,bL,irR,rR,rRi,bR,minL3,minR3,jL,jR)

		for ( j = 2; j <= ny; j++) {
			for (k = 3; k <= nz-1; k++) {

				if (k <= 2 | k >= ny) {

					jL = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

					jR = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);

					iu1 = U1_[i][j][k-1];
					iu2 = U2_[i][j][k-1];
					iu3 = U3_[i][j][k-1];
					iu4 = U4_[i][j][k-1];
					iu5 = U5_[i][j][k-1];

					u1 = U1_[i][j][k];
					u2 = U2_[i][j][k];
					u3 = U3_[i][j][k];
					u4 = U4_[i][j][k];
					u5 = U5_[i][j][k];

					u1i = U1_[i][j][k+1];
					u2i = U2_[i][j][k+1];
					u3i = U3_[i][j][k+1];
					u4i = U4_[i][j][k+1];
					u5i = U5_[i][j][k+1];

					rL = (u1i - u1)/( u1  - iu1);
					rR = (u1 - iu1)/( u1i - u1 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML1[i-1][j-1][k  ] = (u1-0.5*max(0,minL3)*(u1  -iu1))*jL;
					MR1[i-1][j-1][k-1] = (u1-0.5*max(0,minR3)*(u1i - u1))*jR;


					rL = (u2i - u2)/( u2  - iu2);
					rR = (u2 - iu2)/( u2i - u2 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML2[i-1][j-1][k  ] = (u2-0.5*max(0,minL3)*(u2  -iu2))*jL;
					MR2[i-1][j-1][k-1] = (u2-0.5*max(0,minR3)*(u2i - u2))*jR;


					rL = (u3i - u3)/( u3  - iu3);
					rR = (u3 - iu3)/( u3i - u3 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML3[i-1][j-1][k  ] = (u3-0.5*max(0,minL3)*(u3  -iu3))*jL;
					MR3[i-1][j-1][k-1] = (u3-0.5*max(0,minR3)*(u3i - u3))*jR;


					rL = (u4i - u4)/( u4  - iu4);
					rR = (u4 - iu4)/( u4i - u4 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML4[i-1][j-1][k  ] = (u4-0.5*max(0,minL3)*(u4  -iu4))*jL;
					MR4[i-1][j-1][k-1] = (u4-0.5*max(0,minR3)*(u4i - u4))*jR;


					rL = (u5i - u5)/( u5  - iu5);
					rR = (u5 - iu5)/( u5i - u5 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML5[i-1][j-1][k  ] = (u5-0.5*max(0,minL3)*(u5  -iu5))*jL;
					MR5[i-1][j-1][k-1] = (u5-0.5*max(0,minR3)*(u5i - u5))*jR;

				}

				else {

					iiu1 = U1_[i][j][k-2];
					iiu2 = U2_[i][j][k-2];
					iiu3 = U3_[i][j][k-2];
					iiu4 = U4_[i][j][k-2];
					iiu5 = U5_[i][j][k-2];

					iu1 = U1_[i][j][k-1];
					iu2 = U2_[i][j][k-1];
					iu3 = U3_[i][j][k-1];
					iu4 = U4_[i][j][k-1];
					iu5 = U5_[i][j][k-1];

					u1 = U1_[i][j][k];
					u2 = U2_[i][j][k];
					u3 = U3_[i][j][k];
					u4 = U4_[i][j][k];
					u5 = U5_[i][j][k];

					u1i = U1_[i][j][k+1];
					u2i = U2_[i][j][k+1];
					u3i = U3_[i][j][k+1];
					u4i = U4_[i][j][k+1];
					u5i = U5_[i][j][k+1];

					u1ii = U1_[i][j][k+2];
					u2ii = U2_[i][j][k+2];
					u3ii = U3_[i][j][k+2];
					u4ii = U4_[i][j][k+2];
					u5ii = U5_[i][j][k+2];


					jL = 1./(2/J[i][j][k-2]-13/J[i][j][k-1]+47/J[i][j][k]+27/J[i][j][k+1]-3/J[i][j][k+2]);
					jR = 1./(-3/J[i][j][k-2]+27/J[i][j][k-1]+47/J[i][j][k]-13/J[i][j][k+1]+2/J[i][j][k+2]);


					irL  = (u1   - iu1 )/(iu1  - iiu1);
					rL  = (u1i  -  u1 )/( u1  -  iu1);
					rLi = (u1ii -  u1i)/( u1i -   u1); 

					rRi = (u1i  -  u1 )/(u1ii -  u1i);
					rR  = ( u1  -  iu1)/( u1i -   u1);
					irR  = (iu1  - iiu1)/( u1   - iu1);

					bL = 1./30*(-2./irL+11+24*rL-3*rL*rLi);
					bR = 1./30*(-2./rRi+11+24*rR-3*rR*irR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML1[i-1][j-1][k  ] = (u1-0.5*max(0,minL3)*(u1  -iu1))*jL;
					MR1[i-1][j-1][k-1] = (u1-0.5*max(0,minR3)*(u1i - u1))*jR;


					irL  = (u2   - iu2 )/(iu2  - iiu2);
					rL  = (u2i  -  u2 )/( u2  -  iu2);
					rLi = (u2ii -  u2i)/( u2i -   u2); 

					rRi = (u2i  -  u2 )/(u2ii -  u2i);
					rR  = ( u2  -  iu2)/( u2i -   u2);
					irR  = (iu2  - iiu2)/( u2   - iu2);

					bL = 1./30*(-2./irL+11+24*rL-3*rL*rLi);
					bR = 1./30*(-2./rRi+11+24*rR-3*rR*irR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML2[i-1][j-1][k  ] = (u2-0.5*max(0,minL3)*(u2  -iu2))*jL;
					MR2[i-1][j-1][k-1] = (u2-0.5*max(0,minR3)*(u2i - u2))*jR;


					irL  = (u3   - iu3 )/(iu3  - iiu3);
					rL  = (u3i  -  u3 )/( u3  -  iu3);
					rLi = (u3ii -  u3i)/( u3i -   u3); 

					rRi = (u3i  -  u3 )/(u3ii -  u3i);
					rR  = ( u3  -  iu3)/( u3i -   u3);
					irR  = (iu3  - iiu3)/( u3   - iu3);

					bL = 1./30*(-2./irL+11+24*rL-3*rL*rLi);
					bR = 1./30*(-2./rRi+11+24*rR-3*rR*irR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML3[i-1][j-1][k  ] = (u3-0.5*max(0,minL3)*(u3  -iu3))*jL;
					MR3[i-1][j-1][k-1] = (u3-0.5*max(0,minR3)*(u3i - u3))*jR;


					irL  = (u4   - iu4 )/(iu4  - iiu4);
					rL  = (u4i  -  u4 )/( u4  -  iu4);
					rLi = (u4ii -  u4i)/( u4i -   u4); 

					rRi = (u4i  -  u4 )/(u4ii -  u4i);
					rR  = ( u4  -  iu4)/( u4i -   u4);
					irR  = (iu4  - iiu4)/( u4   - iu4);

					bL = 1./30*(-2./irL+11+24*rL-3*rL*rLi);
					bR = 1./30*(-2./rRi+11+24*rR-3*rR*irR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML4[i-1][j-1][k  ] = (u4-0.5*max(0,minL3)*(u4  -iu4))*jL;
					MR4[i-1][j-1][k-1] = (u4-0.5*max(0,minR3)*(u4i - u4))*jR;


					irL  = (u5   - iu5 )/(iu5  - iiu5);
					rL  = (u5i  -  u5 )/( u5  -  iu5);
					rLi = (u5ii -  u5i)/( u5i -   u5); 

					rRi = (u5i  -  u5 )/(u5ii -  u5i);
					rR  = ( u5  -  iu5)/( u5i -   u5);
					irR  = (iu5  - iiu5)/( u5   - iu5);

					bL = 1./30*(-2./irL+11+24*rL-3*rL*rLi);
					bR = 1./30*(-2./rRi+11+24*rR-3*rR*irR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML5[i-1][j-1][k  ] = (u5-0.5*max(0,minL3)*(u5  -iu5))*jL;
					MR5[i-1][j-1][k-1] = (u5-0.5*max(0,minR3)*(u5i - u5))*jR;


				}
			
				
			}
		}
	}
	
	for (i = istart; i <= iend; i++) {
		for (j = 2; j <= ny; j++) {
		
			k = 1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
			
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
				//if (myid == 0 & i == 3 & j == 24)
						//printf("%d\t%f\t%f\t%f\t%f\n",k,ML1[i-1][j-1][k],MR1[i-1][j-1][k],ML2[i-1][j-1][k],MR2[i-1][j-1][k]);
						
			

// ====================================================================================//
		
			/*
			k = 2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
					
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			*/
			
			k = 2; 
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
			
// ====================================================================================//		
			/*
			k = 3; 
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
			k = nz-1;
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			*/
			
			
			
// ====================================================================================//


			/*
			k = nz;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k-1]);
					
			ML1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k-1]);
			ML2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k-1]);
			ML3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k-1]);
			ML4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k-1]);
			ML5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k-1]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k-1]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j][k-1]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j][k-1]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j][k-1]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j][k-1]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j][k-1]);
			*/
			
			k = nz;
			
			temp = 1./(-1/J[i][j][k-1]+5/J[i][j][k]+2/J[i][j][k+1]);

			ML1[i-1][j-1][k] = temp*(-U1_[i][j][k-1]+5*U1_[i][j][k]+2*U1_[i][j][k+1]);
			ML2[i-1][j-1][k] = temp*(-U2_[i][j][k-1]+5*U2_[i][j][k]+2*U2_[i][j][k+1]);
			ML3[i-1][j-1][k] = temp*(-U3_[i][j][k-1]+5*U3_[i][j][k]+2*U3_[i][j][k+1]);
			ML4[i-1][j-1][k] = temp*(-U4_[i][j][k-1]+5*U4_[i][j][k]+2*U4_[i][j][k+1]);
			ML5[i-1][j-1][k] = temp*(-U5_[i][j][k-1]+5*U5_[i][j][k]+2*U5_[i][j][k+1]);
						
			temp = 1./(2/J[i][j][k-1]+5/J[i][j][k]-1/J[i][j][k+1]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j][k-1]+5*U1_[i][j][k]-U1_[i][j][k+1]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j][k-1]+5*U2_[i][j][k]-U2_[i][j][k+1]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j][k-1]+5*U3_[i][j][k]-U3_[i][j][k+1]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j][k-1]+5*U4_[i][j][k]-U4_[i][j][k+1]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j][k-1]+5*U5_[i][j][k]-U5_[i][j][k+1]);
			
			
// ====================================================================================//

			k = nz;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j][k+1]);
				
			MR1[i-1][j-1][k] = temp*(U1_[i][j][k]+U1_[i][j][k+1]);
			MR2[i-1][j-1][k] = temp*(U2_[i][j][k]+U2_[i][j][k+1]);
			MR3[i-1][j-1][k] = temp*(U3_[i][j][k]+U3_[i][j][k+1]);
			MR4[i-1][j-1][k] = temp*(U4_[i][j][k]+U4_[i][j][k+1]);
			MR5[i-1][j-1][k] = temp*(U5_[i][j][k]+U5_[i][j][k+1]);
			
			
				//if (myid == 0 & i == 3 & j == 24)
						//printf("%d\t%f\t%f\t%f\t%f\n",k,ML1[i-1][j-1][k],MR1[i-1][j-1][k],ML2[i-1][j-1][k],MR2[i-1][j-1][k]);
						
			
				
		}
	}
	
	


	

//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ZTx,ZTy,ZTz,\
	_rho,_u,_v,_w,_U,_V,_W,__W,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,W__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_W_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	j,k\
	)

		/*---Z fluxes---*/
		for (i = istart; i <= iend; i++) {
			for (j = 1; j < ny; j++) {
				for (k = 1; k < nzz; k++) {

				xix=xidx_w[i][j][k];
				xiy=xidy_w[i][j][k];
				xiz=xidz_w[i][j][k];
				etx=etdx_w[i][j][k];
				ety=etdy_w[i][j][k];
				etz=etdz_w[i][j][k];          
				ztx=ztdx_w[i][j][k];
				zty=ztdy_w[i][j][k];
				ztz=ztdz_w[i][j][k];
				ZTx=ztx/(sqrt(ztx*ztx+zty*zty+ztz*ztz));
				ZTy=zty/(sqrt(ztx*ztx+zty*zty+ztz*ztz));
				ZTz=ztz/(sqrt(ztx*ztx+zty*zty+ztz*ztz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				_U = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;
				_W = ztx*_u+zty*_v+ztz*_w;

				__W = ZTx*_u+ZTy*_v+ZTz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_T = _P/_rho;
				_C = K*_P/_rho;
				_H = 0.5*_VV+_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				W__ = ZTx*u_+ZTy*v_+ZTz*w_;       

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ = K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);

				/* flux varibale */
				rho = 0.5*(_rho+rho_);
				u = 0.5*(_u+u_);
				v = 0.5*(_v+v_);
				w = 0.5*(_w+w_);

				U = 0.5*(_U+U_);
				V = 0.5*(_V+V_);
				W = 0.5*(_W+W_); 

				_W_ = 0.5*(__W+W__); 

				VV = u*u+v*v+w*w;
				H = 0.5*(_H+H_);
				C = (H-0.5*VV)*(K-1);
				P = rho*C/K;
				T = P/rho;

				/* jump dU */
				dU1 = P_-_P;
				dU2 = u_-_u;
				dU3 = v_-_v;
				dU4 = w_-_w;
				dU5 = T_-_T;

				

				/* preconditioning */
				beta = max(VV/C,e);

				/*W=ztx*u+zty*v+ztz*w;*/
				//S=sqrt(W*W*(beta-1)+4*beta*C*C*(ztx*ztx+zty*zty+ztz*ztz));
				S=sqrt(W*W*(beta-1)*(beta-1)+4*beta*C*(ztx*ztx+zty*zty+ztz*ztz));
				/*_W_=ZTx*u+ZTy*v+ZTz*w;*/
				_S_=sqrt(_W_*_W_*(beta-1)*(beta-1)+4*beta*C*(ZTx*ZTx+ZTy*ZTy+ZTz*ZTz));

				temp = 4*K*_S_*T*beta;
				temp1 = S+W+W*beta;
				temp2 = -S+W+W*beta;
				temp3 = _S_-_W_+_W_*beta;
				temp4 = _S_+_W_-_W_*beta;
				temp5 = _S_*_S_-(beta-1)*(beta-1)*_W_*_W_;
				tempu1 = 2*T*ZTx*K*beta+u*_S_+u*(beta-1)*_W_;
				tempu2 = 2*T*ZTx*K*beta-u*_S_+u*(beta-1)*_W_;
				tempv1 = 2*T*ZTy*K*beta+v*_S_+v*(beta-1)*_W_;
				tempv2 = 2*T*ZTy*K*beta-v*_S_+v*(beta-1)*_W_;
				tempw1 = 2*T*ZTz*K*beta+w*_S_+w*(beta-1)*_W_;
				tempw2 = 2*T*ZTz*K*beta-w*_S_+w*(beta-1)*_W_;
				tempuvw = 2*T*K*beta*(u*ZTx+v*ZTy+w*ZTz);
				temph1 = tempuvw+H*_S_+H*_W_*beta-H*_W_; 
				temph2 = tempuvw-H*_S_+H*_W_*beta-H*_W_;


				d11 = 1/temp*(4*(K-1)*_S_*beta*fabs(W)+temp3*fabs(temp1)+temp4*fabs(temp2));
				d12 = -1/(2*temp)*(ZTx*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d13 = -1/(2*temp)*(ZTy*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d14 = -1/(2*temp)*(ZTz*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d15 = -rho*fabs(W)/T;

				d21 = 1/temp*(u*_S_*(4*(K-1)*beta*fabs(W)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ZTx*K*beta+u*_W_*(beta-1)));
				d22 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ZTy*ZTy+ZTz*ZTz)*fabs(W)+ZTx*(temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1))));
				d23 = 1/(2*temp)*(rho*ZTy*(-8*T*K*_S_*beta*ZTx*fabs(W)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d24 = 1/(2*temp)*(rho*ZTz*(-8*T*K*_S_*beta*ZTx*fabs(W)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d25 = -u*rho*fabs(W)/T;

				d31 = 1/temp*(v*_S_*(4*(K-1)*beta*fabs(W)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ZTy*K*beta+v*_W_*(beta-1)));
				d32 = 1/(2*temp)*(rho*ZTx*(-8*T*K*_S_*beta*ZTy*fabs(W)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d33 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ZTx*ZTx+ZTz*ZTz)*fabs(W)+ZTy*(temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1))));
				d34 = 1/(2*temp)*(rho*ZTz*(-8*T*K*_S_*beta*ZTy*fabs(W)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d35 = -v*rho*fabs(W)/T;

				d41 = 1/temp*(w*_S_*(4*(K-1)*beta*fabs(W)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ZTz*K*beta+w*_W_*(beta-1)));
				d42 = 1/(2*temp)*(rho*ZTx*(-8*T*K*_S_*beta*ZTz*fabs(W)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d43 = 1/(2*temp)*(rho*ZTy*(-8*T*K*_S_*beta*ZTz*fabs(W)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d44 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ZTx*ZTx+ZTy*ZTy)*fabs(W)+ZTz*(temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1))));
				d45 = -w*rho*fabs(W)/T;

				d51 = 1/temp*(_S_*(4*beta*(H*K-H-T*K)*fabs(W)+H*(fabs(temp2)+fabs(temp1)))-(fabs(temp2)-fabs(temp1))*(tempuvw+H*_W_*beta-H*_W_));
				d52 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-v*ZTx*ZTy-w*ZTx*ZTz+u*(ZTy*ZTy+ZTz*ZTz)*fabs(W))+ZTx*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d53 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*ZTx*ZTy-w*ZTy*ZTz+v*(ZTx*ZTx+ZTz*ZTz)*fabs(W))+ZTy*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d54 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*ZTx*ZTz-v*ZTy*ZTz+w*(ZTx*ZTx+ZTy*ZTy)*fabs(W))+ZTz*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d55 = (-H/T+K/(K-1))*rho*fabs(W);

				/* artificial viscosity */
				Fav1 = d11*dU1+d12*dU2+d13*dU3+d14*dU4+d15*dU5;
				Fav2 = d21*dU1+d22*dU2+d23*dU3+d24*dU4+d25*dU5;
				Fav3 = d31*dU1+d32*dU2+d33*dU3+d34*dU4+d35*dU5;
				Fav4 = d41*dU1+d42*dU2+d43*dU3+d44*dU4+d45*dU5;
				Fav5 = d51*dU1+d52*dU2+d53*dU3+d54*dU4+d55*dU5;


				/* inviscid fluxes */
				
				//EpZ[i][j][k] = 0.05;
				
				inFz1[i][j][k] = 0.5*((_rho*_W+rho_*W_-EpZ[i][j][k]*Fav1))/J_w[i][j][k];
				inFz2[i][j][k] = 0.5*((_rho*_u*_W+rho_*u_*W_+ztx*(_P+P_))-EpZ[i][j][k]*Fav2)/J_w[i][j][k];
				inFz3[i][j][k] = 0.5*((_rho*_v*_W+rho_*v_*W_+zty*(_P+P_))-EpZ[i][j][k]*Fav3)/J_w[i][j][k];
				inFz4[i][j][k] = 0.5*((_rho*_w*_W+rho_*w_*W_+ztz*(_P+P_))-EpZ[i][j][k]*Fav4)/J_w[i][j][k];
				inFz5[i][j][k] = 0.5*((_W*(3.5*_P+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_))-EpZ[i][j][k]*Fav5)/J_w[i][j][k];
				
				/*
				inFz1[i][j][k] = 0.5*((_rho*_W+rho_*W_))/J_w[i][j][k];
				inFz2[i][j][k] = 0.5*((_rho*_u*_W+rho_*u_*W_+ztx*(_P+P_)))/J_w[i][j][k];
				inFz3[i][j][k] = 0.5*((_rho*_v*_W+rho_*v_*W_+zty*(_P+P_)))/J_w[i][j][k];
				inFz4[i][j][k] = 0.5*((_rho*_w*_W+rho_*w_*W_+ztz*(_P+P_)))/J_w[i][j][k];
				inFz5[i][j][k] = 0.5*((_W*(3.5*_P+0.5*_rho*_VV)+W_*(3.5*P_+0.5*rho_*VV_)))/J_w[i][j][k];
				*/
				
			}
		}
	}

}