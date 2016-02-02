



#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Flux_X_limiter
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
double (*J_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFx1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFx5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpX)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
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
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]+1;  ////
		else iend = gend[myid]+1;				  ////
//// ============================================ ////



	for (i = istart; i <= iend; i++) {
	
	#pragma omp parallel for private(k,\
		u1ii,u2ii,u3ii,u4ii,u5ii,u1i,u2i,u3i,u4i,u5i,\
		u1,u2,u3,u4,u5,\
		iu1,iu2,iu3,iu4,iu5,iiu1,iiu2,iiu3,iiu4,iiu5,\
		irL,rL,rLi,bL,irR,rR,rRi,bR,minL3,minR3,jL,jR)

		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				
				jL = J_u[i  ][j-1][k-1];
				jR = J_u[i-1][j-1][k-1];


				if ((myid==0 & i <= 3) | (myid==np-1 & i >= gend[np-1])) {

					iu1 = U1_[i-1][j][k];
					iu2 = U2_[i-1][j][k];
					iu3 = U3_[i-1][j][k];
					iu4 = U4_[i-1][j][k];
					iu5 = U5_[i-1][j][k];

					u1 = U1_[i][j][k];
					u2 = U2_[i][j][k];
					u3 = U3_[i][j][k];
					u4 = U4_[i][j][k];
					u5 = U5_[i][j][k];

					u1i = U1_[i+1][j][k];
					u2i = U2_[i+1][j][k];
					u3i = U3_[i+1][j][k];
					u4i = U4_[i+1][j][k];
					u5i = U5_[i+1][j][k];

					rL = (u1i - u1)/( u1  - iu1);
					rR = (u1 - iu1)/( u1i - u1 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML1[i  ][j-1][k-1] = (u1+0.5*max(0,minL3)*(u1  -iu1))*jL;
					MR1[i-1][j-1][k-1] = (u1-0.5*max(0,minR3)*(u1i - u1))*jR;
					
					
					rL = (u2i - u2)/( u2  - iu2);
					rR = (u2 - iu2)/( u2i - u2 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML2[i  ][j-1][k-1] = (u2+0.5*max(0,minL3)*(u2  -iu2))*jL;
					MR2[i-1][j-1][k-1] = (u2-0.5*max(0,minR3)*(u2i - u2))*jR;


					rL = (u3i - u3)/( u3  - iu3);
					rR = (u3 - iu3)/( u3i - u3 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML3[i  ][j-1][k-1] = (u3+0.5*max(0,minL3)*(u3  -iu3))*jL;
					MR3[i-1][j-1][k-1] = (u3-0.5*max(0,minR3)*(u3i - u3))*jR;


					rL = (u4i - u4)/( u4  - iu4);
					rR = (u4 - iu4)/( u4i - u4 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML4[i  ][j-1][k-1] = (u4+0.5*max(0,minL3)*(u4  -iu4))*jL;
					MR4[i-1][j-1][k-1] = (u4-0.5*max(0,minR3)*(u4i - u4))*jR;


					rL = (u5i - u5)/( u5  - iu5);
					rR = (u5 - iu5)/( u5i - u5 );
					
					bL = 1./3*(1+2*rL);
					bR = 1./3*(1+2*rR);

					minL3 = min(2,min(2*rL,bL));
					minR3 = min(2,min(2*rR,bR));

					ML5[i  ][j-1][k-1] = (u5+0.5*max(0,minL3)*(u5  -iu5))*jL;
					MR5[i-1][j-1][k-1] = (u5-0.5*max(0,minR3)*(u5i - u5))*jR;

				}
				else {


					iiu1 = U1_[i-2][j][k];
					iiu2 = U2_[i-2][j][k];
					iiu3 = U3_[i-2][j][k];
					iiu4 = U4_[i-2][j][k];
					iiu5 = U5_[i-2][j][k];

					iu1 = U1_[i-1][j][k];
					iu2 = U2_[i-1][j][k];
					iu3 = U3_[i-1][j][k];
					iu4 = U4_[i-1][j][k];
					iu5 = U5_[i-1][j][k];

					u1 = U1_[i][j][k];
					u2 = U2_[i][j][k];
					u3 = U3_[i][j][k];
					u4 = U4_[i][j][k];
					u5 = U5_[i][j][k];

					u1i = U1_[i+1][j][k];
					u2i = U2_[i+1][j][k];
					u3i = U3_[i+1][j][k];
					u4i = U4_[i+1][j][k];
					u5i = U5_[i+1][j][k];

					u1ii = U1_[i+2][j][k];
					u2ii = U2_[i+2][j][k];
					u3ii = U3_[i+2][j][k];
					u4ii = U4_[i+2][j][k];
					u5ii = U5_[i+2][j][k];

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

					ML1[i  ][j-1][k-1] = (u1+0.5*max(0,minL3)*(u1  -iu1))*jL;
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

					ML2[i  ][j-1][k-1] = (u2+0.5*max(0,minL3)*(u2  -iu2))*jL;
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

					ML3[i  ][j-1][k-1] = (u3+0.5*max(0,minL3)*(u3  -iu3))*jL;
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

					ML4[i  ][j-1][k-1] = (u4+0.5*max(0,minL3)*(u4  -iu4))*jL;
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

					ML5[i  ][j-1][k-1] = (u5+0.5*max(0,minL3)*(u5  -iu5))*jL;
					MR5[i-1][j-1][k-1] = (u5-0.5*max(0,minR3)*(u5i - u5))*jR;

				}


			}
		}
	}

// ============================================================================ //
/**** MUSCL 5th-order-end****/


		
//// ============================================ ////
		istart = 2;							      ////
//// ============================================ ////
		iend = gend[myid];			     		  ////
//// ============================================ ////

#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,XIx,XIy,XIz,\
	_rho,_u,_v,_w,_U,_V,_W,__U,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,U__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_U_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	j,k\
	)

	/*---X fluxes---*/
	for (i = istart; i <= iend; i++) {
		for (j = 1; j < ny; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_u[i][j][k];
				xiy=xidy_u[i][j][k];
				xiz=xidz_u[i][j][k];
				etx=etdx_u[i][j][k];
				ety=etdy_u[i][j][k];
				etz=etdz_u[i][j][k];          
				ztx=ztdx_u[i][j][k];
				zty=ztdy_u[i][j][k];
				ztz=ztdz_u[i][j][k];
				XIx=xix/(sqrt(xix*xix+xiy*xiy+xiz*xiz));
				XIy=xiy/(sqrt(xix*xix+xiy*xiy+xiz*xiz));
				XIz=xiz/(sqrt(xix*xix+xiy*xiy+xiz*xiz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;
				

				_U = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;
				_W = ztx*_u+zty*_v+ztz*_w;

				__U = XIx*_u+XIy*_v+XIz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);

				_T = _P/_rho;
				_C = K*_P/_rho; /**** _C = sqrt(K*_P/_rho); ****/
				_H = 0.5*_VV+_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;
				
				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				U__ = XIx*u_+XIy*v_+XIz*w_;  

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);

				T_ = P_/rho_;
				C_ = K*P_/rho_; /**** C_ = sqrt(K*P_/rho_); ****/
				H_ = 0.5*VV_+C_/(K-1);

				
				
				//if (myid == 1 && i == 3 & k == 24)
					//printf("%d\t%d\t%f\t%f\t%f\t%f\n",j,myid,rho_,_rho,u_,_u);
					
				
				/* flux varibale */
				rho = 0.5*(_rho+rho_);
				u = 0.5*(_u+u_);
				v = 0.5*(_v+v_);
				w = 0.5*(_w+w_);

				U = 0.5*(_U+U_);
				V = 0.5*(_V+V_);
				W = 0.5*(_W+W_); 

				_U_ = 0.5*(__U+U__);

				VV = u*u+v*v+w*w;
				H = 0.5*(_H+H_);
				C = (H-0.5*VV)*(K-1); /**** C = sqrt((H-0.5*VV)*(K-1)); ****/
				P = rho*C/K;
				T = P/rho;

				/* jump dU */
				dU1 = P_-_P;
				dU2 = u_-_u;
				dU3 = v_-_v;
				dU4 = w_-_w;
				dU5 = T_-_T;

				/* preconditioning */
				/*
				if (VV/C/C < e)
					beta = e;
				else
					beta = VV/C/C;


				S=sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*C*(xix*xix+xiy*xiy+xiz*xiz));
				_S_=sqrt(_U_*_U_*(beta-1)*(beta-1)+4*beta*C*C*(XIx*XIx+XIy*XIy+XIz*XIz));
				*/

				beta = max(VV/C,e);



				S=sqrt(U*U*(beta-1)*(beta-1)+4*beta*C*(xix*xix+xiy*xiy+xiz*xiz));
				_S_=sqrt(_U_*_U_*(beta-1)*(beta-1)+4*beta*C*(XIx*XIx+XIy*XIy+XIz*XIz));

				temp = 4*K*_S_*T*beta;
				temp1 = S+U+U*beta;
				temp2 = -S+U+U*beta;
				temp3 = _S_-_U_+_U_*beta;
				temp4 = _S_+_U_-_U_*beta;
				temp5 = _S_*_S_-(beta-1)*(beta-1)*_U_*_U_;
				tempu1 = 2*T*XIx*K*beta+u*_S_+u*(beta-1)*_U_;
				tempu2 = 2*T*XIx*K*beta-u*_S_+u*(beta-1)*_U_;
				tempv1 = 2*T*XIy*K*beta+v*_S_+v*(beta-1)*_U_;
				tempv2 = 2*T*XIy*K*beta-v*_S_+v*(beta-1)*_U_;
				tempw1 = 2*T*XIz*K*beta+w*_S_+w*(beta-1)*_U_;
				tempw2 = 2*T*XIz*K*beta-w*_S_+w*(beta-1)*_U_;
				tempuvw = 2*T*K*beta*(u*XIx+v*XIy+w*XIz);
				temph1 = tempuvw+H*_S_+H*_U_*beta-H*_U_; 
				temph2 = tempuvw-H*_S_+H*_U_*beta-H*_U_;

				d11 = 1/temp*(4*(K-1)*_S_*beta*fabs(U)+temp3*fabs(temp1)+temp4*fabs(temp2));
				d12 = -1/(2*temp)*(XIx*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d13 = -1/(2*temp)*(XIy*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d14 = -1/(2*temp)*(XIz*rho*(fabs(temp2)-fabs(temp1))*temp5);
				d15 = -rho*fabs(U)/T;

				d21 = 1/temp*(u*_S_*(4*(K-1)*beta*fabs(U)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*XIx*K*beta+u*_U_*(beta-1)));
				d22 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(XIy*XIy+XIz*XIz)*fabs(U)+XIx*(temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1))));
				d23 = 1/(2*temp)*(rho*XIy*(-8*T*K*_S_*beta*XIx*fabs(U)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d24 = 1/(2*temp)*(rho*XIz*(-8*T*K*_S_*beta*XIx*fabs(U)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d25 = -u*rho*fabs(U)/T;

				d31 = 1/temp*(v*_S_*(4*(K-1)*beta*fabs(U)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*XIy*K*beta+v*_U_*(beta-1)));
				d32 = 1/(2*temp)*(rho*XIx*(-8*T*K*_S_*beta*XIy*fabs(U)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d33 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(XIx*XIx+XIz*XIz)*fabs(U)+XIy*(temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1))));
				d34 = 1/(2*temp)*(rho*XIz*(-8*T*K*_S_*beta*XIy*fabs(U)+temp3*tempv2*fabs(temp2)+temp4*tempv1*fabs(temp1)));
				d35 = -v*rho*fabs(U)/T;

				d41 = 1/temp*(w*_S_*(4*(K-1)*beta*fabs(U)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*XIz*K*beta+w*_U_*(beta-1)));
				d42 = 1/(2*temp)*(rho*XIx*(-8*T*K*_S_*beta*XIz*fabs(U)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d43 = 1/(2*temp)*(rho*XIy*(-8*T*K*_S_*beta*XIz*fabs(U)+temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1)));
				d44 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(XIx*XIx+XIy*XIy)*fabs(U)+XIz*(temp3*tempw2*fabs(temp2)+temp4*tempw1*fabs(temp1))));
				d45 = -w*rho*fabs(U)/T;

				d51 = 1/temp*(_S_*(4*beta*(H*K-H-T*K)*fabs(U)+H*(fabs(temp2)+fabs(temp1)))-(fabs(temp2)-fabs(temp1))*(tempuvw+H*_U_*beta-H*_U_));
				d52 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-v*XIx*XIy-w*XIx*XIz+u*(XIy*XIy+XIz*XIz)*fabs(U))+XIx*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d53 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*XIx*XIy-w*XIy*XIz+v*(XIx*XIx+XIz*XIz)*fabs(U))+XIy*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d54 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(-u*XIx*XIz-v*XIy*XIz+w*(XIx*XIx+XIy*XIy)*fabs(U))+XIz*(temph2*temp3*fabs(temp2)+temph1*temp4*fabs(temp1))));
				d55 = (-H/T+K/(K-1))*rho*fabs(U);

				/* artificial viscosity */
				Fav1 = d11*dU1+d12*dU2+d13*dU3+d14*dU4+d15*dU5;
				Fav2 = d21*dU1+d22*dU2+d23*dU3+d24*dU4+d25*dU5;
				Fav3 = d31*dU1+d32*dU2+d33*dU3+d34*dU4+d35*dU5;
				Fav4 = d41*dU1+d42*dU2+d43*dU3+d44*dU4+d45*dU5;
				Fav5 = d51*dU1+d52*dU2+d53*dU3+d54*dU4+d55*dU5;

				/* inviscid fluxes */
				
				//EpX[i][j][k] = 0.05;
				
				inFx1[i][j][k] = 0.5*((_rho*_U+rho_*U_)-EpX[i][j][k]*Fav1)/J_u[i][j][k];
				inFx2[i][j][k] = 0.5*((_rho*_u*_U+rho_*u_*U_+xix*(_P+P_))-EpX[i][j][k]*Fav2)/J_u[i][j][k];
				inFx3[i][j][k] = 0.5*((_rho*_v*_U+rho_*v_*U_+xiy*(_P+P_))-EpX[i][j][k]*Fav3)/J_u[i][j][k];
				inFx4[i][j][k] = 0.5*((_rho*_w*_U+rho_*w_*U_+xiz*(_P+P_))-EpX[i][j][k]*Fav4)/J_u[i][j][k];
				inFx5[i][j][k] = 0.5*((_U*(3.5*_P+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_))-EpX[i][j][k]*Fav5)/J_u[i][j][k];
				
				/*
				inFx1[i][j][k] = 0.5*((_rho*_U+rho_*U_))/J_u[i][j][k];
				inFx2[i][j][k] = 0.5*((_rho*_u*_U+rho_*u_*U_+xix*(_P+P_)))/J_u[i][j][k];
				inFx3[i][j][k] = 0.5*((_rho*_v*_U+rho_*v_*U_+xiy*(_P+P_)))/J_u[i][j][k];
				inFx4[i][j][k] = 0.5*((_rho*_w*_U+rho_*w_*U_+xiz*(_P+P_)))/J_u[i][j][k];
				inFx5[i][j][k] = 0.5*((_U*(3.5*_P+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_)))/J_u[i][j][k];
				*/
			}
		}
	}
	

}