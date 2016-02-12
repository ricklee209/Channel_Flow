



#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "Resolution.h"
#include "Pre_selection.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Flux_X
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
	double temp,temp1,temp2,temp3,temp4,temp5,temp6;
	double deltaU, deltaP, Cdiss;
	double beta,S,_S_,insqr;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double Fav1,Fav2,Fav3,Fav4,Fav5;

/**** MUSCL 5th-order ****/
// ============================================================================ //


//// ============================================ ////
		if (myid ==0) istart = 4;		          ////	
		else istart = 2;	             		  ////	
//// ============================================ ////
		if (myid ==nproc-1) iend = gend[myid]-1;  ////
		else iend = gend[myid]+1;				  ////
//// ============================================ ////



	for (i = istart; i <= iend; i++) {
	
	#pragma omp parallel for private(k,temp)
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
				
				temp = 1./(2/J[i-2][j][k]-13/J[i-1][j][k]+47/J[i][j][k]+27/J[i+1][j][k]-3/J[i+2][j][k]);

				ML1[i][j-1][k-1] = temp*(2*U1_[i-2][j][k]-13*U1_[i-1][j][k]+47*U1_[i][j][k]+27*U1_[i+1][j][k]-3*U1_[i+2][j][k]);
				ML2[i][j-1][k-1] = temp*(2*U2_[i-2][j][k]-13*U2_[i-1][j][k]+47*U2_[i][j][k]+27*U2_[i+1][j][k]-3*U2_[i+2][j][k]);
				ML3[i][j-1][k-1] = temp*(2*U3_[i-2][j][k]-13*U3_[i-1][j][k]+47*U3_[i][j][k]+27*U3_[i+1][j][k]-3*U3_[i+2][j][k]);
				ML4[i][j-1][k-1] = temp*(2*U4_[i-2][j][k]-13*U4_[i-1][j][k]+47*U4_[i][j][k]+27*U4_[i+1][j][k]-3*U4_[i+2][j][k]);
				ML5[i][j-1][k-1] = temp*(2*U5_[i-2][j][k]-13*U5_[i-1][j][k]+47*U5_[i][j][k]+27*U5_[i+1][j][k]-3*U5_[i+2][j][k]);

				temp = 1./(-3/J[i-2][j][k]+27/J[i-1][j][k]+47/J[i][j][k]-13/J[i+1][j][k]+2/J[i+2][j][k]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i-2][j][k]+27*U1_[i-1][j][k]+47*U1_[i][j][k]-13*U1_[i+1][j][k]+2*U1_[i+2][j][k]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i-2][j][k]+27*U2_[i-1][j][k]+47*U2_[i][j][k]-13*U2_[i+1][j][k]+2*U2_[i+2][j][k]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i-2][j][k]+27*U3_[i-1][j][k]+47*U3_[i][j][k]-13*U3_[i+1][j][k]+2*U3_[i+2][j][k]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i-2][j][k]+27*U4_[i-1][j][k]+47*U4_[i][j][k]-13*U4_[i+1][j][k]+2*U4_[i+2][j][k]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i-2][j][k]+27*U5_[i-1][j][k]+47*U5_[i][j][k]-13*U5_[i+1][j][k]+2*U5_[i+2][j][k]);
					
				//if (myid == 1 && k == 24 & j == 24)
					//printf("%d\t%f\t%f\t%f\t%f\n",myid,ML1[i-1][j][k-1],MR1[i-1][j-1][k-1],ML2[i-1][j][k-1],MR2[i-1][j-1][k-1]);
				
			}
		}
	}

	if (myid == 0) {
		//#pragma omp parallel for private(k,k_)
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
		
				i = 2;
				
				temp = 1./(1/J[i][j][k]+1/J[i+1][j][k]);
			
				ML1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);
				
				i = 3; 
			
				temp = 1./(-1/J[i-1][j][k]+5/J[i][j][k]+2/J[i+1][j][k]);

				ML1[i][j-1][k-1] = temp*(-U1_[i-1][j][k]+5*U1_[i][j][k]+2*U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(-U2_[i-1][j][k]+5*U2_[i][j][k]+2*U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(-U3_[i-1][j][k]+5*U3_[i][j][k]+2*U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(-U4_[i-1][j][k]+5*U4_[i][j][k]+2*U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(-U5_[i-1][j][k]+5*U5_[i][j][k]+2*U5_[i+1][j][k]);
				
				temp = 1./(2/J[i-1][j][k]+5/J[i][j][k]-1/J[i+1][j][k]);
							
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i-1][j][k]+5*U1_[i][j][k]-U1_[i+1][j][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i-1][j][k]+5*U2_[i][j][k]-U2_[i+1][j][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i-1][j][k]+5*U3_[i][j][k]-U3_[i+1][j][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i-1][j][k]+5*U4_[i][j][k]-U4_[i+1][j][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i-1][j][k]+5*U5_[i][j][k]-U5_[i+1][j][k]);
					
				
				}
			}
		}


	if (myid == nproc-1) {
	
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {
			
				
				i = gend[myid];
							
				temp = 1./(-1/J[i-1][j][k]+5/J[i][j][k]+2/J[i+1][j][k]);

				ML1[i][j-1][k-1] = temp*(-U1_[i-1][j][k]+5*U1_[i][j][k]+2*U1_[i+1][j][k]);
				ML2[i][j-1][k-1] = temp*(-U2_[i-1][j][k]+5*U2_[i][j][k]+2*U2_[i+1][j][k]);
				ML3[i][j-1][k-1] = temp*(-U3_[i-1][j][k]+5*U3_[i][j][k]+2*U3_[i+1][j][k]);
				ML4[i][j-1][k-1] = temp*(-U4_[i-1][j][k]+5*U4_[i][j][k]+2*U4_[i+1][j][k]);
				ML5[i][j-1][k-1] = temp*(-U5_[i-1][j][k]+5*U5_[i][j][k]+2*U5_[i+1][j][k]);
							
				temp = 1./(2/J[i-1][j][k]+5/J[i][j][k]-1/J[i+1][j][k]);
							
				MR1[i-1][j-1][k-1] = temp*(2*U1_[i-1][j][k]+5*U1_[i][j][k]-U1_[i+1][j][k]);
				MR2[i-1][j-1][k-1] = temp*(2*U2_[i-1][j][k]+5*U2_[i][j][k]-U2_[i+1][j][k]);
				MR3[i-1][j-1][k-1] = temp*(2*U3_[i-1][j][k]+5*U3_[i][j][k]-U3_[i+1][j][k]);
				MR4[i-1][j-1][k-1] = temp*(2*U4_[i-1][j][k]+5*U4_[i][j][k]-U4_[i+1][j][k]);
				MR5[i-1][j-1][k-1] = temp*(2*U5_[i-1][j][k]+5*U5_[i][j][k]-U5_[i+1][j][k]);
				
				//if (k == 2 & j == 2)
				//printf("%d\t%f\t%f\t%f\t%f\n",i,ML1[i][j-1][k-1],MR1[i-1][j-1][k-1],ML2[i][j-1][k-1],MR2[i-1][j-1][k-1]);

					
				i = gend[myid];
				
				temp = 1./(1/J[i][j][k]+1/J[i+1][j][k]);
					
				MR1[i][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i+1][j][k]);
				MR2[i][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i+1][j][k]);
				MR3[i][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i+1][j][k]);
				MR4[i][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i+1][j][k]);
				MR5[i][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i+1][j][k]);
				
					
				//if (k == 2 & j == 2)
				//printf("%d\t%f\t%f\t%f\t%f\n",i,ML1[i][j-1][k-1],MR1[i-1][j-1][k-1],ML2[i][j-1][k-1],MR2[i-1][j-1][k-1]);
						

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

	/*---X fluxes---*/
	for (i = istart; i <= iend; i++) {

		
#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,XIx,XIy,XIz,\
	_rho,_u,_v,_w,_U,_V,_W,__U,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,U__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_U_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,\
	temp,temp1,temp2,temp3,temp4,temp5,temp6,\
	deltaU, deltaP, Cdiss,insqr, \
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)

		for (j = 1; j < ny; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_u[i][j][k];
				xiy=xidy_u[i][j][k];
				xiz=xidz_u[i][j][k];

				insqr = 1.0/sqrt(xix*xix+xiy*xiy+xiz*xiz);

				XIx=xix*insqr;
				XIy=xiy*insqr;
				XIz=xiz*insqr;


				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;
				

				_U = xix*_u+xiy*_v+xiz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);

				_C = K*_P/_rho; /**** _C = sqrt(K*_P/_rho); ****/
				_H = 0.5*_VV+_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;
				
				U_ = xix*u_+xiy*v_+xiz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);

				C_ = K*P_/rho_; /**** C_ = sqrt(K*P_/rho_); ****/
				H_ = 0.5*VV_+C_/(K-1);
				

				

				#if ROE != 2

					/* flux varibale */
					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					u = (temp5*_u+temp6*u_)/temp4;
					v = (temp5*_v+temp6*v_)/temp4;;
					w = (temp5*_w+temp6*w_)/temp4;

					U = (temp5*_U+temp6*U_)/temp4;
					V = (temp5*_V+temp6*V_)/temp4;
					W = (temp5*_W+temp6*W_)/temp4;

					VV = u*u+v*v+w*w;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1); /**** C = sqrt((H-0.5*VV)*(K-1)); ****/
					P = rho*C/K;

					/* jump dU */
				
					dU1 = rho_-_rho;
					dU2 = rho_*u_-_rho*_u;
					dU3 = rho_*v_-_rho*_v;
					dU4 = rho_*w_-_rho*_w;
					dU5 = (P_/(K-1)+0.5*rho_*VV_)-(_P/(K-1)+0.5*_rho*_VV);

				#endif


				#if ROE == 1 

					S = sqrt(C);
					temp1 = (P_-_P)/rho/C;
					temp2 = U*insqr/S*(U_-_U)*insqr;

					deltaU = (S-fabs(U*insqr))*temp1+temp2;
					deltaP = U*insqr/S*(P_-_P)+( S - fabs(U*insqr) )*rho*(U_-_U)*insqr;


				#elif ROE == 2

					beta = sqrt(max(VV_/C_,_VV/_C));

					
					temp1 = 0.5*(_u+u_)+0.5*beta*(_u-u_);
					temp2 = 0.5*(_v+v_)+0.5*beta*(_v-v_);
					temp3 = 0.5*(_w+w_)+0.5*beta*(_w-w_);

					u_ = 0.5*(u_+_u)+0.5*beta*(u_-_u);
					v_ = 0.5*(v_+_v)+0.5*beta*(v_-_v);
					w_ = 0.5*(w_+_w)+0.5*beta*(w_-_w);

					_u = temp1;
					_v = temp2;
					_w = temp3;

					temp1 = 0.5*(_U+U_)+0.5*beta*(_U-U_);
					U_ = 0.5*(U_+_U)+0.5*beta*(U_-_U);
					_U = temp1;



					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					u = (temp5*_u+temp6*u_)/temp4;
					v = (temp5*_v+temp6*v_)/temp4;;
					w = (temp5*_w+temp6*w_)/temp4;

					U = (temp5*_U+temp6*U_)/temp4;

					VV = u*u+v*v+w*w;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1); /**** C = sqrt((H-0.5*VV)*(K-1)); ****/
					P = rho*C/K;


					S = sqrt(C);
					temp1 = (P_-_P)/rho/C;
					temp2 = U*insqr/S*(U_-_U)*insqr;

					deltaU = (S-fabs(U*insqr))*temp1+temp2;
					deltaP = U*insqr/S*(P_-_P)+( S - fabs(U*insqr) )*rho*(U_-_U)*insqr;

					/* jump dU */
				
					dU1 = rho_-_rho;
					dU2 = rho_*u_-_rho*_u;
					dU3 = rho_*v_-_rho*_v;
					dU4 = rho_*w_-_rho*_w;
					dU5 = (P_/(K-1)+0.5*rho_*VV_)-(_P/(K-1)+0.5*_rho*_VV);



				#elif ROE == 3

					S = sqrt(C);
					beta = (fabs(u)+fabs(v)+fabs(w))/S;

					temp1 = (P_-_P)/rho/C;
					temp2 = U*insqr/S*beta*(U_-_U)*insqr;

					deltaU = (S-fabs(U)*insqr)*temp1+temp2;
					deltaP = U*insqr/S*(P_-_P)+( beta*S - fabs(U*insqr) )*rho*(U_-_U)*insqr;

					        
					
				#endif


				/* artificial viscosity */

				#ifdef ILES

					Fav1 = EpX[i][j][k]*fabs(U*insqr)*dU1+deltaU*rho;
					Fav2 = EpX[i][j][k]*fabs(U*insqr)*dU2+deltaU*rho*u+deltaP*XIx;
					Fav3 = EpX[i][j][k]*fabs(U*insqr)*dU3+deltaU*rho*v+deltaP*XIy;
					Fav4 = EpX[i][j][k]*fabs(U*insqr)*dU4+deltaU*rho*w+deltaP*XIz;
					Fav5 = EpX[i][j][k]*fabs(U*insqr)*dU5+deltaU*rho*H+deltaP*U;

				#else

					Fav1 = fabs(U*insqr)*dU1+deltaU*rho;
					Fav2 = fabs(U*insqr)*dU2+deltaU*rho*u+deltaP*XIx;
					Fav3 = fabs(U*insqr)*dU3+deltaU*rho*v+deltaP*XIy;
					Fav4 = fabs(U*insqr)*dU4+deltaU*rho*w+deltaP*XIz;
					Fav5 = fabs(U*insqr)*dU5+deltaU*rho*H+deltaP*U;

				#endif
					
				/* inviscid fluxes */
				
				inFx1[i][j][k] = 0.5*((_rho*_U+rho_*U_)-Fav1/insqr)/J_u[i][j][k];
				inFx2[i][j][k] = 0.5*((_rho*_u*_U+rho_*u_*U_+xix*(_P+P_))-Fav2/insqr)/J_u[i][j][k];
				inFx3[i][j][k] = 0.5*((_rho*_v*_U+rho_*v_*U_+xiy*(_P+P_))-Fav3/insqr)/J_u[i][j][k];
				inFx4[i][j][k] = 0.5*((_rho*_w*_U+rho_*w_*U_+xiz*(_P+P_))-Fav4/insqr)/J_u[i][j][k];
				inFx5[i][j][k] = 0.5*((_U*(3.5*_P+0.5*_rho*_VV)+U_*(3.5*P_+0.5*rho_*VV_))-Fav5/insqr)/J_u[i][j][k];
				
			}
		}
	}
	

}