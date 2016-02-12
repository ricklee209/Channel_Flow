



#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "Resolution.h"
#include "Pre_selection.h"

extern int X_np;

#define min(a,b) (((a)<(b))?(a):(b)) 
#define max(a,b) (((a)>(b))?(a):(b)) 


void Flux_Y
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
double (*J_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_v)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*inFy1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*inFy5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
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
	double deltaU, deltaP, Cdiss, insqr;
	double beta,S,_S_;

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
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;						  ////	
//// ============================================ ////
		iend = gend[myid];						  ////
//// ============================================ ////
	
	

	for (i = istart; i <= iend; i++) {
#pragma omp parallel for private(k,temp)
		for ( j = 4; j <= ny-2; j++) {
			for (k = 2; k <= nz; k++) {
			
				temp = 1./(2/J[i][j-2][k]-13/J[i][j-1][k]+47/J[i][j][k]+27/J[i][j+1][k]-3/J[i][j+2][k]);

				ML1[i-1][j][k-1] = temp*(2*U1_[i][j-2][k]-13*U1_[i][j-1][k]+47*U1_[i][j][k]+27*U1_[i][j+1][k]-3*U1_[i][j+2][k]);
				ML2[i-1][j][k-1] = temp*(2*U2_[i][j-2][k]-13*U2_[i][j-1][k]+47*U2_[i][j][k]+27*U2_[i][j+1][k]-3*U2_[i][j+2][k]);
				ML3[i-1][j][k-1] = temp*(2*U3_[i][j-2][k]-13*U3_[i][j-1][k]+47*U3_[i][j][k]+27*U3_[i][j+1][k]-3*U3_[i][j+2][k]);
				ML4[i-1][j][k-1] = temp*(2*U4_[i][j-2][k]-13*U4_[i][j-1][k]+47*U4_[i][j][k]+27*U4_[i][j+1][k]-3*U4_[i][j+2][k]);
				ML5[i-1][j][k-1] = temp*(2*U5_[i][j-2][k]-13*U5_[i][j-1][k]+47*U5_[i][j][k]+27*U5_[i][j+1][k]-3*U5_[i][j+2][k]);

				temp = 1./(-3/J[i][j-2][k]+27/J[i][j-1][k]+47/J[i][j][k]-13/J[i][j+1][k]+2/J[i][j+2][k]);
						
				MR1[i-1][j-1][k-1] = temp*(-3*U1_[i][j-2][k]+27*U1_[i][j-1][k]+47*U1_[i][j][k]-13*U1_[i][j+1][k]+2*U1_[i][j+2][k]);
				MR2[i-1][j-1][k-1] = temp*(-3*U2_[i][j-2][k]+27*U2_[i][j-1][k]+47*U2_[i][j][k]-13*U2_[i][j+1][k]+2*U2_[i][j+2][k]);
				MR3[i-1][j-1][k-1] = temp*(-3*U3_[i][j-2][k]+27*U3_[i][j-1][k]+47*U3_[i][j][k]-13*U3_[i][j+1][k]+2*U3_[i][j+2][k]);
				MR4[i-1][j-1][k-1] = temp*(-3*U4_[i][j-2][k]+27*U4_[i][j-1][k]+47*U4_[i][j][k]-13*U4_[i][j+1][k]+2*U4_[i][j+2][k]);
				MR5[i-1][j-1][k-1] = temp*(-3*U5_[i][j-2][k]+27*U5_[i][j-1][k]+47*U5_[i][j][k]-13*U5_[i][j+1][k]+2*U5_[i][j+2][k]);
					
				
			}
		}
	}
	
	
	for (i = istart; i <= iend; i++) {
		for (k = 2; k <= nz; k++) {
		
			j = 1;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			

// ====================================================================================//
		
			
			j = 2;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
// ====================================================================================//		
			j = 3; 
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//	
		
			j = ny-1;
			
			temp = 1./(-1/J[i][j-1][k]+5/J[i][j][k]+2/J[i][j+1][k]);

			ML1[i-1][j][k-1] = temp*(-U1_[i][j-1][k]+5*U1_[i][j][k]+2*U1_[i][j+1][k]);
			ML2[i-1][j][k-1] = temp*(-U2_[i][j-1][k]+5*U2_[i][j][k]+2*U2_[i][j+1][k]);
			ML3[i-1][j][k-1] = temp*(-U3_[i][j-1][k]+5*U3_[i][j][k]+2*U3_[i][j+1][k]);
			ML4[i-1][j][k-1] = temp*(-U4_[i][j-1][k]+5*U4_[i][j][k]+2*U4_[i][j+1][k]);
			ML5[i-1][j][k-1] = temp*(-U5_[i][j-1][k]+5*U5_[i][j][k]+2*U5_[i][j+1][k]);
						
			temp = 1./(2/J[i][j-1][k]+5/J[i][j][k]-1/J[i][j+1][k]);
						
			MR1[i-1][j-1][k-1] = temp*(2*U1_[i][j-1][k]+5*U1_[i][j][k]-U1_[i][j+1][k]);
			MR2[i-1][j-1][k-1] = temp*(2*U2_[i][j-1][k]+5*U2_[i][j][k]-U2_[i][j+1][k]);
			MR3[i-1][j-1][k-1] = temp*(2*U3_[i][j-1][k]+5*U3_[i][j][k]-U3_[i][j+1][k]);
			MR4[i-1][j-1][k-1] = temp*(2*U4_[i][j-1][k]+5*U4_[i][j][k]-U4_[i][j+1][k]);
			MR5[i-1][j-1][k-1] = temp*(2*U5_[i][j-1][k]+5*U5_[i][j][k]-U5_[i][j+1][k]);
			
// ====================================================================================//



			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
					
			ML1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			ML2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			ML3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			ML4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			ML5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
			temp = 1./(1/J[i][j][k]+1/J[i][j-1][k]);
			
			MR1[i-1][j-1][k-1] = temp*(U1_[i][j][k]+U1_[i][j-1][k]);
			MR2[i-1][j-1][k-1] = temp*(U2_[i][j][k]+U2_[i][j-1][k]);
			MR3[i-1][j-1][k-1] = temp*(U3_[i][j][k]+U3_[i][j-1][k]);
			MR4[i-1][j-1][k-1] = temp*(U4_[i][j][k]+U4_[i][j-1][k]);
			MR5[i-1][j-1][k-1] = temp*(U5_[i][j][k]+U5_[i][j-1][k]);
			
			
// ====================================================================================//

			j = ny;
			
			temp = 1./(1/J[i][j][k]+1/J[i][j+1][k]);
				
			MR1[i-1][j][k-1] = temp*(U1_[i][j][k]+U1_[i][j+1][k]);
			MR2[i-1][j][k-1] = temp*(U2_[i][j][k]+U2_[i][j+1][k]);
			MR3[i-1][j][k-1] = temp*(U3_[i][j][k]+U3_[i][j+1][k]);
			MR4[i-1][j][k-1] = temp*(U4_[i][j][k]+U4_[i][j+1][k]);
			MR5[i-1][j][k-1] = temp*(U5_[i][j][k]+U5_[i][j+1][k]);
			
				
		}
	}
	
	
	



//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

	/*---Y fluxes---*/
	for (i = istart; i <= iend; i++) {

		
#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ETx,ETy,ETz,\
	_rho,_u,_v,_w,_U,_V,_W,__V,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,V__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_V_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,\
	temp,temp1,temp2,temp3,temp4,temp5,temp6,\
	deltaU, deltaP, Cdiss, insqr,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	k\
	)


		for (j = 1; j < nyy; j++) {
			for (k = 1; k < nz; k++) {

				etx=etdx_v[i][j][k];
				ety=etdy_v[i][j][k];
				etz=etdz_v[i][j][k];          

				insqr = 1.0/sqrt(etx*etx+ety*ety+etz*etz);

				ETx=etx*insqr;
				ETy=ety*insqr;
				ETz=etz*insqr;


				
				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				_V = etx*_u+ety*_v+etz*_w;

				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_C = K*_P/_rho;
				_H = 0.5*_VV+_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				V_ = etx*u_+ety*v_+etz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				C_ =K*P_/rho_;
				H_ = 0.5*VV_+C_/(K-1);


				dU1 = rho_-_rho;
				dU2 = rho_*u_-_rho*_u;
				dU3 = rho_*v_-_rho*_v;
				dU4 = rho_*w_-_rho*_w;
				dU5 = (P_/(K-1)+0.5*rho_*VV_)-(_P/(K-1)+0.5*_rho*_VV);



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
					temp2 = V*insqr/S*(V_-_V)*insqr;

					deltaU = (S-fabs(V*insqr))*temp1+temp2;

					deltaP = V*insqr/S*(P_-_P)+( S - fabs(V*insqr) )*rho*(V_-_V)*insqr;

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

					temp1 = 0.5*(_V+V_)+0.5*beta*(_V-V_);
					V_ = 0.5*(V_+_V)+0.5*beta*(V_-_V);
					_V = temp1;



					temp5 = sqrt(_rho);
					temp6 = sqrt(rho_);
					
					temp4 = temp5+temp6;

					rho = sqrt(_rho*rho_);
					u = (temp5*_u+temp6*u_)/temp4;
					v = (temp5*_v+temp6*v_)/temp4;;
					w = (temp5*_w+temp6*w_)/temp4;

					V = (temp5*_V+temp6*V_)/temp4;

					VV = u*u+v*v+w*w;
					H = (temp5*_H+temp6*H_)/temp4;
					C = (H-0.5*VV)*(K-1); /**** C = sqrt((H-0.5*VV)*(K-1)); ****/
					P = rho*C/K;


					S = sqrt(C);
					temp1 = (P_-_P)/rho/C;
					temp2 = V*insqr/S*(V_-_V)*insqr;

					deltaU = (S-fabs(V*insqr))*temp1+temp2;

					deltaP = V*insqr/S*(P_-_P)+( S - fabs(V*insqr) )*rho*(V_-_V)*insqr;
					
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
					temp2 = V*insqr/S*beta*(V_-_V)*insqr;

					deltaU = (S-fabs(V)*insqr)*temp1+temp2;
					deltaP = V*insqr/S*(P_-_P)+( beta*S - fabs(V*insqr) )*rho*(V_-_V)*insqr;


				#endif


					
				/* artificial viscosity */


				#ifdef ILES

					
					Fav1 = EpY[i][j][k]*fabs(V*insqr)*dU1+deltaU*rho;
					Fav2 = EpY[i][j][k]*fabs(V*insqr)*dU2+deltaU*rho*u+deltaP*ETx;
					Fav3 = EpY[i][j][k]*fabs(V*insqr)*dU3+deltaU*rho*v+deltaP*ETy;
					Fav4 = EpY[i][j][k]*fabs(V*insqr)*dU4+deltaU*rho*w+deltaP*ETz;
					Fav5 = EpY[i][j][k]*fabs(V*insqr)*dU5+deltaU*rho*H+deltaP*V*insqr;

				#else
					
					Fav1 = fabs(V*insqr)*dU1+deltaU*rho;
					Fav2 = fabs(V*insqr)*dU2+deltaU*rho*u+deltaP*ETx;
					Fav3 = fabs(V*insqr)*dU3+deltaU*rho*v+deltaP*ETy;
					Fav4 = fabs(V*insqr)*dU4+deltaU*rho*w+deltaP*ETz;
					Fav5 = fabs(V*insqr)*dU5+deltaU*rho*H+deltaP*V*insqr;

				#endif
					
				/* inviscid fluxes */

				
				inFy1[i][j][k] = 0.5*((_rho*_V+rho_*V_)-Fav1/insqr)/J_v[i][j][k];
				inFy2[i][j][k] = 0.5*((_rho*_u*_V+rho_*u_*V_+etx*(_P+P_))-Fav2/insqr)/J_v[i][j][k];
				inFy3[i][j][k] = 0.5*((_rho*_v*_V+rho_*v_*V_+ety*(_P+P_))-Fav3/insqr)/J_v[i][j][k];
				inFy4[i][j][k] = 0.5*((_rho*_w*_V+rho_*w_*V_+etz*(_P+P_))-Fav4/insqr)/J_v[i][j][k];
				inFy5[i][j][k] = 0.5*((_V*(3.5*_P+0.5*_rho*_VV)+V_*(3.5*P_+0.5*rho_*VV_))-Fav5/insqr)/J_v[i][j][k];
				
				if (j==1 | j==ny) {

					inFy1[i][j][k] = 0;
					inFy2[i][j][k] = 0;
					inFy3[i][j][k] = 0.5*(ety*(_P+P_)-Fav3/insqr)/J_v[i][j][k];
					inFy4[i][j][k] = 0;
					inFy5[i][j][k] = 0;

				}
				
				
			}
		}
	}
	

}