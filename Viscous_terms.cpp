




#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>

#include "Resolution.h"

extern int X_np;

void Viscous_terms
(
// ============================================================================ //
int myid,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*vF2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*vF4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*vF5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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


/**** LR = vFx2_2 ****/
double (*LR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** LR = vFx2_2-end ****/


/**** LL = vFx2_1 ****/
double (*LL1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LL2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LL3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LL4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*LL5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** LL = vFx2_1-end ****/


/**** MR = vFy2_2 ****/
double (*MR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** MR = vFy2_2 ****/


/**** ML = vFy2_1 ****/
double (*ML1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ML5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** ML = vFy2_1-end ****/


/**** NR = vFz2_2 ****/
double (*NR1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
/**** NR = vFz2_2-end ****/


/**** NL = vFz2_1 ****/
double (*NL1)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL2)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL3)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL4)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*NL5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
/**** NL = vFz2_1-end ****/
// ============================================================================ //
)

{

	#include "ijk.h"
	#include "Viscous_terms.h"

	#include "MPI_prm.h"
	#include "Mpi_division.h"

	double rho,U,V,W,VV,P,C,T,h,H;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;

	double 	irho,iU,iV,iW,iP,iT,
			jrho,jU,jV,jW,jP,jT,
			krho,kU,kV,kW,kP,kT,

			rhoi,Ui,Vi,Wi,Pi,Ti,
			rhoj,Uj,Vj,Wj,Pj,Tj,
			rhok,Uk,Vk,Wk,Pk,Tk,

			ijrho,ijU,ijV,ijW,ijP,ijT,
			jkrho,jkU,jkV,jkW,jkP,jkT,
			ikrho,ikU,ikV,ikW,ikP,ikT,

			rhoij,Uij,Vij,Wij,Pij,Tij,
			rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,
			rhoik,Uik,Vik,Wik,Pik,Tik,

			irhoj,iUj,iVj,iWj,iPj,iTj,
			irhok,iUk,iVk,iWk,iPk,iTk,

			jrhoi,jUi,jVi,jWi,jPi,jTi,
			jrhok,jUk,jVk,jWk,jPk,jTk,

			krhoj,kUj,kVj,kWj,kPj,kTj,
			krhoi,kUi,kVi,kWi,kPi,kTi,

			du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,
			du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,
			duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,
			dT_dx,dT_dy,dT_dz,

			mu_E,mu_T, Pr_E;

	double a1,a2,a3,a4,a5,a6,a7,
		   b1,b2,b3,b4,b5,b6,b7,
		   c1,c2,c3,c4,c5,c6,c7,
		   d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,
		   e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,
		   f1,f2,f3,f4,f5,f6,f7,f8,f9,f10;

	double Ux,Uy,Uz,
		   Vx,Vy,Vz,
		   Wx,Wy,Wz;

	

	double invXI = 1./(deltaXI);
	double invET = 1./(deltaET);
	double invZT = 1./(deltaZT);

	double inv4XI = 1./(4*deltaXI);
	double inv4ET = 1./(4*deltaET);
	double inv4ZT = 1./(4*deltaZT);
	
//// ============================================ ////
		istart = 3;		             			  ////	
//// ============================================ ////
		iend = gend[myid]+1;					  ////
//// ============================================ ////
		
		for (i = istart; i <= iend; i++) {


		
#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,\
	a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6,c7,\
	d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,\
	f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,\
	rho,U,V,W,VV,P,T,\
	irho,iU,iV,iW,iP,iT,\
	jrho,jU,jV,jW,jP,jT,\
	krho,kU,kV,kW,kP,kT,\
	rhoi,Ui,Vi,Wi,Pi,Ti,\
	rhoj,Uj,Vj,Wj,Pj,Tj,\
	rhok,Uk,Vk,Wk,Pk,Tk,\
	ijrho,ijU,ijV,ijW,ijP,ijT,\
	jkrho,jkU,jkV,jkW,jkP,jkT,\
	ikrho,ikU,ikV,ikW,ikP,ikT,\
	rhoij,Uij,Vij,Wij,Pij,Tij,\
	rhojk,Ujk,Vjk,Wjk,Pjk,Tjk,\
	rhoik,Uik,Vik,Wik,Pik,Tik,\
	irhoj,iUj,iVj,iWj,iPj,iTj,\
	irhok,iUk,iVk,iWk,iPk,iTk,\
	jrhoi,jUi,jVi,jWi,jPi,jTi,\
	jrhok,jUk,jVk,jWk,jPk,jTk,\
	krhoj,kUj,kVj,kWj,kPj,kTj,\
	krhoi,kUi,kVi,kWi,kPi,kTi,\
	Ux,Uy,Uz,\
	Vx,Vy,Vz,\
	Wx,Wy,Wz,\
	du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz,\
	du2_dx,du2_dy,du2_dz,dv2_dx,dv2_dy,dv2_dz,dw2_dx,dw2_dy,dw2_dz,\
	duv_dx,dvw_dx,duw_dx,duv_dy,dvw_dy,duw_dy,duv_dz,dvw_dz,duw_dz,\
	dT_dx,dT_dy,dT_dz,mu_E,Pr_E,\
	_k,__k,k\
	)
	
			for (j = 2; j <= nyy; j++) {
				for (__k = 0,_k = 1, k = 2; k <= nzz; __k++,_k++, k++) {


					// ==== 0 0 0 ==== //
				rho = U1_[i][j][k];
				U = U2_[i][j][k]/rho;
				V = U3_[i][j][k]/rho;
				W = U4_[i][j][k]/rho;     
				VV = U*U+V*V+W*W;
				P = (U5_[i][j][k]-0.5*rho*VV)*(K-1)*J[i][j][k];
				T = P/(rho*R)/J[i][j][k];


				// ==== -1 1 0 ==== //
				irhoj = U1_[i-1][j+1][k];
				iUj = U2_[i-1][j+1][k]/irhoj;
				iVj = U3_[i-1][j+1][k]/irhoj;
				iWj = U4_[i-1][j+1][k]/irhoj;
				iPj = (U5_[i-1][j+1][k]-0.5*irhoj*(iUj*iUj+iVj*iVj+iWj*iWj))*(K-1)*J[i-1][j+1][k];
				iTj = iPj/(irhoj*R)/J[i-1][j+1][k];

				// ==== -1 0 1 ==== //
				irhok = U1_[i-1][j][k+1];
				iUk = U2_[i-1][j][k+1]/irhok;
				iVk = U3_[i-1][j][k+1]/irhok;
				iWk = U4_[i-1][j][k+1]/irhok;
				iPk = (U5_[i-1][j][k+1]-0.5*irhok*(iUk*iUk+iVk*iVk+iWk*iWk))*(K-1)*J[i-1][j][k+1];
				iTk = iPk/(irhok*R)/J[i-1][j][k+1];

				// ==== 1 -1 0 ==== //
				jrhoi = U1_[i+1][j-1][k];
				jUi = U2_[i+1][j-1][k]/jrhoi;
				jVi = U3_[i+1][j-1][k]/jrhoi;
				jWi = U4_[i+1][j-1][k]/jrhoi;
				jPi = (U5_[i+1][j-1][k]-0.5*jrhoi*(jUi*jUi+jVi*jVi+jWi*jWi))*(K-1)*J[i+1][j-1][k];
				jTi = jPi/(jrhoi*R)/J[i+1][j-1][k];

				// ==== 1 0 -1 ==== //
				krhoi = U1_[i+1][j][k-1];
				kUi = U2_[i+1][j][k-1]/krhoi;
				kVi = U3_[i+1][j][k-1]/krhoi;
				kWi = U4_[i+1][j][k-1]/krhoi;
				kPi = (U5_[i+1][j][k-1]-0.5*krhoi*(kUi*kUi+kVi*kVi+kWi*kWi))*(K-1)*J[i+1][j][k-1];
				kTi = kPi/(krhoi*R)/J[i+1][j][k-1];


				// ==== 0 -1 1 ==== //
				jrhok = U1_[i][j-1][k+1];
				jUk = U2_[i][j-1][k+1]/jrhok;
				jVk = U3_[i][j-1][k+1]/jrhok;
				jWk = U4_[i][j-1][k+1]/jrhok;
				jPk = (U5_[i][j-1][k+1]-0.5*jrhok*(jUk*jUk+jVk*jVk+jWk*jWk))*(K-1)*J[i][j-1][k+1];
				jTk = jPk/(jrhok*R)/J[i][j-1][k+1];


				// ==== 0 1 -1 ==== //
				krhoj = U1_[i][j+1][k-1];
				kUj = U2_[i][j+1][k-1]/krhoj;
				kVj = U3_[i][j+1][k-1]/krhoj;
				kWj = U4_[i][j+1][k-1]/krhoj;
				kPj = (U5_[i][j+1][k-1]-0.5*krhoj*(kUj*kUj+kVj*kVj+kWj*kWj))*(K-1)*J[i][j+1][k-1];
				kTj = kPj/(krhoj*R)/J[i][j+1][k-1];

				// ==== -1 0 0 ==== //
				irho = U1_[i-1][j][k];
				iU = U2_[i-1][j][k]/irho;
				iV = U3_[i-1][j][k]/irho;
				iW = U4_[i-1][j][k]/irho;
				iP = (U5_[i-1][j][k]-0.5*irho*(iU*iU+iV*iV+iW*iW))*(K-1)*J[i-1][j][k];
				iT = iP/(irho*R)/J[i-1][j][k];

				// ==== 0 -1 0 ==== //
				jrho = U1_[i][j-1][k];
				jU = U2_[i][j-1][k]/jrho;
				jV = U3_[i][j-1][k]/jrho;
				jW = U4_[i][j-1][k]/jrho;
				jP = (U5_[i][j-1][k]-0.5*jrho*(jU*jU+jV*jV+jW*jW))*(K-1)*J[i][j-1][k];
				jT = jP/(jrho*R)/J[i][j-1][k];

				// ==== 0 0 -1 ==== //
				krho = U1_[i][j][k-1];
				kU = U2_[i][j][k-1]/krho;
				kV = U3_[i][j][k-1]/krho;
				kW = U4_[i][j][k-1]/krho;
				kP = (U5_[i][j][k-1]-0.5*krho*(kU*kU+kV*kV+kW*kW))*(K-1)*J[i][j][k-1];
				kT = kP/(krho*R)/J[i][j][k-1];

				// ==== -1 -1 0 ==== //
				ijrho = U1_[i-1][j-1][k];
				ijU = U2_[i-1][j-1][k]/ijrho;
				ijV = U3_[i-1][j-1][k]/ijrho;
				ijW = U4_[i-1][j-1][k]/ijrho;
				ijP = (U5_[i-1][j-1][k]-0.5*ijrho*(ijU*ijU+ijV*ijV+ijW*ijW))*(K-1)*J[i-1][j-1][k];
				ijT = ijP/(ijrho*R)/J[i-1][j-1][k];


				// ==== 0 -1 -1 ==== //
				jkrho = U1_[i][j-1][k-1];
				jkU = U2_[i][j-1][k-1]/jkrho;
				jkV = U3_[i][j-1][k-1]/jkrho;
				jkW = U4_[i][j-1][k-1]/jkrho;
				jkP = (U5_[i][j-1][k-1]-0.5*jkrho*(jkU*jkU+jkV*jkV+jkW*jkW))*(K-1)*J[i][j-1][k-1];
				jkT = jkP/(jkrho*R)/J[i][j-1][k-1];


				// ==== -1 0 -1 ==== //
				ikrho = U1_[i-1][j][k-1];
				ikU = U2_[i-1][j][k-1]/ikrho;
				ikV = U3_[i-1][j][k-1]/ikrho;
				ikW = U4_[i-1][j][k-1]/ikrho;
				ikP = (U5_[i-1][j][k-1]-0.5*ikrho*(ikU*ikU+ikV*ikV+ikW*ikW))*(K-1)*J[i-1][j][k-1];
				ikT = ikP/(ikrho*R)/J[i-1][j][k-1];
				



				// ==== 1 0 0 ==== //
				rhoi = U1_[i+1][j][k];
				Ui = U2_[i+1][j][k]/rhoi;
				Vi = U3_[i+1][j][k]/rhoi;
				Wi = U4_[i+1][j][k]/rhoi;
				Pi = (U5_[i+1][j][k]-0.5*rhoi*(Ui*Ui+Vi*Vi+Wi*Wi))*(K-1)*J[i+1][j][k];
				Ti = Pi/(rhoi*R)/J[i+1][j][k];

				// ==== 0 1 0 ==== //
				rhoj = U1_[i][j+1][k];
				Uj = U2_[i][j+1][k]/rhoj;
				Vj = U3_[i][j+1][k]/rhoj;
				Wj = U4_[i][j+1][k]/rhoj;
				Pj = (U5_[i][j+1][k]-0.5*rhoj*(Uj*Uj+Vj*Vj+Wj*Wj))*(K-1)*J[i][j+1][k];
				Tj = Pj/(rhoj*R)/J[i][j+1][k];

				// ==== 0 0 1  ==== //
				rhok = U1_[i][j][k+1];
				Uk = U2_[i][j][k+1]/rhok;
				Vk = U3_[i][j][k+1]/rhok;
				Wk = U4_[i][j][k+1]/rhok;
				Pk = (U5_[i][j][k+1]-0.5*rhok*(Uk*Uk+Vk*Vk+Wk*Wk))*(K-1)*J[i][j][k+1];
				Tk = Pk/(rhok*R)/J[i][j][k+1];
				

				/* derivatives of velocity */
				/* X-direction */
				du_dx = (U-iU)*invXI;
				dv_dx = (V-iV)*invXI;
				dw_dx = (W-iW)*invXI;


				du2_dx = (U*U-iU*iU)*invXI;
				dv2_dx = (V*V-iV*iV)*invXI;
				dw2_dx = (W*W-iW*iW)*invXI;


				duv_dx = (U*V-iU*iV)*invXI;
				dvw_dx = (V*W-iV*iW)*invXI;
				duw_dx = (U*W-iU*iW)*invXI;


				dT_dx = (T-iT)*invXI;


				/* Y-direction */
				du_dy = (Uj+iUj-jU-ijU)*inv4ET;
				dv_dy = (Vj+iVj-jV-ijV)*inv4ET;
				dw_dy = (Wj+iWj-jW-ijW)*inv4ET;

				du2_dy = (Uj*Uj+iUj*iUj-jU*jU-ijU*ijU)*inv4ET;
				dv2_dy = (Vj*Vj+iVj*iVj-jV*jV-ijV*ijV)*inv4ET;
				dw2_dy = (Wj*Wj+iWj*iWj-jW*jW-ijW*ijW)*inv4ET;


				duv_dy = (Uj*Vj+iUj*iVj-jU*jV-ijU*ijV)*inv4ET;
				dvw_dy = (Vj*Wj+iVj*iWj-jV*jW-ijV*ijW)*inv4ET;
				duw_dy = (Uj*Wj+iUj*iWj-jU*jW-ijU*ijW)*inv4ET;


				dT_dy = (Tj+iTj-jT-ijT)*inv4ET;


				/* Z-direction */
				du_dz = (Uk+iUk-kU-ikU)*inv4ZT;
				dv_dz = (Vk+iVk-kV-ikV)*inv4ZT;
				dw_dz = (Wk+iWk-kW-ikW)*inv4ZT;

				du2_dz = (Uk*Uk+iUk*iUk-kU*kU-ikU*ikU)*inv4ZT;
				dv2_dz = (Vk*Vk+iVk*iVk-kV*kV-ikV*ikV)*inv4ZT;
				dw2_dz = (Wk*Wk+iWk*iWk-kW*kW-ikW*ikW)*inv4ZT;


				duv_dz = (Uk*Vk+iUk*iVk-kU*kV-ikU*ikV)*inv4ZT;
				dvw_dz = (Vk*Wk+iVk*iWk-kV*kW-ikV*ikW)*inv4ZT;
				duw_dz = (Uk*Wk+iUk*iWk-kU*kW-ikU*ikW)*inv4ZT;


				dT_dz = (Tk+iTk-kT-ikT)*inv4ZT;


				/* viscous*/
				mu_E = mu_L;
				Pr_E = Pr_L;



				xix=xidx_u[i-1][j-1][k-1];
				xiy=xidy_u[i-1][j-1][k-1];
				xiz=xidz_u[i-1][j-1][k-1];
				etx=etdx_u[i-1][j-1][k-1];
				ety=etdy_u[i-1][j-1][k-1];
				etz=etdz_u[i-1][j-1][k-1];          
				ztx=ztdx_u[i-1][j-1][k-1];
				zty=ztdy_u[i-1][j-1][k-1];
				ztz=ztdz_u[i-1][j-1][k-1];


				a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
				a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
				a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
				a4 = xix*xix+xiy*xiy+xiz*xiz;
				a5 = (1./3)*xix*xiy;
				a6 = (1./3)*xiy*xiz;
				a7 = (1./3)*xix*xiz;

				b1 = (4./3)*etx*etx+ety*ety+etz*etz;
				b2 = etx*etx+(4./3)*ety*ety+etz*etz;
				b3 = etx*etx+ety*ety+(4./3)*etz*etz;
				b4 = etx*etx+ety*ety+etz*etz;
				b5 = (1./3)*etx*ety;
				b6 = (1./3)*ety*etz;
				b7 = (1./3)*etx*etz;

				c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
				c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
				c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
				c4 = ztx*ztx+zty*zty+ztz*ztz;
				c5 = (1./3)*ztx*zty;
				c6 = (1./3)*zty*ztz;
				c7 = (1./3)*ztx*ztz;

				d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
				d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
				d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
				d4 = xix*etx+xiy*ety+xiz*etz;
				d5 = xix*ety-(2./3)*xiy*etx;
				d6 = xix*etz-(2./3)*xiz*etx;
				d7 = xiy*etx-(2./3)*xix*ety;
				d8 = xiy*etz-(2./3)*xiz*ety;
				d9 = xiz*etx-(2./3)*xix*etz;
				d10 = xiz*ety-(2./3)*xiy*etz;

				e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
				e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
				e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
				e4 = xix*ztx+xiy*zty+xiz*ztz;
				e5 = xix*zty-(2./3)*xiy*ztx;
				e6 = xix*ztz-(2./3)*xiz*ztx;
				e7 = xiy*ztx-(2./3)*xix*zty;
				e8 = xiy*ztz-(2./3)*xiz*zty;
				e9 = xiz*ztx-(2./3)*xix*ztz;
				e10 = xiz*zty-(2./3)*xiy*ztz;

				f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
				f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
				f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
				f4 = etx*ztx+ety*zty+etz*ztz;
				f5 = etx*zty-(2./3)*ety*ztx;
				f6 = etx*ztz-(2./3)*etz*ztx;
				f7 = ety*ztx-(2./3)*etx*zty;
				f8 = ety*ztz-(2./3)*etz*zty;
				f9 = etz*ztx-(2./3)*etx*ztz;
				f10 = etz*zty-(2./3)*ety*ztz;

				Ux = 0.5*(U+iU);
				Vx = 0.5*(V+iV);
				Wx = 0.5*(W+iW);

				/* X viscid fluxes */
				LL2[i][j][k] = mu_E*(a1*du_dx+a5*dv_dx+a7*dw_dx+
					d1*du_dy+d7*dv_dy+d9*dw_dy+
					e1*du_dz+e7*dv_dz+e9*dw_dz)/J_u[i-1][j-1][k-1];

				LL3[i][j][k] = mu_E*(a5*du_dx+a2*dv_dx+a6*dw_dx+
					d5*du_dy+d2*dv_dy+d10*dw_dy+
					e5*du_dz+e2*dv_dz+e10*dw_dz)/J_u[i-1][j-1][k-1];

				LL4[i][j][k] = mu_E*(a7*du_dx+a6*dv_dx+a3*dw_dx+
					d6*du_dy+d8*dv_dy+d3*dw_dy+
					e6*du_dz+e8*dv_dz+e3*dw_dz)/J_u[i-1][j-1][k-1];

				LL5[i][j][k] = mu_E*(0.5*a1*du2_dx+0.5*a2*dv2_dx+0.5*a3*dw2_dx+a5*duv_dx+a6*dvw_dx+a7*duw_dx+Cv*K*a4*dT_dx/(Pr_E)+
					0.5*d1*du2_dy+0.5*d2*dv2_dy+0.5*d3*dw2_dy+d5*Vx*du_dy+d6*Wx*du_dy+d7*Ux*dv_dy+d8*Wx*dv_dy+d9*Ux*dw_dy+d10*Vx*dw_dy+Cv*K*d4*dT_dy/(Pr_E)+
					0.5*e1*du2_dz+0.5*e2*dv2_dz+0.5*e3*dw2_dz+e5*Vx*du_dz+e6*Wx*du_dz+e7*Ux*dv_dz+e8*Wx*dv_dz+e9*Ux*dw_dz+e10*Vx*dw_dz+Cv*K*e4*dT_dz/(Pr_E))/J_u[i-1][j-1][k-1];


				xix=xidx_v[i-1][j-1][k-1];
				xiy=xidy_v[i-1][j-1][k-1];
				xiz=xidz_v[i-1][j-1][k-1];
				etx=etdx_v[i-1][j-1][k-1];
				ety=etdy_v[i-1][j-1][k-1];
				etz=etdz_v[i-1][j-1][k-1];          
				ztx=ztdx_v[i-1][j-1][k-1];
				zty=ztdy_v[i-1][j-1][k-1];
				ztz=ztdz_v[i-1][j-1][k-1];


				/* derivatives of velocity */
				
				/* X-direction */
				du_dx = (Ui+jUi-iU-ijU)*inv4XI;
				dv_dx = (Vi+jVi-iV-ijV)*inv4XI;
				dw_dx = (Wi+jWi-iW-ijW)*inv4XI;

				du2_dx = (Ui*Ui+jUi*jUi-iU*iU-ijU*ijU)*inv4XI;
				dv2_dx = (Vi*Vi+jVi*jVi-iV*iV-ijV*ijV)*inv4XI;
				dw2_dx = (Wi*Wi+jWi*jWi-iW*iW-ijW*ijW)*inv4XI;


				duv_dx = (Ui*Vi+jUi*jVi-iU*iV-ijU*ijV)*inv4XI;
				dvw_dx = (Vi*Wi+jVi*jWi-iV*iW-ijV*ijW)*inv4XI;
				duw_dx = (Ui*Wi+jUi*jWi-iU*iW-ijU*ijW)*inv4XI;


				dT_dx = (Ti+jTi-iT-ijT)*inv4XI;


				/* Y-direction */
				du_dy = (U-jU)*invET;
				dv_dy = (V-jV)*invET;
				dw_dy = (W-jW)*invET;


				du2_dy = (U*U-jU*jU)*invET;
				dv2_dy = (V*V-jV*jV)*invET;
				dw2_dy = (W*W-jW*jW)*invET;


				duv_dy = (U*V-jU*jV)*invET;
				dvw_dy = (V*W-jV*jW)*invET;
				duw_dy = (U*W-jU*jW)*invET;


				dT_dy = (T-jT)*invET;


				/* Z-direction */
				du_dz = (Uk+jUk-kU-jkU)*inv4ZT;
				dv_dz = (Vk+jVk-kV-jkV)*inv4ZT;
				dw_dz = (Wk+jWk-kW-jkW)*inv4ZT;

				du2_dz = (Uk*Uk+jUk*jUk-kU*kU-jkU*jkU)*inv4ZT;
				dv2_dz = (Vk*Vk+jVk*jVk-kV*kV-jkV*jkV)*inv4ZT;
				dw2_dz = (Wk*Wk+jWk*jWk-kW*kW-jkW*jkW)*inv4ZT;


				duv_dz = (Uk*Vk+jUk*jVk-kU*kV-jkU*jkV)*inv4ZT;
				dvw_dz = (Vk*Wk+jVk*jWk-kV*kW-jkV*jkW)*inv4ZT;
				duw_dz = (Uk*Wk+jUk*jWk-kU*kW-jkU*jkW)*inv4ZT;


				dT_dz = (Tk+jTk-kT-jkT)*inv4ZT;


				a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
				a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
				a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
				a4 = xix*xix+xiy*xiy+xiz*xiz;
				a5 = (1./3)*xix*xiy;
				a6 = (1./3)*xiy*xiz;
				a7 = (1./3)*xix*xiz;

				b1 = (4./3)*etx*etx+ety*ety+etz*etz;
				b2 = etx*etx+(4./3)*ety*ety+etz*etz;
				b3 = etx*etx+ety*ety+(4./3)*etz*etz;
				b4 = etx*etx+ety*ety+etz*etz;
				b5 = (1./3)*etx*ety;
				b6 = (1./3)*ety*etz;
				b7 = (1./3)*etx*etz;

				c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
				c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
				c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
				c4 = ztx*ztx+zty*zty+ztz*ztz;
				c5 = (1./3)*ztx*zty;
				c6 = (1./3)*zty*ztz;
				c7 = (1./3)*ztx*ztz;

				d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
				d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
				d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
				d4 = xix*etx+xiy*ety+xiz*etz;
				d5 = xix*ety-(2./3)*xiy*etx;
				d6 = xix*etz-(2./3)*xiz*etx;
				d7 = xiy*etx-(2./3)*xix*ety;
				d8 = xiy*etz-(2./3)*xiz*ety;
				d9 = xiz*etx-(2./3)*xix*etz;
				d10 = xiz*ety-(2./3)*xiy*etz;

				e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
				e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
				e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
				e4 = xix*ztx+xiy*zty+xiz*ztz;
				e5 = xix*zty-(2./3)*xiy*ztx;
				e6 = xix*ztz-(2./3)*xiz*ztx;
				e7 = xiy*ztx-(2./3)*xix*zty;
				e8 = xiy*ztz-(2./3)*xiz*zty;
				e9 = xiz*ztx-(2./3)*xix*ztz;
				e10 = xiz*zty-(2./3)*xiy*ztz;

				f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
				f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
				f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
				f4 = etx*ztx+ety*zty+etz*ztz;
				f5 = etx*zty-(2./3)*ety*ztx;
				f6 = etx*ztz-(2./3)*etz*ztx;
				f7 = ety*ztx-(2./3)*etx*zty;
				f8 = ety*ztz-(2./3)*etz*zty;
				f9 = etz*ztx-(2./3)*etx*ztz;
				f10 = etz*zty-(2./3)*ety*ztz;

				Uy = 0.5*(U+jU);
				Vy = 0.5*(V+jV);
				Wy = 0.5*(W+jW);

				/* Y viscid fluxes */
				
				ML2[i][j][k] = mu_E*(d1*du_dx+d5*dv_dx+d6*dw_dx+
					b1*du_dy+b5*dv_dy+b7*dw_dy+
					f1*du_dz+f7*dv_dz+f9*dw_dz)/J_v[i-1][j-1][k-1];

				ML3[i][j][k] = mu_E*(d7*du_dx+d2*dv_dx+d8*dw_dx+
					b5*du_dy+b2*dv_dy+b6*dw_dy+
					f5*du_dz+f2*dv_dz+f10*dw_dz)/J_v[i-1][j-1][k-1];

				ML4[i][j][k] = mu_E*(d9*du_dx+d10*dv_dx+d3*dw_dx+
					b7*du_dy+b6*dv_dy+b3*dw_dy+
					f6*du_dz+f8*dv_dz+f3*dw_dz)/J_v[i-1][j-1][k-1];

				ML5[i][j][k] = mu_E*(0.5*b1*du2_dy+0.5*b2*dv2_dy+0.5*b3*dw2_dy+b5*duv_dy+b6*dvw_dy+b7*duw_dy+Cv*K*b4*dT_dy/(Pr_E)+
					0.5*d1*du2_dx+0.5*d2*dv2_dx+0.5*d3*dw2_dx+d5*Uy*dv_dx+d6*Uy*dw_dx+d7*Vy*du_dx+d8*Vy*dw_dx+d9*Wy*du_dx+d10*Wy*dv_dx+Cv*K*d4*dT_dx/(Pr_E)+
					0.5*f1*du2_dz+0.5*f2*dv2_dz+0.5*f3*dw2_dz+f5*Vy*du_dz+f6*Wy*du_dz+f7*Uy*dv_dz+f8*Wy*dv_dz+f9*Uy*dw_dz+f10*Vy*dw_dz+Cv*K*f4*dT_dz/(Pr_E))/J_v[i-1][j-1][k-1];


				/* derivatives of velocity */
				
				/* X-direction */
				du_dx = (Ui+kUi-iU-ikU)*inv4XI;
				dv_dx = (Vi+kVi-iV-ikV)*inv4XI;
				dw_dx = (Wi+kWi-iW-ikW)*inv4XI;

				du2_dx = (Ui*Ui+kUi*kUi-iU*iU-ikU*ikU)*inv4XI;
				dv2_dx = (Vi*Vi+kVi*kVi-iV*iV-ikV*ikV)*inv4XI;
				dw2_dx = (Wi*Wi+kWi*kWi-iW*iW-ikW*ikW)*inv4XI;


				duv_dx = (Ui*Vi+kUi*kVi-iU*iV-ikU*ikV)*inv4XI;
				dvw_dx = (Vi*Wi+kVi*kWi-iV*iW-ikV*ikW)*inv4XI;
				duw_dx = (Ui*Wi+kUi*kWi-iU*iW-ikU*ikW)*inv4XI;


				dT_dx = (Ti+kTi-iT-ikT)*inv4XI;


				/* Y-direction */
				du_dy = (Uj+kUj-jU-jkU)*inv4ET;
				dv_dy = (Vj+kVj-jV-jkV)*inv4ET;
				dw_dy = (Wj+kWj-jW-jkW)*inv4ET;

				du2_dy = (Uj*Uj+kUj*kUj-jU*jU-jkU*jkU)*inv4ET;
				dv2_dy = (Vj*Vj+kVj*kVj-jV*jV-jkV*jkV)*inv4ET;
				dw2_dy = (Wj*Wj+kWj*kWj-jW*jW-jkW*jkW)*inv4ET;


				duv_dy = (Uj*Vj+kUj*iVj-jU*jV-jkU*jkV)*inv4ET;
				dvw_dy = (Vj*Wj+kVj*iWj-jV*jW-jkV*jkW)*inv4ET;
				duw_dy = (Uj*Wj+kUj*iWj-jU*jW-jkU*jkW)*inv4ET;


				dT_dy = (Tj+kTj-jT-jkT)*inv4ET;


				/* Z-direction */
				du_dz = (U-kU)*invZT;
				dv_dz = (V-kV)*invZT;
				dw_dz = (W-kW)*invZT;

				du2_dz = (U*U-kU*kU)*invZT;
				dv2_dz = (V*V-kV*kV)*invZT;
				dw2_dz = (W*W-kW*kW)*invZT;

				duv_dz = (U*V-kU*kV)*invZT;
				dvw_dz = (V*W-kV*kW)*invZT;
				duw_dz = (U*W-kU*kW)*invZT;

				dT_dz = (T-kT)*invZT;


				xix=xidx_w[i-1][j-1][k-1];
				xiy=xidy_w[i-1][j-1][k-1];
				xiz=xidz_w[i-1][j-1][k-1];
				etx=etdx_w[i-1][j-1][k-1];
				ety=etdy_w[i-1][j-1][k-1];
				etz=etdz_w[i-1][j-1][k-1];          
				ztx=ztdx_w[i-1][j-1][k-1];
				zty=ztdy_w[i-1][j-1][k-1];
				ztz=ztdz_w[i-1][j-1][k-1];


				a1 = (4./3)*xix*xix+xiy*xiy+xiz*xiz;
				a2 = xix*xix+(4./3)*xiy*xiy+xiz*xiz;
				a3 = xix*xix+xiy*xiy+(4./3)*xiz*xiz;
				a4 = xix*xix+xiy*xiy+xiz*xiz;
				a5 = (1./3)*xix*xiy;
				a6 = (1./3)*xiy*xiz;
				a7 = (1./3)*xix*xiz;

				b1 = (4./3)*etx*etx+ety*ety+etz*etz;
				b2 = etx*etx+(4./3)*ety*ety+etz*etz;
				b3 = etx*etx+ety*ety+(4./3)*etz*etz;
				b4 = etx*etx+ety*ety+etz*etz;
				b5 = (1./3)*etx*ety;
				b6 = (1./3)*ety*etz;
				b7 = (1./3)*etx*etz;

				c1 = (4./3)*ztx*ztx+zty*zty+ztz*ztz;
				c2 = ztx*ztx+(4./3)*zty*zty+ztz*ztz;
				c3 = ztx*ztx+zty*zty+(4./3)*ztz*ztz;
				c4 = ztx*ztx+zty*zty+ztz*ztz;
				c5 = (1./3)*ztx*zty;
				c6 = (1./3)*zty*ztz;
				c7 = (1./3)*ztx*ztz;

				d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
				d2 = xix*etx+(4./3)*xiy*ety+xiz*etz;
				d3 = xix*etx+xiy*ety+(4./3)*xiz*etz;
				d4 = xix*etx+xiy*ety+xiz*etz;
				d5 = xix*ety-(2./3)*xiy*etx;
				d6 = xix*etz-(2./3)*xiz*etx;
				d7 = xiy*etx-(2./3)*xix*ety;
				d8 = xiy*etz-(2./3)*xiz*ety;
				d9 = xiz*etx-(2./3)*xix*etz;
				d10 = xiz*ety-(2./3)*xiy*etz;

				e1 = (4./3)*xix*ztx+xiy*zty+xiz*ztz;
				e2 = xix*ztx+(4./3)*xiy*zty+xiz*ztz;
				e3 = xix*ztx+xiy*zty+(4./3)*xiz*ztz;
				e4 = xix*ztx+xiy*zty+xiz*ztz;
				e5 = xix*zty-(2./3)*xiy*ztx;
				e6 = xix*ztz-(2./3)*xiz*ztx;
				e7 = xiy*ztx-(2./3)*xix*zty;
				e8 = xiy*ztz-(2./3)*xiz*zty;
				e9 = xiz*ztx-(2./3)*xix*ztz;
				e10 = xiz*zty-(2./3)*xiy*ztz;

				f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
				f2 = etx*ztx+(4./3)*ety*zty+etz*ztz;
				f3 = etx*ztx+ety*zty+(4./3)*etz*ztz;
				f4 = etx*ztx+ety*zty+etz*ztz;
				f5 = etx*zty-(2./3)*ety*ztx;
				f6 = etx*ztz-(2./3)*etz*ztx;
				f7 = ety*ztx-(2./3)*etx*zty;
				f8 = ety*ztz-(2./3)*etz*zty;
				f9 = etz*ztx-(2./3)*etx*ztz;
				f10 = etz*zty-(2./3)*ety*ztz;


				Uz = 0.5*(U+kU);
				Vz = 0.5*(V+kV);
				Wz = 0.5*(W+kW);


				/* Z viscid fluxes */
				NL2[i][j][k] = mu_E*(e1*du_dx+e5*dv_dx+e6*dw_dx+
					f1*du_dy+f5*dv_dy+f6*dw_dy+
					c1*du_dz+c5*dv_dz+c7*dw_dz)/J_w[i-1][j-1][k-1];

				NL3[i][j][k] = mu_E*(e7*du_dx+e2*dv_dx+e8*dw_dx+
					f7*du_dy+f2*dv_dy+f8*dw_dy+
					c5*du_dz+c2*dv_dz+c6*dw_dz)/J_w[i-1][j-1][k-1];

				NL4[i][j][k] = mu_E*(e9*du_dx+e10*dv_dx+e3*dw_dx+
					f9*du_dy+f10*dv_dy+f3*dw_dy+
					c7*du_dz+c6*dv_dz+c3*dw_dz)/J_w[i-1][j-1][k-1];

				NL5[i][j][k] = mu_E*(0.5*c1*du2_dz+0.5*c2*dv2_dz+0.5*c3*dw2_dz+c5*duv_dz+c6*dvw_dz+c7*duw_dz+Cv*K*c4*dT_dz/(Pr_E)+
					0.5*e1*du2_dx+0.5*e2*dv2_dx+0.5*e3*dw2_dx+e5*Uz*dv_dx+e6*Uz*dw_dx+e7*Vz*du_dx+e8*Vz*dw_dx+e9*Wz*du_dx+e10*Wz*dv_dx+Cv*K*e4*dT_dx/(Pr_E)+
					0.5*f1*du2_dy+0.5*f2*dv2_dy+0.5*f3*dw2_dy+f5*Uz*dv_dy+f6*Uz*dw_dy+f7*Vz*du_dy+f8*Vz*dw_dy+f9*Wz*du_dy+f10*Wz*dv_dy+Cv*K*f4*dT_dy/(Pr_E))/J_w[i-1][j-1][k-1];
					
				
					}
				}
			}


// =================== //
	istart = 3;        //
	iend = gend[myid]; //
// =================== //
	
		for (i = istart; i <= iend; i++) {


#pragma omp parallel for private(\
	k,k_\
	)

			for (j = 2; j <= ny; j++) {
				for (k = 2, k_ = 3; k <= nz; k++, k_++) {

					
					vF2[i][j][k] = (LL2[i+1][j][k]-LL2[i][j][k])*invXI+
						(ML2[i][j+1][k]-ML2[i][j][k])*invET+
						(NL2[i][j][k_]-NL2[i][j][k])*invZT;
						
					vF3[i][j][k] = (LL3[i+1][j][k]-LL3[i][j][k])*invXI+
						(ML3[i][j+1][k]-ML3[i][j][k])*invET+
						(NL3[i][j][k_]-NL3[i][j][k])*invZT;

					vF4[i][j][k] = (LL4[i+1][j][k]-LL4[i][j][k])*invXI+
						(ML4[i][j+1][k]-ML4[i][j][k])*invET+
						(NL4[i][j][k_]-NL4[i][j][k])*invZT;

					vF5[i][j][k] = (LL5[i+1][j][k]-LL5[i][j][k])*invXI+
						(ML5[i][j+1][k]-ML5[i][j][k])*invET+
						(NL5[i][j][k_]-NL5[i][j][k])*invZT;
					
				}
			}
		}
		


}