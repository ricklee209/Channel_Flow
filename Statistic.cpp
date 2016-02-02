




#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <omp.h>
#include <math.h>


#include "Resolution.h"

extern int X_np;

void Statistic
(
// ============================================================================ //
int myid,
int step,
int iteration,
int statistic_step,

double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

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

double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
)

{

	
#include "ijk.h"
#include "Viscous_terms.h"

#include "MPI_prm.h"
#include "Mpi_division.h"

	char LESdata[100];
	
	double rho,U,V,W,VV,P,C,T,h,H;
	double u,v,w;
	double temp,temp1,temp2,temp3,temp4,temp5;
	double beta,S,_S_;

	double m11,m12,m13,m14,m15,m22,m23,m24,m25,m33,m34,m35,m44,m45,m55;
	double thedac1,thedac2,thedac3,thedad1,thedad2,thedad3;

	double xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
	double XIx,XIy,XIz,ETx,ETy,ETz,ZTx,ZTy,ZTz;
	double _rho,_u,_v,_w,_U,_V,_W,__U,__V,__W,_VV,_P,_T,_C,_H;
	double rho_,u_,v_,w_,U_,V_,W_,U__,V__,W__,VV_,P_,T_,C_,H_;
	double _U_,_V_,_W_;
	double dU1,dU2,dU3,dU4,dU5;
	double tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2;
	double d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35;
	double d41,d42,d43,d44,d45,d51,d52,d53,d54,d55, Fav1,Fav2,Fav3,Fav4,Fav5;


// ============================================================================================================= //
	
	static double Um[Y_m],Vm[Y_m],Wm[Y_m],UUm[Y_m],VVm[Y_m],WWm[Y_m],UVm[Y_m],
		ULm[Y_m],URm[Y_m],
		VLm[Y_m],VRm[Y_m],
		rhoRm[Y_m],rhoLm[Y_m],
		RUVm[Y_m],SXYm[Y_m],Roe_Dm[Y_m],dRLm[Y_m],
		SXY_up,SXY_bottom;

	// ============================ //

	static double Um0[Y_m],Vm0[Y_m],Wm0[Y_m],UUm0[Y_m],VVm0[Y_m],WWm0[Y_m],UVm0[Y_m],
		ULm0[Y_m],URm0[Y_m],
		VLm0[Y_m],VRm0[Y_m],
		rhoRm0[Y_m],rhoLm0[Y_m],
		RUVm0[Y_m],SXYm0[Y_m],Roe_Dm0[Y_m],dRLm0[Y_m];

	// ============================ //
	
	double (*MUL)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];
	double (*MUR)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];

	double (*MVL)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];
	double (*MVR)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];

	double (*MrhoL)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];
	double (*MrhoR)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];

	double (*MRUV)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];
	double (*MSXY)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];

	double (*MRoe_D)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];
	double (*MdRL)[Y_m][Z_m] = new double[X_m][Y_m][Z_m];

// ============================================================================================================= //


/**** MUSCL 5th-order ****/
// ============================================================================ //
	
//// ============================================ ////
		if (myid ==0) istart = 2;		          ////	
		else istart = 3;						  ////	
//// ============================================ ////
		iend = gend[myid];						  ////
//// ============================================ ////

	#pragma omp parallel for private(\
	j,k,j_,k_,j__,k__,_j,_k,__j,__k,\
	m11,m12,m13,m14,m15,m22,m23,m24,m25,m33,m34,m35,m44,m45,m55,\
	thedac1,thedac2,thedac3,thedad1,thedad2,thedad3)

	for (i = istart; i <= iend; i++) {
		for ( j = 2; j <= ny; j++) {
			for (k = 2; k <= nz; k++) {

				
				m11 = 1./J[i][j-2][k];
				m12 = m11+1./J[i][j-1][k];
				m13 = m12+1./J[i][j][k];
				m14 = m13+1./J[i][j+1][k];
				m15 = m14+1./J[i][j+2][k];

				m22 = 1./J[i][j-1][k];
				m23 = m22+1./J[i][j][k];
				m24 = m23+1./J[i][j+1][k];
				m25 = m24+1./J[i][j+2][k];

				m33 = 1./J[i][j][k];
				m34 = m33+1./J[i][j+1][k];
				m35 = m34+1./J[i][j+2][k];

				m44 = 1./J[i][j+1][k];
				m45 = m44+1./J[i][j+2][k];

				m55 = 1./J[i][j+2][k];

				thedac1 = m13/m15*m23/m25;
				thedac2 = m13/m15*m45/m25*(m25/m14+1);
				thedac3 = m44/m15*m45/m14;

				thedad1 = m22/m15*m12/m25;
				thedad2 = m12/m15*m35/m25*(m25/m14+1);
				thedad3 = m34/m15*m35/m14;

				if (j==2 | j==ny) {

					ML1[i-1][j][k-1] = U1_[i][j][k]*J[i][j][k]+m33/m24*m23/m34*(U1_[i][j+1][k]*J[i][j+1][k]-U1_[i][j][k]*J[i][j][k])-m33/m24*m44/m23*(U1_[i][j-1][k]*J[i][j-1][k]-U1_[i][j][k]*J[i][j][k]);
					ML2[i-1][j][k-1] = U2_[i][j][k]*J[i][j][k]+m33/m24*m23/m34*(U2_[i][j+1][k]*J[i][j+1][k]-U2_[i][j][k]*J[i][j][k])-m33/m24*m44/m23*(U2_[i][j-1][k]*J[i][j-1][k]-U2_[i][j][k]*J[i][j][k]);
					ML3[i-1][j][k-1] = U3_[i][j][k]*J[i][j][k]+m33/m24*m23/m34*(U3_[i][j+1][k]*J[i][j+1][k]-U3_[i][j][k]*J[i][j][k])-m33/m24*m44/m23*(U3_[i][j-1][k]*J[i][j-1][k]-U3_[i][j][k]*J[i][j][k]);
					ML4[i-1][j][k-1] = U4_[i][j][k]*J[i][j][k]+m33/m24*m23/m34*(U4_[i][j+1][k]*J[i][j+1][k]-U4_[i][j][k]*J[i][j][k])-m33/m24*m44/m23*(U4_[i][j-1][k]*J[i][j-1][k]-U4_[i][j][k]*J[i][j][k]);
					ML5[i-1][j][k-1] = U5_[i][j][k]*J[i][j][k]+m33/m24*m23/m34*(U5_[i][j+1][k]*J[i][j+1][k]-U5_[i][j][k]*J[i][j][k])-m33/m24*m44/m23*(U5_[i][j-1][k]*J[i][j-1][k]-U5_[i][j][k]*J[i][j][k]);

					MR1[i-1][j-1][k-1] = U1_[i][j][k]*J[i][j][k]+m33/m24*m34/m23*(U1_[i][j-1][k]*J[i][j-1][k]-U1_[i][j][k]*J[i][j][k])-m33/m24*m22/m34*(U1_[i][j+1][k]*J[i][j+1][k]-U1_[i][j][k]*J[i][j][k]);
					MR2[i-1][j-1][k-1] = U2_[i][j][k]*J[i][j][k]+m33/m24*m34/m23*(U2_[i][j-1][k]*J[i][j-1][k]-U2_[i][j][k]*J[i][j][k])-m33/m24*m22/m34*(U2_[i][j+1][k]*J[i][j+1][k]-U2_[i][j][k]*J[i][j][k]);
					MR3[i-1][j-1][k-1] = U3_[i][j][k]*J[i][j][k]+m33/m24*m34/m23*(U3_[i][j-1][k]*J[i][j-1][k]-U3_[i][j][k]*J[i][j][k])-m33/m24*m22/m34*(U3_[i][j+1][k]*J[i][j+1][k]-U3_[i][j][k]*J[i][j][k]);
					MR4[i-1][j-1][k-1] = U4_[i][j][k]*J[i][j][k]+m33/m24*m34/m23*(U4_[i][j-1][k]*J[i][j-1][k]-U4_[i][j][k]*J[i][j][k])-m33/m24*m22/m34*(U4_[i][j+1][k]*J[i][j+1][k]-U4_[i][j][k]*J[i][j][k]);
					MR5[i-1][j-1][k-1] = U5_[i][j][k]*J[i][j][k]+m33/m24*m34/m23*(U5_[i][j-1][k]*J[i][j-1][k]-U5_[i][j][k]*J[i][j][k])-m33/m24*m22/m34*(U5_[i][j+1][k]*J[i][j+1][k]-U5_[i][j][k]*J[i][j][k]);

				}



				else  {

					ML1[i-1][j][k-1] = thedac1*(U1_[i][j+1][k]*J[i][j+1][k]+m44/m35*m45/m34*(U1_[i][j][k]*J[i][j][k]    -U1_[i][j+1][k]*J[i][j+1][k])    -m44/m35*m33/m45*(U1_[i][j+2][k]*J[i][j+2][k]-U1_[i][j+1][k]*J[i][j+1][k]))+
						thedac2*(U1_[i][j][k]  *J[i][j][k]  +m33/m24*m23/m34*(U1_[i][j+1][k]*J[i][j+1][k]-U1_[i][j][k]*J[i][j][k])        -m33/m24*m44/m23*(U1_[i][j-1][k]*J[i][j-1][k]-U1_[i][j][k]*J[i][j][k]))+
						thedac3*(U1_[i][j-1][k]*J[i][j-1][k]+m33/m12*m23/m13*(U1_[i][j-2][k]*J[i][j-2][k]-U1_[i][j-1][k]*J[i][j-1][k])+(1+m33/m23+m33/m13)*(U1_[i][j][k]*J[i][j][k]    -U1_[i][j-1][k]*J[i][j-1][k]));

					ML2[i-1][j][k-1] = thedac1*(U2_[i][j+1][k]*J[i][j+1][k]+m44/m35*m45/m34*(U2_[i][j][k]*J[i][j][k]    -U2_[i][j+1][k]*J[i][j+1][k])    -m44/m35*m33/m45*(U2_[i][j+2][k]*J[i][j+2][k]-U2_[i][j+1][k]*J[i][j+1][k]))+
						thedac2*(U2_[i][j][k]  *J[i][j][k]  +m33/m24*m23/m34*(U2_[i][j+1][k]*J[i][j+1][k]-U2_[i][j][k]*J[i][j][k])        -m33/m24*m44/m23*(U2_[i][j-1][k]*J[i][j-1][k]-U2_[i][j][k]*J[i][j][k]))+
						thedac3*(U2_[i][j-1][k]*J[i][j-1][k]+m33/m12*m23/m13*(U2_[i][j-2][k]*J[i][j-2][k]-U2_[i][j-1][k]*J[i][j-1][k])+(1+m33/m23+m33/m13)*(U2_[i][j][k]*J[i][j][k]    -U2_[i][j-1][k]*J[i][j-1][k]));

					ML3[i-1][j][k-1] = thedac1*(U3_[i][j+1][k]*J[i][j+1][k]+m44/m35*m45/m34*(U3_[i][j][k]*J[i][j][k]    -U3_[i][j+1][k]*J[i][j+1][k])    -m44/m35*m33/m45*(U3_[i][j+2][k]*J[i][j+2][k]-U3_[i][j+1][k]*J[i][j+1][k]))+
						thedac2*(U3_[i][j][k]  *J[i][j][k]  +m33/m24*m23/m34*(U3_[i][j+1][k]*J[i][j+1][k]-U3_[i][j][k]*J[i][j][k])        -m33/m24*m44/m23*(U3_[i][j-1][k]*J[i][j-1][k]-U3_[i][j][k]*J[i][j][k]))+
						thedac3*(U3_[i][j-1][k]*J[i][j-1][k]+m33/m12*m23/m13*(U3_[i][j-2][k]*J[i][j-2][k]-U3_[i][j-1][k]*J[i][j-1][k])+(1+m33/m23+m33/m13)*(U3_[i][j][k]*J[i][j][k]    -U3_[i][j-1][k]*J[i][j-1][k]));

					ML4[i-1][j][k-1] = thedac1*(U4_[i][j+1][k]*J[i][j+1][k]+m44/m35*m45/m34*(U4_[i][j][k]*J[i][j][k]    -U4_[i][j+1][k]*J[i][j+1][k])    -m44/m35*m33/m45*(U4_[i][j+2][k]*J[i][j+2][k]-U4_[i][j+1][k]*J[i][j+1][k]))+
						thedac2*(U4_[i][j][k]  *J[i][j][k]  +m33/m24*m23/m34*(U4_[i][j+1][k]*J[i][j+1][k]-U4_[i][j][k]*J[i][j][k])        -m33/m24*m44/m23*(U4_[i][j-1][k]*J[i][j-1][k]-U4_[i][j][k]*J[i][j][k]))+
						thedac3*(U4_[i][j-1][k]*J[i][j-1][k]+m33/m12*m23/m13*(U4_[i][j-2][k]*J[i][j-2][k]-U4_[i][j-1][k]*J[i][j-1][k])+(1+m33/m23+m33/m13)*(U4_[i][j][k]*J[i][j][k]    -U4_[i][j-1][k]*J[i][j-1][k]));

					ML5[i-1][j][k-1] = thedac1*(U5_[i][j+1][k]*J[i][j+1][k]+m44/m35*m45/m34*(U5_[i][j][k]*J[i][j][k]    -U5_[i][j+1][k]*J[i][j+1][k])    -m44/m35*m33/m45*(U5_[i][j+2][k]*J[i][j+2][k]-U5_[i][j+1][k]*J[i][j+1][k]))+
						thedac2*(U5_[i][j][k]  *J[i][j][k]  +m33/m24*m23/m34*(U5_[i][j+1][k]*J[i][j+1][k]-U5_[i][j][k]*J[i][j][k])        -m33/m24*m44/m23*(U5_[i][j-1][k]*J[i][j-1][k]-U5_[i][j][k]*J[i][j][k]))+
						thedac3*(U5_[i][j-1][k]*J[i][j-1][k]+m33/m12*m23/m13*(U5_[i][j-2][k]*J[i][j-2][k]-U5_[i][j-1][k]*J[i][j-1][k])+(1+m33/m23+m33/m13)*(U5_[i][j][k]*J[i][j][k]    -U5_[i][j-1][k]*J[i][j-1][k]));


					MR1[i-1][j-1][k-1] = thedad1*(U1_[i][j+1][k]*J[i][j+1][k]+(1+m33/m34+m33/m35)*(U1_[i][j][k]*J[i][j][k]    -U1_[i][j+1][k]*J[i][j+1][k])+m33/m35*m34/m45*(U1_[i][j+2][k]*J[i][j+2][k]-U1_[i][j+1][k]*J[i][j+1][k]))+
						thedad2*(U1_[i][j][k]*J[i][j][k]     +m33/m24*m34/m23*   (U1_[i][j-1][k]*J[i][j-1][k]-U1_[i][j][k]*J[i][j][k])    -m33/m24*m22/m34*(U1_[i][j+1][k]*J[i][j+1][k]-U1_[i][j][k]*J[i][j][k]))+
						thedad3*(U1_[i][j-1][k]*J[i][j-1][k] +m22/m13*m12/m23*   (U1_[i][j][k]*J[i][j][k]    -U1_[i][j-1][k]*J[i][j-1][k])-m22/m13*m33/m12*(U1_[i][j-2][k]*J[i][j-2][k]-U1_[i][j-1][k]*J[i][j-1][k]));

					MR2[i-1][j-1][k-1] = thedad1*(U2_[i][j+1][k]*J[i][j+1][k]+(1+m33/m34+m33/m35)*(U2_[i][j][k]*J[i][j][k]    -U2_[i][j+1][k]*J[i][j+1][k])+m33/m35*m34/m45*(U2_[i][j+2][k]*J[i][j+2][k]-U2_[i][j+1][k]*J[i][j+1][k]))+
						thedad2*(U2_[i][j][k]*J[i][j][k]     +m33/m24*m34/m23*   (U2_[i][j-1][k]*J[i][j-1][k]-U2_[i][j][k]*J[i][j][k])    -m33/m24*m22/m34*(U2_[i][j+1][k]*J[i][j+1][k]-U2_[i][j][k]*J[i][j][k]))+
						thedad3*(U2_[i][j-1][k]*J[i][j-1][k] +m22/m13*m12/m23*   (U2_[i][j][k]*J[i][j][k]    -U2_[i][j-1][k]*J[i][j-1][k])-m22/m13*m33/m12*(U2_[i][j-2][k]*J[i][j-2][k]-U2_[i][j-1][k]*J[i][j-1][k]));

					MR3[i-1][j-1][k-1] = thedad1*(U3_[i][j+1][k]*J[i][j+1][k]+(1+m33/m34+m33/m35)*(U3_[i][j][k]*J[i][j][k]    -U3_[i][j+1][k]*J[i][j+1][k])+m33/m35*m34/m45*(U3_[i][j+2][k]*J[i][j+2][k]-U3_[i][j+1][k]*J[i][j+1][k]))+
						thedad2*(U3_[i][j][k]*J[i][j][k]     +m33/m24*m34/m23*   (U3_[i][j-1][k]*J[i][j-1][k]-U3_[i][j][k]*J[i][j][k])    -m33/m24*m22/m34*(U3_[i][j+1][k]*J[i][j+1][k]-U3_[i][j][k]*J[i][j][k]))+
						thedad3*(U3_[i][j-1][k]*J[i][j-1][k] +m22/m13*m12/m23*   (U3_[i][j][k]*J[i][j][k]    -U3_[i][j-1][k]*J[i][j-1][k])-m22/m13*m33/m12*(U3_[i][j-2][k]*J[i][j-2][k]-U3_[i][j-1][k]*J[i][j-1][k]));

					MR4[i-1][j-1][k-1] = thedad1*(U4_[i][j+1][k]*J[i][j+1][k]+(1+m33/m34+m33/m35)*(U4_[i][j][k]*J[i][j][k]    -U4_[i][j+1][k]*J[i][j+1][k])+m33/m35*m34/m45*(U4_[i][j+2][k]*J[i][j+2][k]-U4_[i][j+1][k]*J[i][j+1][k]))+
						thedad2*(U4_[i][j][k]*J[i][j][k]     +m33/m24*m34/m23*   (U4_[i][j-1][k]*J[i][j-1][k]-U4_[i][j][k]*J[i][j][k])    -m33/m24*m22/m34*(U4_[i][j+1][k]*J[i][j+1][k]-U4_[i][j][k]*J[i][j][k]))+
						thedad3*(U4_[i][j-1][k]*J[i][j-1][k] +m22/m13*m12/m23*   (U4_[i][j][k]*J[i][j][k]    -U4_[i][j-1][k]*J[i][j-1][k])-m22/m13*m33/m12*(U4_[i][j-2][k]*J[i][j-2][k]-U4_[i][j-1][k]*J[i][j-1][k]));

					MR5[i-1][j-1][k-1] = thedad1*(U5_[i][j+1][k]*J[i][j+1][k]+(1+m33/m34+m33/m35)*(U5_[i][j][k]*J[i][j][k]    -U5_[i][j+1][k]*J[i][j+1][k])+m33/m35*m34/m45*(U5_[i][j+2][k]*J[i][j+2][k]-U5_[i][j+1][k]*J[i][j+1][k]))+
						thedad2*(U5_[i][j][k]*J[i][j][k]     +m33/m24*m34/m23*   (U5_[i][j-1][k]*J[i][j-1][k]-U5_[i][j][k]*J[i][j][k])    -m33/m24*m22/m34*(U5_[i][j+1][k]*J[i][j+1][k]-U5_[i][j][k]*J[i][j][k]))+
						thedad3*(U5_[i][j-1][k]*J[i][j-1][k] +m22/m13*m12/m23*   (U5_[i][j][k]*J[i][j][k]    -U5_[i][j-1][k]*J[i][j-1][k])-m22/m13*m33/m12*(U5_[i][j-2][k]*J[i][j-2][k]-U5_[i][j-1][k]*J[i][j-1][k]));

					
				}

			}
		}
	}
	
//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

//#pragma omp parallel for private(k,k_)
	for (i = istart; i <= iend; i++) {
		for (k = 1, k_ = 2; k < nz; k++, k_++) { 

			ML1[i][1][k] = 0.5*(U1_[i+1][1][k_]*J[i+1][1][k_]+U1_[i+1][2][k_]*J[i+1][2][k_]);
			ML2[i][1][k] = 0.5*(U2_[i+1][1][k_]*J[i+1][1][k_]+U2_[i+1][2][k_]*J[i+1][2][k_]);
			ML3[i][1][k] = 0.5*(U3_[i+1][1][k_]*J[i+1][1][k_]+U3_[i+1][2][k_]*J[i+1][2][k_]);
			ML4[i][1][k] = 0.5*(U4_[i+1][1][k_]*J[i+1][1][k_]+U4_[i+1][2][k_]*J[i+1][2][k_]);
			ML5[i][1][k] = 0.5*(U5_[i+1][1][k_]*J[i+1][1][k_]+U5_[i+1][2][k_]*J[i+1][2][k_]);

			MR1[i][ny][k] = 0.5*(U1_[i+1][nyy][k_]*J[i+1][nyy][k_]+U1_[i+1][ny][k_]*J[i+1][ny][k_]);
			MR2[i][ny][k] = 0.5*(U2_[i+1][nyy][k_]*J[i+1][nyy][k_]+U2_[i+1][ny][k_]*J[i+1][ny][k_]);
			MR3[i][ny][k] = 0.5*(U3_[i+1][nyy][k_]*J[i+1][nyy][k_]+U3_[i+1][ny][k_]*J[i+1][ny][k_]);
			MR4[i][ny][k] = 0.5*(U4_[i+1][nyy][k_]*J[i+1][nyy][k_]+U4_[i+1][ny][k_]*J[i+1][ny][k_]);
			MR5[i][ny][k] = 0.5*(U5_[i+1][nyy][k_]*J[i+1][nyy][k_]+U5_[i+1][ny][k_]*J[i+1][ny][k_]);

		}
	}
// ============================================================================ //
/**** MUSCL 5th-order-end****/

	
//// ============================================ ////
			istart = 2;							  ////	
//// ============================================ ////
			iend = gend[myid]-1;	     		  ////
//// ============================================ ////

	#pragma omp parallel for private(\
	xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz,ETx,ETy,ETz,\
	_rho,_u,_v,_w,_U,_V,_W,__V,_VV,_P,_T,_C,_H,\
	rho_,u_,v_,w_,U_,V_,W_,V__,VV_,P_,T_,C_,H_,\
	rho,u,v,w,U,V,W,_V_,VV,H,C,P,T,\
	dU1,dU2,dU3,dU4,dU5,\
	beta,S,_S_,\
	temp,temp1,temp2,temp3,temp4,temp5,tempu1,tempu2,tempv1,tempv2,tempw1,tempw2,tempuvw,temph1,temph2,\
	d11,d12,d13,d14,d15,d21,d22,d23,d24,d25,d31,d32,d33,d34,d35,\
	d41,d42,d43,d44,d45,d51,d52,d53,d54,d55,\
	Fav1,Fav2,Fav3,Fav4,Fav5,\
	j,k\
	)

	/*---Y fluxes---*/
	for (i = istart; i <= iend; i++) {
		for (j = 1; j < nyy; j++) {
			for (k = 1; k < nz; k++) {

				xix=xidx_v[i][j][k];
				xiy=xidy_v[i][j][k];
				xiz=xidz_v[i][j][k];
				etx=etdx_v[i][j][k];
				ety=etdy_v[i][j][k];
				etz=etdz_v[i][j][k];          
				ztx=ztdx_v[i][j][k];
				zty=ztdy_v[i][j][k];
				ztz=ztdz_v[i][j][k];
				ETx=etx/(sqrt(etx*etx+ety*ety+etz*etz));
				ETy=ety/(sqrt(etx*etx+ety*ety+etz*etz));
				ETz=etz/(sqrt(etx*etx+ety*ety+etz*etz));

				/* lefr parameter */
				_rho = ML1[i][j][k];
				_u = ML2[i][j][k]/_rho;
				_v = ML3[i][j][k]/_rho;
				_w = ML4[i][j][k]/_rho;

				_U = xix*_u+xiy*_v+xiz*_w;
				_V = etx*_u+ety*_v+etz*_w;

				_W = ztx*_u+zty*_v+ztz*_w;


				__V = ETx*_u+ETy*_v+ETz*_w;


				_VV = _u*_u+_v*_v+_w*_w;
				_P = (ML5[i][j][k]-0.5*_rho*_VV)*(K-1);
				_T = _P/_rho;
				_C = sqrt(K*_P/_rho);
				_H = 0.5*_VV+_C*_C/(K-1);

				/* right parameter */
				rho_ = MR1[i][j][k];
				u_ = MR2[i][j][k]/rho_;
				v_ = MR3[i][j][k]/rho_;
				w_ = MR4[i][j][k]/rho_;

				U_ = xix*u_+xiy*v_+xiz*w_;
				V_ = etx*u_+ety*v_+etz*w_;
				W_ = ztx*u_+zty*v_+ztz*w_;

				V__ = ETx*u_+ETy*v_+ETz*w_;

				VV_ = u_*u_+v_*v_+w_*w_;
				P_ = (MR5[i][j][k]-0.5*rho_*VV_)*(K-1);
				T_ = P_/rho_;
				C_ = sqrt(K*P_/rho_);
				H_ = 0.5*VV_+C_*C_/(K-1);

				/* flux varibale */
				rho = 0.5*(_rho+rho_);
				u = 0.5*(_u+u_);
				v = 0.5*(_v+v_);
				w = 0.5*(_w+w_);

				U = 0.5*(_U+U_);
				V = 0.5*(_V+V_);
				W = 0.5*(_W+W_); 

				_V_ = 0.5*(__V+V__);

				VV = u*u+v*v+w*w;
				H = 0.5*(_H+H_);
				C = sqrt((H-0.5*VV)*(K-1));
				P = rho*C*C/K;
				T = P/rho;


				/* jump dU */
				dU1 = P_-_P;
				dU2 = u_-_u;
				dU3 = v_-_v;
				dU4 = w_-_w;
				dU5 = T_-_T;

				/* preconditioning */
				if (VV/C/C <= e)
					beta = e;
				else
					beta = VV/C/C;

				/*V=etx*u+ety*v+etz*w;*/
				//S=sqrt(V*V*(beta-1)+4*beta*C*C*(etx*etx+ety*ety+etz*etz));
				S=sqrt(V*V*(beta-1)*(beta-1)+4*beta*C*C*(etx*etx+ety*ety+etz*etz));
				/*_V_=ETx*u+ETy*v+ETz*w;*/
				_S_=sqrt(_V_*_V_*(beta-1)*(beta-1)+4*beta*C*C*(ETx*ETx+ETy*ETy+ETz*ETz));

				temp = 4*K*_S_*T*beta;
				temp1 = S+V+V*beta;
				temp2 = -S+V+V*beta;
				temp3 = _S_-_V_+_V_*beta;
				temp4 = _S_+_V_-_V_*beta;
				temp5 = _S_*_S_-(beta-1)*(beta-1)*_V_*_V_;
				tempu1 = 2*T*ETx*K*beta+u*_S_+u*(beta-1)*_V_;
				tempu2 = 2*T*ETx*K*beta-u*_S_+u*(beta-1)*_V_;
				tempv1 = 2*T*ETy*K*beta+v*_S_+v*(beta-1)*_V_;
				tempv2 = 2*T*ETy*K*beta-v*_S_+v*(beta-1)*_V_;
				tempw1 = 2*T*ETz*K*beta+w*_S_+w*(beta-1)*_V_;
				tempw2 = 2*T*ETz*K*beta-w*_S_+w*(beta-1)*_V_;
				tempuvw = 2*T*K*beta*(u*ETx+v*ETy+w*ETz);
				temph1 = tempuvw+H*_S_+H*_V_*beta-H*_V_; 
				temph2 = tempuvw-H*_S_+H*_V_*beta-H*_V_;


				d21 = 1/temp*(u*_S_*(4*(K-1)*beta*fabs(V)+fabs(temp2)+fabs(temp1))-(fabs(temp2)-fabs(temp1))*(2*T*ETx*K*beta+u*_V_*(beta-1)));
				d22 = 1/(2*temp)*(rho*(8*T*K*_S_*beta*(ETy*ETy+ETz*ETz)*fabs(V)+ETx*(temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1))));
				d23 = 1/(2*temp)*(rho*ETy*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d24 = 1/(2*temp)*(rho*ETz*(-8*T*K*_S_*beta*ETx*fabs(V)+temp3*tempu2*fabs(temp2)+temp4*tempu1*fabs(temp1)));
				d25 = -u*rho*fabs(V)/T;


				Fav2 = d21*dU1+d22*dU2+d23*dU3+d24*dU4+d25*dU5;


				MUL[i][j][k] = _u;
				MUR[i][j][k] = u_;

				MVL[i][j][k] = _v;
				MVR[i][j][k] = v_;

				MrhoL[i][j][k] = _rho;
				MrhoR[i][j][k] = rho_;

				MRUV[i][j][k] = 0.5*(_rho*_u*_v+rho_*u_*v_);

				MRoe_D[i][j][k] = 0.5*EpY[i][j][k]*Fav2;
				MdRL[i][j][k] = dU2;
				
			}
		}
	}




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


//// ============================ ////
	istart = 3;				      ////	
//// ============================ ////
	iend = gend[myid]+1;	      ////
//// ============================ ////

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
	_j,_k,__j,__k,j,k\
	)

	for (i = istart; i <= iend; i++) {
		for (__j = 0,_j = 1, j = 2; j <= nyy; __j++,_j++, j++) {
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

				/* viscous*/
				mu_E = mu_L;
				Pr_E = Pr_L;


				b1 = (4./3)*etx*etx+ety*ety+etz*etz;
				b5 = (1./3)*etx*ety;
				b7 = (1./3)*etx*etz;

				d1 = (4./3)*xix*etx+xiy*ety+xiz*etz;
				d5 = xix*ety-(2./3)*xiy*etx;
				d6 = xix*etz-(2./3)*xiz*etx;
				
				f1 = (4./3)*etx*ztx+ety*zty+etz*ztz;
				f7 = ety*ztx-(2./3)*etx*zty;
				f9 = etz*ztx-(2./3)*etx*ztz;
				

				/* Y viscid fluxes */
				ML2[i][j][k] = mu_E*(d1*du_dx+d5*dv_dx+d6*dw_dx+
					b1*du_dy+b5*dv_dy+b7*dw_dy+
					f1*du_dz+f7*dv_dz+f9*dw_dz)/J_v[i-1][j-1][k-1];


				/**** ES ****/
				MSXY[i][j-1][k] = ML2[i][j][k];
				/**** ES ****/

			}
		}
	}



//// ============================================ ////
		istart = 2;	         					  ////	
//// ============================================ ////
		iend = gend[myid]-1;	    			  ////
//// ============================================ ////

		for (j = 1; j < nyy; j++) {
			for (i = istart; i <= iend; i++) {
				for (k = 1; k < nz; k++) {

					ULm[j] = ULm[j]+MUL[i][j][k];
					URm[j] = URm[j]+MUR[i][j][k];

					VLm[j] = VLm[j]+MVL[i][j][k];
					VRm[j] = VRm[j]+MVR[i][j][k];

					rhoLm[j] = rhoLm[j]+MrhoL[i][j][k];
					rhoRm[j] = rhoRm[j]+MrhoR[i][j][k];

					RUVm[j] = RUVm[j]+MRUV[i][j][k];
					//SXYm[j] = SXYm[j]+MSXY[i][j][k];

					Roe_Dm[j] = Roe_Dm[j]+MRoe_D[i][j][k];
					dRLm[j] = dRLm[j]+MdRL[i][j][k];
					
				}
			}

		}



//// ============================================ ////
		 istart = 3;		     	              ////	
//// ============================================ ////
		iend = gend[myid];				    	  ////
//// ============================================ ////

		for (j = 1; j < nyy; j++) {
			for (i = istart; i <= iend; i++) {
				for (k = 2; k <= nz; k++) {

					SXYm[j] = SXYm[j]+MSXY[i][j][k];

				}
			}
		}



MPI_Comm comm;
comm=MPI_COMM_WORLD;
icount = Y_m;
idest = 0;
MPI_Reduce ((void*)&SXYm, (void*)&SXYm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);	

/**** Shear_stress output ****/
if (myid == 0) {
// =============================================================================================================== //

	FILE *fptr;
	sprintf(LESdata,"Shear_stress.dat");
	fptr = fopen(LESdata,"a");

	fprintf(fptr,"%f\t%f\n",SXYm0[1]/(nx-1)/(nz-1)-SXY_bottom,SXYm0[ny]/(nx-1)/(nz-1)-SXY_up);

	fclose(fptr);

	SXY_bottom = SXYm0[1]/(nx-1)/(nz-1);

	SXY_up = SXYm0[ny]/(nx-1)/(nz-1);

// =============================================================================================================== //
}
//**** Shear_stress output-end ****/


// ====================================== //
	if ( (step%statistic_step) == 0 ) {   //   
// ====================================== //

		MPI_Comm comm;
		comm=MPI_COMM_WORLD;
		
		//MPI_Barrier(MPI_COMM_WORLD);

		icount = Y_m;
		idest = 0;
		
		MPI_Reduce ((void*)&ULm, (void*)&ULm0, icount, MPI_DOUBLE ,MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&URm, (void*)&URm0, icount, MPI_DOUBLE ,MPI_SUM, idest, comm);

		MPI_Reduce ((void*)&VLm, (void*)&VLm0, icount, MPI_DOUBLE ,MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&VRm, (void*)&VRm0, icount, MPI_DOUBLE ,MPI_SUM, idest, comm);

		MPI_Reduce ((void*)&rhoLm, (void*)&rhoLm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&rhoRm, (void*)&rhoRm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);

		MPI_Reduce ((void*)&RUVm, (void*)&RUVm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&SXYm, (void*)&SXYm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);

		MPI_Reduce ((void*)&Roe_Dm, (void*)&Roe_Dm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&dRLm, (void*)&dRLm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		

// ========================= //
		if (myid == 0) {     //
// ========================= //

			double inv = 1./statistic_step/(nx-1)/(nz-1);

			for (j = 1; j <= ny; j++) {

				ULm0[j] = ULm0[j]*inv;
				URm0[j] = URm0[j]*inv;

				VLm0[j] = VLm0[j]*inv;
				VRm0[j] = VRm0[j]*inv;

				rhoLm0[j] = rhoLm0[j]*inv;
				rhoRm0[j] = rhoRm0[j]*inv;

				RUVm0[j] = RUVm0[j]*inv;
				SXYm0[j] = SXYm0[j]*inv;

				Roe_Dm0[j] = Roe_Dm0[j]*inv;
				dRLm0[j] = dRLm0[j]*inv;

			}


			FILE *fptrND;
			sprintf(LESdata,"ND""%0.5d"".dat",step);
			fptrND = fopen(LESdata,"w");

			for (j = 1; j <= ny; j++) {	fprintf(fptrND,"%.16f\n",ULm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",URm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VLm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",VRm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoLm0[j]);	}

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",rhoRm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",RUVm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",SXYm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",Roe_Dm0[j]); }

			for (j = 1; j <= ny; j++) { fprintf(fptrND,"%.16f\n",dRLm0[j]); }

			fclose(fptrND);

// ========================= //
		}					 //
// ========================= //

#pragma omp parallel for private(j)
		for (j = 1; j <= ny; j++) {

			ULm[j] = 0;
			URm[j] = 0;

			VLm[j] = 0;
			VRm[j] = 0;

			rhoLm[j] = 0;
			rhoRm[j] = 0;

			RUVm[j] = 0;
			SXYm[j] = 0;

			Roe_Dm[j] = 0;
			dRLm[j] = 0;

		}

		SXY_bottom = SXY_up = 0;

// ====================================== //
	}                                     //
// ====================================== //



	delete [] MUL;
	delete [] MUR;

	delete [] MVL;
	delete [] MVR;

	delete [] MrhoL;
	delete [] MrhoR;

	delete [] MRUV;
	delete [] MSXY;

	delete [] MRoe_D;
	delete [] MdRL;



/**** mean profile and turbulenct intensties ****/



//// ============================================ ////
		 istart = 3;		     	              ////	
//// ============================================ ////
		iend = gend[myid];				    	  ////
//// ============================================ ////

		for (j = 2; j <= ny; j++) {
			for (i = istart; i <= iend; i++) {
				for (k = 2; k <= nz; k++) {

					U = U2_[i][j][k]/U1_[i][j][k];
					V = U3_[i][j][k]/U1_[i][j][k];
					W = U4_[i][j][k]/U1_[i][j][k];     

					Um[j] = Um[j]+U;
					Vm[j] = Vm[j]+V;
					Wm[j] = Wm[j]+W;

					UUm[j] = UUm[j]+U*U;
					VVm[j] = VVm[j]+V*V;
					WWm[j] = WWm[j]+W*W;

					UVm[j] = UVm[j]+U*V;

				}
			}

		}



// ======================================================== //
	if ( (step%statistic_step) == 0) {					    //   
// ======================================================== //
		
		MPI_Comm comm;
		comm=MPI_COMM_WORLD;

		icount = Y_m;
		idest = 0;
		
		MPI_Reduce ((void*)&Um, (void*)&Um0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&Vm, (void*)&Vm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&Wm, (void*)&Wm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);

		MPI_Reduce ((void*)&UUm, (void*)&UUm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&VVm, (void*)&VVm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		MPI_Reduce ((void*)&WWm, (void*)&WWm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		
		MPI_Reduce ((void*)&UVm, (void*)&UVm0, icount, MPI_DOUBLE, MPI_SUM, idest, comm);
		
// ========================= //
		if (myid == 0) {     //
// ========================= //

			double inv = 1./statistic_step/(nx-1)/(nz-1);

			for (j = 2; j <= ny; j++) {

				Um0[j] = Um0[j]*inv;
				Vm0[j] = Vm0[j]*inv;
				Wm0[j] = Wm0[j]*inv;

				UUm0[j] = UUm0[j]*inv;
				VVm0[j] = VVm0[j]*inv;
				WWm0[j] = WWm0[j]*inv;

				UVm0[j] = UVm0[j]*inv;
				
			
				//printf("%f\t%d\n",Um0[j],myid);
			}


			FILE *fptr;
			sprintf(LESdata,"st""%0.5d"".dat",step);
			fptr = fopen(LESdata,"w");

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Um0[j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Vm0[j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",Wm0[j]); }

			for (j = 2; j <= ny; j++) {	fprintf(fptr,"%.16f\n",UUm0[j]);	}

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",VVm0[j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",WWm0[j]); }

			for (j = 2; j <= ny; j++) { fprintf(fptr,"%.16f\n",UVm0[j]); }
			fclose(fptr);

// ========================= //
		}					 //
// ========================= //

		for (j = 2; j <= ny; j++) {

			Um[j] = Vm[j] = Wm[j] = 0;
			UUm[j] = VVm[j] = WWm[j] = 0;
			UVm[j] = 0;

		}

// ======================================================== //
	}  												        //
// ======================================================== //

/**** mean profile and turbulenct intensties - end ****/


}
