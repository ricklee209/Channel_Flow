



#include "Resolution.h"

extern int X_np;


void Grid
(
// ======================================================== //
int myid,

double (*X_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*Y_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*Z_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ======================================================== //
);


void Initial_condition
(
// ======================================================== //
int myid,

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

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ======================================================== //
);


// ============================================================================ //
void Jacobian
(
int myid,

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
);
// ============================================================================ //




// ============================================================================ //
void Jacobian_u
(
int myid,

double (*J_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_u)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
);
// ============================================================================ //




// ============================================================================ //
void Jacobian_v
(
int myid,

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

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 

);
// ============================================================================ //




// ============================================================================ //
void Jacobian_w
(
int myid,

double (*J_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*xidx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*xidz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*etdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*etdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*ztdx_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdy_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*ztdz_w)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
);
// ============================================================================ //




// ============================================================================ //
void Driven_force
(
int myid,

double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
);
// ============================================================================ //




// ============================================================================ //
void X_boundary_condition
(
double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
);
// ============================================================================ //




// ============================================================================ //
void Y_boundary_condition
(
int myid,
double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
);
// ============================================================================ //




// ============================================================================ //
void Z_boundary_condition
(
int myid,
double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
);
// ============================================================================ //




// ============================================================================ //
void Flux_X
(
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
);
// ============================================================================ //



// ============================================================================ //
void Flux_Y
(
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
);
// ============================================================================ //



// ============================================================================ //
void Flux_Z
(
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
);
// ============================================================================ //






// ============================================================================ //
void Flux_X_limiter
(
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
);
// ============================================================================ //



// ============================================================================ //
void Flux_Y_limiter
(
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
);
// ============================================================================ //



// ============================================================================ //
void Flux_Z_limiter
(
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
);
// ============================================================================ //






void Viscous_terms
// ============================================================================ //
(
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
// ============================================================================ //
);



void Dual_time_stepping
(
// ============================================================================ //
int myid,

double deltaT,

double *er1,
double *er2,
double *er3,
double *er4,
double *er5,


int (*Nijk)[2] = new int[(X_np-3)+ny+nz+7][2],
int (*ijk)[3] = new int[(X_np-6)*Y_out*Z_out*3][3],


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
;




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
;


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
);




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

//double (*EpZ)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]
// ============================================================================ //
);


void Filter
(
// ============================================================================ //
int myid,


double Roe_criterion,
double E_high,
double E_low,


double (*U1_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U2_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U3_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U4_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*U5_)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 

double (*Y_point)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*EpX)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],
double (*EpY)[Y_m][Z_m] = new double[X_np][Y_m][Z_m], 
double (*EpZ)[Y_m][Z_m] = new double[X_np][Y_m][Z_m],

double (*J)[Y_m][Z_m] = new double[X_np][Y_m][Z_m]

// ============================================================================ //
);

void Hyperplane 
(
// =========================================================== //
int myid,

int (*Nijk)[2] = new int[(X_np-3)+ny+nz+7][2],

int (*ijk)[3] = new int[(X_np-6)*Y_out*Z_out*3][3]

// =========================================================== //
);



void RK3NODT
(
// ============================================================================ //
int myid,

int RK,

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
double (*MR5)[Y_m][Z_m] = new double[X_np][Y_m][Z_m] 
/**** MR = RHS-end ****/


// ============================================================================ //
);

