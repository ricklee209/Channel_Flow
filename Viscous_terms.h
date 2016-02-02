#include "prm.h"
     
#define Cv 717.5                      /**** specific constant volume heat capacity J/(kg.k) ****/
#define mu_L 0.0000185
#define Pr_L 0.72
#define R 287

double Cp = Cv*K;                     /**** Cp=Cv*K ****/
double lambda_L = mu_L*Cv*K/Pr_L;
