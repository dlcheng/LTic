#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

/*
Functions:
   1)
   2)
*/

double f_k(double k)              /* defined as T_k^2  * k^(2+ns) * W_k^2 */
{
 double tk = T_k(k);
 double k_ns = pow(k , 2+ns);
 double wk,k8;

 k8 = 8.0 * k;
 tk = tk * tk;

 if(k8 < 1e-3)
   wk = 1.0 / 3.0 * (1.0 - 0.1 * k8 * k8);
 else 
   wk = sin(k8) / k8 / k8 / k8 - cos(k8) / k8 / k8;

 wk = wk * wk;

 return tk * k_ns * wk;
}                                  /* end f_k */

double f_q(double q, void *param)
{
 return f_k(1.0 / q - 1.0) / q / q;
}                                  /* end f_q */

double growth_factor(double z)
{
 double g_square = g_square_factor(z);
 double omega_m_z = Omega_m * pow(1.0 + z, 3) / g_square;
 double omega_v_z = Omega_v / g_square;

 return 1.0 / (1.0 + z) * 2.5 * omega_m_z 
        / (pow(omega_m_z, 4.0/7.0) - omega_v_z + (1.0 + omega_m_z * 0.5) * (1.0 + omega_v_z / 70.0));

}                                /* end growth_factor */

double g_square_factor(double z)
{
  return Omega_m * pow(1.0 + z, 3) + (1.0 - Omega_m - Omega_v) * pow(1.0 + z, 2) + Omega_v;
}                                /* end g_factor */

double T_k(double k)
{
 double tk;

#ifdef BBKS_T_K
 double b1, b2;
 double T_sig = Omega_m * h0 * exp(-1.0 * Omega_b - sqrt(2.0 * h0) * Omega_b / Omega_m);
 double q = k / T_sig;

 double a0 = 2.34;
 double a1 = 3.89;
 double a2 = 16.19;
 double a3 = 5.46;
 double a4 = 6.71;
    
 if(q <= 1e-8)
    tk = 1.0;                    /* for k = 0 return 1 */
 else 
   {
     b1 = 1.0 + a1 * q + a2 * a2 * q * q + a3 * a3 * a3 * q * q * q + a4 * a4 * a4 * a4 * q * q * q * q;
     b1 = 1.0 / pow (b1 , 0.25);
     b2 = log(1.0 + a0 * q)/ a0 / q;
     tk = b1 * b2;
   }
#endif

#ifdef EH_T_K
if(k > 1e-8)
 tk = eh_t_k(k);
else
 tk = 1.0;
#endif

#ifdef WD_ST
 tk *= pow(1+pow(WD_alpha * k, 2.25), -3.08);
#endif

 return tk;

}                                 /* end T_k */

void normalization()
{
int WORKSIZE = 10000;
 double d_growth = growth_factor(0.0);
 d_growth *= d_growth;

 double c = 4.5 * d_growth / Pi_square;

 gsl_function F;
 gsl_integration_workspace *workspace;
 double result,abserr;

 workspace = gsl_integration_workspace_alloc(WORKSIZE);
 F.function = &f_q;

 gsl_integration_qag(&F,0,1,0,1.0e-6, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
 gsl_integration_workspace_free(workspace);

 norm_A = Sigma_8 * Sigma_8 / c / result;
}                                   /* end normalization */

