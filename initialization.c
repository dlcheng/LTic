#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"


/*
Function:
   1) create the data structure 
   2) sampling according to certain distribution
*/

void init_global()
{ 
  init_param_format();
  
#ifdef EH_T_K
  init_eh_t_k();
#endif

#ifdef WD_ST 
  WD_alpha = 0.189 * pow(M_WD, -0.858) * pow(Omega_m/0.26, -0.136) * pow(H0/0.7, 0.692);
#endif

  normalization();
  init_memory();
  init_constants();

  int i,j,k;

  
#pragma omp parallel shared(np1, np2, np3, L_box, pf, C1, C2) private(i, j, k)
{    
  gsl_rng * my_rg;
  int my_id;
  
  my_id = omp_get_thread_num();  
  my_rg = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(my_rg, Ran_seed + 2 * my_id);	            /* start-up seed */
  
  #pragma omp for 
  for(i=0; i<= np1 / 2; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)        
        init_sampling(i,j,k,my_rg);                     /* random_generator is a shared variable, 
                                                           making OpenMP not working here */
    }
    
   gsl_rng_free(my_rg);
}        
  state(" Complete sampling.");

#pragma omp parallel shared(np1, np2, np3, pf) private(i, j, k)
{ 
  #pragma omp for 
  for(i=(np1 / 2 + 1); i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
        {
        pf[i][j][k].Re = pf[(np1-i) % np1][(np2-j) % np2][(np3-k) % np3].Re;
        pf[i][j][k].Im = -1.0 * pf[(np1-i) % np1][(np2-j) % np2][(np3-k) % np3].Im;
        }
     }
}

}                                                           /* end init_global */

void init_param_format()
{
 Z_start = Z_START;
 L_box = L_BOX;
 Ran_seed = RAN_SEED;
 F_nl = F_NL;
 h0 = H0;

 Omega_m = OMEGA_M;
 Omega_v = OMEGA_V;
 Omega_b = OMEGA_B;

 ns = NS;
 Sigma_8 = SIGMA_8;

 np1 = NP1;
 np2 = NP2;
 np3 = NP3;
}                                                             /* init_param_format */

void init_memory()
{
  alloc_3d_array(1);
  alloc_3d_array(2);
  state(" Complete memory allocation for FFT data.");

  alloc_part_data();
  state(" Complete memory allocation for particle data.");
}                                                           /* end init_memory */


void init_constants()
{
 double N_total = np1 * np2 * np3;
 C1 = 1.5 * Omega_m * 1e-6 / 9.0;
 C2 = N_total * sqrt(norm_A / 2.0 / L_box / L_box / L_box);
 
 state(" Complete constants calculation.");
}                                                           /* end init_constants */

void init_sampling(int k1, int k2, int k3, gsl_rng *my_rg_1)
{
 int k_1 = k1 - np1/2;
 int k_2 = k2 - np2/2;
 int k_3 = k3 - np3/2;

 double k = 2.0 * Pi / L_box * sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

 double C3;
 
 if(k_1 == 0 && k_2 == 0 && k_3 == 0)
   {
   pf[k1][k2][k3].Re = 0.0;
   pf[k1][k2][k3].Im = 0.0;                                  /* k=0, the sampling is 0 */
   }
 else
   {
   C3 = pow(k, (ns-4.0)/2.0);
   double sig = C1 * C2 * C3;
   pf[k1][k2][k3].Re = gsl_ran_gaussian(my_rg_1, sig);
   pf[k1][k2][k3].Im = gsl_ran_gaussian(my_rg_1, sig);
   }
}                                                            /* end init_sampling */
 
#ifdef PEN_CORRECTION



#endif
