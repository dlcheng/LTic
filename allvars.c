#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"

 double norm_A;                 /* the normalization factor */
 double mean_phi_square;        /* mean phi_square get from the data */
 GRID ***pf;                    /* pointer of the cubic field */
 GRID ***ptemp;                 /* temporary field used when doing zeldovich 
                                   approximation to generate the velocity field */
 PART_DATA *p_part;             /* pointer of the particle data structure */



 double C1;                     /* 3*omegam*H0*h0/2   */
 double C2;                     /* sqrt(N*N*nor_A/2/V) */


 FILE *fp_g;                    /* output file */

/* parameters*/

 double Z_start;
 double L_box; 
 int Ran_seed;
 double F_nl;
 double h0;

 double Omega_m;
 double Omega_v;
 double Omega_b;

 double ns;
 double Sigma_8;

 int np1;
 int np2;
 int np3;
 
#ifdef WD_ST
double WD_alpha;
#endif

