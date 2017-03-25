#ifndef ALLVAR_H
#define ALLVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_fft_complex.h>

#include "define.h"

typedef struct gadget_head GADGET_HEAD;
typedef struct part_data PART_DATA;
typedef struct grid GRID;

struct gadget_head
{  
   unsigned int Npart[6];
   double Massarr[6];
   double Time;
   double Redshift;
   int FlagSfr;
   int FlagFeedback;
   int Nall[6];
   int  FlagCooling;
   int NumFiles;
   double BoxSize;
   double Omega0;
   double OmegaLambda;
   double HubbleParam;
   int FlagAge;
   int FlagMetals;
   int NallHW[6];
   int Flag_entr_ics;
   char unused[60];          
};

struct part_data
{
  float pos[3];   
                               /* 
                                  pos[0] = x coordinate
                                  pos[1] = y coordinate
                                  pos[2] = z coordinate
                               */
  float vel[3];
                                /*
                                   vel[0] = V_x coordinate
                                   vel[1] = V_y coordinate
                                   vel[2] = V_z coordinate
                                */
};

struct grid
{
  float Re;                      /* real part of the complex number */
  float Im;                      /* imaginary part of the complex number */
};

extern double norm_A;            /* the normalization factor */
extern double mean_phi_square;   /* mean phi_square get from the data */
extern GRID ***pf;               /* pointer of the cubic field */
extern GRID ***ptemp;            /* temporary field used for zeldovich approximation to generate the velocity field */
extern PART_DATA *p_part;        /* pointer of the particle data structure */


extern double C1;                /* 3*omegam*H0*h0/2   */
extern double C2;                /* sqrt(N*N*nor_A/2/V) */

extern FILE *fp_g;               /* output file name */


/* parameters */

extern double Z_start;
extern double L_box; 
extern int Ran_seed;
extern double F_nl;
extern double h0;

extern double Omega_m;
extern double Omega_v;
extern double Omega_b;

extern double ns;
extern double Sigma_8;

extern int np1;
extern int np2;
extern int np3;

#ifdef WD_ST
extern double WD_alpha;
#endif

#endif

