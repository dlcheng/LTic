#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


/* 
Functions:
    1) do the zeldovich approximation 
    2) fill the particle data structure (position and velocity)
*/

void zel_fft_main_proc()
{
  int k1,k2,k3;
  int i,j,k;
  double k_s;
  double Dz = growth_factor(Z_start);
 
#pragma omp parallel shared(np1, np2, np3, L_box, pf, Dz, C1) private(i, j, k, k1, k2, k3, k_s)
{ 
  #pragma omp for
  for(i=0; i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
         {
         k1 = i - np1/2;
         k2 = j - np2/2;
         k3 = k - np3/2;
         k_s = 2.0 * Pi / L_box * sqrt(k1 * k1 + k2 * k2 + k3 * k3);
         
         pf[i][j][k].Re *= Dz * T_k(k_s) / C1;
         pf[i][j][k].Im *= Dz * T_k(k_s) / C1; 
         }
     }
}

 zel_fft_field(1);
 zel_fft_field(2);
 zel_fft_field(3);

}                                 /* end zel_fft_main_proc */



void zel_fft_field(int flag)
{
  int k1,k2,k3;
  int i,j,k;

#pragma omp parallel shared(np1, np2, np3, flag, L_box, pf, ptemp) private(i, j, k, k1, k2, k3)
{ 
  #pragma omp for  
  for(i=0; i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
        {
         k1 = i - np1/2;
         k2 = j - np2/2;
         k3 = k - np3/2;

         if(flag == 1)
           {
           ptemp[i][j][k].Re = -1.0 * pf[i][j][k].Im * 2.0 * Pi * k1 / L_box;
           ptemp[i][j][k].Im = pf[i][j][k].Re * 2.0 * Pi * k1 / L_box;  
           }

         if(flag == 2)
           {
           ptemp[i][j][k].Re = -1.0 * pf[i][j][k].Im * 2.0 * Pi * k2 / L_box;
           ptemp[i][j][k].Im = pf[i][j][k].Re * 2.0 * Pi * k2 / L_box;  
           }

         if(flag == 3)
           {
           ptemp[i][j][k].Re = -1.0 * pf[i][j][k].Im * 2.0 * Pi * k3 / L_box;
           ptemp[i][j][k].Im = pf[i][j][k].Re * 2.0 * Pi * k3 / L_box;  
           }
     
        }
     }
}

  fft_3d(0, ptemp);                  /* from k->x */
  flip_fft_cube(ptemp);              /* from z->F */
  fill_part_data(flag);
}                                    /* end zel_fft_field */

void fill_part_data(int flag)
{
 int i,j,k;
 int m = 0;
 
 double f1 = Omega_m * pow(1.0 + Z_start, 3) / (Omega_m * pow(1.0 + Z_start, 3) + Omega_v);
        f1 = pow(f1, 0.55);                   /* f1 here is very close to 1(0.9999986)
                                                 Omegam(z)^0.6 and Omegam(z)^0.55 are no difference.
                                              */

 double h2 = sqrt(Omega_m * (1.0 + Z_start) + Omega_v / pow(1.0 + Z_start, 2)) * 1e-3 / 3.0;

  for(i=0; i<np1; i++)
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
         {
         if(flag == 1)
           {
           p_part[m].pos[flag -1] = (i * L_box / (double) np1 + ptemp[i][j][k].Re) * 1000.0; /* lengh unit back to kpc/h */                                            
           p_part[m].vel[flag -1] = f1 * h2 * ptemp[i][j][k].Re * 3e5 * sqrt(1.0 + Z_start); /* speed unit km/s */                                           
           }

         if(flag == 2)
           {
           p_part[m].pos[flag -1] = (j * L_box / (double) np2 + ptemp[i][j][k].Re) * 1000.0; 
           p_part[m].vel[flag -1] = f1 * h2 * ptemp[i][j][k].Re * 3e5 * sqrt(1.0 + Z_start);
           }

         if(flag == 3)
           {
           p_part[m].pos[flag -1] = (k * L_box / (double) np3 + ptemp[i][j][k].Re) * 1000.0; 
           p_part[m].vel[flag -1] = f1 * h2 * ptemp[i][j][k].Re * 3e5 * sqrt(1.0 + Z_start);
           }
         m++;
         } 
}                                           /* end fill_part_data */
