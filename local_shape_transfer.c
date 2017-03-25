#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"



/* 
Functions:
    1) do one FFT to real space
    2) do another one to k space again
*/

void local_type_transfer()
{
  fft_3d(0, pf);             /* k->x */
  flip_fft_cube(pf);         /* to the real value of the field */
  local_transfer();          /* apply the equation in real field */
  flip_fft_cube(pf);         /* flip back to z value for FFT */
  fft_3d(1, pf);             /* FFT to k space again */
  
}                            /* end local_type_transfer */


void fft_3d(int flag, GRID ***f)
{                             /* flag = 0 means k->x, 
                                 flag = 1 means x->k 
                                 f is the fft cube
                              */


 int i,j,k;
 double *fft_data_local;

/* using OpenMP to speed up */
#pragma omp parallel shared(np1, np2, np3, f, flag) private(j, k, fft_data_local)
{ 
 fft_data_local = (double *) malloc(2 * np1 * sizeof(double));

 #pragma omp for 
 for(j=0; j<np2; j++)      
  { 
   for(k=0; k<np3; k++)
     {
     copy_to_fft_array(-1,j,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,np1);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,np1);

     copy_from_fft_array(-1,j,k,f,fft_data_local);     
     }
  }

  free(fft_data_local);
}

#pragma omp parallel shared(np1, np2, np3, f, flag) private(k, i, fft_data_local)
{ 
  fft_data_local = (double *) malloc(2 * np1 * sizeof(double));
  
  #pragma omp for 
  for(k=0; k<np3; k++)
  {
   for(i=0; i<np1; i++)
     {
     copy_to_fft_array(i,-1,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,np2);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,np2);

     copy_from_fft_array(i,-1,k,f,fft_data_local);     
     }
   }

   free(fft_data_local);
}

#pragma omp parallel shared(np1, np2, np3, f, flag) private(i, j, fft_data_local)
{	
 fft_data_local = (double *) malloc(2 * np1 * sizeof(double));
 
 #pragma omp for 
 for(i=0; i<np1; i++)
   {
   for(j=0; j<np2; j++)
     {
     copy_to_fft_array(i,j,-1,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,np3);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,np3);

     copy_from_fft_array(i,j,-1,f,fft_data_local);     
     }

   }

   free(fft_data_local);
}

}                             /* end three_d_fft*/

void copy_to_fft_array(int i, int j, int k, GRID ***f, double *fft_data_local)
{
 int m;
 if(i == -1)
   {
    for(m=0; m<np1; m++)
      {
      fft_data_local[2*m] = f[m][j][k].Re;
      fft_data_local[2*m+1] = f[m][j][k].Im;
      }
   }

 if(j == -1)
   {
    for(m=0; m<np2; m++)
      {
      fft_data_local[2*m] = f[i][m][k].Re;
      fft_data_local[2*m+1] = f[i][m][k].Im;
      }
   }

 if(k == -1)
   {
    for(m=0; m<np3; m++)
      {
      fft_data_local[2*m] = f[i][j][m].Re;
      fft_data_local[2*m+1] = f[i][j][m].Im;
      }
   }
}                             /* end copy_to_fft_array */

void copy_from_fft_array(int i, int j, int k, GRID *** f, double * fft_data_local)
{
 int m;
 if(i == -1)
   {
    for(m=0; m<np1; m++)
      {
      f[m][j][k].Re = fft_data_local[2*m];
      f[m][j][k].Im = fft_data_local[2*m+1];
      }
   }

 if(j == -1)
   {
    for(m=0; m<np2; m++)
      {
      f[i][m][k].Re = fft_data_local[2*m];
      f[i][m][k].Im = fft_data_local[2*m+1];
      }
   }

 if(k == -1)
   {
    for(m=0; m<np3; m++)
      {
      f[i][j][m].Re = fft_data_local[2*m];
      f[i][j][m].Im = fft_data_local[2*m+1];
      }
   }
}                             /* end copy_from_fft_array */

void flip_fft_cube(GRID ***f)
{
  int i,j,k;

#pragma omp parallel shared(np1, np2, np3, f) private(i, j, k)
{  
  #pragma omp for
  for(i=0; i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
         {
         if((i+j+k) % 2 == 1)
           f[i][j][k].Re *= -1.0;

//         f[i][j][k].Im = 0.0;
         }
     }
}
}                             /* end flip_fft_cube */

void local_transfer()
{
  int i,j,k;
  double k_min = 2 * Pi / L_box;
  double k_max = sqrt(1) * Pi * np1 / L_box;
  double factor_ns = 1.0 / (ns -1);
  double factor_A = C1 * C1 * norm_A / 2.0 / Pi_square;

#ifdef CUT_OFF
  double k5 = 0.0;
  double mean_phi_square_op_1 = 0.0;
  double mean_phi_square_op_2, mean_phi_square_op_3;
  
#pragma omp parallel shared(np1, np2, np3, pf, L_box, ns, mean_phi_square_op_1, k5) private(i, j, k)
{ 
  double k1, k2, k3, k4;
  double my_mean_phi_square_op_1 = 0.0;
  double my_k5 = 0.0;
  
  #pragma omp for  
  for(i=0; i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
        {
        my_mean_phi_square_op_1 += pf[i][j][k].Re * pf[i][j][k].Re;

        if(i != np1/2 && j != np2/2 && k != np3/2)
          {
          k1 = i - np1/2;
          k2 = j - np2/2;
          k3 = k - np3/2;
          k4 = sqrt(k1 * k1 + k2 * k2 + k3 * k3) * 2 * Pi / L_box;
          my_k5 += pow(k4, ns-4);
	      }
	    }
	 }   
	 
   #pragma omp atomic
   mean_phi_square_op_1 += my_mean_phi_square_op_1;
   #pragma omp atomic
   k5 += my_k5;
}	    
  mean_phi_square_op_1 = mean_phi_square_op_1 / (double) np1 / (double) np2 / (double) np3;	    
  mean_phi_square_op_2 = factor_A * pow(2* Pi / L_box, 3) * k5 / 4 / Pi;
  mean_phi_square_op_3 = -1.0 * factor_A * factor_ns * pow(k_min, ns-1);
#endif	    

 mean_phi_square = factor_ns * factor_A  * (pow(k_max, ns-1) - pow(k_min, ns-1));
 
#ifdef CUT_OFF 
  printf(" mean_phi_square:\n");
  printf(" (0) %.8e with both high and low k cut off\n", mean_phi_square);
  printf(" (1) %.8e average over real space points\n", mean_phi_square_op_1);
  printf(" (2) %.8e average over k space points\n", mean_phi_square_op_2);
  printf(" (3) %.8e with only low k cut off\n", mean_phi_square_op_3);
#endif     

#pragma omp parallel shared(np1, np2, np3, pf, F_nl, mean_phi_square) private(i, j, k)
{  
  #pragma omp for 
  for(i=0; i<np1; i++)
    {
    for(j=0; j<np2; j++)
      for(k=0; k<np3; k++)
         pf[i][j][k].Re += F_nl * (pf[i][j][k].Re * pf[i][j][k].Re - mean_phi_square);    
    }
}
}                              /* end local_transfer */
