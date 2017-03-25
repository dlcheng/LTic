#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"


/* 
Functions:
    1) convert the particle date to gadget file.
    2) do split if initial condition is too large (not finished).
*/

void output_gadget_file()
{

 GADGET_HEAD gh;
 char filename[128];
 char file_path[100];
 
 strcpy(file_path, Outputfold);
 if(file_path[strlen(file_path) - 1] != '/')
  sprintf(file_path,"%s""%s", Outputfold, "/");

 sprintf(filename,"%s%s%.0f%s%.0f%s%.0f%s",file_path,"gadget_fnl_",F_nl,"_z_",Z_start,"_L_",L_box,".dat");
 fp_g = fopen(filename,"wb+");
 rewind(fp_g);
 make_gadget_head(&gh);

 write_gadget_file(&gh);

}                                             /* end output_gadget_file */


void make_gadget_head(GADGET_HEAD *gh)
{
 double GSI = 6.67428e-11;
 double M_S = 1.98892e30;
 double MPC = 3.0857e22;
 double boxsize = L_box;
 double numd = np1;
 double mass_p = 3.* Omega_m * MPC * (boxsize / numd) * (boxsize / numd) * (boxsize / numd) / (8. * Pi * GSI * M_S); 
 int i;
 unsigned int N_total = np1 * np2 * np3;

 
 gh->Npart[0] = 0;                             /* type 0 stands for gas particles */
 gh->Npart[1] = np1 * np2 * np3;               /* this is dark matter particles  */

 for(i=2; i<=5; i++)
   {
   gh->Npart[i] = 0;
   }

 gh->Massarr[0] = 0.;                          /* Type 0 gas particle has no mass now */
 gh->Massarr[1] = mass_p;

 for(i=2; i<=5; i++)
   {
   gh->Massarr[i]=0.;
   }

 gh->Time = 1.0 / (1.0 + Z_start);             /* the startup time can be set in Gadget2 parameter file */
 gh->Redshift = Z_start;
 gh->FlagSfr = 0;
 gh->FlagFeedback = 0;
 gh->Nall[0] = 0;                              
 gh->Nall[1] = N_total;

 for(i=2; i<=5; i++)
   {
   gh->Nall[i] = 0;
   }

 gh->FlagCooling = 0;
 gh->NumFiles = 1;                              /* may used when split the initial files */
                             
 gh->BoxSize = L_box * 1000.0 ;          
 gh->Omega0 = Omega_m;
 gh->OmegaLambda = Omega_v;
 gh->HubbleParam = h0;                          /* H0 = h0 (100km/s/Mpc) */
 gh->FlagAge = 0;
 gh->FlagMetals = 0;
 
 for(i=0;i<=5;i++)
   {
   gh->NallHW[i]=0;
   }

 gh->Flag_entr_ics = 0;

}                                                /*end make_gadget_head */

void write_gadget_file(GADGET_HEAD * gh)
{
 unsigned int block_size;
 unsigned int N_total = np1 * np2 * np3;
 unsigned int i;

 block_size = 256;                                  /* the head */
 fwrite(&block_size, sizeof(block_size), 1, fp_g);
 fwrite(gh, 256, 1, fp_g); 
 fwrite(&block_size, sizeof(block_size), 1, fp_g);

 block_size = N_total * 3 * sizeof(float);          /* the position */
 fwrite(&block_size, sizeof(block_size), 1, fp_g);

 for(i=0; i<N_total; i++)
   fwrite(p_part[i].pos, 3 * sizeof(float), 1, fp_g);

 fwrite(&block_size, sizeof(block_size), 1, fp_g);

 fwrite(&block_size, sizeof(block_size), 1, fp_g);  /* the velocity */

 for(i=0; i<N_total; i++)
   fwrite(p_part[i].vel, 3 * sizeof(float), 1, fp_g);   

 fwrite(&block_size, sizeof(block_size), 1, fp_g);

 block_size = N_total * sizeof(int);                /* the ID */
 fwrite(&block_size, sizeof(block_size), 1, fp_g);
 
 for(i=0; i<N_total; i++)                           /* Very important, the particle ID starts from 0, 
                                                       this is linked to the halo shape analysis code */
   fwrite(&i, sizeof(int), 1, fp_g);

 fwrite(&block_size, sizeof(block_size), 1, fp_g);

}                                                   /* end write_gadget_fiel */ 

