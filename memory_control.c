#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*
Functions:
   1) create 3-D array
   2) create particle data array
   3) free 3-D array
   4) free particle data array
*/

void alloc_3d_array(int flag)
{
 int i,j;

 if(flag == 1)
  {
  pf = (GRID ***) malloc(np1 * sizeof(GRID **));

  for(i=0; i<np1 ; i++)
    pf[i] = (GRID **) malloc(np2 * sizeof(GRID *));
   
  for(i=0; i<np1 ; i++)
    for(j=0; j<np2 ; j++)
       pf[i][j] = (GRID *) malloc(np3 * sizeof(GRID));
  }
 else
  {
  ptemp = (GRID ***) malloc(np1 * sizeof(GRID **));

  for(i=0; i<np1 ; i++)
    ptemp[i] = (GRID **) malloc(np2 * sizeof(GRID *));
   
  for(i=0; i<np1 ; i++)
    for(j=0; j<np2 ; j++)
       ptemp[i][j] = (GRID *) malloc(np3 * sizeof(GRID));
  }
}                              /* end create_3d_array */

void free_3d_array(GRID *** array)
{
 int i,j;
  
 for(i=0; i<np1; i++)
   for(j=0; j<np2 ; j++)
      free(array[i][j]);

 for(i=0; i<np1; i++)
   free(array[i]);

 free(array);
}                               /* end free_3d_array */

void alloc_part_data()
{
 if(( p_part = (PART_DATA *) malloc(np1 * np2 * np3 * sizeof(PART_DATA)) ) == NULL)
   warn_and_end("Fail to allocate memory to store the particle data.");
}                               /* end create_particle_data */


void free_part_data()
{
 free(p_part);
}                               /* end free_particle_data */


void free_all()
{
 free_part_data();
 free_3d_array(ptemp);
 free_3d_array(pf);
#ifndef CUT_OFF
 fclose(fp_g);
#endif
}                               /* end free_all */
