#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"

int main()
{
 omp_set_num_threads(N_thread);
 
 state("Start initializing ...");
 init_global();  
 state("Complete \n");
 
 state("Start local type transforming ...");
 local_type_transfer();
 state("Complete\n");
 
#ifndef CUT_OFF
 state("Start Zeldo'vich approximation ...");
 zel_fft_main_proc();
 state("Complete\n");

 state("Begin to write output file ...");
 output_gadget_file();
 state("Compelete\n");
#endif

 free_all();

 return 0;
}                           /* end main */


