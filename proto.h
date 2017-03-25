#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void init_global();
void init_param_format();
void init_memory();
void init_constants();
void init_sampling(int k1, int k2, int k3, gsl_rng *my_rg);

void normalization();
double f_k(double k);
double f_q(double q, void *param);
double growth_factor(double z);
double g_square_factor(double z);
double T_k(double k);
#ifdef EH_T_K
void init_eh_t_k();
double eh_t_k(double k_0);
double tf_press_less(double q, double a , double b);
#endif

void local_type_transfer();
void fft_3d(int flag, GRID ***f);
void copy_to_fft_array(int i, int j, int k, GRID ***f,double *fft_data_local);
void copy_from_fft_array(int i, int j, int k, GRID *** f, double *fft_data_local);
void flip_fft_cube(GRID ***f);
void local_transfer();

void alloc_3d_array(int flag);
void free_3d_array(GRID *** array);
void alloc_part_data();
void free_part_data();
void free_all();

void zel_fft_main_proc();
void zel_fft_field(int flag);
void fill_part_data(int flag);

void output_gadget_file();
void make_gadget_head(GADGET_HEAD *gh);
void write_gadget_file(GADGET_HEAD * gh);

void warn_and_end(char *s);
void state(char *s);
