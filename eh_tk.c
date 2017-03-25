#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"

#ifdef EH_T_K
static double  y,
               omhh,          /* The density of CDM and baryons, in units of critical dens,
                                 multiplied by the square of the Hubble constant, in units
                                 of 100 km/s/Mpc */
               obhh,
               Tcmb,          /* The temperature of the CMB in Kelvin, 2.728(4) by COBE  */
               theta_cmb,     /* Tcmb in units of 2.7 K */
               z_equality,    /* Redshift of matter-radiation equality, really 1+z */
               k_equality,    /* Scale of equality, in Mpc^-1 */
               z_drag,        /* Redshift of drag epoch */
               R_drag,        /* Photon-baryon ratio at drag epoch */
               R_equality,    /* Photon-baryon ratio at equality epoch */
               sound_horizon, /* Sound horizon at drag epoch, in Mpc */
               k_silk,        /* Silk damping scale, in Mpc^-1 */
               alpha_c,       /* CDM suppression */
               beta_c,        /* CDM log shift */
               alpha_b,       /* Baryon suppression */
               beta_b,        /* Baryon envelope shift */
               f_baryon,      /* The fraction of baryons to CDM */
               beta_node;

void init_eh_t_k()
{
  f_baryon = Omega_b / Omega_m;
  omhh = Omega_m * H0 * H0;
  obhh = omhh * f_baryon;
  Tcmb = T_CMB;
  theta_cmb = Tcmb / 2.7;
 
/* main variables */ 
  z_equality = 2.50e4 * omhh * pow(theta_cmb, -4.0);                       /* different form Mpgrafic */
  k_equality = 0.0746 * omhh * pow(theta_cmb, -2.0); 
  z_drag = 0.313 * pow(omhh, -0.419) * (1.0 + 0.607 * pow(omhh, 0.674));
  z_drag = 1e0 + z_drag * pow(obhh, (0.238 * pow(omhh, 0.223)));
  z_drag = 1291e0 * pow(omhh, 0.251) / (1e0 + 0.659 * pow(omhh, 0.828)) * z_drag;
  R_drag = 31.5 * obhh * pow(theta_cmb, -4.0) *1000e0 / (z_drag);          /* different form Mpgrafic */
  R_equality = 31.5 * obhh * pow(theta_cmb, -4.0) * 1000e0 / (z_equality); /* different form Mpgrafic */
  sound_horizon = 2.0 / 3.0 / k_equality * sqrt(6.0/R_equality) * 
                  log((sqrt(1.0+ R_drag)+sqrt(R_drag+R_equality)) / (1.0+sqrt(R_equality)));
  k_silk = 1.6 * pow(obhh, 0.52) * pow(omhh, 0.73) * (1e0 + pow(10.4 * omhh, -0.95));
  alpha_c = pow(46.9*omhh, 0.670) * (1e0+pow(32.1*omhh, -0.532));
  alpha_c = pow(alpha_c, -1.0 * f_baryon);
  alpha_c = alpha_c * pow(pow(12.0*omhh, 0.424)*(1e0 + pow(45.0*omhh, -0.582)), pow(-1.0 * f_baryon, 3.0));
  beta_c = 0.944 / (1+pow(458.0*omhh, -0.708));
  beta_c = 1.0 + beta_c * (pow(1.0-f_baryon, pow(0.395*omhh, -0.0266)) - 1e0);
  beta_c = 1.0 / beta_c;
  y = (1e0 + z_equality)/(1e0 + z_drag);
  alpha_b = y * (-6.0 * sqrt(1.0+y)+(2.0+3.0*y) * log((sqrt(1.0+y)+1.0) / (sqrt(1.0+y)-1.0)));
  alpha_b = 2.07 * k_equality * sound_horizon * pow(1.0+R_drag, -0.75) * alpha_b;
  beta_b = 0.5 + f_baryon + (3.0-2.0*f_baryon)* sqrt(pow(17.2*omhh, 2.0) + 1e0);
  beta_node = 8.41 * pow(omhh, 0.435);	
	
}

double eh_t_k(double k_0)                                                   /* k here is in unit of h*Mpc^-1 */
{
  double q, k, ks;
  double tf_cdm, s_tilde, tf_baryon, tf_full;
  
  k = k_0 * H0;            /* changed from k_0/H0 */
  q = k / 13.41 / k_equality;
  ks = k * sound_horizon;

  tf_cdm = 1.0 / (1.0 + pow(ks/5.4, 4.0));
  tf_cdm = tf_cdm * tf_press_less(q, 1.0, beta_c) + (1.0-tf_cdm) * tf_press_less(q, alpha_c, beta_c);
  s_tilde = sound_horizon / pow(1.0 + pow(beta_node/ks, 3.0), 1.0/3.0); 
  tf_baryon = tf_press_less(q,1.0,1.0) / (1.0 + pow(ks/5.2, 2.0));
  tf_baryon = tf_baryon + alpha_b / (1.0+ pow(beta_b/ks, 3.0)) * exp(-1.0 * pow(k/k_silk, 1.4));
  tf_baryon = tf_baryon * (sin(k * s_tilde) / (k * s_tilde));
  tf_full = f_baryon * tf_baryon + (1-f_baryon) * tf_cdm;

  return tf_full;

}

double tf_press_less(double q, double a , double b)
{
  double result;
  double C;
  
  C = 14.2 / a + 386.0 / (1.0 + 69.9 * pow(q, 1.08));
  result = log(exp(1.0) + 1.8*b*q);
  result = result / (result  + C * q * q);	
	
  return result;
}
#endif
