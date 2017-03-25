/* 
Attention:
(1) the code uses natural unit in the calculation, 1 length = 1 Mpc/h
    while the Gadget uints for output:
                         1 length = 1kpc/h
                         1 velocity = 1 km/s, peculiar velocity / sqrt(a)
                         1 mass = 10^10 M_sun /h
(2) the ID of the partilces shoud start from 0, to be in consistent with the halo shape code.

*/
#define N_thread   10                          /* number of cores used */
/* PARAMETERS */
#define Z_START    100                         /* redshift where the initial condition is made */
#define RAN_SEED   1988100                     /* the seed of random number generator          */
#define F_NL       0                           /* F_nl */
#define L_BOX      256                         /* box size in unit of Mpc/h     */  

/* NOW ONLY ALLOW np1=np2=np3=2^n  */ 
#define NP1        256                         /* gird number in x direction */
#define NP2        256                         /* grid number in y direction */
#define NP3        256                         /* grid number in z direction */

/* MAKE SURE YOUR HAVE THIS FOLDER */                   
#define Outputfold "./Result/"    

/* COSMOLOGICAL PARAMETERS */
#define H0         0.7                         /* H(today) = H0 * 100km/s/Mpc */
#define T_CMB      2.728                       /* CMB average temperature in unit of Kelvin */
#define OMEGA_M    0.3                         /* the matter energy density in unit of critical density today */
#define OMEGA_B    0.024 / H0 / H0             /* the baryon energy density in unit of critical density today */
#define OMEGA_V    0.7                         /* the vaccum energy density in unit of cirtical density today */
#define NS         0.96                        /* the power law index of primordial power spectrum */
#define SIGMA_8    0.8                         /* fluctuation of matter excess fraction in sphere of 8 Mpc/h today */
#ifdef WD_ST
#define M_WD       0.5                         /* mass of sterile neutrino, fitting valid for range [0.3,15] in unit 
                                                  of keV */
#endif
                                   
/* CONSTANTS */
#define Pi         3.1415926535897932384626433832795028842            
#define Pi_square  9.86960440108935861883449099987615113531           
                
