#include "udf.h"

/* Initializing the fields using UDFs */
DEFINE_INIT(tv_init,d)
{
        cell_t c        ;
        Thread  *t      ;

/* Loop over all cell threads in the domain */
  thread_loop_c(t,d)
        {
           /* loop over all cells */
                begin_c_loop_all(c,t)
                {
                C_UDMI(c,t,0) = 1200.0;
                }
                end_c_loop_all(c,t)
        }
}

/* Adjust UDFs for computing the vibrational
 * temperature and manipulating the field
 * during ongoing simulations */
DEFINE_ADJUST(tv_comp,d)
{
        Thread *t;
        cell_t c ;
        real enerV;
        real thetaV =3395.0;
        real MW = 28.0134;
        real Ru = 8314.0;
        real Tv;

        /* Loop over all the thread */
        thread_loop_c (t,d)
        {
                begin_c_loop (c,t)
                {
                enerV = C_UDSI(c,t,0);
                Tv = thetaV*pow(log((thetaV*Ru/MW/(enerV))+1),-1);
                C_UDMI(c,t,1) = Tv;
                }
                end_c_loop (c,t)
        }
}

/* Diffusion coefficent used in the UDS equation 
 * computed using UDFs using Hilfelshders
 * formulation for temperatrure dependent viscosity */
DEFINE_DIFFUSIVITY(mun2,c,t,i)
{
        real mu;
        real MW = 28.0134;
        real CsR;
        real temp;
        real A_22 = -7.6303990E-3 ;
        real B_22 = 1.6878089E-1  ;
        real C_22 = -1.4004234E0 ;
        real D_22 = 2.1427708E3   ;

        temp  = C_T(c,t);
        CsR   = D_22*pow(temp,((A_22*log(temp)+B_22)*log(temp)+C_22));
        mu    = 2.6693E-6*3.14125*pow(MW*temp,0.5)/CsR ;
        return mu;
}

