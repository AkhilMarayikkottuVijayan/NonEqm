#include "udf.h"

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
