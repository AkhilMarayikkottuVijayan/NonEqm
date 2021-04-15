#include "udf.h"

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

