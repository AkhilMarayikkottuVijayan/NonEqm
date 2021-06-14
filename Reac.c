#include "udf.h"

DEFINE_VR_RATE(vol_reac_rate,c,t,r,wk,yk,rate,rr_t)
{
  real ci, prod;
  int i;

/* Calculate Arrhenius reaction rate  */

  prod = 1.;

/*for(i = 0; i < r->n_reactants; i++)
    {
       ci     = C_R(c,t) * yk[r->reactant[i]] / wk[r->reactant[i]];
       prod  *= pow(ci, r->exp_reactant[i]);
    }
  
   *rate = r->A * exp( - r->E / (UNIVERSAL_GAS_CONSTANT * C_T(c,t)))*pow(C_T(c,t), r->b) * prod;
  *rr_t = *rate;*/

  /* No "return..;" value. */
}

