#include "udf.h"

DEFINE_VR_RATE(vol_reac_rate,c,t,r,wk,yk,rate,rr_t)
{
  real ci, prod;
  int i;

/* Calculate Arrhenius reaction rate  */

  prod = 1.;

/*  for(i = 0; i < r->n_reactants; i++)
    {
       ci     = C_R(c,t) * yk[r->reactant[i]] / wk[r->reactant[i]];
       prod  *= pow(ci, r->exp_reactant[i]);
       printf ("yk:  %f \n", yk[r->reactant[i]]);
       printf ("wk:  %f \n", wk[r->reactant[i]]);
    }*/
 
   if (!strcmp(r->name, "N2 + N2 <-> 2N + N2"))
      {
      /* printf ("Reactant : %f \n ", yk[r->reactant[0]]);
       printf ("Reaction rate : %f \n", *rate);
       printf ("Reaction name: %s\n", r->name);
       printf ("A : %f \n", r->A);
       printf ("Activation energy : %f \n", r->E);
       printf ("Pre-exponent : %f \n", r->b);*/
      } 
   else if (!strcmp(r->name, "O2 + O2 <-> 2O + O2"))
      {
      /* printf ("Reactant : %f \n ", yk[r->reactant[0]]);
       printf ("Reaction rate : %f \n", *rate);
       printf ("Reaction name: %s\n", r->name);
       printf ("A : %f \n", r->A);
       printf ("Activation energy : %f \n", r->E);
       printf ("Pre-exponent : %f \n", r->b);*/
      } 
   else if (!strcmp(r->name, "O2 + N2 <-> 2O + N2"))
      {
      /* printf ("Reactant : %f \n ", yk[r->reactant[0]]);
 *        printf ("Reaction rate : %f \n", *rate);
       printf ("Reaction name: %s\n", r->name);
       printf ("A : %f \n", r->A);
       printf ("Activation energy : %f \n", r->E);
       printf ("Pre-exponent : %f \n", r->b);*/
      } 

  /* *rate = r->A * exp( - r->E / (UNIVERSAL_GAS_CONSTANT * C_T(c,t)))*pow(C_T(c,t), r->b) * prod;
  *rr_t = *rate;*/
  /* No "return..;" value. */
}

