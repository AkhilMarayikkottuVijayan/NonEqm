#include "udf.h"

DEFINE_UDS_FLUX(massflux_plus_diffusiveflux,f,t,i)
{
    cell_t c0, c1 = -1;
    Thread *t0, *t1 = NULL;
    real NV_VEC(psi_vec), NV_VEC(A), flux = 0.0;
    c0 = F_C0(f,t);
    t0 = F_C0_THREAD(f,t);
    F_AREA(A, f, t);
    
 if (BOUNDARY_FACE_THREAD_P(t)) /*Most face values will be available*/
    {
      real dens;

      if (NNULLP(THREAD_STORAGE(t,SV_DENSITY)))
        dens = F_R(f,t);   
      else
        dens = C_R(c0,t0); 

      NV_DS(psi_vec,  =, F_U(f,t), F_V(f,t), F_W(f,t), *, dens);

      flux = NV_DOT(psi_vec, A); /* flux through Face */
    }
  else
    {
      c1 = F_C1(f,t);       /* Get cell on other side of face */
      t1 = F_C1_THREAD(f,t); 

      NV_DS(psi_vec,  =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
      NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

      flux = NV_DOT(psi_vec, A)/2.0; /* Average flux through face */
    }

  /* ANSYS FLUENT will multiply the returned value by phi_f (the scalar's
     value at the face) to get the "complete'' advective term.  */

  return flux;
}
