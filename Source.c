#include "udf.h"

DEFINE_SOURCE(vt_source,c,t,dS,eqn)
{
real	vt,msr,MW,theta;
real	xr,tau_sr,Asr,temp;	
real	pres,e,C_temp,tau_vs;
real	density, Ys,tempV;
real	evsT,evsTV;

e	= 2.718281828	;
MW	= 28.0134	;
theta	= 3.395E+3	;
msr	= MW*0.5	;
temp	= C_T(c,t)	;
pres	= C_P(c,t)	;
density	= C_R(c,t)	;
Ys	= C_YI(c,t,1)	;
Asr	= 1.16E-3*pow(msr,0.5)*pow(theta,4/3);
C_temp	= Asr*(pow(temp,-1/3)-0.015*pow(msr,0.25))-18.42;
tau_sr	= (1/pres)*pow(e,C_temp)/pres;

/* Since single species N2 ,
 * Only one avergaring term */
tau_vs	= tau_sr	;
tempV	= C_UDMI(c,t,1)	;
evsT	= theta*Ru/(MW*(pow(e,theta/temp)-1));
evsTV   = theta*Ru/(MW*(pow(e,theta/tempV)-1));

/*  Landau-Teller formula */
vt 	= density*Ys*(evsT-evsTV)/tau_vs;

return vt;
}
