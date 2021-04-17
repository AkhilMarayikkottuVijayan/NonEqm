#include "udf.h"

DEFINE_SOURCE(vt_source,c,t,dS,eqn)
{
real	vt,msr,MW,theta;
real	xr,tau_sr,Asr,temp;	
real	pres,e,C_temp,tau_vs;
real	density, Ys,tempV;
real	evsT,evsTV,Ru,NA;
real	tau_cs,sig_s,C_s,n_s;

e	= 2.718281828	;
NA	= 6.0221476E23	;
MW	= 30.0061	;
theta	= 2.817E+3	;
Ru	= 8314.0	;
msr	= MW*0.5	;
temp	= C_T(c,t)	;
pres	= C_P(c,t)	;
density	= C_R(c,t)	;
Ys	= C_YI(c,t,1)	;
Asr	= 1.16E-3*pow(msr,0.5)*pow(theta,4/3);
C_temp	= Asr*(pow(temp,-1/3)-0.015*pow(msr,0.25))-18.42;
tau_sr	= (1/pres)*pow(e,C_temp)/pres;

/* Parks correction */

sig_s	= 1E-20*pow(50000/temp,2);
C_s	= pow(8*Ru*temp/(3.14125*MW),0.5);
n_s	= NA*density*Ts/MW;
tau_cs	= pow(sig_s*C_s*n_s,-1);

/* Since we have only one diatomic species at this 
 * point, we have to just add the two relaxation 
 * time steps to get tau_vs */
tau_vs 	= tau_sr+tau_cs;

tempV	= C_UDMI(c,t,0)	;
evsT	= theta*Ru/(MW*(pow(e,theta/temp)-1));
evsTV   = theta*Ru/(MW*(pow(e,theta/tempV)-1));

/*  Landau-Teller formula */
vt 	= density*Ys*(evsT-evsTV)/tau_vs;

return vt;
}
