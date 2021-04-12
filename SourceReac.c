#include "udf.h"

DEFINE_SOURCE(n2_dis_src1,c,t,dS,eqn)
{
real	A,Tcf,Tcb,b,Ed,Ru;
real	kf,kb,e,temp,tempV;
real	Keq,Yn,Yn2,rho;
real	rhoN2,rhoN,Mn2,Mn;
real	wn2_1,theta,evs,src1;

e	= 2.718281828	;
Ru	= 8314.0	;
A	= 7.0E+18	;
b	= -1.6		;
Ed	= 9.4E+8	;
Mn2	= 28.0134	;
Mn	= 14.0067	;
theta	= 3395.0	;

rho	= C_R(c,t)	;
temp	= C_T(c,t)	;
tempV	= C_UDMI(c,t,1) ;
Yn2	= C_YI(c,t,1)	;
Yn	= C_YI(c,t,2)	;
/* forward reaction rate */
Tcf	= pow(temp*tempV,0.5);
kf	= A*pow(Tcf,b)*pow(e,
	  -1*Ed/(Ru*Tcf));

/* backward reaction rate */
Tcb 	= temp;
Keq	=         //Parks eqn
kb	= kf/Keq;

/* species production */
rhoN2	= Yn2*rho;
rhoN	= Yn*rho;
wn2_1	= (kf*pow(rhoN2/Mn2,2)-
	  kb*(rhoN2/Mn2)*pow(rhoN/Mn,2));

/* vibrational energy from specis production */
evs	= theta*Ru/(Mn2*
	  (pow(e,theta/tempV)-1));
src1	= evs*wn2_1;

return src1;
}

DEFINE_SOURCE(n2_dis_src2,c,t,dS,eqn)
{
real    A,Tcf,Tcb,b,Ed,Ru;
real    kf,kb,e,temp,tempV;
real    Keq,Yn,Yn2,rho;
real    rhoN2,rhoN,Mn2,Mn;
real    wn2_2,evs,theta,src2;

e       = 2.718281828   ;
Ru      = 8314.0        ;
A       = 3.0E+19       ;
b       = -1.6          ;
Ed      = 9.4E+8        ;
Mn2     = 28.0134       ;
Mn      = 14.0067       ;
theta	= 3395.0	;

rho     = C_R(c,t)      ;
temp    = C_T(c,t)      ;
tempV   = C_UDMI(c,t,1) ;
Yn2     = C_YI(c,t,1)   ;
Yn      = C_YI(c,t,2)   ;
/* forward reaction rate */
Tcf     = pow(temp*tempV,0.5);
kf      = A*pow(Tcf,b)*pow(e,
          -1*Ed/(Ru*Tcf));

/* backward reaction rate */
Tcb     = temp;
Keq     =         //Parks eqn
kb      = kf/Keq;

/* species production */
rhoN2   = Yn2*rho;
rhoN    = Yn*rho;
wn2_2   = (kf*(rhoN2/Mn2)*(rhoN/Mn)-
          kb*pow(rhoN/Mn,3));
/* vibrational energy from specis production */
evs     = theta*Ru/(Mn2*
          (pow(e,theta/tempV)-1));
src2    = evs*wn2_2;

return src2;
}
~  
