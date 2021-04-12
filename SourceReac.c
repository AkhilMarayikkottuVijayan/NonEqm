#include "udf.h"

DEFINE_SOURCE(n2_dis_src1,c,t,dS,eqn)
{
real	A,Tcf,Tcb,b,Ed,Ru;
real	kf,kb,e,temp,tempV;
real	Keq;

e	= 2.718281828	;
Ru	= 8314.0	;
A	= 7.0E+18	;
b	= -1.6		;
Ed	= 9.4E+8	;

temp	= C_T(c,t)	;
tempV	= C_UDMI(c,t,1) ; 
/* forward reaction rate */
Tcf	= pow(temp*tempV,0.5);
kf	= A*pow(Tcf,b)*pow(e,
	  -1*Ed/(Ru*Tcf));

/* backward reaction rate */
Tcb 	= temp;
kb	= kf/Keq;
}
