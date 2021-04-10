#include "udf.h"

DEFINE_PROPERTY(evn2,c,t)
{
	real enerV;
	real thetaV =3395.0;
	real MW = 28.0134;
	real Ru = 8.314;
	real e = 2.717;
	real temp;

        temp  = C_T(c,t);
	enerV = thetaV*Ru/MW/(pow(e,(thetaV/temp))-1);
	return enerV;
}

DEFINE_DIFFUSIVITY(mun2,c,t,i)
{
	real mu;
	real MW = 28.0134;
	real CsR;
	real temp;
        real A_22 = -7.6303990E-13 ;
        real B_22 = 1.6878089E-11  ;
        real C_22 = -1.4004234E-10 ;
        real D_22 = 2.1427708E-7   ;

        temp  = C_T(c,t);
	CsR   = D_22*temp*((A_22*log(temp)+B_22)*log(temp)+C_22);
        mu    = 2.6693E-6*3.14125*pow(MW*temp,0.5)/CsR ;
	return mu;
}
