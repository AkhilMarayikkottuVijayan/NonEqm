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