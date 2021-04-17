#include "udf.h"

DEFINE_DIFFUSIVITY(mun2,c,t,i)
{
	real mu;
	real MW = 28.0134;
	real CsR;
	real temp;
        real A_22 = -7.6303990E-3 ;
        real B_22 = 1.6878089E-1  ;
        real C_22 = -1.4004234E0 ;
        real D_22 = 2.1427708E3   ;

        temp  = C_T(c,t);
	CsR   = D_22*pow(temp,((A_22*log(temp)+B_22)*log(temp)+C_22));
        mu    = 2.6693E-6*3.14125*pow(MW*temp,0.5)/CsR ;
	return mu;
}

DEFINE_DIFFUSIVITY(muo2,c,t,i)
{
        real mu;
        real MW = 31.9988;
        real CsR;
        real temp;
        real A_22 = -6.2931612E-3 ;
        real B_22 = 1.4624645E-1  ;
        real C_22 = -1.3006927E0 ;
        real D_22 = 1.8066892E+3   ;

        temp  = C_T(c,t);
        CsR   = D_22*pow(temp,((A_22*log(temp)+B_22)*log(temp)+C_22));
        mu    = 2.6693E-6*3.14125*pow(MW*temp,0.5)/CsR ;
        return mu;

DEFINE_DIFFUSIVITY(muo2,c,t,i)
{
        real mu;
        real MW = 30.0061;
        real CsR;
        real temp;
        real A_22 = -7.4942466E-3 ;
        real B_22 = 1.6626193E-1  ;
        real C_22 = -1.4107027E0 ;
        real D_22 = 2.3097604E+3   ;

        temp  = C_T(c,t);
        CsR   = D_22*pow(temp,((A_22*log(temp)+B_22)*log(temp)+C_22));
        mu    = 2.6693E-6*3.14125*pow(MW*temp,0.5)/CsR ;
        return mu;

