#include "udf.h"

/* Initializing the fields using UDFs */
DEFINE_INIT(tv_init,d)
{
        cell_t c        ;
        Thread  *t      ;
        real temp,theta,Ru,e;
        real MW ;
 
        theta	= 3395.0;
        e	= 2.718281;
        Ru	= 8314.0;
        MW      = 28.134; 
/* Loop over all cell threads in the domain */
  thread_loop_c(t,d)
        {
           /* loop over all cells */
                begin_c_loop_all(c,t)
                {
                temp	      = C_T(c,t);
                C_UDMI(c,t,0) = temp;
                C_UDSI(c,t,0) = theta*Ru/(MW*(pow(e,theta/temp)-1));
                }
                end_c_loop_all(c,t)
        }
}

/* Adjust UDFs for computing the vibrational
 * temperature and manipulating the field
 * during ongoing simulations */
DEFINE_ADJUST(tv_comp,d)
{
        Thread *t;
        cell_t c ;
        real enerV;
        real thetaV =3395.0;
        real MW = 28.0134;
        real Ru = 8314.0;
        real Tv;

        /* Loop over all the thread */
        thread_loop_c (t,d)
        {
                begin_c_loop (c,t)
                {
                enerV = C_UDSI(c,t,0);
                Tv = thetaV*pow(log((thetaV*Ru/(MW*enerV))+1),-1);
                C_UDMI(c,t,0) = Tv;
               /* printf ("Pressure	: %f \n", C_P(c,t));
                printf ("Temperature    : %f \n", C_T(c,t));
                printf ("Density        : %f \n", C_R(c,t));
                printf ("Y0             : %f \n", C_YI(c,t,0));
                printf ("Y1             : %f \n", C_YI(c,t,1));
                printf ("Y2             : %f \n", C_YI(c,t,2));
                printf ("Y3             : %f \n", C_YI(c,t,3));
                printf ("Y4             : %f \n", C_YI(c,t,4));*/
                //printf ("Vibrational temperature : %f \n", C_UDMI(c,t,0));
                }
                end_c_loop (c,t)
        }
}

/* Diffusion coefficent used in the UDS equation 
 * computed using UDFs using Hilfelshders
 * formulation for temperatrure dependent viscosity */
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
        C_UDMI(c,t,2) = mu;
        return mu;
}

/* Vibrational-Translational relaxation source 
 * terms computed using Landaue-Teller model
 * is coded up as a source UDF               */
DEFINE_SOURCE(vt_source,c,t,dS,eqn)
{
real    vt,msr,MW,theta;
real    xr,tau_sr,Asr,temp;
real    pres,e,C_temp,tau_vs;
real    density, Ys,tempV;
real    evsT,evsTV,Ru,NA;
real    tau_cs,sig_s,C_s,n_s;

e       = 2.718281828   ;
NA      = 6.0221476E23  ;
MW      = 28.0134       ;
theta   = 3.395E+3      ;
Ru      = 8314.0        ;
msr     = MW*0.5        ;
temp    = C_T(c,t)      ;
pres    = C_P(c,t)      ;
density = C_R(c,t)      ;
Ys      = C_YI(c,t,1)   ;
Asr     = 1.16E-3*pow(msr,0.5)*pow(theta,4/3);
C_temp  = Asr*(pow(temp,-1/3)-0.015*pow(msr,0.25))-18.42;
tau_sr  = (101325.0/pres)*pow(e,C_temp);

/* Parks correction */

sig_s   = 1E-20*pow(50000/temp,2);
C_s     = pow(8*Ru*temp/(3.14125*MW),0.5);
n_s     = NA*density*Ys/MW;
tau_cs  = pow(sig_s*C_s*n_s,-1);

/* Since we have only one diatomic species at this 
 *  * point, we have to just add the two relaxation 
 *   * time steps to get tau_vs */
tau_vs  = tau_sr+tau_cs;

tempV   = C_UDMI(c,t,0) ;
evsT    = theta*Ru/(MW*(pow(e,theta/temp)-1));
//evsTV   = theta*Ru/(MW*(pow(e,theta/tempV)-1));
evsTV   = C_UDSI(c,t,0);
/*  Landau-Teller formula */
vt      = density*Ys*(evsT-evsTV)/tau_vs;
//printf ("Tau relaxation source term : %f \n", tau_vs);
C_UDMI(c,t,1) = vt;
C_UDMI(c,t,3) = tau_vs;
//printf ("LT - VT : %f \n ", vt);
//vt =0.0;
return vt;
}

DEFINE_SOURCE(e_vt_source,c,t,dS,eqn)
{
real    evt;
evt     = -1*C_UDMI(c,t,1);
//printf ("evt : %f \n", evt);
//evt     = 0.0;
return evt;
}

/*
/* Reaction source terms into the energy equation
 * and the UDS equation for non-equilibrium */
/*DEFINE_SOURCE(n2_dis_src1,c,t,dS,eqn)
{
real    A,Tcf,Tcb,b,Ed,Ru;
real    kf,kb,e,temp,tempV;
real    Keq,Yn,Yn2,rho;
real    rhoN2,rhoN,Mn2,Mn;
real    wn2_1,theta,evs,src1;
real	A1,A2,A3,A4,A5;
real	ndens,NA,expon;

e       = 2.718281828   ;
Ru      = 8314.0        ;
A       = 7.0E+18       ;
b       = -1.6          ;
Ed      = 9.4E+8        ;
Mn2     = 28.0134       ;
Mn      = 14.0067       ;
theta   = 3395.0        ;
NA	= 6.022E+23	;

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
ndens	= rho*1000.0*NA/Mn2;
if ( ndens <= 1E+14)
{A1=3.4907;A2=0.8313;A3=4.0978;A4=-12.728;A5=0.0748};
if ( 1E+14 < ndens <= 1E+15)
{A1=2.0723;A2=1.3897;A3=2.0617;A4=-11.828;A5=0.1510};
if ( 1E+15 < ndens <= 1E+16)
{A1=1.6060;A2=1.5732;A3=1.3923;A4=-11.533;A5=-0.0045};
if ( 1E+16 < ndens <= 1E+17)
{A1=1.5351;A2=1.6061;A3=1.2993;A4=-11.494;A5=-0.0069};
if ( 1E+17 < ndens <= 1E+18)
{A1=1.4766;A2=1.6291;A3=1.2153;A4=-11.457;A5=-0.0094};
if ( 1E+18 < ndens)
{A1=1.4766;A2=1.6291;A3=1.2153;A4=-11.457;A5=-0.0094};

expon	= A1*temp/10000+A2+A3*log(10000/temp)+
          A4*10000/temp+A5*pow(10000/temp,2);
Keq     = pow(e,expon);        //Parks eqn
kb      = kf/Keq;

/* species production */
rhoN2   = Yn2*rho;
rhoN    = Yn*rho;
wn2_1   = (kf*pow(rhoN2/Mn2,2)-
          kb*(rhoN2/Mn2)*pow(rhoN/Mn,2));

/* vibrational energy from specis production */
evs     = theta*Ru/(Mn2*
          (pow(e,theta/tempV)-1));
src1    = evs*wn2_1;

return src1;
}*/

