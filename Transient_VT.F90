FUNCTION Asr(M_sr, theta_vs)
    REAL    :: Asr
    REAL    :: M_sr,theta_va

    Asr = 1.16E-3*SQRT(M_sr)*theta_vs**(4/3)
    RETURN
END FUNCTION

FUNCTION Cs(T,M_ws)
    REAL ::Cs
    REAL :: T,M_ws
    
    Cs = SQRT(8*8314*T/3.141/M_ws)
    RETURN
END FUNCTION

FUNCTION Sigma_s(T)
    REAL :: Sigma_s
    REAL :: T
 
    Sigma_s = 1E-20*(50000*T**(-1))**2
END FUNCTION

FUNCTION eV(theta,MW,T)
   REAL ::eV
   REAL ::theta,MW,T
   REAL :: Ru,e
   Ru = 8314
   e  = 2.718281828
   eV = theta*Ru/(MW*(e**(theta/T)-1))
END FUNCTION

FUNCTION fTR(MW,SV,Dens)
   REAL :: fTR
   REAL :: MW,SV,Dens
   REAL :: dof, Ru,cp
   Ru = 8314
   dof = 7

   cp = 0.5*Ru*dof/MW

   fTR = SV/Dens/cp
END FUNCTION

FUNCTION fTV(MW,t,Tv,SV,Dens)
  REAL :: fTV
  REAL :: MW,t,Tv,SV,Dens
  REAL ::Ru, A,B,e
  Ru = 8314
  e = 2.718281828
  
  A = (SV*(e**(t/Tv)-1)**2)/(e**(t/Tv)*Dens)
  B = MW*Tv**2/(Ru*t**2) 
  
  fTV = A*B
END FUNCTION

PROGRAM VIBRATIONAL_RELAXATION
IMPLICIT NONE

INTEGER :: I,N
REAL	:: Asr,Cs,Sigma_s,eV,fTR,fTV
REAL	:: Temp,Pres,MW,thetaV,NA
REAL	:: Dens, m_sr,e,ns,c_s,A_sr
REAL	:: sig_s,tau_sr,tau_cs,tauV
REAL	:: dt,Tv,Tv0,Ttr,Sv
REAL    :: Time

! Millikan-White correlation

NA   =  6.0221409E+23 ! Avagadros number
e    = 2.718281828

WRITE (*,*) "Molecular weight ?"
READ  (*,*) MW
WRITE (*,*) "Pressure in atm ?"
READ  (*,*) Pres
WRITE (*,*) "Initial Translational Temperature ?"
READ  (*,*) Temp
WRITE (*,*) "Initial Vibrational Temperature ?"
READ  (*,*) Tv0
WRITE (*,*) "Vibrational temperature ?"
READ  (*,*) thetaV
WRITE (*,*) "Numb"
READ  (*,*) N

Ttr = Temp
Tv  = Tv0
dt = thetaV/N
Time = 0.0

DO 100 I = 1,N

Time    = Time+dt
Dens	= Pres*1.01325E+5*MW/8314/Ttr

m_sr	= 0.5*MW
A_sr	= Asr(m_sr, thetaV)
tau_sr	= (1/Pres)*e**(A_sr*(Ttr**(-1/3)&
          &-0.015*m_sr**(1/4))-18.43)

! Parks correction

ns	= NA*Dens/MW
c_s	= Cs(Temp,MW)
sig_s	= Sigma_s(Ttr)
tau_cs	= (sig_s*c_s*ns)**(-1)

! Net vibrational relaxation time

tauV	= tau_sr+tau_cs


! Source term !
Sv =  Dens*(eV(thetaV,MW,Ttr)-eV(thetaV,MW,Tv))/tauV

! Temperature updation
Ttr = Ttr-fTR(MW,Sv,Dens)*dt
Tv  = Tv+fTV(MW,thetaV,Tv,Sv,Dens)*dt
WRITE (2,*) Time,Ttr,Tv
100 CONTINUE
END PROGRAM
