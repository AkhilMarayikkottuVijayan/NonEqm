
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

PROGRAM VIBRATIONAL_RELAXATION
IMPLICIT NONE

REAL	:: Asr,Cs,Sigma_s
REAL	:: Temp,Pres,MW,thetaV,NA
REAL	:: Dens, m_sr,e,ns,c_s,A_sr
REAL	:: sig_s,tau_sr,tau_cs,tauV
! Millikan-White correlation

NA   =  6.0221409E+23 ! Avagadros number
e    = 2.718281828

WRITE (*,*) "Molecular weight ?"
READ  (*,*) MW
WRITE (*,*) "Pressure in atm ?"
READ  (*,*) Pres
WRITE (*,*) "Temperature ?"
READ  (*,*) Temp
WRITE (*,*) "Vibrational temperature ?"
READ  (*,*) thetaV


Dens	= Pres*1.01325E+5*MW/8314/Temp

m_sr	= 0.5*MW
A_sr	= Asr(m_sr, thetaV)
tau_sr	= (1/Pres)*e**(A_sr*(Temp**(-1/3)&
          &-0.015*m_sr**(1/4))-18.43)

! Parks correction

ns	= NA*Dens/MW
c_s	= Cs(Temp,MW)
sig_s	= Sigma_s(Temp)
tau_cs	= (sig_s*c_s*ns)**(-1)

! Net vibrational relaxation time

tauV	= tau_sr+tau_cs

WRITE (*,*) "Vibrational relaxation time :   ", tauV
END PROGRAM
