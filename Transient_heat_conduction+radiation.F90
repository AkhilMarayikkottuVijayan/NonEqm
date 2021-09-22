PROGRAM STEMP
IMPLICIT NONE
INTEGER		:: I,J,L
DOUBLE PRECISION :: Q,t,X,Y,T1,T2,rho,Cp
DOUBLE PRECISION :: Ti,k,sig,eps,dX,r1,r2,PI,dT
DOUBLE PRECISION :: Tw1,Tw2,Qc,Qr,R_slab,R_shell
DOUBLE PRECISION :: Tn(100),ddX,dtC,d,error
OPEN (99,FILE="Trajectory_values.dat",FORM="FORMATTED")
OPEN (2,FILE="Transient_rad_plus_cond.dat",FORM="FORMATTED")
OPEN (3,FILE="TempProf.dat",FORM="FORMATTED")
READ(99,*) 

PI      = 3.14125
rho     = 2500.0
Cp      = 850.0
T1      = 300.0
T2      = 300.0
dX      = 0.01
ddX     = dX/100
eps     = 0.5
r1      = 0.1
r2      = 0.1-0.01
sig     = 5.67E-8
k       = 0.084
Ti      = 300.0
Tn(:)   = 300.0
d       = 4.94E-8/ddX**2 
error   = 1
dtC     = 0.025/d

! Read this outside the loop transient !
DO 1 I = 1,130501
READ(99,*) Q,t,X,Y
!DO 1 I = 1,130501
!______________________________________!
!
!               FTCS
!______________________________________!
DO 369 J = 1,2
75  Qc      = Q-eps*sig*Tn(2)**4
    Tn(1)   = Tn(3)+2*ddX*Qc/k
    Tn(2)   = Tn(2)+d*dtC*(Tn(1)-2*Tn(2)+Tn(3))
    error   = (Qc-Q+eps*sig*Tn(2)**4)/Qc
!    WRITE(*,*) "error",error,Tn(1),Tn(2),Tn(3),dtC
IF(error .ge. 0.5) GO TO 75
! Inner thermal heat conduction loop !
DO 55 L = 3,99

Tn(L) = Tn(L)+d*dtC*(Tn(L-1)-2*Tn(L)+Tn(L+1))

55 CONTINUE
369  CONTINUE
IF ( I .EQ. 120000) THEN
DO 58 L = 3,100
WRITE(3,*) Tn(L)
58 CONTINUE
END IF
WRITE(2,*) t,Tn(2)
1 CONTINUE

END PROGRAM
