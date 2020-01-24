module mod_direct_aux
implicit none

!--internal 
double precision::DELTA
double precision::PT
double precision::Q2FAC
double precision::Q2FRAG
double precision,parameter::NC=3d0
double precision,parameter::CF=4d0/3d0
double precision,parameter::PI=3.14159265359d0
double precision::NF;
contains

SUBROUTINE SETUP(PT_,Q2FAC_,Q2FRAG_,NF_,DELTA_)
double precision::DELTA_,PT_,Q2FAC_,Q2FRAG_,NF_
DELTA=DELTA_
PT=PT_
Q2FAC=Q2FAC_
Q2FRAG=Q2FRAG_
NF=NF_
END SUBROUTINE
SUBROUTINE SIGLO(SIGQQB,SIGQG,SIGGQ,V)
double precision::SIGGQ
double precision::SIGQG
double precision::SIGQQB
double precision::V
SIGQQB=2.D0*CF/NC*PI*(V**2+(1.D0-V)**2)
SIGQG =1.D0/NC*PI*V*(1.D0+(1.D0-V)**2)
SIGGQ =1.D0/NC*PI*(1.D0-V)*(1.D0+V**2)
!print*,'SIGQQB=',SIGQQB
!print*,'SIGQG=',SIGQG 
!print*,'SIGGQ=',SIGGQ 
RETURN
END SUBROUTINE
SUBROUTINE DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMSS)
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::SHAT
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
SHAT=PT**2/V/W/(1.D0-V)
V1=1.D0-V
V2=2.D0-V
V3=3.D0-V
V4=4.D0-V
V5=5.D0-V
X=1.D0-V*W
Y=1.D0-V+V*W
LV=DLOG(V)
L1V=DLOG(1.D0-V)
LW=DLOG(W)
IF(W.NE.1.) THEN
	L1W=DLOG(1.D0-W)
ELSE
	L1W=0.D0
ENDIF
LLVW=DLOG(1.D0-V*W)
L1VVW=DLOG(1.D0-V+V*W)
LMU=0.D0
LMS=DLOG(Q2FAC/SHAT)
LMSS=DLOG(Q2FRAG/SHAT)
RETURN
END SUBROUTINE
SUBROUTINE DEFIN1(V,W,V1,Y,BFAC)
double precision::BFAC
double precision::DELTA
double precision::LTOT
double precision::V
double precision::V1
double precision::W
double precision::Y
V1=1.D0-V
Y=1.D0-V+V*W
LTOT=DLOG(V**2*(1.D0-W)**2*PT**2*DELTA**2/Q2FRAG)
BFAC=Y**2+LTOT*(1.D0+V**2*(1.D0-W)**2)
!print*,'LTOT=',LTOT
!print*,'BFAC=',BFAC
RETURN
END SUBROUTINE
FUNCTION FQQB1(V)
double precision::B0
double precision::C1
double precision::CA
double precision::FQQB1
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQQB
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::X
double precision::Y
CALL DEFIN(V,1.D0,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,&
&LMSS)
B0=11.D0/6.D0*NC-NF/3.D0
TQQB=1.D0-2.D0*V+2.D0*V**2
C1= - CF/NC*NF*TQQB*(5.D0 - 3.D0*LV)/9.D0 +CF**2/NC*(-(7.D0- PI**2)*TQ&
&QB - 2.D0*TQQB*L1V*LV +L1V*V*(2.D0+ V) + LV*V1*(2.D0+ V1)+LV**2*(3.D0*&
&V**2 + 2.D0*V1)+L1V**2*(2.D0*V + 3.D0*V1**2)) +CF*((67.D0 - 6.D0*PI**2&
&)*TQQB/18.D0 - L1V*V*V1 -LV*(11.D0 - 16.D0*V*V1)/6.D0- L1V**2*(2.D0*V&
&+ 3.D0*V1**2)/2.D0 -LV**2*V*V2/2.D0)
CA=-(CF**2/NC*TQQB*LMS*(3.D0 - 2.D0*L1V + 2.D0*LV))
FQQB1=C1+CA
!print*,'FQQB1=',FQQB1
RETURN
END FUNCTION
FUNCTION FQQB1MU(V)
double precision::B0
double precision::FQQB1MU
double precision::TQQB
double precision::V
B0=11.D0/6.D0*NC-NF/3.D0
TQQB=1.D0-2.D0*V+2.D0*V**2
FQQB1MU=CF/NC*TQQB*B0
!print*,'FQQB1MU=',FQQB1MU
RETURN
END FUNCTION
FUNCTION FQQB2(V,W)
double precision::C2
double precision::CB
double precision::FQQB2
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQQB
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMSS)
TQQB=1.D0-2.D0*V+2.D0*V**2
C2=TQQB*(CF/NC*NF/3.D0 - CF*(11.D0 - 12.D0*L1V)/6.D0 -4.D0*CF**2/NC*(L1V - LV))
CB=-4.D0*CF**2/NC*LMS*(1.D0 - 2.D0*V*V1)
FQQB2=C2+CB
!print*,'FQQB2=',FQQB2
RETURN
END FUNCTION
FUNCTION FQQB3(V,W)
double precision::C3
double precision::FQQB3
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQQB
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMSS)
TQQB=1.D0-2.D0*V+2.D0*V**2
C3=2.D0*CF/NC*(4.D0*CF - NC)*TQQB
FQQB3=C3
!print*,'FQQB3=',FQQB3
RETURN
END FUNCTION
FUNCTION FQG1(V)
double precision::B0
double precision::C1
double precision::CA
double precision::FQG1
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQG
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::X
double precision::Y
CALL DEFIN(V,1.D0,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,&
&LMSS)
B0=11.D0/6.D0*NC-NF/3.D0
TQG=2.D0-2.D0*V+V**2
C1= - LV**2*V*(V**2 - 2.D0*V1)/4.D0+L1V*V*V1/2.D0-LV*V*V1 + PI**2*V*(1&
&.D0 + V)*V1/4.D0+ L1V**2*V*(1.D0 + V)*V1/4.D0 -L1V*LV*V*(1.D0 + V)*V1/&
&2.D0 +CF/NC*(-7.D0*TQG*V/4.D0 + L1V*V*(1.D0 + 2.D0*V)/2.D0 +PI**2*V*(1&
&.D0 - 4.D0*V + 5.D0*V**2)/6.D0-LV*V*(3.D0*V**2 - 2.D0*V1)/4.D0+LV**2*V&
&*(3.D0*V**2 + 2.D0*V1)/2.D0+ L1V**2*V*(V**2 + V1**2)/2.D0 -L1V*LV*V*(V&
&**2 + V1**2))
CA=TQG*LMS*(-3.D0*CF/NC/4.D0+NF/(6.D0*NC)-(11.D0-12.D0*L1V+12.D0*LV)/1&
&2.D0)*V
FQG1=C1+CA
!print*,'FQG1=',FQG1
RETURN
END FUNCTION
FUNCTION FQG1MU(V)
double precision::B0
double precision::FQG1MU
double precision::TQG
double precision::V
B0=11.D0/6.D0*NC-NF/3.D0
TQG=2.D0-2.D0*V+V**2
FQG1MU=TQG*B0*V/(2.D0*NC)
!print*,'FQG1MU=',FQG1MU
RETURN
END FUNCTION
FUNCTION FQG2(V,W)
double precision::C2
double precision::CB
double precision::FQG2
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQG
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQG=2.D0-2.D0*V+V**2
C2=TQG*(-L1V - CF/NC*(3.D0 - 4.D0*LV)/4.D0 + LV)*V
CB=-(CF + NC)*TQG*LMS*V/NC
FQG2=C2+CB
!print*,'FQG2=',FQG2
RETURN
END FUNCTION
FUNCTION FQG3(V,W)
double precision::C3
double precision::FQG3
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQG
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQG=2.D0-2.D0*V+V**2
C3=(CF + 2.D0*NC)*TQG*V/NC
FQG3=C3
!print*,'FQG3=',FQG3
RETURN
END FUNCTION
FUNCTION FGQ1(V)
double precision::B0
double precision::C1
double precision::CA
double precision::FGQ1
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQGC
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::X
double precision::Y
CALL DEFIN(V,1.D0,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,&
&LMSS)
B0=11.D0/6.D0*NC-NF/3.D0
TQGC=1.D0+V**2
C1=( - (L1V*V*V1) + LV*V*V1/2.D0- L1V**2*(1 - 4.D0*V + V**2)*V1/4.D0 +&
&PI**2*V*V2*V1/4.D0 - L1V*LV*V*V2*V1/2.D0+ LV**2*V*V2*V1/4.D0) +CF/NC*(&
&-7.D0*TQGC*V1/4.D0 + 2.D0*L1V*V*V1+ LV*(3.D0 - 4.D0*V - 3.D0*V**2)*V1/&
&4.D0 +PI**2*(2.D0 - 6.D0*V + 5.D0*V**2)*V1/6.D0 + L1V**2*V1**3 +LV**2*&
&V1*(3.D0*V**2 + 2.D0*V1)/2.D0- L1V*LV*V1*(V**2 + V1**2))
CA=TQGC*LMS*(-11.D0/12.D0 + NF/(6.D0*NC)- CF/NC*(3.D0 - 4.D0*L1V + 4.D&
&0*LV)/4.D0)*V1
FGQ1=C1+CA
!print*,'FGQ1=',FGQ1
RETURN
END FUNCTION
FUNCTION FGQ1MU(V)
double precision::B0
double precision::FGQ1MU
double precision::TQGC
double precision::V
B0=11.D0/6.D0*NC-NF/3.D0
TQGC=1.D0+V**2
FGQ1MU=TQGC*B0*(1.D0-V)/(2.D0*NC)
!print*,'FGQ1MU=',FGQ1MU
RETURN
END FUNCTION
FUNCTION FGQ2(V,W)
double precision::C2
double precision::CB
double precision::FGQ2
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQGC
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQGC=1.D0+V**2
C2=TQGC*(-L1V - CF/NC*(3.D0 - 4.D0*LV)/4.D0 + LV)*V1
CB=-(CF + NC)*TQGC*LMS*V1/NC
FGQ2=C2+CB
!print*,'FGQ2=',FGQ2
RETURN
END FUNCTION
FUNCTION FGQ3(V,W)
double precision::C3
double precision::FGQ3
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQGC
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQGC=1.D0+V**2
C3=(CF + 2.D0*NC)*TQGC*V1/NC
FGQ3=C3
!print*,'FGQ3=',FGQ3
RETURN
END FUNCTION
FUNCTION FQQB(V,W)
double precision::C10
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C7
double precision::C8
double precision::C9
double precision::CC
double precision::FQQB
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQQB
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQQB=1.D0-2.D0*V+2.D0*V**2
CC=CF**2/NC*LMS*((1.D0 - 2.D0*V)**2 - (1.D0 - 2.D0*V)*V/X +V*V1/X**2 +&
&(3.D0*V**2 + V1**2)*W)
C4=-CC*LV/LMS
C5=CF**2/NC*L1V*(-2.D0*(2.D0 - X) + 2.D0*(1.D0 + V1**2)/X) +CF*L1V*((1&
&.D0 + V)*V1 - (1.D0 + V1**2)/X + V*V1*W)
C6=CF*LW*(1.D0 - V**2 - (1.D0 + V1**2)/X + V*V1*W) +2.D0*CF**2/NC*LW*(&
&V1**2 + V**2*W**2)/X
C7=CF*TQQB*LW/(1.D0 - W)
C8=CF*L1W*(-1.D0 - 2.D0*V*V1 + (1.D0 + V1**2)/X- (1.D0 - 2.D0*V)*V*W)&
&+CF**2/NC*L1W*(1.D0 + 8.D0*V*V1 - V*V1/X**2- (3.D0*v + 4.D0*v1**2)/X -&
&(1.D0 - 4.D0*V + 8.D0*V**2)*W)
C9=-4.D0*CF**2/NC*LLVW*V*(V1 - V*W) + CF*LLVW*V*(V2 - V*W)
C10=CF/NC*(4.D0*CF - NC)*TQQB*(L1V-LLVW)/(1.D0-W)
C13=CF*NF*V*(V2 - V*W)/(3.D0*NC) + CF**2/NC*(-4.D0*V*V1 + V*V1/X**2 +&
&V*V1/X + (1.D0 + V - 4.D0*V**2)*W) +CF*(-(1.D0 - 4.D0*V)*V/(2*X) + V*V&
&1/(2.D0*X**2)- (2.D0 - 11.D0*V1**2)/6.D0+(3.D0 - 12.D0*V + 11.D0*V**2)&
&*W/6.D0)
FQQB=CC+C4+C5+C6+C7+C8+C9+C10+C13

!print*,'TQQB=',TQQB
!print*,'LW=',LW
!print*,'W=',W

!print*,'CC=',CC
!print*,'C4=',C4
!print*,'C5=',C5
!print*,'C6=',C6
!print*,'C7=',C7
!print*,'C8=',C8
!print*,'C9=',C9
!print*,'C10=',C10
!print*,'C13=',C13
!print*,'FQQB=',FQQB
RETURN
END FUNCTION
FUNCTION FQG(V,W)
double precision::C10
double precision::C11
double precision::C12
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C7
double precision::C8
double precision::C9
double precision::CC
double precision::CD
double precision::FQG
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQG
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQG=2.D0-2.D0*V+V**2
CC=CF/NC*LMS*(-(V*V1**2/X**3) + V*V1*(1.D0 + 2.D0*V1)/X**2 -V*(1.D0 +&
&2.D0*V1)**2/(2.D0*X) +V*(V**2 + 6.D0*V1)/2.D0 + V*(2.D0+ V**2)*W/2.D0)&
&+LMS*(-V + V*(4.D0- 3.D0*V + V**2)/X+ V*V1**2/X**3-V*V1*V2/X**2-V**2*&
&(2.D0*V**2 + 3.D0*V1)*W/V1 + V**3.D0*(1.D0 + V)*W**2/V1- V**4*W**3/V1)
CD=-LMSS*V*(1.D0 + V**2*(1.D0 - W)**2)*(CF/NC*V1**2 + Y*V*W)*(V1**2 +&
&2.D0*Y*V*W)/(2.D0*Y**3*V1)
C4=CF/NC*LV*V*((3.D0*V**2 - 10.D0*V1)/2.D0+ (3.D0- V)**2*V1/(2.D0*Y)+&
&V1**2/X**3 -(3.D0 - V)*V1**2/Y**2 + V1**3/Y**3- V1*(1.D0 + 2.D0*V1)/X*&
&*2 +(1.D0 + 2.D0*V1)**2/(2.D0*X) -(11.D0*V**2 + 2.D0*V1)*W/2.D0+4.D0*V&
&**2*W**2)+LV*V*(-(4.D0 - 3.D0*V + V**2)*V1 + (4.D0 - V)*V1**2/Y- V1**3&
&/X**3 -V1**3/Y**2 - V1*V2/X + V1**2*V2/X**2 +V*(7.D0+ V**2)*W/2.D0 -4.&
&D0*V**2*W**2 + 2.D0*V**3*W**3)/V1
C5=-CF/NC*L1V*V**3*(1.D0 + W) - L1V*V*((1.D0 + V)*V1**2 - 4.D0*V*W*(2.&
&D0*X - V*V1+ Y*V*W))/(2.D0*V1)
C6=-CF/NC*LW*V**3*(1.D0 + 4.D0*W - 2.D0*W**2) +LW*V*(2.D0*TQG/X + V**2&
&- 2.D0*V1- V*W*(3.D0 - 4.D0*V + V*W))/2.D0
C7=(2.D0*CF - NC)*LW*V*(3.D0*V**2+ 2.D0*V1)/(2.D0*NC*(1.D0 - W))
C8=C4*L1W/LV-CF/NC*L1W*V**3*(1.D0 - 2.D0*W + 2.D0*W**2) -L1W*V*(TQG/X&
&- (4.D0 - 3.D0*V)*V*W/2.D0 - V**2*W**2/2.D0)
C9=CF/NC*LLVW*V**3*(3.D0 - 2.D0*W)*W + LLVW*V*((5.D0 - 4.D0*V + V**2)*&
&V1 - V*(8.D0 - 3.D0*V*V1)*W+ V**2*W**2*(7.D0 + V - 4.D0*V*W))/(2.D0*V1&
&)
C10=(L1V-LLVW)/(1.D0-W)*V*(5.D0 - 4.D0*V + V**2)/2.D0 +CF/NC*(L1V-LLVW&
&)/(1.D0-W)*V*(V**2 + V1**2)
C11=CF/NC*L1VVW*V*(-4.D0*V1 + (3.D0 - V)**2*V1/Y -2.D0*(3.D0 - V)*V1**&
&2/Y**2 +2.D0*V1**3/Y**3 + V*(2.D0*X + V)*W) +L1VVW*V*(-11.D0/2.D0 + V&
&+ 2.D0*(4.D0 - V)*V1/Y -2.D0*V1**2/Y**2 + 3.D0*V*W/2.D0 +  V**2*W**2/2&
&.D0)
C12=-L1VVW/(1.D0-W)*(1.D0 - 2.D0*V)*V/2.D0 -CF/NC*L1VVW/(1.D0-W)*V*(1.&
&D0 + V**2)
C13=V*W*(-TQG*V1/2.D0 - 2.D0*V*V1**2/Y**2 -(9.D0 - 4.D0*V)*V*V1**2/(2.&
&D0*X**2) +3.D0*V*V1**3/X**3 + V*V1*(1.D0 + 3.D0*V1)/(2.D0*X) +V*V1*V2/&
&Y + Y*V**2.D0*W)/V1+ CF/NC*(-V*(22.D0 - 35.D0*V + 11.D0*V**2)/(2.D0*X)&
&+(23.D0 - 16.D0*V)*V*V1/(2.D0*X**2)- V*(13.D0 - 3.D0*V - 3.D0*V**2)*V&
&1/(2.D0*Y) -4.D0*V*V1**2/X**3 + (17.D0 - 5.D0*V)*V*V1**2/(2.D0*Y**2)-&
&3.D0*V*V1**3/Y**3 +3.D0*V*(V**2 + 10.D0*V1)/4.D0 + V*(4.D0 - 7.D0*V**2&
&)*W/4.D0)
FQG=CC+CD+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13
!print*,'FQG=',FQG
RETURN
END FUNCTION
FUNCTION FQGC(V,W)
double precision::BFAC
double precision::DD
double precision::FQGC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=V*BFAC*(CF*V1**2 + NC*Y*V*W)*(Y**2 + V**2*W**2)/(2.D0*NC*Y**3*V1)
FQGC=DD
!print*,'FQGC=',FQGC
RETURN
END FUNCTION
FUNCTION FGQ(V,W)  
double precision::FGQ
double precision::HELP1
double precision::HELP2
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::TQGC
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
TQGC=1.D0+V**2
HELP1=TQGC*(L1V-LLVW)/(1.D0-W)*V1 - TQGC*LW*V1/(1.D0 - W) +(-3.D0*CF +&
&4.D0*CF*L1W + 8.D0*NC*L1W - 4.D0*NC*LLVW+ 4.D0*CF*LV + 4.D0*NC*LV -4.&
&D0*NC*LW)*V*(1.D0 - V + V**2 - V*W + V**2*W+ V**2*W**2)/(4.D0*NC)
HELP2= - (CF + NC)*LMS*V*(1.D0 - V + V**2 - V*W+ V**2*W + V**2*W**2)/N&
&C
FGQ=HELP1+HELP2+V/(1.D0-V*W)*FQG(1.D0-V*W,(1.D0-V)/(1.D0-V*W))
!print*,'FGQ=',FGQ
RETURN
END FUNCTION
FUNCTION FGQC(V,W)
double precision::FGQC
double precision::V
double precision::W
FGQC=V/(1.D0-V*W)*FQGC(1.D0-V*W,(1.D0-V)/(1.D0-V*W))
!print*,'FGQC=',FGQC
RETURN
END FUNCTION
FUNCTION FGG(V,W)
double precision::C11
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C8
double precision::C9
double precision::CC
double precision::CD
double precision::FGG
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
CC=1.D0/NC*LMS*V*(-(7.D0 - 8.D0*V + 3.D0*V**2)/2.D0- 2.D0*V1**2/X**2+&
&2.D0*V1*V2/X +(2.D0 - 3.D0*V + 2.D0*V**2)*W- (3.D0*V**2 + 4.D0*V1)*W**&
&2/2.D0)
CD=-LMSS*V*(X**2 + V2**2 - 4.D0*V2*V1/Y + 4.D0*V1**2/Y**2)*(1.D0/NC -&
&(1.D0 - V)*V*W/(CF*Y**2))/2.D0
C4=1.D0/CF*LV*V**2*V1*(3.D0 - 4.D0*V1**2/Y**4 + 4.D0*V1/(Y*V2)- 2.D0*(&
&1.D0 + V1**2)/(X*V2) +4.D0*V1*V2/Y**3 - 2.D0*V2**2/Y**2)*W/2.D0 +2.D0*&
&1.D0/NC*LV*V*(V1**2/X**2 + V1**2/Y**2- 2.D0*V1**2/(X*V2) -2.D0*V1**2/(&
&Y*V2) + (V**2 + V1)*(1.D0 - W + W**2))
C5=1.D0/CF*L1V*V**2*V1*W*(1.D0 - 2.D0*W+ 2.D0*(1.D0 + V1**2)*W/(X*Y))/&
&2.D0+1.D0/NC*L1V*V*((4.D0 - 3.D0*V + V**3)/V- 2.D0*V1*(1.D0 + V1**2)/(&
&X*Y*V)- 2.D0*V*(1.D0 - W)*W)
C6=1.D0/CF*LW*V*V1*(1.D0 - 3.D0*V1+ 2.D0*V1*(1.D0 + V1**2)/(X*Y))*W/2.&
&D0 +1.D0/NC*LW*V*(-2.D0*V1*(1.D0 - W)+ 2.D0*V1*(1.D0 + V1**2)*(1.D0 -&
&W)/(X*Y) +  V**2*W**2)
C8=1.D0/CF*L1W*V*V1*(1.D0 - 2.D0*V*V1**2/Y**4- (1.D0 + V1**2)/(X*V2)+&
&2.D0*V*V1*V2/Y**3 -V*V2**2/Y**2 - V1*(2.D0*V1 - V*V2)/(Y*V2))*W +1.D0/&
&NC*L1W*V*(2.D0*V1**2/X**2 + 2.D0*V1**2/Y**2+ 3.D0*(1.D0 + V1**2) -2.D0&
&*(2.D0 - V**2)*V1/(Y*V*V2) +2.D0*V1**2*(2.D0*V1 - V*V2)/(X*V*V2) -4.D0&
&*(V**2 + V1)*W + (3.D0*V**2 + 2.D0*V1)*W**2)
C9=-1.D0/NC*LLVW*V*(V**2 + V1**2 + V*W*(2.D0*V1 + V*W))
C11=1.D0/CF*L1VVW*V**2*V1*(-3.D0 - 4.D0*V1**2/Y**4 + 2.D0*V2/Y+ 4.D0*V&
&1*V2/Y**3 -    2.D0*V2**2/Y**2)*W- 4.D0*1.D0/NC*L1VVW*V**3*V1*(1.D0 -&
&W)*W/Y**2
C13=1.D0/CF*V*V1*W*(-2.D0*V/Y + 2.D0*(7.D0 - 6.D0*V)*V/Y**2- 4.D0*V*V1&
&/X**2+ 12.D0*V*V1**2/Y**4 +V2 + 2.D0*V*V2/X - 12.D0*V*V1*V2/Y**3 - 2.D&
&0*V2*W)/2.D0 +1.D0/NC*V*(-(3.D0 - V)*V1 - 4.D0*V1**2/X**2 + 4.D0*V1*V2&
&/X +2.D0*(2.D0 - 3.D0*V + 2.D0*V**2)*W - V2**2*W**2)/2.D0
FGG=CC+CD+C4+C5+C6+C8+C9+C11+C13
!print*,'FGG=',FGG
RETURN
END FUNCTION
FUNCTION FGGC(V,W)
double precision::BFAC
double precision::DD
double precision::FGGC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=V*BFAC*(CF*Y**2 - NC*V*V1*W)*(V1**2 + V**2*W**2)/(2.D0*CF*NC*Y**4)
FGGC=DD
!print*,'FGGC=',FGGC
RETURN
END FUNCTION
FUNCTION FQQ(V,W)
double precision::C11
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C8
double precision::C9
double precision::CC
double precision::CD
double precision::FQQ
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
CC=CF/NC*LMS*(-(V1**2/X**2) - V1*(5.D0 - 3.D0*V + 2.D0*V**2*V1)+ 2.D0*&
&V1*V2/X -2.D0*V1*(V**3 - V1**2)*W-(4.D0*V**2-2.D0*V**3*V1+V1**2)*W**2&
&+2.D0*(2 - Y)*V**3*W**3)/(2.D0*V1*W)
CD=LMSS*(1.D0 + V**2*(1.D0 - W)**2)*(CF/NC**2*Y**2*V*V1*W -CF/NC*(V1**&
&4 + V*V1**3*W + V**2*V1**2*W**2 + Y*V**3*W**3))/(Y**2*(1.D0 - V)*W)
C4=CF/NC**2*LV*V*(-(1.D0 + V)**2 + 4.D0*V1/(Y*V2)- 2.D0*(1.D0 + V1**2)&
&/(X*V2) +2.D0*V*(1.D0 + V)*W - V**2*W**2) +CF/NC*LV*((43.D0 - 33.D0*V&
&+ 12.D0*V**2 - 4.D0*V**3)*V1+ V1**2/X**2 -8.D0*(5.D0 - 4.D0*V + V**2)*&
&V1**2/Y + 4.D0*V1**4/Y**2- 2.D0*V1*V2/X -2.D0*V1*(1.D0 + 11.D0*V + V**&
&2*(1.D0 + 2.D0*V1))*W +(22.D0*V**2 - 18.D0*V**3 + 4.D0*V**4 + V1**2)*W&
&**2- 4.D0*(2.D0 - Y)*V**3*W**3)/(2.D0*V1*W)
C5=CF/NC**2*L1V*(V1 - 4.D0*V1/(Y*V2) + V*(1.D0 + V1**2)/(X*V2)- V**3*(&
&1.D0 - W)**2 -V*W) + CF/NC*L1V*V*(2.D0*(1.D0 - 2.D0/Y)*V1**2+ V*(3.D0&
&+ V**2)*W -2.D0*(2.D0 - Y)*V**2*W**2)/V1
C6=CF/NC*LW*(2.D0*(3.D0 + V**2)*V1 - 4.D0*V1**2/Y -2.D0*V*(1.D0 + V)*V&
&2*W +(3.D0 - V)*V**2*W**2)/W +CF/NC**2*LW*(-(1.D0 + V + V**3) + 4.D0*V&
&1/(Y*V2) -V*(1.D0 + V1**2)/(X*V2) +V*(1.D0 + 2.D0*V**2)*W - V**3*W**2)
C8=C4*L1W/LV-CF/NC*L1W*(X - V)*V*(X + V)/Y -CF/NC**2*L1W*(1.D0 + V**2*&
&(1.D0 - W)**2)*(V1**2- V*(1.D0 + V**2*(1.D0 - W))*W)/(X*Y)
C9=CF/NC*LLVW*V*((1.D0 + V)*V1 - V*(4.D0 - V*V1)*W+ 2.D0*(2.D0 - Y)*V*&
&*2*W**2)/V1
C11=4.D0*CF/NC*L1VVW*V**3*V1*(1.D0 - W)*W/Y**2
C13=-CF/NC**2*Y**2*V + CF/NC*(-(27.D0-30.D0*V+2.D0*V**2+ 2.D0*V**3)*V1&
&- V1**2/X**2 +8.D0*(5.D0 - v)*v1**3/Y-16.D0*v1**4/Y**2+V1*(1.D0 + 5.D0&
&*V1)/X +V1*(1.D0+2.D0*(3.D0-V)*V*V1)*W-V1*(1.D0 - 4.D0*V*(1.D0 + V)*V1&
&)*W**2 +2.D0*Y*V**3*W**3)/(2.D0*V1*W)
FQQ=CC+CD+C4+C5+C6+C8+C9+C11+C13
!print*,'FQQ=',FQQ
RETURN
END FUNCTION
FUNCTION FQQC(V,W)
double precision::BFAC
double precision::DD
double precision::FQQC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=CF*BFAC*( - (Y**2*V*V1*W) +NC*(Y**4 - 3.D0*Y**2*V*V1*W + V**2*V1**2&
&*W**2))/(NC**2*Y**2*V1*W)
FQQC=DD
!print*,'FQQC=',FQQC
RETURN
END FUNCTION
FUNCTION FQB(V,W)
double precision::C10
double precision::C11
double precision::C12
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C7
double precision::C8
double precision::C9
double precision::CC
double precision::CD
double precision::FQB
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
CC=CF/NC*LMS*(-V1**2/X**2 - V1*(5.D0 - 3.D0*V + 2.D0*V**2*V1)+ 2.D0*V1&
&*V2/X -2.D0*V1*(V**3 - V1**2)*W - (4.D0*V**2 - 2.D0*V**3*V1+ V1**2)*W*&
&*2 +2.D0*(2.D0 - Y)*V**3*W**3)/(2.D0*V1*W)
CD=LMSS*(1.D0 + V**2*(1.D0 - W)**2)*(-CF/NC**2*Y**2*V*V1*W*(V1**2 - V*&
&V1*W + V**2*W**2) -CF/NC*(V1**6 + 3.D0*V*V1**5*W + 5.D0*V**2*V1**4*W**&
&2 +4.D0*V**3*V1**3*W**3 + 5.D0*V**4*V1**2*W**4 +3.D0*V**5*V1*W**5 + V*&
&*6*W**6))/(Y**4*V1*W)
C4=CF/NC**2*LV*V*(-2.D0 - 2.D0*(1.D0 + V - V**2)*V1/Y - V1**2+ 6.D0*V1&
&**2/Y**2 + V*(5.D0 + V)*W - 3.D0*V**2*W**2) +CF/NC*LV*((11.D0 - 53.D0*&
&V + 40.D0*V**2 - 12.D0*V**3)*V1+ V1**2/X**2 -4.D0*(10.D0 - 27.D0*V + 1&
&6.D0*V**2 - 3.D0*V**3)*V1**2/Y +4.D0*V3*V5*V1**4/Y**2- 8.D0*V4*V1**5/Y&
&**3+8.D0*V1**6/Y**4 - 2.D0*V1*V2/X - 2.D0*(1.D0 - 5.D0*V- 7.D0*V**2)*V&
&1*W +(2.D0*V**3 + 6.D0*V**4 + V1**2)*W**2- 4.D0*(2.D0 - Y)*V**3*W**3)/&
&(2.D0*V1*W)
C5=CF/NC**2*L1V*V*(4.D0 - 3.D0*V- 2.D0*V1*(1.D0 + V**2 + 2.D0*V1)/Y- V&
&*(X + 2.D0*V1)*W) +CF/NC*L1V*V*(-2.D0*V1**2 + 4.D0*V1**2/Y + V*(3.D0 +&
&V**2)*W -2.D0*(2.D0 - Y)*V**2*W**2)/V1
C6=CF/NC**2*LW*V*(-(5.D0 - 4.D0*V + V**2)+2.D0*V1*(1.D0+V**2 + 2.D0*V1&
&)/Y+ 2.D0*(1.D0 + X)*V*W) +CF/NC*LW*(4.D0*V1**2/Y - 2.D0*(1.D0 + V)*V1&
&**2+2.D0*V*(V**2 + V2)*W -V**2*(1.D0 + V)*W**2)/W
C7=-CF/NC**2*LW*(1.D0 + V**2)/(1.D0 - W)
C8=C4*L1W/LV+L1W*(CF/NC*(2.D0 - Y)*V*(Y - 2.D0*V*W)/Y -CF/NC**2*V**2*(&
&Y - 2.D0*V*W)*(2.D0*V - W - 3.D0*V*W+ 2.D0*V*W**2)/Y)
C9=CF/NC**2*LLVW*V*(V2 - (3.D0 - Y)*V*W) -CF/NC*LLVW*V*((1.D0+ V)*V1 +&
&2.D0*X*V*W+ V**2*W*(1.D0 + X - Y + 2.D0*V*W**2))/V1
C10=-CF/NC**2*(L1V-LLVW)/(1.D0-W)*(1.D0 + V1**2)
C11=2.D0*CF/NC*L1VVW*V*(-2.D0*(4.D0 - 3.D0*V + V**2)- 2.D0*V3**2*V1**2&
&/Y**2 +4.D0*V3*V1**3/Y**3 - 4.D0*V1**4/Y**4+ 4.D0*V1*V2**2/Y + V3*V*W)&
&+ CF/NC**2*L1VVW*V*(15.D0 - 13.D0*V + 3.D0*V**2+ 12.D0*V1**2/Y**2 -12.&
&D0*V1*V2/Y - 7.D0*V*W + 3.D0*V**2*W**2)
C12=CF/NC**2*L1VVW/(1.D0-W)*(3.D0 - 2.D0*V*V1)
C13=CF/NC**2*(-3.D0 - (5.D0 - 2.D0*V)*V/X + V*V1/X**2+ 4.D0*V3*V*V1/Y&
&+ 2.D0*V*V1**2 -8.D0*V*V1**2/Y**2 + (V**2 + 2.D0*V**3 + V1**2)*W+ 2.D0&
&*V**3*W**2)/2.D0 +CF/NC*(-(31.D0 - 54.D0*V + 22.D0*V**2+ 2.D0*V**3)*V1&
&- V1**2/X**2 +4.D0*(27.D0 - 19.D0*V)*V1**3/Y- 4.D0*(37.D0 - 18.D0*V)*&
&V1**4/Y**2 +24.D0*V4*V1**5/Y**3 - 24.D0*V1**6/Y**4 + V1*V2/X +V1*(1.D0&
&+ 2.D0*V3*V*V1)*W-V1*(1.D0 - 4.D0*V**2*V1)*W**2 +2.D0*Y*V**3*W**3)/(2&
&.D0*V1*W)
FQB=CC+CD+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13
!print*,'FQB=',FQB
RETURN
END FUNCTION
FUNCTION FQBC(V,W)
double precision::BFAC
double precision::DD
double precision::FQBC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=CF*BFAC*(2.D0*NC*V**2*V1**2*W**2*(Y**2 - V*V1*W) +Y**2*(Y**2 - 3.D0&
&*V*V1*W)*(NC*Y**2 + V*V1*W))/(NC**2*Y**4*V1*W)
FQBC=DD
!print*,'FQBC=',FQBC
RETURN
END FUNCTION
FUNCTION FQBS(V,W)
double precision::C11
double precision::C13
double precision::C4
double precision::C8
double precision::CCC
double precision::CD
double precision::FQBS
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
CCC=CF/NC*V**2*V1*(-1.D0 - 4.D0*V1**2/Y**4+ 2.D0*V2/Y + 4.D0*V1*V2/Y**&
&3 - 2.D0*V2**2/Y**2)*W
CD=CCC*LMSS
C4=-CCC*LV
C8=-CCC*L1W
C11=-2.D0*CCC*L1VVW
C13=2.D0*CF/NC*V**3*V1*(1.D0 - W)*W*(-1.D0 + 6.D0*V*V1*W/Y**2)/Y**2
FQBS=CD+C4+C8+C11+C13
!print*,'FQBS=',FQBS
RETURN
END FUNCTION
FUNCTION FQBSC(V,W)
double precision::BFAC
double precision::DD
double precision::FQBSC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=CF*V**2*V1*BFAC*W*(V1**2 + V**2*W**2)/(NC*Y**4)
FQBSC=DD
!print*,'FQBSC=',FQBSC
RETURN
END FUNCTION
FUNCTION FQS(V,W,CQ,CQS)
double precision::C11
double precision::C13
double precision::C4
double precision::C5
double precision::C6
double precision::C8
double precision::C9
double precision::CC
double precision::CD
double precision::CQ
double precision::CQS
double precision::FQS
double precision::L1V
double precision::L1VVW
double precision::L1W
double precision::LLVW
double precision::LMS
double precision::LMSS
double precision::LMU
double precision::LV
double precision::LW
double precision::V
double precision::V1
double precision::V2
double precision::V3
double precision::V4
double precision::V5
double precision::W
double precision::X
double precision::Y
CALL DEFIN(V,W,V1,V2,V3,V4,V5,X,Y,LV,L1V,LW,L1W,LLVW,L1VVW,LMU,LMS,LMS&
&S)
CC=CF/NC*CQ**2*LMS*V**2*W*(-(3.D0 + V**2)/2.D0 + V1/X- V1**2/(2.D0*X**&
&2) + (2.D0 - Y)*V*W)/V1- CF/NC*CQS**2*LMS*(1.D0 + V**2)*V1*(2.D0- 2.D0&
&*W + W**2)/(2.D0*W)
CD=CF/NC*LMSS*(1.D0 + V**2*(1.D0 - W)**2)*(-(CQ**2*V**2*W**2*(V1**2 +&
&2.D0*Y*V*W)) -CQS**2*V1**2*(2.D0*V1**2 + V*W*(2.D0*Y - V*W)))/(2.D0*Y*&
&*2*(1.D0 - V)*W)
C4=2.D0*CF/NC*CQ*CQS*LV*V*(-(5.D0+ V)+ (7.D0 - 4.D0*V + V**2)/Y+ 3.D0*&
&V*W) +CF/NC*CQ**2*LV*V**2*W*(-2.D0*V1/X - 2.D0*(3.D0 - V)*V1/Y+ V1**2/&
&X**2 +2.D0*V1**2/Y**2 + 2.D0*(4.D0 - V*V1)- 4.D0*(2.D0 - Y)*V*W)/(2.D0&
&*V1) +CF/NC*CQS**2*LV*V1*(2.D0*(2.D0*V**2 + V2)+ 2.D0*V*V1**2*(1.D0 -&
&W)/Y**2 -(1.D0 + 2.D0*V**2)*(2.D0 - W)*W)/(2.D0*W)
C5=2.D0*CF/NC*CQ*CQS*(1.D0 - 2.D0/Y)*L1V*V*V1 +CF/NC*CQ**2*L1V*V**2*W*&
&(3.D0 + V**2 - 2.D0*(2.D0 - Y)*V*W)/V1
C6=-2.D0*CF/NC*CQ*CQS*(2.D0 - Y)*LW*V**2*W/Y +CF/NC*CQS**2*LW*V1*(2.D0&
&+ V**2*(2.D0 - 2.D0*W + W**2))/W
C8=C4*L1W/LV-CF/NC*CQ*CQS*(2.D0 - Y)*L1W*V*(V1 - V*W)/Y
C9=CF/NC*CQ*CQS*(2.D0 - Y)*LLVW*V - CF/NC*CQ**2*LLVW*V**2*W*(3.D0 + V*&
&*2 - 2.D0*(2.D0 - Y)*V*W)/V1
C11=2.D0*CF/NC*CQ*CQS*(2.D0 - Y)*L1VVW*V +2.D0*CF/NC*CQS**2*L1VVW*V*V1&
&*(-V/Y - V1/Y**2) +2.D0*CF/NC*CQ**2*L1VVW*V**2*(1.D0 - V3/Y + V1/Y**2)&
&*W
C13=CF/NC*CQ*CQS*(2.D0 - Y)*V*(1.D0/X - 2.D0*V/Y**2)*V1*W +CF/NC*CQ**2&
&*V**2*W*((1.D0 - 2.D0*V)*V1 - V1**2/X**2- 4.D0*V1**2/Y**2 +V1*V2/X + 2&
&.D0*V1*V2/Y + 2.D0*Y*V*W)/(2.D0*V1) +CF/NC*CQS**2*V1*(2.D0*V4*V1/Y - 4&
&.D0*V1**2/Y**2- 2.D0*(1.D0 + V*V1) +(1.D0 + 3.D0*V - 2.D0*V**2)*W -(1.&
&D0 - 2.D0*V)*(1.D0 + V)*W**2)/(2.D0*W)
FQS=CC+CD+C4+C5+C6+C8+C9+C11+C13
!print*,'FQS=',FQS
RETURN
END FUNCTION
FUNCTION FQSC(V,W,CQ,CQS)
double precision::BFAC
double precision::CQ
double precision::CQS
double precision::DD
double precision::FQSC
double precision::V
double precision::V1
double precision::W
double precision::Y
CALL DEFIN1(V,W,V1,Y,BFAC)
DD=CF*BFAC*(CQS**2*V1**2*(Y**2 + V1**2)+ CQ**2*V**2*W**2*(Y**2 + V**2*&
&W**2))/(2.D0*NC*Y**2*V1*W)
FQSC=DD
!print*,'FQSC=',FQSC
RETURN
END FUNCTION
end module

module mod_fragment_aux
implicit none

double precision::DELTA
double precision::PT
double precision::Q2FAC
double precision::Q2FRAG
integer,parameter::JMAR=2 !MSBAR
double precision::Nf

double precision,parameter::NC=3d0
double precision,parameter::CF=4d0/3d0
double precision,parameter::PI=3.14159265359d0

double precision,parameter::N =3d0
double precision,parameter::AL=1d0
double precision,parameter::CQ=1d0

double precision::VC
double precision::V1
double precision::V2
double precision::V3
double precision::V4

double precision::GTR
double precision::S
double precision::GV,GW
integer::J0

contains

SUBROUTINE SET_J0(J0_)
integer::J0_
J0=J0_+1
END SUBROUTINE
SUBROUTINE SETUP(PT_,Q2FAC_,Q2FRAG_,NF_,DELTA_,S_,GV_,GW_)
double precision::DELTA_,PT_,Q2FAC_,Q2FRAG_,NF_,S_,GV_,GW_
DELTA=DELTA_
PT=PT_
Q2FAC=Q2FAC_
Q2FRAG=Q2FRAG_
NF=NF_
GTR=NF/2d0
S=S_
GV=GV_
GW=GW_
VC=NC*NC-1d0
V1=VC*VC/NC
V2=VC/NC
V3=(NC**4-1d0)/2d0/NC**2
V4=VC**2/2d0/NC**2
END SUBROUTINE
SUBROUTINE FBOR(V,SHD,F0)
double precision::F0
double precision::SHD
double precision::V
DIMENSION F0(16)
F0(1)=PI*CF/N/SHD/V/(1.-V)*(V**2+1.)/(1.-V)**2
F0(2)=0.D0
F0(3)=PI*CF/N/SHD/V/(1.-V)*(V**2+1.)/(1.-V)**2
F0(4)=0.D0
F0(5)=PI*CF/N/SHD/V/(1.-V)*(2.*V**2-2.*V+1.)
F0(6)=PI*2.*CF/(N**2)/SHD/V/(1.-V)*(N*V**4-2.*N*V**3+4.*N*V**2+V**2-(3&
&.*N+1.)*V+N)/V**2/(1.-V)**2
F0(7)=0.D0
F0(8)=0.D0
F0(9)=0.D0
F0(10)=0.D0
F0(11)=PI*2.*CF/(N**2)/SHD/V/(1.-V)*(N*V**4-(3.*N+1.)*V**3+(4.*N+1.)*V&
&**2-2.*N*V+N)/(1.-V)**2
F0(12)=PI*CF/(N**2)/SHD/V/(1.-V)*(2.*V**2-2.*V+1.)*(2.*N**2*V**2-2.*N*&
&*2*V+N**2-1.)/V/(1.-V)
F0(13)=PI/(2.*N**2)/SHD/V/(1.-V)*(V**2+1.)*((N**2-1.)*V**2+2.*V+(N**2-&
&1.))/V/(1.-V)**2
F0(14)=PI/(2.*N**2)/SHD/V/(1.-V)*(V**2-2.*V+2.)*((N**2-1.)*V**2-2.*N**&
&2*V+2.*N**2)/V**2/(1.-V)
F0(15)=PI*(4.*N**2)/VC/SHD/V/(1.-V)*(3.-V*(1.-V)+V/(1.-V)**2+(1.-V)/V**2)
F0(16)=PI/(2.*N)/VC/SHD/V/(1.-V)*(V**2+(1.-V)**2)*(2.*N**2*(V**2-V)+N**2-1.)/V/(1.-V)
RETURN
END SUBROUTINE
SUBROUTINE FMU(V,SHD,FDELMU)
double precision::FDELMU
double precision::SHD
double precision::V
integer::J0
DIMENSION FDELMU(16)
DO J0=1,16
FDELMU(J0)=AVDELMU(V,SHD,J0)/(1.D0-V)/SHD
ENDDO
RETURN
END SUBROUTINE
FUNCTION AVDELMU(V,S,J0)
double precision::AVDELMU 
integer::J0
double precision::S
double precision::V
IF (J0.EQ.1) THEN
AVDELMU = 12*(V**2+1)*(4*GTR*V4+4*GTR*V3-11*V2-11*V1)/((V-1)**2*V)/18.&
&D0/(8.D0*N**2)
ELSE IF (J0.EQ.2) THEN
AVDELMU =0.D0
ELSE IF (J0.EQ.3) THEN
AVDELMU = 12*(V**2+1)*(4*GTR*V4+4*GTR*V3-11*V2-11*V1)/((V-1)**2*V)/18.&
&D0/(8.D0*N**2)
ELSE IF (J0.EQ.4) THEN
AVDELMU =0.D0
ELSE IF (J0.EQ.5) THEN
AVDELMU = 12*(2*V**2-2*V+1)*(4*GTR*V4+4*GTR*V3-11*V2-11*V1)/V/18.D0/(8&
&.D0*N**2)
ELSE IF (J0.EQ.6) THEN
AVDELMU = 12*(4*GTR*(V**4-2*V**3+4*V**2-3*V+1)*V4-11*(V-1)*V*V4+4*GTR*&
&(V**4-2*V**3+4*V**2-3*V+1)*V3-11*(V-1)*V*V3-11*(V**4-2*V**3+4*V**2-3*V&
&+1)*V2+4*GTR*(V-1)*V*V2-11*(V**4-2*V**3+4*V**2-3*V+1)*V1)/((V-1)**2*V*&
&*3)/9.D0/(8.D0*N**2)
ELSE IF (J0.EQ.7) THEN
AVDELMU =0.D0
ELSE IF (J0.EQ.8) THEN
AVDELMU = 0.D0
ELSE IF (J0.EQ.9) THEN
AVDELMU = 0.D0
ELSE IF (J0.EQ.10) THEN
AVDELMU = 0.D0
ELSE IF (J0.EQ.11) THEN
AVDELMU = 12*(4*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V4+11*(V-1)*V**2*V4+4*G&
&TR*(V**4-3*V**3+4*V**2-2*V+1)*V3+11*(V-1)*V**2*V3-11*(V**4-3*V**3+4*V*&
&*2-2*V+1)*V2-4*GTR*(V-1)*V**2*V2-11*(V**4-3*V**3+4*V**2-2*V+1)*V1)/((V&
&-1)**2*V)/9.D0/(8.D0*N**2)
ELSE IF (J0.EQ.12) THEN
AVDELMU = (132*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC-48*GTR*N**3*(V-1)**4&
&*(2*V**2-2*V+1)**2*VC-132*N**2*(V-1)**4*(2*V**2-2*V+1)*VC+48*GTR*N*(V-&
&1)**4*(2*V**2-2*V+1)*VC)/(N**2*(V-1)**5*V**2)/18.D0/(8.D0*N**2)
ELSE IF (J0.EQ.13) THEN
AVDELMU = (-132*N**4*(V-1)**4*(V**2+1)**2*VC+48*GTR*N**3*(V-1)**4*(V**&
&2+1)**2*VC+132*N**2*(V-1)**6*(V**2+1)*VC-48*GTR*N*(V-1)**6*(V**2+1)*VC&
&)/(N**2*(V-1)**6*V**2)/18.D0/(8.D0*N*VC)
ELSE IF (J0.EQ.14) THEN
AVDELMU = (44*N**4*(V-1)**4*(V**2-2*V+2)**2*VC-16*GTR*N**3*(V-1)**4*(V&
&**2-2*V+2)**2*VC-44*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC+16*GTR*N*(V-1)*&
&*4*V**2*(V**2-2*V+2)*VC)/(N**2*(V-1)**5*V**3)/6.D0/(8.D0*N*VC)
ELSE IF (J0.EQ.15) THEN
AVDELMU =(192*GTR*N**2*(V-1)**4*(V**2-V+1)**3*VC-528*N**3*(V-1)**4*(V*&
&*2-V+1)**3*VC)/((V-1)**6*V**3)/9.D0/(8.D0*VC**2)
ELSE IF (J0.EQ.16) THEN
AVDELMU = (44*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC-16*GTR*N**3*(V-1)**4*&
&(2*V**2-2*V+1)**2*VC-44*N**2*(V-1)**4*(2*V**2-2*V+1)*VC+16*GTR*N*(V-1)&
&**4*(2*V**2-2*V+1)*VC)/(N**2*(V-1)**5*V**2)/6.D0/(8.D0*VC**2)
ENDIF
RETURN
END FUNCTION
FUNCTION FDEL1(V,X3)
double precision::BX1
double precision::BX2
double precision::FDEL1
double precision::FKEL
double precision::SHD
double precision::V
double precision::X3
BX1=GV*GW/V/X3
BX2=(1.-GV)/(1.-V)/X3
SHD=BX1*BX2*S
FKEL=AVDEL(V,SHD)
FDEL1=FKEL/SHD
RETURN
END FUNCTION
FUNCTION FDEL2(V,X3)
double precision::BX1
double precision::BX2
double precision::FDEL2
double precision::FKELC
double precision::SHD
double precision::UN
double precision::V
double precision::X3
BX1=GV*GW/V/X3
BX2=(1.-GV)/(1.-V)/X3
SHD=BX1*BX2*S
UN=1.D0
FKELC=(AVDEL(1.-V,SHD)+DLOG(V/(1.-V))*AVWPL(UN,1.-V,SHD)+.5*AVLO(UN,1.&
&-V,SHD)*(DLOG((1.-V)/V))**2)*(1.-V)/V
FDEL2=FKELC/SHD
RETURN
END FUNCTION
FUNCTION FVWPL1(W,V,X3)
double precision::FVWPL1
double precision::RVWPL
double precision::SH
double precision::V
double precision::W
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
SH=X1*X2*S
RVWPL=AVWPL(W,V,SH)
FVWPL1=RVWPL/SH
RETURN
END FUNCTION
FUNCTION FVWPL2(W,V,X3)
double precision::FVWPL2
double precision::RVWPLC
double precision::SH
double precision::V
double precision::VX
double precision::W
double precision::WX
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
VX=1.-V*W
WX=(1.-V)/(1.-V*W)
SH=X1*X2*S
RVWPLC=(AVWPL(WX,VX,SH)+AVLO(WX,VX,SH)*DLOG(V/VX))*VX/V
FVWPL2=RVWPLC/SH
RETURN
END FUNCTION
FUNCTION FVLO1(W,V,X3)
double precision::FVLO1
double precision::RVWLO
double precision::SH
double precision::V
double precision::W
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
SH=X1*X2*S
RVWLO=AVLO(W,V,SH)
FVLO1=RVWLO/SH
RETURN
END FUNCTION
FUNCTION FVLO2(W,V,X3)
double precision::FVLO2
double precision::RVWLOC
double precision::SH
double precision::V
double precision::VX
double precision::W
double precision::WX
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
VX=1.-V*W
WX=(1.-V)/(1.-V*W)
SH=X1*X2*S
RVWLOC=AVLO(WX,VX,SH)*VX/V
FVLO2=RVWLOC/SH
RETURN
END FUNCTION
FUNCTION FRESC1(W,V,X3)
double precision::FRESC1
double precision::RRESC
double precision::SH
double precision::V
double precision::W
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
SH=X1*X2*S
RRESC=(STRUV(W,V,X3,SH)+AVGO(W,V))
FRESC1=RRESC/SH
RETURN
END FUNCTION
FUNCTION FRESC2(W,V,X3)
double precision::FRESC2
double precision::RRESCC
double precision::SH
double precision::V
double precision::VX
double precision::W
double precision::WX
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
VX=1.-V*W
WX=(1.-V)/(1.-V*W)
SH=X1*X2*S
RRESCC=STRUV(WX,VX,X3,SH)+AVGO(WX,VX)
FRESC2=RRESCC/SH
RETURN
END FUNCTION
FUNCTION FRESCC1(W,V,X3)
double precision::FRESCC1
double precision::RRESC
double precision::SH
double precision::V
double precision::W
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
SH=X1*X2*S
RRESC=STRUVC(W,V,X3,SH)
FRESCC1=RRESC/SH
RETURN
END FUNCTION
FUNCTION FRESCC2(W,V,X3)
double precision::FRESCC2
double precision::RRESCC
double precision::SH
double precision::V
double precision::VX
double precision::W
double precision::WX
double precision::X1
double precision::X2
double precision::X3
X1=GV*GW/V/W/X3
X2=(1.-GV)/(1.-V)/X3
VX=1.-V*W
WX=(1.-V)/(1.-V*W)
SH=X1*X2*S
RRESCC=STRUVC(WX,VX,X3,SH)*V/VX
FRESCC2=RRESCC/SH
RETURN
END FUNCTION
FUNCTION FQQD(X)
double precision::FQQD
double precision::PI2
double precision::X
PI2=PI**2
FQQD=4./3.*(-9./2.+2.*PI2/3.)
RETURN
END FUNCTION
FUNCTION FQQW(X)
double precision::FQQW
double precision::X
FQQW=4./3.*(-3./2.+2.*(1.+X**2)*LOG(X)+(1.-X)**2*3./2.)
RETURN
END FUNCTION
FUNCTION FQQL(X)
double precision::FQQL
double precision::X
FQQL=4./3.*(1.+X**2)
RETURN
END FUNCTION
FUNCTION FQGL(X)
double precision::FQGL
double precision::X
FQGL=0.
RETURN
END FUNCTION
FUNCTION FQGD(X)
double precision::FQGD
double precision::X
FQGD=0.
RETURN
END FUNCTION
FUNCTION FQGW(X)
double precision::FQGW
double precision::X
FQGW=0.
RETURN
END FUNCTION
FUNCTION FGGL(X)
double precision::FGGL
double precision::X
FGGL=0.
RETURN
END FUNCTION
FUNCTION FGGD(X)
double precision::FGGD
double precision::X
FGGD=0.
RETURN
END FUNCTION
FUNCTION FGGW(X)
double precision::FGGW
double precision::X
FGGW=0.
RETURN
END FUNCTION
FUNCTION FGQL(X)
double precision::FGQL
double precision::X
FGQL=0.
RETURN
END FUNCTION
FUNCTION FGQD(X)
double precision::FGQD
double precision::X
FGQD=0.
RETURN
END FUNCTION
FUNCTION FGQW(X)
double precision::FGQW
double precision::X
FGQW=0.
RETURN
END FUNCTION
FUNCTION CGQD(X)
double precision::CGQD
double precision::PI2
double precision::X
PI2=PI**2
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGQD=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGQD=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGQD=-4./3.D0*(9./2.D0+PI2/3.D0)
ENDIF
RETURN
END FUNCTION
FUNCTION CGQW(X)
double precision::CGQW
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGQW=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGQW=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGQW=4./3.D0*(-3./2.D0-(1.+X**2)*DLOG(X)+(1.-X)*(3.+2.*X))
ENDIF
RETURN
END FUNCTION
FUNCTION CGQL(X)
double precision::CGQL
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGQL=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGQL=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGQL=4./3.D0*(1.+X**2)
ENDIF
RETURN
END FUNCTION
FUNCTION CQQD(X)
double precision::CQQD
double precision::PI2
double precision::X
PI2=PI**2
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQQD=-(9./2.D0+PI2/3.D0)*CF
ELSE IF (JMAR.EQ.1) THEN
CQQD=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQQD=0.D0
ENDIF
RETURN
END FUNCTION
FUNCTION CQQW(X)
double precision::CQQW
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQQW=CF*(-3./2.D0+(3.+2.*X)*(1.-X)-(1.+X**2)*DLOG(X))
ELSE IF (JMAR.EQ.1) THEN
CQQW=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQQW=0.D0
ENDIF
RETURN
END FUNCTION
FUNCTION CQQL(X)
double precision::CQQL
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQQL=4./3.D0*(1.+X**2)
ELSE IF (JMAR.EQ.1) THEN
CQQL=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQQL=0.D0
ENDIF
RETURN
END FUNCTION
FUNCTION CQGD(X)
double precision::CQGD
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQGD=0.D0
ELSE IF (JMAR.EQ.1) THEN
CQGD=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQGD=0.D0
ENDIF
RETURN
END FUNCTION
FUNCTION CQGW(X)
double precision::CQGW
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQGW=0.D0
ELSE IF (JMAR.EQ.1) THEN
CQGW=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQGW=-(1.-X)/2.D0*(-(X**2+(1.-X)**2)*DLOG(X)+8.*X*(1.-X)-1.)
ENDIF
RETURN
END FUNCTION
FUNCTION CQGL(X)
double precision::CQGL
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CQGL=0.D0
ELSE IF (JMAR.EQ.1) THEN
CQGL=0.D0
ELSE IF (JMAR.EQ.3) THEN
CQGL=-(1.-X)/2.D0*(X**2+(1.-X)**2)
ENDIF
RETURN
END FUNCTION
FUNCTION CGGD(X)
double precision::CGGD
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGGD=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGGD=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGGD=0.D0
ENDIF
RETURN
END FUNCTION
FUNCTION CGGW(X)
double precision::CGGW
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGGW=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGGW=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGGW=2.*NF*(1.-X)/2.*(-(X**2+(1.-X)**2)*DLOG(X)+8.*X*(1.-X)-1.)
ENDIF
RETURN
END FUNCTION
FUNCTION CGGL(X)
double precision::CGGL
double precision::X
IF (JMAR.EQ.0.OR.JMAR.EQ.2) THEN
CGGL=0.D0
ELSE IF (JMAR.EQ.1) THEN
CGGL=0.D0
ELSE IF (JMAR.EQ.3) THEN
CGGL=2.*NF*(1.-X)/2.D0*(X**2+(1.-X)**2)
ENDIF
RETURN
END FUNCTION
FUNCTION AVWPL(W,V,S)
double precision::AVWPL 
double precision::M
double precision::MP
double precision::S
double precision::V
double precision::VY
double precision::VZ
double precision::W
double precision::WPLUS
double precision::ZWPL
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
WPLUS=1.
IF (J0.EQ.1) THEN
AVWPL =(-2*LOG(1-V)*(V-1)*(V**2+1)*(4*V2+V1)-2*(V-1)*(V**2+1)*LOG(V)*(&
&4*V2-5*V1)+2*LOG(S/MP**2)*(V-1)*(V**2+1)*V1+4*LOG(S/M**2)*(V-1)*(V**2+&
&1)*V1+3*(V-1)*(V**2+1)*V1)*WPLUS/((V-1)**3*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=A0(V,S)/V/W*HQQW(W,V,1)+A0(V*W,S)/(1.-V)*HQQW(W,V,2)+E0(W*V,S)*HG&
&QW(W,V,2)/(1.-V)*N/VC+A0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.2) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))*N/VC+HF&
&GQW(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.3) THEN
AVWPL =(2*(V-1)*(V**2+1)*LOG(V)*(8*V2+V1)-2*LOG(1-V)*(V-1)*(V**2+1)*(4&
&*V2+V1)+2*LOG(S/MP**2)*(V-1)*(V**2+1)*V1+4*LOG(S/M**2)*(V-1)*(V**2+1)*&
&V1+3*(V-1)*(V**2+1)*V1)*WPLUS/((V-1)**3*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=A0(V,S)/V/W*HQQW(W,V,1)+A0(V*W,S)/(1.-V)*HQQW(W,V,2)+E0(W*V,S)*HG&
&QW(W,V,2)/(1.-V)*N/VC+A0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.4) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))*N/VC+HF&
&GQW(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.5) THEN
AVWPL =(2*(V-1)*(2*V**2-2*V+1)*LOG(V)*(8*V2+V1)-2*LOG(1-V)*(V-1)*(2*V*&
&*2-2*V+1)*(4*V2-3*V1)+2*LOG(S/MP**2)*(V-1)*(2*V**2-2*V+1)*V1+4*LOG(S/M&
&**2)*(V-1)*(2*V**2-2*V+1)*V1+3*(V-1)*(2*V**2-2*V+1)*V1)*WPLUS/((V-1)*V&
&)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=A2(1.-V,S)/V/W*HQQW(W,V,1)+A2(1.-V*W,S)/(1.-V)*HQQW(W,V,2)+A2(1.-&
&VZ,S)*HFQQW(W,V)/(1.-V+V*W)+D1(VZ,S)*HFQGW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.6) THEN
AVWPL =(4*(V-1)*LOG(V)*(6*(V-1)*V*V4-4*(V-1)*V*V3-4*(V**4-2*V**3+4*V**&
&2-3*V+1)*V2+(3*V**4-2*V**3+6*V**2-3*V+1)*V1)+4*LOG(1-V)*(V-1)*(2*(V-1)&
&*V*V4-4*(V-1)*V*V3-4*(V**4-2*V**3+4*V**2-3*V+1)*V2+(V**4-6*V**3+10*V**&
&2-9*V+3)*V1)+4*LOG(S/MP**2)*(V-1)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*&
&V+1)*V1)+8*LOG(S/M**2)*(V-1)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*&
&V1)+6*(V-1)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1))*WPLUS/((V-1)&
&**3*V**3)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=B0(V,S)/V/W*HQQW(W,V,1)+B0(V*W,S)/(1.-V)*HQQW(W,V,2)+E0(1.-V,S)/W&
&/V*N/VC*HGQW(W,V,1)+E0(W*V,S)/(1.-V)*N/VC*HGQW(W,V,2)+B0(VZ,S)*HFQQW(W&
&,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.7) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=(E0(V,S)/V/W*HGQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQW(W,V,2))*N/VC+HF&
&GQW(W,V)/(1.-V+V*W)*B0(VZ,S)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.8) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=N/VC*D1(V,S)/V/W*HGQW(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)/N*VC&
&+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGW(W,V,2)/(1.-V)*A2(1.-V*W,S&
&)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.9) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=N/VC*D1(VY,S)/V/W*HGQW(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)/N*V&
&C+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGW(W,V,2)/(1.-V)*A2(V*W,S)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.10) THEN
AVWPL =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=N/VC*D1(VY,S)/V/W*HGQW(W,V,1)+D0(1.-V*W,S)/(1.-V)*HQGW(W,V,2)/N*V&
&C+HFQGW(W,V)/(1.-V+V*W)*E0(1.-VZ,S)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.11) THEN
AVWPL =-2*(2*LOG(V)*(2*(V-1)*V**2*V4+8*(V-1)*V**2*V3-8*(V**4-3*V**3+4*&
&V**2-2*V+1)*V2+(-V**4+3*V**3-4*V**2+2*V-1)*V1)+2*LOG(1-V)*(2*(V-1)*V**&
&2*V4-4*(V-1)*V**2*V3+4*(V**4-3*V**3+4*V**2-2*V+1)*V2+(-3*V**4+9*V**3-1&
&0*V**2+6*V-1)*V1)+2*LOG(S/MP**2)*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2&
&+2*V-1)*V1)+4*LOG(S/M**2)*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2+2*V-1)&
&*V1)+3*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1))*WPLUS/((V-1)*&
&*2*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=D0(V,S)/V/W*HQQW(W,V,1)+D0(V*W,S)/(1.-V)*HQQW(W,V,2)+E0(W*V,S)/(1&
&.-V)*N/VC*HGQW(W,V,2)+D0(VZ,S)*HFQQW(W,V)/(1.-V+V*W)+D1(VZ,S)*HFQGW(W,&
&V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.12) THEN
AVWPL =(LOG(V)*(6*N**4*(V-1)*(2*V**2-2*V+1)*(2*CQ*(2*V**2-2*V+1)-10*V*&
&*2+14*V-7)*VC+12*N**2*(V-1)*(V**2-V-CQ+4)*(2*V**2-2*V+1)*VC-6*(V-1)*(2&
&*V**2-2*V+1)*VC)+LOG(1-V)*(-6*N**4*(V-1)*(2*V**2-2*V+1)*(2*V**2+2*V-1)&
&*VC-12*N**2*(V-1)*(V**2-V-1)*(2*V**2-2*V+1)*VC+6*(V-1)*(2*V**2-2*V+1)*&
&VC)+LOG(S/M**2)*(-12*N**4*(V-1)*(2*V**2-2*V+1)**2*VC+24*N**2*(V-1)*(V*&
&*2-V+1)*(2*V**2-2*V+1)*VC-12*(V-1)*(2*V**2-2*V+1)*VC)+LOG(S/MP**2)*(12&
&*N**2*(V-1)*(2*V**2-2*V+1)*VC-12*N**4*(V-1)*(2*V**2-2*V+1)**2*VC)+N**2&
&*(V-1)*(2*V**2-2*V+1)*(18*V**2-18*V+7)*VC+2*N**4*(V-1)*(2*V**2-2*V+1)*&
&*2*VC-4*GTR*N**3*(V-1)*(2*V**2-2*V+1)**2*VC+4*GTR*N*(V-1)*(2*V**2-2*V+&
&1)*VC-9*(V-1)*(2*V**2-2*V+1)*VC)*WPLUS/(N**2*(V-1)**2*V**2)/3.D0
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=D1(V,S)/V/W*HQQW(W,V,1)+D1(V*W,S)/(1.-V)*HQQW(W,V,2)+2.*(2.*GTR-1&
&.)*A2(VZ,S)/(1.-V+V*W)*HFGQW(W,V)+E0(1.-V*W,S)/(1.-V)*N/VC*HGQW(W,V,2)&
&+N/VC*E0(V,S)/W/V*HGQW(W,V,1)+D1(VZ,S)*HFGGW(W,V)/(1.-V+V*W)+(D0(1.-VZ&
&,S)+D0(VZ,S))*HFGQW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.13) THEN
AVWPL = (LOG(1-V)*(-12*N**2*(V-1)*(V**2+1)*(CQ*(V-1)**2-2*(V**2-3*V+1)&
&)*VC+12*(CQ-1)*N**4*(V-1)*(V**2+1)**2*VC+12*(V-1)**3*(V**2+1)*VC)+LOG(&
&V)*(-6*N**4*(V-1)*(V**2+1)*(2*CQ*(V**2+1)-3*V**2-7)*VC+12*N**2*(V-1)*(&
&-4*V**2+7*V+CQ*(V-1)**2-4)*(V**2+1)*VC+6*(V-1)**3*(V**2+1)*VC)+LOG(S/M&
&**2)*(-12*N**2*(V-1)*(V**2+1)*(2*V**2-3*V+2)*VC+18*N**4*(V-1)*(V**2+1)&
&**2*VC+6*(V-1)**3*(V**2+1)*VC)+LOG(S/MP**2)*(-12*N**2*(V-1)*(V**2+1)*(&
&V**2-V+1)*VC+6*N**4*(V-1)*(V**2+1)**2*VC+6*(V-1)**3*(V**2+1)*VC)-N**2*&
&(V-1)*(V**2+1)*(7*V**2+4*V+7)*VC-2*N**4*(V-1)*(V**2+1)**2*VC+4*GTR*N**&
&3*(V-1)*(V**2+1)**2*VC-4*GTR*N*(V-1)**3*(V**2+1)*VC+9*(V-1)**3*(V**2+1&
&)*VC)*WPLUS/(N**2*(V-1)**3*V**2)/3.D0
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=E0(V,S)/V/W*HQQW(W,V,1)+E0(V*W,S)/(1.-V)*HGGW(W,V,2)+2.*(2.*GTR-1&
&.)*VC/N*A0(V*W,S)/(1.-V)*HQGW(W,V,2)+(B0(V*W,S)+D0(V*W,S))/(1.-V)*VC/N&
&*HQGW(W,V,2)+N/VC*D1(V,S)/W/V*HGQW(W,V,1)+E0(VZ,S)*HFQQW(W,V)/(1.-V+V*&
&W)+E0(1.-VZ,S)*HFQGW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.14) THEN
AVWPL = (LOG(1-V)*(4*CQ*N**2*(V-1)**2*V**2*(V**2-2*V+2)*VC-4*N**4*(V-1&
&)**2*(V**2-2*V+2)*(CQ*(V**2-2*V+2)-2*(V-1)**2)*VC)+LOG(V)*(-4*N**2*(V-&
&1)**2*(V**2-2*V+2)*(2*CQ*V**2-3*V**2+V-1)*VC+2*(4*CQ-9)*N**4*(V-1)**2*&
&(V**2-2*V+2)**2*VC-2*(V-1)**2*V**2*(V**2-2*V+2)*VC)+LOG(S/M**2)*(4*N**&
&2*(V-1)**2*(V**2-2*V+2)*(2*V**2-V+1)*VC-6*N**4*(V-1)**2*(V**2-2*V+2)**&
&2*VC-2*(V-1)**2*V**2*(V**2-2*V+2)*VC)+LOG(S/MP**2)*(4*N**2*(V-1)**2*V*&
&*2*(V**2-2*V+2)*VC-4*N**4*(V-1)**2*(V**2-2*V+2)**2*VC))*WPLUS/(N**2*(V&
&-1)**3*V**3)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=E0(1.-V,S)/V/W*HQQW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGGW(W,V,2)+N/VC*F&
&2(V,S)/W/V*HGQW(W,V,1)+VC/N*D1(V*W,S)/(1.-V)*HQGW(W,V,2)+E0(VZ,S)*HFGQ&
&W(W,V)/(1.-V+V*W)+E0(1.-VZ,S)*HFGGW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.15) THEN
AVWPL = (-96*N**3*(V-1)*(V**2-V+1)**2*(2*CQ*(V**2-V+1)-4*V**2+5*V-5)*L&
&OG(V)*VC+96*N**3*LOG(1-V)*(V-1)*(V**2-V+1)**2*(CQ*(V**2-V+1)-(V-1)**2)&
&*VC+96*N**3*LOG(S/MP**2)*(V-1)*(V**2-V+1)**3*VC+192*N**3*LOG(S/M**2)*(&
&V-1)*(V**2-V+1)**3*VC-88*N**3*(V-1)*(V**2-V+1)**3*VC+32*GTR*N**2*(V-1)&
&*(V**2-V+1)**3*VC)*WPLUS/((V-1)**3*V**3)/3.D0
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=F2(V,S)/V/W*HGGW(W,V,1)+F2(V*W,S)/(1.-V)*HGGW(W,V,2)+4.*GTR*VC/N*&
&(E0(VY,S)/W/V*HQGW(W,V,1)+E0(V*W,S)/(1.-V)*HQGW(W,V,2))+F2(VZ,S)*HFGGW&
&(W,V)/(1.-V+V*W)+4.*GTR*D1(VZ,S)*HFGQW(W,V)/(1.-V+V*W)
AVWPL=AVWPL+ZWPL
ELSE IF (J0.EQ.16) THEN
AVWPL = (LOG(V)*(4*N**4*(V-1)**2*(2*V**2-2*V+1)*(CQ*(2*V**2-2*V+1)-2*(&
&3*V**2-4*V+2))*VC-4*(CQ-4)*N**2*(V-1)**2*(2*V**2-2*V+1)*VC)+LOG(1-V)*(&
&4*CQ*N**2*(V-1)**2*(2*V**2-2*V+1)*VC-4*N**4*(V-1)**2*(2*V**2-2*V+1)*(C&
&Q*(2*V**2-2*V+1)-2*(V-1)**2)*VC)+LOG(S/MP**2)*(-2*N**4*(V-1)**2*(2*V**&
&2-2*V+1)**2*VC+4*N**2*(V-1)**2*(V**2-V+1)*(2*V**2-2*V+1)*VC-2*(V-1)**2&
&*(2*V**2-2*V+1)*VC)+LOG(S/M**2)*(8*N**2*(V-1)**2*(2*V**2-2*V+1)*VC-8*N&
&**4*(V-1)**2*(2*V**2-2*V+1)**2*VC))*WPLUS/(N**2*(V-1)**3*V**2)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZWPL=D1(V,S)/V/W*HGGW(W,V,1)+D1(V*W,S)/(1.-V)*HGGW(W,V,2)+VC/N*(E0(V,S&
&)/W/V*HQGW(W,V,1)+E0(1.-V*W,S)/(1.-V)*HQGW(W,V,2))+D1(VZ,S)*HFQQW(W,V)&
&/(1.-V+V*W)+F2(VZ,S)*HFQGW(W,V)/(1.-V+V*W)
AVWPL=ZWPL+AVWPL
ENDIF
RETURN
END FUNCTION
FUNCTION AVDEL(V,S)
double precision::AVDEL 
double precision::DELW
double precision::M
double precision::MP
double precision::S
double precision::V
double precision::ZDEL
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
DELW=1.
IF (J0.EQ.1) THEN
AVDEL =DELW*(3*LOG(1-V)*(V-1)**4*(V**2+1)*(16*GTR*(V**2+1)*V4+16*GTR*(&
&V**2+1)*V3-4*(17*V**2-3*V+8)*V2+(-17*V**2+12*V-29)*V1)-80*GTR*(V-1)**4&
&*(V**2+1)**2*(V4+V3)+2*(V-1)**4*(V**2+1)*((9*PI**2*V**2+170*V**2+27* P&
&I**2+170)*V2+(24* PI**2*V**2+179*V**2+6* PI**2+179)*V1)+36*LOG(1-V)*(V&
&-1)**4*(V**2+1)*LOG(V)*((3*V**2+1)*V2+(-5*V**2-3)*V1)+18*LOG(1-V)**2*(&
&V-1)**4*(V**2+1)*((V**2+7)*V2+(3*V**2+1)*V1)+9*(V-1)**4*(V**2+1)*LOG(V&
&)*(4*(V-1)*V2+(3*V**2-4*V+7)*V1)-18*(V-1)**4*(V**2+1)*(9*V**2+7)*LOG(V&
&)**2*(V2-V1)+LOG(S/M**2)*(36*(V-1)**4*(V**2+1)**2*LOG(V)*V1-36*LOG(1-V&
&)*(V-1)**4*(V**2+1)**2*V1+54*(V-1)**4*(V**2+1)**2*V1)+LOG(S/MP**2)*(36&
&*(V-1)**4*(V**2+1)**2*LOG(V)*V1+27*(V-1)**4*(V**2+1)**2*V1))/((V-1)**6&
&*V*(V**2+1))/18.D0
ZDEL=A0(V,S)/V*HQQD(V,1)+A0(V,S)/(1.-V)*HQQD(V,2)+E0(V,S)*HGQD(V,2)/(1&
&.-V)*N/VC+A0(V,S)*HFQQD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.2) THEN
AVDEL =1.D-30
ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))*N/VC+HFGQD(V)*(&
&A0(V,S)+A0(1.-V,S))
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.3) THEN
AVDEL =DELW*(3*LOG(1-V)*(V-1)**4*(V**2+1)*(16*GTR*(V**2+1)*V4+16*GTR*(&
&V**2+1)*V3-4*(8*V**2-3*V+17)*V2+(-29*V**2+12*V-17)*V1)-80*GTR*(V-1)**4&
&*(V**2+1)**2*(V4+V3)+2*(V-1)**4*(V**2+1)*(2*(18* PI**2*V**2+85*V**2+85&
&)*V2+(15* PI**2+179)*(V**2+1)*V1)+18*LOG(1-V)**2*(V-1)**4*(V**2+1)*((7&
&*V**2+1)*V2+(V**2+3)*V1)-9*(V-1)**4*(V**2+1)*LOG(V)*(8*(V-1)*V2-3*(V**&
&2+1)*V1)+36*(V-1)**4*(V**2+1)*(9*V**2+7)*LOG(V)**2*V2-144*LOG(1-V)*(V-&
&1)**4*(V**2+1)*(3*V**2+2)*LOG(V)*V2+LOG(S/M**2)*(36*(V-1)**4*(V**2+1)*&
&*2*LOG(V)*V1-36*LOG(1-V)*(V-1)**4*(V**2+1)**2*V1+54*(V-1)**4*(V**2+1)*&
&*2*V1)+LOG(S/MP**2)*(36*(V-1)**4*(V**2+1)**2*LOG(V)*V1+27*(V-1)**4*(V*&
&*2+1)**2*V1))/((V-1)**6*V*(V**2+1))/18.D0
ZDEL=A0(V,S)/V*HQQD(V,1)+A0(V,S)/(1.-V)*HQQD(V,2)+E0(V,S)*HGQD(V,2)/(1&
&.-V)*N/VC+A0(V,S)*HFQQD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.4) THEN
AVDEL =1.D-30
ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))*N/VC+HFGQD(V)*(&
&A0(V,S)+A0(1.-V,S))
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.5) THEN
AVDEL =DELW*( -80*GTR*(V-1)**3*(2*V**2-2*V+1)**2*(V4+V3)-9*LOG(1-V)*(V&
&-1)**3*(2*V**2-2*V+1)*(4*V*V2+(6*V**2-10*V+3)*V1)-9*(V-1)**3*(2*V**2-2&
&*V+1)*LOG(V)*(8*(V-1)*V2-3*(2*V**2-2*V+1)*V1)-2*(V-1)**3*(2*V**2-2*V+1&
&)**2*(2*(9* PI**2-85)*V2+(3* PI**2-179)*V1)-72*LOG(1-V)*(V-1)**3*(2*V*&
&*2-2*V+1)**2*LOG(V)*(3*V2-V1)+18*LOG(1-V)**2*(V-1)**3*(2*V-1)*(2*V**2-&
&2*V+1)*(V2-V1)+36*(V-1)**3*(2*V**2-2*V+1)*(16*V**2-14*V+7)*LOG(V)**2*V&
&2+LOG(S/M**2)*(36*(V-1)**3*(2*V**2-2*V+1)**2*LOG(V)*V1-36*LOG(1-V)*(V-&
&1)**3*(2*V**2-2*V+1)**2*V1+54*(V-1)**3*(2*V**2-2*V+1)**2*V1)+LOG(S/MP*&
&*2)*(36*(V-1)**3*(2*V**2-2*V+1)**2*LOG(V)*V1+27*(V-1)**3*(2*V**2-2*V+1&
&)**2*V1))/((V-1)**3*V*(2*V**2-2*V+1))/18.D0
ZDEL=A2(1.-V,S)/V*HQQD(V,1)+A2(1.-V,S)/(1.-V)*HQQD(V,2)+A2(1.-V,S)*HFQ&
&QD(V)+D1(V,S)*HFQGD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.6) THEN
AVDEL =DELW*(LOG(S/M**2)*(36*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*(2*&
&(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)-36*LOG(1-V)*(V-1)**4*(V**2+1&
&)*(V**2-2*V+2)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)+54*(V-1)**&
&4*(V**2+1)*(V**2-2*V+2)*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1))+&
&LOG(S/MP**2)*(36*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*(2*(V-1)*V*V4+(&
&V**4-2*V**3+4*V**2-3*V+1)*V1)+27*(V-1)**4*(V**2+1)*(V**2-2*V+2)*(2*(V-&
&1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1))-80*GTR*(V-1)**4*(V**2+1)*(V**2&
&-2*V+2)*((V**4-2*V**3+4*V**2-3*V+1)*V4+(V**4-2*V**3+4*V**2-3*V+1)*V3+(&
&V-1)*V*V2)+(V-1)**4*(V**2+1)*(V**2-2*V+2)*((V-1)*V*(18* PI**2*V**2-18*&
&PI**2*V+15* PI**2+376)*V4-(V-1)*V*(18* PI**2*V**2-18* PI**2*V-45* PI*&
&*2-340)*V3+2*(9*PI**2*V**4+170*V**4-18* PI**2*V**3-340*V**3+54* PI**2*&
&V**2+680*V**2-45* PI**2*V-510*V+18* PI**2+170)*V2+2*(24* PI**2*V**4+17&
&9*V**4-48* PI**2*V**3-358*V**3+78* PI**2*V**2+716*V**2-54* PI**2*V-537&
&*V+15* PI**2+179)*V1)+9*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)**2*((V-1&
&)*V*(2*V**2-2*V+15)*V4-(V-1)*V*(2*V**2-2*V+11)*V3-2*(8*V**4-14*V**3+25&
&*V**2-15*V+4)*V2+2*(6*V**4-6*V**3+13*V**2-7*V+2)*V1)-18*LOG(1-V)*(V-1)&
&**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*((V-1)*V*(2*V**2-2*V+3)*V4-(V-1)*V*(2&
&*V**2-2*V+3)*V3-2*(V**2-V+2)*(3*V**2-3*V+1)*V2+2*V*(3*V**3-2*V**2+4*V-&
&1)*V1)+9*LOG(1-V)**2*(V-1)**4*V*(V**2+1)*(V**2-2*V+2)*((V-1)*(2*V**2-2&
&*V-1)*V4-(V-1)*(2*V**2-2*V-5)*V3+2*(2*V**2+V+1)*V2+2*(2*V**3-2*V**2+3*&
&V-1)*V1)+3*(V-1)**4*(V**2+1)*(V**2-2*V+2)*LOG(V)*(8*GTR*(V-1)**2*(V**2&
&-2*V+2)*V4-2*(V-1)*V*(6*V-19)*V4+8*GTR*(V-1)**2*(V**2-2*V+2)*V3+4*(V-1&
&)*V*(3*V-7)*V3-2*(V-1)*(17*V**3-51*V**2+53*V-22)*V2+8*GTR*(V-1)*V*V2+(&
&5*V**4-14*V**3+26*V**2-9*V+1)*V1)+3*LOG(1-V)*(V-1)**4*(V**2+1)*(V**2-2&
&*V+2)*(8*GTR*V**2*(V**2+1)*V4+2*(V-1)*V*(6*V-5)*V4+8*GTR*V**2*(V**2+1)&
&*V3-4*(V-1)*V*(3*V+4)*V3-2*V*(17*V**3+2*V+3)*V2+8*GTR*(V-1)*V*V2+(-13*&
&V**4+30*V**3-58*V**2+33*V-9)*V1))/((V-1)**6*V**3*(V**2+1)*(V**2-2*V+2)&
&)/9.D0
ZDEL=B0(V,S)/V*HQQD(V,1)+B0(V,S)/(1.-V)*HQQD(V,2)+E0(1.-V,S)/V*N/VC*HG&
&QD(V,1)+E0(V,S)/(1.-V)*N/VC*HGQD(V,2)+B0(V,S)*HFQQD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.7) THEN
AVDEL =1.D-30
ZDEL=(E0(V,S)/V*HGQD(V,1)+E0(1.-V,S)/(1.-V)*HGQD(V,2))*N/VC+HFGQD(V)*B&
&0(V,S)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.8) THEN
AVDEL = 1.D-30
ZDEL=N/VC*D1(V,S)/V*HGQD(V,1)+A0(1.-V,S)/(1.-V)*HQGD(V,2)/N*VC+HFQGD(V&
&)*E0(1.-V,S)+VC/N*HQGD(V,2)/(1.-V)*A2(1.-V,S)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.9) THEN
AVDEL = 1.D-30
ZDEL=N/VC*D1(1.-V,S)/V*HGQD(V,1)+A0(1.-V,S)/(1.-V)*HQGD(V,2)/N*VC+HFQG&
&D(V)*E0(1.-V,S)+VC/N*HQGD(V,2)/(1.-V)*A2(V,S)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.10) THEN
AVDEL = 1.D-30
ZDEL=N/VC*D1(1.-V,S)/V*HGQD(V,1)+D0(1.-V,S)/(1.-V)*HQGD(V,2)/N*VC+HFQG&
&D(V)*E0(1.-V,S)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.11) THEN
AVDEL =-DELW*(2*(40*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V4+4*(3* PI**2+47)*&
&(V-1)*V**2*V4+40*GTR*(V**4-3*V**3+4*V**2-2*V+1)*V3+170*(V-1)*V**2*V3+(&
&18* PI**2*V**4-170*V**4-54* PI**2*V**3+510*V**3+45*PI**2*V**2-680*V**2&
&-36* PI**2*V+340*V+9* PI**2-170)*V2-40*GTR*(V-1)*V**2*V2+(3* PI**2*V**&
&4-179*V**4-9* PI**2*V**3+537*V**3+3*PI**2*V**2-716*V**2-6* PI**2*V+358&
&*V-6* PI**2-179)*V1)-9*LOG(1-V)**2*((V-1)*(V**2+2*V-2)*V4-(V-1)*(5*V**&
&2+2*V-2)*V3+2*V*(V**2+V+2)*V2-2*(V**3-3*V**2+2*V-2)*V1)-3*LOG(1-V)*(8*&
&GTR*(V**2+1)*V4+2*(V-1)*V*(5*V-6)*V4+8*GTR*(V**2+1)*V3+4*(V-1)*V*(4*V+&
&3)*V3-2*(3*V**3+2*V**2+17)*V2-8*GTR*(V-1)*V**2*V2+(-9*V**4+33*V**3-58*&
&V**2+30*V-13)*V1)+9*LOG(V)*(6*(V-1)*V**2*V4+4*(V-1)*(V**2-2*V+2)*V2-3*&
&(V**4-3*V**3+4*V**2-2*V+1)*V1)+36*LOG(1-V)*LOG(V)*(2*(V-1)*V**2*V4-8*(&
&V-1)*V**2*V3+(6*V**4-18*V**3+27*V**2-12*V+7)*V2-(V-1)**2*(2*V**2-2*V+1&
&)*V1)+9*LOG(S/MP**2)*(4*LOG(V)+3)*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**&
&2+2*V-1)*V1)+18*LOG(S/M**2)*(2*LOG(V)-2*LOG(1-V)+3)*(2*(V-1)*V**2*V4+(&
&-V**4+3*V**3-4*V**2+2*V-1)*V1)+36*LOG(V)**2*(8*(V-1)*V**2*V3+(-8*V**4+&
&23*V**3-30*V**2+14*V-7)*V2))/((V-1)**2*V)/9.D0
ZDEL=D0(V,S)/V*HQQD(V,1)+D0(V,S)/(1.-V)*HQQD(V,2)+E0(V,S)/(1.-V)*N/VC*&
&HGQD(V,2)+D0(V,S)*HFQQD(V)+D1(V,S)*HFQGD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.12) THEN
AVDEL =DELW*(LOG(S/M**2)*(LOG(1-V)*(36*N**4*(V-1)**4*(2*V**2-2*V+1)**2&
&*VC-72*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+1)*VC+36*(V-1)**4*(2*V**2-&
&2*V+1)*VC)+LOG(V)*(-36*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC+72*N**2*(V-1&
&)**4*(V**2-V+1)*(2*V**2-2*V+1)*VC-36*(V-1)**4*(2*V**2-2*V+1)*VC)-54*N*&
&*4*(V-1)**4*(2*V**2-2*V+1)**2*VC+108*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-&
&2*V+1)*VC-54*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(S/MP**2)*(LOG(V)*(72*N**2&
&*(V-1)**4*(2*V**2-2*V+1)*VC-72*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC)-66*&
&N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC+24*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)&
&**2*VC+66*N**2*(V-1)**4*(2*V**2-2*V+1)*VC-24*GTR*N*(V-1)**4*(2*V**2-2*&
&V+1)*VC)+LOG(V)**2*(18*N**4*(V-1)**4*(2*CQ*(2*V**2-2*V+1)**2-24*V**4+6&
&0*V**3-65*V**2+36*V-9)*VC-18*N**2*(V-1)**4*(V**3+2*CQ*(2*V**2-2*V+1)-2&
&0*V**2+21*V-10)*VC-18*(V-1)**4*(V**2+1)*VC)+LOG(V)*(3*N**4*(V-1)**4*(5&
&2*V**4-74*V**3+26*V**2+14*V-5)*VC+6*N**2*(V-1)**4*(18*V**4-39*V**3+53*&
&V**2-50*V+16)*VC-9*(V-1)**4*(8*V**2-14*V+9)*VC-24*GTR*N**3*(V-1)**4*(2&
&*V**2-2*V+1)**2*VC+24*GTR*N*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(1-V)*(9*N*&
&*4*(V-1)**4*(12*V**4-34*V**3+28*V**2-12*V+3)*VC-18*N**2*(V-1)**4*(6*V*&
&*4-13*V**3+8*V**2-7*V+3)*VC+9*(V-1)**4*(4*V**2-10*V+3)*VC)+LOG(1-V)**2&
&*(18*N**4*(V-1)**4*V*(4*V**2-5*V+2)*VC+18*N**2*(V-1)**5*(V**2-2*V+2)*V&
&C-18*(V-1)**4*(V**2-2*V+2)*VC)+LOG(1-V)*LOG(V)*(-72*N**4*(V-1)**4*(2*V&
&-1)*(2*V**2-2*V+1)*VC-72*N**2*(V-1)**4*(2*V**2-2*V+1)*VC)+2*N**4*(V-1)&
&**4*(3*(4*PI**2-51)*CQ*(2*V**2-2*V+1)**2-36* PI**2*V**4-376*V**4+72* P&
&I**2*V**3+752*V**3-72* PI**2*V**2-725*V**2+36* PI**2*V+349*V-9* PI**2-&
&85)*VC+2*N**2*(V-1)**4*(60* PI**2*V**4+36*V**4-120* PI**2*V**3-72*V**3&
&-3*(4* PI**2-51)*CQ*(2*V**2-2*V+1)+138* PI**2*V**2+269*V**2-78* PI**2*&
&V-233*V+24* PI**2+103)*VC+(63*CQ+40)*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)*&
&*2*VC-(63*CQ+40)*GTR*N*(V-1)**4*(2*V**2-2*V+1)*VC-6*(5* PI**2+6)*(V-1)&
&**4*(2*V**2-2*V+1)*VC)/(N**2*(V-1)**5*V**2)/18.D0
ZDEL=D1(V,S)/V*HQQD(V,1)+D1(V,S)/(1.-V)*HQQD(V,2)+2.*(2.*GTR-1.)*A2(V,&
&S)*HFGQD(V)+E0(1.-V,S)/(1.-V)*N/VC*HGQD(V,2)+N/VC*E0(V,S)/V*HGQD(V,1)+&
&D1(V,S)*HFGGD(V)+(D0(1.-V,S)+D0(V,S))*HFGQD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.13) THEN
AVDEL = DELW*(LOG(S/MP**2)*(LOG(V)*(-72*N**2*(V-1)**4*(V**2+1)*(V**2-V&
&+1)*VC+36*N**4*(V-1)**4*(V**2+1)**2*VC+36*(V-1)**6*(V**2+1)*VC)-54*N**&
&2*(V-1)**4*(V**2+1)*(V**2-V+1)*VC+27*N**4*(V-1)**4*(V**2+1)**2*VC+27*(&
&V-1)**6*(V**2+1)*VC)+LOG(S/M**2)*(LOG(V)*(72*N**4*(V-1)**4*(V**2+1)**2&
&*VC-72*N**2*(V-1)**6*(V**2+1)*VC)+LOG(1-V)*(72*N**2*(V-1)**6*(V**2+1)*&
&VC-72*N**4*(V-1)**4*(V**2+1)**2*VC)-6*N**2*(V-1)**4*(V**2+1)*(20*V**2-&
&31*V+20)*VC+93*N**4*(V-1)**4*(V**2+1)**2*VC-24*GTR*N**3*(V-1)**4*(V**2&
&+1)**2*VC+24*GTR*N*(V-1)**6*(V**2+1)*VC+27*(V-1)**6*(V**2+1)*VC)+LOG(1&
&-V)*LOG(V)*(36*N**4*(V-1)**4*(2*CQ*(V**2+1)**2-V**4-2*V**3-5*V**2-4)*V&
&C-36*N**2*(V-1)**4*(-2*V**4+6*V**3+2*CQ*(V-1)**2*(V**2+1)-5*V**2+9*V-4&
&)*VC-36*(V-2)*(V-1)**6*V*VC)+LOG(V)**2*(-18*N**4*(V-1)**4*(2*CQ*(V**2+&
&1)**2-2*V**4-2*V**3-11*V**2-9)*VC+18*N**2*(V-1)**5*(-8*V**3+2*CQ*(V-1)&
&*(V**2+1)+8*V**2-9*V+10)*VC+18*(V-1)**6*(2*V**2-2*V+1)*VC)+LOG(V)*(-3*&
&N**4*(V-1)**4*(13*V**4+38*V**2+6*V-5)*VC+6*N**2*(V-1)**4*(2*V**4-19*V*&
&*3+V**2+14*V-16)*VC+9*(V-1)**6*(3*V**2-4*V+9)*VC+24*GTR*N**3*(V-1)**4*&
&(V**2+1)**2*VC-24*GTR*N*(V-1)**6*(V**2+1)*VC)+LOG(1-V)**2*(-18*N**4*(V&
&-1)**4*(V**2+1)*(2*CQ*(V**2+1)-(V+1)**2)*VC+18*N**2*(V-1)**6*(2*CQ*(V*&
&*2+1)-V)*VC+18*(V-1)**6*(3*V**2-4*V+3)*VC)+LOG(1-V)*(18*N**4*(V-1)**4*&
&V*(V**2+10*V+1)*VC-18*N**2*(V-1)**4*V*(V**2+10*V+1)*VC+72*(V-1)**6*V*V&
&C)+2*N**4*(V-1)**4*(6*CQ*(PI**2+3)*(V**2+1)**2+9*PI**2*V**4+85*V**4+18&
&*PI**2*V**3+9*V**3+9*PI**2*V**2+188*V**2+9*V+85)*VC-2*N**2*(V-1)**4*(-&
&12*PI**2*V**4+103*V**4+18*PI**2*V**3-179*V**3+6*CQ*(PI**2+3)*(V-1)**2*&
&(V**2+1)-15*PI**2*V**2+188*V**2-9*PI**2*V-179*V+6*PI**2+103)*VC+6*(V-1&
&)**6*(5*PI**2*V**2+6*V**2-6*PI**2*V+2*PI**2+6)*VC-10*(3*CQ+4)*GTR*N**3&
&*(V-1)**4*(V**2+1)**2*VC+10*(3*CQ+4)*GTR*N*(V-1)**6*(V**2+1)*VC)/(N**2&
&*(V-1)**6*V**2)/18.D0
ZDEL=E0(V,S)/V*HQQD(V,1)+E0(V,S)/(1.-V)*HGGD(V,2)+2.*(2.*GTR-1.)*VC/N*&
&A0(V,S)/(1.-V)*HQGD(V,2)+(B0(V,S)+D0(V,S))/(1.-V)*VC/N*HQGD(V,2)+N/VC*&
&D1(V,S)/V*HGQD(V,1)+E0(V,S)*HFQQD(V)+E0(1.-V,S)*HFQGD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.14) THEN
AVDEL = DELW*(LOG(S/M**2)*(LOG(1-V)*(24*N**4*(V-1)**4*(V**2-2*V+2)**2*&
&VC-24*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC)+LOG(V)*(24*N**2*(V-1)**4*V**&
&2*(V**2-2*V+2)*VC-24*N**4*(V-1)**4*(V**2-2*V+2)**2*VC)+2*N**2*(V-1)**4&
&*(V**2-2*V+2)*(20*V**2-9*V+9)*VC-31*N**4*(V-1)**4*(V**2-2*V+2)**2*VC+8&
&*GTR*N**3*(V-1)**4*(V**2-2*V+2)**2*VC-8*GTR*N*(V-1)**4*V**2*(V**2-2*V+&
&2)*VC-9*(V-1)**4*V**2*(V**2-2*V+2)*VC)+LOG(S/MP**2)*(LOG(V)*(24*N**2*(&
&V-1)**4*V**2*(V**2-2*V+2)*VC-24*N**4*(V-1)**4*(V**2-2*V+2)**2*VC)-22*N&
&**4*(V-1)**4*(V**2-2*V+2)**2*VC+8*GTR*N**3*(V-1)**4*(V**2-2*V+2)**2*VC&
&+22*N**2*(V-1)**4*V**2*(V**2-2*V+2)*VC-8*GTR*N*(V-1)**4*V**2*(V**2-2*V&
&+2)*VC)+LOG(1-V)**2*(6*N**4*(V-1)**4*(2*CQ*(V**2-2*V+2)**2-2*V**4+10*V&
&**3-21*V**2+20*V-8)*VC-6*N**2*(V-1)**4*V*(2*CQ*V*(V**2-2*V+2)-V-1)*VC-&
&6*(V-1)**4*V**2*(2*V**2-2*V+1)*VC)+LOG(V)*(-6*N**2*(V-1)**4*(3*V**4-8*&
&V**3+2*V**2+12*V-6)*VC+3*N**4*(V-1)**4*(3*V**4-10*V**3-2*V**2+24*V-12)&
&*VC+3*(V-1)**4*V**2*(3*V**2+2*V-2)*VC)+LOG(1-V)*LOG(V)*(12*N**2*(V-1)*&
&*4*V*((V-1)*(2*V**2-2*V+1)+2*CQ*V*(V**2-2*V+2))*VC-12*N**4*(V-1)**4*(2&
&*CQ*(V**2-2*V+2)**2-4*V**4+18*V**3-37*V**2+36*V-16)*VC+12*(V-1)**4*V**&
&2*(2*V**2-2*V+1)*VC)+LOG(V)**2*(-6*N**2*(V-1)**4*V**2*(4*CQ*(V**2-2*V+&
&2)-2*V**2+5*V-5)*VC+6*N**4*(V-1)**4*(V**2-2*V+2)*(4*CQ*(V**2-2*V+2)-11&
&*V**2+24*V-24)*VC-6*(V-1)**4*V**2*(3*V**2-2*V+2)*VC)+LOG(1-V)*(-6*N**2&
&*(V-1)**4*V*(2*V**2-7*V-1)*VC-6*(V-1)**4*V**2*(2*V+1)*VC+6*N**4*(V-1)*&
&*4*V*(2*V-5)*VC)+2*N**4*(V-1)**4*(CQ*(2*PI**2-57)*(V**2-2*V+2)**2-3*(2&
&*PI**2*V**4+V**4-10*PI**2*V**3-5*V**3+21*PI**2*V**2+13*V**2-20*PI**2*V&
&-16*V+8*PI**2+8))*VC-2*N**2*(V-1)**4*V*(CQ*(2*PI**2-57)*V*(V**2-2*V+2)&
&-3*(2*V**3-5*V**2+PI**2*V+5*V+PI**2))*VC-6*(V-1)**4*V**2*(2*PI**2*V**2&
&+V**2-2*PI**2*V-2*V+PI**2+2)*VC+31*CQ*GTR*N**3*(V-1)**4*(V**2-2*V+2)**&
&2*VC-31*CQ*GTR*N*(V-1)**4*V**2*(V**2-2*V+2)*VC)/(N**2*(V-1)**5*V**3)/6&
&.D0
ZDEL=E0(1.-V,S)/V*HQQD(V,1)+E0(1.-V,S)/(1.-V)*HGGD(V,2)+N/VC*F2(V,S)/V&
&*HGQD(V,1)+VC/N*D1(V,S)/(1.-V)*HQGD(V,2)+E0(V,S)*HFGQD(V)+E0(1.-V,S)*H&
&FGGD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.15) THEN
AVDEL = DELW*(LOG(S/M**2)*(288*N**3*(V-1)**4*(V**2-V+1)**3*LOG(V)*VC-2&
&88*N**3*LOG(1-V)*(V-1)**4*(V**2-V+1)**3*VC+528*N**3*(V-1)**4*(V**2-V+1&
&)**3*VC-192*GTR*N**2*(V-1)**4*(V**2-V+1)**3*VC)+LOG(S/MP**2)*(288*N**3&
&*(V-1)**4*(V**2-V+1)**3*LOG(V)*VC+264*N**3*(V-1)**4*(V**2-V+1)**3*VC-9&
&6*GTR*N**2*(V-1)**4*(V**2-V+1)**3*VC)+LOG(1-V)*LOG(V)*(72*N**3*(V-1)**&
&4*(4*CQ*(V**2-V+1)**3-4*V**6+16*V**5-37*V**4+50*V**3-47*V**2+26*V-8)*V&
&C+72*GTR*N**2*(V-1)**5*V*(2*V**2-2*V+1)*VC)+LOG(V)**2*(-72*N**3*(V-1)*&
&*4*(4*CQ*(V**2-V+1)**3-8*V**6+27*V**5-57*V**4+72*V**3-66*V**2+36*V-12)&
&*VC-36*GTR*N**2*(V-1)**5*V**2*(V**2+V-1)*VC)+LOG(1-V)**2*(36*GTR*N**2*&
&(V-1)**6*V*(V**2-3*V+1)*VC-72*N**3*(V-1)**4*(2*CQ*(V**2-V+1)**3-2*V**6&
&+7*V**5-14*V**4+16*V**3-14*V**2+7*V-2)*VC)+LOG(V)*(24*GTR*N**2*(V-1)**&
&4*(V**2-V+1)*(4*V**4-3*V**3-V**2+8*V-4)*VC-24*N**3*(V-1)**4*(V**2-V+1)&
&*(11*V**4-15*V**3+4*V**2+22*V-11)*VC)+LOG(1-V)*(24*N**3*(V-1)**4*V*(V*&
&*2-V+1)*(7*V**2+8*V+7)*VC-24*GTR*N**2*(V-1)**4*V*(V**2-V+1)*(5*V**2-2*&
&V+5)*VC)+4*N**3*(V-1)**4*(378*CQ*(V**2-V+1)**3+24*PI**2*V**6-134*V**6-&
&72*PI**2*V**5+402*V**5+153*PI**2*V**4-831*V**4-186*PI**2*V**3+992*V**3&
&+171*PI**2*V**2-831*V**2-90*PI**2*V+402*V+24*PI**2-134)*VC-4*GTR*N**2*&
&(V-1)**4*(123*CQ*(V**2-V+1)**3-40*V**6+120*V**5+18*PI**2*V**4-294*V**4&
&-36*PI**2*V**3+388*V**3+27*PI**2*V**2-294*V**2-9*PI**2*V+120*V-40)*VC)&
&/((V-1)**6*V**3)/9.D0
ZDEL=F2(V,S)/V*HGGD(V,1)+F2(V,S)/(1.-V)*HGGD(V,2)+VC*(E0(1.-V,S)/V*HQG&
&D(V,1)+E0(V,S)/(1.-V)*HQGD(V,2))*4.*GTR/N+F2(V,S)*HFGGD(V)+4.*GTR*D1(V&
&,S)*HFGQD(V)
AVDEL=AVDEL+ZDEL
ELSE IF (J0.EQ.16) THEN
AVDEL = DELW*(LOG(S/M**2)*(LOG(1-V)*(24*N**4*(V-1)**4*(2*V**2-2*V+1)**&
&2*VC-24*N**2*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(V)*(24*N**2*(V-1)**4*(2*V&
&**2-2*V+1)*VC-24*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC)-44*N**4*(V-1)**4*&
&(2*V**2-2*V+1)**2*VC+16*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC+44*N**2&
&*(V-1)**4*(2*V**2-2*V+1)*VC-16*GTR*N*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(S&
&/MP**2)*(LOG(V)*(-12*N**4*(V-1)**4*(2*V**2-2*V+1)**2*VC+24*N**2*(V-1)*&
&*4*(V**2-V+1)*(2*V**2-2*V+1)*VC-12*(V-1)**4*(2*V**2-2*V+1)*VC)-9*N**4*&
&(V-1)**4*(2*V**2-2*V+1)**2*VC+18*N**2*(V-1)**4*(V**2-V+1)*(2*V**2-2*V+&
&1)*VC-9*(V-1)**4*(2*V**2-2*V+1)*VC)+LOG(1-V)**2*(6*N**4*(V-1)**4*(2*CQ&
&*(2*V**2-2*V+1)**2-8*V**4+20*V**3-21*V**2+10*V-2)*VC+6*N**2*(V-1)**4*(&
&V**2*(V+1)-2*CQ*(2*V**2-2*V+1))*VC-6*(V-1)**4*(V**2-2*V+2)*VC)+LOG(V)*&
&*2*(6*N**4*(V-1)**4*(2*CQ*(2*V**2-2*V+1)**2-24*V**4+60*V**3-65*V**2+36&
&*V-9)*VC-6*N**2*(V-1)**4*(V**3+2*CQ*(2*V**2-2*V+1)-20*V**2+21*V-10)*VC&
&-6*(V-1)**4*(V**2+1)*VC)+LOG(1-V)*LOG(V)*(24*(CQ-2)*N**2*(V-1)**4*(2*V&
&**2-2*V+1)*VC-24*N**4*(V-1)**4*(2*V**2-2*V+1)*(CQ*(2*V**2-2*V+1)-2*(V-&
&1)**2)*VC)+LOG(1-V)*(6*N**2*(V-1)**4*V*(V**2+7*V-2)*VC-6*N**4*(V-1)**4&
&*V**2*(5*V-2)*VC-6*(V-1)**4*V*(V+2)*VC)+LOG(V)*(-6*N**2*(V-1)**5*(V**2&
&-9*V+6)*VC+6*N**4*(V-1)**6*(5*V-3)*VC-6*(V-3)*(V-1)**5*VC)-2*N**4*(V-1&
&)**4*(3*(8*V**4-16*V**3+13*V**2-5*V+1)+4*CQ*(PI**2+3)*(2*V**2-2*V+1)**&
&2)*VC+2*N**2*(V-1)**4*(3*(5*V**2-5*V+2)+4*CQ*(PI**2+3)*(2*V**2-2*V+1))&
&*VC+20*CQ*GTR*N**3*(V-1)**4*(2*V**2-2*V+1)**2*VC-20*CQ*GTR*N*(V-1)**4*&
&(2*V**2-2*V+1)*VC-6*(V-1)**4*(2*V**2-2*V+1)*VC)/(N**2*(V-1)**5*V**2)/6&
&.D0
ZDEL=D1(V,S)/V*HGGD(V,1)+D1(V,S)/(1.-V)*HGGD(V,2)+VC*(E0(V,S)/V*HQGD(V&
&,1)+E0(1.-V,S)/(1.-V)*HQGD(V,2))/N+D1(V,S)*HFQQD(V)+F2(V,S)*HFQGD(V)
AVDEL=AVDEL+ZDEL
ENDIF
RETURN
END FUNCTION
FUNCTION AVLO(W,V,S)
double precision::AVLO 
double precision::LOPLUS
double precision::M
double precision::MP
double precision::S
double precision::V
double precision::VY
double precision::VZ
double precision::W
double precision::ZVLO
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LOPLUS=1.
IF (J0.EQ.1) THEN
AVLO = 4*LOPLUS*(V**2+1)*V1/((V-1)**2*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=A0(V,S)/V/W*HQQL(W,V,1)+A0(V*W,S)/(1.-V)*HQQL(W,V,2)+E0(W*V,S)*HG&
&QL(W,V,2)/(1.-V)*N/VC+A0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.2) THEN
AVLO =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))*N/VC+HF&
&GQL(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.3) THEN
AVLO = 4*LOPLUS*(V**2+1)*V1/((V-1)**2*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=A0(V,S)/V/W*HQQL(W,V,1)+A0(V*W,S)/(1.-V)*HQQL(W,V,2)+E0(W*V,S)*HG&
&QL(W,V,2)/(1.-V)*N/VC+A0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.4) THEN
AVLO =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))*N/VC+HF&
&GQL(W,V)/(1.-V+V*W)*(A0(VZ,S)+A0(1.-VZ,S))
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.5) THEN
AVLO = 4*LOPLUS*(2*V**2-2*V+1)*V1/V
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=A2(1.-V,S)/V/W*HQQL(W,V,1)+A2(1.-V*W,S)/(1.-V)*HQQL(W,V,2)+A2(1.-&
&VZ,S)*HFQQL(W,V)/(1.-V+V*W)+D1(VZ,S)*HFQGL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.6) THEN
AVLO =8*LOPLUS*(2*(V-1)*V*V4+(V**4-2*V**3+4*V**2-3*V+1)*V1)/((V-1)**2*&
&V**3)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=B0(V,S)/V/W*HQQL(W,V,1)+B0(V*W,S)/(1.-V)*HQQL(W,V,2)+E0(1.-V,S)/W&
&/V*N/VC*HGQL(W,V,1)+E0(W*V,S)/(1.-V)*N/VC*HGQL(W,V,2)+B0(VZ,S)*HFQQL(W&
&,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.7) THEN
AVLO =1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=(E0(V,S)/V/W*HGQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGQL(W,V,2))*N/VC+HF&
&GQL(W,V)/(1.-V+V*W)*B0(VZ,S)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.8) THEN
AVLO = 1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=N/VC*D1(V,S)/V/W*HGQL(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)/N*VC&
&+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGL(W,V,2)/(1.-V)*A2(1.-V*W,S&
&)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.9) THEN
AVLO = 1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=N/VC*D1(VY,S)/V/W*HGQL(W,V,1)+A0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)/N*V&
&C+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)+VC/N*HQGL(W,V,2)/(1.-V)*A2(V*W,S)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.10) THEN
AVLO = 1.D-30
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=N/VC*D1(VY,S)/V/W*HGQL(W,V,1)+D0(1.-V*W,S)/(1.-V)*HQGL(W,V,2)/N*V&
&C+HFQGL(W,V)/(1.-V+V*W)*E0(1.-VZ,S)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.11) THEN
AVLO =-8*LOPLUS*(2*(V-1)*V**2*V4+(-V**4+3*V**3-4*V**2+2*V-1)*V1)/((V-1&
&)**2*V)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=D0(V,S)/V/W*HQQL(W,V,1)+D0(V*W,S)/(1.-V)*HQQL(W,V,2)+E0(W*V,S)/(1&
&.-V)*N/VC*HGQL(W,V,2)+D0(VZ,S)*HFQQL(W,V)/(1.-V+V*W)+D1(VZ,S)*HFQGL(W,&
&V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.12) THEN
AVLO =LOPLUS*(4*N**2*(V-1)*(2*V**2-2*V+1)*(2*V**2-2*V-CQ+3)*VC+4*(CQ-2&
&)*N**4*(V-1)*(2*V**2-2*V+1)**2*VC-4*(V-1)*(2*V**2-2*V+1)*VC)/(N**2*(V-&
&1)**2*V**2)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=D1(V,S)/V/W*HQQL(W,V,1)+D1(V*W,S)/(1.-V)*HQQL(W,V,2)+2.*(2.*GTR-1&
&.)*A2(VZ,S)/(1.-V+V*W)*HFGQL(W,V)+E0(1.-V*W,S)/(1.-V)*N/VC*HGQL(W,V,2)&
&+N/VC*E0(V,S)/W/V*HGQL(W,V,1)+D1(VZ,S)*HFGGL(W,V)/(1.-V+V*W)+(D0(1.-VZ&
&,S)+D0(VZ,S))*HFGQL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.13) THEN
AVLO = LOPLUS*(-4*(CQ-2)*N**4*(V-1)*(V**2+1)**2*VC+4*N**2*(V-1)*(-3*V*&
&*2+4*V+CQ*(V-1)**2-3)*(V**2+1)*VC+4*(V-1)**3*(V**2+1)*VC)/(N**2*(V-1)*&
&*3*V**2)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=E0(V,S)/V/W*HQQL(W,V,1)+E0(V*W,S)/(1.-V)*HGGL(W,V,2)+2.*(2.*GTR-1&
&.)*VC/N*A0(V*W,S)/(1.-V)*HQGL(W,V,2)+(B0(V*W,S)+D0(V*W,S))/(1.-V)*VC/N&
&*HQGL(W,V,2)+N/VC*D1(V,S)/W/V*HGQL(W,V,1)+E0(VZ,S)*HFQQL(W,V)/(1.-V+V*&
&W)+E0(1.-VZ,S)*HFQGL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.14) THEN
AVLO = LOPLUS*(8*(CQ-2)*N**2*(V-1)**2*(V**2-2*V+2)**2*VC-8*(CQ-2)*(V-1&
&)**2*V**2*(V**2-2*V+2)*VC)/((V-1)**3*V**3)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=E0(1.-V,S)/V/W*HQQL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HGGL(W,V,2)+N/VC*F&
&2(V,S)/W/V*HGQL(W,V,1)+VC/N*D1(V*W,S)/(1.-V)*HQGL(W,V,2)+E0(VZ,S)*HFGQ&
&L(W,V)/(1.-V+V*W)+E0(1.-VZ,S)*HFGGL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.15) THEN
AVLO = -32*(3*CQ-5)*LOPLUS*N**3*(V**2-V+1)**3*VC/((V-1)**2*V**3)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=F2(V,S)/V/W*HGGL(W,V,1)+F2(V*W,S)/(1.-V)*HGGL(W,V,2)+4.*GTR*VC/N*&
&(E0(VY,S)/W/V*HQGL(W,V,1)+E0(V*W,S)/(1.-V)*HQGL(W,V,2))+F2(VZ,S)*HFGGL&
&(W,V)/(1.-V+V*W)+4.*GTR*D1(VZ,S)*HFGQL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ELSE IF (J0.EQ.16) THEN
AVLO = LOPLUS*(8*(CQ-2)*N**2*(V-1)**2*(2*V**2-2*V+1)**2*VC-8*(CQ-2)*(V&
&-1)**2*(2*V**2-2*V+1)*VC)/((V-1)**3*V**2)
VZ=V*W/(1.-V+V*W)
VY=1.-V
ZVLO=D1(V,S)/V/W*HGGL(W,V,1)+D1(V*W,S)/(1.-V)*HGGL(W,V,2)+VC/N*(E0(V,S&
&)/W/V*HQGL(W,V,1)+E0(1.-V*W,S)/(1.-V)*HQGL(W,V,2))+D1(VZ,S)*HFQQL(W,V)&
&/(1.-V+V*W)+F2(VZ,S)*HFQGL(W,V)/(1.-V+V*W)
AVLO=AVLO+ZVLO
ENDIF
RETURN
END FUNCTION
FUNCTION AVGO(W,V)
double precision::AVGO 
double precision::LOPLUS
double precision::LOTVW
double precision::LOVW
double precision::LOW
double precision::M
double precision::MP
double precision::V
double precision::W
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LOVW=LOG((1.-V*W)/(1.-V))/(1.-W)
LOTVW=LOG(1.-V+V*W)/(1.-W)
LOW=LOG(W)/(1.-W)
LOPLUS=1.
IF (J0.EQ.1) THEN
AVGO =-2*((5*LOW*V**2+LOVW*V**2-2*LOTVW*V**2+3*LOW+7*LOVW+2*LOTVW)*V2+&
&(-6*LOW*V**2+5*LOVW*V**2+2*LOTVW*V**2-4*LOW+3*LOVW+2*LOTVW)*V1)/((V-1)&
&**2*V)
ELSE IF (J0.EQ.2) THEN
AVGO =1.D-30
ELSE IF (J0.EQ.3) THEN
AVGO =2*((10*LOW*V**2-7*LOVW*V**2-LOTVW*V**2+6*LOW-LOVW+LOTVW)*V2+(LOW&
&*V**2-3*LOVW*V**2-LOTVW*V**2+LOW-5*LOVW-3*LOTVW)*V1)/((V-1)**2*V)
ELSE IF (J0.EQ.4) THEN
AVGO =1.D-30
ELSE IF (J0.EQ.5) THEN
AVGO =2*((16*LOW*V**2-8*LOTVW*V**2-12*LOW*V-2*LOVW*V+2*LOTVW*V+6*LOW+L&
&OVW-LOTVW)*V2+(2*LOW*V**2-4*LOVW*V**2-8*LOTVW*V**2-2*LOW*V+6*LOVW*V+10&
&*LOTVW*V+LOW-3*LOVW-5*LOTVW)*V1)/V
ELSE IF (J0.EQ.6) THEN
AVGO =2*((V-1)*V*(2*LOW*V**2-2*LOVW*V**2-2*LOW*V+2*LOVW*V+11*LOW-7*LOV&
&W-8*LOTVW)*V4-(V-1)*V*(2*LOW*V**2-2*LOVW*V**2-2*LOW*V+2*LOVW*V+3*LOW+5&
&*LOVW)*V3-2*V*(4*LOW*V**3-2*LOTVW*V**3-6*LOW*V**2+2*LOVW*V**2+4*LOTVW*&
&V**2+9*LOW*V+LOVW*V-4*LOTVW*V-3*LOW+LOVW+2*LOTVW)*V2+2*(5*LOW*V**4-4*L&
&OVW*V**4-2*LOTVW*V**4-8*LOW*V**3+6*LOVW*V**3+4*LOTVW*V**3+15*LOW*V**2-&
&11*LOVW*V**2-8*LOTVW*V**2-10*LOW*V+7*LOVW*V+6*LOTVW*V+3*LOW-2*LOVW-2*L&
&OTVW)*V1)/((V-1)**2*V**3)
ELSE IF (J0.EQ.7) THEN
AVGO =1.D-30
ELSE IF (J0.EQ.8) THEN
AVGO = 1.D-30
ELSE IF (J0.EQ.9) THEN
AVGO = 1.D-30
ELSE IF (J0.EQ.10) THEN
AVGO = 1.D-30
ELSE IF (J0.EQ.11) THEN
AVGO =-2*((V-1)*(4*LOW*V**2-7*LOVW*V**2-7*LOTVW*V**2+2*LOVW*V+2*LOTVW*&
&V-2*LOVW-2*LOTVW)*V4+(V-1)*(16*LOW*V**2-5*LOVW*V**2-5*LOTVW*V**2-2*LOV&
&W*V-2*LOTVW*V+2*LOVW+2*LOTVW)*V3-2*(8*LOW*V**4-4*LOTVW*V**4-22*LOW*V**&
&3-LOVW*V**3+9*LOTVW*V**3+28*LOW*V**2-LOVW*V**2-7*LOTVW*V**2-12*LOW*V-2&
&*LOVW*V+2*LOTVW*V+6*LOW)*V2-2*(LOW*V**4-2*LOVW*V**4-4*LOTVW*V**4-3*LOW&
&*V**3+7*LOVW*V**3+13*LOTVW*V**3+4*LOW*V**2-11*LOVW*V**2-17*LOTVW*V**2-&
&2*LOW*V+6*LOVW*V+10*LOTVW*V+LOW-4*LOVW-4*LOTVW)*V1)/((V-1)**2*V)
ELSE IF (J0.EQ.12) THEN
AVGO =(-2*N**4*(V-1)*(-4*CQ*LOTVW*(2*V**2-2*V+1)**2+12*LOW*V**4-8*LOVW&
&*V**4-28*LOW*V**3+20*LOVW*V**3+29*LOW*V**2-21*LOVW*V**2+2*LOTVW*V**2-1&
&6*LOW*V+10*LOVW*V-2*LOTVW*V+4*LOW-2*LOVW+LOTVW)*VC+2*N**2*(V-1)*(4*LOW&
&*V**4-8*LOVW*V**4-9*LOW*V**3+15*LOVW*V**3-4*CQ*LOTVW*(2*V**2-2*V+1)+18&
&*LOW*V**2-17*LOVW*V**2+3*LOTVW*V**2-15*LOW*V+8*LOVW*V-3*LOTVW*V+6*LOW-&
&2*LOVW+2*LOTVW)*VC-2*(V-1)*(3*LOW*V**2-5*LOVW*V**2-2*LOTVW*V**2-2*LOW*&
&V+6*LOVW*V+2*LOTVW*V+2*LOW-4*LOVW-3*LOTVW)*VC)/(N**2*(V-1)**2*V**2)
ELSE IF (J0.EQ.13) THEN
AVGO = (-2*N**2*(V-1)*(4*LOW*V**4-2*LOVW*V**4-4*LOTVW*V**4-6*LOW*V**3+&
&3*LOVW*V**3+3*LOTVW*V**3+9*LOW*V**2-2*LOVW*V**2-5*LOTVW*V**2-9*LOW*V+3&
&*LOVW*V+6*LOW-2*LOVW-2*LOTVW)*VC+2*N**4*(V-1)*(LOW*V**4-3*LOVW*V**4-LO&
&TVW*V**4+2*LOW*V**3-2*LOVW*V**3+5*LOW*V**2-6*LOVW*V**2-3*LOTVW*V**2-2*&
&LOVW*V-2*LOTVW*V+4*LOW-3*LOVW-2*LOTVW)*VC+2*(V-1)**3*(3*LOW*V**2-3*LOV&
&W*V**2-3*LOTVW*V**2-2*LOW*V+4*LOVW*V+2*LOTVW*V+2*LOW-3*LOVW-4*LOTVW)*V&
&C)/(N**2*(V-1)**3*V**2)
ELSE IF (J0.EQ.14) THEN
AVGO = (-2*N**4*(V-1)**2*(-4*CQ*LOTVW*(V**2-2*V+2)**2+4*LOW*V**4-4*LOV&
&W*V**4+LOTVW*V**4-18*LOW*V**3+18*LOVW*V**3-4*LOTVW*V**3+38*LOW*V**2-37&
&*LOVW*V**2+7*LOTVW*V**2-40*LOW*V+36*LOVW*V-4*LOTVW*V+20*LOW-16*LOVW)*V&
&C-2*N**2*(V-1)**2*(2*LOVW*V**4+3*LOW*V**3-4*LOVW*V**3-LOTVW*V**3+4*CQ*&
&LOTVW*V**2*(V**2-2*V+2)-7*LOW*V**2+5*LOVW*V**2+8*LOW*V+LOVW*V-LOTVW*V-&
&4*LOW)*VC-2*(V-1)**2*V**2*(4*LOW*V**2-2*LOVW*V**2-LOTVW*V**2-4*LOW*V+2&
&*LOVW*V+4*LOW-LOVW-LOTVW)*VC)/(N**2*(V-1)**3*V**3)
ELSE IF (J0.EQ.15) THEN
AVGO = (16*N**3*(V-1)*(2*CQ*(LOW-2*LOTVW)*(V**2-V+1)**3+2*LOW*V**6-4*L&
&OVW*V**6-7*LOW*V**5+13*LOVW*V**5+15*LOW*V**4-26*LOVW*V**4+2*LOTVW*V**4&
&-20*LOW*V**3+30*LOVW*V**3-4*LOTVW*V**3+20*LOW*V**2-26*LOVW*V**2+3*LOTV&
&W*V**2-12*LOW*V+13*LOVW*V-LOTVW*V+4*LOW-4*LOVW)*VC-8*GTR*N**2*(V-1)**2&
&*V*(LOW*V**3+LOVW*V**3+LOW*V**2-4*LOVW*V**2-LOTVW*V**2-LOW*V+4*LOVW*V+&
&LOTVW*V-LOVW+LOTVW)*VC)/((V-1)**3*V**3)
ELSE IF (J0.EQ.16) THEN
AVGO = (-2*N**4*(V-1)**2*(2*CQ*LOW*(2*V**2-2*V+1)**2+8*LOW*V**4-16*LOV&
&W*V**4-16*LOTVW*V**4-20*LOW*V**3+36*LOVW*V**3+32*LOTVW*V**3+21*LOW*V**&
&2-37*LOVW*V**2-30*LOTVW*V**2-12*LOW*V+18*LOVW*V+14*LOTVW*V+3*LOW-4*LOV&
&W-3*LOTVW)*VC-2*N**2*(V-1)**2*(LOW*V**3+LOVW*V**3-2*CQ*LOW*(2*V**2-2*V&
&+1)-8*LOW*V**2+5*LOVW*V**2+5*LOTVW*V**2+9*LOW*V-4*LOVW*V-5*LOTVW*V-4*L&
&OW+2*LOVW+2*LOTVW)*VC-2*(V-1)**2*(LOW*V**2-LOVW*V**2-2*LOTVW*V**2+2*LO&
&VW*V+2*LOTVW*V+LOW-2*LOVW-3*LOTVW)*VC)/(N**2*(V-1)**3*V**2)
ENDIF
RETURN
END FUNCTION
FUNCTION STRUV(W,V,X3,S)
double precision::S
double precision::STRUV
double precision::V
double precision::W
double precision::X3
IF (J0.EQ.1) THEN
STRUV=STRUV1(W,V,X3,S)
ELSE IF (J0.EQ.2) THEN
STRUV=STRUV2(W,V,X3,S)
ELSE IF (J0.EQ.3) THEN
STRUV=STRUV3(W,V,X3,S)
ELSE IF (J0.EQ.4) THEN
STRUV=STRUV4(W,V,X3,S)
ELSE IF (J0.EQ.5) THEN
STRUV=STRUV5(W,V,X3,S)
ELSE IF (J0.EQ.6) THEN
STRUV=STRUV6(W,V,X3,S)
ELSE IF (J0.EQ.7) THEN
STRUV=STRUV7(W,V,X3,S)
ELSE IF (J0.EQ.8) THEN
STRUV=STRUV8(W,V,X3,S)
ELSE IF (J0.EQ.9) THEN
STRUV=STRUV9(W,V,X3,S)
ELSE IF (J0.EQ.10) THEN
STRUV=STRUV10(W,V,X3,S)
ELSE IF (J0.EQ.11) THEN
STRUV=STRUV11(W,V,X3,S)
ELSE IF (J0.EQ.12) THEN
STRUV=STRUV12(W,V,X3,S)
ELSE IF (J0.EQ.13) THEN
STRUV=STRUV13(W,V,X3,S)
ELSE IF (J0.EQ.14) THEN
STRUV=STRUV14(W,V,X3,S)
ELSE IF (J0.EQ.15) THEN
STRUV=STRUV15(W,V,X3,S)
ELSE IF (J0.EQ.16) THEN
STRUV=STRUV16(W,V,X3,S)
ENDIF
RETURN
END FUNCTION
FUNCTION STRUV1(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV1
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**4*(V**2+6*V-1)*W&
&**4+2*V**4*(V+7)*W**3-2*V**2*(3*V**2+V+7)*W**2+2*V*(V**2+5)*W-V**2-1)+&
&CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*(V-&
&1)*V*V2*W*(3*V**3*W**3+V**2*W**2-V*(2*V+5)*W+3)+2*CQ*V*V2*W*(V**2*W**2&
&+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
LV1 = -LOG(1-V)*(V1*(4*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(2*V**2+9*V+1)&
&*W**4-V**3*(7*V**2+11*V+10)*W**3+V**2*(8*V**2+13*V+15)*W**2-V*(7*V**2+&
&3*V+10)*W+2*(V**2+1))+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**&
&4+2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V*&
&*2+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)*(2*V**2*W**2+V**2+4*V-3))/((V-1)**2&
&*V*W*(V*W-1)**3)
LV = LOG(V)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**3*(V**3+13*V**2-7*V&
&+1)*W**4+V**2*(5*V**3+27*V**2-3*V+3)*W**3-V*(13*V**3+7*V**2+27*V+3)*W*&
&*2+(7*V**3-V**2+21*V+1)*W-2*(V**2+1))-2*V*V2*W*(2*V**4*W**4-4*V**3*(2*&
&V-1)*W**3+V**2*(3*V**2+10*V-1)*W**2-4*V*(V**2+V+3)*W+3*V**2-2*V+9)+CQ*&
&V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*CQ*V*V&
&2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W&
&-1)**3)
LVW = 2*(V1*(2*V**3*W**3-2*(V-2)*V**2*W**2+V*(V**2+V+4)*W-V**2-1)+V*V2&
&*W*(4*V*W-V+5))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(1-V*W)&
&/((V-1)**2*V*W*(V*W-1)**4)
LTVW = 2*(V1*(2*V**2*W**2-2*(V-2)*V*W+V**2-V+4)+(1-V)*V2)*(V**4*W**4-4&
&*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(V*W-V+1)/((V-1)**2*(V*W-1)**4)
LW = -(V1*(V*(11*V**2-6*V+1)*W**2-(5*V**3+3*V**2+11*V+1)*W+3*(V**2+1))&
&+2*V*V2*W*(2*V**2*W**2-5*(V-1)*V*W+3*V**2+V+6))*(V**3*W**3-3*V**2*W**2&
&+3*V*W-1)*LOG(W)/((V-1)**2*V*W*(V*W-1)**4)
CVC = (3*(V-1)*V1*(V*W-1)*(6*V**7*W**7-4*V**6*(3*V+1)*W**6+V**4*(7*V**&
&3+6*V**2-5*V-2)*W**5-V**3*(V**4-V**3+13*V**2-37*V-4)*W**4+2*V**3*(7*V*&
&*2-22*V-2)*W**3-2*V*(7*V**4-25*V**3+17*V**2-3*V+2)*W**2-(8*V**4-27*V**&
&3+26*V**2-11*V-2)*W-(V-1)*(V**2+4*V+1))+12*(V-1)*V2*W*(V*W-1)*(2*V**6*&
&W**5-V**4*(6*V**2-3*V+1)*W**4+V**3*(6*V**3-3*V**2-V+2)*W**3-V**3*(4*V*&
&*3-16*V**2+27*V-11)*W**2-V*(4*V**4-10*V**3+3*V**2-V+2)*W+(V-1)*(2*V**3&
&-6*V**2+3*V-1))-6*AL*(V-1)*V1*(V*W-1)*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)&
&*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)-8*CQ*(V-1)**2*V1*(V*W-V+1)*(V**2*&
&W**2-2*V*W+1)*(V**4*W**4+2*V**2*W**2+1)-12*AL*(V-1)*V*V2*W*(V*W-1)*(V*&
&W-V+1)*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-16*CQ*(V-1)**2*V&
&*V2*W*(V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-2*V*W+1))/((V-1)**3*V*W*(V*W-&
&1)**4*(V*W-V+1))/6.D0
LM = LOG(S/M**2)*(-V1*(2*V**6*W**6-2*V**5*(V+1)*W**5+V**3*(V**3+2*V**2&
&+4*V+1)*W**4-V**2*(2*V**3+10*V**2+5*V+3)*W**3+V*(5*V**3+8*V**2+10*V+3)&
&*W**2-(4*V**3+2*V**2+9*V+1)*W+2*(V**2+1))-2*V*V2*W*(V**2*W**2+1)*(2*V*&
&*2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
LMP = -LOG(S/MP**2)*V1*(2*V**3*W**3-2*V**2*(2*V-3)*W**2+V*(3*V**2-8*V+&
&9)*W-(V-1)*(V**2-V+4))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)/((V&
&-1)**2*(V*W-1)**4*(V*W-V+1))
STRUV1=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV2(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV2
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V**4+&
&24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)*W**5+2*V**&
&2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V**2*(6*V**6+2*V**5&
&+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+12*V**7-17*V**6+54*V**5-5&
&*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)*(4*V**6-2*V**5+7*V**4+6*V**3+6*&
&V**2-4*V-1)*W+4*(V-1)**2*(V**2+1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*&
&(3*V+1)*W**7-V**4*(16*V**4+2*V**3+17*V**2+1)*W**6+2*V**4*(8*V**4-11*V*&
&*3+25*V**2+4*V+2)*W**5-2*V**2*(8*V**6-19*V**5+20*V**4+33*V**3-10*V**2+&
&V-1)*W**4+2*V**2*(6*V**6-9*V**5-11*V**4+63*V**3-16*V**2-6*V-3)*W**3-(4&
&*V**8+8*V**7-49*V**6+88*V**5+13*V**4-46*V**3-V**2-2*V+1)*W**2+2*(V-1)*&
&(4*V**6-6*V**5+4*V**4+11*V**3+9*V**2-5*V-1)*W-4*(V-1)**2*(V**2+1)*(V**&
&2-2*V+3))+2*V2*(V*W-1)*(V*W-V+1)*(V**3*(6*V**2-3*V+1)*W**4-V**3*(8*V**&
&2+9*V-5)*W**3+V*(6*V**4+9*V**3+15*V**2-13*V-1)*W**2-V*(9*V**3+V**2+11*&
&V-13)*W+4*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)&
&*V*W+(V-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2&
&+1)))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LV1 = LOG(1-V)*(2*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**4*W**4-2*V**3&
&*(V+2)*W**3+V**2*(V**2+2*V+9)*W**2-2*V*(V**2+2*V+3)*W+2*(V+1)*(V**2+1)&
&)+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4*V**3*W&
&**3+8*V**2*W**2-8*V*W+4)-4*V2*(V*W-V-1)*(V**2*W**2-2*V*W+1)*(V**3*W**3&
&-V**2*(V+1)*W**2+V*(3*V**2-2*V+1)*W-(V-1)*(V**2+1))-2*CQ*V2*(V*W-1)*(V&
&*W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**&
&2*V**2*W**2*(V*W-1)**2*(V*W-V+1))
LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*(5*V&
&+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6*V**3+46*V*&
&*2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-2*V-1)*W**3+2*V**2&
&*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W**2-2*(V-1)*V*(2*V**5-3*V&
&**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*(V**2+1)*(V**2-2*V+3))-V1*(-4*V*&
&*8*W**8+4*V**7*(3*V+2)*W**7-V**4*(16*V**4+12*V**3+23*V**2+1)*W**6+2*V*&
&*4*(8*V**4-8*V**3+37*V**2+3*V+2)*W**5-2*V**2*(8*V**6-22*V**5+35*V**4+3&
&7*V**3-14*V**2+V-1)*W**4+2*V**2*(6*V**6-14*V**5-5*V**4+74*V**3-20*V**2&
&-10*V-3)*W**3-(4*V**8+4*V**7-57*V**6+114*V**5+7*V**4-52*V**3-3*V**2-2*&
&V+1)*W**2+2*(V-1)*(4*V**6-10*V**5+7*V**4+10*V**3+14*V**2-8*V-1)*W-4*(V&
&-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3&
&*(11*V**2-2*V+1)*W**4+2*V**3*(4*V**2+11*V-5)*W**3-V*(V**4+16*V**3+24*V&
&**2-20*V-1)*W**2-2*V*(V**4-4*V**3-2*V**2-8*V+9)*W+2*(V-3)*(V-1)*(V**2+&
&1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**&
&2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V&
&*W-V+1)**2)
LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V&
&**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+3*V+12)*W**2-V*(V+3)**2*W+4*(V**&
&2+1))+V2*(-4*V**3*W**3+V**2*(5*V+7)*W**2-V*(3*V**2+6*V+7)*W+4*(V**2+1)&
&))*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LTVW = -4*V1*(V**2*W**2-2*V*W+1)*(CQ*(V**2*W**2-2*V**2*W+V**2+1)*(V**4&
&*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)-(V-1)*&
&V*W*(V**3*W**3-V**2*(5*V-3)*W**2+(V-1)*V*(5*V+1)*W-(V-1)**2*(V+1)))*LO&
&G(V*W-V+1)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LW = -(V**2*W**2-2*V*W+1)*(2*V1*(V*W-V+1)*(2*V**3*W**3-V**2*(V**2+7)*W&
&**2+2*V*(V**3-V**2+3*V+3)*W-2*(V**2+1)*(V**2-V+2))+4*V2*(V*W-V-1)*(V**&
&3*W**3-V**2*(3*V-1)*W**2+V*(V**2+2*V-1)*W-(V-1)*(V**2+1))-2*CQ*V*(V**2&
&+1)*V2*(V*W-V+1)*(W**2-2*W+2)-CQ*(V**2+1)**2*V1*(V*W-V+1)*(W**2-2*W+2)&
&)*LOG(W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1))
CVC = (3*V1*(8*V**7*W**7+V**4*(2*V**4-24*V**3-13*V**2+5*V-2)*W**6-V**4&
&*(4*V**4-32*V**3-46*V**2+17*V-7)*W**5+2*V**2*(V**6-12*V**5-36*V**4+7*V&
&**3+5*V**2-7*V+2)*W**4+2*V**2*(8*V**5+19*V**4+22*V**3-27*V**2+11*V-5)*&
&W**3-(8*V**7+17*V**6+3*V**5+40*V**4-66*V**3+21*V**2-9*V+2)*W**2+(V-1)*&
&(16*V**5-11*V**4+34*V**3+2*V**2-6*V-3)*W-8*(V-1)**2*V*(V**2-V+2))+3*AL&
&*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(V+5)*W**5+V*&
&*2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V+1&
&)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+V*&
&*4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))+12*V*V2*(V*W-1)*(V**5*&
&W**6-V**3*(4*V**2-V+1)*W**5+V**2*(5*V**3+2*V**2+3*V-1)*W**4-V*(5*V**4-&
&2*V**3+16*V**2-8*V-1)*W**3+(4*V**5-3*V**4+12*V**3-4*V**2-6*V+1)*W**2-(&
&V-1)*(V**4+2*V**3+4*V**2+4*V-3)*W+(V-1)**2*(V**2+3))+4*CQ*V1*(V*W-1)*(&
&V*W-V+1)*(V**5*(V+2)*W**5-V**4*(V**2+5*V-3)*W**4+V**2*(V**4+7*V**3-V**&
&2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2-5*V+1)*W**2+(V-1)*(4*V**4-6*V**&
&3+13*V**2-4*V+1)*W-(V-1)**2*(3*V**2-6*V+7))-6*AL*V2*(V*W-1)*(V*W-V-1)*&
&(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*&
&(V**2+3)*W-2*(V**2+1))-8*CQ*V2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*&
&V*W+(V-1)**2)*((V-1)*V**2*W**2-V*(V**2+2*V-1)*W+2*(V-1)))/((V-1)**2*V*&
&*2*W**2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*V&
&**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10&
&*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*&
&W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))-2*V2*&
&(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V&
&**2+1)))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LMP = 2*LOG(S/MP**2)*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V**2+1&
&)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)&
&/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
STRUV2=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV3(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV3
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**4*(V**2+2*V+3)*W&
&**4+2*V**3*(V**2+5*V+2)*W**3-2*V**2*(3*V**2+3*V+5)*W**2+2*V*(V**2+2*V+&
&3)*W-V**2-1)+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*&
&W**2+1)-2*(V-1)*V*V2*W*(3*V**3*W**3-7*V**2*W**2+V*(2*V-1)*W+3)+2*CQ*V*&
&V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*&
&W-1)**3)
LV1 = -LOG(1-V)*(V1*(4*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(2*V**2+9*V+1)&
&*W**4-V**3*(7*V**2+3*V+18)*W**3+V**2*(8*V**2-3*V+31)*W**2-V*(7*V**2-5*&
&V+18)*W+2*(V**2+1))+CQ*V1*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+&
&2*V**2*W**2+1)+2*CQ*V*V2*W*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2&
&+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)*(2*V**2*W**2+V**2-8*V+9))/((V-1)**2*V&
&*W*(V*W-1)**3)
LV = LOG(V)*(V1*(-2*V**6*W**6+2*V**5*(V+1)*W**5-V**3*(V**3+V**2+5*V+1)&
&*W**4+V**2*(V**3+7*V**2+5*V+3)*W**3-V*(5*V**3+3*V**2+7*V+3)*W**2+(3*V*&
&*3+3*V**2+5*V+1)*W-2*(V**2+1))-2*V*V2*W*(2*V**4*W**4+2*V**3*(5*V-7)*W*&
&*3-V**2*(3*V**2+20*V-11)*W**2+2*V*(4*V**2+V+9)*W-3*V**2+4*V-15)+CQ*V1*&
&(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*CQ*V*V2*W&
&*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)&
&**3)
LVW = 2*(V1*(2*V**3*W**3-2*(V-2)*V**2*W**2+V*(V**2-V+6)*W-V**2-1)+V*V2&
&*W*(4*V*W+5*V-1))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(1-V*&
&W)/((V-1)**2*V*W*(V*W-1)**4)
LTVW = 2*(V1*(2*V**2*W**2-2*(V-2)*V*W+V**2-2*V+5)+2*(V-1)*V2)*(V**4*W*&
&*4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(V*W-V+1)/((V-1)**2*(V*W-1)**4)
LW = -(V1*(V*(V**2+4*V+1)*W**2-(V+1)*(V**2+1)*W+3*(V**2+1))+2*V*V2*W*(&
&2*V**2*W**2+10*(V-1)*V*W-3*V**2-2*V-9))*(V**3*W**3-3*V**2*W**2+3*V*W-1&
&)*LOG(W)/((V-1)**2*V*W*(V*W-1)**4)
CVC = (3*(V-1)*V1*(V*W-1)*(6*V**7*W**7-4*V**6*(3*V+1)*W**6+V**4*(7*V**&
&3+2*V**2+3*V-6)*W**5-V**3*(V**4-9*V**3+21*V**2-29*V-12)*W**4-2*V**3*(4&
&*V**3-19*V**2+34*V-2)*W**3-2*V*(7*V**4-21*V**3+13*V**2-7*V+6)*W**2+(V+&
&2)*(3*V**2-4*V+3)*W-(V-1)*(V**2+4*V+1))+12*(V-1)*V2*W*(V*W-1)*(2*V**6*&
&W**5-V**4*(3*V**2+3*V-2)*W**4+V**3*(3*V**2+5*V-4)*W**3+V**3*(2*V**3-2*&
&V**2-9*V+5)*W**2-V*(4*V**4-16*V**3+9*V**2+5*V-4)*W-2*(V-1)*(2*V-1)*(V*&
&*2-V+1))-6*AL*(V-1)*V1*(V*W-1)*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2&
&+1)*(V**4*W**4+2*V**2*W**2+1)-8*CQ*(V-1)**2*V1*(V*W-V+1)*(V**2*W**2-2*&
&V*W+1)*(V**4*W**4+2*V**2*W**2+1)-12*AL*(V-1)*V*V2*W*(V*W-1)*(V*W-V+1)*&
&(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-16*CQ*(V-1)**2*V*V2*W*(&
&V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-2*V*W+1))/((V-1)**3*V*W*(V*W-1)**4*(&
&V*W-V+1))/6.D0
LM = LOG(S/M**2)*(-V1*(2*V**6*W**6-2*V**5*(V+1)*W**5+V**3*(V**3+2*V**2&
&+4*V+1)*W**4-V**2*(2*V**3+10*V**2+5*V+3)*W**3+V*(5*V**3+8*V**2+10*V+3)&
&*W**2-(4*V**3+2*V**2+9*V+1)*W+2*(V**2+1))-2*V*V2*W*(V**2*W**2+1)*(2*V*&
&*2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3)
LMP = -LOG(S/MP**2)*V1*(2*V**3*W**3-2*V**2*(2*V-3)*W**2+V*(3*V**2-8*V+&
&9)*W-(V-1)*(V**2-V+4))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)/((V&
&-1)**2*(V*W-1)**4*(V*W-V+1))
STRUV3=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV4(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV4
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V**4+&
&24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)*W**5+2*V**&
&2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V**2*(6*V**6+2*V**5&
&+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+12*V**7-17*V**6+54*V**5-5&
&*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)*(4*V**6-2*V**5+7*V**4+6*V**3+6*&
&V**2-4*V-1)*W+4*(V-1)**2*(V**2+1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*&
&(3*V+1)*W**7-V**4*(16*V**4+12*V**3+7*V**2+1)*W**6+2*V**4*(8*V**4+2*V**&
&3+19*V**2-3*V+2)*W**5-2*V**2*(8*V**6-6*V**5+27*V**4+13*V**3-10*V**2+V-&
&1)*W**4+2*V**2*(6*V**6-4*V**5+5*V**4+46*V**3-24*V**2-2*V-3)*W**3-(4*V*&
&*8+8*V**7-29*V**6+74*V**5+7*V**4-56*V**3+9*V**2-2*V+1)*W**2+2*(V-1)*(4&
&*V**6-6*V**5+9*V**4+6*V**3+12*V**2-8*V-1)*W-4*(V-1)**2*(V**2+1)*(V**2-&
&2*V+3))-2*V2*(V*W-1)*(V*W-V+1)*(V**3*(9*V**2-12*V-1)*W**4-4*V**3*(4*V*&
&*2-3*V-4)*W**3+V*(9*V**4-24*V**2-2*V+1)*W**2-2*V*(3*V**3-8*V**2-V+2)*W&
&-4*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V&
&-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1)))/(&
&(V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LV1 = LOG(1-V)*(2*V1*(V**2*W**2-2*V*W+1)*(2*V**5*W**5-2*V**4*(2*V+1)*W&
&**4+V**3*(3*V**2+4*V+5)*W**3-V**2*(V**3+5*V**2+7*V-1)*W**2+2*V*(3*V**3&
&+V**2+V-1)*W-2*(V-1)*(V+1)*(V**2+1))+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*&
&(V+1)*W+V**2+1)*(V**4*W**4-4*V**3*W**3+8*V**2*W**2-8*V*W+4)-4*V2*(V*W-&
&V-1)*(V**2*W**2-2*V*W+1)*(V**3*W**3-V**2*(V+1)*W**2+2*V*(2*V-1)*W-(V-1&
&)*(V**2+1))-2*CQ*V2*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*W**2&
&-2*V*(V+1)*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1))
LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*(5*V&
&+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6*V**3+46*V*&
&*2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-2*V-1)*W**3+2*V**2&
&*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W**2-2*(V-1)*V*(2*V**5-3*V&
&**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*(V**2+1)*(V**2-2*V+3))-V1*(-4*V*&
&*8*W**8+4*V**7*(3*V+2)*W**7-V**4*(16*V**4+24*V**3+11*V**2+1)*W**6+2*V*&
&*4*(8*V**4+6*V**3+33*V**2-7*V+2)*W**5-2*V**2*(8*V**6-10*V**5+47*V**4+1&
&5*V**3-16*V**2+V-1)*W**4+2*V**2*(6*V**6-10*V**5+11*V**4+64*V**3-36*V**&
&2-4*V-3)*W**3-(4*V**8+4*V**7-41*V**6+106*V**5+7*V**4-76*V**3+13*V**2-2&
&*V+1)*W**2+2*(V-1)*(4*V**6-10*V**5+11*V**4+6*V**3+18*V**2-12*V-1)*W-4*&
&(V-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5+V*&
&*3*(7*V**2-16*V-1)*W**4-4*V**3*(4*V**2-4*V-5)*W**3+V*(11*V**4-4*V**3-2&
&4*V**2-4*V+1)*W**2-2*V*(V**4+2*V**3-8*V**2-2*V+3)*W+2*(V-3)*(V-1)*(V**&
&2+1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-1)&
&**2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*&
&(V*W-V+1)**2)
LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V&
&**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+4*V+11)*W**2-2*V*(V**2+3*V+4)*W+&
&4*(V**2+1))-2*V2*(2*V**3*W**3-V**2*(V+5)*W**2+V*(3*V+5)*W-2*(V**2+1)))&
&*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LTVW = -4*(V**2*W**2-2*V*W+1)*(CQ*V1*(V**2*W**2-2*V**2*W+V**2+1)*(V**4&
&*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)-3*(V-1&
&)*V*V2*W*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)+2*(V-1)**2*V**3*V1&
&*(W-1)*W**2)*LOG(V*W-V+1)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LW = -(V**2*W**2-2*V*W+1)*(2*V1*(2*V**4*W**4-V**3*(V+1)*(V+3)*W**3+V**&
&2*(3*V**3-V**2+13*V-3)*W**2-2*V*(2*V-1)*(V**3-V**2+3*V+1)*W+2*(V-1)*(V&
&**2+1)*(V**2-V+2))+4*V2*(V*W-V-1)*(V**3*W**3-2*V**2*W**2+V*(V**2+2*V-1&
&)*W-(V-1)*(V**2+1))-2*CQ*V*(V**2+1)*V2*(V*W-V+1)*(W**2-2*W+2)-CQ*(V**2&
&+1)**2*V1*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/((V-1)**2*V**2*W**2*(V*W-1)**&
&2*(V*W-V+1))
CVC = (3*AL*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(V&
&+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+&
&20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2&
&*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))+6*V*V2*(V*W&
&-1)*(2*V**5*W**6-V**3*(V+1)*(5*V-1)*W**5+V**2*(7*V**3+13*V**2-3*V+1)*W&
&**4-V*(7*V**4+20*V**3-4*V**2-4*V+1)*W**3+(5*V**5+9*V**4+12*V**3-20*V**&
&2+3*V-1)*W**2-2*(V-1)*(V**4+2*V**3+4*V**2+4*V-3)*W+2*(V-1)**2*(V**2+3)&
&)+3*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(8*V**5*W**5+V**2*(2*V**4-10*V&
&**3-25*V**2+3*V-2)*W**4+V*(2*V**4+28*V**3+33*V**2-3*V+4)*W**3-(8*V**5+&
&3*V**4+37*V**3+27*V**2+3*V+2)*W**2+(16*V**4-11*V**3+39*V**2+9*V+3)*W-8&
&*V*(V**2-V+2))+4*CQ*V1*(V*W-1)*(V*W-V+1)*(V**5*(V+2)*W**5-V**4*(V**2+5&
&*V-3)*W**4+V**2*(V**4+7*V**3-V**2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2&
&-5*V+1)*W**2+(V-1)*(4*V**4-6*V**3+13*V**2-4*V+1)*W-(V-1)**2*(3*V**2-6*&
&V+7))-6*AL*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V*&
&*3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1))-8*CQ*V2*(V**2*W&
&**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*((V-1)*V**2*W**2-V*(V**2&
&+2*V-1)*W+2*(V-1)))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*V&
&**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10&
&*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*&
&W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))-2*V2*&
&(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V&
&**2+1)))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LMP = 2*LOG(S/MP**2)*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V**2+1&
&)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)&
&/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
STRUV4=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV5(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV5
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V1*(-2*V**6*W**6+2*V**&
&5*(2*V-1)*W**5-V**4*(6*V**2-8*V+3)*W**4+2*(V-1)*V**3*(8*V**2-9*V+2)*W*&
&*3-2*(V-1)**2*V**2*(11*V**2-13*V+5)*W**2+2*(V-1)**3*V*(6*V**2-8*V+3)*W&
&-(V-1)**4*(2*V**2-2*V+1))+CQ*V1*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+&
&1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)-2*(V-1)*V*V2*W*(3*V**3*W*&
&*3-7*(V-1)*V**2*W**2+(V-1)*V*(V+1)*W+3*(V-1)**3)+2*CQ*(V-1)*V*V2*W*(V*&
&*2*W**2+(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V*W&
&*(V*W-1)**3*(V*W-V+1)**3)
LV1 = -LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-&
&1)**2*V*W-(V-1)**3)*(V1*(6*V**4*W**4-12*V**3*W**3+V**2*(6*V**2-4*V+7)*&
&W**2-V*(20*V**3-48*V**2+50*V-19)*W+(V-1)*(20*V**2-33*V+18))+2*(V-1)*V2&
&*(V*W-1)*(2*V**2*W**2-5*V*W+10*V**2-17*V+10))/((V-1)*(V*W-1)**3*(V*W-V&
&+1)**4)
LV = LOG(V)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V1*(-2*V**6*W**6+2*V**5*(&
&2*V-1)*W**5-V**3*(8*V**3-14*V**2+8*V-1)*W**4+(V-1)*V**2*(16*V**3-26*V*&
&*2+14*V-3)*W**3-(V-1)**2*V*(18*V**3-26*V**2+16*V-3)*W**2+(V-1)**3*(12*&
&V**3-16*V**2+8*V-1)*W-2*(V-1)**4*(2*V**2-2*V+1))-2*(V-1)*V*V2*W*(2*V**&
&4*W**4-2*V**3*(2*V-7)*W**3-V**2*(12*V**2+2*V-11)*W**2+2*(V-1)*V*(14*V*&
&*2-19*V+9)*W-(V-1)**2*(14*V**2-26*V+15))+CQ*V1*(2*V**2*W**2-2*V*(2*V-1&
&)*W+2*V**2-2*V+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)+2*CQ*(V-1)&
&*V*V2*W*(V**2*W**2+(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))&
&/((V-1)*V*W*(V*W-1)**3*(V*W-V+1)**3)
LVW = 2*(V1*(2*V**2*W**2+2*(V-2)*V*W+4*V**2-8*V+5)+2*(V-1)*V2)*(V**3*W&
&**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*&
&W**2-4*(V-1)**3*V*W+(V-1)**4)*LOG(1-V*W)/((V-1)*(V*W-1)**3*(V*W-V+1)**&
&4)
LTVW = 2*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(CQ*V1*(2*V**2*W**2-2*V*(2*V-&
&1)*W+2*V**2-2*V+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)+(V-1)*V*V&
&2*W*(9*V**3*W**3-V**2*(8*V**2+3*V-9)*W**2+(V-1)**2*V*(16*V+3)*W-(V-1)*&
&*2*(8*V**2-7*V+1))+(V-1)*V*V1*W*(5*V**3*W**3-V**2*(8*V**2-9*V+5)*W**2+&
&(V-1)*V*(16*V**2-25*V+13)*W-(V-1)**2*(8*V**2-11*V+7))+2*CQ*(V-1)*V*V2*&
&W*(V**2*W**2+(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))*LOG(V&
&*W-V+1)/((V-1)*V*W*(V*W-1)**3*(V*W-V+1)**3)
LW = -(V1*(V*(6*V**2-6*V+1)*W**2-(2*V-1)*(2*V**2-2*V+1)*W+3*(V-1)*(2*V&
&**2-2*V+1))+2*V*V2*W*(2*V**2*W**2+10*V*W-14*V**2+20*V-9))*(V**3*W**3-3&
&*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)*&
&*3)*LOG(W)/(V*W*(V*W-1)**3*(V*W-V+1)**4)
CVC = (V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(6*V**7*W**7-4*V**6*(2*V+1)*W*&
&*6-(V-2)*V**4*(18*V**2-14*V+3)*W**5+V**3*(52*V**4-150*V**3+162*V**2-77&
&*V+12)*W**4-2*(V-1)*V**3*(21*V**3-43*V**2+22*V+2)*W**3+2*(V-1)**2*V*(2&
&*V**4+12*V**3-32*V**2+17*V-6)*W**2+(V-1)**3*(6*V**4-16*V**3+19*V**2-11&
&*V+6)*W-(V-1)**4*(6*V**2-6*V+1))+4*(V-1)*V2*W*(V*W-V+1)*(V**2*W**2-2*V&
&*W+1)*(2*V**6*W**5-V**4*(4*V**2+V-2)*W**4+(V-1)*V**3*(4*V**2+3*V-4)*W*&
&*3-V**3*(4*V**3-5*V**2-6*V+5)*W**2+(V-1)*V*(2*V**4+V**3-11*V+4)*W-2*(V&
&-1)**2*(V+1)*(V**2-V+1)))/((V-1)*V*W*(V*W-1)**3*(V*W-V+1)**4)/2.D0
LM = -LOG(S/M**2)*V1*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)&
&**3)*(2*V**4*W**4+2*(V-3)*V**3*W**3+V*(6*V**3-14*V**2+12*V-1)*W**2-(2*&
&V**4+2*V**3-8*V**2+6*V-1)*W+(V-1)*(2*V**2-2*V+1))/((V-1)*V*W*(V*W-1)*(&
&V*W-V+1)**3)
LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V1*(2*V**6*W**6-2*V**5*(2*V-1&
&)*W**5+V**4*(6*V**2-9*V+4)*W**4-(V-1)*V**3*(12*V**2-15*V+4)*W**3+(V-1)&
&**2*V**2*(14*V**2-19*V+7)*W**2-(V-1)**3*V*(8*V**2-13*V+6)*W+(V-1)**4*(&
&2*V**2-2*V+1))+2*(V-1)*V*V2*W*(V**2*W**2+(V-1)**2)*(2*V**2*W**2-2*V*(2&
&*V-1)*W+2*V**2-2*V+1))/((V-1)*V*W*(V*W-1)**2*(V*W-V+1)**3)
STRUV5=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV6(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV6
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(CQ*V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**3*(2*V**5-5*V*&
&*4+17*V**3-16*V**2+12*V-4)*W**5-V**2*(2*V**6-7*V**5+13*V**4+8*V**3-24*&
&V**2+28*V-12)*W**4+V*(2*V**7-4*V**6-V**5+25*V**4-32*V**3+16*V**2+12*V-&
&12)*W**3-(6*V**7-24*V**6+43*V**5-27*V**4-10*V**3+32*V**2-12*V-4)*W**2+&
&(6*V**6-28*V**5+63*V**4-72*V**3+41*V**2-8)*W-2*(V-1)*(V**2-2*V+2)**2)+&
&V1*(-2*V**8*W**7+2*V**7*(V+1)*W**6-V**3*(2*V**5+3*V**4+7*V**3-14*V**2+&
&12*V-4)*W**5+V**2*(2*V**6+V**5+5*V**4+30*V**3-46*V**2+36*V-12)*W**4-V*&
&(2*V**7-4*V**6+21*V**5-39*V**4+80*V**3-62*V**2+36*V-12)*W**3+(6*V**7-2&
&4*V**6+63*V**5-99*V**4+108*V**3-50*V**2+12*V-4)*W**2-V*(6*V**5-28*V**4&
&+69*V**3-96*V**2+79*V-28)*W+2*(V-1)*(V**2-2*V+2)**2)+2*CQ*V2*(2*V**7*W&
&**6-V**3*(3*V**4-2*V**3+7*V**2-6*V+2)*W**5+V**2*(3*V**5-5*V**4+5*V**3+&
&9*V**2-14*V+6)*W**4-V*(2*V**6-2*V**5-5*V**4+20*V**3-11*V**2-6*V+6)*W**&
&3+(6*V**6-17*V**5+19*V**4+3*V**3-17*V**2+6*V+2)*W**2-2*(V-1)**2*(3*V+1&
&)*(V**2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+2))+2*(V-1)*V2*W*(V**3*(3*V**3-2&
&*V**2+4*V-2)*W**4-V**2*(3*V**4-14*V**3+8*V**2+8*V-6)*W**3+V*(7*V**4-44&
&*V**3+36*V**2-6)*W**2-(9*V**4-42*V**3+40*V**2-8*V-2)*W+(V-1)*(3*V**2-1&
&0*V+4))+2*(V-1)*V*V4*W*(V**2*W**2-2*V*W+1)*(2*(V-1)*V**3*W**3-4*V**4*W&
&**2+V*(2*V**3+2*V**2-2*V+7)*W-4*V**2+7*V-5)+2*(V-1)*V*V3*W*(V**2*W**2-&
&2*V*W+1)*(2*V**3*W**3-2*(V-1)*V**2*W**2+V*(2*V**2+1)*W-2*V**3+6*V**2-9&
&*V+3))/((V-1)**2*V**3*W**2*(V*W-1)**3)
LV1 = -LOG(1-V)*(2*V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**5*(V**3+7*V**2&
&-3*V+1)*W**5-V**3*(9*V**4-9*V**3+35*V**2-31*V+10)*W**4+V**2*(23*V**4-6&
&0*V**3+130*V**2-105*V+30)*W**3-V*(28*V**4-87*V**3+158*V**2-119*V+30)*W&
&**2+(15*V**4-48*V**3+77*V**2-52*V+10)*W-3*(V-1)*(V**2-2*V+2))+2*V2*(V*&
&*2*W**2-2*V*W+1)*(2*V**5*W**4-3*(V-1)**2*V**3*W**3+V*(6*V**4-11*V**3+2&
&0*V**2-21*V+8)*W**2-(V-1)**2*(2*V**3+V**2-4*V+8)*W+2*(V-1)**2*(V**2-2*&
&V+2))+2*(V-1)*V*V4*W*(V*W-1)*(2*V**5*W**4-2*(V-1)*V**3*(2*V+1)*W**3+V*&
&*2*(2*V**3+V-14)*W**2-V*(2*V**3-6*V**2+11*V-22)*W-2*(2*V**2-4*V+5))+CQ&
&*V**2*V1*W*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+&
&2*CQ*V**3*V2*W**2*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)+2*(V-&
&1)*V*V3*W*(V*(2*V**2-3*V+4)*W-2*(V-1)*(V**2-2*V+2))*(V**2*W**2-2*V*W+1&
&))/((V-1)**2*V**3*W**2*(V*W-1)**3)
LV = LOG(V)*(V1*(-2*V**8*W**7+2*V**7*(V+1)*W**6-V**3*(2*V**5+9*V**4+V*&
&*3-8*V**2+6*V-2)*W**5+V**2*(2*V**6-V**5+35*V**4-4*V**3-6*V**2+12*V-6)*&
&W**4-V*(2*V**7-4*V**6+13*V**5+23*V**4+22*V**2-6)*W**3+(6*V**7-24*V**6+&
&53*V**5-41*V**4+26*V**3+22*V**2-12*V-2)*W**2-(6*V**6-28*V**5+65*V**4-7&
&6*V**3+49*V**2-6*V-6)*W+2*(V-1)*(V**2-2*V+2)**2)-2*V2*(2*V**7*W**6-V**&
&3*(9*V**4-8*V**3+7*V**2-6*V+2)*W**5+V**2*(5*V**5-3*V**4+25*V**3-15*V**&
&2-6*V+6)*W**4+V*(2*V**6-18*V**5+55*V**4-126*V**3+91*V**2-18*V-6)*W**3-&
&(6*V**6-33*V**5+91*V**4-155*V**3+113*V**2-30*V-2)*W**2+2*(V-1)**2*(3*V&
&**3-7*V**2+12*V-6)*W-2*(V-1)**2*(V**2-2*V+2))+CQ*V**2*V1*W*(2*V**2*W**&
&2-2*V*(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+4*(V-1)*V*V3*W*(V**2*W&
&**2-2*V*W+1)*(V**3*W**2-4*V*W-(V-2)*(V**2-2*V+2))-4*(V-1)*V*V4*W*(V**2&
&*W**2-2*V*W+1)*(V*(V**2+1)*W**2-(5*V+1)*W-V**3+4*V**2-6*V+5)+2*CQ*V**3&
&*V2*W**2*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V**&
&3*W**2*(V*W-1)**3)
LVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(2*V**5*W**4-2&
&*(V-2)*V**4*W**3+V**3*(V**2+2*V+3)*W**2-(3*V**4-8*V**3+17*V**2-14*V+4)&
&*W+2*(V-1)*(V**2-2*V+2))+V*V2*W*(4*V**3*W**2-(V-5)*V**2*W-2*(V-1)**2)+&
&(V-1)*V*V4*W*(4*V**2*W-4*V**2+7*V-7)+(V-1)*V*V3*W*(2*V**2*W-2*V**2+3*V&
&-5))*LOG(1-V*W)/((V-1)**2*V**3*W**2*(V*W-1)**4)
LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(2*V**5*W**4-&
&2*(V-2)*V**4*W**3+V**3*(2*V**2-3*V+5)*W**2-(V-1)**2*(2*V**3-V**2-4*V+4&
&)*W+2*(V-1)**3*(V**2-2*V+2))+4*(V-1)*V*V4*W*(V**3*W**2-V**2*(2*V-3)*W+&
&(V-1)*(V**2-2*V+2))-(V-1)*V*V2*W*(V**2*W-(V-4)*(V-1)))*LOG(V*W-V+1)/((&
&V-1)**2*V**3*W**2*(V*W-1)**4)
LW = -(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V2*(2*V**5*W**4-(V-1)*V**3*(6&
&*V-1)*W**3+V**2*(3*V**3+2*V**2+4*V+1)*W**2-2*(V-1)**2*V*(V**2-2*V+2)*W&
&+2*(V-1)**2*(V**2-2*V+2))+2*V1*(V*(V**5+4*V**4-3*V**3+3*V**2-3*V+1)*W*&
&*3-(2*V**6-2*V**5+2*V**4+17*V**3-12*V**2+2*V+1)*W**2+(2*V**6-8*V**5+20&
&*V**4-30*V**3+37*V**2-23*V+5)*W-2*(V-1)*(V**2-2*V+2)**2)+2*(V-1)*V*V4*&
&W*(2*(V-1)*V**3*W**3-2*V*(2*V**3-3*V**2-1)*W**2+(2*V**4-2*V**3+V**2-10&
&*V-2)*W-2*(V**3-4*V**2+6*V-6))+2*(V-1)*V**2*V3*W**2*(2*V**2*W**2-2*(V-&
&1)*V*W+V+4)-2*CQ*(V-1)**2*(V**2-2*V+2)*V2*(V*W-1)*(W**2-2*W+2)+CQ*(V-1&
&)*(V**2-2*V+2)**2*V1*(V*W-1)*(W**2-2*W+2))*LOG(W)/((V-1)**2*V**3*W**2*&
&(V*W-1)**4)
CVC = (3*(V-1)*V1*(V*W-1)*(6*V**9*W**8-4*V**8*(3*V+1)*W**7+V**4*(8*V**&
&5+V**4+19*V**3-62*V**2+64*V-24)*W**6-(V-2)*V**3*(8*V**5-22*V**4+75*V**&
&3-35*V**2-22*V+24)*W**5+V**3*(12*V**6-43*V**5+55*V**4+62*V**3-344*V**2&
&+368*V-144)*W**4-V*(6*V**8+4*V**7-93*V**6+245*V**5-304*V**4+90*V**3+20&
&0*V**2-200*V+48)*W**3+(18*V**8-60*V**7-3*V**6+289*V**5-551*V**4+503*V*&
&*3-198*V**2-16*V+24)*W**2-(V-1)*(18*V**6-66*V**5+59*V**4+56*V**3-119*V&
&**2+100*V-36)*W+2*(V-1)**2*(3*V**4-10*V**3+8*V**2+2))+12*(V-1)*V2*(V*W&
&-1)*(2*V**8*W**7-V**4*(10*V**4-23*V**3+35*V**2-24*V+6)*W**6+V**3*(12*V&
&**5-25*V**4+13*V**3+28*V**2-36*V+12)*W**5-V**3*(10*V**5-25*V**4-11*V**&
&3+98*V**2-98*V+30)*W**4+V*(2*V**7+3*V**6-27*V**5+24*V**4+62*V**3-116*V&
&**2+66*V-12)*W**3-(V-1)*(6*V**6-10*V**5-3*V**4+22*V**3-7*V**2-12*V+6)*&
&W**2+(V-1)**2*(V**2-2*V+2)*(6*V**2-V-3)*W-2*(V-1)**3*(V**2-2*V+2))-6*A&
&L*(V-1)*V1*(V*W-1)*(V*W-V+1)*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**3*(2*V*&
&*5-5*V**4+17*V**3-16*V**2+12*V-4)*W**5-V**2*(2*V**6-7*V**5+13*V**4+8*V&
&**3-24*V**2+28*V-12)*W**4+V*(2*V**7-4*V**6-V**5+25*V**4-32*V**3+16*V**&
&2+12*V-12)*W**3-(6*V**7-24*V**6+43*V**5-27*V**4-10*V**3+32*V**2-12*V-4&
&)*W**2+(6*V**6-28*V**5+63*V**4-72*V**3+41*V**2-8)*W-2*(V-1)*(V**2-2*V+&
&2)**2)-12*AL*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**7*W**6-V**3*(3*V**4-2*V*&
&*3+7*V**2-6*V+2)*W**5+V**2*(3*V**5-5*V**4+5*V**3+9*V**2-14*V+6)*W**4-V&
&*(2*V**6-2*V**5-5*V**4+20*V**3-11*V**2-6*V+6)*W**3+(6*V**6-17*V**5+19*&
&V**4+3*V**3-17*V**2+6*V+2)*W**2-2*(V-1)**2*(3*V+1)*(V**2-2*V+2)*W+2*(V&
&-1)**2*(V**2-2*V+2))-8*CQ*(V-1)**2*V1*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*&
&(V**6*W**4+V**2*(V**4-4*V**3+10*V**2-8*V+4)*W**2-2*V*(V**2-2*V+2)**2*W&
&+V**4-4*V**3+9*V**2-8*V+4)+12*(V-1)**2*V*V4*W*(V*W-V+1)*(V**2*W**2-2*V&
&*W+1)*(V**5*W**4-V**2*(2*V**3+2*V**2+V+4)*W**3+V*(V**4+4*V**3-6*V**2+7&
&*V+8)*W**2-(2*V**4+3*V**3-16*V**2+11*V+4)*W+V**3+5*V**2-13*V+5)-16*CQ*&
&(V-1)**2*V2*W*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**5*W**3-(V-1)*V**2*(V**&
&2-2*V+2)*W**2+V*(2*V**3-5*V**2+8*V-4)*W-(V-1)*(V**2-2*V+2))-12*(V-1)**&
&3*V*V3*W**2*(V*W-V+1)*(V**2*(4*V-3)*W**2-2*V*(2*V**2+4*V-3)*W+V**2+4*V&
&-3)*(V**2*W**2-2*V*W+1))/((V-1)**3*V**3*W**2*(V*W-1)**4*(V*W-V+1))/6.D&
&0
LM = LOG(S/M**2)*(-V1*(2*V**8*W**7-2*V**7*(V+1)*W**6+V**3*(2*V**5-V**4&
&+11*V**3-8*V**2+6*V-2)*W**5-V**2*(2*V-1)*(V**5-V**4+8*V**3+10*V**2-4*V&
&+6)*W**4+V*(2*V**7-4*V**6+9*V**5+11*V**4+20*V**3-18*V**2+12*V-6)*W**3-&
&(6*V**7-24*V**6+57*V**5-63*V**4+68*V**3-26*V**2-2)*W**2+(6*V**6-28*V**&
&5+71*V**4-96*V**3+83*V**2-30*V-2)*W-2*(V-1)*(V**2-2*V+2)*(V**2-2*V+3))&
&-2*V2*(2*V**7*W**6-V**3*(3*V**4-2*V**3+7*V**2-6*V+2)*W**5+V**2*(3*V**5&
&-5*V**4+5*V**3+9*V**2-14*V+6)*W**4-V*(2*V**6-2*V**5-5*V**4+20*V**3-11*&
&V**2-6*V+6)*W**3+(6*V**6-17*V**5+19*V**4+3*V**3-17*V**2+6*V+2)*W**2-2*&
&(V-1)**2*(3*V+1)*(V**2-2*V+2)*W+2*(V-1)**2*(V**2-2*V+2))-4*(V-1)*V*V4*&
&W*(V*W-1)*(V**2*W**3-V*(3*V+2)*W**2-(V**2-6*V-1)*W-V**2+2*V-3))/((V-1)&
&**2*V**3*W**2*(V*W-1)**3)
LMP = -2*LOG(S/MP**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*&
&(V**6*W**5-V**5*(2*V-3)*W**4+V**4*(2*V**2-5*V+5)*W**3-(V-1)*V*(V+1)*(2&
&*V-1)*(V**2-2*V+2)*W**2+(V-1)**3*(2*V**3-3*V**2+2)*W-(V-1)**4*(V**2-2*&
&V+2))+2*(V-1)*V*V4*W*(V*W-V+1)*(V**3*W**2-V**2*(2*V-3)*W+(V-1)*(V**2-2&
&*V+2)))/((V-1)**2*V**3*W**2*(V*W-1)**4*(V*W-V+1))
STRUV6=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV7(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV7
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+V**4*(16*V**4+&
&24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V**2-V+2)*W**5+2*V**&
&2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*V**2*(6*V**6+2*V**5&
&+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+12*V**7-17*V**6+54*V**5-5&
&*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)*(4*V**6-2*V**5+7*V**4+6*V**3+6*&
&V**2-4*V-1)*W+4*(V-1)**2*(V**2+1)*(V**2-V+2))-V1*(-4*V**8*W**8+4*V**7*&
&(3*V+1)*W**7-V**4*(16*V**4+2*V**3+17*V**2+1)*W**6+2*V**4*(8*V**4-11*V*&
&*3+25*V**2+4*V+2)*W**5-2*V**2*(8*V**6-19*V**5+20*V**4+33*V**3-10*V**2+&
&V-1)*W**4+2*V**2*(6*V**6-9*V**5-11*V**4+63*V**3-16*V**2-6*V-3)*W**3-(4&
&*V**8+8*V**7-49*V**6+88*V**5+13*V**4-46*V**3-V**2-2*V+1)*W**2+2*(V-1)*&
&(4*V**6-6*V**5+4*V**4+11*V**3+9*V**2-5*V-1)*W-4*(V-1)**2*(V**2+1)*(V**&
&2-2*V+3))+2*V2*(V*W-1)*(V*W-V+1)*(V**3*(6*V**2-3*V+1)*W**4-V**3*(8*V**&
&2+9*V-5)*W**3+V*(6*V**4+9*V**3+15*V**2-13*V-1)*W**2-V*(9*V**3+V**2+11*&
&V-13)*W+4*(V-1)*(V**2+1))+2*CQ*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)&
&*V*W+(V-1)**2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2&
&+1))-4*CQ*(V-1)*V*V4*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1&
&)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*V4*W*(V*W-1)*(V*W-V-1)*(V**2&
&*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*V*V3*W&
&*(V*W-1)*(V*W-V+1)*(V*W**2-W+V-1)*(V**2*W**2-2*V**2*W+V**2+1))/((V-1)*&
&*2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LV1 = LOG(1-V)*(2*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**4*W**4-2*V**3&
&*(V+2)*W**3+V**2*(V**2+2*V+9)*W**2-2*V*(V**2+2*V+3)*W+2*(V+1)*(V**2+1)&
&)+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4*V**3*W&
&**3+8*V**2*W**2-8*V*W+4)-4*V2*(V*W-V-1)*(V**2*W**2-2*V*W+1)*(V**3*W**3&
&-V**2*(V+1)*W**2+V*(3*V**2-2*V+1)*W-(V-1)*(V**2+1))-2*CQ*V2*(V*W-1)*(V&
&*W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)+4*(V-1)*V&
&**2*V4*W**2*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)**2*V&
&3*W*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*V**2*W+V**2+1))/((V-1)**2*V**2*W**2&
&*(V*W-1)**2*(V*W-V+1))
LV = LOG(V)*(-CQ*V1*(4*V**8*W**8-4*V**7*(3*V+2)*W**7+3*V**6*(V+1)*(5*V&
&+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(9*V**4+6*V**3+46*V*&
&*2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2-2*V-1)*W**3+2*V**2&
&*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W**2-2*(V-1)*V*(2*V**5-3*V&
&**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*(V**2+1)*(V**2-2*V+3))-V1*(-4*V*&
&*8*W**8+4*V**7*(3*V+2)*W**7-V**4*(16*V**4+12*V**3+23*V**2+1)*W**6+2*V*&
&*4*(8*V**4-8*V**3+37*V**2+3*V+2)*W**5-2*V**2*(8*V**6-22*V**5+35*V**4+3&
&7*V**3-14*V**2+V-1)*W**4+2*V**2*(6*V**6-14*V**5-5*V**4+74*V**3-20*V**2&
&-10*V-3)*W**3-(4*V**8+4*V**7-57*V**6+114*V**5+7*V**4-52*V**3-3*V**2-2*&
&V+1)*W**2+2*(V-1)*(4*V**6-10*V**5+7*V**4+10*V**3+14*V**2-8*V-1)*W-4*(V&
&-1)**2*(V**2+1)*(V**2-3*V+4))-2*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3&
&*(11*V**2-2*V+1)*W**4+2*V**3*(4*V**2+11*V-5)*W**3-V*(V**4+16*V**3+24*V&
&**2-20*V-1)*W**2-2*V*(V**4-4*V**3-2*V**2-8*V+9)*W+2*(V-3)*(V-1)*(V**2+&
&1))+2*CQ*V2*(V*W-1)*(V**2*W**2-2*V*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**&
&2)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-4*CQ*(V-1)*V*V4*W*(V**2*W**2-2*V*W&
&+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-&
&1)*V*V4*W*(V*W-2)*(V*W-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-&
&2*V**2*W+V**2+1)+4*(V-1)*V*V3*W*(V*W-1)*(V*W-V+1)*(V*W+V-1)*(V**2*W**2&
&-2*V**2*W+V**2+1))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LVW = -2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V&
&**4*W**4-2*V**3*(V+3)*W**3+V**2*(V**2+3*V+12)*W**2-V*(V+3)**2*W+4*(V**&
&2+1))+V2*(-4*V**3*W**3+V**2*(5*V+7)*W**2-V*(3*V**2+6*V+7)*W+4*(V**2+1)&
&))*LOG(1-V*W)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LTVW = -4*(V**2*W**2-2*V*W+1)*(CQ*V1*(V**2*W**2-2*V**2*W+V**2+1)*(V**4&
&*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)-(V-1)*&
&V*V1*W*(V**3*W**3-V**2*(5*V-3)*W**2+(V-1)*V*(5*V+1)*W-(V-1)**2*(V+1))+&
&2*CQ*(V-1)*V*V4*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W&
&+V**2+1))*LOG(V*W-V+1)/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LW = (-2*V1*(V*W-1)*(V*W-V+1)*(2*V**3*W**3-V**2*(V**2+7)*W**2+2*V*(V**&
&3-V**2+3*V+3)*W-2*(V**2+1)*(V**2-V+2))-4*V2*(V*W-1)*(V*W-V-1)*(V**3*W*&
&*3-V**2*(3*V-1)*W**2+V*(V**2+2*V-1)*W-(V-1)*(V**2+1))+4*(V-1)**2*V4*W*&
&(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V**2*W+V**2+1)-4*(V-1)*V*V3*W**2*(V*W-V&
&-1)*(V**2*W**2-2*V**2*W+V**2+1)+2*CQ*V*(V**2+1)*V2*(V*W-1)*(V*W-V+1)*(&
&W**2-2*W+2)+CQ*(V**2+1)**2*V1*(V*W-1)*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/(&
&(V-1)**2*V**2*W**2*(V*W-1)*(V*W-V+1))
CVC = (3*V1*(8*V**7*W**7+V**4*(2*V**4-24*V**3-13*V**2+5*V-2)*W**6-V**4&
&*(4*V**4-32*V**3-46*V**2+17*V-7)*W**5+2*V**2*(V**6-12*V**5-36*V**4+7*V&
&**3+5*V**2-7*V+2)*W**4+2*V**2*(8*V**5+19*V**4+22*V**3-27*V**2+11*V-5)*&
&W**3-(8*V**7+17*V**6+3*V**5+40*V**4-66*V**3+21*V**2-9*V+2)*W**2+(V-1)*&
&(16*V**5-11*V**4+34*V**3+2*V**2-6*V-3)*W-8*(V-1)**2*V*(V**2-V+2))+3*AL&
&*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(V+5)*W**5+V*&
&*2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V+1&
&)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+V*&
&*4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))+12*V*V2*(V*W-1)*(V**5*&
&W**6-V**3*(4*V**2-V+1)*W**5+V**2*(5*V**3+2*V**2+3*V-1)*W**4-V*(5*V**4-&
&2*V**3+16*V**2-8*V-1)*W**3+(4*V**5-3*V**4+12*V**3-4*V**2-6*V+1)*W**2-(&
&V-1)*(V**4+2*V**3+4*V**2+4*V-3)*W+(V-1)**2*(V**2+3))+4*CQ*V1*(V*W-1)*(&
&V*W-V+1)*(V**5*(V+2)*W**5-V**4*(V**2+5*V-3)*W**4+V**2*(V**4+7*V**3-V**&
&2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2-5*V+1)*W**2+(V-1)*(4*V**4-6*V**&
&3+13*V**2-4*V+1)*W-(V-1)**2*(3*V**2-6*V+7))+12*(V-1)*V*V4*W*(V**2*W**2&
&-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3&
&*V*W+(V-1)**4)-6*AL*V2*(V*W-1)*(V*W-V-1)*(V**2*W**2-2*(V-1)*V*W+(V-1)*&
&*2)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V**2+1))+24*CQ*&
&(V-1)*V*V4*W*(V**2*W**2-2*V*W+1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)*&
&*2*V*W-(V-1)**3)-8*CQ*V2*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V&
&-1)**2)*((V-1)*V**2*W**2-V*(V**2+2*V-1)*W+2*(V-1)))/((V-1)**2*V**2*W**&
&2*(V*W-1)**2*(V*W-V+1)**2)/3.D0
LM = LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*V&
&**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10&
&*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*&
&W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))-2*V2*&
&(V*W-1)*(V*W-V-1)*(2*V**3*W**3-V*(V**2+4*V+1)*W**2+2*V*(V**2+3)*W-2*(V&
&**2+1)))/((V-1)**2*V**2*W**2*(V*W-1)**2*(V*W-V+1)**2)
LMP = 2*LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+V**2+1)*(&
&V1*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4&
&)+2*(V-1)*V*V4*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2))/((V-1)**2*V**2*W**2&
&*(V*W-1)**2*(V*W-V+1)**2)
STRUV7=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV8(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV8
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*&
&W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V**4-59*V**3&
&+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4+25*V**3-16*V**2+5&
&*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V**2+20*V-3)*W*&
&*3+(16*V**8+40*V**7-204*V**6+280*V**5-150*V**4+12*V**3+15*V**2-6*V+1)*&
&W**2-2*(V-1)*(16*V**6-28*V**5+6*V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)&
&**2*(2*V**2-3*V+2)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**&
&7-V**4*(36*V**4-30*V**3+13*V**2-4*V+1)*W**6+2*V**4*(28*V**4-39*V**3+22&
&*V**2-5*V+2)*W**5-2*V**2*(32*V**6-46*V**5+V**4+37*V**3-20*V**2+5*V-1)*&
&W**4+2*V**2*(24*V**6-20*V**5-66*V**4+130*V**3-79*V**2+20*V-3)*W**3-(16&
&*V**8+32*V**7-204*V**6+292*V**5-138*V**4-12*V**3+23*V**2-6*V+1)*W**2+2&
&*(V-1)*(16*V**6-32*V**5+4*V**4+46*V**3-43*V**2+14*V-1)*W-4*(V-1)**2*(2&
&*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(V**3*(4*V**&
&2-14*V+1)*W**4-4*V**3*(3*V**2-11*V+4)*W**3+V*(16*V**4-50*V**3+24*V**2+&
&2*V-1)*W**2-2*(V-1)*V*(4*V**3-4*V**2-5*V+2)*W+4*(V-1)**2*(2*V**2-2*V+1&
&))-2*CQ*(V-1)*V2*(V*W-2*V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**&
&3-V*(6*V**2-6*V+1)*W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**2-2*V+1)))/&
&((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)&
&*(CQ*V1*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W**3+V**2*W*&
&*2-V*W+1)+2*V*V1*(W-1)*(V*W-1)*(V**4*W**4-V**3*(2*V-1)*W**3+V**2*(2*V*&
&*2-3*V+2)*W**2-V**2*W+2*V**2-2*V+1)-2*(V-1)*V2*(V*W-1)*(V**4*W**4-2*V*&
&*3*(V+1)*W**3+V**3*(2*V+3)*W**2-V*(8*V**2-9*V+5)*W+2*V**2-2*V+1))/((V-&
&1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*W**&
&7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+1)*W**5+V**&
&4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**2+&
&6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V**2+3)*W**2-2*(&
&V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V+5)*W+2*(V-1)**2*(2*V**2-4*V&
&+3)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4&
&-50*V**3+17*V**2-4*V+1)*W**6+2*V**4*(42*V**4-59*V**3+24*V**2-V+2)*W**5&
&-2*V**2*(44*V**6-64*V**5-9*V**4+59*V**3-26*V**2+5*V-1)*W**4+2*V**2*(28&
&*V**6-22*V**5-98*V**4+180*V**3-101*V**2+22*V-3)*W**3-(16*V**8+40*V**7-&
&260*V**6+380*V**5-178*V**4-16*V**3+27*V**2-6*V+1)*W**2+2*(V-1)*(16*V**&
&6-36*V**5+2*V**4+62*V**3-57*V**2+18*V-1)*W-4*(V-1)**2*(2*V**2-5*V+4)*(&
&2*V**2-2*V+1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3*(10*V**2&
&-18*V+1)*W**4+4*V**3*(5*V**2-14*V+5)*W**3-V*(20*V**4-60*V**3+30*V**2-1&
&)*W**2+2*V*(4*V**4-8*V**3-4*V**2+10*V-3)*W-2*(V-1)*(2*V-3)*(2*V**2-2*V&
&+1))-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*&
&W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V**2*W*&
&*2*(V*W-1)**4*(V*W-V+1)**2)
LVW = 4*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2&
&*W**2-4*V*W+1)*(V1*(V**4*W**4-V**3*(2*V-1)*W**3+V**2*(2*V**2-4*V+3)*W*&
&*2+V*(2*V**2-4*V+1)*W+2*V**2-2*V+1)-3*(V-1)*V*V2*W*(V*W-2*V+1))*LOG(1-&
&V*W)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1*(2*V**2*W*&
&*2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+8*(V-1)**2&
&*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W*&
&*2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)-2*&
&(V-1)*V*V2*W*(V*W-V+1)*(V**2*W**2-V*(2*V-1)*W+(V-1)**2)+(V-1)**2*V**2*&
&V1*W*((2*V-1)*W-2*(V-1)))*LOG(V*W-V+1)/((V-1)*V**2*W**2*(V*W-1)**4*(V*&
&W-V+1)**2)
LW = (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*&
&(2*V1*(2*(V-1)*V**4*W**4-V**3*(2*V-1)*(4*V-3)*W**3+V**2*(12*V**3-16*V*&
&*2+4*V+3)*W**2-2*V*(V+1)*(4*V**3-8*V**2+6*V-1)*W+2*(2*V**2-3*V+2)*(2*V&
&**2-2*V+1))+4*(V-1)*V2*(V*W-2*V+1)*(V**3*W**3-2*(V-1)*V**2*W**2+V*(2*V&
&**2-1)*W-2*V**2+2*V-1)-2*CQ*(V-1)*V*(2*V**2-2*V+1)*V2*(V*W-1)*(W**2-2*&
&W+2)-CQ*(2*V**2-2*V+1)**2*V1*(V*W-1)*(W**2-2*W+2))*LOG(W)/((V-1)*V**2*&
&W**2*(V*W-1)**4*(V*W-V+1)**2)
CVC = (3*AL*V1*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2&
&*V**6*W**6-2*V**5*(2*V+1)*W**5+V**2*(8*V**3-4*V**2+4*V-1)*W**4+2*V*(2*&
&V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W**3-(8*V**6-16*V**4+16*V**3-10*V**2+1&
&)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V**2+2*V-1)*W-4*(V-1)*V*(2*V**2-2*V+1)&
&)-6*(V-1)*V*V2*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**5*W**6-V**3*(2*V-1)&
&*(4*V+1)*W**5+V**2*(18*V**3-10*V**2-1)*W**4-V*(20*V**4-4*V**3-10*V**2+&
&1)*W**3+(8*V**5+20*V**4-40*V**3+18*V**2-2*V+1)*W**2-2*(8*V**4-10*V**3-&
&2*V**2+8*V-3)*W+2*(V-1)*(4*V**2-6*V+3))-3*V1*(V**4*W**4-4*V**3*W**3+6*&
&V**2*W**2-4*V*W+1)*(8*(V-1)*V**5*W**5-V**2*(32*V**4-59*V**3+28*V**2-5*&
&V+2)*W**4+(V-1)*V*(64*V**4-101*V**3+48*V**2-13*V+4)*W**3-(V-1)*(80*V**&
&5-180*V**4+156*V**3-59*V**2+13*V-2)*W**2+(V-1)**2*(56*V**4-106*V**3+84&
&*V**2-21*V+3)*W-8*(V-1)**3*V*(2*V**2-3*V+2))-6*AL*(V-1)*V*(2*V**2-2*V+&
&1)*V2*(W**2-2*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*&
&W**3+6*V**2*W**2-4*V*W+1)-8*CQ*(V-1)*V*(2*V**2-2*V+1)*V2*W*(V**2*W**2-&
&2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)-4*CQ&
&*(2*V**2-2*V+1)**2*V1*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*&
&V**3*W**3+6*V**2*W**2-4*V*W+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)*&
&*2)/3.D0
LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*&
&V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2*V*(4*V**5-&
&2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V**2+1)*W**2-2*(4*V&
&**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**2-2*V+1))+2*(V-1)*&
&V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*V*W+1))/((V-1)*V**2*W**2&
&*(V*W-1)**2*(V*W-V+1)**2)
LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)&
&**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V1*(V**2*W**2-2*(V-1)*V&
&*W+2*(V-1)**2)-2*(V-1)*V2*(V*W-V+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-&
&V+1)**2)
STRUV8=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV9(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV9
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*&
&W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V**4-59*V**3&
&+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4+25*V**3-16*V**2+5&
&*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V**2+20*V-3)*W*&
&*3+(16*V**8+40*V**7-204*V**6+280*V**5-150*V**4+12*V**3+15*V**2-6*V+1)*&
&W**2-2*(V-1)*(16*V**6-28*V**5+6*V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)&
&**2*(2*V**2-3*V+2)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**&
&7-V**4*(36*V**4-40*V**3+23*V**2-4*V+1)*W**6+2*V**4*(28*V**4-59*V**3+49&
&*V**2-12*V+2)*W**5-2*V**2*(32*V**6-79*V**5+54*V**4+17*V**3-20*V**2+5*V&
&-1)*W**4+2*V**2*(24*V**6-46*V**5-23*V**4+121*V**3-91*V**2+24*V-3)*W**3&
&-(16*V**8+16*V**7-196*V**6+354*V**5-232*V**4+38*V**3+13*V**2-6*V+1)*W*&
&*2+2*(V-1)*(16*V**6-40*V**5+26*V**4+23*V**3-31*V**2+11*V-1)*W-4*(V-1)*&
&*2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(V**3*(&
&4*V**2+V+1)*W**4-V**3*(12*V**2+V-5)*W**3+V*(16*V**4+4*V**3-30*V**2+17*&
&V-1)*W**2-(V-1)*V*(8*V**3+16*V**2-28*V+13)*W+4*(V-1)**2*(2*V**2-2*V+1)&
&)-2*CQ*(V-1)*V2*(V*W-2*V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**3&
&-V*(6*V**2-6*V+1)*W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**2-2*V+1)))/(&
&(V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)&
&*(CQ*V1*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W**3+V**2*W*&
&*2-V*W+1)+V*V1*(V**2*W**2-2*V*W+1)*(2*V**3*W**4-2*V**2*(3*V-2)*W**3+V*&
&(8*V**2-11*V+5)*W**2-(2*V-1)*(2*V**2-3*V+3)*W+2*(2*V**2-2*V+1))-(V-1)*&
&V2*(V*W-1)*(V**2*W**2+2*V*W+1)*(2*V**2*W**2-V*(4*V-1)*W+2*(2*V**2-2*V+&
&1)))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*W**&
&7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+1)*W**5+V**&
&4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**2+&
&6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V**2+3)*W**2-2*(&
&V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V+5)*W+2*(V-1)**2*(2*V**2-4*V&
&+3)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4&
&-62*V**3+29*V**2-4*V+1)*W**6+2*V**4*(42*V**4-83*V**3+58*V**2-11*V+2)*W&
&**5-2*V**2*(44*V**6-102*V**5+57*V**4+29*V**3-24*V**2+5*V-1)*W**4+2*V**&
&2*(28*V**6-50*V**5-48*V**4+166*V**3-115*V**2+28*V-3)*W**3-(16*V**8+24*&
&V**7-252*V**6+452*V**5-298*V**4+56*V**3+11*V**2-6*V+1)*W**2+2*(V-1)*(1&
&6*V**6-44*V**5+26*V**4+34*V**3-41*V**2+14*V-1)*W-4*(V-1)**2*(2*V**2-5*&
&V+4)*(2*V**2-2*V+1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3*(1&
&0*V**2+1)*W**4+2*V**3*(10*V**2-V-5)*W**3-V*(20*V**4-42*V**2+24*V-1)*W*&
&*2+2*V*(4*V**4+4*V**3-28*V**2+28*V-9)*W-2*(V-1)*(2*V-3)*(2*V**2-2*V+1)&
&)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2&
&*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*&
&(V*W-1)**4*(V*W-V+1)**2)
LVW = 4*V1*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-V**3*(2*V-1)*W*&
&*3+V**2*(2*V**2-5*V+4)*W**2+V*(4*V**2-7*V+2)*W+2*V**2-2*V+1)*(V**4*W**&
&4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*LOG(1-V*W)/((V-1)*V**2*W**2*(V*W-1)&
&**4*(V*W-V+1)**2)
LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1*(2*V**2*W*&
&*2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+8*(V-1)**2&
&*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-(V-1)*V*V1*W*(V**3*W**3-V**2*(4*&
&V-3)*W**2+(V-1)*V*(3*V-2)*W+(V-1)**2)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W*&
&*2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)+(V&
&-1)*V*V2*W*(V*W-V+1)*(V**2*W**2-V*(5*V-4)*W+(V-1)*(4*V-1)))*LOG(V*W-V+&
&1)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LW = (V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*&
&(2*V1*(V*W-1)*(2*(V-1)*V**3*W**3-V**2*(8*V**2-14*V+7)*W**2+2*V*(6*V**3&
&-14*V**2+12*V-3)*W-2*(2*V**2-3*V+2)*(2*V**2-2*V+1))+4*(V-1)*V2*(V*W-2*&
&V+1)*(V**3*W**3-V**2*(2*V+1)*W**2+V*(2*V**2-1)*W-2*V**2+2*V-1)-2*CQ*(V&
&-1)*V*(2*V**2-2*V+1)*V2*(V*W-1)*(W**2-2*W+2)-CQ*(2*V**2-2*V+1)**2*V1*(&
&V*W-1)*(W**2-2*W+2))*LOG(W)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
CVC = (-3*V1*(V**2*W**2-2*V*W+1)*(8*(V-1)*V**7*W**7-V**4*(32*V**4-43*V&
&**3+10*V**2-3*V+2)*W**6+V**4*(64*V**4-101*V**3+37*V**2-11*V+7)*W**5-2*&
&V**2*(40*V**6-66*V**5+25*V**4-3*V**3+5*V-2)*W**4+2*(V-1)*V**2*(28*V**5&
&-V**4-43*V**3+33*V**2-14*V+5)*W**3-(V-1)*(16*V**7+56*V**6-156*V**5+124&
&*V**4-26*V**3-9*V**2+5*V-2)*W**2+(V-1)**2*(32*V**5-24*V**4-26*V**3+52*&
&V**2-21*V+3)*W-8*(V-1)**3*V*(2*V**2-3*V+2))+3*AL*V1*(V**2*W**2-2*V*W+1&
&)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V*&
&*2*(8*V**3-4*V**2+4*V-1)*W**4+2*V*(2*V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W&
&**3-(8*V**6-16*V**4+16*V**3-10*V**2+1)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V&
&**2+2*V-1)*W-4*(V-1)*V*(2*V**2-2*V+1))-12*(V-1)*V*V2*(V*W-V+1)*(V**2*W&
&**2-2*V*W+1)*(V**5*W**6-V**3*(4*V**2-V+1)*W**5+V**2*(9*V**3-5*V**2+1)*&
&W**4-V*(10*V**4-2*V**3-14*V**2+12*V-1)*W**3+(4*V**5+10*V**4-26*V**3+18&
&*V**2-V-1)*W**2-(8*V**4-10*V**3-2*V**2+8*V-3)*W+(V-1)*(4*V**2-6*V+3))-&
&6*AL*(V-1)*V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*(V-1)*V*W+(V-&
&1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)-8*CQ*(V-1)*V*(2*V**&
&2-2*V+1)*V2*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+&
&6*V**2*W**2-4*V*W+1)-4*CQ*(2*V**2-2*V+1)**2*V1*W*(V**2*W**2-2*(V-1)*V*&
&W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1))/((V-1)*V**2*W&
&**2*(V*W-1)**4*(V*W-V+1)**2)/3.D0
LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*&
&V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2*V*(4*V**5-&
&2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V**2+1)*W**2-2*(4*V&
&**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**2-2*V+1))+2*(V-1)*&
&V*(2*V**2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*V*W+1))/((V-1)*V**2*W**2&
&*(V*W-1)**2*(V*W-V+1)**2)
LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)&
&**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V1*(V**2*W**2-2*(V-1)*V&
&*W+2*(V-1)**2)-2*(V-1)*V2*(V*W-V+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-&
&V+1)**2)
STRUV9=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV10(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV10
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*&
&W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V**4-59*V**3&
&+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4+25*V**3-16*V**2+5&
&*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V**2+20*V-3)*W*&
&*3+(16*V**8+40*V**7-204*V**6+280*V**5-150*V**4+12*V**3+15*V**2-6*V+1)*&
&W**2-2*(V-1)*(16*V**6-28*V**5+6*V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-1)&
&**2*(2*V**2-3*V+2)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(4*V-1)*W**&
&7-V**4*(36*V**4-40*V**3+23*V**2-4*V+1)*W**6+2*V**4*(28*V**4-59*V**3+49&
&*V**2-12*V+2)*W**5-2*V**2*(32*V**6-79*V**5+54*V**4+17*V**3-20*V**2+5*V&
&-1)*W**4+2*V**2*(24*V**6-46*V**5-23*V**4+121*V**3-91*V**2+24*V-3)*W**3&
&-(16*V**8+16*V**7-196*V**6+354*V**5-232*V**4+38*V**3+13*V**2-6*V+1)*W*&
&*2+2*(V-1)*(16*V**6-40*V**5+26*V**4+23*V**3-31*V**2+11*V-1)*W-4*(V-1)*&
&*2*(2*V**2-4*V+3)*(2*V**2-2*V+1))-2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(V**3*(&
&4*V**2+V+1)*W**4-V**3*(12*V**2+V-5)*W**3+V*(16*V**4+4*V**3-30*V**2+17*&
&V-1)*W**2-(V-1)*V*(8*V**3+16*V**2-28*V+13)*W+4*(V-1)**2*(2*V**2-2*V+1)&
&)-2*CQ*(V-1)*V2*(V*W-2*V+1)*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**3&
&-V*(6*V**2-6*V+1)*W**2+2*V*(4*V**2-6*V+3)*W-2*(V-1)*(2*V**2-2*V+1))+4*&
&CQ*V*V4*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W&
&**2-2*V**2*W+2*V**2-2*V+1)-4*(V-1)*V4*W*(V*W-2*V+1)*(V*W-V+1)*(V**2*W*&
&*2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-4*(V-1)*V*V3*W*(V*W-1)*(&
&V*W-V+1)*(V*W**2-(V-1)*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1))/((V-1)*&
&V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LV1 = -2*LOG(1-V)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)&
&*(CQ*V1*(V*W-V+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W*&
&*3+V**2*W**2-V*W+1)+V*V1*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(2*V**3*W**4-2*&
&V**2*(3*V-2)*W**3+V*(8*V**2-11*V+5)*W**2-(2*V-1)*(2*V**2-3*V+3)*W+2*(2&
&*V**2-2*V+1))-(V-1)*V2*(V*W-1)*(V*W-V+1)*(V**2*W**2+2*V*W+1)*(2*V**2*W&
&**2-V*(4*V-1)*W+2*(2*V**2-2*V+1))+2*CQ*V*V4*W*(V*W-V+1)*(V**2*W**2-2*V&
&*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+2*V*V4*W*(V*W-V+1)*(V**2*W**2-&
&2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-2*(V-1)*V3*W*(V*W-1)*(V*(V+&
&1)*W-V+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)*&
&*4*(V*W-V+1)**3)
LV = LOG(V)*(V**2*W**2-2*V*W+1)*(CQ*V1*(4*V**8*W**8-4*V**7*(5*V-2)*W**&
&7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+1)*W**5+V**&
&4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**2+&
&6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V**2+3)*W**2-2*(&
&V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V+5)*W+2*(V-1)**2*(2*V**2-4*V&
&+3)*(2*V**2-2*V+1))+V1*(-4*V**8*W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4&
&-62*V**3+29*V**2-4*V+1)*W**6+2*V**4*(42*V**4-83*V**3+58*V**2-11*V+2)*W&
&**5-2*V**2*(44*V**6-102*V**5+57*V**4+29*V**3-24*V**2+5*V-1)*W**4+2*V**&
&2*(28*V**6-50*V**5-48*V**4+166*V**3-115*V**2+28*V-3)*W**3-(16*V**8+24*&
&V**7-252*V**6+452*V**5-298*V**4+56*V**3+11*V**2-6*V+1)*W**2+2*(V-1)*(1&
&6*V**6-44*V**5+26*V**4+34*V**3-41*V**2+14*V-1)*W-4*(V-1)**2*(2*V**2-5*&
&V+4)*(2*V**2-2*V+1))+2*(V-1)*V2*(V*W-1)*(V*W-V+1)*(2*V**5*W**5-V**3*(1&
&0*V**2+1)*W**4+2*V**3*(10*V**2-V-5)*W**3-V*(20*V**4-42*V**2+24*V-1)*W*&
&*2+2*V*(4*V**4+4*V**3-28*V**2+28*V-9)*W-2*(V-1)*(2*V-3)*(2*V**2-2*V+1)&
&)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2&
&*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)+4*CQ*V*V4*W*(V**2*&
&W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+2*V&
&**2-2*V+1)-4*V*V4*W*(V*W-2*(V-1))*(V*W-V+1)*(V**2*W**2-2*V*W+1)*(V**2*&
&W**2-2*V**2*W+2*V**2-2*V+1)-4*(V-1)*V*V3*W*(V*W-1)*(V*W+1)*(V*W-V+1)*(&
&V**2*W**2-2*V**2*W+2*V**2-2*V+1))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1&
&)**2)
LVW = 4*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W*&
&*4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(V**4*W**4-V**3*(2*V-1)*W**3+V&
&**2*(2*V**2-5*V+4)*W**2+V*(4*V**2-7*V+2)*W+2*V**2-2*V+1)+2*V*V4*W*(V**&
&2*W**2-2*V**2*W+2*V**2-2*V+1))*LOG(1-V*W)/((V-1)*V**2*W**2*(V*W-1)**4*&
&(V*W-V+1)**3)
LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*V1*(2*V**2*W*&
&*2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+8*(V-1)**2&
&*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-(V-1)*V*V1*W*(V**3*W**3-V**2*(4*&
&V-3)*W**2+(V-1)*V*(3*V-2)*W+(V-1)**2)-2*CQ*(V-1)*V2*(V*W-V+1)*(V**2*W*&
&*2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)+(V&
&-1)*V*V2*W*(V*W-V+1)*(V**2*W**2-V*(5*V-4)*W+(V-1)*(4*V-1)))*LOG(V*W-V+&
&1)/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
LW = -(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)&
&*(-2*V1*(V*W-1)*(V*W-V+1)*(2*(V-1)*V**3*W**3-V**2*(8*V**2-14*V+7)*W**2&
&+2*V*(6*V**3-14*V**2+12*V-3)*W-2*(2*V**2-3*V+2)*(2*V**2-2*V+1))-4*(V-1&
&)*V2*(V*W-2*V+1)*(V*W-V+1)*(V**3*W**3-V**2*(2*V+1)*W**2+V*(2*V**2-1)*W&
&-2*V**2+2*V-1)+4*V4*W*(V*W-1)*(V*W-V+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V&
&+1)-4*(V-1)*V*V3*W**2*(V*W-2*V+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+2*&
&CQ*(V-1)*V*(2*V**2-2*V+1)*V2*(V*W-1)*(V*W-V+1)*(W**2-2*W+2)+CQ*(2*V**2&
&-2*V+1)**2*V1*(V*W-1)*(V*W-V+1)*(W**2-2*W+2))*LOG(W)/((V-1)*V**2*W**2*&
&(V*W-1)**4*(V*W-V+1)**3)
CVC = (AL*V1*(6*V**11*W**11-30*V**11*W**10+(54*V**11+30*V**10-30*V**9+&
&12*V**8-3*V**7)*W**9+(-18*V**11-144*V**10+120*V**9-36*V**8+3*V**7+3*V*&
&*6)*W**8+(-84*V**11+258*V**10-96*V**9-108*V**8+129*V**7-60*V**6+9*V**5&
&)*W**7+(144*V**11-132*V**10-306*V**9+540*V**8-375*V**7+111*V**6+9*V**5&
&-9*V**4)*W**6+(-96*V**11-180*V**10+744*V**9-678*V**8+126*V**7+270*V**6&
&-264*V**5+90*V**4-9*V**3)*W**5+(24*V**11+264*V**10-480*V**9-120*V**8+8&
&34*V**7-900*V**6+474*V**5-78*V**4-27*V**3+9*V**2)*W**4+(-96*V**10-96*V&
&**9+816*V**8-1140*V**7+642*V**6+24*V**5-300*V**4+195*V**3-48*V**2+3*V)&
&*W**3+(144*V**9-336*V**8+60*V**7+540*V**6-774*V**5+528*V**4-177*V**3+3&
&*V**2+15*V-3)*W**2+(-96*V**8+384*V**7-612*V**6+480*V**5-150*V**4-60*V*&
&*3+84*V**2-36*V+6)*W+24*V**7-120*V**6+252*V**5-288*V**4+192*V**3-72*V*&
&*2+12*V)+V2*((12*V**10-12*V**11)*W**10+(72*V**11-84*V**10+24*V**9-12*V&
&**8)*W**9+(-216*V**11+276*V**10-60*V**9-12*V**8+12*V**7)*W**8+(384*V**&
&11-420*V**10-264*V**9+516*V**8-252*V**7+36*V**6)*W**7+(-396*V**11+72*V&
&**10+1416*V**9-1752*V**8+684*V**7+12*V**6-36*V**5)*W**6+(216*V**11+552&
&*V**10-2256*V**9+1764*V**8+516*V**7-1224*V**6+468*V**5-36*V**4)*W**5+(&
&-48*V**11-600*V**10+1188*V**9+840*V**8-3576*V**7+3132*V**6-936*V**5-36&
&*V**4+36*V**3)*W**4+(192*V**10+240*V**9-2256*V**8+3384*V**7-1428*V**6-&
&864*V**5+996*V**4-276*V**3+12*V**2)*W**3+(-288*V**9+720*V**8+108*V**7-&
&2160*V**6+2928*V**5-1680*V**4+360*V**3+24*V**2-12*V)*W**2+(192*V**8-84&
&0*V**7+1416*V**6-1032*V**5+84*V**4+348*V**3-204*V**2+36*V)*W-48*V**7+2&
&64*V**6-612*V**5+768*V**4-552*V**3+216*V**2-36*V)+V1*((24*V**10-24*V**&
&11)*W**10+(120*V**11-129*V**10+6*V**9-9*V**8+6*V**7)*W**9+(-288*V**11+&
&288*V**10+60*V**9-12*V**8-18*V**7-6*V**6)*W**8+(432*V**11-315*V**10-37&
&2*V**9+201*V**8-36*V**7+72*V**6-18*V**5)*W**7+(-408*V**11-54*V**10+120&
&0*V**9-888*V**8+318*V**7-180*V**6+18*V**5+18*V**4)*W**6+(216*V**11+594&
&*V**10-1902*V**9+1245*V**8-48*V**7-243*V**6+240*V**5-126*V**4+18*V**3)&
&*W**5+(-48*V**11-600*V**10+1056*V**9+600*V**8-2340*V**7+2040*V**6-912*&
&V**5+204*V**4+18*V**3-18*V**2)*W**4+(192*V**10+240*V**9-2052*V**8+3048&
&*V**7-1653*V**6-132*V**5+591*V**4-300*V**3+72*V**2-6*V)*W**3+(-288*V**&
&9+720*V**8-24*V**7-1794*V**6+2628*V**5-1704*V**4+510*V**3-36*V**2-18*V&
&+6)*W**2+(192*V**8-840*V**7+1458*V**6-1146*V**5+147*V**4+432*V**3-324*&
&V**2+90*V-9)*W-48*V**7+264*V**6-624*V**5+816*V**4-624*V**3+264*V**2-48&
&*V)+AL*V4*(12*V**10*W**10+(-60*V**10-12*V**9)*W**9+(132*V**10+72*V**9-&
&24*V**8)*W**8+(-156*V**10-204*V**9+96*V**8+24*V**7)*W**7+(96*V**10+336&
&*V**9-180*V**8-72*V**7)*W**6+(-24*V**10-288*V**9+60*V**8+252*V**7-72*V&
&**6)*W**5+(96*V**9+192*V**8-480*V**7+252*V**6-72*V**5+24*V**4)*W**4+(-&
&144*V**8+192*V**7+60*V**6-180*V**5+96*V**4-24*V**3)*W**3+(96*V**7-288*&
&V**6+336*V**5-204*V**4+72*V**3-12*V**2)*W**2+(-24*V**6+96*V**5-156*V**&
&4+132*V**3-60*V**2+12*V)*W)+V4*(-12*V**10*W**10+(36*V**10+36*V**9)*W**&
&9+(-36*V**10-144*V**9)*W**8+(12*V**10+180*V**9+144*V**8-96*V**7)*W**7+&
&(-72*V**9-324*V**8+144*V**7+72*V**6)*W**6+(180*V**8+180*V**7-360*V**6+&
&72*V**5)*W**5+(-240*V**7+180*V**6+144*V**5-96*V**4)*W**4+(180*V**6-324&
&*V**5+144*V**4)*W**3+(-72*V**5+180*V**4-144*V**3+36*V**2)*W**2+(12*V**&
&4-36*V**3+36*V**2-12*V)*W)+AL*V2*((-12*V**11+24*V**10-18*V**9+6*V**8)*&
&W**9+(60*V**11-108*V**10+66*V**9-12*V**8-6*V**7)*W**8+(-132*V**11+168*&
&V**10+30*V**9-150*V**8+102*V**7-18*V**6)*W**7+(156*V**11-36*V**10-426*&
&V**9+516*V**8-228*V**7+18*V**5)*W**6+(-96*V**11-216*V**10+744*V**9-492&
&*V**8-156*V**7+360*V**6-162*V**5+18*V**4)*W**5+(24*V**11+264*V**10-444&
&*V**9-264*V**8+960*V**7-756*V**6+198*V**5+36*V**4-18*V**3)*W**4+(-96*V&
&**10-96*V**9+816*V**8-1044*V**7+360*V**6+270*V**5-294*V**4+90*V**3-6*V&
&**2)*W**3+(144*V**9-336*V**8+36*V**7+588*V**6-738*V**5+372*V**4-48*V**&
&3-24*V**2+6*V)*W**2+(-96*V**8+384*V**7-600*V**6+432*V**5-84*V**4-84*V*&
&*3+60*V**2-12*V)*W+24*V**7-120*V**6+252*V**5-288*V**4+192*V**3-72*V**2&
&+12*V)+CQ*V2*((-16*V**11+32*V**10-24*V**9+8*V**8)*W**8+(48*V**11-80*V*&
&*10+40*V**9-8*V**7)*W**7+(-48*V**11+168*V**9-216*V**8+120*V**7-24*V**6&
&)*W**6+(16*V**11+112*V**10-312*V**9+256*V**8-48*V**7-48*V**6+24*V**5)*&
&W**5+(-64*V**10+32*V**9+288*V**8-544*V**7+432*V**6-168*V**5+24*V**4)*W&
&**4+(96*V**9-288*V**8+288*V**7-48*V**6-120*V**5+96*V**4-24*V**3)*W**3+&
&(-64*V**8+272*V**7-480*V**6+456*V**5-248*V**4+72*V**3-8*V**2)*W**2+(16&
&*V**7-80*V**6+168*V**5-192*V**4+128*V**3-48*V**2+8*V)*W)+CQ*V1*((-16*V&
&**11+32*V**10-32*V**9+16*V**8-4*V**7)*W**8+(48*V**11-80*V**10+64*V**9-&
&16*V**8-4*V**7+4*V**6)*W**7+(-48*V**11+144*V**9-240*V**8+180*V**7-72*V&
&**6+12*V**5)*W**6+(16*V**11+112*V**10-304*V**9+320*V**8-140*V**7-12*V*&
&*6+36*V**5-12*V**4)*W**5+(-64*V**10+32*V**9+256*V**8-560*V**7+560*V**6&
&-312*V**5+96*V**4-12*V**3)*W**4+(96*V**9-288*V**8+336*V**7-144*V**6-72&
&*V**5+120*V**4-60*V**3+12*V**2)*W**3+(-64*V**8+272*V**7-512*V**6+560*V&
&**5-384*V**4+164*V**3-40*V**2+4*V)*W**2+(16*V**7-80*V**6+176*V**5-224*&
&V**4+180*V**3-92*V**2+28*V-4)*W))/((V-1)*V**2*W**2*(V*W-1)**4*(V*W-V+1&
&)**3)/3.D0
LM = -LOG(S/M**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V1*(2*V**6*W**6-2*&
&V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2*V*(4*V**5-&
&2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V**2+1)*W**2-2*(4*V&
&**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**2-2*V+1))+4*V*V4*W&
&*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+2*(V-1)*V*(2*V*&
&*2-2*V+1)*V2*(W**2-2*W+2)*(V**2*W**2-2*V*W+1))/((V-1)*V**2*W**2*(V*W-1&
&)**2*(V*W-V+1)**2)
LMP = -LOG(S/MP**2)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+2*(V-1)&
&**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V1*(V**2*W**2-2*(V-1)*V&
&*W+2*(V-1)**2)-2*(V-1)*V2*(V*W-V+1))/((V-1)*V**2*W**2*(V*W-1)**2*(V*W-&
&V+1)**2)
STRUV10=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV11(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV11
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(CQ*V1*(2*V**10*W**9-4*V**8*(V**2+2*V-2)*W**8+3*V**8*(2&
&*V**2+3*V-3)*W**7-V**6*(8*V**4+V**3-15*V**2+28*V-14)*W**6+3*V**6*(2*V*&
&*4+V**3-10*V**2+18*V-9)*W**5-V**4*(4*V**6-3*V**5-7*V**4+26*V**3-28*V**&
&2+18*V-6)*W**4+V**4*(2*V**6-17*V**4+45*V**3-50*V**2+33*V-11)*W**3-(V-1&
&)*V**2*(6*V**6-18*V**5+17*V**4+5*V**2-6*V+2)*W**2+(V-1)**2*V**2*(6*V**&
&4-20*V**3+25*V**2-10*V+5)*W-2*(V-1)**3*(V**4-3*V**3+4*V**2-2*V+1))+V1*&
&(-2*V**10*W**9+4*V**8*(V**2+2*V-2)*W**8-V**8*(6*V**2+11*V-11)*W**7+V**&
&6*(16*V**4-21*V**3+19*V**2+4*V-2)*W**6-V**6*(22*V**4-39*V**3+24*V**2+3&
&0*V-15)*W**5+3*V**4*(4*V**6+V**5-27*V**4+62*V**3-56*V**2+30*V-10)*W**4&
&-V**4*(2*V**6+24*V**5-93*V**4+143*V**3-84*V**2+15*V-5)*W**3+(V-1)*V**2&
&*(6*V**6+6*V**5-65*V**4+128*V**3-89*V**2+30*V-10)*W**2-(V-1)**2*V**2*(&
&6*V**4-12*V**3+V**2+22*V-11)*W+2*(V-1)**3*(V**4-3*V**3+4*V**2-2*V+1))+&
&2*CQ*V**2*V2*W*(2*V**6*(V**2-2*V+2)*W**7-4*V**6*(V**2-V+1)*W**6+2*V**4&
&*(2*V**4-V**3+4*V**2-6*V+3)*W**5-V**4*(4*V**4-6*V**3+13*V**2-14*V+7)*W&
&**4+2*V**2*(V**6+V**5-8*V**4+12*V**3-V**2-6*V+2)*W**3-(V-1)*V**2*(6*V*&
&*4-11*V**3-V**2+24*V-12)*W**2+2*(V-1)**2*(3*V**4-7*V**3+8*V**2-2*V+1)*&
&W-(V-1)**3*(2*V**2-3*V+3))-2*(V-1)*V**2*V2*W*(3*V**6*W**6-V**4*(7*V**2&
&+4*V-4)*W**5+V**4*(V**2+29*V-29)*W**4+3*V**2*(V**4-6*V**3-6*V**2+24*V-&
&12)*W**3-(V-1)*V**2*(9*V**2-43*V+43)*W**2+(V-1)**2*(7*V**2-16*V+16)*W-&
&3*(V-1)**3)+2*(V-1)*V**2*V4*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W&
&+(V-1)**2)*(2*V**2*W**4-V**2*W**3-4*(V**2-2*V+2)*W**2+2*(V-1))+2*(V-1)&
&*V**2*V3*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**&
&2*W**3-3*V**2*W**2+4*(V-1)*W-2*(V**2-2*V+2)))/((V-1)**2*V*W*(V*W-1)**3&
&*(V*W-V+1)**3)
LV1 = -LOG(1-V)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*&
&(2*V1*(V**7*(3*V-1)*W**7-4*V**6*(4*V-3)*W**6+V**5*(3*V**3-2*V**2+25*V-&
&20)*W**5-V**4*(10*V**4-27*V**3+46*V**2-20*V+5)*W**4+V**3*(30*V**4-98*V&
&**3+147*V**2-92*V+31)*W**3-V**2*(30*V**4-103*V**3+151*V**2-89*V+21)*W*&
&*2+V*(10*V**4-33*V**3+47*V**2-23*V+1)*W-(V-1)*(V**2+1))-2*(V-1)*V**2*V&
&4*(V*W-1)*(6*V**3*W**5-4*V**2*W**4-V*(7*V**2+2*V+2)*W**3+(22*V**2-15*V&
&+4)*W**2-(2*V**2+11*V-11)*W+2*(V-1))+CQ*V1*(V*W-V+1)*(2*V**2*W**2-2*V*&
&(V+1)*W+V**2+1)*(V**4*W**4+2*V**2*W**2+1)+2*V*V2*W*(V**2*W**2-2*V*W+1)&
&*(2*V**3*(V**2-2*V+2)*W**3-(V-1)*V**2*(7*V-5)*W**2+V*(10*V**4-37*V**3+&
&60*V**2-55*V+24)*W-(V-1)*(10*V**3-26*V**2+19*V-1))+2*CQ*V*V2*W*(V*W-V+&
&1)*(V**2*W**2+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)-2*(V-1)*V**2*V3*W*(2&
&*(V-1)*V*W**2+(3*V**2-4*V+4)*W+(V-1)*(2*V-5))*(V**2*W**2-2*V*W+1))/((V&
&-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**4)
LV = LOG(V)*(CQ*V1*(2*V**10*W**9-4*V**8*(V**2+2*V-2)*W**8+3*V**8*(2*V*&
&*2+3*V-3)*W**7-V**6*(8*V**4+V**3-15*V**2+28*V-14)*W**6+3*V**6*(2*V**4+&
&V**3-10*V**2+18*V-9)*W**5-V**4*(4*V**6-3*V**5-7*V**4+26*V**3-28*V**2+1&
&8*V-6)*W**4+V**4*(2*V**6-17*V**4+45*V**3-50*V**2+33*V-11)*W**3-(V-1)*V&
&**2*(6*V**6-18*V**5+17*V**4+5*V**2-6*V+2)*W**2+(V-1)**2*V**2*(6*V**4-2&
&0*V**3+25*V**2-10*V+5)*W-2*(V-1)**3*(V**4-3*V**3+4*V**2-2*V+1))+V1*(-2&
&*V**10*W**9+4*V**8*(V**2+2*V-2)*W**8-V**6*(8*V**4+3*V**3-V**2-4*V+2)*W&
&**7+V**6*(16*V**4-23*V**3+17*V**2+12*V-6)*W**6-3*V**4*(6*V**6-9*V**5-2&
&*V**4+24*V**3-17*V**2+6*V-2)*W**5+V**4*(12*V**6-9*V**5-47*V**4+136*V**&
&3-128*V**2+72*V-24)*W**4-V**2*(4*V**8+12*V**7-75*V**6+137*V**5-90*V**4&
&+9*V**3+25*V**2-24*V+6)*W**3+(V-1)*V**2*(12*V**6-24*V**5-7*V**4+72*V**&
&3-61*V**2+30*V-10)*W**2-(V-1)**2*(12*V**6-36*V**5+41*V**4-8*V**3-V**2+&
&6*V-2)*W+4*(V-1)**3*(V**4-3*V**3+4*V**2-2*V+1))-2*V**2*V2*W*(2*V**6*(V&
&**2-2*V+2)*W**7-4*(V-2)**2*V**6*W**6-2*V**4*(6*V**4-17*V**3+42*V**2-50&
&*V+25)*W**5+V**4*(28*V**4-86*V**3+175*V**2-178*V+89)*W**4-2*V**2*(7*V*&
&*6+V**5-54*V**4+164*V**3-227*V**2+174*V-58)*W**3+(V-1)*V**2*(42*V**4-1&
&17*V**3+173*V**2-112*V+56)*W**2-6*(V-1)**2*(7*V**4-21*V**3+30*V**2-18*&
&V+9)*W+(V-1)**3*(14*V**2-37*V+37))+2*CQ*V**2*V2*W*(2*V**6*(V**2-2*V+2)&
&*W**7-4*V**6*(V**2-V+1)*W**6+2*V**4*(2*V**4-V**3+4*V**2-6*V+3)*W**5-V*&
&*4*(4*V**4-6*V**3+13*V**2-14*V+7)*W**4+2*V**2*(V**6+V**5-8*V**4+12*V**&
&3-V**2-6*V+2)*W**3-(V-1)*V**2*(6*V**4-11*V**3-V**2+24*V-12)*W**2+2*(V-&
&1)**2*(3*V**4-7*V**3+8*V**2-2*V+1)*W-(V-1)**3*(2*V**2-3*V+3))+4*(V-1)*&
&V**2*V3*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(4*V**2&
&*W**3-2*V**2*W**2-2*(3*V**2-4*V+4)*W-V**2+8*V-8)-4*(V-1)*V**2*V4*(V**2&
&*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**3-4*(V-1)*W**&
&2+3*(V-1)*W-V+1))/((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**3)
LVW = 2*(V1*(2*V**4*W**3+2*(V-2)**2*V**2*W**2+V*(4*V**3-11*V**2+12*V+1&
&)*W-V**2-1)+V*V2*W*(4*V*W+2*V**2+V+1)-(V-1)*V**2*V4*W*(4*W+3)-(V-1)*V*&
&*2*V3*W*(4*W+3))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**&
&4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*LOG(&
&1-V*W)/((V-1)**2*V*W*(V*W-1)**4*(V*W-V+1)**4)
LTVW = 2*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(2*V*V1*W*(V**5*W&
&**5-V**4*(4*V-5)*W**4+V**3*(9*V**2-21*V+14)*W**3-(V-1)*V**2*(4*V**3-3*&
&V**2-7*V+12)*W**2+(V-1)**2*V*(8*V**3-18*V**2+13*V+3)*W-(V-1)**3*(4*V**&
&3-9*V**2+8*V-1))+CQ*(V-1)*V1*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*&
&(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)-(V-1)*V**2*V4*W*(4*W+3)*(V**&
&3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)-(V-1)*V**2*V3*W*(4*W&
&+3)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)+(V-1)*V*V2*W&
&*(V**3*(9*V-7)*W**3-(V-1)*V**2*(8*V**2+3*V-3)*W**2+(V-1)**2*V*(16*V**2&
&-13*V+3)*W-(V-1)**3*(8*V**2-7*V+3))+2*CQ*(V-1)**2*V*V2*W*(V**2*W**2+(V&
&-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1))*LOG(V*W-V+1)/((V-1)*&
&*2*V*W*(V*W-1)**4*(V*W-V+1)**3)
LW = -2*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3&
&*(V-1)**2*V*W-(V-1)**3)*((V-1)*V**2*V4*(2*V**2*W**4-3*V**2*W**3+4*V**2&
&*W**2-2*(2*V**2-V+1)*W+4*(V-1))+V1*(V**2*(3*V**4-9*V**3+10*V**2-2*V+1)&
&*W**3-2*V**2*(V**4-2*V**3+V**2+2*V-1)*W**2+(3*V**6-10*V**5+13*V**4-5*V&
&**3+3*V-1)*W-3*(V-1)*(V**4-3*V**3+4*V**2-2*V+1))+V**2*V2*W*(2*V**2*(V*&
&*2-2*V+2)*W**3+8*(V-1)*V**2*W**2-2*(7*V**4-24*V**3+43*V**2-38*V+19)*W+&
&(V-1)*(14*V**2-31*V+31))-(V-1)*V**2*V3*W*(6*V**2*W**3-V**2*W**2-12*(V*&
&*2-V+1)*W+12*(V-1)))*LOG(W)/((V-1)**2*V*W*(V*W-1)**4*(V*W-V+1)**4)
CVC = (3*(V-1)*V1*(V*W-1)*(V*W-V+1)*(6*V**10*W**9-8*V**8*(V**2+4*V-3)*&
&W**8-V**6*(18*V**4-121*V**3+109*V**2-24*V+12)*W**7+V**6*(52*V**4-193*V&
&**3+207*V**2-100*V+62)*W**6-V**4*(42*V**6-75*V**5-88*V**4+282*V**3-183&
&*V**2+92*V-36)*W**5+V**4*(4*V**6+95*V**5-385*V**4+462*V**3-24*V**2-202&
&*V+54)*W**4+V**2*(6*V**8-48*V**7+45*V**6+293*V**5-876*V**4+965*V**3-49&
&5*V**2+152*V-36)*W**3-(V-1)*V**2*(18*V**6-102*V**5+229*V**4-256*V**3+1&
&49*V**2-38*V+18)*W**2+(V-1)**2*(18*V**6-76*V**5+131*V**4-114*V**3+75*V&
&**2-28*V+12)*W-2*(V-1)**3*(3*V**4-9*V**3+10*V**2-2*V+1))+12*(V-1)*V2*W&
&*(V*W-1)*(V*W-V+1)*(2*V**8*(V**2-2*V+2)*W**7-V**6*(4*V**4-3*V**3-V**2+&
&8*V-4)*W**6+V**6*(4*V**4-V**3-5*V**2+12*V-6)*W**5-V**4*(4*V**6-5*V**5-&
&5*V**4+8*V**3+26*V**2-36*V+12)*W**4+V**4*(V**2+2*V-2)*(2*V**4-V**3-17*&
&V**2+36*V-18)*W**3-(V-1)*V**2*(6*V**6-10*V**5-7*V**4+22*V**3+19*V**2-3&
&6*V+12)*W**2+(V-1)**2*V**2*(6*V**4-13*V**3+27*V**2-28*V+14)*W-2*(V-1)*&
&*3*(V**2-V+1)*(V**2+2*V-2))-12*(V-1)**2*V**2*V4*(V**2*W**2-2*V*W+1)*(V&
&**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(2*V**3*W**5-2*V**&
&2*(2*V+3)*W**4+2*V*(V**2+9*V-3)*W**3+(V**3-14*V**2+6*V+2)*W**2-(2*V**2&
&-9*V+8)*W+V-1)-6*AL*(V-1)*V1*(V*W-1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*&
&(V**4*W**4+2*V**2*W**2+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2&
&*W**2-4*(V-1)**3*V*W+(V-1)**4)-8*CQ*(V-1)**2*V1*(V**2*W**2-2*V*W+1)*(V&
&**4*W**4+2*V**2*W**2+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W&
&**2-4*(V-1)**3*V*W+(V-1)**4)-12*AL*(V-1)*V*V2*W*(V*W-1)*(V**2*W**2+1)*&
&(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)*&
&*2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-16*CQ*(V-1)**2*V*V2*W*(V**2*W**2&
&+1)*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W&
&**2-4*(V-1)**3*V*W+(V-1)**4)+12*(V-1)**3*V**2*V3*W*(V**2*W**2-2*V*W+1)&
&*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**2*W**3-6*V**2*W**2+(V**2+14*V-&
&14)*W-4*(V-1)))/((V-1)**3*V*W*(V*W-1)**4*(V*W-V+1)**4)/6.D0
LM = -LOG(S/M**2)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3&
&)*(V1*(2*V**7*W**6+2*V**5*(V**2-7*V+4)*W**5+V**3*(6*V**4-23*V**3+46*V*&
&*2-23*V+2)*W**4-2*V**2*(V**5+6*V**4-25*V**3+42*V**2-17*V+3)*W**3+2*V*(&
&3*V**5-14*V**3+29*V**2-8*V+3)*W**2-2*(3*V**5-6*V**4+4*V**3+4*V**2+2*V+&
&1)*W+2*V**4-6*V**3+9*V**2-4*V+3)-4*(V-1)*V**2*V4*(V*W-1)*(V**2*W**4+2*&
&(V-1)*V*W**3-(V**2+6*V-2)*W**2+(2*V+3)*W-1)+2*V*V2*W*(V**2*W**2+1)*(2*&
&V**2*W**2-2*V*(V+1)*W+V**2+1))/((V-1)**2*V*W*(V*W-1)**3*(V*W-V+1)**3)
LMP = -LOG(S/MP**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V1*(2&
&*V**7*W**6-2*V**5*(2*V**2+V-4)*W**5+V**4*(6*V**3-2*V**2-19*V+19)*W**4-&
&4*(V-1)*V**3*(3*V**3-4*V**2-2*V+6)*W**3+2*(V-1)**2*V**2*(7*V**3-14*V**&
&2+8*V+5)*W**2-2*(V-1)**3*V*(4*V**3-10*V**2+9*V-1)*W+(V-1)**5*(2*V**2-2&
&*V+1))-4*(V-1)*V**2*V4*W*(V*W-V+1)*(V**2*W**3-(V-2)*V*W**2-2*(V-1)*(V+&
&1)*W+2*(V-1)**2)+2*(V-1)**2*V*V2*W*(V**2*W**2+(V-1)**2)*(2*V**2*W**2-2&
&*V*(2*V-1)*W+2*V**2-2*V+1))/((V-1)**2*V*W*(V*W-1)**4*(V*W-V+1)**3)
STRUV11=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV12(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV12
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)*&
&*2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V**7*(3*V+2)*&
&W**7+V**4*(16*V**4+24*V**3+11*V**2+1)*W**6-2*V**4*(8*V**4+12*V**3+21*V&
&**2-V+2)*W**5+2*V**2*(8*V**6+4*V**5+29*V**4+9*V**3-6*V**2+V-1)*W**4-2*&
&V**2*(6*V**6+2*V**5+11*V**4+30*V**3-16*V**2-2*V-3)*W**3+(4*V**8+12*V**&
&7-17*V**6+54*V**5-5*V**4-32*V**3+V**2-2*V+1)*W**2-2*(V-1)*(4*V**6-2*V*&
&*5+7*V**4+6*V**3+6*V**2-4*V-1)*W+4*(V-1)**2*(V**2+1)*(V**2-V+2))-4*V**&
&8*W**8+4*V**7*(3*V+1)*W**7-V**4*(16*V**4+12*V**3+7*V**2+1)*W**6+2*V**4&
&*(8*V**4+2*V**3+19*V**2-3*V+2)*W**5-2*V**2*(8*V**6-6*V**5+27*V**4+13*V&
&**3-10*V**2+V-1)*W**4+2*V**2*(6*V**6-4*V**5+5*V**4+46*V**3-24*V**2-2*V&
&-3)*W**3-(4*V**8+8*V**7-29*V**6+74*V**5+7*V**4-56*V**3+9*V**2-2*V+1)*W&
&**2+2*(V-1)*(4*V**6-6*V**5+9*V**4+6*V**3+12*V**2-8*V-1)*W-4*(V-1)**2*(&
&V**2+1)*(V**2-2*V+3))+N*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1&
&)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V**7*(3*V+1&
&)*W**7+V**4*(16*V**4+10*V**3+3*V**2-2*V+1)*W**6-2*V**4*(8*V**4+V**3+7*&
&V**2-6*V+2)*W**5+2*V**2*(8*V**6-7*V**5+10*V**4-8*V**3-V**2+3*V-1)*W**4&
&-2*(V-1)*V**2*(6*V**5+V**4-4*V**3+8*V**2-8*V+3)*W**3+(V-1)**2*(4*V**6+&
&16*V**5-11*V**4+6*V**3-2*V**2-2*V+1)*W**2-2*(V-1)**3*(4*V**4+2*V**3+V*&
&*2-1)*W+4*(V-1)**4*(V**2+1))-4*V**8*W**8+4*V**7*(3*V+1)*W**7-V**4*(16*&
&V**4+30*V**3-17*V**2-2*V+1)*W**6+2*V**4*(8*V**4+27*V**3-5*V**2-20*V+2)&
&*W**5-2*V**2*(8*V**6+19*V**5+24*V**4-48*V**3-V**2+3*V-1)*W**4+2*(V-1)*&
&V**2*(6*V**5+11*V**4+38*V**3+16*V**2-16*V+3)*W**3-(V-1)*(4*V**7+12*V**&
&6+13*V**5+29*V**4-8*V**3-20*V**2+3*V-1)*W**2+2*(V-1)**2*(V**2+1)*(4*V*&
&*3-2*V**2+5*V+1)*W-4*(V-1)**4*(V**2+1))+2*N**4*(V-1)*VC*W*(V**2*W**2-2&
&*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**8*W**7-V**7&
&*(6*V-5)*W**6+V**6*(10*V**2-16*V+11)*W**5-V**4*(30*V**4-67*V**3+67*V**&
&2-31*V+4)*W**4+2*(V-1)*V**3*(31*V**4-73*V**3+79*V**2-43*V+8)*W**3-2*CQ&
&*(V**2*W-(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(V-1)&
&*V*(2*V**2-2*V+1)*W+(V-1)**2)-(V-1)**2*V**2*(58*V**4-151*V**3+180*V**2&
&-101*V+24)*W**2+(V-1)**3*V*(22*V**4-62*V**3+81*V**2-53*V+16)*W-(V-1)**&
&4*(V**2-V+1)*(2*V**2-3*V+4))-2*N**2*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V*&
&*4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4&
&)*(2*V**7*W**6-2*V**5*(V**2-V-1)*W**5-CQ*(V*W-1)*(V**5*W**4-(V-1)*V**4&
&*W**3-2*V**2*(3*V**2-5*V+1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4)*W+(V-1)**&
&2*(V**3-2*V**2+5*V-2))+V**4*(4*V**3-13*V**2-V+3)*W**4-V**3*(12*V**4-42&
&*V**3+29*V**2-9*V+4)*W**3+V**2*(10*V**5-35*V**4+23*V**3-V**2+3*V-4)*W*&
&*2-(V-1)*V*(2*V**5-27*V**3+32*V**2-19*V+6)*W+(V-1)**2*(2*V**4-10*V**3+&
&11*V**2-12*V+5))-2*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V&
&**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**5-V**4&
&*(4*V**2+V-1)*W**4+2*V**3*(V+2)*(3*V**2-3*V+1)*W**3+CQ*V*(V*W-1)*(V**2&
&*W**2+(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)-2*V**2*(2*V**4+4&
&*V**3-5*V**2+2*V-2)*W**2+(V-1)*V**2*(V**3+5*V**2-V+5)*W-(V-1)**2*(V**3&
&+V**2-V+1))-4*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**&
&2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*&
&W**3-3*V**2*W**2+3*V*W-1)+4*(CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W&
&**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**&
&2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)&
&**3*(V*W-V+1)**6)
LV1 = LOG(1-V)*(N**4*(V-1)*VC*W*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+1&
&0*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(1&
&2*V**7*W**6-36*V**6*W**5+4*V**5*(3*V**2-3*V+13)*W**4-2*V**3*(12*V**4-1&
&5*V**3+16*V**2+11*V+1)*W**3+V**2*(50*V**4-98*V**3+94*V**2-21*V+2)*W**2&
&-V*(28*V**4-62*V**3+56*V**2-15*V-2)*W+(V-1)*(2*V**3-3*V**2+2))-2*N**2*&
&(V-1)*VC*W*(V*W-1)*(6*V**6*W**5+V**5*(2*V-19)*W**4+V**3*(8*V**3-21*V**&
&2+39*V-1)*W**3-V**2*(15*V**3-28*V**2+36*V+3)*W**2+V*(8*V**3-9*V**2+9*V&
&+7)*W-V**3-V**2+2*V-3)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W&
&**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**&
&6)+N**3*VC*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-&
&10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(2*(V**2*W**2-2*V*W+1)*&
&(2*V**5*W**5-2*V**4*(2*V+1)*W**4+V**3*(3*V**2+4*V+5)*W**3-V**2*(V**3+5&
&*V**2+7*V-1)*W**2+2*V*(3*V**3+V**2+V-1)*W-2*(V-1)*(V+1)*(V**2+1))+CQ*(&
&V*W-V+1)*(2*V**2*W**2-2*V*(V+1)*W+V**2+1)*(V**4*W**4-4*V**3*W**3+8*V**&
&2*W**2-8*V*W+4))+(V-1)*VC*W*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V&
&-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(2*V**&
&6*W**5+2*V**4*(V**2-3*V+1)*W**4-2*V**3*(2*V**3+3*V**2-7*V+3)*W**3+V**2&
&*(2*V-1)*(2*V**3+4*V**2+3*V-6)*W**2-V*(8*V**4-8*V**2+7*V-6)*W+(V-1)*(V&
&+2)*(4*V**2-5*V+4))-N*V*VC*W*(V*W-1)*(2*(V**2*W**2-2*V*W+1)*(2*V**4*W*&
&*4-4*V**4*W**3+V**2*(3*V**2+1)*W**2-(V-1)*V*(V**2+4*V-1)*W+4*(V-1)**2*&
&(V+1))+CQ*V*W*(V*W-V+1)*(V**2*W**2-2*V*W+2)*(2*V**2*W**2-2*V*(V+1)*W+V&
&**2+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3&
&*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)*&
&*3*(V*W-V+1)**6)
LV = LOG(V)*(-N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*&
&V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4*V**8*W**8-4*V**7*(3*V+2)*W**&
&7+3*V**6*(V+1)*(5*V+3)*W**6-2*V**5*(6*V**3+12*V**2+17*V-1)*W**5+V**4*(&
&9*V**4+6*V**3+46*V**2+14*V-15)*W**4-2*V**3*(3*V+5)*(V**4-2*V**3+6*V**2&
&-2*V-1)*W**3+2*V**2*(V**6+2*V**5-7*V**4+20*V**3+3*V**2-18*V+3)*W**2-2*&
&(V-1)*V*(2*V**5-3*V**4+4*V**3+4*V**2+6*V-5)*W+2*(V-1)**2*(V**2+1)*(V**&
&2-2*V+3))-4*V**8*W**8+4*V**7*(3*V+2)*W**7-V**4*(16*V**4+24*V**3+11*V**&
&2+1)*W**6+2*V**4*(8*V**4+6*V**3+33*V**2-7*V+2)*W**5-2*V**2*(8*V**6-10*&
&V**5+47*V**4+15*V**3-16*V**2+V-1)*W**4+2*V**2*(6*V**6-10*V**5+11*V**4+&
&64*V**3-36*V**2-4*V-3)*W**3-(4*V**8+4*V**7-41*V**6+106*V**5+7*V**4-76*&
&V**3+13*V**2-2*V+1)*W**2+2*(V-1)*(4*V**6-10*V**5+11*V**4+6*V**3+18*V**&
&2-12*V-1)*W-4*(V-1)**2*(V**2+1)*(V**2-3*V+4))+N*VC*(V*W-1)*(V**4*W**4-&
&4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(4&
&*V**8*W**8-4*V**7*(3*V+1)*W**7+V**6*(15*V**2+12*V+1)*W**6-2*V**5*(6*V*&
&*3+5*V**2+3*V-2)*W**5+V**4*(9*V**4-2*V**3+12*V**2-10*V-1)*W**4-2*(V-1)&
&*V**3*(3*V**4+V**3-V**2+5*V-2)*W**3+2*(V-1)**2*V**2*(V**4+4*V**3-3*V**&
&2+2)*W**2-2*(V-1)**3*V*(2*V**3+V**2-1)*W+2*(V-1)**4*(V**2+1))-4*V**8*W&
&**8+4*V**7*(3*V+1)*W**7-V**4*(2*V-1)*(8*V**3+21*V**2-1)*W**6+2*V**4*(8&
&*V**4+29*V**3-V**2-26*V+2)*W**5-2*V**2*(8*V**6+17*V**5+34*V**4-52*V**3&
&-5*V**2+3*V-1)*W**4+2*(V-1)*V**2*(6*V**5+9*V**4+36*V**3+28*V**2-20*V+3&
&)*W**3-(V-1)*(4*V**7+12*V**6+5*V**5+33*V**4+8*V**3-32*V**2+3*V-1)*W**2&
&+2*(V-1)**2*(4*V**5-2*V**4+7*V**3-V**2+7*V+1)*W-4*(V-1)**4*(V**2+1))+N&
&**4*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2&
&+3*V*W-1)*(4*V**8*W**7-4*V**7*(3*V-2)*W**6+V**4*(24*V**4-32*V**3+20*V*&
&*2-4*V+1)*W**5-V**3*(80*V**5-196*V**4+212*V**3-120*V**2+33*V-4)*W**4+2&
&*(V-1)*V**2*(82*V**5-222*V**4+268*V**3-162*V**2+41*V-3)*W**3-4*CQ*(V**&
&2*W-(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(V-1)*V*(2&
&*V**2-2*V+1)*W+(V-1)**2)-2*(V-1)**2*V*(78*V**5-228*V**4+296*V**3-180*V&
&**2+49*V-2)*W**2+(V-1)**3*(64*V**5-192*V**4+264*V**3-176*V**2+57*V-1)*&
&W-(V-1)**4*(8*V**4-20*V**3+32*V**2-24*V+13))-2*N**2*(V-1)*VC*W*(V**2*W&
&**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)&
&**3*V*W+(V-1)**4)*(2*V**7*W**6-V**6*(6*V-7)*W**5-CQ*(V*W-1)*(V**5*W**4&
&-(V-1)*V**4*W**3-2*V**2*(3*V**2-5*V+1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4&
&)*W+(V-1)**2*(V**3-2*V**2+5*V-2))+V**3*(2*V-1)*(5*V**3-9*V**2+V-1)*W**&
&4-(V-1)*V**2*(14*V**4-26*V**3+10*V**2-11*V+1)*W**3+V*(12*V**6-35*V**5+&
&32*V**4-13*V**3+8*V**2-4*V-1)*W**2-(V-1)**2*(4*V**5+V**4-16*V**3+14*V*&
&*2-10*V+1)*W+(V-1)**2*(4*V**4-13*V**3+17*V**2-16*V+7))-(V-1)*VC*W*(V**&
&3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V*&
&*2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(6*V**5*W**4-V**2*(14*V**3-2*V+1)*W**&
&3+2*CQ*V*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+V*&
&(12*V**4+6*V**3-18*V**2+7*V-2)*W**2-(V-1)*(2*V**4+10*V**3-2*V**2+5*V-1&
&)*W-(V-1)**2*(2*V**3-2*V**2+4*V-1))-4*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC&
&*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2&
&-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)+4*(CQ-1)*GTR*N*(V-1)&
&**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)&
&*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1))/(N**2*(V&
&-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**6)
LVW = (-2*N**3*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**4*W**4-2*V**3*&
&(V+3)*W**3+V**2*(V**2+4*V+11)*W**2-2*V*(V**2+3*V+4)*W+4*(V**2+1))*(V**&
&6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+1&
&5*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)-2*N**4*(V-1)*VC*W*(V**3*&
&W**3-3*V**2*W**2+3*V*W-1)*(4*V**4*W**3+V**3*(4*V-9)*W**2+V**2*(8*V**2-&
&16*V+13)*W-V**2+2*V-2)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W&
&**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**&
&6)+2*N**2*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**4*W**3+2*V*&
&*3*(2*V-3)*W**2+V**2*(8*V**2-17*V+12)*W-(V-2)*(V-1))*(V**6*W**6-6*(V-1&
&)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V*&
&*2*W**2-6*(V-1)**5*V*W+(V-1)**6)+2*N*V*VC*W*(V**3*W**3-3*V**2*W**2+3*V&
&*W-1)*(2*V**3*W**3-2*V**2*(V+1)*W**2+V*(V+1)**2*W-2*(V-1)*(V+1))*(V**6&
&*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15&
&*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)+2*(V-1)*VC*W*(V**3*W**2-V&
&**2*(V+4)*W+2*(V**2-V+2))*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**6*W**6-6&
&*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)*&
&*4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6))*LOG(1-V*W)/(N**2*(V-1)**2*V**2*&
&W**2*(V*W-1)**3*(V*W-V+1)**6)
LTVW = (-2*N**4*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3&
&-3*V**2*W**2+3*V*W-1)*(V**7*W**6+2*V**6*(4*V-5)*W**5-V**4*(9*V**3-12*V&
&**2+2*V-1)*W**4+4*(V-1)*V**3*(3*V**2-1)*W**3+4*CQ*(V**2*W-(V-1)**2)*(V&
&**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*(V-1)*V*(2*V**2-2*V+1)*W+(&
&V-1)**2)-(V-1)**2*V**2*(9*V**3-6*V**2+7*V-6)*W**2+2*(V-1)**3*V*(4*V**3&
&-5*V**2+7*V-2)*W+(V-1)**4*(V+1)*(V**2-V+1))+2*N**2*(V-1)*VC*W*(V**3*W*&
&*3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W&
&**2-4*(V-1)**3*V*W+(V-1)**4)*(2*CQ*(V**5*W**4-(V-1)*V**4*W**3-2*V**2*(&
&3*V**2-5*V+1)*W**2-(V-1)*V*(V**3-8*V**2+9*V-4)*W+(V-1)**2*(V**3-2*V**2&
&+5*V-2))+2*V**5*W**4-V**4*(9*V-5)*W**3+2*V**2*(7*V**3-6*V**2+1)*W**2-(&
&V-1)*V*(9*V**3-2*V**2+V+4)*W+2*(V-1)**2*(V+1)*(V**2-V+1))+4*N*VC*(V**3&
&*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**&
&2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*W**2-2*V**2*W+V**2+1)*(V**4*&
&W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)+(V-1)*V&
&*W*(3*V**3*W**3-V**2*(7*V-1)*W**2+(V-1)*V*(7*V+3)*W-3*(V-1)**2*(V+1)))&
&-4*N**3*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**&
&3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*W**2-2*V**2*&
&W+V**2+1)*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(&
&V-1)**4)+2*(V-1)**2*V**3*(W-1)*W**2)-2*(V-1)*VC*W*(V**3*W**3-3*V**2*W*&
&*2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)*&
&*3*V*W+(V-1)**4)*(-3*V**5*W**4+V**4*(6*V+1)*W**3+2*CQ*V*(V**2*W**2+(V-&
&1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)-V**2*(6*V**3-3*V**2+2*V+3&
&)*W**2+(V-1)*V*(6*V**3-3*V**2-5*V+6)*W-(V-1)**2*(V+1)*(3*V**2-2*V+3))-&
&8*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**&
&2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*&
&W**2+3*V*W-1)+8*(CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2&
&)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W&
&**3-3*V**2*W**2+3*V*W-1))*LOG(V*W-V+1)/(N**2*(V-1)**2*V**2*W**2*(V*W-1&
&)**3*(V*W-V+1)**6)
LW = (2*N**2*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**3*(4*V**2-5*V+2)*W**3-&
&V*(2*V**4-2*V**3+5*V**2-2*V+1)*W**2+(6*V**5-11*V**4+12*V**3-5*V**2+4*V&
&+1)*W-(2*V**2-V+2)*(3*V**2-5*V+3))*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-&
&1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5&
&*V*W+(V-1)**6)-N**3*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*(2*V**4*W**4&
&-V**3*(V+1)*(V+3)*W**3+V**2*(3*V**3-V**2+13*V-3)*W**2-2*V*(2*V-1)*(V**&
&3-V**2+3*V+1)*W+2*(V-1)*(V**2+1)*(V**2-V+2))-CQ*(V**2+1)**2*(V*W-V+1)*&
&(W**2-2*W+2))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V&
&-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-(V-1)*VC*W*(V**2*W**2-2*V*W+&
&1)*(4*V**5*W**4-V**2*(4*V**3+10*V**2+1)*W**3+2*V**2*(V**3+5*V**2+2*V+1&
&)*W**2+(2*V**5-8*V**4+V**2-5*V+1)*W-(V-1)*(2*V**3-6*V**2+8*V-5))*(V**5&
&*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*&
&(V-1)**4*V*W-(V-1)**5)-N*(V-1)*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*(&
&V**3*(V+3)*W**3-V**2*(3*V**2+5)*W**2+2*(V-1)*V*(2*V**2-V+1)*W-2*(V-1)*&
&*2*(V**2+1))+CQ*(V-1)*(V**2+1)*(V*W-V+1)*(W**2-2*W+2))*(V**5*W**5-5*(V&
&-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V&
&*W-(V-1)**5)-N**4*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**4*W&
&**3-V*(12*V**4-14*V**3+6*V**2-4*V+1)*W**2+(24*V**5-68*V**4+84*V**3-50*&
&V**2+12*V-1)*W-(V-1)*(12*V**4-24*V**3+28*V**2-18*V+7))*(V**5*W**5-5*(V&
&-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V&
&*W-(V-1)**5))*LOG(W)/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**6)
CVC = (-N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-&
&1)**2)*(6*V**9*W**8-6*V**8*(V+1)*W**7-2*V**6*(45*V**3-88*V**2+42*V-3)*&
&W**6+V**5*(250*V**4-678*V**3+688*V**2-291*V+16)*W**5-V**4*(234*V**5-75&
&4*V**4+884*V**3-376*V**2+5*V-6)*W**4+2*(V-1)*V**3*(33*V**5-24*V**4-138&
&*V**3+181*V**2-46*V+8)*W**3+2*(V-1)**2*V**2*(7*V**5-66*V**4+129*V**3-1&
&07*V**2+83*V-7)*W**2-(V-1)**3*V**2*(6*V**4-4*V**3-38*V**2+27*V-19)*W+(&
&V-1)**4*(6*V**4-12*V**3+8*V**2+9*V+2))+N**2*(V-1)*VC*W*(V*W-1)*(V**4*W&
&**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(6&
&*V**8*W**7+6*(V-5)*V**7*W**6-V**4*(84*V**4-174*V**3+79*V**2-15*V+6)*W*&
&*5+V**4*(120*V**4-204*V**3+69*V**2+19)*W**4-V**2*(42*V**6+54*V**5-303*&
&V**4+296*V**3-117*V**2+42*V-12)*W**3-V**2*(6*V**6-102*V**5+245*V**4-23&
&2*V**3+156*V**2-95*V+26)*W**2+(V-1)*(12*V**6-66*V**5+121*V**4-103*V**3&
&+71*V**2-21*V+6)*W-(V-1)**2*(6*V**4-6*V**3-16*V**2+27*V-7))-3*N*VC*(V*&
&W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+&
&(V-1)**4)*(4*V**7*W**7+V**4*(2*V**4-16*V**3+3*V**2+V-2)*W**6-V**4*(4*V&
&**4-20*V**3-6*V**2+17*V-7)*W**5+2*V**2*(V**6-6*V**5-2*V**4+7*V**3-V**2&
&-3*V+2)*W**4+2*(V-1)*V**2*(4*V**4+3*V**3-11*V**2+2*V+5)*W**3-(V-1)**2*&
&(4*V**5+13*V**4-11*V**3-11*V**2-V+2)*W**2+(V-1)**3*(8*V**3+V**2-3)*W-4&
&*(V-1)**4*V)+3*(V-1)*VC*W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)&
&**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**7*W**6-(V-2)*(V-1)*V**4*(&
&2*V-1)*W**5-V**4*(19*V**2-17*V+7)*W**4-V**2*(2*V**5-29*V**4+30*V**3+V*&
&*2-18*V+4)*W**3+V**2*(2*V**5-13*V**4-5*V**3+50*V**2-47*V+10)*W**2-(V-1&
&)*(4*V**5-21*V**4+17*V**3+9*V**2-9*V+2)*W+(V-1)**2*(2*V**3-8*V**2+6*V-&
&3))+3*AL*N**3*VC*(V*W-1)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1)**2*V**4&
&*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)&
&**6)*(2*V**6*W**6-2*V**5*(V+5)*W**5+V**2*(2*V**4+8*V**3+27*V**2+1)*W**&
&4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V+1)*W**3+(2*V**6+4*V**5+13*V**4+24&
&*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+V**4+8*V**3+6*V**2+10*V+1)*W+2*(V*&
&*2+1)*(V**2+3))-3*AL*N*VC*(V*W-1)*(V**6*W**6-6*(V-1)*V**5*W**5+15*(V-1&
&)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)**4*V**2*W**2-6*(V-1)**5*&
&V*W+(V-1)**6)*(2*V**6*W**6-2*V**5*(V+3)*W**5+V**2*(2*V**4+2*V**3+11*V*&
&*2-2*V+1)*W**4-2*V*(V**5+2*V**3+3*V**2-V+1)*W**3+(V**2+1)*(2*V**4-3*V*&
&*2+2*V+1)*W**2-2*(V-1)**2*(2*V+1)*(V**2+1)*W+2*(V-1)**2*(V**2+1))+3*N*&
&*3*VC*(V*W-1)*(8*V**5*W**5+V**2*(2*V**4-10*V**3-25*V**2+3*V-2)*W**4+V*&
&(2*V**4+28*V**3+33*V**2-3*V+4)*W**3-(8*V**5+3*V**4+37*V**3+27*V**2+3*V&
&+2)*W**2+(16*V**4-11*V**3+39*V**2+9*V+3)*W-8*V*(V**2-V+2))*(V**6*W**6-&
&6*(V-1)*V**5*W**5+15*(V-1)**2*V**4*W**4-20*(V-1)**3*V**3*W**3+15*(V-1)&
&**4*V**2*W**2-6*(V-1)**5*V*W+(V-1)**6)-4*GTR*N*(V-1)*VC*W*(V**2*W**2-2&
&*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**6*W**5+V**4*(&
&V**2-2*V-1)*W**4+2*(V-1)*V**3*(14*V**2-15*V+2)*W**3-2*(V-1)**2*V**2*(1&
&4*V**2+V+3)*W**2-(V-1)**3*V*(V**2+V-4)*W-(V-1)**4*(V**2+1))+4*CQ*N**3*&
&VC*(V**2*W**2-2*V*W+1)*(V**5*(V+2)*W**5-V**4*(V**2+5*V-3)*W**4+V**2*(V&
&**4+7*V**3-V**2-4*V+1)*W**3-V**2*(V**4+3*V**3+4*V**2-5*V+1)*W**2+(V-1)&
&*(4*V**4-6*V**3+13*V**2-4*V+1)*W-(V-1)**2*(3*V**2-6*V+7))*(V**5*W**5-5&
&*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**&
&4*V*W-(V-1)**5)-4*CQ*N*VC*(V**2*W**2-2*V*W+1)*(V**5*(V+2)*W**5-V**4*(V&
&**2+3*V-1)*W**4+(V-1)*V**2*(V**3+4*V**2+V-1)*W**3-(V-1)**2*V**2*(V**2+&
&3*V+3)*W**2+(V-1)**3*(4*V**2+1)*W-3*(V-1)**4)*(V**5*W**5-5*(V-1)*V**4*&
&W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)*&
&*5)-12*CQ*N**2*(V-1)*V*VC*W*(V**2*W**2+(V-1)*V*W+(V-1)**2)*(V**3*W**3-&
&3*V**2*W**2+3*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**&
&3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+12*CQ*(V-1)*V*VC*W*(V&
&**2*W**2+(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**5*W**&
&5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1&
&)**4*V*W-(V-1)**5)-4*GTR*N**3*(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)*&
&*2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**6*(2*V-3)*W**5+V**4*(4*V**4-12&
&*V**3+6*V**2+3*V+1)*W**4-(V-1)*V**3*(12*V**4-32*V**3+45*V**2-37*V+4)*W&
&**3+(V-1)**2*V**2*(12*V**4-36*V**3+53*V**2-11*V+6)*W**2-2*(V-1)**3*V*(&
&2*V**4-7*V**3+5*V**2-6*V+2)*W+(V-1)**4*(V**2+1))+24*CQ*GTR*N**3*(V-1)*&
&*2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(&
&V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)-24*CQ*GTR*N*(V-1)&
&**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*&
&(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3))/(N**2*(V-1)**2*&
&V**2*W**2*(V*W-1)**3*(V*W-V+1)**6)/3.D0
LM = LOG(S/M**2)*(N**3*VC*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)&
&**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**6-2*V**5*(V+5)*W**5+&
&V**2*(2*V**4+8*V**3+27*V**2+1)*W**4-2*V*(V**5+3*V**4+10*V**3+20*V**2+V&
&+1)*W**3+(2*V**6+4*V**5+13*V**4+24*V**3+36*V**2+4*V+1)*W**2-2*(2*V**5+&
&V**4+8*V**3+6*V**2+10*V+1)*W+2*(V**2+1)*(V**2+3))-N*VC*(V*W-1)*(V**4*W&
&**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2&
&*V**6*W**6-2*V**5*(V+3)*W**5+V**2*(2*V**4+2*V**3+11*V**2-2*V+1)*W**4-2&
&*V*(V**5+2*V**3+3*V**2-V+1)*W**3+(V**2+1)*(2*V**4-3*V**2+2*V+1)*W**2-2&
&*(V-1)**2*(2*V+1)*(V**2+1)*W+2*(V-1)**2*(V**2+1))+N**4*(V-1)*VC*W*(V*W&
&-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(&
&V-1)**4)*(4*V**6*W**5+4*(V-4)*V**5*W**4+V**2*(12*V**4-32*V**3+40*V**2-&
&4*V+1)*W**3-V*(4*V**5+16*V**4-40*V**3+40*V**2-5*V+2)*W**2+(8*V**5-4*V*&
&*4-4*V**3+7*V**2+2*V+1)*W-4*V**4+8*V**3-9*V**2+6*V-3)-2*N**2*(V-1)*VC*&
&W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3&
&*V*W+(V-1)**4)*(2*V**6*W**5+2*(V-4)*V**5*W**4+V**2*(6*V**4-16*V**3+22*&
&V**2-3*V+1)*W**3-V*(2*V**5+8*V**4-18*V**3+22*V**2-3*V+2)*W**2+(4*V**5-&
&2*V**4+V**3+2*V**2+3*V+1)*W-2*V**4+4*V**3-6*V**2+5*V-3)+(V-1)*VC*W*(V*&
&W-1)*(V**2*(4*V**2-2*V+1)*W**3-V*(4*V**3+4*V**2-V+2)*W**2+(6*V**3-3*V*&
&*2+4*V+1)*W-3*V**2+4*V-3)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2&
&*W**2-4*(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V&
&*W-V+1)**4)
LMP = LOG(S/MP**2)*(4*N**4*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*&
&(V**8*W**7-V**7*(3*V-2)*W**6+V**6*(5*V**2-8*V+5)*W**5-V**4*(11*V**4-26&
&*V**3+25*V**2-10*V+1)*W**4+(V-1)**2*V**3*(19*V**3-26*V**2+19*V-4)*W**3&
&-(V-1)**2*V**2*(17*V**4-44*V**3+49*V**2-26*V+6)*W**2+(V-1)**4*V*(7*V**&
&3-12*V**2+11*V-4)*W-(V-1)**4*(V**2-V+1)**2)-2*N**2*(V-1)*VC*W*(V**2*W*&
&*2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**6*W**5-&
&V**5*(6*V-5)*W**4+V**4*(8*V**2-13*V+7)*W**3-2*(V-1)*V**2*(4*V**3-6*V**&
&2+6*V-1)*W**2+(V-1)*V*(6*V**4-15*V**3+20*V**2-13*V+4)*W-(V-1)**4*(2*V*&
&*2-V+2))+2*N**3*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*&
&W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**4-(V-1)*V**3*W**3+(&
&V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)-2*N*VC*(V**2*W**2-2*(V-1)*V*W&
&+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)&
&*(V**4*W**4-(V-1)*V**3*W**3+(V-1)**2*V**2*W**2-(V-1)**3*V*W+(V-1)**4)+&
&2*(V-1)*V*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2+(V-1)*V*W+(&
&V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)+4&
&*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*V**2&
&*W+V**2+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)-4*GTR*N*(V-1)**2*V**2*VC*W*&
&*2*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)*(V**3*W**3-3*V**2*&
&W**2+3*V*W-1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**4)
STRUV12=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV13(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV13
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**&
&4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*V**8*W**8-4*V**7*(5*V-2&
&)*W**7+V**4*(52*V**4-50*V**3+17*V**2-4*V+1)*W**6-2*V**4*(42*V**4-59*V*&
&*3+30*V**2-7*V+2)*W**5+2*V**2*(44*V**6-64*V**5+15*V**4+25*V**3-16*V**2&
&+5*V-1)*W**4-2*V**2*(28*V**6-22*V**5-60*V**4+114*V**3-71*V**2+20*V-3)*&
&W**3+(16*V**8+40*V**7-204*V**6+280*V**5-150*V**4+12*V**3+15*V**2-6*V+1&
&)*W**2-2*(V-1)*(16*V**6-28*V**5+6*V**4+30*V**3-29*V**2+10*V-1)*W+4*(V-&
&1)**2*(2*V**2-3*V+2)*(2*V**2-2*V+1))-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V&
&**4*(36*V**4-30*V**3+13*V**2-4*V+1)*W**6+2*V**4*(28*V**4-39*V**3+22*V*&
&*2-5*V+2)*W**5-2*V**2*(32*V**6-46*V**5+V**4+37*V**3-20*V**2+5*V-1)*W**&
&4+2*V**2*(24*V**6-20*V**5-66*V**4+130*V**3-79*V**2+20*V-3)*W**3-(16*V*&
&*8+32*V**7-204*V**6+292*V**5-138*V**4-12*V**3+23*V**2-6*V+1)*W**2+2*(V&
&-1)*(16*V**6-32*V**5+4*V**4+46*V**3-43*V**2+14*V-1)*W-4*(V-1)**2*(2*V*&
&*2-4*V+3)*(2*V**2-2*V+1))-N*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**&
&2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*V**8*W**8-4*V**7&
&*(4*V-1)*W**7+V**4*(28*V**4-14*V**3+3*V**2-2*V+1)*W**6-2*V**4*(12*V**4&
&-5*V**3+V**2-2*V+2)*W**5+2*V**2*(4*V**6+6*V**5-5*V**4+2*V**3-V**2+3*V-&
&1)*W**4-2*V**2*(6*V**5+2*V**3-6*V**2+7*V-3)*W**3+(12*V**6-10*V**4+2*V*&
&*3+3*V**2-4*V+1)*W**2-2*(V-1)*(6*V**4-5*V**2+4*V-1)*W+4*(V-1)**2*(2*V*&
&*2-2*V+1))-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V**4*(28*V**4+6*V**3-17*V**&
&2-2*V+1)*W**6+2*V**4*(12*V**4+35*V**3-53*V**2+12*V+2)*W**5-2*V**2*(4*V&
&**6+72*V**5-111*V**4+42*V**3-V**2+3*V-1)*W**4+2*V**2*(58*V**5-86*V**4+&
&20*V**3+18*V**2-V-3)*W**3-(32*V**7-4*V**6-124*V**5+178*V**4-98*V**3+23&
&*V**2-4*V+1)*W**2+2*(V-1)*(2*V**2-2*V+1)*(8*V**3-11*V**2+8*V-1)*W-4*(V&
&-1)**2*(2*V**2-2*V+1))-2*N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**&
&4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V&
&**8*W**7-V**7*(V+5)*W**6+V**6*(5*V**2-6*V+11)*W**5-V**4*(3*V**4+10*V**&
&3-2*V**2+15*V+4)*W**4+2*V**3*(2*V**4+12*V**3-2*V**2+11*V+8)*W**3-2*CQ*&
&(V**2*W-1)*(V**2*W**2-2*V*W+1)*((V-1)**2*V**2*W**2-2*V*(V**2+1)*W+(V-1&
&)**2)-V**2*(10*V**4-2*V**3+21*V**2+5*V+24)*W**2+V*(4*V**4-5*V**3+18*V*&
&*2-11*V+16)*W-(V**2-V+1)*(3*V**2-5*V+4))+2*N**2*(V-1)*VC*W*(V**3*W**3-&
&3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V*&
&*2*W**2-4*V*W+1)*(2*V**7*W**6+2*V**5*(V**2-3*V+1)*W**5-CQ*(V-1)*(V*W-V&
&+1)*(V**5*W**4-V**4*W**3+2*(V-1)*V**2*(V**2-3*V-1)*W**2+V*(2*V**3-2*V*&
&*2+3*V-4)*W+2*V**3-2*V**2-V+2)-V**4*(7*V**3-6*V**2-8*V+3)*W**4+V**3*(6&
&*V**4+5*V**3-26*V**2+7*V-4)*W**3-V**2*(4*V**5+2*V**3-29*V**2+17*V-4)*W&
&**2+V*(6*V**5-4*V**4-15*V**3+16*V**2-11*V+6)*W-(V-1)*(4*V**4-4*V**3-5*&
&V**2+8*V-5))+2*(V-1)**2*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V&
&*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**6*W**5-V*&
&*4*(4*V**2+V-1)*W**4+2*V**3*(3*V-2)*(V**2+V+1)*W**3+CQ*V*(V*W-V+1)*(V*&
&*2*W**2+V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)-2*V**2*(V**4+8*V**3-1&
&1*V**2+6*V-2)*W**2+V**2*(10*V**3-18*V**2+14*V-5)*W-(V-1)*(2*V**3-2*V**&
&2+2*V-1))+4*(CQ-1)*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*&
&W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**&
&3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-4*(CQ-1)*GTR*N*(V&
&-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**&
&2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*&
&(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)*&
&*4)
LV1 = LOG(1-V)*(2*N*(V-1)**2*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V&
&*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W&
&+(V-1)**4)*(2*(V*W-1)*(V**5*W**5-2*V**5*W**4+2*V**3*(V**2-2*V+2)*W**3+&
&2*V**3*(2*V-3)*W**2-2*V*(2*V-3)*(2*V**2-2*V+1)*W-2*V**2+2*V-1)+CQ*(V**&
&2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-V**3*W**3+V**2*W**2-V*W+1))-2&
&*N**4*(V-1)*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**&
&4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*((V*&
&W-1)*(4*V**5*W**4-4*(V-1)*V**4*W**3+V**3*(5*V**2-6*V+9)*W**2-V*(V+1)*(&
&V**3-V**2+3*V+1)*W+(V+1)**2*(V**2-V+1))+2*CQ*(V**2*W-1)*((V-1)**2*V**2&
&*W**2-2*V*(V**2+1)*W+(V-1)**2))+2*N**2*(V-1)*VC*W*(V**3*W**3-3*(V-1)*V&
&**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4&
&*V*W+1)*((V**2*W**2-2*V*W+1)*(4*V**5*W**4-7*(V-1)*V**4*W**3+V**2*(4*V*&
&*3-8*V**2+13*V-1)*W**2+V*(4*V**4-15*V**3+9*V**2-4*V+2)*W-(V-1)*(4*V**4&
&-14*V**3+14*V**2-9*V+1))+CQ*(V-1)*(V*W-V+1)*(V**5*W**4-V**4*W**3+2*(V-&
&1)*V**2*(V**2-3*V-1)*W**2+V*(2*V**3-2*V**2+3*V-4)*W+2*V**3-2*V**2-V+2)&
&)-2*(V-1)**2*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3&
&)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*((V*W-1)*(V**5*W**4-(V-1&
&)**2*V**3*W**3-(V-1)*V**2*(5*V**2-5*V+6)*W**2+V*(4*V**4-7*V**3+2*V**2+&
&V-1)*W-(V-1)*(4*V**3-12*V**2+11*V-4))+CQ*V*(V*W-V+1)*(V**2*W**2+V*W+1)&
&*(V**2*W**2-2*V**2*W+2*V**2-2*V+1))-2*N**3*(V-1)**2*VC*(V**4*W**4-4*V*&
&*3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V&
&**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V**2*W**2-2*V**2*W+2*V**2-2*V+1&
&)*(V**4*W**4-V**3*W**3+V**2*W**2-V*W+1)+2*V*(W-1)*(V*W-1)*(V**4*W**4-V&
&**3*(2*V-1)*W**3+V**2*(2*V**2-3*V+2)*W**2-V**2*W+2*V**2-2*V+1))-4*CQ*G&
&TR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*&
&W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V*&
&*2*W**2-4*(V-1)**3*V*W+(V-1)**4)+4*CQ*GTR*N*(V-1)**2*V**2*VC*W**2*(V**&
&2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*&
&W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4))/&
&(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)
LV = LOG(V)*(N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W&
&**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*V**8*W**8-4*V**7*(5*V-2)*W&
&**7+3*V**6*(2*V-1)*(8*V-3)*W**6-2*V**5*(34*V**3-43*V**2+14*V+1)*W**5+V&
&**4*(60*V**4-80*V**3-2*V**2+46*V-15)*W**4-2*V**3*(8*V-5)*(2*V**4-6*V**&
&2+6*V-1)*W**3+2*V**2*(4*V**6+12*V**5-64*V**4+88*V**3-42*V**2+3)*W**2-2&
&*(V-1)*V*(8*V**5-16*V**4+2*V**3+22*V**2-19*V+5)*W+2*(V-1)**2*(2*V**2-4&
&*V+3)*(2*V**2-2*V+1))-4*V**8*W**8+4*V**7*(5*V-2)*W**7-V**4*(52*V**4-50&
&*V**3+17*V**2-4*V+1)*W**6+2*V**4*(42*V**4-59*V**3+24*V**2-V+2)*W**5-2*&
&V**2*(44*V**6-64*V**5-9*V**4+59*V**3-26*V**2+5*V-1)*W**4+2*V**2*(28*V*&
&*6-22*V**5-98*V**4+180*V**3-101*V**2+22*V-3)*W**3-(16*V**8+40*V**7-260&
&*V**6+380*V**5-178*V**4-16*V**3+27*V**2-6*V+1)*W**2+2*(V-1)*(16*V**6-3&
&6*V**5+2*V**4+62*V**3-57*V**2+18*V-1)*W-4*(V-1)**2*(2*V**2-5*V+4)*(2*V&
&**2-2*V+1))-N*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-&
&4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(CQ*(4*V**8*W**8-4*V**7*(4*V-1)*W**7+&
&V**6*(28*V**2-14*V+1)*W**6-2*V**5*(12*V**3-5*V**2-3*V+2)*W**5+V**4*(8*&
&V**4+12*V**3-24*V**2+14*V-1)*W**4-2*V**3*(6*V**4-6*V**3+2*V**2+3*V-2)*&
&W**3+2*V**2*(4*V**4-6*V**3+9*V**2-8*V+2)*W**2-2*(V-1)*V*(2*V**3+2*V**2&
&-3*V+1)*W+2*(V-1)**2*(2*V**2-2*V+1))-4*V**8*W**8+4*V**7*(4*V-1)*W**7-V&
&**4*(V+1)*(28*V**3-18*V**2-3*V+1)*W**6+2*V**4*(12*V**4+43*V**3-67*V**2&
&+18*V+2)*W**5-2*V**2*(4*V**6+82*V**5-137*V**4+62*V**3-5*V**2+3*V-1)*W*&
&*4+2*V**2*(62*V**5-100*V**4+30*V**3+22*V**2-5*V-3)*W**3-(32*V**7-4*V**&
&6-144*V**5+230*V**4-142*V**3+35*V**2-4*V+1)*W**2+2*(V-1)*(16*V**5-42*V&
&**4+56*V**3-37*V**2+12*V-1)*W-4*(V-1)**2*(2*V**2-2*V+1))-N**4*(V-1)*VC&
&*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*&
&*2-4*(V-1)**3*V*W+(V-1)**4)*(4*V**8*W**7-4*V**7*(V+2)*W**6+V**4*(9*V**&
&4+14*V**2+1)*W**5-V**3*(5*V**5+20*V**4+10*V**3+28*V**2+13*V+4)*W**4+2*&
&V**2*(4*V**5+23*V**4-2*V**3+28*V**2+26*V+3)*W**3-4*CQ*(V**2*W-1)*(V**2&
&*W**2-2*V*W+1)*((V-1)**2*V**2*W**2-2*V*(V**2+1)*W+(V-1)**2)-2*V*(13*V*&
&*5-10*V**4+30*V**3+4*V**2+39*V+2)*W**2+(16*V**5-31*V**4+68*V**3-42*V**&
&2+52*V+1)*W-9*V**4+24*V**3-38*V**2+28*V-13)+2*N**2*(V-1)*VC*W*(V**3*W*&
&*3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6&
&*V**2*W**2-4*V*W+1)*(2*V**7*W**6+(V-7)*V**6*W**5-CQ*(V-1)*(V*W-V+1)*(V&
&**5*W**4-V**4*W**3+2*(V-1)*V**2*(V**2-3*V-1)*W**2+V*(2*V**3-2*V**2+3*V&
&-4)*W+2*V**3-2*V**2-V+2)-V**3*(V+1)*(4*V**3-10*V**2+2*V-1)*W**4+V**2*(&
&12*V**4-35*V**3+17*V**2-7*V-1)*W**3-V*(V**6-4*V**5+14*V**4-41*V**3+27*&
&V**2-10*V+1)*W**2+(6*V**5-24*V**4+24*V**3-16*V**2+5*V+1)*W-(V-1)*(V**4&
&+V**3-11*V**2+12*V-7))+(V-1)**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**&
&2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**&
&3*V*W+(V-1)**4)*(6*V**5*W**4-V**2*(13*V**3+V**2+V-1)*W**3+2*CQ*V*(V**2&
&*W**2+V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+V*(5*V**4+17*V**3-9*V**&
&2+V-2)*W**2-(14*V**4-17*V**3+7*V**2-V-1)*W-3*V**3+3*V**2-V-1)+4*(CQ-1)&
&*GTR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**&
&2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*&
&V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)-4*(CQ-1)*GTR*N*(V-1)**2*V**2*VC*W**&
&2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*&
&(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)&
&**4))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)
LVW = (4*N**3*(V-1)**2*VC*(V**4*W**4-V**3*(2*V-1)*W**3+V**2*(2*V**2-4*&
&V+3)*W**2+V*(2*V**2-4*V+1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**&
&3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5&
&+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-4*N*(V-1)**2*VC*(V**4&
&*W**4-V**3*(2*V-1)*W**3+V**3*(2*V-1)*W**2-V*(4*V**2-5*V+2)*W+2*V**2-2*&
&V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+&
&(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**&
&2-6*V*W+1)-2*N**2*(V-1)*VC*W*(4*V**4*W**3-4*V**3*W**2-V**2*(3*V-7)*W+2&
&*(2*V**4-5*V**3+3*V**2-V-1))*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V&
&**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-&
&20*V**3*W**3+15*V**2*W**2-6*V*W+1)+2*N**4*(V-1)*VC*W*(4*V**4*W**3-V**3&
&*(3*V-7)*W**2+2*V**2*(3*V**2-4*V+7)*W-(V**2-V+1)*(2*V**2-V+3))*(V**4*W&
&**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V&
&**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-2&
&*(V-1)**2*VC*W*(V**3*W**2-V**2*(3*V-1)*W+4*V**3-8*V**2+8*V-3)*(V**4*W*&
&*4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V*&
&*6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1))*L&
&OG(1-V*W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)
LTVW = (2*N**3*(V-1)**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(CQ*(2*V**&
&2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+8*(V-1&
&)**2*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)+(V-1)**2*V**2*W*((2*V-1)*W-2&
&*(V-1)))*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2&
&-6*V*W+1)-2*N**2*(V-1)*VC*W*(4*V**4*W**3-2*(V-3)*V**3*W**2+V**2*(3*V**&
&2-7*V+12)*W+(V-2)*(V-1)**2)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V*&
&*2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-2&
&0*V**3*W**3+15*V**2*W**2-6*V*W+1)+2*N**4*(V-1)*VC*W*(4*V**4*W**3-V**3*&
&(5*V-9)*W**2+V**2*(5*V**2-10*V+13)*W-(V-1)**2*(V**2-2*V+2))*(V**4*W**4&
&-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6&
&*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-2*(V&
&-1)**2*VC*W*(V**3*W**2-V**2*(5*V-4)*W+2*(V-1)*(2*V**2-3*V+2))*(V**4*W*&
&*4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V*&
&*6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-2*&
&N*(V-1)**2*V*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*((V-1)*(2*V**3*W**3&
&-2*V**2*(3*V-2)*W**2+(V-1)*V*(8*V-5)*W-2*(V-1)**2*(2*V-1))+CQ*V*W*(V**&
&2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*W**2-2*V*(2*V-1)*W+2*V**2-2*V+1&
&))*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W&
&+1))*LOG(V*W-V+1)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)
LW = (-2*N**2*(V-1)*VC*W*(V**3*(V**2+V+2)*W**3-V*(4*V**4-6*V**3+5*V**2&
&-2*V+1)*W**2+(7*V**5-19*V**4+31*V**3-21*V**2+9*V-1)*W-(V-1)*(V**2-V+3)&
&*(3*V**2-3*V+2))*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)&
&*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1&
&)+N**3*(V-1)**2*VC*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4&
&*(V-1)**3*V*W+(V-1)**4)*(2*(2*(V-1)*V**4*W**4-V**3*(2*V-1)*(4*V-3)*W**&
&3+V**2*(12*V**3-16*V**2+4*V+3)*W**2-2*V*(V+1)*(4*V**3-8*V**2+6*V-1)*W+&
&2*(2*V**2-3*V+2)*(2*V**2-2*V+1))-CQ*(2*V**2-2*V+1)**2*(V*W-1)*(W**2-2*&
&W+2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)+(V-1)*&
&*2*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(4*V**5*&
&W**4-V**2*(15*V**3-13*V**2+3*V-1)*W**3+2*V**2*(9*V**3-12*V**2+5*V-1)*W&
&**2-(9*V**5-20*V**4+17*V**3-9*V**2+1)*W+(V-1)*(V**3-5*V**2+7*V-5))*(V*&
&*5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)+N*(V-1)**2*VC*(&
&2*(V**3*(4*V-3)*W**3-V**2*(8*V**2-10*V+5)*W**2+2*V*(2*V**2-V+1)*W-2*(2&
&*V**2-2*V+1))+CQ*(2*V**2-2*V+1)*(V*W-1)*(W**2-2*W+2))*(V**4*W**4-4*(V-&
&1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**5*W**5-&
&5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)+N**4*(V-1)*VC*W*(2*(V-1&
&)*V**4*W**3-V*(V**4+10*V**3+1)*W**2+(V**5+7*V**4-4*V**3+12*V**2+7*V+1)&
&*W-5*V**4+6*V**3-16*V**2+10*V-7)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)*&
&*2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V**3*W&
&**3-10*V**2*W**2+5*V*W-1))*LOG(W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*&
&(V*W-V+1)**4)
CVC = (N**4*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1&
&)**2)*(6*V**10*W**9-6*V**9*(3*V+2)*W**8+2*V**7*(10*V**3+2*V**2-6*V-3)*&
&W**7-V**6*(23*V**4-97*V**3+145*V**2-167*V-10)*W**6+2*V**5*(14*V**5-113&
&*V**4+297*V**3-378*V**2+106*V+5)*W**5-V**4*(13*V**6-159*V**5+668*V**4-&
&1072*V**3+475*V**2+11*V+22)*W**4-2*V**3*(14*V**6-164*V**5+452*V**4-474&
&*V**3+226*V**2-62*V+1)*W**3-V**2*(78*V**6-406*V**5+835*V**4-877*V**3+5&
&03*V**2-113*V-14)*W**2-2*(V-1)*V*(14*V**5-56*V**4+79*V**3-60*V**2+18*V&
&-1)*W-(V-1)**2*(13*V**4-33*V**3+35*V**2-11*V+2))-N**2*(V-1)*VC*W*(V**2&
&*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1&
&)*(6*V**8*W**7-6*V**7*(4*V-1)*W**6+V**4*(20*V**4+11*V**3-70*V**2+9*V-6&
&)*W**5+V**4*(4*V**4-70*V**3+207*V**2-64*V+19)*W**4-V**2*(2*V**6-42*V**&
&5+111*V**4+76*V**3-123*V**2+30*V-12)*W**3-V**2*(4*V**6+15*V**5-51*V**4&
&-46*V**3+119*V**2-61*V+26)*W**2+(V-1)*(20*V**6-58*V**5+82*V**4-91*V**3&
&+68*V**2-15*V+6)*W-(V-1)**2*(4*V**4-3*V**3-V**2+13*V-7))+3*N*(V-1)**2*&
&VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2&
&-4*V*W+1)*(4*(V-1)*V**7*W**7-V**4*(12*V**4-15*V**3+6*V**2-7*V+2)*W**6+&
&V**4*(12*V**4-9*V**3-3*V**2-11*V+7)*W**5-2*V**2*(2*V**6+4*V**5-13*V**4&
&+13*V**3-14*V**2+9*V-2)*W**4+2*(V-1)*V**2*(3*V**4-7*V**3+25*V**2-22*V+&
&5)*W**3+(V-1)*(4*V**5-36*V**4+30*V**3+5*V**2-9*V+2)*W**2+(V-1)**2*(6*V&
&**3+8*V**2-9*V+3)*W-4*(V-1)**3*V)+3*AL*N**3*(V-1)**2*VC*(V**4*W**4-4*V&
&**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*&
&V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**6-2*V**5*(2*V+1)*W**5+V*&
&*2*(8*V**3-4*V**2+4*V-1)*W**4+2*V*(2*V**2-V+1)*(2*V**3-2*V**2-2*V+1)*W&
&**3-(8*V**6-16*V**4+16*V**3-10*V**2+1)*W**2+2*(V-1)*(8*V**4-4*V**3+2*V&
&**2+2*V-1)*W-4*(V-1)*V*(2*V**2-2*V+1))-3*N**3*(V-1)**2*VC*(V**2*W**2-2&
&*(V-1)*V*W+(V-1)**2)*(8*(V-1)*V**5*W**5-V**2*(32*V**4-59*V**3+28*V**2-&
&5*V+2)*W**4+(V-1)*V*(64*V**4-101*V**3+48*V**2-13*V+4)*W**3-(V-1)*(80*V&
&**5-180*V**4+156*V**3-59*V**2+13*V-2)*W**2+(V-1)**2*(56*V**4-106*V**3+&
&84*V**2-21*V+3)*W-8*(V-1)**3*V*(2*V**2-3*V+2))*(V**6*W**6-6*V**5*W**5+&
&15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-4*CQ*N**3*(V-1)**2*(2*&
&V**2-2*V+1)**2*VC*W*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-&
&4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W&
&**3+15*V**2*W**2-6*V*W+1)+4*CQ*N*(V-1)**2*(2*V**2-2*V+1)*VC*W*(V**4*W*&
&*4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V*&
&*6*W**6-6*V**5*W**5+15*V**4*W**4-20*V**3*W**3+15*V**2*W**2-6*V*W+1)-4*&
&GTR*N**3*(V-1)*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6&
&*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*((V-3)*(V-1)*V**6*W**5-V*&
&*4*(2*V**4-13*V**3+21*V**2-7*V+1)*W**4-V**3*(8*V**4-37*V**3+42*V**2-21&
&*V-4)*W**3-V**2*(24*V**4-61*V**3+56*V**2-13*V+6)*W**2-2*V*(4*V**4-7*V*&
&*3+V**2+2*V-2)*W-(V-1)**2*(2*V**2-2*V+1))-3*(V-1)**2*VC*W*(V**3*W**3-3&
&*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**4*W**4-4*V**3*W**3+6*V**&
&2*W**2-4*V*W+1)*(2*V**6*W**5+V**3*(2*V**3-V**2-V-2)*W**4-(V-1)*V**2*(7&
&*V**3+8*V+2)*W**3+V*(3*V**5-8*V**4+20*V**3-13*V**2-6*V+2)*W**2+(2*V**5&
&-7*V**4-V**3+12*V**2-2*V-2)*W+(V-1)*(3*V**3-7*V**2+5*V-3))-3*AL*N*(V-1&
&)**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-&
&1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*V**6*W**&
&5-2*V**5*(2*V+1)*W**4+V**2*(4*V**4+2*V**2+2*V-1)*W**3-2*V*(2*V**4-2*V*&
&*3+2*V**2+V-1)*W**2-(4*V**3-8*V**2+2*V+1)*W+2*(V-1)*(2*V**2-1))+4*GTR*&
&N*(V-1)**3*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-&
&1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**6*W**5-V**4*(2*V**2-4*V+1&
&)*W**4+2*V**3*(V**2+11*V+2)*W**3-2*V**2*(18*V**2-7*V+3)*W**2+V*(2*V**2&
&-7*V+4)*W-2*V**2+2*V-1)-6*AL*N**2*(V-1)**2*V*VC*W*(V**2*W**2+V*W+1)*(V&
&**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V&
&*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W&
&+(V-1)**4)+6*AL*(V-1)**2*V*VC*W*(V**2*W**2+V*W+1)*(V**2*W**2-2*V**2*W+&
&2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4&
&*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)+12*AL*G&
&TR*N**3*(V-1)**2*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*&
&W**2-2*V**2*W+2*V**2-2*V+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V*&
&*2*W**2-4*(V-1)**3*V*W+(V-1)**4)-12*AL*GTR*N*(V-1)**2*V**2*VC*W**2*(V*&
&*2*W**2+1)*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(V**4&
&*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4))&
&/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**6*(V*W-V+1)**4)/3.D0
LM = LOG(S/M**2)*(-N**4*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(4*V**8*&
&W**7-4*V**7*(V+2)*W**6+V**4*(9*V**4-8*V**3+22*V**2+1)*W**5-V**3*(5*V**&
&5+12*V**4+6*V**3+32*V**2+5*V+4)*W**4+2*V**2*(2*V**5+19*V**4+4*V**3+20*&
&V**2+10*V+3)*W**3-2*V*(7*V**5+2*V**4+20*V**3+8*V**2+15*V+2)*W**2+(4*V*&
&*5+9*V**4+16*V**3-2*V**2+20*V+1)*W-5*V**4+8*V**3-14*V**2+8*V-5)-N**3*(&
&V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*W&
&**6-2*V**5*(2*V+1)*W**5+V**2*(8*V**4-8*V**3+12*V**2-4*V+1)*W**4-2*V*(4&
&*V**5-2*V**4+6*V**2-3*V+1)*W**3+(8*V**6-8*V**4+16*V**3-2*V**2+1)*W**2-&
&2*(4*V**2-2*V+1)*(2*V**3-2*V**2+V+1)*W+4*(V**2-V+1)*(2*V**2-2*V+1))+N*&
&(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**6*&
&W**6-2*V**5*(2*V+1)*W**5+V**2*(4*V**4+6*V**2-2*V+1)*W**4-2*V*(2*V**4+2&
&*V**3+2*V**2-V+1)*W**3+(2*V+1)*(4*V**3+1)*W**2-2*(6*V**3-2*V**2+V+1)*W&
&+4*(2*V**2-2*V+1))+2*N**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*&
&V*W+(V-1)**2)*(2*V**6*W**5-V**5*(V+5)*W**4+V**2*(3*V**4-2*V**3+9*V**2-&
&V+1)*W**3-V*(V**5+7*V**4-6*V**3+9*V**2+V+2)*W**2+(4*V**5-3*V**4+8*V**3&
&-3*V**2+5*V+1)*W-V**4+V**3-5*V**2+4*V-3)-(V-1)*VC*W*(V**2*W**2-2*V*W+1&
&)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(2*V**5*W**4-V**2*(3*V**3-V**2-V+1)&
&*W**3+(V-1)*V*(3*V**3-6*V**2-V-2)*W**2+(6*V**4-9*V**3+3*V**2-V-1)*W+3*&
&V**3-3*V**2+V+1)-4*GTR*N**3*(V-1)*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**&
&2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+4*GTR*N*(V-1&
&)*V**2*VC*W**2*(V**2*W**2+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W*&
&*2-2*V**2*W+2*V**2-2*V+1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**4*(V*W-V+&
&1)**2)
LMP = LOG(S/MP**2)*(-N**4*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*&
&W+1)*(4*V**6*W**5-4*V**5*(3*V-4)*W**4+8*V**4*(2*V**2-5*V+4)*W**3-2*(V-&
&1)*V**2*(7*V**3-15*V**2+17*V-1)*W**2+(V-1)**2*V*(7*V**3-14*V**2+19*V-4&
&)*W-(V-1)**4*(V**2-2*V+2))+2*N**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W&
&**2-4*V*W+1)*(2*V**6*W**5-2*V**5*(3*V-4)*W**4+V**4*(9*V**2-22*V+17)*W*&
&*3-(V-1)*V**2*(10*V**3-21*V**2+21*V-2)*W**2+(V-1)**2*V*(6*V**3-13*V**2&
&+15*V-4)*W-(V-1)**4*(V**2-2*V+2))-N**3*(V-1)*VC*(2*V**2*W**2-2*V*(2*V-&
&1)*W+2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W&
&**4-4*(V-1)*V**3*W**3+8*(V-1)**2*V**2*W**2-8*(V-1)**3*V*W+4*(V-1)**4)-&
&(V-1)**2*VC*W*(2*V**4*W**3-2*V**2*(3*V**2-3*V+1)*W**2+(V-1)*V*(5*V**2-&
&7*V+4)*W-(V-1)**2*(V**2-2*V+2))*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V&
&*W+1)+N*(V-1)*V**2*VC*W**2*(V**2*W**2-2*(V-1)*V*W+2*(V-1)**2)*(2*V**2*&
&W**2-2*V*(2*V-1)*W+2*V**2-2*V+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*&
&V*W+1))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)**4*(V*W-V+1)**2)
STRUV13=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV14(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV14
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*&
&V*W+(V-1)**2)*(CQ*(4*V**12*W**10-4*V**10*(5*V**2-2*V+2)*W**9+4*V**6*(1&
&3*V**6-13*V**5+19*V**4-14*V**3+12*V**2-6*V+2)*W**8-4*V**6*(23*V**6-37*&
&V**5+67*V**4-62*V**3+36*V**2-6*V+2)*W**7+V**4*(120*V**8-227*V**7+440*V&
&**6-378*V**5+45*V**4+240*V**3-192*V**2+96*V-24)*W**6-V**6*(112*V**6-17&
&8*V**5+275*V**4-24*V**3-413*V**2+510*V-170)*W**5+V**2*(64*V**10-7*V**9&
&-194*V**8+750*V**7-1245*V**6+1020*V**5-228*V**4-240*V**3+240*V**2-120*&
&V+24)*W**4-V**2*(16*V**10+96*V**9-373*V**8+838*V**7-1161*V**6+1004*V**&
&5-596*V**4+368*V**3-272*V**2+120*V-24)*W**3+(V-1)*(48*V**10-48*V**9-29&
&*V**8+325*V**7-662*V**6+793*V**5-563*V**4+208*V**3+8*V**2-40*V+8)*W**2&
&-(V-1)**2*(48*V**8-128*V**7+275*V**6-392*V**5+457*V**4-358*V**3+194*V*&
&*2-64*V+16)*W+16*(V-1)**3*(V**6-3*V**5+7*V**4-11*V**3+13*V**2-9*V+3))-&
&8*V**12*W**10+20*V**10*(2*V**2-V+1)*W**9-4*V**6*(26*V**6-31*V**5+34*V*&
&*4-8*V**3+9*V**2-6*V+2)*W**8+2*V**6*(92*V**6-171*V**5+184*V**4-14*V**3&
&-23*V**2+36*V-12)*W**7-V**4*(240*V**8-529*V**7+492*V**6+166*V**5-337*V&
&**4+372*V**3-236*V**2+96*V-24)*W**6+V**4*(224*V**8-444*V**7+113*V**6+7&
&86*V**5-607*V**4-12*V**3+452*V**2-384*V+96)*W**5-V**2*(128*V**10-73*V*&
&*9-726*V**8+1934*V**7-1699*V**6+552*V**5+432*V**4-672*V**3+348*V**2-12&
&0*V+24)*W**4+V**2*(32*V**10+176*V**9-951*V**8+1822*V**7-1601*V**6+976*&
&V**5-932*V**4+1240*V**3-1210*V**2+600*V-120)*W**3-(V-1)*(96*V**10-160*&
&V**9-191*V**8+1091*V**7-1554*V**6+1303*V**5-565*V**4+64*V**3+44*V**2-4&
&0*V+8)*W**2+(V-1)**2*(96*V**8-336*V**7+621*V**6-650*V**5+573*V**4-432*&
&V**3+368*V**2-192*V+48)*W-32*(V-1)**3*(V**2-2*V+2)*(V**2-V+1)**2)+2*N*&
&*2*(V-1)*V**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2&
&)*(CQ*(2*V**10*W**9-2*V**8*(5*V**2-2*V+2)*W**8+V**8*(22*V**2-13*V+13)*&
&W**7-V**6*(26*V**4-11*V**3+15*V**2-8*V+4)*W**6+2*V**6*(8*V**4+5*V**3+3&
&*V**2-16*V+8)*W**5-V**4*(4*V**6+21*V**5+10*V**4-72*V**3+61*V**2-30*V+1&
&0)*W**4+(V-1)*V**4*(9*V**4+32*V**3-35*V**2+6*V-3)*W**3-(V-1)**2*V**2*(&
&13*V**4+15*V**3-19*V**2+8*V-4)*W**2+(V-1)**3*V**2*(9*V**2+2*V-2)*W-(V-&
&1)**4*(3*V**2-2*V+2))-(V*W-1)*(V*W-V+1)*(4*V**8*W**7-4*V**6*(4*V**2-3*&
&V+3)*W**6+V**6*(28*V**2-31*V+31)*W**5-V**4*(24*V**4-15*V**3-V**2+32*V-&
&16)*W**4+V**4*(8*V**4+20*V**3-41*V**2+42*V-21)*W**3-2*(V-1)*V**2*(8*V*&
&*4+5*V**3-16*V**2+22*V-11)*W**2+4*(V-1)**2*V**2*(4*V**2-3*V+3)*W-2*(V-&
&1)**3*(2*V**2-V+1)))-(V-1)**2*V**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2&
&-2*(V-1)*V*W+(V-1)**2)*(4*V**8*W**8+CQ*(2*V**8*W**7-2*V**6*(V+2)*(3*V-&
&2)*W**6+V**6*(7*V**2+27*V-27)*W**5-V**4*(4*V**4+37*V**3-43*V**2+12*V-6&
&)*W**4+V**4*(V**4+25*V**3-27*V**2+4*V-2)*W**3-(V-1)*V**4*(7*V**2+10*V-&
&10)*W**2+(V-1)**2*V**2*(7*V**2+V-1)*W-(V-1)**3*(3*V**2-2*V+2))-18*V**8&
&*W**7+4*V**6*(8*V**2+5*V-5)*W**6-V**6*(29*V**2+63*V-63)*W**5+V**4*(14*&
&V**4+79*V**3-73*V**2-12*V+6)*W**4-V**4*(3*V**4+49*V**3-27*V**2-44*V+22&
&)*W**3+(V-1)*V**2*(V**2+2*V-2)*(13*V**2+2*V-2)*W**2-(V-1)**2*V**2*(13*&
&V**2+5*V-5)*W+(V-1)**3*(5*V**2-2*V+2)))/(N**2*(V-1)**3*V**3*W**2*(V*W-&
&1)**5*(V*W-V+1)**5)
LV1 = LOG(1-V)*(N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V&
&**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*(V-1)*(V*W-&
&V+1)*(4*V**8*W**7-8*V**7*(V+1)*W**6+4*V**6*(2*V**2+3*V+2)*W**5-4*V**3*&
&(4*V**4+3*V**3+5*V**2-8*V+4)*W**4+V**2*(16*V**4+4*V**3+69*V**2-112*V+4&
&8)*W**3-2*V*(4*V**4-V**3+50*V**2-72*V+24)*W**2+(2*V**4-2*V**3+65*V**2-&
&80*V+16)*W-16*(V-1))+2*(V**2*W**2-2*V*W+1)*(4*V**8*W**6-16*V**6*(V**2-&
&V+1)*W**5+V**5*(28*V**3-47*V**2+48*V+7)*W**4-2*V**2*(12*V**6-26*V**5+3&
&6*V**4-3*V**3-4*V**2+13*V-4)*W**3+V**3*(8*V**5-11*V**4+29*V**3+7*V+7)*&
&W**2-(8*V**7-3*V**6+7*V**5+17*V**4-37*V**3+50*V**2-34*V+8)*W+8*(V-1)*(&
&V+1)*(V**2-V+1)**2))-2*N**2*(V-1)*V*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**&
&4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(CQ*&
&(V-1)*V*(V*W-V+1)*(2*V**6*W**6-4*V**5*(V+1)*W**5+V**4*(4*V**2+6*V+5)*W&
&**4-V**3*(8*V**2+8*V+3)*W**3+2*V**2*(5*V**2+V+2)*W**2-V*(6*V**2-2*V+3)&
&*W+2*V**2-2*V+1)+(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**6*W**4-4*V**4*(&
&3*V**2-3*V+2)*W**3+V**3*(16*V**3-33*V**2+37*V-12)*W**2-(V-1)*V*(8*V**4&
&-14*V**3+16*V**2-V+1)*W+(V-1)**4*(V+1)))+(V-1)**2*V**2*VC*W*(V**2*W**2&
&-2*V*W+1)*(2*(V**2*W-(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1)+CQ*(V**&
&2*W**2-2*V**2*W+2*V**2-2*V+1)*(2*V**2*W**2-2*V*W+1))*(V**5*W**5-5*(V-1&
&)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W&
&-(V-1)**5))/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LV = LOG(V)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W&
&+(V-1)**2)*(CQ*(4*V**12*W**10-4*V**10*(5*V**2-2*V+2)*W**9+4*V**10*(11*&
&V**2-7*V+7)*W**8-4*V**6*(13*V**6-7*V**5+7*V**4+8*V**3-24*V**2+24*V-8)*&
&W**7+V**6*(32*V**6+13*V**5+8*V**4+22*V**3-171*V**2+192*V-64)*W**6-V**4&
&*(8*V**8+38*V**7+35*V**6-160*V**5+211*V**4-426*V**3+590*V**2-384*V+96)&
&*W**5+(V-1)*V**4*(17*V**6+87*V**5-107*V**4+136*V**3-308*V**2+288*V-96)&
&*W**4-(V-1)**2*V**2*(35*V**6+84*V**5-188*V**4+304*V**3-392*V**2+288*V-&
&96)*W**3+(V-1)**3*V**4*(59*V**2-45*V+45)*W**2-(V-1)**4*(51*V**4-98*V**&
&3+130*V**2-64*V+32)*W+16*(V-1)**5*(V**2-2*V+2))-8*V**12*W**10+16*V**10&
&*(3*V**2-2*V+2)*W**9-V**6*(144*V**6-199*V**5+219*V**4-44*V**3+32*V**2-&
&12*V+4)*W**8+4*V**6*(68*V**6-131*V**5+135*V**4+2*V**3-26*V**2+30*V-10)&
&*W**7-V**4*(344*V**8-741*V**7+573*V**6+460*V**5-552*V**4+420*V**3-196*&
&V**2+48*V-12)*W**6+2*V**4*(144*V**8-265*V**7-62*V**6+729*V**5-486*V**4&
&-39*V**3+321*V**2-264*V+66)*W**5-V**2*(144*V**10-34*V**9-1130*V**8+267&
&1*V**7-2073*V**6+537*V**5+437*V**4-600*V**3+240*V**2-60*V+12)*W**4+2*V&
&**2*(16*V**10+104*V**9-569*V**8+1019*V**7-730*V**6+331*V**5-437*V**4+7&
&12*V**3-718*V**2+360*V-72)*W**3-2*(V-1)*(48*V**10-80*V**9-165*V**8+730&
&*V**7-997*V**6+846*V**5-422*V**4+108*V**3-12*V**2-10*V+2)*W**2+2*(V-1)&
&**2*(48*V**8-184*V**7+338*V**6-351*V**5+309*V**4-233*V**3+199*V**2-104&
&*V+26)*W-16*(V-1)**3*(V**2-V+1)**2*(2*V**2-5*V+5))+2*N**2*(V-1)*VC*W*(&
&V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(CQ*V**2*(2*V**10*&
&W**9-2*V**8*(5*V**2-2*V+2)*W**8+V**8*(22*V**2-13*V+13)*W**7-V**6*(26*V&
&**4-11*V**3+15*V**2-8*V+4)*W**6+2*V**6*(8*V**4+5*V**3+3*V**2-16*V+8)*W&
&**5-V**4*(4*V**6+21*V**5+10*V**4-72*V**3+61*V**2-30*V+10)*W**4+(V-1)*V&
&**4*(9*V**4+32*V**3-35*V**2+6*V-3)*W**3-(V-1)**2*V**2*(13*V**4+15*V**3&
&-19*V**2+8*V-4)*W**2+(V-1)**3*V**2*(9*V**2+2*V-2)*W-(V-1)**4*(3*V**2-2&
&*V+2))-(V*W-1)*(V*W-V+1)*(4*V**10*W**7-4*V**8*(4*V**2-3*V+3)*W**6+V**4&
&*(28*V**6-30*V**5+27*V**4+8*V**3-9*V**2+6*V-2)*W**5-V**4*(24*V**6-15*V&
&**5-2*V**4+40*V**3-35*V**2+18*V-6)*W**4+V**2*(8*V**8+19*V**7-35*V**6+3&
&2*V**5-12*V**4-16*V**3+24*V**2-16*V+4)*W**3-(V-1)*V**2*(16*V**6+13*V**&
&5-37*V**4+56*V**3-48*V**2+24*V-8)*W**2+(V-1)**2*(18*V**6-13*V**5+14*V*&
&*4-5*V**2+6*V-2)*W-(V-1)**3*(6*V**4-5*V**3+7*V**2-4*V+2)))-(V-1)**2*V*&
&*2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(8*V**8*W&
&**8+CQ*(2*V**8*W**7-2*V**6*(V+2)*(3*V-2)*W**6+V**6*(7*V**2+27*V-27)*W*&
&*5-V**4*(4*V**4+37*V**3-43*V**2+12*V-6)*W**4+V**4*(V**4+25*V**3-27*V**&
&2+4*V-2)*W**3-(V-1)*V**4*(7*V**2+10*V-10)*W**2+(V-1)**2*V**2*(7*V**2+V&
&-1)*W-(V-1)**3*(3*V**2-2*V+2))-V**6*(35*V**2-2*V+2)*W**7+2*V**6*(31*V*&
&*2+12*V-12)*W**6-3*V**4*(19*V**4+30*V**3-32*V**2+4*V-2)*W**5+2*V**6*(1&
&4*V**2+61*V-61)*W**4-V**2*(6*V**6+80*V**5-49*V**4-68*V**3+49*V**2-18*V&
&+6)*W**3+2*(V-1)*V**2*(11*V**4+23*V**3-25*V**2+4*V-2)*W**2-2*(V-1)**2*&
&(11*V**4+2*V**3-3*V**2+2*V-1)*W+4*(V-1)**3*(2*V**2-V+1)))/(N**2*(V-1)*&
&*3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LVW = (-2*N**4*(V-1)*VC*(4*V**6*W**4-8*V**4*(V**2-V+2)*W**3+V**3*(8*V*&
&*3-17*V**2+26*V+15)*W**2+(6*V**5-24*V**4-11*V**3+33*V**2-52*V+16)*W+16&
&*(V**2-V+1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W&
&-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**&
&2*W**2+5*(V-1)**4*V*W-(V-1)**5)+2*N**2*(V-1)*V*VC*W*(4*V**5*W**3-2*V**&
&3*(4*V**2-5*V+5)*W**2+V**3*(8*V**2-15*V+15)*W-(V-1)*(2*V**3-V**2-V+1))&
&*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-&
&5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)*&
&*4*V*W-(V-1)**5)+2*(V-1)**2*V**2*VC*W*(2*V**2*W**2-3*V**2*W+2*V**2-2*V&
&+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W*&
&*5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-&
&1)**4*V*W-(V-1)**5))*LOG(1-V*W)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V&
&*W-V+1)**5)
LTVW = (-2*N**4*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**5*W**5-5&
&*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(CQ*(4*V**8*W**7-8*V**7*&
&(2*V-1)*W**6+4*V**6*(7*V**2-7*V+2)*W**5-4*(V-1)*V**3*(8*V**4-5*V**3+5*&
&V**2-8*V+4)*W**4+(V-1)**2*V**2*(25*V**4+2*V**3+21*V**2-80*V+48)*W**3-2&
&*(V-1)**3*V*(5*V**4+21*V**3-22*V**2-24*V+24)*W**2+(V-1)**4*(V**4+48*V*&
&*3-79*V**2+16*V+16)*W-16*(V-1)**7)-(V-1)*V*W*(V**5*W**4-V**3*(3*V**3+3&
&*V**2-4)*W**3+(V-1)*V**2*(2*V**3+5*V**2+15*V-12)*W**2+(V-1)**2*V*(5*V*&
&*3-3*V**2-16*V+12)*W-4*(V-1)**3*(V+1)*(V**2-V+1)))+2*N**2*(V-1)*V*VC*W&
&*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-&
&10*V**2*W**2+5*V*W-1)*(2*CQ*V*(2*V**6*W**6-4*V**5*(2*V-1)*W**5+V**4*(1&
&5*V**2-16*V+5)*W**4-(V-1)*V**3*(19*V**2-14*V+3)*W**3+2*(V-1)**2*V**2*(&
&8*V**2-5*V+2)*W**2-(V-1)**3*V*(7*V**2-4*V+3)*W+(V-1)**4*(V**2+1))+(V-1&
&)*(V*W-V+1)*(2*V**5*W**4-V**4*(5*V-6)*W**3+V**2*(11*V**3-14*V**2-4*V+1&
&)*W**2-(V-1)*V*(9*V**3-7*V**2+2)*W+(V-1)**2*(V+1)*(V**2-V+1)))+2*(V-1)&
&**2*V**2*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V*W*(2*V**4*W**4-V**3*&
&(7*V-4)*W**3+V**2*(10*V**2-13*V+5)*W**2-(V-1)*V*(7*V**2-7*V+2)*W+(V-1)&
&**2*(2*V**2-V+1))-CQ*(V-1)*(V**2*W**2-2*V**2*W+V**2+1)*(2*V**2*W**2-2*&
&(V-1)*V*W+(V-1)**2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+&
&5*V*W-1))*LOG(V*W-V+1)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**&
&5)
LW = (-(V-1)**2*V**2*VC*W*(4*V**2*W**2-(9*V**2-2*V+2)*W+5*V**2-6*V+6)*&
&(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5&
&*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**&
&4*V*W-(V-1)**5)-N**4*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W&
&+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(&
&V-1)**4)*(4*V**6*(2*V**2-3*V+3)*W**5-V**2*(40*V**6-89*V**5+97*V**4-12*&
&V**3-4*V**2+12*V-4)*W**4+2*V**2*(40*V**6-104*V**5+121*V**4-14*V**3-43*&
&V**2+60*V-20)*W**3-8*CQ*(V**2-V+1)**3*(V*W-1)*(V*W-V+1)*(W**2-2*W+2)-(&
&80*V**8-223*V**7+308*V**6-174*V**5+101*V**4-28*V**3+28*V**2-16*V+4)*W*&
&*2+(32*V**8-64*V**7+33*V**6+92*V**5-85*V**4-54*V**3+186*V**2-144*V+36)&
&*W-16*(V-1)*(V**2-V+1)**2*(2*V**2-3*V+3))-2*N**2*(V-1)**2*VC*W*(V**2*(&
&V**2-2*V+2)*(V**2-V+1)*W**3-(V-1)*V**2*(V**2+4*V-4)*W**2+(V**2+2*V-2)*&
&(V**4-V**3+2*V**2-2*V+1)*W-(V-1)*(3*V**4-4*V**3+6*V**2-4*V+2))*(V**4*W&
&**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V&
&-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4))*LOG(W)/(N**2*(V-1)**3*V**3*&
&W**2*(V*W-1)**5*(V*W-V+1)**5)
CVC = (3*N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1&
&)**2)*(8*V**12*W**10-8*V**10*(3*V**2+2*V-2)*W**9+V**6*(24*V**6+91*V**5&
&-141*V**4+148*V**3-194*V**2+144*V-48)*W**8+2*V**6*(12*V**6-125*V**5+23&
&0*V**4-306*V**3+393*V**2-288*V+96)*W**7-V**4*(128*V**8-494*V**7+927*V*&
&*6-1036*V**5+799*V**4+66*V**3-694*V**2+576*V-144)*W**6+V**4*(192*V**8-&
&522*V**7+773*V**6-186*V**5-1129*V**4+2676*V**3-2908*V**2+1728*V-432)*W&
&**5-V**2*(128*V**10-107*V**9-422*V**8+1971*V**7-3646*V**6+4107*V**5-24&
&61*V**4+72*V**3+1062*V**2-720*V+144)*W**4+V**2*(32*V**10+192*V**9-891*&
&V**8+1964*V**7-2315*V**6+1082*V**5+1366*V**4-3208*V**3+2962*V**2-1440*&
&V+288)*W**3-(V-1)*(96*V**10-96*V**9-217*V**8+1119*V**7-1982*V**6+2191*&
&V**5-1393*V**4+280*V**3+290*V**2-240*V+48)*W**2+(V-1)**2*(96*V**8-256*&
&V**7+427*V**6-302*V**5+3*V**4+312*V**3-328*V**2+192*V-48)*W-32*(V-1)**&
&3*(V**2-V+1)**3)-6*AL*N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(&
&V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*&
&V*W-(V-1)**5)*(4*(V-1)*V**8*W**7-8*(V-1)*V**7*(V+1)*W**6+4*V**3*(6*V**&
&5-11*V**4+11*V**3-12*V**2+6*V-2)*W**5+4*V**2*(4*V**7-10*V**6+10*V**5+5&
&*V**4-15*V**3+24*V**2-14*V+6)*W**4-V*(16*V**8-40*V**6+128*V**5-113*V**&
&4+81*V**3+16*V**2-24*V+24)*W**3+2*(24*V**8-48*V**7+72*V**6-29*V**5-5*V&
&**4+46*V**3-24*V**2+12*V+4)*W**2-(48*V**7-128*V**6+238*V**5-236*V**4+1&
&73*V**3-47*V**2+16)*W+16*(V**2-V+1)**3)-6*N**2*(V-1)*V**2*VC*W*(V**3*W&
&**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(&
&V-1)**3)*(4*V**8*W**7-8*V**6*(V**2+V-1)*W**6+V**4*(4*V**4+39*V**3-60*V&
&**2+42*V-21)*W**5-3*(V-1)*V**4*(15*V**2-11*V+11)*W**4+3*(V-1)*V**2*(7*&
&V**4-6*V**3-8*V**2+28*V-14)*W**3-(V-1)*V**2*(7*V**4-18*V**3+12*V**2+12&
&*V-6)*W**2+(V-1)**3*(10*V**2-21*V+21)*W-(V-1)**3*(3*V**2-5*V+5))+3*(V-&
&1)**2*V**2*VC*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(&
&V**6*(11*V**2-4*V+4)*W**7-2*V**6*(15*V**2-11*V+11)*W**6+V**4*(30*V**4-&
&23*V**3+11*V**2+24*V-12)*W**5-V**4*(14*V**4-23*V**3+89*V**2-132*V+66)*&
&W**4+V**2*(3*V**6-23*V**5+144*V**4-254*V**3+157*V**2-36*V+12)*W**3+(V-&
&1)*V**2*(5*V**4-44*V**3+54*V**2-20*V+10)*W**2+(V-1)**2*(5*V**4-5*V**3+&
&V**2+8*V-4)*W-(V-1)**3*(V**2-2*V+2))+12*AL*N**2*(V-1)**2*V**2*VC*W*(V*&
&*2*W**2-2*V*W+1)*(V**2*W**2-V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(&
&2*V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**&
&3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-6*AL*(V-1)**2*V**2*VC&
&*W*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)*(2*V**2*W**2-&
&2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**&
&3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+4*CQ*N**4*(V-1)*VC*W*(V**3*W**3-3&
&*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(12*V**6*W**4-24*(V-1)*V**5*&
&W**3+8*V**2*(2*V**6-6*V**5+15*V**4-20*V**3+15*V**2-6*V+2)*W**2-4*(V-1)&
&*V*(8*V**6-24*V**5+51*V**4-62*V**3+51*V**2-24*V+8)*W+(V-1)**2*(16*V**6&
&-48*V**5+99*V**4-118*V**3+99*V**2-48*V+16))*(V**5*W**5-5*V**4*W**4+10*&
&V**3*W**3-10*V**2*W**2+5*V*W-1)-24*CQ*N**2*(V-1)*V**2*VC*W*(V**2*W**2-&
&(V-1)*V*W+(V-1)**2)*(2*V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W**3-3*(V&
&-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V**5*W**5-5*V**4*W**4+10*V**3*&
&W**3-10*V**2*W**2+5*V*W-1)+12*CQ*(V-1)**3*V**2*VC*W*(2*V**2*W**2-2*(V-&
&1)*V*W+(V-1)**2)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)&
&*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1))/(N**2*(V-1&
&)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)/6.D0
LM = LOG(S/M**2)*(N**4*VC*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-&
&(V-1)**3)*(4*V**9*W**7-4*V**7*(2*V**2+V+4)*W**6+V**3*(16*V**6-15*V**5+&
&55*V**4+20*V**3+32*V**2-12*V+4)*W**5-V**2*(16*V**7-3*V**6+26*V**5+77*V&
&**4+24*V**3+108*V**2-40*V+12)*W**4+V*(16*V**8+11*V**6+80*V**5+54*V**4+&
&23*V**3+164*V**2-48*V+12)*W**3-(48*V**8-96*V**7+195*V**6-128*V**5+213*&
&V**4-96*V**3+164*V**2-24*V+4)*W**2+(48*V**7-128*V**6+275*V**5-297*V**4&
&+319*V**3-177*V**2+108*V-4)*W-16*(V**2-V+1)**2*(V**2-V+2))-2*N**2*VC*W&
&*(V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(2*V**8&
&*W**5-4*V**6*(V**2+1)*W**4+V**2*(4*V**6+V**5-V**4+15*V**3-9*V**2+6*V-2&
&)*W**3-V*(V**2-V+1)*(7*V**4+5*V**3-2*V**2+6*V-4)*W**2+(10*V**6-19*V**5&
&+22*V**4-12*V**3+3*V**2+2*V-2)*W-(V-1)*(3*V**4-5*V**3+6*V**2-4*V+2))+(&
&V-1)*V**2*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(&
&V**3*(V**2+2)*W**4-V**2*(V**3+5*V**2-2*V+6)*W**3+V*(7*V**3-3*V**2+3*V+&
&6)*W**2-(7*V**3-7*V**2+6*V+2)*W+3*V**2-4*V+3))/(N**2*(V-1)**2*V**3*W**&
&2*(V*W-1)**3*(V*W-V+1)**3)
LMP = LOG(S/MP**2)*(N**4*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(4*V**9*W*&
&*7-4*V**7*(7*V**2-9*V+4)*W**6+4*V**6*(23*V**3-56*V**2+51*V-16)*W**5-4*&
&(V-1)*V**3*(45*V**5-110*V**4+101*V**3-27*V**2-8*V+4)*W**4+(V-1)**2*V**&
&2*(224*V**5-535*V**4+498*V**3-123*V**2-80*V+48)*W**3-2*(V-1)**3*V*(88*&
&V**5-203*V**4+205*V**3-70*V**2-24*V+24)*W**2+(V-1)**4*(80*V**5-175*V**&
&4+208*V**3-111*V**2+16*V+16)*W-16*(V-1)**5*(V**2-V+1)**2)-2*N**2*V**2*&
&VC*W*(V*W-V+1)*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**6*W**5-4*V**4*(2*&
&V**2-2*V+1)*W**4+V**3*(2*V-1)*(7*V**2-10*V+7)*W**3-6*(V-1)*V**2*(2*V-1&
&)*(V**2-V+1)*W**2+2*(V-1)**2*V*(2*V**3-V**2+2*V+1)*W-(V-1)**3*(V**2+1)&
&)+(V-1)**2*V**2*VC*W*(V**2*W**2-2*V**2*W+V**2+1)*(2*V**2*W**2-2*(V-1)*&
&V*W+(V-1)**2)*(V**3*W**3-3*V**2*W**2+3*V*W-1))/(N**2*(V-1)**2*V**3*W**&
&2*(V*W-1)**3*(V*W-V+1)**3)
STRUV14=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV15(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV15
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-16*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*(4*V**14*W**13-2*V&
&**12*(10*V**2+V-1)*W**12+V**8*(52*V**6-V**5+7*V**4-14*V**3+12*V**2-6*V&
&+2)*W**11-V**8*(92*V**6-36*V**5+71*V**4-80*V**3+65*V**2-30*V+10)*W**10&
&+V**6*(110*V**8-59*V**7+75*V**6-22*V**5-22*V**4+62*V**3-58*V**2+32*V-8&
&)*W**9-V**6*(90*V**8-14*V**7-115*V**6+350*V**5-437*V**4+404*V**3-284*V&
&**2+128*V-32)*W**8+V**4*(62*V**10+V**9-166*V**8+340*V**7-183*V**6-30*V&
&**5+122*V**4-168*V**3+132*V**2-60*V+12)*W**7-2*CQ*(V**2*W-V+1)*(V**2*W&
&**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**6*(V**2-2*V+2)*W**6-&
&V**4*(V**4-2*V**3+4*V**2-4*V+2)*W**5-2*V**4*(V**4-4*V**3+5*V**2-2*V+1)&
&*W**4+2*V**2*(V**6-4*V**5+7*V**4-8*V**3+9*V**2-6*V+2)*W**3-2*(V-1)*V**&
&2*(V**4-4*V**3+5*V**2-2*V+1)*W**2-(V-1)**2*(V**4-2*V**3+4*V**2-4*V+2)*&
&W+(V-1)**3*(V**2-2*V+2))-V**4*(42*V**10-14*V**9-88*V**8+45*V**7+481*V*&
&*6-937*V**5+975*V**4-784*V**3+466*V**2-180*V+36)*W**6+V**2*(20*V**12+3&
&1*V**11-214*V**10+344*V**9-127*V**8-2*V**7-166*V**6+328*V**5-370*V**4+&
&280*V**3-144*V**2+48*V-8)*W**5-V**2*(4*V**12+50*V**11-190*V**10+241*V*&
&*9-22*V**8-71*V**7-221*V**6+592*V**5-739*V**4+570*V**3-290*V**2+96*V-1&
&6)*W**4+(V-1)*(16*V**12+14*V**11-206*V**10+535*V**9-696*V**8+679*V**7-&
&569*V**6+436*V**5-301*V**4+150*V**3-52*V**2+12*V-2)*W**3-(V-1)**2*(24*&
&V**10-58*V**9+32*V**8+109*V**7-202*V**6+189*V**5-77*V**4+15*V**2-10*V+&
&2)*W**2+(V-1)**3*(16*V**8-58*V**7+113*V**6-128*V**5+115*V**4-78*V**3+5&
&4*V**2-24*V+6)*W-4*(V-1)**4*(V**2-2*V+2)*(V**2-V+1)**2)-8*GTR*N**2*(V-&
&1)**2*VC*W*(V*W-1)*(V*W-V+1)*(CQ*(2*V**8*(2*V**4-3*V**3+5*V**2-4*V+2)*&
&W**10-V**8*(19*V**4-30*V**3+50*V**2-40*V+20)*W**9+V**6*(42*V**6-65*V**&
&5+83*V**4-20*V**3-30*V**2+48*V-16)*W**8-2*V**6*(27*V**6-34*V**5+10*V**&
&4+80*V**3-120*V**2+96*V-32)*W**7+2*V**4*(21*V**8-13*V**7-33*V**6+126*V&
&**5-136*V**4+54*V**3+38*V**2-48*V+12)*W**6-V**4*(19*V**8+16*V**7-66*V*&
&*6+80*V**5+82*V**4-348*V**3+452*V**2-288*V+72)*W**5+V**2*(V**2+4*V-4)*&
&(4*V**8+9*V**7-37*V**6+70*V**5-66*V**4+26*V**3+10*V**2-16*V+4)*W**4-2*&
&(V-1)*V**2*(5*V**8+15*V**7-49*V**6+68*V**5-18*V**4-64*V**3+96*V**2-64*&
&V+16)*W**3+2*(V-1)**2*(8*V**8-4*V**7-2*V**6+21*V**5-31*V**4+19*V**3+3*&
&V**2-8*V+2)*W**2-(V-1)**3*(10*V**6-5*V**5+7*V**4-10*V**2+12*V-4)*W+(V-&
&1)**4*(4*V**4-5*V**3+7*V**2-4*V+2))-2*(V**8*(3*V**4-3*V**3+5*V**2-4*V+&
&2)*W**10-2*V**8*(8*V**4-7*V**3+12*V**2-10*V+5)*W**9+V**6*(38*V**6-24*V&
&**5+33*V**4-10*V**3-15*V**2+24*V-8)*W**8-V**6*(50*V**6-3*V**5-19*V**4+&
&76*V**3-118*V**2+96*V-32)*W**7+V**4*(38*V**8+43*V**7-75*V**6+98*V**5-1&
&22*V**4+54*V**3+38*V**2-48*V+12)*W**6-V**4*(16*V**8+60*V**7-43*V**6-44&
&*V**5+83*V**4-174*V**3+226*V**2-144*V+36)*W**5+V**2*(3*V**10+36*V**9+2&
&4*V**8-162*V**7+216*V**6-238*V**5+182*V**4-40*V**3-50*V**2+40*V-8)*W**&
&4-(V-1)*V**2*(9*V**8+47*V**7-58*V**6+20*V**5+11*V**4-70*V**3+98*V**2-6&
&4*V+16)*W**3+(V-1)**2*(14*V**8+16*V**7-21*V**6+19*V**5-30*V**4+19*V**3&
&+3*V**2-8*V+2)*W**2-(V-1)**3*(9*V**6+3*V**5-3*V**4+2*V**3-6*V**2+6*V-2&
&)*W+(V-1)**4*(3*V**4-2*V**3+3*V**2-2*V+1)))+8*(CQ-1)*GTR*N**4*(V-1)**2&
&*VC*W*(V*W-1)*(V*W-V+1)*(2*V**8*(V**4-2*V**3+4*V**2-4*V+2)*W**10-V**8*&
&(9*V**4-20*V**3+40*V**2-40*V+20)*W**9+V**6*(20*V**6-47*V**5+73*V**4-36&
&*V**3-22*V**2+48*V-16)*W**8-2*V**6*(13*V**6-28*V**5+20*V**4+48*V**3-10&
&4*V**2+96*V-32)*W**7+2*V**4*(10*V**8-17*V**7-5*V**6+84*V**5-130*V**4+7&
&2*V**3+32*V**2-48*V+12)*W**6-V**4*(9*V**8-8*V**7-10*V**6+52*V**5+6*V**&
&4-240*V**3+416*V**2-288*V+72)*W**5+V**2*(2*V**10+5*V**9+3*V**8-78*V**7&
&+262*V**6-442*V**5+390*V**4-112*V**3-92*V**2+80*V-16)*W**4-2*(V-1)*V**&
&2*(2*V**8+7*V**7-35*V**6+64*V**5-36*V**4-40*V**3+88*V**2-64*V+16)*W**3&
&+2*(V-1)**2*(4*V**8-7*V**7+6*V**6+12*V**5-29*V**4+22*V**3+2*V**2-8*V+2&
&)*W**2-(V-1)**3*(4*V**6-3*V**5+7*V**4-4*V**3-8*V**2+12*V-4)*W+(V-1)**4&
&*(2*V**4-3*V**3+5*V**2-4*V+2))+8*GTR*(V-1)**2*V**2*VC*W*(V**3*W**3-3*V&
&**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3&
&)*(2*CQ*(V**4*(V**2-V+1)*W**6-3*V**4*(V**2-V+1)*W**5+2*V**2*(2*V**4-2*&
&V**3+V**2+2*V-1)*W**4-V**2*(3*V**4-3*V**3-V**2+8*V-4)*W**3+(V**6-V**4+&
&V**3+2*V**2-3*V+1)*W**2-(V-1)*(V**4+V**3-2*V**2+2*V-1)*W+(V-1)**2*(V**&
&2-V+1))-V**4*(3*V**2-2*V+2)*W**6+2*V**4*(5*V**2-4*V+4)*W**5-2*V**2*(V*&
&*2-V+1)*(7*V**2+2*V-2)*W**4+2*V**2*(5*V**4-2*V**3-2*V**2+8*V-4)*W**3-(&
&V**2+V-1)**2*(3*V**2-2*V+2)*W**2+2*(V-1)*V**2*(V+1)*(2*V-1)*W-(V-1)**2&
&*(3*V**2-2*V+2)))/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LV1 = LOG(1-V)*(-16*N**5*(V-1)*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**&
&4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)&
&*((V*W-1)*(4*V**8*W**7-10*V**8*W**6+2*V**5*(8*V**3-5*V**2+6*V+1)*W**5-&
&V**5*(14*V**3-15*V**2+18*V+5)*W**4+V**2*(6*V**6-6*V**5+6*V**4+15*V**3-&
&8*V**2+5*V-2)*W**3-V**3*(2*V**5-3*V**4+9*V**3-2*V**2+5*V-1)*W**2+(2*V-&
&1)*(V**2-V+1)*(V**4+V**3+V**2+3*V-2)*W-2*(V-1)*(V+1)*(V**2-V+1)**2)+2*&
&CQ*((V-1)*W-1)*(V*W-V+1)*(V*W+V-1)*((V-1)*V*W+1)*(V**2*W-1)*(V**2*W-V+&
&1))-4*GTR*N**4*(V-1)**2*V*VC*W*(V*W-1)*(CQ*V*(V**2*W**2-2*V**2*W+2*V**&
&2-2*V+1)*(V**4*W**4+2*V**2*W**2+1)-2*(V**2*W-(V-1)**2)*(V**4*W**4-4*V*&
&*3*W**3+6*V**2*W**2-4*V*W+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2&
&*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+8*GTR*N**2*(&
&V-1)**2*V*VC*W*(V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W*&
&*2-4*(V-1)**3*V*W+(V-1)**4)*((V**3*W**3-3*V**2*W**2+3*V*W-1)*(2*V**5*W&
&**4-V**4*(3*V+2)*W**3+V**2*(5*V**3-3*V**2+5*V-2)*W**2-2*V**5*W+(V-1)*(&
&V+1)*(2*V**2-3*V+2))+CQ*V*(V*W-V+1)*(V**2*W**2+1)*(V**2*W**2-V*W+1)*(V&
&**2*W**2-2*V**2*W+2*V**2-2*V+1))-4*GTR*(V-1)**2*V**2*VC*W*(V**3*W**3-3&
&*V**2*W**2+3*V*W-1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-&
&4*(V-1)**3*V*W+(V-1)**4)*(2*(V*W-1)*(2*V**3*W**4-4*V**3*W**3+V*(V**3+2&
&*V**2-V+2)*W**2-V*(V**3-V+2)*W+(V-1)**3)+CQ*(V*W-V+1)*(V**2*W**2+1)*(V&
&**2*W**2-2*V**2*W+2*V**2-2*V+1)))/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*&
&(V*W-V+1)**5)
LV = LOG(V)*(-16*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*(4*V**14*W**13-22*V**&
&14*W**12+2*V**8*(31*V**6-5*V**5+8*V**4-7*V**3+6*V**2-3*V+1)*W**11-V**8&
&*(116*V**6-53*V**5+83*V**4-70*V**3+60*V**2-30*V+10)*W**10+2*V**6*(74*V&
&**8-41*V**7+40*V**6+7*V**5-20*V**4+31*V**3-29*V**2+16*V-4)*W**9-V**6*(&
&132*V**8-34*V**7-137*V**6+432*V**5-473*V**4+398*V**3-282*V**2+128*V-32&
&)*W**8+V**4*(92*V**10+8*V**9-261*V**8+514*V**7-265*V**6-36*V**5+124*V*&
&*4-168*V**3+132*V**2-60*V+12)*W**7-2*CQ*(V**2*W-V+1)*(V**2*W**2-2*V*W+&
&1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**6*(V**2-2*V+2)*W**6-V**4*(V**4&
&-2*V**3+4*V**2-4*V+2)*W**5-2*V**4*(V**4-4*V**3+5*V**2-2*V+1)*W**4+2*V*&
&*2*(V**6-4*V**5+7*V**4-8*V**3+9*V**2-6*V+2)*W**3-2*(V-1)*V**2*(V**4-4*&
&V**3+5*V**2-2*V+1)*W**2-(V-1)**2*(V**4-2*V**3+4*V**2-4*V+2)*W+(V-1)**3&
&*(V**2-2*V+2))-V**4*(54*V**10+19*V**9-228*V**8+241*V**7+430*V**6-999*V&
&**5+1005*V**4-792*V**3+468*V**2-180*V+36)*W**6+V**2*(22*V**12+56*V**11&
&-293*V**10+388*V**9+23*V**8-242*V**7-30*V**6+280*V**5-358*V**4+280*V**&
&3-144*V**2+48*V-8)*W**5-V**2*(4*V**12+56*V**11-198*V**10+171*V**9+214*&
&V**8-359*V**7-41*V**6+508*V**5-703*V**4+560*V**3-288*V**2+96*V-16)*W**&
&4+(V-1)*(16*V**12+18*V**11-242*V**10+597*V**9-719*V**8+661*V**7-549*V*&
&*6+424*V**5-298*V**4+150*V**3-52*V**2+12*V-2)*W**3-2*(V-1)**2*(12*V**1&
&0-31*V**9+11*V**8+76*V**7-133*V**6+127*V**5-61*V**4+10*V**3+5*V**2-5*V&
&+1)*W**2+2*(V-1)**3*(V**2-V+1)*(8*V**6-24*V**5+31*V**4-17*V**3+16*V**2&
&-9*V+3)*W-2*(V-1)**4*(V**2-V+1)**2*(2*V**2-5*V+5))+4*GTR*N**4*(V-1)**2&
&*VC*W*(V*W-1)*(V*W-V+1)*(CQ*V**2*(2*V**10*W**10-8*V**10*W**9+V**8*(19*&
&V**2-18*V+18)*W**8-4*V**8*(7*V**2-12*V+12)*W**7+4*V**6*(6*V**4-14*V**3&
&+25*V**2-22*V+11)*W**6-4*V**6*(3*V**4-9*V**3+29*V**2-40*V+20)*W**5+V**&
&4*(3*V**6-6*V**5+56*V**4-144*V**3+182*V**2-132*V+44)*W**4-4*(V-1)*V**4&
&*(V**4+5*V**3-17*V**2+24*V-12)*W**3+2*(V-1)**2*V**2*(5*V**4-8*V**3+17*&
&V**2-18*V+9)*W**2-4*(V-1)**3*V**2*(V**2+2*V-2)*W+(V-1)**4*(3*V**2-2*V+&
&2))-2*(2*V**8*(V**4-2*V**3+4*V**2-4*V+2)*W**10-V**8*(9*V**4-20*V**3+40&
&*V**2-40*V+20)*W**9+V**6*(20*V**6-47*V**5+73*V**4-36*V**3-22*V**2+48*V&
&-16)*W**8-2*V**6*(13*V**6-28*V**5+20*V**4+48*V**3-104*V**2+96*V-32)*W*&
&*7+2*V**4*(10*V**8-17*V**7-5*V**6+84*V**5-130*V**4+72*V**3+32*V**2-48*&
&V+12)*W**6-V**4*(9*V**8-8*V**7-10*V**6+52*V**5+6*V**4-240*V**3+416*V**&
&2-288*V+72)*W**5+V**2*(2*V**10+5*V**9+3*V**8-78*V**7+262*V**6-442*V**5&
&+390*V**4-112*V**3-92*V**2+80*V-16)*W**4-2*(V-1)*V**2*(2*V**8+7*V**7-3&
&5*V**6+64*V**5-36*V**4-40*V**3+88*V**2-64*V+16)*W**3+2*(V-1)**2*(4*V**&
&8-7*V**7+6*V**6+12*V**5-29*V**4+22*V**3+2*V**2-8*V+2)*W**2-(V-1)**3*(4&
&*V**6-3*V**5+7*V**4-4*V**3-8*V**2+12*V-4)*W+(V-1)**4*(2*V**4-3*V**3+5*&
&V**2-4*V+2)))-8*GTR*N**2*(V-1)**2*VC*W*(V*W-1)*(V*W-V+1)*(CQ*V**2*(2*V&
&**10*W**10-9*V**10*W**9+V**8*(21*V**2-10*V+10)*W**8-2*V**8*(15*V**2-14&
&*V+14)*W**7+2*V**6*(13*V**4-15*V**3+25*V**2-20*V+10)*W**6-V**6*(13*V**&
&4-10*V**3+44*V**2-68*V+34)*W**5+V**4*(3*V**6+8*V**5+16*V**4-68*V**3+84&
&*V**2-60*V+20)*W**4-2*(V-1)*V**4*(3*V**4+9*V**3-19*V**2+20*V-10)*W**3+&
&2*(V-1)**2*V**2*(5*V**4-V**3+6*V**2-10*V+5)*W**2-(V-1)**3*V**2*(6*V**2&
&+5*V-5)*W+(V-1)**4*(3*V**2-2*V+2))-2*V**8*(2*V**4-3*V**3+5*V**2-4*V+2)&
&*W**10+V**8*(23*V**4-30*V**3+50*V**2-40*V+20)*W**9-V**6*(58*V**6-61*V*&
&*5+79*V**4-20*V**3-30*V**2+48*V-16)*W**8+8*V**6*(10*V**6-5*V**5-V**4+2&
&0*V**3-30*V**2+24*V-8)*W**7-4*V**4*(16*V**8+9*V**7-29*V**6+57*V**5-65*&
&V**4+27*V**3+19*V**2-24*V+6)*W**6+V**4*(29*V**8+80*V**7-82*V**6-16*V**&
&5+130*V**4-348*V**3+452*V**2-288*V+72)*W**5-V**2*(6*V**10+59*V**9+15*V&
&**8-234*V**7+392*V**6-482*V**5+366*V**4-80*V**3-100*V**2+80*V-16)*W**4&
&+2*(V-1)*V**2*(9*V**8+36*V**7-56*V**6+40*V**5-4*V**4-64*V**3+96*V**2-6&
&4*V+16)*W**3-2*(V-1)**2*(14*V**8+7*V**7-11*V**6+17*V**5-29*V**4+19*V**&
&3+3*V**2-8*V+2)*W**2+(V-1)**3*(18*V**6-V**5+3*V**4-10*V**2+12*V-4)*W-(&
&V-1)**4*(6*V**4-5*V**3+7*V**2-4*V+2))+4*GTR*(V-1)**2*V**2*VC*W*(V**3*W&
&**3-3*V**2*W**2+3*V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(&
&V-1)**3)*(CQ*(2*V**6*W**6-6*V**6*W**5+3*V**4*(3*V**2-2*V+2)*W**4-4*V**&
&4*(2*V**2-3*V+3)*W**3+V**2*(3*V**4-4*V**3+10*V**2-12*V+6)*W**2-2*(V-1)&
&*V**2*(V**2+3*V-3)*W+(V-1)**2*(3*V**2-2*V+2))-4*(V**4*(V**2-V+1)*W**6-&
&3*V**4*(V**2-V+1)*W**5+2*V**2*(V+1)*(2*V-1)*(V**2-V+1)*W**4-V**2*(V+2)&
&*(3*V-2)*(V**2-V+1)*W**3+(V**2+V-1)*(V**4+2*V**3-3*V**2+2*V-1)*W**2-(V&
&-1)*(2*V**4+3*V**3-4*V**2+2*V-1)*W+(V-1)**2*(2*V**2-V+1))))/(N**2*(V-1&
&)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LVW = (16*N**5*(V-1)*VC*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W*&
&*2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1&
&)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(4*V**6*W**5-4*V**5*(V+1)*W**4&
&+V**4*(8*V**2-11*V+19)*W**3-2*V**3*(V**2+V+4)*W**2-(V**2-V+1)*(2*V**3-&
&5*V**2-9*V+4)*W-4*(V**2-V+1)**2)-16*GTR*N**2*(V-1)**2*V**2*VC*W*(V**2*&
&W**2-2*V**2*W+2*V**2-2*V+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**&
&2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*&
&(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+8*GTR*(V-1)**2*V**2*VC*W*(&
&V**2*W**2-2*(V-1)*V*W+2*V**2-2*V+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**&
&3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*&
&W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-8*GTR*N**4*(V-1)**&
&2*V*VC*W*(V**2*W-(V-1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2&
&*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(&
&V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))*LOG(1-V*W)/(N**2*(V-1)**3*&
&V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LTVW = (16*N**5*(V-1)*VC*(V*W-V+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3&
&-10*V**2*W**2+5*V*W-1)*((V-1)*V*W*(V**7*W**6+2*(V-2)*V**6*W**5-V**4*(9&
&*V**3-16*V**2+6*V-1)*W**4+4*(V-1)*V**3*(3*V**3-3*V**2+2*V-1)*W**3-(V-1&
&)**2*V**2*(9*V**3-6*V**2+11*V-6)*W**2+2*(V-1)**3*V*(V**3-V**2+4*V-2)*W&
&+(V-1)**4*(V+1)*(V**2-V+1))+4*CQ*(W-V+1)*(V*W+1)*(V*W+(V-1)**2)*(V**2*&
&W-(V-1)**2)*(V**2*W-V+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2))+8*GTR*N**4*&
&(V-1)**2*V*VC*W*(V*W-V+1)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*&
&W**2+5*V*W-1)*(-V**7*W**6+V**6*(V+1)*W**5+CQ*V*(V**2*W**2-2*V**2*W+V**&
&2+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)+V**4*(V**3-4*V**2+2*V-1&
&)*W**4-2*(V-2)*(V-1)**3*V**3*W**3+(V-1)**2*V**2*(V**3-6*V**2+7*V-6)*W*&
&*2+(V-1)**4*V*(V**2+V-4)*W-(V-1)**4*(V+1)*(V**2-V+1))+8*GTR*(V-1)**2*V&
&**2*VC*W*(CQ*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*V**2*W+V**2+1)+4*(V-1)*&
&V**2*(W-1)*W)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-(V-1)**3)*(V&
&**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)-16*GTR*N**2*(V&
&-1)**2*V**2*VC*W*(V*W-V+1)*(V**2*W**2-2*V**2*W+V**2+1)*(CQ*(V**2*W**2+&
&(V-1)**2)*(V**2*W**2-(V-1)*V*W+(V-1)**2)-(V-1)*V*W*(3*V**2*W**2-4*(V-1&
&)*V*W+3*(V-1)**2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*&
&V*W-1))*LOG(V*W-V+1)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)**5)
LW = (16*N**5*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V*&
&*4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4&
&)*(2*V**6*(V**2-V+1)*W**6-V**6*(8*V**2-13*V+13)*W**5+V**4*(14*V**4-29*&
&V**3+30*V**2-2*V+1)*W**4-2*V**2*(7*V**6-15*V**5+12*V**4+8*V**3-9*V**2+&
&6*V-2)*W**3+V**2*(10*V**6-24*V**5+26*V**4-3*V**3-V**2+3*V-1)*W**2-4*(V&
&**2-V+1)*(V**6-V**5-V**4+3*V**3+V**2-3*V+1)*W+2*(V-1)*(V**2-V+1)**2*(2&
&*V**2-3*V+3))-4*GTR*N**4*(V-1)**2*VC*W*(CQ*(V**2-2*V+2)**2*(2*W**2-2*W&
&+1)-2*V**2*(V**2*W-V+1))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W&
&**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-&
&1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+8*GTR*N**2*(V-1)**2*VC*W*(V**&
&4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6&
&*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**2*(2*V**4*W**4-V**2*(&
&5*V**2-2*V+2)*W**3+3*V**2*(V**2+V-1)*W**2-2*(V**4+V**2-2*V+1)*W+(V-1)*&
&(2*V**2-V+1))+CQ*(V**2-2*V+2)*(V**2-V+1)*(V*W-1)*(V*W-V+1)*(2*W**2-2*W&
&+1))-4*GTR*(V-1)**2*V**2*VC*W*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W&
&+1)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(&
&V-1)**4)*(2*(V**4*W**4-V**2*(V**2+2*V-2)*W**3+5*(V-1)*V**2*W**2-2*(V-1&
&)*(2*V**2-V+1)*W+2*(V-1)*(V**2-V+1))+CQ*(V**2-2*V+2)*(V*W-1)*(V*W-V+1)&
&*(2*W**2-2*W+1)))*LOG(W)/(N**2*(V-1)**3*V**3*W**2*(V*W-1)**5*(V*W-V+1)&
&**5)
CVC = (-8*N**5*(V-1)*VC*(V*W-1)*(V*W-V+1)*((V-1)*V**8*(17*V**4-50*V**3&
&+90*V**2-80*V+40)*W**10-2*(V-1)*V**8*(34*V**4-114*V**3+213*V**2-198*V+&
&99)*W**9+V**6*(11*V**8+103*V**7-571*V**6+1145*V**5-935*V**4-13*V**3+75&
&1*V**2-640*V+160)*W**8-V**6*(33*V**8+71*V**7-602*V**6+806*V**5+869*V**&
&4-3296*V**3+4048*V**2-2528*V+632)*W**7+V**4*(33*V**10+81*V**9-642*V**8&
&+796*V**7+1121*V**6-3554*V**5+3350*V**4-416*V**3-1696*V**2+1200*V-240)&
&*W**6-V**4*(11*V**10+123*V**9-651*V**8+1137*V**7-719*V**6+743*V**5-330&
&9*V**4+6872*V**3-7028*V**2+3540*V-708)*W**5+(V-1)*V**2*(61*V**10-124*V&
&**9+13*V**8-168*V**7+1657*V**6-3402*V**5+3178*V**4-792*V**3-1002*V**2+&
&800*V-160)*W**4-(V-1)**2*V**2*(84*V**8-166*V**7-183*V**6+818*V**5-397*&
&V**4-888*V**3+1752*V**2-1248*V+312)*W**3+(V-1)**3*(86*V**8-285*V**7+38&
&2*V**6-24*V**5-373*V**4+350*V**3+70*V**2-160*V+40)*W**2-(V-1)**4*(29*V&
&**6-47*V**5+87*V**4-42*V**3-74*V**2+114*V-38)*W+17*(V-1)**5*(V**2-V+1)&
&**2)+8*GTR*N**4*(V-1)*VC*(V*W-1)*(V*W-V+1)*((V-1)*V**8*(7*V**4-34*V**3&
&+78*V**2-88*V+44)*W**10-(V-1)*V**8*(31*V**4-156*V**3+372*V**2-432*V+21&
&6)*W**9+V**6*(4*V**8+53*V**7-368*V**6+922*V**5-1015*V**4+172*V**3+764*&
&V**2-704*V+176)*W**8-2*V**6*(6*V**8+23*V**7-194*V**6+314*V**5+257*V**4&
&-1460*V**3+2092*V**2-1376*V+344)*W**7+V**4*(12*V**10+45*V**9-381*V**8+&
&602*V**7+694*V**6-3226*V**5+3670*V**4-640*V**3-1820*V**2+1320*V-264)*W&
&**6-V**4*(4*V**10+51*V**9-351*V**8+876*V**7-952*V**6+892*V**5-3060*V**&
&4+6976*V**3-7504*V**2+3840*V-768)*W**5+(V-1)*V**2*(23*V**10-95*V**9+28&
&7*V**8-834*V**7+2174*V**6-3702*V**5+3362*V**4-768*V**3-1128*V**2+880*V&
&-176)*W**4-2*(V-1)**2*V**2*(9*V**8+20*V**7-213*V**6+422*V**5-133*V**4-&
&564*V**3+972*V**2-672*V+168)*W**3+(V-1)**3*(34*V**8-138*V**7+173*V**6+&
&72*V**5-347*V**4+250*V**3+122*V**2-176*V+44)*W**2+(V-1)**4*(2*V**6-29*&
&V**5+33*V**4-48*V**3+124*V**2-120*V+40)*W+(V-1)**5*(V**2-V+1)*(7*V**2-&
&4*V+4))+24*GTR*N**2*(V-1)**2*VC*(V*W-1)*(V*W-V+1)*(2*(V-1)*V**8*(5*V**&
&2-6*V+6)*W**10+V**8*(V**4-46*V**3+106*V**2-120*V+60)*W**9-V**6*(4*V**6&
&-95*V**5+175*V**4-112*V**3-64*V**2+144*V-48)*W**8+6*V**6*(V**6-18*V**5&
&+14*V**4+40*V**3-100*V**2+96*V-32)*W**7-2*V**4*(2*V**8-49*V**7+35*V**6&
&+142*V**5-320*V**4+198*V**3+102*V**2-144*V+36)*W**6+V**4*(V**8-80*V**7&
&+166*V**6-136*V**5+194*V**4-756*V**3+1260*V**2-864*V+216)*W**5+(V-1)*V&
&**2*(41*V**8-88*V**7+218*V**6-436*V**5+610*V**4-336*V**3-112*V**2+192*&
&V-48)*W**4-2*(V-1)*V**2*(5*V**8-11*V**7+49*V**6-80*V**5+2*V**4+180*V**&
&3-284*V**2+192*V-48)*W**3+2*(V-1)**2*(4*V**8+2*V**7-6*V**6-11*V**5+47*&
&V**4-33*V**3-17*V**2+24*V-6)*W**2-(V-1)**3*(10*V**6-29*V**5+35*V**4-24&
&*V**3+42*V**2-36*V+12)*W-(V-1)**5*V**2)-24*AL*GTR*N**2*(V-1)**2*VC*(V*&
&W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V*&
&*2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(V**4*(3*V**4-6*V**3+10*V**2-8*V+4)*W&
&**6-V**3*(4*V**5+3*V**4-14*V**3+32*V**2-28*V+16)*W**5+V**2*(3*V**6+5*V&
&**5-4*V**4+30*V**2-32*V+24)*W**4-2*V*(3*V**6+V**5-3*V**4+10*V**3-4*V+8&
&)*W**3+(10*V**6-12*V**5+11*V**4+10*V**3-10*V**2+8*V+4)*W**2-(6*V**5-10&
&*V**4+15*V**3-6*V**2+4)*W+3*V**4-5*V**3+6*V**2-4*V+2)+12*AL*GTR*N**4*(&
&V-1)**2*VC*(V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-&
&10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(V**4*(3*V**4-8*V**3+16&
&*V**2-16*V+8)*W**6-4*V**3*(V**5-4*V**3+12*V**2-14*V+8)*W**5+V**2*(3*V*&
&*6+2*V**5-9*V**4+8*V**3+36*V**2-64*V+48)*W**4-4*V*(V**6-2*V**4+8*V**3-&
&4*V**2-4*V+8)*W**3+(10*V**6-20*V**5+21*V**4+8*V**3-24*V**2+16*V+8)*W**&
&2-4*(V**5-3*V**4+6*V**3-4*V**2+2)*W+3*V**4-6*V**3+9*V**2-8*V+4)-24*GTR&
&*(V-1)**2*V**2*VC*(V**2*W**2+V**2-2*V+2)*(V**5*W**5-5*V**4*W**4+10*V**&
&3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*&
&V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)+12*AL*GTR*(V-&
&1)**2*V**2*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**2*(3*V**2-4*V+4)*W**&
&4-4*V*(V**3-V+2)*W**3+(3*V**4-2*V**2+4*V+4)*W**2-2*(V**3+2)*W+3*V**2-4&
&*V+3)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V&
&**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-24*CQ*GTR*(V-1)**2*V**2*VC*(V**2*W**&
&2+(V-1)**2)*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)*&
&*3*V*W+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*&
&W-1)-24*CQ*GTR*N**4*(V-1)**2*V**2*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*&
&(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)*(V**5*W**5-5*V**4*W**4+10*V*&
&*3*W**3-10*V**2*W**2+5*V*W-1)+48*CQ*GTR*N**2*(V-1)**2*V**2*VC*(V**2*W*&
&*2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-(V-1)*V*W+(V-&
&1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1))/(N**&
&2*(V-1)**3*V**3*W*(V*W-1)**5*(V*W-V+1)**5)/3.D0
LM = LOG(S/M**2)*(-32*N**5*VC*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*&
&V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**10*W**9-V**9*(V+5)*W**8+V**4*(3&
&*V**6-2*V**5+21*V**4-7*V**3+6*V**2-3*V+1)*W**7-V**3*(V**7+9*V**6-9*V**&
&5+47*V**4-22*V**3+21*V**2-11*V+4)*W**6+V**2*(V**8+18*V**6-12*V**5+53*V&
&**4-18*V**3+24*V**2-14*V+6)*W**5-V*(V**9+V**8+20*V**6-3*V**5+27*V**4+1&
&0*V**3+6*V**2-6*V+4)*W**4+(4*V**9-6*V**8+13*V**7+2*V**6+15*V**5-2*V**4&
&+25*V**3-6*V**2+V+1)*W**3-(6*V**8-14*V**7+30*V**6-27*V**5+34*V**4-17*V&
&**3+18*V**2-3*V+1)*W**2+V*(V**2-V+1)*(4*V**4-7*V**3+13*V**2-7*V+8)*W-(&
&V**2-V+1)**2*(V**2-V+2))+8*GTR*N**2*(V-1)*VC*W*(V**4*W**4-4*(V-1)*V**3&
&*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**4*(3*V**4-6*V*&
&*3+10*V**2-8*V+4)*W**6-V**3*(4*V**5+3*V**4-14*V**3+32*V**2-28*V+16)*W*&
&*5+V**2*(3*V**6+5*V**5-4*V**4+30*V**2-32*V+24)*W**4-2*V*(3*V**6+V**5-3&
&*V**4+10*V**3-4*V+8)*W**3+(10*V**6-12*V**5+11*V**4+10*V**3-10*V**2+8*V&
&+4)*W**2-(6*V**5-10*V**4+15*V**3-6*V**2+4)*W+3*V**4-5*V**3+6*V**2-4*V+&
&2)-4*GTR*N**4*(V-1)*VC*W*(V**4*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*&
&W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**4*(3*V**4-8*V**3+16*V**2-16*V+8)*W**&
&6-4*V**3*(V**5-4*V**3+12*V**2-14*V+8)*W**5+V**2*(3*V**6+2*V**5-9*V**4+&
&8*V**3+36*V**2-64*V+48)*W**4-4*V*(V**6-2*V**4+8*V**3-4*V**2-4*V+8)*W**&
&3+(10*V**6-20*V**5+21*V**4+8*V**3-24*V**2+16*V+8)*W**2-4*(V**5-3*V**4+&
&6*V**3-4*V**2+2)*W+3*V**4-6*V**3+9*V**2-8*V+4)-4*GTR*(V-1)*V**2*VC*W*(&
&V**2*W**2-2*V*W+1)*(V**2*(3*V**2-4*V+4)*W**4-4*V*(V**3-V+2)*W**3+(3*V*&
&*4-2*V**2+4*V+4)*W**2-2*(V**3+2)*W+3*V**2-4*V+3)*(V**4*W**4-4*(V-1)*V*&
&*3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4))/(N**2*(V-1)**2*&
&V**3*W**2*(V*W-1)**4*(V*W-V+1)**4)
LMP = LOG(S/MP**2)*(-32*N**5*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V&
&*W+1)*(V**10*W**9-V**9*(6*V-5)*W**8+V**8*(18*V**2-31*V+15)*W**7-(V-1)*&
&V**7*(35*V**2-57*V+30)*W**6+(V-1)*V**4*(47*V**5-120*V**4+114*V**3-41*V&
&**2-2*V+1)*W**5-(V-1)**3*V**3*(45*V**4-66*V**3+39*V**2+3*V-4)*W**4+(V-&
&1)**3*V**2*(32*V**5-75*V**4+72*V**3-24*V**2-8*V+6)*W**3-(V-1)**4*V*(17&
&*V**5-38*V**4+40*V**3-16*V**2-2*V+4)*W**2+(V-1)**5*(V**2-V+1)*(6*V**3-&
&7*V**2+3*V+1)*W-(V-1)**6*(V**2-V+1)**2)-4*GTR*N**4*(V-1)*V**2*VC*W*(V*&
&*2*W**2-2*V**2*W+V**2+1)*(V**4*W**4+2*(V-1)**2*V**2*W**2+(V-1)**4)*(V*&
&*4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)+8*GTR*N**2*(V-1)*V**2*VC*W*(V&
&**2*W**2+(V-1)**2)*(V**2*W**2-(V-1)*V*W+(V-1)**2)*(V**2*W**2-2*V**2*W+&
&V**2+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)-4*GTR*(V-1)*V**2*V&
&C*W*(V**2*W**2+(V-1)**2)*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**2*W**2-2&
&*V**2*W+V**2+1)*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1))/(N**2*(V-&
&1)**2*V**3*W**2*(V*W-1)**4*(V*W-V+1)**4)
STRUV15=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUV16(W,V,X3,S)
double precision::CVC 
double precision::LM 
double precision::LMP 
double precision::LTVW 
double precision::LV 
double precision::LV1 
double precision::LVW 
double precision::LW 
double precision::LW1 
double precision::M
double precision::MP
double precision::S
double precision::STRUV16
double precision::V
double precision::W
double precision::X3
M=DSQRT(Q2FAC)
MP=DSQRT(Q2FRAG)
LW1 = LOG(1-W)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*&
&V*W+(V-1)**2)*(CQ*(16*V**11*W**11-16*V**10*(5*V-1)*W**10+V**6*(224*V**&
&5-141*V**4+23*V**3+4*V**2+2)*W**9-V**6*(384*V**5-323*V**4-28*V**3+79*V&
&**2+8)*W**8+V**4*(416*V**7-280*V**6-431*V**5+470*V**4-88*V**3+5*V**2+6&
&*V-6)*W**7-V**4*(320*V**7-120*V**6-770*V**5+829*V**4-170*V**3-37*V**2+&
&18*V-18)*W**6+V**2*(192*V**9-35*V**8-679*V**7+670*V**6+129*V**5-358*V*&
&*4+122*V**3-15*V**2-12*V+6)*W**5-V**2*(80*V**9+83*V**8-710*V**7+903*V*&
&*6-244*V**5-203*V**4+98*V**3+9*V**2-24*V+12)*W**4+(V-1)*(16*V**10+144*&
&V**9-413*V**8+265*V**7+223*V**6-298*V**5+67*V**4+29*V**3-15*V**2-4*V+2&
&)*W**3-(V-1)**2*(48*V**8-293*V**6+522*V**5-320*V**4+60*V**3+7*V**2+2*V&
&-2)*W**2+(V-1)**3*(12*V**2-13*V+5)*(4*V**4-5*V**3-V**2+5*V+1)*W-4*(V-1&
&)**4*(2*V**2-4*V+3)*(2*V**2-2*V+1))-32*V**11*W**11+16*V**9*(8*V**2+V-1&
&)*W**10-V**6*(288*V**5-51*V**4+13*V**3-44*V**2+20*V-2)*W**9+V**6*(448*&
&V**5-335*V**4+262*V**3-243*V**2+116*V-8)*W**8-V**4*(448*V**7-398*V**6-&
&13*V**5+294*V**4-276*V**3+189*V**2-66*V+6)*W**7+V**4*(320*V**7-206*V**&
&6-732*V**5+1539*V**4-1246*V**3+621*V**2-210*V+18)*W**6-V**2*(224*V**9-&
&239*V**8-489*V**7+952*V**6-29*V**5-768*V**4+590*V**3-267*V**2+72*V-6)*&
&W**5+V**2*(128*V**9-83*V**8-666*V**7+1167*V**6-82*V**5-1091*V**4+948*V&
&**3-425*V**2+124*V-12)*W**4-(V-1)*(32*V**10+208*V**9-815*V**8+903*V**7&
&-147*V**6-96*V**5-197*V**4+229*V**3-107*V**2+24*V-2)*W**3+(V-1)**2*(96&
&*V**8-64*V**7-439*V**6+938*V**5-648*V**4+172*V**3+V**2-10*V+2)*W**2-(V&
&-1)**3*(96*V**6-240*V**5+199*V**4+26*V**3-76*V**2+18*V+9)*W+8*(V-1)**4&
&*(2*V**2-4*V+3)*(2*V**2-2*V+1))+2*N**2*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V&
&**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2&
&+5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V**3*(3*V**4-3*V**3+4*V**2-2*V+2)*W**6-&
&V**2*(4*V**5+5*V**4-5*V**3+10*V**2-4*V+6)*W**5+V*(3*V**6+7*V**5+4*V**4&
&-3*V**3+9*V**2+6)*W**4-(9*V**6+V**5+7*V**4-8*V**3+7*V**2+4*V+2)*W**3+(&
&13*V**5-13*V**4+17*V**3-16*V**2+7*V+2)*W**2-(9*V**4-13*V**3+12*V**2-9*&
&V+3)*W+2*(V-1)*(2*V**2-2*V+1))-(V*W-1)*(4*V**3*(V**3+V**2-V+1)*W**5-2*&
&V**2*(3*V**4+2*V**3+5*V**2-2*V+4)*W**4+V*(4*V**5+8*V**4-2*V**3+5*V**2+&
&9*V+4)*W**3-V*(16*V**4-22*V**3+23*V**2-15*V+14)*W**2+(16*V**4-28*V**3+&
&25*V**2-14*V+5)*W-4*(V-1)*(2*V**2-2*V+1)))-(V-1)**2*VC*W*(V**2*W**2-2*&
&V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*&
&V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V**3*(3*V**3-2*V**2+2*V-2)*W**&
&5-2*V**2*(2*V**4+3*V**3-2*V**2+2*V-3)*W**4+V*(3*V**5+7*V**4+4*V**3-V**&
&2-6)*W**3-(7*V**5+3*V**4+V**3+V**2-4*V-2)*W**2+(7*V**4-5*V**3+3*V**2-V&
&-2)*W-(V-1)*(V**2+1))-V**3*(5*V**3-2*V**2+2*V-2)*W**5+2*V**2*(4*V**4+6&
&*V**3-2*V**2+2*V-3)*W**4-V*(5*V**5+19*V**4+12*V**3-V**2-6)*W**3+(13*V*&
&*5+15*V**4+9*V**3+V**2-4*V-2)*W**2-(13*V**4-V**3+9*V**2-V-2)*W+(3*V-1)&
&*(V**2+1)))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
LV1 = LOG(1-V)*(N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V&
&**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(2*(V**2*W**2-2&
&*V*W+1)*(16*V**7*W**7-40*V**7*W**6+2*V**4*(32*V**3-19*V**2+7*V+4)*W**5&
&-V**4*(56*V**3-57*V**2+23*V+14)*W**4+V**2*(24*V**5-20*V**4-9*V**3+31*V&
&**2-4*V-6)*W**3-V**2*(8*V**5-11*V**4+17*V**3-21*V**2+23*V-12)*W**2+(V-&
&1)*(8*V**5-3*V**4+3*V**2-2*V+2)*W-4*(V-1)**2*V*(2*V**2-2*V+1))+CQ*(V-1&
&)*(V*W-V+1)*(V**6*(16*V-15)*W**6+2*V**5*(8*V**2-41*V+30)*W**5-V**3*(46&
&*V**3-150*V**2+87*V+4)*W**4+4*V**2*(10*V**3-30*V**2+12*V+3)*W**3+4*V*(&
&8*V**2+V-3)*W**2-4*(4*V**3-2*V**2+3*V-1)*W+4*(2*V**2-2*V+1)))-2*N**2*(&
&V-1)*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V&
&**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V-1)*(V**&
&6*W**6-V**5*(2*V+3)*W**5+2*V**3*(V**3+2*V**2+3*V-1)*W**4-V**2*(6*V**3+&
&4*V**2+9*V-6)*W**3+V*(10*V**3-2*V**2+11*V-6)*W**2-2*(4*V**3-2*V**2+3*V&
&-1)*W+2*(2*V**2-2*V+1))+(2*V**2*(5*V-1)*W**3-V**2*(7*V+1)*W**2+V*(V**3&
&+4*V**2-3*V+2)*W-4*(V-1)*(2*V**2-2*V+1))*(V**3*W**3-3*V**2*W**2+3*V*W-&
&1))+(V-1)**2*VC*W*(V**2*W**2-2*V*W+1)*(2*(V-1)*(V**2-2*V+2)*(V**3*W**3&
&-3*V**2*W**2+3*V*W-1)+CQ*V**2*W*(V**2*W**2-2*V*W+2)*(V**2*W**2-2*V**2*&
&W+2*V**2-2*V+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10&
&*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**2*(V-1)**3*V**2*W**2&
&*(V*W-1)**5*(V*W-V+1)**5)
LV = LOG(V)*(-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W&
&+(V-1)**2)*(CQ*(16*V**11*W**11-16*V**10*(5*V-1)*W**10+V**9*(224*V**2-1&
&43*V+23)*W**9-V**8*(384*V**3-331*V**2-28*V+63)*W**8+V**6*(416*V**5-293&
&*V**4-437*V**3+450*V**2-100*V+4)*W**7-V**6*(320*V**5-131*V**4-788*V**3&
&+825*V**2-206*V-12)*W**6+(V-1)*V**4*(192*V**6+152*V**5-548*V**4+127*V*&
&*3+226*V**2-113*V+12)*W**5-(V-1)*V**4*(80*V**6+162*V**5-560*V**4+341*V&
&**3+97*V**2-119*V+15)*W**4+4*(V-1)*V**2*(4*V**8+36*V**7-104*V**6+64*V*&
&*5+56*V**4-78*V**3+19*V**2+7*V-3)*W**3-4*(V-1)**2*V**2*(12*V**6-74*V**&
&4+130*V**3-81*V**2+14*V+2)*W**2+4*(V-1)**3*(12*V**6-28*V**5+18*V**4+12&
&*V**3-15*V**2+3*V+1)*W-4*(V-1)**4*(2*V**2-4*V+3)*(2*V**2-2*V+1))-32*V*&
&*11*W**11+144*V**11*W**10-V**6*(352*V**5-109*V**4+71*V**3-44*V**2+20*V&
&-2)*W**9+V**6*(576*V**5-435*V**4+316*V**3-185*V**2+80*V-8)*W**8-V**4*(&
&624*V**7-552*V**6+5*V**5+446*V**4-380*V**3+189*V**2-66*V+6)*W**7+V**4*&
&(480*V**7-340*V**6-858*V**5+1885*V**4-1470*V**3+627*V**2-198*V+18)*W**&
&6-V**2*(304*V**9-219*V**8-927*V**7+1698*V**6-539*V**5-674*V**4+606*V**&
&3-267*V**2+72*V-6)*W**5+V**2*(144*V**9+7*V**8-1102*V**7+1815*V**6-368*&
&V**5-1315*V**4+1234*V**3-539*V**2+144*V-12)*W**4-(V-1)*(32*V**10+240*V&
&**9-913*V**8+873*V**7+243*V**6-586*V**5+47*V**4+181*V**3-107*V**2+24*V&
&-2)*W**3+(V-1)**2*(96*V**8-64*V**7-553*V**6+1186*V**5-808*V**4+176*V**&
&3+35*V**2-22*V+2)*W**2-(V-1)**3*(96*V**6-272*V**5+225*V**4+64*V**3-126&
&*V**2+36*V+9)*W+8*(V-1)**4*(2*V**2-5*V+4)*(2*V**2-2*V+1))+2*N**2*(V-1)&
&*VC*(V**2*W**2-2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*&
&W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(CQ*(V-1)*(V**6*W*&
&*6-V**5*(2*V+3)*W**5+2*V**3*(V**3+2*V**2+3*V-1)*W**4-V**2*(6*V**3+4*V*&
&*2+9*V-6)*W**3+V*(10*V**3-2*V**2+11*V-6)*W**2-2*(4*V**3-2*V**2+3*V-1)*&
&W+2*(2*V**2-2*V+1))-(V*W-1)*(V**3*(3*V**3+5*V**2-4*V+4)*W**5-2*V**2*(2&
&*V**4+3*V**3+5*V**2-2*V+4)*W**4+V*(3*V**2+3*V+1)*(V**3+2*V**2-3*V+4)*W&
&**3-V*(16*V**4-20*V**3+21*V**2-15*V+14)*W**2+(17*V**4-29*V**3+26*V**2-&
&15*V+5)*W-4*(V-1)*(2*V**2-2*V+1)))-(V-1)**2*VC*W*(V**2*W**2-2*V*W+1)*(&
&-V**3*(3*V**3-2*V**2+2*V-2)*W**5+2*V**2*(2*V**4+3*V**3-2*V**2+2*V-3)*W&
&**4-V*(3*V**5+7*V**4+4*V**3-V**2-6)*W**3+CQ*V**2*W*(V**2*W**2-2*V*W+2)&
&*(V**2*W**2-2*V**2*W+2*V**2-2*V+1)+(7*V**5+3*V**4+V**3+V**2-4*V-2)*W**&
&2-(7*V**4-5*V**3+3*V**2-V-2)*W+(V-1)*(V**2+1))*(V**5*W**5-5*(V-1)*V**4&
&*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)&
&**5))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
LVW = (-2*N**4*(V-1)*VC*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W*&
&*2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1&
&)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(16*V**5*W**5-16*V**4*(V+1)*W*&
&*4+V**3*(32*V**2-43*V+43)*W**3-V**2*(10*V**2-7*V+13)*W**2-(6*V**4-21*V&
&**3+12*V**2+3*V-4)*W-4*(V-1)*(2*V**2-2*V+1))+2*N**2*(V-1)*VC*(V**3*(V+&
&7)*W**3-V**2*(2*V**2+V+5)*W**2+(2*V**4+V**3-V+2)*W-4*(V-1)*(2*V**2-2*V&
&+1))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5*W&
&**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V&
&-1)**4*V*W-(V-1)**5)-2*(V-1)**2*VC*W*(V**3*W**2-2*V**3*W+(2*V-1)*(V**2&
&-V+2))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V**5&
&*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*&
&(V-1)**4*V*W-(V-1)**5))*LOG(1-V*W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5&
&*(V*W-V+1)**5)
LTVW = (-2*N**4*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**5*W**5-5&
&*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(8*CQ*(2*V**2*W**2-2*V*(&
&2*V-1)*W+2*V**2-2*V+1)*(V**6*W**6-3*(V-1)*V**5*W**5+6*(V-1)**2*V**4*W*&
&*4-7*(V-1)**3*V**3*W**3+6*(V-1)**4*V**2*W**2-3*(V-1)**5*V*W+(V-1)**6)+&
&(V-1)*W*(5*V**6*W**5+V**5*(5*V-7)*W**4-(V-1)*V**3*(16*V**3-8*V**2+8*V-&
&3)*W**3+(V-1)**2*V**2*(32*V**3-40*V**2+36*V-9)*W**2-(V-1)**3*V*(16*V**&
&3-27*V**2+27*V-9)*W-(V-1)**4*(5*V**2-3*V+3)))+2*N**2*(V-1)**2*VC*W*(V*&
&*3*W**2-(V-3)*V**2*W+(V-2)*(V-1)**2)*(V**5*W**5-5*V**4*W**4+10*V**3*W*&
&*3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3&
&*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)-2*(V-1)**2*VC*W*(&
&V**3*W**2+V**2*W+(V-1)*(V**2-2*V+3))*(V**5*W**5-5*V**4*W**4+10*V**3*W*&
&*3-10*V**2*W**2+5*V*W-1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3&
&*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5))*LOG(V*W-V+1)/(N*&
&*2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
LW = (-N**4*(V-1)*VC*(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**4&
&*W**4-4*(V-1)*V**3*W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*&
&(2*(8*V**5*(V**2-V+1)*W**6-V**5*(32*V**2-51*V+43)*W**5+2*V**3*(28*V**4&
&-56*V**3+47*V**2-4*V-1)*W**4-V**2*(56*V**5-117*V**4+67*V**3+45*V**2-38&
&*V+3)*W**3+V*(40*V**6-94*V**5+54*V**4+55*V**3-74*V**2+29*V-6)*W**2-(V-&
&1)*(16*V**6-16*V**5-30*V**4+73*V**3-45*V**2+7*V+3)*W+4*(V-1)**2*(2*V**&
&2-3*V+2)*(2*V**2-2*V+1))-CQ*(V**2+1)**2*W*(V*W-1)*(V*W-V+1)*(2*W**2-2*&
&W+1))+(V-1)**2*VC*W*(2*(V**3*W**2+V**2*W+V-1)+CQ*(V-1)*(V**2+1)*(2*W**&
&2-2*W+1))*(V**5*W**5-5*V**4*W**4+10*V**3*W**3-10*V**2*W**2+5*V*W-1)*(V&
&**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2&
&+5*(V-1)**4*V*W-(V-1)**5)-2*N**2*(V-1)*VC*(V**4*(V+7)*W**4-V**3*(V+3)*&
&*2*W**3+CQ*(V**2+1)*(V**2-V+1)*W*(V*W-1)*(2*W**2-2*W+1)+V*(8*V**3+V**2&
&-V+4)*W**2-(10*V**4-15*V**3+16*V**2-11*V+4)*W+4*(V-1)*(2*V**2-2*V+1))*&
&(V**4*W**4-4*V**3*W**3+6*V**2*W**2-4*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W*&
&*4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5&
&))*LOG(W)/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*(V*W-V+1)**5)
CVC = (-N**4*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V**2*W**2-2*(V-1)*V*W+(V-1)&
&**2)*(32*V**11*W**11-32*V**10*(5*V-1)*W**10+V**6*(416*V**5-245*V**4+59&
&*V**3-78*V**2+52*V-12)*W**9-V**6*(704*V**5-641*V**4+96*V**3-177*V**2+1&
&92*V-46)*W**8+V**4*(832*V**7-906*V**6-165*V**5+228*V**4-123*V**3+314*V&
&**2-192*V+36)*W**7-V**4*(704*V**7-714*V**6-964*V**5+1507*V**4-1102*V**&
&3+965*V**2-522*V+102)*W**6+V**2*(416*V**9-161*V**8-1817*V**7+2436*V**6&
&-1023*V**5+48*V**4+341*V**3-456*V**2+228*V-36)*W**5-V**2*(160*V**9+235&
&*V**8-1864*V**7+2185*V**6-180*V**5-1327*V**4+1396*V**3-943*V**2+396*V-&
&66)*W**4+(V-1)*(32*V**10+288*V**9-793*V**8+13*V**7+1340*V**6-1476*V**5&
&+683*V**4-49*V**3-134*V**2+76*V-12)*W**3-(V-1)**2*(96*V**8-667*V**6+99&
&2*V**5-496*V**4-62*V**3+125*V**2-46*V+10)*W**2+(V-1)**3*(96*V**6-224*V&
&**5+53*V**4+220*V**3-233*V**2+66*V-10)*W-8*(V-1)**4*(4*V**4-12*V**3+12&
&*V**2-4*V-1))-2*N**2*(V-1)*VC*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)*(V**3*W&
&**3-3*V**2*W**2+3*V*W-1)*((V-1)*V**6*(3*V**2-6*V+10)*W**8-(V-1)*V**5*(&
&10*V**3-24*V**2+41*V-10)*W**7+V**4*(15*V**5-60*V**4+125*V**3-77*V**2-1&
&9*V+20)*W**6-V**3*(15*V**6-80*V**5+198*V**4-119*V**3-80*V**2+94*V-20)*&
&W**5+V**2*(10*V**7-55*V**6+143*V**5-52*V**4-145*V**3+123*V**2-10*V-10)&
&*W**4-(V-1)*V*(3*V**7-15*V**6+62*V**5-17*V**4-30*V**3-2*V**2+25*V-10)*&
&W**3+(V-1)**2*V*(23*V**4-27*V**3+54*V**2-39*V+13)*W**2-(V-1)**3*(7*V**&
&4-18*V**3+39*V**2-20*V+8)*W+4*(V-1)**4)-2*AL*N**4*(V-1)*VC*(V**2*W**2-&
&2*V*W+1)*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**&
&3*V**2*W**2+5*(V-1)**4*V*W-(V-1)**5)*(V**3*(3*V**4-V**3+4*V**2+2)*W**6&
&-2*V**2*(2*V**5+4*V**4+6*V**2+V+3)*W**5+V*(3*V**6+10*V**5+11*V**4+3*V*&
&*3+13*V**2+6*V+6)*W**4-(11*V**6+6*V**5+12*V**4+7*V**2+6*V+2)*W**3+(19*&
&V**5-14*V**4+18*V**3-8*V**2+3*V+2)*W**2-(17*V**4-24*V**3+18*V**2-8*V+1&
&)*W+4*(V-1)*(2*V**2-2*V+1))+4*AL*N**2*(V-1)*VC*(V**2*W**2-2*V*W+1)*(V*&
&*5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V**2*W**2+&
&5*(V-1)**4*V*W-(V-1)**5)*(V**3*(3*V**4-3*V**3+4*V**2-2*V+2)*W**6-V**2*&
&(4*V**5+5*V**4-5*V**3+10*V**2-4*V+6)*W**5+V*(3*V**6+7*V**5+4*V**4-V**3&
&+7*V**2+6)*W**4-(9*V**6+V**5+5*V**4+V**2+4*V+2)*W**3+(13*V**5-13*V**4+&
&13*V**3-6*V**2+V+2)*W**2-(9*V**4-13*V**3+10*V**2-5*V+1)*W+2*(V-1)*(2*V&
&**2-2*V+1))+(V-1)**2*VC*W*(V**2*W**2-2*V*W+1)*(V**4*W**4-4*(V-1)*V**3*&
&W**3+6*(V-1)**2*V**2*W**2-4*(V-1)**3*V*W+(V-1)**4)*(V**4*(V**3+2*V**2-&
&12*V+12)*W**6-V**3*(V+1)*(V**3+6*V**2-26*V+24)*W**5+V**3*(V**4+11*V**3&
&+6*V**2-56*V+56)*W**4-V*(V**6+5*V**5+20*V**4-22*V**3-22*V**2+60*V-24)*&
&W**3+(5*V**6-23*V**5+78*V**4-85*V**3+40*V**2+4*V-12)*W**2+(V-1)*(5*V**&
&4-8*V**3-29*V**2+24*V-10)*W+(V-1)**2*(3*V**2-2*V+6))-2*AL*(V-1)**2*VC*&
&W*(V**2*W**2-2*V*W+1)*(V**3*(3*V**3-2*V**2+2*V-2)*W**5-2*V**2*(2*V**4+&
&3*V**3-2*V**2+2*V-3)*W**4+V*(3*V**5+7*V**4+4*V**3-V**2-6)*W**3-(7*V**5&
&+3*V**4+V**3+V**2-4*V-2)*W**2+(7*V**4-5*V**3+3*V**2-V-2)*W-(V-1)*(V**2&
&+1))*(V**5*W**5-5*(V-1)*V**4*W**4+10*(V-1)**2*V**3*W**3-10*(V-1)**3*V*&
&*2*W**2+5*(V-1)**4*V*W-(V-1)**5))/(N**2*(V-1)**3*V**2*W**2*(V*W-1)**5*&
&(V*W-V+1)**5)/2.D0
LM = LOG(S/M**2)*(N**4*VC*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V*W-&
&(V-1)**3)*(16*V**8*W**8-16*V**7*(V+4)*W**7+V**3*(48*V**5-45*V**4+207*V&
&**3-44*V**2+20*V-2)*W**6-2*V**2*(8*V**6+50*V**5-76*V**4+176*V**3-56*V*&
&*2+29*V-3)*W**5+V*(16*V**7-13*V**6+138*V**5-149*V**4+263*V**3-71*V**2+&
&54*V-6)*W**4-(16*V**8-37*V**6+134*V**5-72*V**4+64*V**3+19*V**2+14*V-2)&
&*W**3+(48*V**7-96*V**6+83*V**5+6*V**4-10*V**3-4*V**2+23*V-2)*W**2-(48*&
&V**6-128*V**5+165*V**4-112*V**3+42*V**2-8*V+1)*W+8*(V-1)*(V**2-V+1)*(2&
&*V**2-2*V+1))-2*N**2*VC*(V*W-1)*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**&
&2*V*W-(V-1)**3)*(V**3*(3*V**3+5*V**2-4*V+4)*W**5-2*V**2*(2*V**4+5*V**3&
&+3*V**2-2*V+4)*W**4+V*(3*V**5+11*V**4+2*V**3+3*V**2+5*V+4)*W**3-V*(14*&
&V**4-10*V**3+9*V**2-3*V+6)*W**2+(3*V**2-3*V+1)*(5*V**2-2*V+1)*W-4*(V-1&
&)*(2*V**2-2*V+1))+(V-1)*VC*W*(V**3*W**3-3*(V-1)*V**2*W**2+3*(V-1)**2*V&
&*W-(V-1)**3)*(V**3*(3*V**3-2*V**2+2*V-2)*W**5-2*V**2*(2*V**4+3*V**3-2*&
&V**2+2*V-3)*W**4+V*(3*V**5+7*V**4+4*V**3-V**2-6)*W**3-(7*V**5+3*V**4+V&
&**3+V**2-4*V-2)*W**2+(7*V**4-5*V**3+3*V**2-V-2)*W-(V-1)*(V**2+1)))/(N*&
&*2*(V-1)**2*V**2*W**2*(V*W-1)**3*(V*W-V+1)**3)
LMP = LOG(S/MP**2)*(N**4*VC*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(16*V**8*W&
&**8-16*V**7*(5*V-4)*W**7+V**6*(208*V**2-351*V+151)*W**6-(V-1)*V**5*(35&
&2*V**2-559*V+230)*W**5+2*(V-1)**2*V**3*(204*V**3-307*V**2+123*V+1)*W**&
&4-2*(V-1)**3*V**2*(168*V**3-235*V**2+90*V+3)*W**3+(V-1)**4*V*(200*V**3&
&-255*V**2+95*V+6)*W**2-(V-1)**5*(80*V**3-95*V**2+38*V+2)*W+8*(V-1)**6*&
&(2*V**2-2*V+1))-2*N**2*(V-1)*VC*W*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**&
&6*W**5-V**5*(2*V-3)*W**4-(V-1)*V**3*(4*V**3-7*V**2+10*V-2)*W**3+(V-1)*&
&*2*V**2*(8*V**3-15*V**2+17*V-6)*W**2-(V-1)**3*V*(4*V**3-10*V**2+13*V-6&
&)*W-(V-1)**4*(V**2-2*V+2))+(V-1)*VC*W*(V**2*W**2-2*(V-1)*V*W+(V-1)**2)&
&*(V**3*W**3-3*V**2*W**2+3*V*W-1)*(V**4*W**3-(V-2)*V**3*W**2+(V-1)*V*(V&
&**2-5*V+2)*W-(V-1)**2*(V**2-2*V+2)))/(N**2*(V-1)**2*V**2*W**2*(V*W-1)*&
&*3*(V*W-V+1)**3)
STRUV16=LV1+LW1+LV+LW+LVW+LTVW+CVC+LM+LMP
RETURN
END FUNCTION
FUNCTION STRUVC(W,V,X3,S)
double precision::FAC1
double precision::FAC2
double precision::FAC3
double precision::LTOT
double precision::PGG
double precision::PGQ
double precision::PQG
double precision::PQQ
double precision::S
double precision::STRUVC
double precision::T
double precision::U
double precision::V
double precision::W
double precision::X3
double precision::Y
double precision::Z
LTOT=DLOG(V**2*(1.D0-W)**2*PT**2*DELTA**2/Q2FRAG/X3**2)
Z=1.D0-V+V*W
Y=V*W/(1.D0-V+V*W)
T=-S*(1.D0-Y)
U=-S*Y
FAC1=-4.D0*N**2/(1.D0-V)/W/(1.D0-V+V*W)
FAC2=2.D0*CF*FAC1
FAC3=2.D0*CF*FAC2
PQQ=(-LTOT)*CF*(1.D0+Z**2)/(1.D0-Z) - CF*(1.D0-Z)
PQG=(-LTOT)*1.D0/2.D0*(Z**2+(1.D0-Z)**2) - Z*(1.D0-Z)
PGQ=(-LTOT)*CF*(1.D0+(1.D0-Z)**2)/Z - CF*Z
PGG=(-LTOT)*2.D0*N*( (1.D0-Z)/Z+Z/(1.D0-Z)+Z*(1.D0-Z) )
IF (J0.EQ.1) THEN
STRUVC=FAC1*( CF/N*(S**2+U**2)/T**2 * PQQ )
ELSE IF (J0.EQ.2) THEN
STRUVC=FAC1*( CF/N*(S**2+U**2)/T**2 * PGQ +CF/N*(S**2+T**2)/U**2 * PGQ&
&)
ELSE IF (J0.EQ.3) THEN
STRUVC=FAC1*( CF/N*(S**2+U**2)/T**2 * PQQ )
ELSE IF (J0.EQ.4) THEN
STRUVC=FAC1*( CF/N*(S**2+U**2)/T**2 * PGQ +CF/N*(S**2+T**2)/U**2 * PGQ&
&)
ELSE IF (J0.EQ.5) THEN
STRUVC=FAC1*( CF/N*(T**2+U**2)/S**2 * PQQ +2.D0*CF*(T**2+U**2)*(CF/N/T&
&/U-1.D0/S**2) * PQG )
ELSE IF (J0.EQ.6) THEN
STRUVC=FAC1*( CF/N*(T**4+U**4+S**2*(T**2+U**2) -2.D0/N*U*T*S**2)/U**2/&
&T**2 * PQQ )
ELSE IF (J0.EQ.7) THEN
STRUVC=FAC1*( CF/N*(T**4+U**4+S**2*(T**2+U**2) -2.D0/N*U*T*S**2)/U**2/&
&T**2 * PGQ )
ELSE IF (J0.EQ.8) THEN
STRUVC=FAC2*( (1.D0/U**2-CF/N/S/T)*(S**2+T**2) * PQG )
ELSE IF (J0.EQ.9) THEN
STRUVC=FAC2*( (1.D0/U**2-CF/N/S/T)*(S**2+T**2) * PQG )
ELSE IF (J0.EQ.10) THEN
STRUVC=FAC2*( (1.D0/U**2-CF/N/S/T)*(S**2+T**2) * PQG )
ELSE IF (J0.EQ.11) THEN
STRUVC=FAC1*( 2.D0*CF*(T**2+U**2)*(CF/N/T/U-1.D0/S**2) * PQG +CF/N*(S*&
&*4+T**4+U**2*(S**2+T**2) -2.D0/N*S*T*U**2)/S**2/T**2 * PQQ )
ELSE IF (J0.EQ.12) THEN
STRUVC=FAC1*( 2.D0*CF*(T**2+U**2)*(CF/N/T/U-1.D0/S**2) * PGG +CF/N*(S*&
&*4+T**4+U**2*(S**2+T**2) -2.D0/N*S*T*U**2)/S**2/T**2 * PGQ +CF/N*(S**4&
&+U**4+T**2*(S**2+U**2) -2.D0/N*S*U*T**2)/S**2/U**2 * PGQ +CF/N*(T**2+U&
&**2)/S**2*2.D0*(NF-1.D0) * PGQ )
ELSE IF (J0.EQ.13) THEN
STRUVC=FAC2*( (1.D0/U**2-CF/N/S/T)*(S**2+T**2) * PQG +(1.D0/T**2-CF/N/&
&S/U)*(S**2+U**2) * PQQ  )
ELSE IF (J0.EQ.14) THEN
STRUVC=FAC2*( (1.D0/U**2-CF/N/S/T)*(S**2+T**2) * PGG +(1.D0/T**2-CF/N/&
&S/U)*(S**2+U**2) * PGQ  )
ELSE IF (J0.EQ.15) THEN
STRUVC=FAC3*( N/2.D0/CF*(S**4+T**4+U**4)*(S**2+T**2+U**2)/S**2/T**2/U*&
&*2 * PGG +1.D0/2.D0/CF*(CF/N/T/U-1.D0/S**2)*(T**2+U**2) * PGQ *2.D0*NF&
&)
ELSE IF (J0.EQ.16) THEN
STRUVC=FAC3*( N/2.D0/CF*(S**4+T**4+U**4)*(S**2+T**2+U**2)/S**2/T**2/U*&
&*2 * PQG +1.D0/2.D0/CF*(CF/N/T/U-1.D0/S**2)*(T**2+U**2) * PQQ )
ENDIF
RETURN
END FUNCTION
FUNCTION HQQD(V,J)
double precision::HQQD
integer::J
double precision::UN
double precision::V
UN=1.D0
IF(J.EQ.1) THEN
HQQD=CQQD(UN)
ELSE
HQQD=(1.-V)/V*(CQQD(UN)+DLOG(V/(1.-V))*CQQW(UN)+DLOG((1.-V)/V)**2/2.*C&
&QQL(UN))
ENDIF
RETURN
END FUNCTION
FUNCTION HGQD(V,J)
double precision::HGQD
integer::J
double precision::UN
double precision::V
UN=1.D0
IF(J.EQ.1) THEN
HGQD=CGQD(UN)
ELSE
HGQD=(1.-V)/V*(CGQD(UN)+DLOG(V/(1.-V))*CGQW(UN)+DLOG((1.-V)/V)**2/2.*C&
&GQL(UN))
ENDIF
RETURN
END FUNCTION
FUNCTION HQGD(V,J)
double precision::HQGD
integer::J
double precision::UN
double precision::V
UN=1.D0
IF(J.EQ.1) THEN
HQGD=CQGD(UN)
ELSE
HQGD=(1.-V)/V*(CQGD(UN)+DLOG(V/(1.-V))*CQGW(UN)+DLOG((1.-V)/V)**2/2.*C&
&QGL(UN))
ENDIF
RETURN
END FUNCTION
FUNCTION HGGD(V,J)
double precision::HGGD
integer::J
double precision::UN
double precision::V
UN=1.D0
IF(J.EQ.1) THEN
HGGD=CGGD(UN)
ELSE
HGGD=(1.-V)/V*(CGGD(UN)+DLOG(V/(1.-V))*CGGW(UN)+DLOG((1.-V)/V)**2/2.*C&
&GGL(UN))
ENDIF
RETURN
END FUNCTION
FUNCTION HQQW(W,V,J)
double precision::HQQW
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HQQW=CQQW(W)
ELSE
HQQW=(1.-V*W)/V*(CQQW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CQQL((1.-V)/(1&
&.-V*W)))
ENDIF
RETURN
END FUNCTION
FUNCTION HQGW(W,V,J)
double precision::HQGW
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HQGW=CQGW(W)
ELSE
HQGW=(1.-V*W)/V*(CQGW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CQGL((1.-V)/(1&
&.-V*W)))
ENDIF
RETURN
END FUNCTION
FUNCTION HGGW(W,V,J)
double precision::HGGW
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HGGW=CGGW(W)
ELSE
HGGW=(1.-V*W)/V*(CGGW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CGGL((1.-V)/(1&
&.-V*W)))
ENDIF
RETURN
END FUNCTION
FUNCTION HGQW(W,V,J)
double precision::HGQW
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HGQW=CGQW(W)
ELSE
HGQW=(1.-V*W)/V*(CGQW((1.-V)/(1.-V*W))+DLOG(V/(1.-V*W))*CGQL((1.-V)/(1&
&.-V*W)))
ENDIF
RETURN
END FUNCTION
FUNCTION HGQL(W,V,J)
integer::J
double precision::HGQL
double precision::V
double precision::W
IF(J.EQ.1) THEN
HGQL=CGQL(W)
ELSE
HGQL=(1.-V*W)/V*CGQL((1.-V)/(1.-V*W))
ENDIF
RETURN
END FUNCTION
FUNCTION HGGL(W,V,J)
double precision::HGGL
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HGGL=CGGL(W)
ELSE
HGGL=(1.-V*W)/V*CGGL((1.-V)/(1.-V*W))
ENDIF
RETURN
END FUNCTION
FUNCTION HQQL(W,V,J)
double precision::HQQL
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HQQL=CQQL(W)
ELSE
HQQL=(1.-V*W)/V*CQQL((1.-V)/(1.-V*W))
ENDIF
RETURN
END FUNCTION
FUNCTION HQGL(W,V,J)
double precision::HQGL
integer::J
double precision::V
double precision::W
IF(J.EQ.1) THEN
HQGL=CQGL(W)
ELSE
HQGL=(1.-V*W)/V*CQGL((1.-V)/(1.-V*W))
ENDIF
RETURN
END FUNCTION
FUNCTION HFQQL(W,V)
double precision::HFQQL
double precision::V
double precision::W
HFQQL=FQQL(1.-V+V*W)/V
RETURN
END FUNCTION
FUNCTION HFQQW(W,V)
double precision::HFQQW
double precision::V
double precision::W
HFQQW=1./V*(FQQW(1.-V+V*W)+DLOG(V)*FQQL(1.-V+V*W))
RETURN
END FUNCTION
FUNCTION HFQQD(V)
double precision::HFQQD
double precision::UN
double precision::V
UN=1.D0
HFQQD=1./V*(FQQD(UN)+DLOG(V)*FQQW(UN)+DLOG(V)**2/2.*FQQL(UN))
RETURN
END FUNCTION
FUNCTION HFQGL(W,V)
double precision::HFQGL
double precision::V
double precision::W
HFQGL=FQGL(1.-V+V*W)/V
RETURN
END FUNCTION
FUNCTION HFQGW(W,V)
double precision::HFQGW
double precision::V
double precision::W
HFQGW=1./V*(FQGW(1.-V+V*W)+DLOG(V)*FQGL(1.-V+V*W))
RETURN
END FUNCTION
FUNCTION HFQGD(V)
double precision::HFQGD
double precision::UN
double precision::V
UN=1.D0
HFQGD=1./V*(FQGD(UN)+DLOG(V)*FQGW(UN)+DLOG(V)**2/2.*FQGL(UN))
RETURN
END FUNCTION
FUNCTION HFGGL(W,V)
double precision::HFGGL
double precision::V
double precision::W
HFGGL=FGGL(1.-V+V*W)/V
RETURN
END FUNCTION
FUNCTION HFGGW(W,V)
double precision::HFGGW
double precision::V
double precision::W
HFGGW=1./V*(FGGW(1.-V+V*W)+DLOG(V)*FGGL(1.-V+V*W))
RETURN
END FUNCTION
FUNCTION HFGGD(V)
double precision::HFGGD
double precision::UN
double precision::V
UN=1.D0
HFGGD=1./V*(FGGD(UN)+DLOG(V)*FGGW(UN)+DLOG(V)**2/2.*FGGL(UN))
RETURN
END FUNCTION
FUNCTION HFGQL(W,V)
double precision::HFGQL
double precision::V
double precision::W
HFGQL=FGQL(1.-V+V*W)/V
RETURN
END FUNCTION
FUNCTION HFGQW(W,V)
double precision::HFGQW
double precision::V
double precision::W
HFGQW=1./V*(FGQW(1.-V+V*W)+DLOG(V)*FGQL(1.-V+V*W))
RETURN
END FUNCTION
FUNCTION HFGQD(V)
double precision::HFGQD
double precision::UN
double precision::V
UN=1.D0
HFGQD=1./V*(FGQD(UN)+DLOG(V)*FGQW(UN)+DLOG(V)**2/2.*FGQL(UN))
RETURN
END FUNCTION
FUNCTION A(S,T,U)
double precision::A
double precision::S
double precision::T
double precision::U
double precision::VC
VC=(N**2-1.)
A=2.*VC*(S**2+U**2)/T**2
RETURN
END FUNCTION
FUNCTION B(S,T,U)
double precision::B
double precision::N
double precision::S
double precision::T
double precision::U
double precision::VC
N=3.D0
VC=(N**2-1.)
B=-4.*VC/N*S**2/U/T
RETURN
END FUNCTION
FUNCTION C(S,T,U)
double precision::C
double precision::N
double precision::S
double precision::T
double precision::U
double precision::VC
N=3.D0
VC=(N**2-1.)
C=2.*VC/N*(VC/U/T-2.*N**2/S**2)*(T**2+U**2)
RETURN
END FUNCTION
FUNCTION D(S,T,U)
double precision::D
double precision::N
double precision::S
double precision::T
double precision::U
double precision::VC
N=3.D0
VC=(N**2-1.)
D=16.*VC*N**2*(3.-U*T/S**2-U*S/T**2-S*T/U**2)
RETURN
END FUNCTION
FUNCTION A0(X,S)
double precision::A0
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
A0=A(S,T,U)
RETURN
END FUNCTION
FUNCTION A2(X,S)
double precision::A2
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
A2=A(T,S,U)
RETURN
END FUNCTION
FUNCTION B0(X,S)
double precision::B0
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
B0=A(S,T,U)+A(S,U,T)+B(S,T,U)
RETURN
END FUNCTION
FUNCTION D0(X,S)
double precision::D0
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
D0=A(U,T,S)+A(U,S,T)+B(U,T,S)
RETURN
END FUNCTION
FUNCTION D1(X,S)
double precision::D1
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
D1=C(S,T,U)
RETURN
END FUNCTION
FUNCTION E0(X,S)
double precision::E0
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
E0=-C(T,S,U)
RETURN
END FUNCTION
FUNCTION F2(X,S)
double precision::F2
double precision::S
double precision::T
double precision::U
double precision::X
T=-S*(1.-X)
U=-S*X
F2=D(S,T,U)
RETURN
END FUNCTION
SUBROUTINE STRU(XUHA,XUBHA,XDHA,XDBHA,XSHA,XCHA,XGPROA,XUHB,XUBHB,XDHB&
&,XDBHB,XSHB,XCHB,XGPROB,XDUP,XDUBP,XDDP,XDDBP,XDSP,XDSBP,XDCP,XDCBP,XD&
&GP,GPPV,GPPC)
double precision::GPPC
double precision::GPPV
double precision::XCHA
double precision::XCHB
double precision::XDBHA
double precision::XDBHB
double precision::XDCBP
double precision::XDCP
double precision::XDDBP
double precision::XDDP
double precision::XDGP
double precision::XDHA
double precision::XDHB
double precision::XDSBP
double precision::XDSP
double precision::XDUBP
double precision::XDUP
double precision::XGPROA
double precision::XGPROB
double precision::XSHA
double precision::XSHB
double precision::XUBHA
double precision::XUBHB
double precision::XUHA
double precision::XUHB
DIMENSION GPPV(16),GPPC(16)
GPPV(1)=&
& XUHA*(XDHB+XSHB+XCHB)*XDUP&
&+XDHA*(XUHB+XSHB+XCHB)*XDDP&
&+XSHA*(XUHB+XDHB+XCHB)*XDSP&
&+XCHA*(XUHB+XDHB+XSHB)*XDCP&
&+XUBHA*(XDBHB+XSHB+XCHB)*XDUBP&
&+XDBHA*(XUBHB+XSHB+XCHB)*XDDBP&
&+XSHA*(XUBHB+XDBHB+XCHB)*XDSP&
&+XCHA*(XUBHB+XDBHB+XSHB)*XDCP


GPPC(1)=XUHB*(XDHA+XSHA+XCHA)*XDUP+XDHB*(XUHA+XSHA+XCHA)*XDDP+XSHB*(XU&
&HA+XDHA+XCHA)*XDSP+XCHB*(XUHA+XDHA+XSHA)*XDCP+XUBHB*(XDBHA+XSHA+XCHA)*&
&XDUBP+XDBHB*(XUBHA+XSHA+XCHA)*XDDBP+XSHB*(XUBHA+XDBHA+XCHA)*XDSP+XCHB*&
&(XUBHA+XDBHA+XSHA)*XDCP
GPPV(2)=(XUHA*(XDHB+XSHB+XCHB)+XDHA*(XSHB+XCHB)+XSHA*XCHB+XUBHA*(XDBH&
&B+XSHB+XCHB)+XDBHA*(XSHB+XCHB)+XSHA*XCHB)*XDGP
GPPC(2)=(XUHB*(XDHA+XSHA+XCHA)+XDHB*(XSHA+XCHA)+XSHB*XCHA+XUBHB*(XDBH&
&A+XSHA+XCHA)+XDBHB*(XSHA+XCHA)+XSHB*XCHA)*XDGP
GPPV(3)=XUHA*(XDBHB+XSHB+XCHB)*XDUP+XDHA*(XUBHB+XSHB+XCHB)*XDDP+XSHA*(&
&XUBHB+XDBHB+XCHB)*XDSP+XCHA*(XUBHB+XDBHB+XSHB)*XDCP+XUBHA*(XDHB+XSHB+X&
&CHB)*XDUBP+XDBHA*(XUHB+XSHB+XCHB)*XDDBP+XSHA*(XUHB+XDHB+XCHB)*XDSP+XCH&
&A*(XUHB+XDHB+XSHB)*XDCP
GPPC(3)=XUHB*(XDBHA+XSHA+XCHA)*XDUP+XDHB*(XUBHA+XSHA+XCHA)*XDDP+XSHB*(&
&XUBHA+XDBHA+XCHA)*XDSP+XCHB*(XUBHA+XDBHA+XSHA)*XDCP+XUBHB*(XDHA+XSHA+X&
&CHA)*XDUBP+XDBHB*(XUHA+XSHA+XCHA)*XDDBP+XSHB*(XUHA+XDHA+XCHA)*XDSP+XCH&
&B*(XUHA+XDHA+XSHA)*XDCP
GPPV(4)=(XUHA*(XDBHB+XSHB+XCHB)+XDHA*(XUBHB+XSHB+XCHB)+XSHA*(XUBHB+XDB&
&HB+XCHB)+XCHA*(XUBHB+XDBHB+XSHB))*XDGP
GPPC(4)=(XUHB*(XDBHA+XSHA+XCHA)+XDHB*(XUBHA+XSHA+XCHA)+XSHB*(XUBHA+XDB&
&HA+XCHA)+XCHB*(XUBHA+XDBHA+XSHA))*XDGP
GPPV(5)=(XDHA*XDBHB+XSHA*XSHB+XCHA*XCHB)*XDUP+(XUHA*XUBHB+XSHA*XSHB+XC&
&HA*XCHB)*XDDP+(XUHA*XUBHB+XDHA*XDBHB+XCHA*XCHB)*XDSP+(XUHA*XUBHB+XDHA*&
&XDBHB+XSHA*XSHB)*XDCP+(XDBHA*XDHB+XSHA*XSHB+XCHA*XCHB)*XDUBP+(XUBHA*XU&
&HB+XSHA*XSHB+XCHA*XCHB)*XDDBP+(XUBHA*XUHB+XDBHA*XDHB+XCHA*XCHB)*XDSP+(&
&XUBHA*XUHB+XDBHA*XDHB+XSHA*XSHB)*XDCP
GPPC(5)=(XDHB*XDBHA+XSHB*XSHA+XCHB*XCHA)*XDUP+(XUHB*XUBHA+XSHB*XSHA+XC&
&HB*XCHA)*XDDP+(XUHB*XUBHA+XDHB*XDBHA+XCHB*XCHA)*XDSP+(XUHB*XUBHA+XDHB*&
&XDBHA+XSHB*XSHA)*XDCP+(XDBHB*XDHA+XSHB*XSHA+XCHB*XCHA)*XDUBP+(XUBHB*XU&
&HA+XSHB*XSHA+XCHB*XCHA)*XDDBP+(XUBHB*XUHA+XDBHB*XDHA+XCHB*XCHA)*XDSP+(&
&XUBHB*XUHA+XDBHB*XDHA+XSHB*XSHA)*XDCP
GPPV(6)=XUHA*XUHB*XDUP+XDHA*XDHB*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP+XU&
&BHA*XUBHB*XDUBP+XDBHA*XDBHB*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
GPPV(6)=GPPV(6)/2.D0
GPPC(6)=GPPV(6)
GPPV(7)=(XUHA*XUHB+XDHA*XDHB+XSHA*XSHB+XCHA*XCHB+XUBHA*XUBHB+XDBHA*XDB&
&HB+XSHA*XSHB+XCHA*XCHB)*XDGP
GPPV(7)=GPPV(7)/2.D0
GPPC(7)=GPPV(7)
GPPV(8)=((XDHA+XSHA+XCHA)*XDUP+(XUHA+XSHA+XCHA)*XDDP+(XUHA+XDHA+XCHA)*&
&XDSP+(XUHA+XDHA+XSHA)*XDCP+(XDBHA+XSHA+XCHA)*XDUP+(XUBHA+XSHA+XCHA)*XD&
&DP+(XUBHA+XDBHA+XCHA)*XDSP+(XUBHA+XDBHA+XSHA)*XDCP)*XGPROB
GPPC(8)=((XDHB+XSHB+XCHB)*XDUP+(XUHB+XSHB+XCHB)*XDDP+(XUHB+XDHB+XCHB)*&
&XDSP+(XUHB+XDHB+XSHB)*XDCP+(XDBHB+XSHB+XCHB)*XDUP+(XUBHB+XSHB+XCHB)*XD&
&DP+(XUBHB+XDBHB+XCHB)*XDSP+(XUBHB+XDBHB+XSHB)*XDCP)*XGPROA
GPPV(9)=GPPV(8)
GPPC(9)=GPPC(8)
GPPV(10)=(XUHA*XDUBP+XDHA*XDDBP+XSHA*XDSP+XCHA*XDCP+XUBHA*XDUP+XDBHA*X&
&DDP+XSHA*XDSP+XCHA*XDCP)*XGPROB
GPPC(10)=(XUHB*XDUBP+XDHB*XDDBP+XSHB*XDSP+XCHB*XDCP+XUBHB*XDUP+XDBHB*X&
&DDP+XSHB*XDSP+XCHB*XDCP)*XGPROA
GPPV(11)=XUHA*XUBHB*XDUP+XDHA*XDBHB*XDDP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP&
&+XUBHA*XUHB*XDUBP+XDBHA*XDHB*XDDBP+XSHA*XSHB*XDSP+XCHA*XCHB*XDCP
GPPC(11)=XUHB*XUBHA*XDUP+XDHB*XDBHA*XDDP+XSHB*XSHA*XDSP+XCHB*XCHA*XDCP&
&+XUBHB*XUHA*XDUBP+XDBHB*XDHA*XDDBP+XSHB*XSHA*XDSP+XCHB*XCHA*XDCP
GPPV(12)=(XUHA*XUBHB+XDHA*XDBHB+XSHA*XSHB+XCHA*XCHB)*XDGP
GPPC(12)=(XUHB*XUBHA+XDHB*XDBHA+XSHB*XSHA+XCHB*XCHA)*XDGP
GPPV(13)=(XUHA*XDUP+XUBHA*XDUBP+XDHA*XDDP+XDBHA*XDDBP+2.*XSHA*XDSP+2.*&
&XCHA*XDCP)*XGPROB
GPPC(13)=(XUHB*XDUP+XUBHB*XDUBP+XDHB*XDDP+XDBHB*XDDBP+2.*XSHB*XDSP+2.*&
&XCHB*XDCP)*XGPROA
GPPV(14)=(XUHA+XUBHA+XDHA+XDBHA+2.*XSHA+2.*XCHA)*XGPROB*XDGP
GPPC(14)=(XUHB+XUBHB+XDHB+XDBHB+2.*XSHB+2.*XCHB)*XGPROA*XDGP
GPPV(15)=XGPROA*XGPROB*XDGP/2.
GPPC(15)=GPPV(15)
GPPV(16)=XGPROA*XGPROB*(XDUP+XDDP+XDSP+XDCP)
GPPC(16)=GPPV(16)
RETURN
END SUBROUTINE

end module




