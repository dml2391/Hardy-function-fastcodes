!  Routine for calculating the Hardy function for Z(t) for large t, typically t>10^15.
!
!  In this version, the Main Sum of the Riemann-Siegel Formula (MSRSF) making up Z(t),
!  is formulated in terms of quartic Gauss sums (GS), rather than cubic Gauss sums. Since the former are comprised of many
!  more summands than the latter, the number of quartic GS are required to represent the MSRSF is vastly reduced, which in
!  turn implies that the number of operations to calculate Z(t) should also reduced. Unfortunately the iterative proceedure
!  designed to reduce an individual cubic GS down to a modest kernel sum of perhaps a few hundred summands, works far less
!  efficiently for quartic GS. Typically when the hierarchical reduction scheme has reduced the original quartic GS down to a kernel
!  of some thousand summands, 'technical difficulties' arise, preventing any further reductions in the length of the kernel.
!  The 'technical difficulties' arise when the third sum coefficient, phi3, gets too large with respect to length of the kernel L and
!  the second coefficient phi2. In particular abs(3*phi3*L)~2*phi2. Phi3 always increases with respect to phi2 as one moves down
!  the hierarchical chain, but for cubic sums the relative increase is slower, meaning that 2*phi2>>3*phi3*L all way down the chain and
!  the original sum reduced can be reduced  to a suitably small kernel. If these technical difficulties can be overcome then one
!  could create an algorithm based upon mth order GS which would require only O(t^(1/(m+1)) operations for a single Z(t) computation. 
!
!  It is possible to move one step further down the chain, to a kernel somewhat smaller than a few thousand summands. However,
!  the phase of this reduced kernel is made up of saddle points that cannot be represented in terms of the coefficients of phi1,phi2,
!  phi3...phiq in a closed form. So they have to be evaluated numerically. So whilst on paper the reduced kernel may have around a
!  thousand summands, its evaluation is considerably more expensive, because of the extra numerical evalation of each saddle associated
!  with an individual summand (and sometimes a high value summand involves evalating two separate "double" saddle points). This extra
!  numerical saddle point evaluation is implemented here, since the benefits of reducing the kernel, tend to outweigh the extra operations
!  needed to calculate a (real) polynomial root within a fixed range. 
!
!  These 'technical difficulties' become more acute as t increases. So whilst in theory a quartic sum algorithm should allow Z(t)
!  to be evalated in O(t^(1/5)) operations - as mentioned above - in practice this increases back up to O(t^(1/4)*(log(t))^3) operations,
!  very similar to the cubic case. In practical computations, an epsilon setting around 0.0002 (as compared to 0.005 for the cubic case),
!  is sufficient to ensure that this code runs significantly faster for t values up to around 10^28. Beyond that, the time savings dwindle
!  and the performance of the two algorithms become pretty similar. Indeed one is better off using the cubic code, since it tends to be a
!  bit more accurate.

!  Full details of the various subroutines are given in the "HardyCode-Detailed Documentation.pdf", file which can be downloaded from
!  the website. The coding comments included here are designed to be read in conjugation with routine descriptions to found in the
!  aforementioned file.


module PARI
  use ISO_C_BINDING, only : C_LONG, C_DOUBLE, C_PTR
  interface

     subroutine pari_init(parisize, maxprime) bind(C,name='pari_init')
       import C_LONG
       integer(kind=C_LONG), VALUE              :: parisize
       integer(kind=C_LONG), VALUE              :: maxprime
     end subroutine pari_init
     !
     subroutine pari_close() bind(C,name='pari_close')
     end subroutine pari_close
     !
     type(C_PTR) function dbltor( r ) bind(C,name='dbltor')
       import C_DOUBLE, C_PTR
       real(kind=C_DOUBLE), VALUE  :: r
     end function dbltor
     !
     real(kind=C_DOUBLE) function rtodbl( x ) bind(C,name='rtodbl')
       import C_DOUBLE, C_PTR
       type(C_PTR), VALUE :: x
     end function rtodbl
     !
     type(C_PTR) function gsqr( x ) bind(C,name='gsqr')
       import C_PTR
       type(C_PTR), VALUE :: x
     end function gsqr
     !
     type(C_PTR) function gfloor( x ) bind(C,name='gfloor')
       import C_PTR
       type(C_PTR), VALUE :: x
     end function gfloor
     !
     type(C_PTR) function gceil( x ) bind(C,name='gceil')
       import C_PTR
       type(C_PTR), VALUE :: x
     end function gceil
     !
     integer(kind=C_LONG) function sizedigit( x ) bind(C,name='sizedigit')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
     end function sizedigit
     !
     integer function gtos( x ) bind(C,name='gtos')
       import C_PTR
       type(C_PTR), VALUE :: x
     end function gtos
     !
     type(C_PTR) function gfrac( x ) bind(C,name='gfrac')
       import C_PTR
       type(C_PTR), VALUE :: x
     end function gfrac
     !
     type(C_PTR) function gadd( x , y) bind(C,name='gadd')
       import C_PTR
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: y
     end function gadd
     !
     type(C_PTR) function gsub( x , y) bind(C,name='gsub')
       import C_PTR
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: y
     end function gsub     
     !
     type(C_PTR) function gmul( x , y) bind(C,name='gmul')
       import C_PTR
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: y
     end function gmul
     !
     type(C_PTR) function gdiv( x , y) bind(C,name='gdiv')
       import C_PTR
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: y
     end function gdiv
     !
     type(C_PTR) function gpow(x, n, prec) bind(C,name='gpow')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: n
       integer(kind=C_LONG), VALUE :: prec
     end function gpow
     !
     type(C_PTR) function gprec( x , d) bind(C,name='gprec')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: d
     end function gprec
     !
     type(C_PTR) function gmod( x , y) bind(C,name='gmod')
       import C_PTR
       type(C_PTR), VALUE :: x
       type(C_PTR), VALUE :: y
     end function gmod
     !
     type(C_PTR) function glog( x , prec) bind(C,name='glog')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function glog
     !
     type(C_PTR) function gatan( x , prec) bind(C,name='gatan')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function gatan
     !
     type(C_PTR) function gerfc( x , prec) bind(C,name='gerfc')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function gerfc
     !
     type(C_PTR) function gexp( x , prec) bind(C,name='gexp')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function gexp
     !
     type(C_PTR) function glngamma( x , prec) bind(C,name='glngamma')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function glngamma
     !
     type(C_PTR) function gsqrt( x , prec) bind(C,name='gsqrt')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function gsqrt
     !
     type(C_PTR) function gimag( x ) bind(C,name='gimag')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
     end function gimag
     !
     type(C_PTR) function greal( x ) bind(C,name='greal')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
     end function greal
     !
     type(C_PTR) function stoi(x) bind(C,name='stoi')
       import C_PTR, C_LONG
       integer(kind=C_LONG), VALUE :: x
     end function stoi
     !
     type(C_PTR) function stor(x , prec) bind(C,name='stor')
       import C_PTR, C_LONG
       integer(kind=C_LONG), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function stor
     !
     type(C_PTR) function gpsi( x , prec) bind(C,name='gpsi')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function gpsi
     !
     integer(kind=C_LONG) function itos(x) bind(C,name='itos')
       import C_PTR, C_LONG
       type(C_PTR), VALUE :: x
     end function itos
     !
     type(C_PTR) function Pi2n(x, prec) bind(C,name='Pi2n')
       import C_PTR, C_LONG
       integer(kind=C_LONG), VALUE :: x
       integer(kind=C_LONG), VALUE :: prec
     end function Pi2n
     !
     type(C_PTR) function vecsum(x) bind(C,name='vecsum')
       import C_PTR, C_LONG
       type(C_PTR), DIMENSION(:)  :: x
       !  integer(kind=C_LONG), VALUE :: x
     end function vecsum
     !
      type(C_PTR) function gcos( x , prec) bind(C,name='gcos')
      import C_PTR, C_LONG
      type(C_PTR), VALUE :: x
      integer(kind=C_LONG), VALUE :: prec
    end function gcos
    !
    type(C_PTR) function gerepilecopy(av, z) bind(C,name='gerepilecopy')
      import C_PTR, C_LONG
      integer(kind=C_LONG), VALUE :: av
      type(C_PTR), VALUE :: z
    end function gerepilecopy
    !
    subroutine set_avma( av ) bind(C,name='set_avma')
      import C_LONG
      integer(kind=C_LONG), VALUE  :: av
    end subroutine set_avma
    !
    integer(kind=C_LONG) function get_avma( ) bind(C,name='get_avma')
      import C_LONG
    end function get_avma
  end interface
end module PARI


program compute_hardy_fast

  use mpi
  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI

  implicit none

  integer :: COMM_SIZE,COMM,IERR,RANK

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16)

  integer       :: ip,icj(0:50),MIT,totblock,irsp(7,2),i,mc,m1(50),pblocks,mmax,istartblock
  integer(dp1)  :: NC,jsum,RSN,L(0:50),KSET,K,KLP,KCUT(3:4)
  real(dp)      :: t,a,C,CE,p,zet(3:13),sp,xr(50),tim(0:350)
  real(dp)      :: et,ht,Y,F1,LC,zsum,tsum,rsphase,rszsum,t1,t2,rsfrac,ain,yphase,GAM,tpp,tpm
  real(dp)      :: phicoeff(0:15,50),fracL(0:50),see(15),xpbloc,fudge,GAMCOF(6)
  real(dp)      :: raecutoff,RN1,rae,raestartofblock,RN1startofblock
  real(dp)      :: rszsumtot,rszsum1,rstime
  complex(kind=16)   :: coeff(7,0:6),zlarcoeff(-1:6),zsmalcoeff(-1:6)
  complex(kind=16)   :: c1,c2,c3,c4,epi4,errcon,c5,c6,c7
  character(len=45)  :: YT
  integer(kind=8)    :: N_INTERNAL, N_MAX, NIS

  namelist /inputs/ et, YT


  COMMON/PARMS/ coeff,zlarcoeff,zsmalcoeff
  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS2/ GAM,ZET
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  COMMON/PARS5/ m1
  COMMON/PARS6/ t,a,ain,yphase,et,ht,Y,K,KLP,ip,mc
  COMMON/PARS7/ see
  COMMON/PARS8/ GAMCOF
  COMMON/PARS9/ YT
  
  !       First initialise a number of standard constants for use in the calculation of Z(t).
  !       HIGH PRECISION IS ABSOLUTELY ESSENTIAL FOR THE CALCULATIONS. The potential
  !       for round-off errors is ever present. Hence the necessity for real(33) precision values.
  !       The vagaries of the GNU compiler mean that these standard constants have to be set in weird 
  !       ways to achieve the Zero point 32 decimal accuracy.


  integer       :: iblock,IFACT,m
  integer(dp1)  :: M2,N1,ib,ial,aenums,MT,MTM,MTO
  integer(dp1)  :: mchange(3:4)
  real(dp)      :: t66(3:4),b,term,X,pcsum,e,t6,sca1,sca2
  real(dp)      :: rmc1(3:4),rmc2(3:4),rmc3(3:4),rmc4(3:4),rmc5(3:4),rmc6(3:4),zsum1,zsum2,RN4

  CALL MPI_Init(IERR)
  COMM = MPI_COMM_WORLD
  CALL  MPI_COMM_RANK(COMM, RANK, IERR)
  CALL  MPI_COMM_SIZE(COMM, COMM_SIZE, IERR)

  
!  Set p=pi, tpp=2*pi and tpm=-2*pi.  

  C=1.0
  p=4*ATAN(C)
  mmax=10
  tpp=2*p
  tpm=-tpp

!       Standard constants sqrt(pi) and exp(pi*i/4).

  C=1.0
  sp=sqrt(p)
  epi4=(1.0,1.0)/sqrt(2*C)

!       More standard constants relating to the zeta function and gamma function for real values.
!       Double precision is fine for the values of zeta(3-9).

  GAM=0.57721566490153286060D0

  zet(3)=1.20205690315959428540d0
  zet(5)=1.03692775514336992633d0
  zet(7)=1.00834927738192282684d0
  zet(9)=1.00200839282608221442d0
  zet(11)=1.0004941886041194646d0
  zet(13)=1.0001227133475784891d0


!     !! Coefficients used to calculate the complex error function. If you can make the
!     intrinsic error function erf work for complex arguments all the below is unnecessary

  C1=((1.0,0.0)*cos(1.0d0)+(0.0,1.0)*sin(1.0d0))*sqrt(2.0/p) 
  C2=((1.0,0.0)*cos(4.0d0)+(0.0,1.0)*sin(4.0d0))*sqrt(2.0/p)        
  C3=((1.0,0.0)*cos(2.25d0)+(0.0,1.0)*sin(2.25d0))*sqrt(2.0/p)
  C4=((1.0,0.0)*COS(6.25d0)+(0.0,1.0)*SIN(6.25d0))*sqrt(2.0/p)
  C5=((1.0,0.0)*cos(1.5625d0)+(0.0,1.0)*sin(1.5625d0))*sqrt(2.0/p)        
  C6=((1.0,0.0)*cos(3.0625d0)+(0.0,1.0)*sin(3.0625d0))*sqrt(2.0/p)
  C7=((1.0,0.0)*COS(5.0625d0)+(0.0,1.0)*SIN(5.0625d0))*sqrt(2.0/p)
  
!    The  error function coefficients start here.

!     COEFF(1,..) WILL APPLY FOR Z=1

  COEFF(1,0)=(0.969264211944215930381490d0,-0.474147636640994245161680d0)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1)

  COEFF(1,1)=C1*(1.0,-1.0)
  COEFF(1,2)=C1*(1.0,1.0)
  COEFF(1,3)=-C1*(1.0,-3.0)/3.0d0
  COEFF(1,4)=-C1*(5.0,-1.0)/6.0d0
  COEFF(1,5)=-C1*(11.0,13.0)/30.0d0
  COEFF(1,6)=C1*(9.0d0,-31.0d0)/90.0d0

!     NOW DO SAME FOR ERF(EXP(-i*PI/4)*2)

  COEFF(2,0)=(1.010311712025489491642626d0,0.273925759463539899021137d0)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2)

  COEFF(2,1)=C2*(1.0,-1.0)
  COEFF(2,2)=C2*(2.0,2.0)
  COEFF(2,3)=C2*(-7.0,9.0)/3.0d0
  COEFF(2,4)=-C2*(11.0,5.0)/3.0d0
  COEFF(2,5)=C2*(13.0,-109.0)/30.0d0
  COEFF(2,6)=C2*(129.0,-31.0)/45.0d0

  COEFF(3,0)=(1.338389640116239225341021d0,-0.096501782737190712909169d0)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.5)

  COEFF(3,1)=C3*(1.0,-1.0)
  COEFF(3,2)=C3*(1.5,1.5)
  COEFF(3,3)=C3*(-7.0,11.0)/6.0d0
  COEFF(3,4)=C3*(-15.0,-3.0)/8.0d0
  COEFF(3,5)=C3*(-13.0,-59.0)/40.0d0
  COEFF(3,6)=C3*(67.0,-53.0)/80.0d0

  COEFF(4,0)=(8.264692402664926098554228d-1,-1.394623184589031908854228d-1)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2.5)

  COEFF(4,1)=C4*(1.0,-1.0)
  COEFF(4,2)=C4*(2.5,2.5)
  COEFF(4,3)=C4*(-23.0,27.0)/6.0d0
  COEFF(4,4)=C4*(-155.0,-95.0)/24.0d0
  COEFF(4,5)=C4*(313.0,-913.0)/120.0d0
  COEFF(4,6)=C4*(1065.0,65.0)/144.0d0

  COEFF(5,0)=(1.215497305362876593269559d0,-0.344267547866857083642667d0)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.25)

  COEFF(5,1)=C5*(1.0,-1.0)
  COEFF(5,2)=C5*(1.25,1.25)
  COEFF(5,3)=C5*(-17.0,33.0)/24.0d0
  COEFF(5,4)=C5*(-245.0,-5.0)/192.0d0
  COEFF(5,5)=C5*(-767.0,-1633.0)/1920.0d0
  COEFF(5,6)=C5*(1665.0,-2335.0)/4608.0d0

  COEFF(6,0)=(1.260051027257371974124166d0,1.664741404414525369657245d-1)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.75)

  COEFF(6,1)=C6*(1.0,-1.0)
  COEFF(6,2)=C6*(1.75,1.75)
  COEFF(6,3)=C6*(-41.0,57.0)/24.0d0
  COEFF(6,4)=C6*(-511.0,-175.0)/192.0d0
  COEFF(6,5)=C6*(-143.0,-4561.0)/1920.0d0
  COEFF(6,6)=C6*(37527.0,-17353.0)/23040.0d0

  COEFF(7,0)=(7.873277503318070969692577d-1,1.234468235979886882660918d-1)

!     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2.25)

  COEFF(7,1)=C7*(1.0,-1.0)
  COEFF(7,2)=C7*(2.25,2.25)
  COEFF(7,3)=C7*(-73.0,89.0)/24.0d0
  COEFF(7,4)=C7*(-315.0,-171.0)/64.0d0
  COEFF(7,5)=C7*(827.0,-3419.0)/640.0d0
  COEFF(7,6)=C7*(12081.0,-879.0)/2560.0d0


!     NEXT COMPUTE THE COEFFICIENTS RELEVANT TO THE ASYMPTOTIC APPROXIMATION
!     OF ERF(Z*EXP(-i*PI/4)) WHEN ABS(Z)>>1.

  ZLARCOEFF(-1)=(1.0,1.0)/SQRT(2*p)
  ZLARCOEFF(0)=(1.0,0.0)
  DO I=1,6
     ZLARCOEFF(I)=(2*I-1)*ZLARCOEFF(I-1)/(0.0,2.0)
  ENDDO

!     NEXT COMPUTE THE COEFFICIENTS RELEVANT TO THE APPROXIMATION
!     OF ERF(Z*EXP(-i*PI/4)) FOR SMALL ABS(Z)<<1

  ZSMALCOEFF(-1)=(1.0,-1.0)/SQRT(0.5*p)
  ZSMALCOEFF(0)=(1.0,0.0)
  C1=ZSMALCOEFF(0)
  IFACT=1
  DO I=1,6
     IFACT=IFACT*I
     C1=C1*(0.0,1.0)
     ZSMALCOEFF(I)=C1/((2*I+1.0)*IFACT*1.0)
  ENDDO

!     !! All of the above code starting from !! can be excluded if one is able to access the intrinsic built in 
!     error function of the GNU fortran compiler to compute erf(Z) with Z=x*exp(+/-I*pi/4) and x>0 to at least double precision.

!     All ready.

!     Now this version of the code computes the first part of the Hardy function Z(t) using sequences of generalised quartic Gaussian sums added
!     together. This replaces the earlier version which computed the Hardy function Z(t) using a sequences of cubic Gaussian sums.
!     Theoretically the operational count necessary to compute Z(t) will scale as O((epsilon*t)^(1/5)). Here the value of mc denotes the
!     degree of the initial generalised Gaussian sums (mc=2 -quadratic, mc=3 - cubic, mc= 4 quartic etc.)  
!     For the mc=2 the quadratic case the operational count is around (et*t)^(1/3)*(log(t)^2.33), with mc=3 the cubic case
!     this goes down to O((epsilon*t)^(1/4))). Unfortunately the 'technical difficulties' mentioned above mean that the theoretical
!     operational count for the quartic case is unattainable and this code works much like the cubic algorithm.

!     The code requires two inputs t obviously, and epsilon=et in the code which specifies the error tolerance. These are read in from
!     a data file called "inputs15.nml". Obviously the source code must have access to this file.  

  open(3,file='inputs15.nml',status='old')
  read(unit = 3, nml=inputs)
  close(3,status='keep')


!     Output data pertaining to the relevant calculation is sent to a file called "zeta15v6resa". Again this file must exist and be accessible
!     by the source code.
   
  open(2,file='zeta15v6resa',status='old')

  write(2,4004) 'et= ',et
  write(2,4005) 't= ',YT
  READ(YT(1:1)  , '(I1)'  )  irsp(7,1)
  READ(YT(3:10)  , '(I8)'  )  irsp(1,1)
  READ(YT(11:18)  , '(I8)'  )  irsp(2,1)
  READ(YT(19:26)  , '(I8)'  )  irsp(3,1)
  READ(YT(27:34)  , '(I8)'  )  irsp(4,1)
  READ(YT(35:42)  , '(I8)'  )  irsp(5,1)
  READ(YT(44:45)  , '(I2)'  )  irsp(6,1)  
  
  CE=10
  C=0.0
  do i=5,1,-1
     C=C+irsp(i,1)/(CE**(8*i-1))
  enddo

  C=C/(CE**1)
  CE=1.0
  C=irsp(7,1)*CE+C
  CE=10.0
  t=(CE**(irsp(6,1)))*C
  write(2,4008) 'fortran t= ',t
  
!   The t value is stored as both a fortran variable and a C variable. The latter conversion is the function of the routine TGRABBER.

!   Now mc=4 which means that for the calculation Z(t) is expressed in terms of quartic Gauss sums.
  
  mc=4
  ip=4
  KCUT(3)=400
  KLP=40
  istartblock=1
  
!  For the first few blocks (block numbers 1 to 5 or 6) cubic Gauss sums must be used and these values of KCUT and KLP only apply in these
!  instances. It is recommended that these values are left as they are.

!  The code now sets values for KCUT (for block 7 onwards) and the total number of blocks denoted by totblock. There are two choices.
!  For inexperienced users just set KSET=0 and the code will generate suitable values of KCUT and totblock for you from the inputs of t
!  and epsilon. If you are more experienced then fix KSET something other than zero and you can input your own values for these parameters to
!  suit your needs. Varying KCUT and totblock produces minor changes in run times but typically nothing major.
  
  KSET=0
  IF (KSET.EQ.0) THEN
     KCUT(4)=MAX0(100,FLOOR(((log(0.0002)/log(et))**0.20)*(60*log(t)/log(10.0)-560.0),dp1))
     totblock=MAX0(25,FLOOR((((log(et)/log(0.0002))**0.20)*(3.545*log(t))-61.8)*0.75,dp1))
  ELSE
     KCUT(4)=850
     totblock=127
  ENDIF 
  K=KCUT(4)  
  C=100.0
  Y=111/C
  CALL pari_phases(yphase,rsphase,a)
  
! a=sqrt(8*t/pi) and is in practice just as important as t, since most of the parameters and calculations
! depend on sqrt(t), rather than t itself. It needs to be calculated to great accuracy.
  
  ain=1/a
  ht=log(et)/log(t)

  raecutoff=exp((log(et)+mc*log(t))/(mc+1.0d0))/(sp*log(t))
  pblocks=14
  xpbloc=(1-pblocks/log(t))/(mc+1.0)
  fudge=(et/(t**(1/(mc+2.0))))**xpbloc
  F1=1+fudge*(1-0.25/(Y*Y))/(2.0**(3*(mc-2)/(2*(mc+1.0d0))))
  LC=((mc-1)*log(t)/(2*mc+2.0d0)-log(log(t))+log((et**(1/(mc+1.0d0)))/(Y*sqrt(8.0d0))))/log(F1)

  write(2,3001) 'Current t value= ', t
  write(2,3001) 'Current a value= ', a
  write(2,3002) 'Relative error scale= ',et,'Cut off integer for QG sums= ',K
  write(2,3003) 'Max alpha cut off value= ',raecutoff,'Estimate of number of blocks needed= ',totblock
  write(2,3009) 'Riemann-Siegel phase= ',rsphase
  write(2,3009) 'pc phase            = ',yphase
  write(2,*) '******************************************************************************************'

!     errcon is a constant needed to complete the evauation of the additive term in routine q. It depends upon the 
!     set value of ip.

  errcon=(1.0,0.0)
  do i=1,ip
     errcon=errcon*(0.0,-1.0)
  enddo
  errcon=errcon*exp(-ip*1.0)*(ip**(ip*1.0))*(0.0,1.0)/(p**(ip+1.5))

!   Initialisation complete.
  
!   Begin the calculation of Z(t). The MSRSF is split into two parts;

!      ZP(t)+RSremainder(t).
  
!   Z(t) involves exactly sqrt(t/(2*pi)) summands. ZP(t) also involves O(sqrt(t)) summands (slightly fewer than Z(t)),
!   whilst RSremainder(t) involves only NC=O(t^(1/4)) summands.
  
!   ZP(t) is calculated in the variable zsum. This is, by far, the largest part of the MSRSF and it is
!   this term that is expressed in terms of quartic Gaussian sums. It is calculated in a
!   series of calculation blocks, numbered zero,one,two up to totblock.

!   The final calculation of involves finding RSremainder(t). Since the number of summands is tiny compared to both Z(t)
!   and ZP(t) it can be evaluated directly, term by term, in the form it appears in MSRSF. It amounts, in effect, to one final calculation block.
!   RSremainder(t) is calculated in the variable rszsum.


  call cpu_time(t1)

  !   The first calculation is that of Block zero. This block is special, because unlike the other blocks the terms are
  !   not amalgamated together to form Gauss sums but are calculated directly one by one. In practice each one of these terms
  !   corresponds to the result of adding O(t^(1/4)) summands making up Z(t) together, a colossal efficiency saving.
  
  C=1.0
  t6=((sqrt(p/(8*C))*a)**(1/(3*C)))
  call start(a,t6,yphase,zsum,M2)
  write(2,2007) 'Block 0. Initial alpha value= ',M2
  ib=2*floor(t6**(1.2d0),dp1)  
  if (mod(ib,2).gt.0.5) then
     ib=ib+1
  endif
  N1=M2+ib
  RN1=real(N1,dp)
  c=sqrt(8/a)        
  do ial=M2,N1,2
     call ter(ial,ain,yphase,t,term)
     zsum=zsum+term
  enddo
  zsum=c*zsum
  call cpu_time(t2)
  tim(0)=t2-t1
  write(2,2003) 'zsum after this block= ',zsum
  write(2,2012) 'alpha value at end of this block= ',N1,'Cpu time needed for this block= ',tim(0)
  write(2,*) '----------------------------------------------------------------------'

! Calculation of block zero complete. Now move on to the Gauss sum blocks 1,2,...,totblock.
! In this code the early blocks involve cubic m=3 sums. Then quickly followed by quartic m=4 sums.
! Because two orders of Gauss sums are utilised, two sets of sum parameters must be initialised.
  
  X=exp(((1/(mc+2.0d0)-ht)/(mc+1.0d0)))

  do m=3,mc
   rmc1(m)=1/(m+1.0)
   rmc2(m)=(m-1.0)*rmc1(m)
   rmc3(m)=(2*m-1)*rmc1(m)/2.0
   rmc4(m)=(3*(m-3)/2.0)*rmc1(m)
   rmc5(m)=1/(m-1.0)
   rmc6(m)=((4*m-5.0)/2.0)*rmc1(m)
   t66(m)=(t**(0.5*rmc2(m)))
  enddo
  mchange(3)=M2+2*floor((t6**(1.5d0)),dp1)
  m=3
  sca1=1.0
  sca2=1.0
  zsum1=0.0
  zsum2=0.0

! If a restart run is necessary activate code below. This facility allows one to restart a calculation
! that stopped prematurely for some reason, without having to start all over again from Block zero.
  
! zsum=4913.412791572dp
! write(2,2003) 'zsum after this block= ',zsum

! The value zsum reached at the last completed block calculation. Listed in zeta15v6resa or the output file used.

! RN1=123456789.0dp
! write(2,2018) 'Alpha value restarts at= ',RN1
! The number that alpha has reached at the end of the last completed block, an odd integer but input as a odd real.

! aenums=7219310dp1
! The number of sums in the last completed block.

! istartblock=120  
! The number of the next block to be calculated.
  
! Finally if a restart is being set up and aenums has reached a constant value (usual case for restarts) then comment out
! the four code lines below, beginning  "MTO=MT" and "if (MT.gt.(1.5*MTO)) then MT=floor(1.5*MTO,dp1) endif". Obviously,
! PUT THIS CODE BACK IN after the restart has finished. N.B. if aenums has not reached a constant, usually around block 30
! or so, then one is better off doing another complete rerun from Block zero.  
  
  CALL pari_init(10000000_8, 2_8)  

! Main loop, incorporating the Gauss sum blocks
  
  do iblock=istartblock,totblock          

     if(RANK .eq. 0) then

        write(6,*) 'iblock= ',iblock

     endif

     call cpu_time(t1)
     if (RN1.lt.(Y*a)) then
        aenums=floor((X**iblock)*ib/2.0,dp1)
        e=RN1/a
        MT=floor(sca1*((((et/p)**rmc5(m))*a)**rmc2(m))*((e-1.0)**rmc3(m))/(2.0**rmc4(m)),dp1)         
     else
        MTO=MT
        b=a/(2.0*RN1)
        MT=floor(sca2*(et**rmc1(m))*t66(m)*(1/b-b)/(sp*(2.0**rmc6(m))),dp1)
        if (MT.gt.(1.5*MTO)) then
         MT=floor(1.5*MTO,dp1)
        endif
     endif
     if (mod(MT,2).lt.0.25) then
        MT=MT-1
     endif

     write(2,2002) 'Block number= ',iblock,' Number of Gsums= ',aenums,'Length of Gsums= ',MT,'Type of Gsums= ',m   
     MTM=(MT-1)/2
     rae=RN1+2.0+real(MT,dp)

!    MTM is the length of (number of summands) of the Gauss sums calculated in Block No. iblock.
!    rae is the first pivot integer (even integer but stored as a fortran real, since it will exceed the largest size of a fortran long integer)
!    of the first GS of Block No. iblock. aenums is the number of individual GS (all of length MTM) making up Block No. iblock.
     
     if (m.eq.3) then
      K=KCUT(3)
     else
      KCUT(3)=MTM/2
      K=min0(KCUT(3),KCUT(4))
     endif


     raestartofblock=RN1+2.0+real(MT,dp)
     RN1startofblock=RN1

     zsum1=0.0      

!    Parallelise the calculation of all the individual GS making up Block No. iblock. 

     do jsum=(rank+1),aenums,COMM_SIZE  


        RN4=2.0*(real(MT,dp)+1)*(real(jsum-1,dp))
        rae=raestartofblock  + RN4            
        RN1=RN1startofblock  + RN4

!     rae is the pivot integer of each GS. RN1 marks the boundary between individual GS within the Block itself.
!     At the end of a Block calculation, RN1 forms the boundary point between the end of Block number iblock and the start of
!     Block number iblock+1.

        call alphasum(MTM,rae,pcsum,m)   
        zsum1=zsum1+pcsum


     enddo

     RN1=RN1startofblock + (2*(real(MT,dp)+1.0)*real(aenums,dp))


     !CALL MPI_ALLREDUCE(zsum1,zsum2,1,MPI_REAL16,MPI_SUM,COMM,IERR)
     ! The global sum will only make sense in process 0
     CALL MPI_REDUCE(zsum1,zsum2,1,MPI_REAL16,MPI_SUM,0,COMM,IERR)


     zsum=zsum+zsum2

     call cpu_time(t2)
     tim(iblock)=t2-t1

     if(RANK.EQ.0) then
          write(2,2003) 'zsum after this block= ',zsum
     endif
        write(2,2008) 'alpha value at end of this block= ',RN1,'Cpu time needed for this block= ',tim(iblock)
        write(2,*) '----------------------------------------------------------------------'

!   Calculation of Block No. iblock complete. Now check if the overall calculation has progressed sufficiently to
!   change from smaller cubic GS to larger quartic GS. This usually occurs very quickly certainly by Block 5 or so.

    if (m.lt.mc) then
     if (RN1.gt.real(mchange(m),dp)) then
       m=m+1
       sca1=0.70
       sca2=0.40
      endif
     endif
  enddo


  CALL pari_close()  
  
  raecutoff=RN1
  tsum=tim(0)
  do i=istartblock,totblock
     tsum=tsum+tim(i)
  enddo
  if(RANK.EQ.0) then
     write(2,3004) 'Final value of ZP=', zsum
  endif
  write(2,3005) 'Total cpu time needed for this calculation= ',tsum

! All the contributions from blocks numbered from 1 to totblock are now done, completing the calculation of ZP(t).
! Now the overall calculation of Z(t) is finalised by computing the very much smaller RSremainder(t) term.

  call cpu_time(t1)
  rszsum=0.0
  rszsum1=0.0
  rszsumtot=0.0
  CE=4.0
  C=1.0
  CE=raecutoff*(1-sqrt(1-(a/(raecutoff*C))**2))/CE
  NC=aint(CE,dp1)
  write(2,3006) 'NC cut off integer of RS part of sum= ',NC

  CE=1.0

     
  NIS = 0
  N_INTERNAL = 100000
  N_MAX = floor( NC / (N_INTERNAL*1.0) )

  DO jsum=(RANK+1),N_MAX,COMM_SIZE

     CALL pari_calc(rszsum,jsum,rsphase,N_INTERNAL,NIS)

  ENDDO

  rszsum1 = 0.0

  !CALL MPI_ALLREDUCE(RSZSUM,RSZSUM1,1,MPI_REAL16,MPI_SUM,COMM,IERR)
  ! Global sum will only be available at process 0.
  CALL MPI_REDUCE(RSZSUM,RSZSUM1,1,MPI_REAL16,MPI_SUM,0,COMM,IERR)

  rszsumtot=rszsumtot+rszsum1

! Additional bit of code to do the final bit of calculation N_MAX*N_INTERNAL+1 to NC.
     
  rszsum1 = 0.0
  NIS=N_MAX*N_INTERNAL+1
  N_INTERNAL=NC
  jsum=1
  CALL pari_calc(rszsum1,jsum,rsphase,N_INTERNAL,NIS)

  rszsumtot=rszsumtot+rszsum1
     

  RSZSUM=RSZSUMTOT

  CE=sqrt(t/(2*p))
  RSN=aint(CE,dp1)
  rsfrac=CE-RSN
  CE=sqrt(1/CE)
  C=cos(2*p*(rsfrac*(rsfrac-1.0)-0.0625))/cos(2*p*rsfrac)
  rszsum=2*rszsum+((-1)**(RSN-1))*CE*C

!   Rszsum is the remainder contribution to MSRSF of Z(t), plus an additional O(t^(-0.25)) term. This extra term was first formulated by Riemann and
!   subsequently re-discovered by Siegel, 73 years later, amongst Riemann's unpublished papers stored in the Gottingen University Library.

  call cpu_time(t2)
  rstime=t2-t1
  if(RANK.EQ.0) then
     write(2,3004) 'Final value of RS contribution to sum= ', rszsum  
  endif                                        
  write(2,3005) 'Total cpu time needed for this calculation= ',rstime
  write(2,3010) 'Riemann-Siegel NT value= ',RSN
  write(2,3009) 'Riemann-Sieg NT frac= ',rsfrac
  write(2,3007) 'Value of RS contribution to sum= ', rszsum                                          
  write(2,3005) 'Total cpu time needed for this calculation= ',(t2-t1) 
  write(2,3010) 'Riemann-Siegel NT value= ',RSN
  write(2,3009) 'Riemann-Sieg NT frac= ',rsfrac
  tsum=tsum+rstime
  zsum=zsum+rszsum

  !     Estimate for Z(t)=ZP(t)+RSremainder(t)

  write(2,3008) 'Grand total of Hardy function Z(t)= ',zsum                                          
  write(2,3005) 'Total cpu time needed for this calculation= ',tsum
  close(2,status='keep')

  !     All Done.

  STOP

2002 FORMAT(A14,1X,I4,A18,1X,I14,1X,A17,1X,I14,1X,A15,I1)
2003 FORMAT(A23,2X,F15.9)
2007 FORMAT(A31,1X,I17)
2008 FORMAT(A34,2X,F24.1,1X,A32,2X,E13.6)
2012 FORMAT(A34,2X,I19,1X,A32,2X,E13.6)
2018 FORMAT(A24,1X,F24.1)  


3001 FORMAT(A17,E41.34)
3002 FORMAT(A23,F14.8,2X,A28,I4)
3003 FORMAT(A25,E15.21,2X,A37,I3)
3004 FORMAT(A18,1X,F17.10)
3005 FORMAT(A43,1X,E13.6)
3006 FORMAT(A38,1X,I21)
3007 FORMAT(A34,1X,F17.10)
3008 FORMAT(A37,1X,F17.10)
3009 FORMAT(A22,1X,F29.23)
3010 FORMAT(A25,1X,I21)

4004 FORMAT(A4,F7.5)
4005 FORMAT(A3,A45)
4008 FORMAT(A12,E47.40)
  
  call MPI_Finalize(IERR)

end program compute_hardy_fast





SUBROUTINE ERF(Z,N,ER)

  !     COMPUTES ERF(EXP(-i*PI/4)*Z) FOR Z REAL.

  !     IT USES JUST N TERMS ONLY.
  !     USES TAYLOR SERIES IF ABS(Z) NEAR 1, POWER SERIES IF ABS(Z)<<1 AND
  !     ASYMPTOTIC SERIES IF ABS(Z)>>1. THE COEFFICIENTS OF THESE SERIES ARE
  !     STORED IN CEF,ZSMAL AND ZLAR RESPECTIVELY. SIX TAYLOR SERIES ARE USED
  !     CENTRED ON Z=1,1.25,1.5,1.75,2,2.25 and 2.5. BELOW ABS(Z)<0.8 THE POWER SERIES
  !     ABOUT Z=0 IS EMPLOYED AND ABOVE ABS(Z)>2.625 THE ASYMPTOTIC SERIES IS USED

  implicit none

  integer, parameter :: dp = selected_real_kind(33) 

  integer       :: N,I
  real(dp)      :: Z,P,E,ZA,SP,TPP,TPM
  complex(kind=16)   :: ER,CEF(7,0:6),ZLAR(-1:6),ZSMAL(-1:6),V,EPI4

  COMMON/PARMS/ CEF,ZLAR,ZSMAL
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4

  IF (N.GT.6) THEN
     WRITE(6,*) 'INCREASE NUMBER OF COEFFICIENTS IN SUBROUTINE ERF.'
     STOP
  ENDIF

  !     DO ALL MAIN CALCULATIONS ASSUMING Z>0.

  ZA=ABS(Z)
  IF (ZA.LT.1E-14) THEN
     ER=ZSMAL(-1)*ZA
     GOTO 1000
  ENDIF
  ER=(0.0,0.0)
  IF (ZA.GE.0.8.AND.ZA.LT.1.125) THEN

     !     ZA CLOSE TO 1, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1)
     !     TO GET ACCURATE SOLUTION

     E=ZA-1.0
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(1,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(1,0)
  ELSE IF (ZA.GE.1.125.AND.ZA.LT.1.375) THEN

     !     ZA CLOSE TO 1.25, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.25)
     !     TO GET ACCURATE SOLUTION

     E=ZA-1.25
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(5,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(5,0)
  ELSE IF (ZA.GE.1.375.AND.ZA.LT.1.625) THEN

     !     ZA CLOSE TO 1.5, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.5)
     !     TO GET ACCURATE SOLUTION


     E=ZA-1.5
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(3,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(3,0)
  ELSE IF (ZA.GE.1.625.AND.ZA.LT.1.875) THEN

     !     ZA CLOSE TO 1.75, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.75)
     !     TO GET ACCURATE SOLUTION


     E=ZA-1.75
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(6,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(6,0)
  ELSE IF (ZA.GE.1.875.AND.ZA.LT.2.125) THEN

     !     ZA CLOSE TO 2.0, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.0
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(2,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(2,0) 
  ELSE IF (ZA.GE.2.125.AND.ZA.LT.2.375) THEN

     !     ZA CLOSE TO 2.25, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2.25)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.25
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(7,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(7,0)

  ELSE IF (ZA.GE.2.375.AND.ZA.LT.2.625) THEN

     !     ZA CLOSE TO 2.5, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2.5)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.5
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(4,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(4,0)     
  ELSE IF (ZA.GE.2.625) THEN

     !    ZA LARGE, USE ASMPTOTIC APPROXIMATION. ACTUALLY WE COMPUTE
     !    ERFC(EXP(-i*PI/4)*Z) FIRST AND THEN USE ERF(Z)=1-ERFC(Z).
     !    THE APPROXIMATION IS VALID WHEN ABS(ARG(EXP(-i*PI/4)*Z))=PI/4,
     !    WHICH APPLIES WHEN Z>0, IS LESS THAN 3*PI/4. THE CASE Z<0 IS DEALT
     !    WITH AT THE END.

     E=ZA*ZA
     V=((1.0,0.0)*COS(E)+(0.0,1.0)*SIN(E))/ZA
     DO I=0,N
        ER=ER+ZLAR(I)/(E**I)
     ENDDO

     !    CONVERT FROM ERFC TO ERF

     ER=(1.0D0,0.0D0)-ZLAR(-1)*V*ER
  ELSE

     !    ZA SMALL<1, USE THE SERIES SOLUTION WHICH CONVERGES RAPIDLY.

     E=ZA*ZA
     DO I=0,N
        ER=ER+ZSMAL(I)*(E**I)
     ENDDO
     ER=ZA*ZSMAL(-1)*ER
  ENDIF

  !     NOW IF Z<0 ORIGINALLY USE ERF(-Z)=-ERF(Z)

1000 IF (Z.LT.0.0D0) THEN
     ER=-ER
  ENDIF
  RETURN
END SUBROUTINE ERF


SUBROUTINE LGAM(GA,XX)

  !     COMPUTES THE LOGARITHM OF THE GAMMA FUNCTION FOR REAL GA>0.
  !     IN THIS APPLICATION ONLY CALLED WHEN GA>12 AND FACTORIAL FUNCTION GETS LARGE VERY FAST.

  implicit none

  integer, parameter :: dp = selected_real_kind(33) 

  integer       :: I
  real(dp)      :: GA,XX,P,SP,TPP,TPM,GAMCOF(6),XXX,Y,TMP,SUM
  complex(kind=16)  :: EPI4

  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  COMMON/PARS8/ GAMCOF

  IF (ABS(GA).LT.1D-14) THEN
     XX=0.0
     RETURN
  ENDIF

  XXX=GA
  Y=XXX
  TMP=XXX+5.5D0
  TMP=(XXX+0.5D0)*LOG(TMP)-TMP
  SUM=1.000000000190015D0
  DO I=1,6
     Y=Y+1.0D0
     SUM=SUM+GAMCOF(I)/Y
  ENDDO

  XX=TMP+LOG(SQRT(TPP)*SUM/XXX)

  RETURN
END SUBROUTINE LGAM


! Subroutine ter calculates the individual terms of Block No. Zero.

subroutine ter(k,ain,yphase,t,term)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)  
  integer, parameter :: dp1 = selected_int_kind(16)

  integer(dp1)  :: k
  real(dp)      :: ain,yphase,t,term,s,v,b

  s=k*ain
  b=sqrt(s*s-1)
  v=t*(b*(b-s)+log(s+b))+yphase
  term=cos(v)/sqrt(b)
  return
END subroutine ter

! Subroutine start sets the starting point of the calculation (the beginning of Block No. zero)
! and occasionally the transition term if a=sqrt(8*t/pi) is within O(t^(-1/6)) of an odd integer.

subroutine start(a,t6,yphase,zsum,M2)

  implicit none

  integer, parameter :: dp = selected_real_kind(33) 
  integer, parameter :: dp1 = selected_int_kind(16) 

  integer       :: itra,i
  integer(dp1)  :: M2
  real(dp)      :: a,t6,yphase,zsum,g,gm,p,sp,tpp,tpm,bot,top,t3,C,cc,AA,wb
  real(dp)      :: abb(64),wee(64),ri
  complex(kind=16) :: EPI4,sum,BB,expon,func,AAC
       
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  
  itra=0
  if (mod(floor(a,dp1)+2,2).eq.0) then
     M2=floor(a,dp1)+1
  else
     M2=floor(a,dp1)+2
  endif
  g=(M2-a)*t6
  gm=(a-(M2-2))*t6
  zsum=0.0d0
  if (g.lt.3.2) then
     itra=1
     M2=M2+2
     C=1.0
     sum=(0.0,0.0)
     t3=t6*t6
     cc=tpp**0.25
     bot=-2*cc*sqrt(g)/t3
     AA=cc*(sqrt(2*g)*t3*t3*0.25-t3*g*cc)*0.5
     top=4/sqrt(AA)
     wb=t3*((0.25*t3-1.5*cc*sqrt(2*g))*t3+6.375*cc*cc*g)/(6*C)
     BB=EPI4*wb
        
     call gauleg(bot,top,abb,wee)

     do i=1,64
      expon=-(abb(i)**2)*(AA+abb(i)*BB)
      func=exp(expon)
      sum=sum+wee(i)*func
     enddo
      
     write(2,6001) '   Initial int = ',real(sum),aimag(sum)
        
     wb=yphase+0.5*sqrt(tpp)*g*t3-sqrt(p)*cc*(g**1.5)/(3*C)
     expon=(0.0,1.0)*wb
     ri=2*real(exp(expon)*sum)
     zsum=(M2-2)*ri/(sqrt(2*a)*(2**(1.5)))
     write(2,6000) '  Initial zsum = ',zsum
  endif
  if (gm.lt.3.2) then
     itra=-1
     M2=M2
     C=1.0
     sum=(0.0,0.0)
     t3=t6*t6
     cc=tpp**0.25
     bot=0.0
     AA=0.5*(sqrt(tpp)*gm*t3)
     AAC=AA+(0.0,1.0)*0.5*(cc*sqrt(2*gm)*t3*t3*0.25-1.75*(cc**3)*(gm**1.5)/sqrt(2*C))
     BB=(0.0,1.0)*1.5*cc*sqrt(2*gm)*t3*t3/(6*C)
     BB=(t3*(0.25*t3*t3-6.375*cc*cc*gm))/(6*C)-BB
     BB=EPI4*BB
     top=3/(real(BB)**(1/(3*C)))
                
     call gauleg(bot,top,abb,wee)

     do i=1,64
      expon=-(abb(i)**2)*(AAC+abb(i)*BB)
      func=exp(expon)
      sum=sum+wee(i)*func
     enddo
     write(2,6001) '   Initial int = ',real(sum),aimag(sum)

     wb=yphase-0.5*sqrt(tpp)*gm*t3
     expon=(0.0,1.0)*wb
     ri=2*real(exp(expon)*sum*exp(-sqrt(p)*cc*(gm**1.5)/(3*C)))
     zsum=(M2-2)*ri/(sqrt(2*a)*(2**(1.5)))
     write(2,6000) '  Initial zsum = ',zsum      
  endif
  return

6000   FORMAT(A17,2X,E38.31)       
6001   FORMAT(A17,2X,E38.31,2X,E38.31)  
END subroutine start

! alphasum parameterises the various quartic sums based on the pivot integer rae.

subroutine alphasum(MTM,rae,pcsum,m)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)  
  integer, parameter :: dp1 = selected_int_kind(16) 

  integer       :: m,ip,mc,i
  integer(dp1)  :: MTM,K,KLP,jj
  real(dp)      :: rae,pcsum,p,sp,t,a,et,e,s,amp,s1,pc,xx,wm,wp,initialcoeff(4)
  real(dp)      :: rs,ht,ain,yphase,Y,pc2,con2,tpp,tpm,initialcoeffcc(4)
  real(dp)      :: cc,cc1,tse,Ppc(3:4),con3
  complex(kind=16)   :: c1,c2,pcgsum(1:2),EPI4,pha

  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  COMMON/PARS6/ t,a,ain,yphase,et,ht,Y,K,KLP,ip,mc


6927 e=(rae*ain)**2
  s=sqrt(e-1.0)
  pc=2*(e+sqrt(e)*s)-1
  pc2=(pc*pc-1)
  s1=a/sqrt(pc)
  amp=sqrt(8/(a*s))
  xx=1/(pc-1)
  rs=p*pc/(2*pc2)
  con2=((sqrt(pc)*xx)**3)/(a)
  Ppc(4)=pc+1.0
  con3=pc*pc*(xx**5)*2*Ppc(4)/(a*a)

  ! first zeroth order phases, just I*wp and I*wm will do.

  CALL pari_tse(rae,tse)

  wp=-tse-p*(1-2*(xx+s1))/8.0
  wm=wp-p*s1/2.0

  if (m.eq.3) then

  ! find parameters for a cubic sum
     
   wp=wp+p*con2/3.0
   wm=wm-p*con2/3.0

  ! parameters phi1

   initialcoeff(1)=xx/2.0+s1/4.0+con2
   initialcoeffcc(1)=initialcoeff(1)-s1/2.0-2*con2

  !  parameters phi2

   initialcoeff(2)=xx/2.0+2*con2
   initialcoeffcc(2)=initialcoeff(2)-4*con2

  ! parametes phi3
   
   initialcoeff(3)=4*con2/3.0
   initialcoeffcc(3)=-initialcoeff(3)

  else if (m.eq.4) then

  ! parameters for a quartic sum
     
   wp=wp+p*con2/3.0+p*con3/4.0
   wm=wm-p*con2/3.0+p*con3/4.0

! parameters phi1

   initialcoeff(1)=xx/2.0+s1/4.0+con2+con3
   initialcoeffcc(1)=initialcoeff(1)-s1/2.0-2*con2

! parameters phi2

   initialcoeff(2)=xx/2.0+2*con2+3*con3
   initialcoeffcc(2)=initialcoeff(2)-4*con2

! parameters phi3

   initialcoeff(3)=4*(con2+3*con3)/3.0
   initialcoeffcc(3)=initialcoeff(3)-8*con2/3.0

! parameters phi4

   initialcoeff(4)=2*con3
   initialcoeffcc(4)=initialcoeff(4)
  endif

  initialcoeff(1)=initialcoeff(1)-anint(initialcoeff(1),dp)
  initialcoeffcc(1)=initialcoeffcc(1)-anint(initialcoeffcc(1),dp)


  !      If sum length is short - i.e. MTM<KCUT, then compute each sum directly.

  if (MTM.lt.K) then
     pcgsum(1)=(1.0,0.0)
     pcgsum(2)=pcgsum(1)
     do jj=1,MTM
        cc=jj*initialcoeff(m)
        cc1=jj*initialcoeffcc(m)
        do i=m-1,1,-1
           cc=(cc+initialcoeff(i))*jj
           cc1=(cc1+initialcoeffcc(i))*jj
        enddo
        pha=(0.0,1.0)*tpp*cc
        pcgsum(1)=pcgsum(1)+exp(pha)
        pha=(0.0,1.0)*tpp*cc1
        pcgsum(2)=pcgsum(2)+exp(pha)
     enddo
  else

  ! compute the sums by finding the appropriate hiearchical chain, compute the kernel sum and move back up the chain
     
     call computarrs(m,initialcoeff,MTM,K,KLP)
     call mgausssum(m,ip,pcgsum(1))
     call computarrs(m,initialcoeffcc,MTM,K,KLP)
     call mgausssum(m,ip,pcgsum(2))

  endif

! Multiply each sum by its constant phase, add, take real parts and multiply by amplitude.

  c1=(0.0,1.0)*wp
  c1=exp(c1)
  c2=(0.0,1.0)*wm
  c2=exp(c2)
  pcsum=amp*(real(c1*pcgsum(1)+c2*pcgsum(2)))


  if (abs(pcsum).gt.0.5) then
   if (KLP.eq.40) then
    KLP=2*KLP
    goto 6927
   else if (KLP.eq.80) then
    KLP=2*KLP
    goto 6927
   else
    write(2,*) 'error in the sum set by rae= ',rae
    write(2,*) pcgsum(1)
    write(2,*) pcgsum(2)
    stop
   endif
  endif
  if (KLP.gt.60) then
    KLP=40
   endif

  return
END subroutine alphasum



!     THE NEXT ROUTINE COMPUTES THE PSI FUNCTION FOR REAL X

SUBROUTINE PSI(X,PS)

  IMPLICIT NONE

  integer, parameter :: dp = selected_real_kind(33) 

  real(dp) X,PS,P,GAM,ZET(3:13),XX,Z,SUM,CA,Q,SP,TPP,TPM
  complex(kind=16)   :: EPI4
  INTEGER  :: I

  COMMON/PARS2/ GAM,ZET
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4

  XX=ABS(X)

  IF (XX.LT.1e-12) THEN
     WRITE(6,*) 'PSI FUNCTION UNDEFINED AT X=0'
     STOP
  ENDIF

  IF (X.LT.0.0.AND.ABS(X-NINT(X)).LT.1e-12) THEN
     WRITE(6,*) 'PSI FUNCTION UNDEFINED AT NEGATIVE INTEGERS'
     STOP
  ENDIF

  !      SOME AWKWARD SPECIAL VALUES PSI(1)=-GAM, PSI(2), PSI(3), PSI(0.5).

  IF (ABS(XX-1.0).LT.1e-12) THEN
     PS=-GAM
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  IF (ABS(XX-2.0).LT.1e-12) THEN
     PS=0.422784335098467d0
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  IF (ABS(XX-3.0).LT.1e-12) THEN
     PS=0.922784335098467d0
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  !      PSI(1/2)=-GAM-2*LOG(2)

  IF (ABS(XX-0.5).LT.1e-12) THEN
     PS=-GAM-2*LOG(2.0)
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  !      WORK WITH POSITIVE X VALUES VIZ. XX=ABS(X) FIRST. WORK OUT PSI(XX)

  IF (XX.LT.4.0) THEN

     !   magnitude of xx<4 use series solution valid for small xx. For positive xx=abs(x)

     IF (XX.LT.0.5) THEN
        Z=XX
        Q=-1/Z
     ELSE IF (XX.GT.0.5.AND.XX.LE.1.5) THEN
        Z=XX-1.0
        Q=0.0
     ELSE IF (XX.GT.1.5.AND.XX.LE.2.5) THEN
        Z=XX-2.0
        Q=1/(Z+1.0)
     ELSE IF (XX.GT.2.5.AND.XX.LE.3.5) THEN
        Z=XX-3.0   
        Q=1/(Z+1.0)+1/(Z+2.0)
     ELSE
        Z=XX-4.0
        Q=1/(Z+1.0)+1/(Z+2.0)+1/(Z+3.0)
     ENDIF

     !   series for psi(1+z) equation 5.7.5. psi(2+z)=psi(1+z)+1/(1+z); psi(3+z)=psi(1+z)+1/(1+z)+1/(2+z); with z an element [-0.5,0.5].

     SUM=0.0
     DO I=1,6
        SUM=SUM+(ZET(2*I+1)-1.0)*(Z**(2*I))
     ENDDO
     IF (ABS((Z-0.5)-NINT(Z-0.5)).LT.1e-12) THEN
        CA=0.0
     ELSE
        CA=0.5*P/TAN(P*Z)
     ENDIF
     PS=1/(2*Z)-CA+1/(Z*Z-1.0)+1-GAM-SUM+Q
  ELSE

     ! magnitude of xx>4 asymptotic series converges fast

     PS=LOG(XX)-0.5/XX-1/(12.0*XX*XX)+1/(120.0*(XX**4))-1/(252.0*(XX**6))+1/(240.0*(XX**8))
  ENDIF

  !  if x<0 then use reflection formula

4020 IF (X.LT.0.0)THEN
     IF (ABS((XX-0.5)-NINT(XX-0.5)).LT.1e-12) THEN
        CA=0.0
     ELSE
        CA=P/TAN(P*XX)
     ENDIF
     PS=PS+1/XX+CA
  ENDIF
  RETURN
END SUBROUTINE PSI

! See documentation for Scheme

subroutine scheme(m,j,small)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16)    

  integer       :: m,j,icj(0:50),MIT,i,ki,is,mmax
  integer(dp1)  :: L(0:50)
  real(dp)      :: small,phicoeff(0:15,50),xr(50),fracL(0:50),see(15),y
  real(dp)      :: rfact1,con1,rfact2,con2,rfact3,con3,con4,p(18),q(18),r(18)

  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS7/ see

  y=2*phicoeff(2,j)
  see(1)=1/(2*y)
  see(2)=-phicoeff(3,j)/(y**3)
  ki=1
  rfact1=0.5
  rfact2=1.0
  rfact3=2.0
  do i=3,m
     if (i.eq.3) then
        small=9*(phicoeff(3,j)**2)/(2*(y**5))
        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.4) then
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**7)
        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.5) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**9)

        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**7)   

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.6) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**11)

        rfact2=rfact2*6*(i-2)*(2*i-5)/((i+1)*1.0*(i-4))
        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**9)

        con3=(phicoeff(3,j)**(i-6))*(9*phicoeff(3,j)*phicoeff(6,j)+10*(i-5)*phicoeff(4,j)*phicoeff(5,j) &
             +(i-6)*(i-5)*16*(phicoeff(4,j)**3)/(9*phicoeff(3,j)))
        small=small-rfact3*con3/(y**8)

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.7) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**13)

        rfact2=rfact2*6*(i-2)*(2*i-5)/((i+1)*1.0*(i-4))
        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**11)

        rfact3=rfact3*6*(i-3)*(2*i-5)/((i+1)*1.0*(i-4))
        con3=(phicoeff(3,j)**(i-6))*(9*phicoeff(3,j)*phicoeff(6,j)+10*(i-5)*phicoeff(4,j)*phicoeff(5,j) &
             +(i-6)*(i-5)*16*(phicoeff(4,j)**3)/(9*phicoeff(3,j)))
        small=small-rfact3*con3/(y**10)

        small=small+(42*phicoeff(3,j)*phicoeff(7,j)+48*phicoeff(4,j)*phicoeff(6,j)+25*(phicoeff(5,j)**2)) &
             /(2*(y**9))

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.8) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**15)

        rfact2=rfact2*6*(i-2)*(2*i-5)/((i+1)*1.0*(i-4))
        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**13)

        rfact3=rfact3*6*(i-3)*(2*i-5)/((i+1)*1.0*(i-4))
        con3=(phicoeff(3,j)**(i-6))*(9*phicoeff(3,j)*phicoeff(6,j)+10*(i-5)*phicoeff(4,j)*phicoeff(5,j) &
             +(i-6)*(i-5)*16*(phicoeff(4,j)**3)/(9*phicoeff(3,j)))
        small=small-rfact3*con3/(y**12)

        small=small+5*(63*(phicoeff(3,j)**2)*phicoeff(7,j)+144*phicoeff(3,j)*phicoeff(4,j)*phicoeff(6,j) &
             +75*phicoeff(3,j)*(phicoeff(5,j)**2)+80*(phicoeff(4,j)**2)*phicoeff(5,j))/(y**11)

        small=small-2*(12*phicoeff(3,j)*phicoeff(8,j)+14*phicoeff(4,j)*phicoeff(7,j)+15*phicoeff(5,j)*phicoeff(6,j))/(y**10)

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.9) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**17)

        rfact2=rfact2*6*(i-2)*(2*i-5)/((i+1)*1.0*(i-4))
        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**15)

        rfact3=rfact3*6*(i-3)*(2*i-5)/((i+1)*1.0*(i-4))
        con3=(phicoeff(3,j)**(i-6))*(9*phicoeff(3,j)*phicoeff(6,j)+10*(i-5)*phicoeff(4,j)*phicoeff(5,j) &
             +(i-6)*(i-5)*16*(phicoeff(4,j)**3)/(9*phicoeff(3,j)))
        small=small-rfact3*con3/(y**14)

        small=small+11*(378*(phicoeff(3,j)**3)*phicoeff(7,j)+1296*(phicoeff(3,j)**2)*phicoeff(4,j)*phicoeff(6,j)+675* &
             ((phicoeff(3,j)*phicoeff(5,j))**2)+1440*phicoeff(3,j)*(phicoeff(4,j)**2)*phicoeff(5,j)+128*(phicoeff(4,j)**4))/(y**13)

        small=small-22*(18*(phicoeff(3,j)**2)*phicoeff(8,j)+42*phicoeff(3,j)*phicoeff(4,j)*phicoeff(7,j)+45*phicoeff(3,j)* &
             phicoeff(5,j)*phicoeff(6,j)+24*(phicoeff(4,j)**2)*phicoeff(6,j)+25*phicoeff(4,j)*(phicoeff(5,j)**2))/(y**12)

        small=small+(27*phicoeff(3,j)*phicoeff(9,j)+32*phicoeff(4,j)*phicoeff(8,j)+35*phicoeff(5,j)*phicoeff(7,j)+18* &
             (phicoeff(6,j)**2))/(y**11)

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)
     else if (i.eq.10) then       
        rfact1=rfact1*(2*i-3)*6/((i+1)*1.0)
        con1=(phicoeff(3,j)**(i-3))*(9*(phicoeff(3,j)**2)-2*(i-2)*y*phicoeff(4,j))
        small=rfact1*con1/(y**19)

        rfact2=rfact2*6*(i-2)*(2*i-5)/((i+1)*1.0*(i-4))
        con2=(phicoeff(3,j)**(i-5))*(15*phicoeff(3,j)*phicoeff(5,j)+8*(i-4)*(phicoeff(4,j)**2))
        small=small+rfact2*con2/(y**17)

        rfact3=rfact3*6*(i-3)*(2*i-5)/((i+1)*1.0*(i-4))
        con3=(phicoeff(3,j)**(i-6))*(9*phicoeff(3,j)*phicoeff(6,j)+10*(i-5)*phicoeff(4,j)*phicoeff(5,j) &
             +(i-6)*(i-5)*16*(phicoeff(4,j)**3)/(9*phicoeff(3,j)))
        small=small-rfact3*con3/(y**16)

        con4=273*(189*(phicoeff(3,j)**3)*phicoeff(7,j)+864*(phicoeff(3,j)**2)*phicoeff(4,j)*phicoeff(6,j)+450* &
             ((phicoeff(3,j)*phicoeff(5,j))**2)+1440*phicoeff(3,j)*(phicoeff(4,j)**2)*phicoeff(5,j)+256*(phicoeff(4,j)**4))
        small=small+con4*phicoeff(3,j)/(y**15)

        con4=108*(phicoeff(3,j)**3)*phicoeff(8,j)+378*(phicoeff(3,j)**2)*phicoeff(4,j)*phicoeff(7,j)
        con4=con4+405*(phicoeff(3,j)**2)*phicoeff(5,j)*phicoeff(6,j)+432*(phicoeff(4,j)**2)*phicoeff(3,j)*phicoeff(6,j)
        con4=con4+450*(phicoeff(5,j)**2)*phicoeff(3,j)*phicoeff(4,j)+160*(phicoeff(4,j)**3)*phicoeff(5,j)
        small=small-52*con4/(y**14)

        con4=243*(phicoeff(3,j)**2)*phicoeff(9,j)+576*phicoeff(3,j)*phicoeff(4,j)*phicoeff(8,j)
        con4=con4+630*phicoeff(3,j)*phicoeff(5,j)*phicoeff(7,j)+324*(phicoeff(6,j)**2)*phicoeff(3,j)
        con4=con4+336*(phicoeff(4,j)**2)*phicoeff(7,j)+720*phicoeff(4,j)*phicoeff(5,j)*phicoeff(6,j)+125*(phicoeff(5,j)**3)
        small=small+2*con4/(y**13)

        con4=(30*phicoeff(3,j)*phicoeff(10,j)+36*phicoeff(4,j)*phicoeff(9,j)+40*phicoeff(5,j)*phicoeff(8,j)  &
             +42*phicoeff(6,j)*phicoeff(7,j))
        small=small-con4/(y**12)

        if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
        else 
           see(i)=ki*small
        endif
        ki=ki*(-1)

       else if (i.eq.11) then
          con1=(24.5*(phicoeff(7,j)**2)+33*phicoeff(11,j)*phicoeff(3,j)+40*phicoeff(10,j)*phicoeff(4,j)+45*phicoeff(5,j)  &
           *phicoeff(9,j)+48*phicoeff(6,j)*phicoeff(8,j))/(y**13)
          con2=-(phicoeff(3,j)*(585*phicoeff(3,j)*phicoeff(10,j)+1404*phicoeff(4,j)*phicoeff(9,j)+1560*phicoeff(5,j)*phicoeff(8,j)&
           +1638*phicoeff(6,j)*phicoeff(7,j))+832*(phicoeff(4,j)**2)*phicoeff(8,j)+1820*phicoeff(4,j)*phicoeff(5,j)*phicoeff(7,j))
          con2=(con2-936*(phicoeff(6,j)**2)*phicoeff(4,j)-975*(phicoeff(5,j)**2)*phicoeff(6,j))/(y**14)
          con3=((phicoeff(3,j)**2)*(7371*phicoeff(3,j)*phicoeff(9,j)+26208*phicoeff(4,j)*phicoeff(8,j)+28665*phicoeff(5,j)* &
           phicoeff(7,j)+14742*phicoeff(6,j)*phicoeff(6,j)))
          con3=(con3+(phicoeff(4,j)**2)*(30576*phicoeff(3,j)*phicoeff(7,j)+11648*phicoeff(4,j)*phicoeff(6,j) &
           +18200*phicoeff(5,j)*phicoeff(5,j)))
          con3=(con3+phicoeff(3,j)*phicoeff(5,j)*(65520*phicoeff(4,j)*phicoeff(6,j)+11375*phicoeff(5,j)*phicoeff(5,j)))/(y**15)
          con4=-(73710*(phicoeff(3,j)**4)*phicoeff(8,j)+(phicoeff(3,j)**3)*(343980*phicoeff(4,j)*phicoeff(7,j) &
           +368550*phicoeff(5,j)*phicoeff(6,j)))
          con4=(con4-(phicoeff(3,j)**2)*phicoeff(4,j)*(589680*phicoeff(4,j)*phicoeff(6,j)+614250*phicoeff(5,j)*phicoeff(5,j)))
          con4=(con4-(phicoeff(4,j)**3)*(436800*phicoeff(5,j)*phicoeff(3,j)+23296*phicoeff(4,j)*phicoeff(4,j)))/(y**16)

          small=con1+con2+con3+con4

          con1=((phicoeff(3,j)**4)*(619164*phicoeff(3,j)*phicoeff(7,j)+3538080*phicoeff(4,j)*phicoeff(6,j) &
           +1842750*phicoeff(5,j)*phicoeff(5,j))) 
          con1=(con1+((phicoeff(3,j)*phicoeff(4,j))**2)*(7862400*phicoeff(3,j)*phicoeff(5,j)+2096640*phicoeff(4,j)*phicoeff(4,j))) &
           /(y**17)
          con2=-((phicoeff(3,j)**5)*(4511052*phicoeff(3,j)*phicoeff(6,j)+30073680*phicoeff(4,j)*phicoeff(5,j)))
          con2=(con2-26732160*(phicoeff(3,j)**4)*(phicoeff(4,j)**3))/(y**18)         
          con3=((phicoeff(3,j)**6)*(28999620*phicoeff(3,j)*phicoeff(5,j)+108265248*phicoeff(4,j)*phicoeff(4,j)))/(y**19)
          con4=-(165297834*(phicoeff(3,j)**8)*phicoeff(4,j))/(y**20)+82648917*(phicoeff(3,j)**10)/(y**21)

          small=small+con1+con2+con3+con4

          if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
          else 
           see(i)=ki*small
          endif
          ki=ki*(-1)
       else if (i.eq.12) then
          con1=(50*phicoeff(10,j)*phicoeff(5,j)+44*phicoeff(11,j)*phicoeff(4,j)+36*phicoeff(12,j)*phicoeff(3,j)+54*phicoeff(6,j) &
           *phicoeff(9,j)+56*phicoeff(7,j)*phicoeff(8,j))/(y**14)
          con2=(phicoeff(3,j)*(1680*phicoeff(10,j)*phicoeff(4,j)+693*phicoeff(11,j)*phicoeff(3,j)+1890*phicoeff(5,j)*phicoeff(9,j) &
           +2016*phicoeff(6,j)*phicoeff(8,j)+1029*(phicoeff(7,j)**2))+phicoeff(4,j)*(1008*phicoeff(4,j)*phicoeff(9,j)+2240* &
           phicoeff(5,j)*phicoeff(8,j)+2352*phicoeff(6,j)*phicoeff(7,j))+phicoeff(5,j)*(1225*phicoeff(5,j)*phicoeff(7,j)+1260 &
           *(phicoeff(6,j)**2)))/(y**15)

          con3=((phicoeff(3,j)**2)*(9450*phicoeff(3,j)*phicoeff(10,j)+34020*phicoeff(4,j)*phicoeff(9,j)+37800*phicoeff(5,j)* & 
           phicoeff(8,j)+39690*phicoeff(6,j)*phicoeff(7,j))+(phicoeff(3,j)*phicoeff(4,j))*(40320*phicoeff(4,j)*phicoeff(8,j) &
           +88200*phicoeff(5,j)*phicoeff(7,j)+45360*(phicoeff(6,j)**2))+(phicoeff(5,j)*phicoeff(6,j))*(47250*phicoeff(3,j)* &
           phicoeff(5,j)+50400*(phicoeff(4,j)**2))+15680*(phicoeff(4,j)**3)*phicoeff(7,j)+17500*phicoeff(4,j)* &
           (phicoeff(5,j)**3))/(y**16)
          con4=(102060*(phicoeff(3,j)**4)*phicoeff(9,j)+(phicoeff(3,j)**3)*(483840*phicoeff(4,j)*phicoeff(8,j)+529200* &
           phicoeff(5,j)*phicoeff(7,j)+272160*(phicoeff(6,j)**2))+(phicoeff(3,j)**2)*(846720*(phicoeff(4,j)**2)*phicoeff(7,j) &
           +1814400*phicoeff(4,j)*phicoeff(5,j)*phicoeff(6,j)+315000*(phicoeff(5,j)**3))+(phicoeff(4,j)**2)*phicoeff(3,j)* &
           (645120*phicoeff(4,j)*phicoeff(6,j)+1008000*(phicoeff(5,j)**2))+179200*(phicoeff(4,j)**4)*phicoeff(5,j))/(y**17)

          small=-con1+con2-con3+con4
 
          con1=(925344*(phicoeff(3,j)**5)*phicoeff(8,j)+(phicoeff(3,j)**4)*(5397840*phicoeff(4,j)*phicoeff(7,j)+5783400* &
           phicoeff(5,j)*phicoeff(6,j))+(phicoeff(3,j)**3)*phicoeff(4,j)*(12337920*phicoeff(4,j)*phicoeff(6,j)+12852000* &
           (phicoeff(5,j)**2))+(phicoeff(4,j)**3)*phicoeff(3,j)*(13708800*phicoeff(3,j)*phicoeff(5,j)+1462272*(phicoeff(4,j)**2))) &
           /(y**18)
          con2=((phicoeff(3,j)**5)*(7287084*phicoeff(3,j)*phicoeff(7,j)+49968576*phicoeff(4,j)*phicoeff(6,j)+26025300* &
          (phicoeff(5,j)**2))+(phicoeff(3,j)**3)*(phicoeff(4,j)**2)*(138801600*phicoeff(3,j)*phicoeff(5,j)+49351680* &
          (phicoeff(4,j)**2)))/(y**19)
          con3=(phicoeff(3,j)**5)*(50860872*(phicoeff(3,j)**2)*phicoeff(6,j)+395584560*phicoeff(3,j)*phicoeff(4,j)*phicoeff(5,j)+ &
           421956864*(phicoeff(4,j)**3))/(y**20)
          con4=((phicoeff(3,j)**7)*(317880450*phicoeff(3,j)*phicoeff(5,j)+1356289920*(phicoeff(4,j)**2)))/(y**21)

          small=small-con1+con2-con3+con4-((phicoeff(3,j)**9)*(1780130520*phicoeff(4,j)*y-801058734*(phicoeff(3,j)**2)))/(y**23)

          if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
          else 
           see(i)=ki*small
          endif
          ki=ki*(-1)

       else if (i.eq.13) then
          p(1)=phicoeff(3,j)
          q(1)=phicoeff(4,j)
          do is=2,12
           p(is)=p(is-1)*phicoeff(3,j)
           if (is.lt.7) then
            q(is)=q(is-1)*phicoeff(4,j)
           endif
          enddo
           
          con1=(60*phicoeff(10,j)*phicoeff(6,j)+55*phicoeff(11,j)*phicoeff(5,j)+48*phicoeff(12,j)*q(1)+39*phicoeff(13,j)*p(1)+63* &
           phicoeff(7,j)*phicoeff(9,j)+32*(phicoeff(8,j)**2))/(y**15)

          con2=(phicoeff(10,j)*(2250*p(1)*phicoeff(5,j)+1200*q(2))+p(1)*(1980*phicoeff(11,j)*q(1)+810*phicoeff(12,j)*p(1)+2430* &
           phicoeff(6,j)*phicoeff(9,j)+2520*phicoeff(7,j)*phicoeff(8,j))+q(1)*(2700*phicoeff(5,j)*phicoeff(9,j)+2880*phicoeff(6,j) &
           *phicoeff(8,j)+1470*(phicoeff(7,j)**2)))
          con2=(con2+phicoeff(5,j)*(1500*phicoeff(5,j)*phicoeff(8,j)+3150*phicoeff(6,j)*phicoeff(7,j))+540*(phicoeff(6,j)**3)) &
           /(y**16)

          con3=(20*p(2)*(2160*phicoeff(10,j)*q(1)+594*phicoeff(11,j)*p(1)+2430*phicoeff(5,j)*phicoeff(9,j)+2592*phicoeff(6,j)* &
           phicoeff(8,j)+1323*(phicoeff(7,j)**2))+5760*p(1)*q(1)*(9*q(1)*phicoeff(9,j)+20*phicoeff(5,j)*phicoeff(8,j)+21* &
           phicoeff(6,j)*phicoeff(7,j)))
          con3=(con3+1800*p(1)*phicoeff(5,j)*(35*phicoeff(5,j)*phicoeff(7,j)+36*(phicoeff(6,j)**2))+640*q(2)*(32*q(1)* &
           phicoeff(8,j)+105*phicoeff(5,j)*phicoeff(7,j)+54*(phicoeff(6,j)**2))+250*(phicoeff(5,j)**2)*(288*q(1)*phicoeff(6,j) &
           +25*(phicoeff(5,j)**2)))/(y**17)
          con4=(9180*p(3)*(15*p(1)*phicoeff(10,j)+72*q(1)*phicoeff(9,j)+80*phicoeff(5,j)*phicoeff(8,j)+84*phicoeff(6,j)* &
           phicoeff(7,j))+73440*p(2)*q(1)*(16*q(1)*phicoeff(8,j)+35*phicoeff(5,j)*phicoeff(7,j)+18*(phicoeff(6,j)**2)))
          con4=(con4+51000*p(1)*(phicoeff(5,j)**2)*(27*p(1)*phicoeff(6,j)+20*q(1)*phicoeff(5,j))+65280*p(1)*q(2)* &
           (14*q(1)*phicoeff(7,j)+45*phicoeff(5,j)*phicoeff(6,j))+21760*q(3)*(12*q(1)*phicoeff(6,j)+25*(phicoeff(5,j)**2))) &
           /(y**18)
           
          small=con1-con2+con3-con4
          
          con1=(49572*p(4)*(27*p(1)*phicoeff(9,j)+160*q(1)*phicoeff(8,j)+175*phicoeff(5,j)*phicoeff(7,j)+90*(phicoeff(6,j)**2)) &
           +1321920*p(2)*q(2)*(14*p(1)*phicoeff(7,j)+16*q(1)*phicoeff(6,j)+25*(phicoeff(5,j)**2)))
          con1=(con1+275400*p(3)*phicoeff(5,j)*(144*q(1)*phicoeff(6,j)+25*(phicoeff(5,j)**2))+52224*q(4)*(225*p(1)*phicoeff(5,j) &
           +8*q(2)))/(y**19)
          con2=(5651208*p(5)*(2*p(1)*phicoeff(8,j)+14*q(1)*phicoeff(7,j)+15*phicoeff(5,j)*phicoeff(6,j))+9418680*p(4)*q(1)* &
           (24*q(1)* phicoeff(6,j)+25*(phicoeff(5,j)**2))+13395456*p(2)*q(3)*(25*p(1)*phicoeff(5,j)+4*q(2)))/(y**20)
          con3=(14128020*p(6)*(6*p(1)*phicoeff(7,j)+48*q(1)*phicoeff(6,j)+25*(phicoeff(5,j)**2))+251164800*p(4)*q(2)*(9*p(1)* &
           phicoeff(5,j)+4*q(2)))/(y**21)
          con4=(63576090*p(7)*(9*p(1)*phicoeff(6,j)+80*q(1)*phicoeff(5,j))+6329352960.0*p(6)*q(3))/(y**22)

          small=small+con1-con2+con3-con4

          con1=(699336990*p(8)*(5*p(1)*phicoeff(5,j)+24*q(2)))/(y**23)
          con2=(877350042*p(10)*(22*q(1)*y-9*p(2)))/(y**25)

          small=small+con1-con2

          if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
          else 
           see(i)=ki*small
          endif
          ki=ki*(-1)

       else if (i.eq.14) then
          p(1)=phicoeff(3,j)
          q(1)=phicoeff(4,j)
          do is=2,12
           p(is)=p(is-1)*phicoeff(3,j)
           if (is.lt.7) then
            q(is)=q(is-1)*phicoeff(4,j)
           endif
          enddo
          r(1)=phicoeff(5,j)
          r(2)=r(1)**2

          con1=(70*phicoeff(10,j)*phicoeff(7,j)+66*phicoeff(11,j)*phicoeff(6,j)+60*phicoeff(12,j)*r(1) &
           +52*phicoeff(13,j)*q(1)+42*phicoeff(14,j)*p(1)+72*phicoeff(8,j)*phicoeff(9,j))/(y**16)
          con2=(320*phicoeff(10,j)*(9*p(1)*phicoeff(6,j)+10*q(1)*r(1))+176*phicoeff(11,j)*(8*q(2)+15* &
           p(1)*r(1))+72*p(1)*(13*p(1)*phicoeff(13,j)+32*p(1)*q(1)))
          con2=(con2+48*p(1)*(32*(phicoeff(8,j)**2)+63*phicoeff(7,j)*phicoeff(9,j))+128*q(1)*( &
           27*phicoeff(6,j)*phicoeff(9,j)+28*phicoeff(7,j)*phicoeff(8,j))+120*r(1)*(15*r(1)*phicoeff(9,j) &
           +32*phicoeff(6,j)*phicoeff(8,j))+56*phicoeff(7,j)*(35*r(1)*phicoeff(7,j)+36*(phicoeff(6,j)**2)))
          con2=con2/(y**17)
          con3=(2448*p(2)*(25*r(1)*phicoeff(10,j)+22*q(1)*phicoeff(11,j)+6*p(1)*phicoeff(12,j)+27*phicoeff &
           (6,j)*phicoeff(9,j)+28*phicoeff(7,j)*phicoeff(8,j)))
          con3=con3+(1632*p(1)*q(1)*(40*q(1)*phicoeff(10,j)+90*r(1)*phicoeff(9,j)+96*phicoeff(6,j)*phicoeff &
           (8,j)+49*(phicoeff(7,j)**2))+8160*p(1)*r(1)*(10*r(1)*phicoeff(8,j)+21*phicoeff(6,j)*phicoeff(7,j)))     
          con3=con3+(4352*q(2)*(20*r(1)*phicoeff(8,j)+21*phicoeff(6,j)*phicoeff(7,j))+2720*q(1)*r(1)*(35*r(1) &
           *phicoeff(7,j)+36*(phicoeff(6,j)**2)))+29376*p(1)*(phicoeff(6,j)**3)+26112*q(3)*phicoeff(9,1)+ &
           34000*r(2)*r(1)*phicoeff(6,j)
          con3=con3/(y**18) 
          con4=(5508*p(3)*(160*q(1)*phicoeff(10,j)+33*p(1)*phicoeff(11,j)+180*r(1)*phicoeff(9,j)+192*phicoeff(6,j) &
           *phicoeff(8,j)+98*(phicoeff(7,j)**2))+176256*p(2)*q(1)*(9*q(1)*phicoeff(9,j)+20*r(1)*phicoeff(8,j)+21* &
           phicoeff(6,j)*phicoeff(7,j)))
          con4=con4+(55080*p(2)*r(1)*(35*r(1)*phicoeff(7,j)+36*(phicoeff(6,j)**2))+39168*p(1)*q(2)*(32*q(1)*phicoeff &
           (8,j)+105*r(1)*phicoeff(7,j)+54*(phicoeff(6,j)**2)))
          con4=(con4+15300*p(1)*r(2)*(288*q(1)*phicoeff(6,j)+25*r(2))+6528*q(2)*(56*q(2)*phicoeff(7,j)+240*q(1)*r(1)* &
           phicoeff(6,j)+125*r(2)*r(1)))/(y**19)
          small=con1-con2+con3-con4
          con1=(627912*p(4)*(3*p(1)*phicoeff(10,j)+18*q(1)*phicoeff(9,j)+20*r(1)*phicoeff(8,j)+21*phicoeff(6,j)* &
           phicoeff(7,j))+1674432*p(3)*q(1)*(16*q(1)*phicoeff(8,j)+35*r(1)*phicoeff(7,j)+18*(phicoeff(6,j)**2)))
          con1=con1+(139536*p(2)*(225*p(1)*r(2)*phicoeff(6,j)+224*q(3)*phicoeff(7,j)+720*q(2)*r(1)*phicoeff(6,j)) &
           +1162800*p(1)*q(1)*r(2)*(3*p(1)*r(1)+32*q(2))+1984512.0*q(4)*(9*p(1)*phicoeff(6,j)+2*p(1)*phicoeff(5,j)))
          con1=con1/(y**20)
          con2=(1883736.0*p(5)*(9*p(1)*phicoeff(9,j)+64*q(1)*phicoeff(8,j)+70*r(1)*phicoeff(7,1)+36*(phicoeff(6,j)**2)) &
           +1046520.0*p(4)*(336*q(2)*phicoeff(7,j)+720*q(1)*r(1)*phicoeff(6,j)+125*r(2)*r(1)))
          con2=(con2+33488640.0*p(3)*q(2)*(16*q(1)*phicoeff(6,j)+25*r(2))+1984512.0*q(4)*p(1)*(225*p(1)*r(1)+16*q(2))) &
           /(y**21)
          con3=(11302416.0*p(6)*(12*p(1)*phicoeff(8,j)+98*q(1)*phicoeff(7,j)+105*r(1)*phicoeff(6,j))+158233824.0*p(5) &
           *q(1)*(24*q(1)*phicoeff(6,j)+25*r(2))+93768192.0*p(3)*q(3)*(75*p(1)*r(1)+16*q(2)))/(y**22)
          con4=(46622466.0*p(7)*(21*p(1)*phicoeff(7,j)+192*q(1)*phicoeff(6,j)+100*r(2))+2320762752.0*p(5)*q(2)*(15* &
           p(1)*r(1)+8*q(2)))/(y**23)
          small=small+con1-con2+con3-con4
    
          con1=(714877812.0*p(7)*(9*p(2)*phicoeff(6,j)+90*p(1)*q(1)*r(1)+128*q(3)))/(y**24)
          con2=(12867800616.0*p(9)*(3*p(1)*r(1)+16*q(2)))/(y**25)
          con3=(26320501260.0*p(11)*(8*q(1)*y-3*p(2)))/(y**27)
          small=con1-con2+con3

          if (i.lt.m) then
           see(i)=ki*small-phicoeff(i+1,j)/(y**(i+1))
          else 
           see(i)=ki*small
          endif
          ki=ki*(-1)
     endif
  enddo

  return
end subroutine scheme

! See documentation for phifactors

subroutine phifactors(m,j)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16)    

  integer       :: m,j,icj(0:50),MIT,i,nk,n,if,mmax
  integer(dp1)  :: L(0:50)
  real(dp)      :: phicoeff(0:15,50),xr(50),fracL(0:50),see(15),x,y,sum,y1,y2,GA

  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS7/ see

  x=-phicoeff(1,j)
  do i=0,m
     if (i.eq.0.or.i.eq.1) then
        nk=2
     else
        nk=i
     endif
     sum=0.0d0

     !        for large m uses the fact that m!=gamma(m+1) and works out log(m!)=log(gamma(m+1)) and then
     !        combines logs to get combinatorial function n!/((n-i)!i!).

     if (m.gt.12) then
        ga=(i+1)*1.0
        call lgam(ga,y)
        do n=nk,m
           ga=(n+1)*1.0
           call lgam(ga,y1)
           ga=(n-i+1)*1.0
           call lgam(ga,y2)
           sum=sum+exp(y1-y2-y)*(x**(n-i))*see(n-1)
        enddo
     else

        !        for small m just work out factorials long hand

        y=1.0
        do if=2,i
           y=y*if
        enddo
        do n=nk,m
           y1=1.0
           do if=2,n
              y1=y1*if
           enddo
           y2=1.0
           do if=2,(n-i)
              y2=y2*if
           enddo
           sum=sum+y1*(x**(n-i))*see(n-1)/(y*y2)
        enddo
     endif
     phicoeff(i,j+1)=sum
     if (i.gt.2.and.abs(phicoeff(i,j+1)).gt.0.5) then
        phicoeff(i,j+1)=phicoeff(i,j+1)-anint(phicoeff(i,j+1),dp1)
     endif
  enddo

  return
end subroutine phifactors


subroutine computecnewton(n,nit,con1,tol,c,c2,numb)

  !       Computes value of saddle point c for large dimensional schemes.
  !       Specifically it evaluate c such that
  !
  !       f(c)=(m-phicoeff(1,nit))-Sum[q=2..n] of q*phicoeff(q,nit)*(c^(q-1))=0
  !
  !       for given integer value. Uses standard 1D Newton-Raphson to a tolerance tol>0
  !       con1=(m-phicoeff(1)), n=scheme value only used for m>=6. On entry c is a guess of the answer.
  !       On exit c=final answer.

  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16) 

  integer       :: n,nit,numb,icj(0:50),MIT,itermax,i,j,mmax,itry
  integer(dp1)  :: L(0:50)
  real(dp)      :: con1,tol,c,c2,cold,sum1,sum2
  real(dp)      :: phicoeff(0:15,50),xr(50),fracL(0:50)


  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax

  cold=c
  if (cold.lt.0.0) then
   cold=1000.0
  endif
  itermax=12
  itry=0
 2600  do i=1,itermax
        sum1=n*phicoeff(n,nit)*cold
        sum2=-n*(n-1)*phicoeff(n,nit)
        do j=(n-1),2,-1
         sum1=(sum1+j*phicoeff(j,nit))*cold
         sum2=(sum2*cold-j*(j-1)*phicoeff(j,nit))
        enddo
        sum1=con1-sum1
        c=cold-sum1/sum2
        if (abs(c-cold).lt.tol) then
         if (c.lt.L(nit-1)) then
          if (numb.eq.1) then
           return
          else

! In this instance numb=2 and we want to find two saddle points. First one,
! the smallest already found, now find second one.

           goto 2620
          endif
         else

! converged but onto wrong saddle c>L(nit-1), try again

          goto 2610
         endif
        else
         cold=c
        endif
       enddo
       sum1=con1
       do j=2,n
        sum1=sum1-j*phicoeff(j,nit)*(c**(j-1))
       enddo
       if (abs(sum1).lt.tol) then
        if (numb.eq.1) then
         return
        else
         goto 2620
        endif
       endif
2610   itry=itry+1
       if (itry.eq.1.and.con1.lt.1.0) then
        cold=2.95*con1
        goto 2600
       else if (itry.eq.1.and.con1.gt.1.0) then
        cold=2*con1
        goto 2600
       else
        write(6,*) 'no convergence in computecnewton after 12 iterations'
        write(6,*) L(nit),L(nit-1)
        write(6,*) 'NIT and n= ',NIT,n
        write(6,*) 'con1= ',con1
        write(6,*) 'c= ',c
        do j=1,n
         write(6,*) j,phicoeff(j,nit)
        enddo
        stop
       endif

!  now try to find second saddle c2 such that 0<c<c2<L(nit-1)

2620   cold=0.5*(L(nit-1)-1.0+c)
       do i=1,itermax
        sum1=n*phicoeff(n,nit)*cold
        sum2=-n*(n-1)*phicoeff(n,nit)
        do j=(n-1),2,-1
         sum1=(sum1+j*phicoeff(j,nit))*cold
         sum2=(sum2*cold-j*(j-1)*phicoeff(j,nit))
        enddo
        sum1=con1-sum1

        c2=cold-sum1/sum2
        if (abs(c2-cold).lt.tol) then
         if ((c2.lt.L(nit-1)).and.(c2.gt.c)) then
          return
         else
          write(6,*) 'second saddle converged but greater than L(nit-1)'
          write(6,*)  L(nit),L(nit-1)
          write(6,*) 'con1= ',con1
          write(6,*) 'two sads= ',c,c2
          do j=1,n
           write(6,*) j,phicoeff(j,nit)
          enddo
          stop
         endif
        else 
         cold=c2
        endif
       enddo
       stop


end subroutine computecnewton

! subroutine mgaussum computes the kernel sum and iterates this value up the hierarchical chain to find the original quartic sum

subroutine mgausssum(m,ip,csum)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)  
  integer, parameter :: dp1 = selected_int_kind(16)

  integer       :: m,ip,MIT,icj(0:50),k,j,m1(50),mmax,nit,istart,numb,jj
  integer(dp1)  :: L(0:50),nosubdivide,Lnosub,i
  real(dp)      :: p,sp,tpp,tpm,sum,tol,c,gc,c2,rminorcor3,bb
  real(dp)      :: phicoeff(0:15,50),xr(50),fracL(0:50),x,y,rminorcor,denn,rminorcor2,h1,h2,h3
  complex(kind=16)   :: csum,EPI4,fn,pp,c1,qq

  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  COMMON/PARS5/ m1
  COMMON/DIVISION/ nosubdivide,Lnosub

  numb=1
  tol=1d-7

! if frac(phi2) is practically zero, so L(1)=2*frac(phi2)*L(0)+3*phi3*(L(0)^2)
! is close to zero then no subdivision is possible and the
! initial sum must be computed longhand. But this is an extremely rare case.
  
  if (MIT.eq.1) then
     csum=(1.0,0.0)
     do i=1,L(0)
        c1=i*phicoeff(3,1)
        do jj=2,1,-1
           c1=(c1+phicoeff(jj,1))*i
        enddo
        c1=(0.0,1.0)*tpm*c1
        csum=csum+exp(c1)
     enddo
     if (icj(1).eq.-1) then
        csum=conjg(csum)
     endif
     return
  endif

! But usually MIT>1 and subdivision is possible. First find kernel sum.
  
  if (nosubdivide.eq.0) then
   if (L(MIT-1).eq.0.or.L(MIT-1).eq.(-1)) then
     csum=(1.0,0.0)
   else
     if (MIT.gt.1) then
        h1=3*phicoeff(3,MIT-1)/(xr(MIT-1)**2)
        h3=h1*h1
        if (m1(MIT-1).gt.3) then
           h2=phicoeff(4,MIT-1)/(xr(MIT-1)**3)
        else
           h2=0.0
        endif   
        denn=1-phicoeff(1,MIT-1)*(h1-phicoeff(1,MIT-1)*(6*h2-1.5*h3))
        rminorcor2=1.5*h3-6*h2
        rminorcor=h1+3*phicoeff(1,MIT-1)*(h3-4*h2)
        if (m1(MIT-1).lt.5) then
           rminorcor3=-22*h1*h2
        else
           rminorcor3=0.5*phicoeff(5,MIT-1)/(xr(MIT-1)**4)-22*h1*h2
        endif   
     endif

     ! Coefficients necessary to find the amplitude of each term in the kernel sum.
     ! Now find kernel sum itself.
     
     csum=(0.0,0.0)
     do i=0,L(MIT-1)
        y=i*1.0
        fn=(0.0,1.0)*tpm*y*(phicoeff(1,MIT)+y*(phicoeff(2,MIT)+y*phicoeff(3,MIT)))
        if (m1(MIT).gt.3) then
         sum=0.0
         do k=m1(MIT),4,-1
          sum=y*(phicoeff(k,MIT)+sum)
         enddo
         fn=fn+(0.0,1.0)*tpm*(y**3)*sum
        endif
        if (MIT.gt.1) then
           csum=csum+exp(fn)/(denn+i*(rminorcor-i*(rminorcor2-i*rminorcor3)))
        else   
           csum=csum+exp(fn)
        endif   
     enddo
     if (MIT.gt.1.and.phicoeff(1,MIT-1).gt.0.0) then
        csum=csum-1.0
     endif
   endif
   if (icj(MIT).eq.(-1)) then
     csum=conjg(csum)
   endif
  else

!    nosubdivide=1 initiates double saddle kernel sum


   csum=(0.0,0.0)
   nit=MIT-1

!   This is the MIT sum but the coefficients for MIT are not found - instead you work with saddle points evaluated
!   from the MIT-1 coefficient equation. Hence nit=MIT-1

   if (MIT.ge.3.and.L(nit).lt.7) then
    MIT=MIT-1
    nit=MIT-1
    Lnosub=L(nit)
   endif
   if (phicoeff(1,nit).gt.0.0) then
    istart=1
   else
    istart=0
   endif
   do i=istart,L(nit)
    y=i*1.0-phicoeff(1,nit)
    if (phicoeff(3,nit).gt.0.0) then
     c=phicoeff(2,nit)*(sqrt(1+3*phicoeff(3,nit)*y/(phicoeff(2,nit)**2))-1)/(3*phicoeff(3,nit))
    else
     c=i*0.834*L(nit-1)/(L(nit)*1.0)
    endif
    if (Lnosub.gt.0.and.i.gt.Lnosub) then
     numb=2
    endif

!   Compute the appropriate saddle point or points numerically.    
    
    call computecnewton(m1(nit),nit,y,tol,c,c2,numb)
    gc=0.0
    denn=0.0
    do j=m1(nit),1,-1
       gc=c*(phicoeff(j,nit)+gc)
       if (j.ge.3) then
          denn=denn+j*(j-1)*phicoeff(j,nit)*(c**(j-2))/2.0
       endif
    enddo
    gc=i*c-gc
    fn=(0.0,1.0)*tpm*gc
    pp=1/sqrt(abs(1+denn/phicoeff(2,nit)))
    
!  Note equivalent to rminorcor correction to pp here

    if (nit.gt.2) then
       bb=1+0.75*phicoeff(3,nit-1)*(c-phicoeff(1,nit-1))/(phicoeff(2,nit-1)**2)
    else
       bb=1.0
    endif

    csum=csum+exp(fn)*pp/bb
    if (numb.eq.2) then
     gc=0.0
     denn=0.0  
     do j=m1(nit),1,-1
        gc=c2*(phicoeff(j,nit)+gc)
        if (j.ge.3) then
           denn=denn+j*(j-1)*phicoeff(j,nit)*(c2**(j-2))/2.0
        endif   
     enddo
     gc=i*c2-gc
     fn=(0.0,1.0)*tpm*gc-(0.0,1.0)*p/2.0
     pp=1/sqrt(abs(1+denn/phicoeff(2,nit)))

     if (nit.gt.2) then
        bb=1+0.75*phicoeff(3,nit-1)*(c2-phicoeff(1,nit-1))/(phicoeff(2,nit-1)**2)
     else
        bb=1.0
     endif
     
     csum=csum+exp(fn)*pp/bb
    endif
   enddo

   c=2*phicoeff(2,nit)-nint(2*phicoeff(2,nit))
   if (c.gt.0.0) then
    csum=conjg(csum)
   endif
  endif

! Appropriate kernel sum evaluated. Now recursively evaluate those values in hierarchical up to the actual sum (csum) required.

  do k=(MIT-1),1,-1
   x=abs(xr(k))
   if (k.eq.(MIT-1).and.nosubdivide.eq.1) then
    pp=(0.0,0.0)
   else  
    pp=(0.0,1.0)*tpp*phicoeff(0,k+1)
   endif  
   c1=exp(pp)/(EPI4*sqrt(x))

   call q(m,k,ip,qq,Lnosub)
   Lnosub=-100
   csum=c1*csum+qq
   if (icj(k).eq.-1) then
    csum=conjg(csum)
   endif
   if (k.gt.1.and.phicoeff(1,k-1).gt.0.0) then
    csum=csum-1.0
   endif
  enddo

  return

end subroutine mgausssum

! routine computarrs sets up the appropriate hierarchical chain starting from the initial quartic sum

subroutine computarrs(m,initialcoeff,N,KCUT,KLP)

  implicit none

  integer, parameter :: dp = selected_real_kind(33)  
  integer, parameter :: dp1 = selected_int_kind(16)


  integer       :: m,icj(0:50),ix,j,jj,MIT,m1(50),ichange,mmax,negL
  integer(dp1)  :: N,KCUT,KLP,L(0:50),nosubdivide,Lnosub
  real(dp)      :: initialcoeff(3),phicoeff(0:15,50),xr(50),fracL(0:50),x0,sum
  real(dp)      :: etaa,small,smalle

  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS5/ m1
  COMMON/DIVISION/ nosubdivide,Lnosub


  nosubdivide=0
  m1(1)=m
  if (N.lt.1) then
     write(2,*) 'N is too small for computarrs'
     return
  endif

  !   Initial setup and first iteration

  L(0)=N     
  icj(0)=1
  negL=0

  x0=2*initialcoeff(2)
  ix=nint(x0)
  xr(1)=x0-ix
  phicoeff(2,1)=xr(1)/2.0
  if (mod(ix,2).eq.1) then
     if (initialcoeff(1).lt.0.0) then
        phicoeff(1,1)=initialcoeff(1)+0.5
     else
        phicoeff(1,1)=initialcoeff(1)-0.5
     endif
  else
     phicoeff(1,1)=initialcoeff(1)
  endif
  do jj=3,m
     phicoeff(jj,1)=initialcoeff(jj)
  enddo
  if (xr(1).gt.0.0) then
     icj(1)=-1
  else
     icj(1)=1
     phicoeff(1,1)=-phicoeff(1,1)
     xr(1)=abs(xr(1))
     phicoeff(2,1)=abs(phicoeff(2,1))
     do jj=3,m
        phicoeff(jj,1)=-phicoeff(jj,1)
     enddo
  endif
  sum=m*phicoeff(m,1)
  etaa=real(L(0),dp)
  do jj=m-1,1,-1
   sum=sum*etaa+jj*phicoeff(jj,1)
  enddo


  L(1)=floor(sum,dp1)
  fracL(1)=sum-L(1)*1.0
  if (phicoeff(1,1).gt.0.0) then
     etaa=1.0-phicoeff(1,1)
  else
     etaa=1.0+phicoeff(1,1)
  endif

  !     Extra check. If sum very close to 1 on first iteration this is likely to lead to significant errors
  !     i.e. sum>1 but fracL(1)<<1 then setting L(1)=1 is unlikely to be good. So instead set it to L(1)=0
  !     and just do a full sum with no iterations to be on the safe side.

  if (L(1).eq.1.and.abs(kcut*fracL(1)).lt.1.0) then
     L(1)=0
  endif
  if ((L(1).le.0)) then
     MIT=1
     return
  endif

  !  Main iteration loop-maximum length of chain is set to 50: in practice 10 is usually sufficient.

  do j=2,50
     ichange=0
     m1(j)=m1(j-1)

985  if (ichange.eq.0) then

!   if counter ichange=0 we are trying to reduce proto-generator sum without changing its order. So we start by computing
!   the small factor which checks the validity of the scheme and a new set of phi coefficients in anticipation of it
!   being successful

      call scheme(m1(j-1),j-1,small)
      call phifactors(m1(j-1),j-1)
     else

!    but if ichange=1 we are trying to reduce the proto-generator sum by moving to a higher order. This means
!    we have to set phi(m1(j))=0 in the previous coeffecients and then compute a new set of phi factors including phi(m1(j))
!    for the next set. So the call to phifactors must come before the call to scheme. NB m1(j)>m1(j-1). 

      call phifactors(m1(j),j-1)
      call scheme(m1(j),j-1,small) 
     endif

     ix=nint(2*phicoeff(2,j))
     xr(j)=2*phicoeff(2,j)-ix
     phicoeff(2,j)=0.5*xr(j)
     phicoeff(1,j)=phicoeff(1,j)-anint(phicoeff(1,j),dp1)
     if (mod(abs(ix),2).eq.1) then
        if (phicoeff(1,j).lt.0.0) then
           phicoeff(1,j)=phicoeff(1,j)+0.5
        else
           phicoeff(1,j)=phicoeff(1,j)-0.5
        endif
     endif
     if (xr(j).gt.0.0) then
        icj(j)=-1
     else
        xr(j)=-xr(j)
        icj(j)=1
        do jj=1,m1(j)
           phicoeff(jj,j)=-phicoeff(jj,j)
        enddo
     endif
     if (L(j-1).le.0) then
        L(j)=-2
     else
      sum=m1(j)*phicoeff(m1(j),j)
      etaa=real(L(j-1),dp)
      do jj=m1(j)-1,1,-1
       sum=sum*etaa+jj*phicoeff(jj,j)
      enddo
      L(j)=floor(sum,dp1)
      fracL(j)=sum-real(L(j),dp)
      if (phicoeff(1,j).gt.0.0) then
       etaa=1.0-phicoeff(1,j)
      else
       etaa=1.0+phicoeff(1,j)
      endif
     endif

! The code above readies things for a next reduction stage before the initial one is confirmed.

     x0=(L(j-1)-phicoeff(1,j-1))
     smalle=(x0**(m1(j)+1))*small


! smalle is the parameter that determines if this new sum length L(j-1) is appropriate for this order sum. If smalle is
! too big then the order must be raised by 1 and the calculation repeated until (if or when) the order exceeds mmax.

     if (abs(smalle).gt.0.001) then
        m1(j)=m1(j)+1

!   Here we reach a critical point. In order to find the next generation of the sum the order is increased, one step at a time.

        if (m1(j).gt.mmax.or.negL.eq.1) then
         small=(3*L(j-2)*(phicoeff(3,j-1)+2*phicoeff(4,j-1)*L(j-2))/phicoeff(2,j-1))

!   potential double saddle exit here

         if ((small.lt.0.0.and.abs(small).lt.1.00)) then
          nosubdivide=1
          MIT=j
          Lnosub=L(j-1)
          return
         else if ((small.gt.0.0.and.abs(small).lt.1.30)) then
          nosubdivide=1
          MIT=j
          Lnosub=L(j-1)
          return
         endif
         MIT=j-1
         return
        endif

!  Increase value of m and then go back and try to achieve a successful reduction with a higher order.

        phicoeff(m1(j),j-1)=0.0
        ichange=1
        goto 985
     endif

! After any successful reduction of a proto-generator sum ichange is reset to zero.

         if (ichange.eq.1) then
          ichange=0
         endif
   
         etaa=2*phicoeff(2,j)+6*phicoeff(3,j)*L(j-1)

!    Here is where double saddles code is set.

         if (etaa.lt.0.0) then
          nosubdivide=1
          MIT=j+1
          Lnosub=L(j)
          etaa=8*phicoeff(4,j)*phicoeff(2,j)/(3*(phicoeff(3,j)**2))
          if (etaa.lt.1.0) then
           etaa=(sqrt(1-etaa)-1)*phicoeff(3,j)/(4*phicoeff(4,j))
           sum=m1(j)*phicoeff(m1(j),j)
           do jj=m1(j)-1,1,-1
            sum=sum*etaa+jj*phicoeff(jj,j)
           enddo
           L(j)=max0(L(j),floor(sum,dp1))
           if (Lnosub.gt.KCUT.and.(L(j).lt.L(j-1))) then
            return
           else
            MIT=j
            Lnosub=L(j-1)
            return
           endif
          else
           MIT=j
           Lnosub=L(j-1)
           return
          endif
         endif

!   Otherwise exit as before


     if ((L(j).le.KCUT.and.L(j-1).le.KCUT).or.m1(j).eq.mmax) then
        MIT=j
        if (L(MIT-1).lt.KLP.and.MIT.ge.3) then
           MIT=MIT-1
        endif
        return
     endif
     if (L(j).le.0.0) then
        negL=1
     endif
     if (L(j).gt.L(j-1)) then
        MIT=j
        return
     endif
  enddo
  write(2,*) 'no convergence after 50 iterations'

end subroutine computarrs

! See documentation for more details of the additive term qq which is calculated in routine q.

subroutine q(m,nit,ip,qq,Lnosub)
 
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: dp1 = selected_int_kind(16) 

  integer       :: m,nit,ip,j,ip1,icj(0:50),MIT,m1(50),jbot,mmax,numb
  integer(dp1)  :: Lnosub,L(0:50),i
  real(dp)      :: p,sp,tpp,tpm,e3,sum,ps1,ps2
  real(dp)      :: phicoeff(0:15,50),xr(50),fracL(0:50),sx,con1,con2,con3,z,wm,c,beta,c2
  real(dp)      :: xbeta,c1,ecor,gc,sav1,sav2,tol,gcddashmt
  complex(kind=16)   :: qq,epi4,t1,t2,t3,t4,t5,endpoint
  complex(kind=16)   :: cr1,cr2,cr3,fn,sum0,sum1

  COMMON/PARS1/ phicoeff,fracL,xr,L,icj,MIT,mmax
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4
  COMMON/PARS5/ m1

  !      COMPUTES THE ADDITIVE TERM QQ IN THE RECURSIVE SERIES WHICH MAKES UP THE HIERARCHICAL CHAIN LINKING THE KERNEL SUM TO THE 
  !      VALUE OF THE QUARTIC SUM WHICH IS TO BE COMPUTED. SPECIFICALLY, EACH LINK IN THE CHAIN SATISFIES THE FOLLOWING EQUATION

  !       GENERAL GAUSS SUM LENGTH L(NIT-1) = MULTIPLICATIVE FACTOR*(NEW GENERAL GAUSS SUM LENGTH L(NIT)+QQ)

  !      WHERE L(NIT)~ABS(xr(nit))*L(NIT-1) WITH ABS(x(nit))<0.5. RECURSIVELY APPLIED IT CAN BE USED TO ESTIMATE 
  !      A VERY LONG GAUSS SUM (TRILLIONS UPON TRILLIONS OF SUMMANDS) FROM A RELATIVELY SHORT ONE (HUNDREDS OF SUMMANDS)
  !      L(0)>L(1)>L(2)....>L(MIT). NB the first sum is of length L(0) and coefficients phicoeff(1 to m,1), the nit sum
  !      is of length L(nit-1) and coefficients phicoeff(1 to m, nit) etc.

  !      Inputs
  !       m: classification of general gauss sum, m=3 cubic, m=4 quartic etc
  !       nmax: used to calulate wm factor. Defined to be largest integer
  !              less than (m+1-q)/l   m>=q>=2,  1<=l<=q. For m=3 nmax=3.
  !       nit: current iteration value, lies between 1 and MIT.
  !       ip: >=1 endpoint factor. Certain integrals derived from terms near the beginning and end of each
  !           Gauss sum are not asymptotically small in the way as the same integrals derived from the main
  !           bulk of the terms. These integrals near the endpoints must be computed and summed. The sums
  !           converge but relatively slowly, so ip must be large enough to ensure the sums have converged
  !           close enough to their actual values. I have found ip~3-4 to be adequate.
  !      Output
  !       qq:  The correction factor to the recursive scheme. It is made up of 5 sub-terms I have denoted
  !       t1-5. These terms are derived from integrals that appear in the Euler-Maclaurin summation formula
  !       which forms the crux of the recursive scheme.

  !       xcube=(xr(nit)^3)
  !       e1=-phicoeff(1,nit)/xr(nit)-3*(phicoeff(1,nit)^2)*phicoeff(3,nit)/xcube
  !       e2=1/(2*xr(nit))+3*phicoeff(1,nit)*phicoeff(3,nit)/xcube
  !       e3=-phicoeff(3,nit)/xcube

  !      e1,e2,e3 are the first 3 coefficients of the shorter Gauss sum on the RHS of the equation
  !      derived from the corresponding coefficients of the longer sum phicoeff(1-3,nit). These are used
  !      in the calculation of term t5.

  !       fn=(0.0,1.0)*tpp*(phicoeff(1,nit)^2)*(0.5/xr(nit)-phicoeff(1,nit)*e3)
  !       con5=exp(fn)

  !      con5 is actually the muliplicative factor in the recursive scheme. Also used for t5


       tol=1e-7
       sx=sqrt(xr(nit))
       sum=0.0
       gcddashmt=0.0
       numb=1
      
       e3=real(L(nit-1),dp)
       do i=1,m1(nit)
        sum=sum+phicoeff(i,nit)*(e3**i)
        if (i.ge.3) then
          gcddashmt=gcddashmt+i*(i-1)*phicoeff(i,nit)*(e3**(i-2))
        endif        
       enddo
       gcddashmt=0.5*(gcddashmt+2*phicoeff(2,nit))
       fn=(0.0,1.0)*tpm*sum
       endpoint=exp(fn)
       t3=0.5*(1.0+endpoint)

!      t3 easiest term to calculate, given by 0.5*(terms at k=0 and k=L(nit-1)) that appear in E-M summation formula.

       con1=1.0-fracL(nit)
       con2=real(L(nit),dp)+fracL(nit)+1.0
       call psi(con2,ps2)
       call psi(con1,ps1)
       sav1=ps1
       con2=ps2-ps1

!      ignore psi(2,z) correction here, which is included in MAPLE version

       call psi(phicoeff(1,nit)+1.0,ps2)
       call psi(L(nit)+1.0-phicoeff(1,nit),ps1)
       sav2=ps1
       con3=ps2-ps1
       t1=(0.0,1.0)*(conjg(endpoint)*con2-con3)/tpp

!      t1 is the contribution from the int[0..L(nit-1)] psi(x)*f'(x)*exp(2*p*I*f(x))dx that doesn't
!      involve saddles. PSi(x)=-1 sum[k=1 to infinity)]sin(2*p*k*x)/(p*k)

       z=sp*con1/sx
       call erf(z,6,cr1)
       cr2=1.0-cr1
       fn=-(0.0,1.0)*p*(con1**2)/xr(nit)
       cr3=p*exp(fn)*cr2/(sx*EPI4)-1/con1-1/(real(L(nit),dp)+1.0)
       t2=(0.0,1.0)*conjg(endpoint)*cr3/tpp

!      t2 acts as a correction to t1. t1 is derived using an asymptotic expansion for the erfc(z), with z large.
!      Normally this is good enough, but if phicoeff(1,nit) is very near zero, or fracL(nit) very near zero it
!      means that z is not really large at the first and last points to apply the asymptotic expansion.
!      Effectively t2 knocks out the inaccurate asmptotic terms used for erfc(z) in t1 and replaces them with their
!      exact erfc(z) values. Usually abs(t2) is relatively small, indicating that the asymptotic terms are pretty good.

       if (phicoeff(1,nit).gt.0.0) then
        con1=phicoeff(1,nit)
       else
        con1=phicoeff(1,nit)+1.0
       endif
       z=sp*con1/sx
       call erf(z,6,cr1)
       cr2=1.0-cr1
       fn=-(0.0,1.0)*p*(con1**2)/xr(nit)
       cr3=(0.0,1.0)*exp(fn)*cr2/(2*sx*EPI4)
       if (phicoeff(1,nit).gt.0.0) then
        t4=-(0.0,1.0)*conjg(endpoint)/(tpp*(L(nit)+fracL(nit)))+cr3
       else
        t4=(0.0,1.0)*(phicoeff(1,nit)/(con1)-1.0)/tpp+cr3
       endif       

!      If phicoeff(1,nit)>0, there is an integral that does not contain a saddle which must be estimated.
!      This is the first correction made in t4. In this case the new sum starts at k=1,2,...L(nit).
!      If phicoeff(1,nit)<0, the integral has a saddle at k=0 effectively. This saddle term acts as the first point
!      of the new shortened sum k=0,1... L(nit). In addition to this zeroth saddle, there is a correction that is 
!      included in t4. Overall the abs(t4) is usually small compared to t1,t3 and t5.

!      Next is t5. This is the most difficult part of the correction term to find.
!      Essentially these are second order terms associated with integrals over k with saddles.
!      The contribution at the saddle is the first order correction and included in the new shortened sum.
!      Most second order terms are tiny in comparison with their first order counterparts and can be ignored.
!      But those close to the endpoints k=0,1,2  and k=L(nit),L(nit)-1,... are close to being the same size
!      as the first order terms (because the saddle is situated very close to an endpoint of the integration
!      range. So at both endpoints these second order tems have computed and sums. The value of ip sets how
!      many are summed at each endpoint. Unfortunately convergence is slow so potentially ip can be large.
!      I have got quite good results for ip=3 or 4.


       sum0=(0.0,0.0)
       jbot=ceiling(phicoeff(1,nit))
       
       con1=jbot-phicoeff(1,nit)
       call psi(con1,ps1)
       sav2=sav2-ps1
       con2=jbot-(L(nit)+fracL(nit))
       call psi(con2,ps2)
       if (fracL(nit)/(2*sqrt(abs(gcddashmt))).gt.0.1) then
        sav1=ps2-sav1
       else
        con1=-fracL(nit)
        call psi(con1,ps1)
        sav1=ps2-ps1
       endif
       t5=(0.0,1.0)*(sav2+endpoint*sav1)/tpp

       ip1=ip
       if (ip1.gt.(0.5*L(nit))) then
        ip1=max1(0.2*L(nit),1.0)
       endif
       do i=jbot,ip1

!     Compute the value of the saddle point c. This has an approximate value
!
!      (i-phi1)/(2*phi2)-Sum[n=3 to nmax] n*((i-phi1)^(n-1))*phin/(2*phi2)^n +smaller terms
!
!     Unfortunately for large i the smaller terms, whilst small are not negligible. And they become
!     inceasingly complex. It is possible to show that whilst in theory nmax=infinity one can cut if
!     at a value nmax<=(mc+1-q)/l+2 at which point any higher order terms give a negligible contribution
!     to the value of g(c).
!     However, in turns out that as the length of the sum reduces g(c) becomes more complex and
!     the value of g(c appropriate for I=L(nit)) goes from third order, to 4th order, to fifth order etc.
!     This is determined by the paramter small found in the routine scheme. So as the scheme complexity increases
!     so does the complexity associated with c.
!     Rather than try and take this explcitly into account (as HAS to be done in routine scheme), here it is sufficient to
!     compute c numerically to some specified tolerance. So for schemes beyond m(nit)>5 a Newton-Raphson method is
!     used to compute c as in computenewton. So first a guess for c is computed. For schemes m(nit)=3,4 or 5 this is sufficiently accurate.
!     For higher order schemes the guess is refined using Newton-Raphson.

        con1=real(i,dp)-phicoeff(1,nit)
        wm=3*phicoeff(3,nit)*(con1**2)/(xr(nit)**3)
        if (m1(nit).ge.4) then
         wm=wm+(4*phicoeff(4,nit)-18*(phicoeff(3,nit)**2)/xr(nit))*(con1**3)/(xr(nit)**4)
        endif
        if (m1(nit).ge.5) then
         wm=wm+(5*(xr(nit)**2)*phicoeff(5,nit)+135*(phicoeff(3,nit)**3)-60*xr(nit)*phicoeff(3,nit) &
          *phicoeff(4,nit))*(con1**4)/(xr(nit)**7)
        endif
        c=con1/xr(nit)-wm
        if (m1(nit).gt.5) then
         numb=1
         call computecnewton(m1(nit),nit,con1,tol,c,c2,numb)
        endif
        sav1=sp*sx*c*sqrt(2.0)
        sav1=erfc(sav1)
        gc=i*c
        do j=1,m1(nit)
         gc=gc-phicoeff(j,nit)*(c**j)
        enddo 
        cr1=(0.0,1.0)*tpp*gc
        cr1=exp(cr1)
        sav2=exp(tpm*con1*c)
        cr2=(0.0,1.0)*phicoeff(2,nit)/p
        cr3=(1/con1+cr2*(1+tpp*c*con1*(1.0+p*con1*c))/(con1**3))
        sum0=sum0+(-cr1*sav1/(2*sx*EPI4)-(0.0,1.0)*(sav2*cr3-cr2/(con1**3))/tpp)
       enddo


!      These are the ip terms for k=possibly 0,1,2,..ip for the lower endpoint. Basically there are two contributions
!      sum and sum2. The integral is split into I4+I5+I6. I6=0, I5 contains the saddle. Sum is the I4 terms.
!      sum2 gives the second order terms in I5 (the main terms being the saddles used in the next Gaussian sum
!      iteration). The terms making up sum and sum2 are comparable in magnitude. The main saddle terms are
!      considerably larger.

       sum1=(0.0,0.0)
       c1=L(nit)+fracL(nit)

       ecor=0.0
       
       do i=L(nit),L(nit)-ip1+1,-1
        if (Lnosub.gt.0.and.i.gt.Lnosub) then
         numb=2
        else
         numb=1
        endif
        con1=real(i,dp)-phicoeff(1,nit)
        wm=3*phicoeff(3,nit)*(con1**2)/(xr(nit)**3)
        if (m1(nit).ge.4) then
         wm=wm+(4*phicoeff(4,nit)-18*(phicoeff(3,nit)**2)/xr(nit))*(con1**3)/(xr(nit)**4)
        endif
        if (m1(nit).ge.5) then
         wm=wm+(5*(xr(nit)**2)*phicoeff(5,nit)+135*(phicoeff(3,nit)**3)-60*xr(nit)*phicoeff(3,nit) &
          *phicoeff(4,nit))*(con1**4)/(xr(nit)**7)
        endif
        c=con1/xr(nit)-wm
        if (m1(nit).gt.3) then
         call computecnewton(m1(nit),nit,con1,tol,c,c2,numb)
        endif
        if (numb.eq.2) then
         c=c2
        endif
        beta=c/real(L(nit-1),dp)
        xbeta=L(nit-1)*(1-beta)
        con2=sqrt(abs(xr(nit)/2.0+3*phicoeff(3,nit)*c))
        sav1=2*sp*con2*xbeta
        sav1=erfc(sav1)
        gc=i*c
        do j=1,m1(nit)
         gc=gc-phicoeff(j,nit)*(c**j)
        enddo
        if ((c1-i)/(2*sqrt(abs(gcddashmt))).gt.0.1) then
         cr1=(0.0,1.0)*tpp*gc
         cr1=exp(cr1)
         cr2=(0.0,1.0)*(gcddashmt)/p
         sav2=exp(tpm*xbeta*(c1-i))
         cr3=(1/(c1-i)+cr2*(1+tpp*xbeta*(c1-i)*(1.0+p*(c1-i)*xbeta))/((c1-i)**3))
         if (numb.eq.1) then
          sum1=sum1+(-sav1)*cr1/(2*(sqrt(2.0)*con2)*EPI4)
         else
          sum1=sum1+(-sav1)*cr1*EPI4/(2*(sqrt(2.0)*con2))
         endif
         sum1=sum1-((0.0,1.0)*endpoint*(sav2*cr3-cr2/((c1-i)**3))/tpp)
        else
         con1=sp*(c1-i)/sqrt(2*abs(gcddashmt))
         call erf(con1,6,cr2)
         cr2=conjg(cr2)
         cr3=epi4*exp((0.0,1.0)*con1*con1)*(1.0-cr2)/sqrt(8*abs(gcddashmt))
         if (gcddashmt.lt.0.0) then
          cr3=conjg(cr3)
         endif
         sum1=sum1+(0.0,1.0)*endpoint*cr3
        endif      
       enddo
     
       t5=t5+(SUM0+SUM1)

!      These are the ip terms for k=L(nit),L(nit)-1 for the upper endpoint. Basically there are two contributions
!      sum1 and sum3. The integral is split into I1+I2+I3. I1=0, I2 contains the saddle. Sum1 is the I3 terms.
!      sum3 gives the second order terms in I2 (the main terms being the saddles used in the next Gaussian sum
!      iteration). The terms making up sum1 and sum3 are comparable in magnitude. The main saddle terms are
!      considerably larger. Note the exp(I*ecor) correction. One expectS this to be very close to one but for these large
!      k values it about exp(I*0.05) which means it is worth including. All these terms make up t5.


       qq=conjg(t1+t2+t4)+t3+t5

      END subroutine q             

      !     Routine used to set up a Gauss-Legendre numerical integration scheme. Only used in routine Start
      !     and only very occasionally when sqrt(8*t/pi) approximates an odd integer.
      
       SUBROUTINE gauleg(x1,x2,ab,we)

       implicit none

       integer, parameter :: dp = selected_real_kind(33)
       integer, parameter :: dp1 = selected_int_kind(16)
       
       real(dp) :: x1,x2,ab(64),we(64),eps,p1,p2,p3,pp,xl,xm,z,z1,p,sp,tpp,tpm
       integer :: i,j,m,n
       complex(kind=16) :: EPI4

       COMMON/PARS4/ P,SP,TPP,TPM,EPI4

       eps=3.0E-14
       n=64

       m=(n+1)/2
       xm=0.5*(x2+x1)
       xl=0.5*(x2-x1)

       do i=1,m
        z=cos(p*(i-0.25)/(n+0.5))

756     continue

        p1=1.0
        p2=0.0
        do j=1,n
         p3=p2
         p2=p1
         p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/(j*1.0)
        enddo

        pp=n*(z*p1-p2)/(z*z-1.0)
        z1=z
        z=z1-p1/pp
        if (abs(z-z1).gt.eps) goto 756
        ab(i)=xm-xl*z
        ab(n+1-i)=xm+xl*z
        we(i)=2*xl/((1.0-z*z)*pp*pp)
        we(n+1-i)=we(i)
       enddo
       return
       end subroutine gauleg

! PARI Library routines       

subroutine pari_calc(rzsum,iteration,rsphase,N,NIS)





  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none


  integer, parameter   :: dp = selected_real_kind(33)
  integer(kind=C_LONG) :: prec   = 5 ! 152 decimal digits

  type(C_PTR)          :: t, rsphase_pari
  type(C_PTR)          :: u,v,w

  real(dp), INTENT(OUT):: rzsum
  real(dp), INTENT(IN) :: rsphase

  integer(kind=C_LONG) :: I
  integer(kind=C_LONG) :: av
  integer*8 :: I_START, I_END, iteration, N,  NIS



  CALL pari_init(10000000_8, 2_8)


  av = get_avma()

  CALL t_grabber(t)
  CALL real2pari(rsphase,rsphase_pari)


  if (NIS.eq.0) then
     I_START = 1 + (iteration-1) * N
  else
     I_START = NIS
  endif   
  I_END   =         iteration * N

  w = stor(0_8,prec)

  DO I=I_START,I_END



     u =  gmul(t, glog( stor(I,prec) ,prec )  )  
     u =  gsub(u , rsphase_pari )
     u =  gcos( u, prec         )
     v =  gsqrt( stor(I,prec) , prec)
     u =  gdiv(u,v)
     w =  gadd( u , w )


     if( modulo(I,1000) .EQ. 0) then

        w = gerepilecopy(av,w)
        CALL t_grabber(t)
        CALL real2pari(rsphase,rsphase_pari)

     endif

  ENDDO

  rzsum = rzsum + rtodbl(w)


  CALL pari_close()

END subroutine pari_calc




subroutine t_grabber(t)


  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none

  integer(kind=C_LONG) :: prec   = 5 ! 152 decimal digits

  type(C_PTR)            :: u1, u2, u3, u4, u5, u6
  type(C_PTR)            :: t

  integer(kind=C_LONG)   ::  c1, c2, c3, c4, c5, n 
  CHARACTER(LEN=45)      ::   Y

  COMMON/PARS9/ Y

  !       CALL pari_init(10000000_8, 2_8)

!  Y = '4.6709141854660972368505489032740000000000e28'  ! keep decimal digits 40 long
  ! Y = '3.0123456789987654321000000000000000000000e32'  ! example input.


  READ(Y(1:1)  , '(I1)'  )  c1            ! first digit and 10-digit pieces of t
  READ(Y(3:12)  , '(I10)'  )  c2
  READ(Y(13:22)  , '(I10)'  )  c3
  READ(Y(23:32)  , '(I10)'  )  c4
  READ(Y(33:42)  , '(I10)'  )  c5
  READ(Y(44:45)  , '(I10)'  )  n


  u1 = gpow( stoi(10_8), stoi(  0_8 ), prec)    ! construct t in blocks of 10 digits
  u2 = gpow( stoi(10_8), stoi(  -10_8 ), prec)
  u3 = gpow( stoi(10_8), stoi(  -20_8 ), prec)
  u4 = gpow( stoi(10_8), stoi(  -30_8 ), prec)
  u5 = gpow( stoi(10_8), stoi(  -40_8 ), prec)
  u6 = gpow( stoi(10_8), stoi(  n ), prec)


  t = gadd( gadd( gadd( gadd(                           &
       gmul(stor(c1,prec) , u1)        ,           &
       gmul(stor(c2,prec) , u2)      ) ,           &
       gmul(stor(c3,prec) , u3)      ) ,           &
       gmul(stor(c4,prec) , u4)      ) ,           &
       gmul(stor(c5,prec) , u5)      )
  
  !  constructed high precision t value.

  t = gmul(t,u6)


END subroutine t_grabber




subroutine real2pari(r,r_pari)


  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none

  integer, parameter   :: dp = selected_real_kind(33)
  integer(kind=C_LONG) :: prec   = 5 ! 152 decimal digits

  type(C_PTR)            :: u0, u1, u2, u3, u4, u5
  type(C_PTR)            :: r_pari

  integer(kind=C_LONG)   :: c0, c1, c2, c3, c4, c5
  CHARACTER(LEN=70)      ::   Y

  real(dp)               ::  r


     WRITE(Y, '(F68.45)') r


    READ(Y(3:12)   , '(I10)'  )  c0
    READ(Y(13:22)  , '(I10)'  )  c1            ! first digit and 10-digit pieces of r
    READ(Y(24:33)  , '(I10)'  )  c2
    READ(Y(34:43)  , '(I10)'  )  c3
    READ(Y(44:53)  , '(I10)'  )  c4
    READ(Y(54:63)  , '(I10)'  )  c5



  u0 = gpow( stoi(10_8), stoi(   10_8 ), prec)
  u1 = gpow( stoi(10_8), stoi(    0_8 ), prec)    ! construct t in blocks of 10 digits
  u2 = gpow( stoi(10_8), stoi(  -10_8 ), prec)
  u3 = gpow( stoi(10_8), stoi(  -20_8 ), prec)
  u4 = gpow( stoi(10_8), stoi(  -30_8 ), prec)
  u5 = gpow( stoi(10_8), stoi(  -40_8 ), prec)


    r_pari = gadd( gadd( gadd( gadd( gadd(         &
       gmul(stor(c0,prec) , u0)        ,           &
       gmul(stor(c1,prec) , u1)      ) ,           &
       gmul(stor(c2,prec) , u2)      ) ,           &
       gmul(stor(c3,prec) , u3)      ) ,           &
       gmul(stor(c4,prec) , u4)      ) ,           &
       gmul(stor(c5,prec) , u5)      )
    
  !  constructed high precision t value.


END subroutine real2pari





subroutine pari_erfc(z,rz,iz)





  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none


  integer, parameter   :: dp = selected_real_kind(33)
  integer(kind=C_LONG) :: prec   = 5 ! 152 decimal digits

  type(C_PTR)          :: z_pari
  type(C_PTR)          :: u, S

  real(dp)             :: z, rz, iz



  CALL pari_init(10000000_8, 2_8)

  CALL real2pari(z,z_pari)   ! I THINK THIS IS A PLACE AN ERROR MIGHT POP UP


  S = gatan( stor(1_8 , prec) ,prec)        ! PI/4
  u = gsqrt( stoi(-1_8) ,prec)              ! imag(1)
  S = gmul ( S, u      )                     ! I*pi/4
  S = gmul ( S, z_pari )                     ! Z*I*pi/4
  S = gexp ( S, prec   )                     ! exp( Z*I*PI/4)
  S = gerfc( S, prec   )                     ! erfc( above )

  rz = rtodbl( greal(S) )
  iz = rtodbl( gimag(S) )


  CALL pari_close()

END subroutine pari_erfc





subroutine pari_tse(rae,tse)

        ! Routine to calculate:   tse = (t*SE) mod 2pi

  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none


  integer, parameter   :: dp = selected_real_kind(33)
  integer(kind=C_LONG) :: prec   = 5 ! 152 decimal digits change 5
  integer(kind=C_LONG) :: av         ! change 6

  type(C_PTR)          :: t, rae_pari, SE_pari
  type(C_PTR)          :: u, S, e, pc, v

  real(dp), INTENT(IN) :: rae
  real(dp), INTENT(OUT):: tse
  real(dp)             :: rae1
  real                 :: frac_part

  
  av = get_avma()    !change 7
  
  CALL real2pari(rae,rae_pari)   
  CALL t_grabber(t)



  frac_part = rtodbl( gfrac(rae_pari) )

  ! Manual rounding of rae

  if(frac_part.ge.0.5) then

          rae1 = ceiling(rae,dp)

  else

          rae1 = floor(rae,dp)

  endif


  CALL real2pari(rae1,rae_pari)


  S = stor(8_8,prec)
  u = Pi2n(1_8, prec)
  u = gdiv(u, stor(2_8,prec) )
  S = gmul(t,S)
  S = gdiv(u,S)


  S = gsqrt(S,prec)    ! We are up to S = 1/ sqrt( t*8 / pi ) - which is ain

  e = gmul(rae_pari, S) ! this is 'e'


  u = gmul(e,e)
  u = gmul( stor(2_8,prec) , u )

  v = gmul(e,e)
  v = gdiv( stor(1_8,prec), v )
  v = gsub(stor(1_8,prec),v)
  v = gsqrt(v,prec)
  v = gadd( stor(1_8,prec) , v )

  S = gmul(u,v)
  pc = gsub(S , stor(1_8,prec) )


  u = glog(pc,prec)


  v = gdiv( stor(1_8,prec), pc)
  S = gadd( u , v )
  S = gadd( S , stor(1_8,prec) )
  SE_pari = gdiv( S , stor(2_8,prec) )

  S = gmul(t,SE_pari)
  S = gmod(S, Pi2n(1_8, prec))

  tse = rtodbl(S)

  CALL set_avma(av)   ! change 8


END subroutine pari_tse


subroutine pari_phases(yphase,rsphase,a)

        ! Routine to calculate:   yphase     = (t+pi/8) mod 2pi
        !                         theta(t)   = IM( log(Gamma(1/4+i*t/2))  -t*log(pi)/2 )
        !                         rsphase(t) = theta(t) mod 2pi

  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none



  integer, parameter   :: dp = selected_real_kind(33)
  integer(kind=C_LONG) :: prec   = 25 ! 152 decimal digits


  type(C_PTR)          :: t,u,v,w,x,y
  type(C_PTR)          :: yphase_pari, theta_pari, rsphase_pari, a_pari

  real(dp), INTENT(OUT):: yphase, rsphase, a

  yphase = 0.0
  rsphase= 0.0


  CALL pari_init(10000000_8, 2_8)


    CALL t_grabber(t)


    v           = Pi2n(1_8, prec) ! 2pi
    w           = gsqrt( stoi(-1_8) ,prec) ! imag(1)
    y           = gdiv(v, stor(2_8,prec)) ! pi
    u           = gdiv(y, stor(8_8,prec)) ! pi/8
    u           = gadd(t,u)               ! t + pi/8
    yphase_pari = gmod(u, v)              ! t + pi/8 mod 2pi

    CALL rtoquad(yphase_pari,yphase)
    
    u           = gdiv(t, stor(2_8,prec))                       ! t/2
    u           = gmul(w, u)                                    ! i*t/2
    x           = gdiv( stor(1_8,prec), stor(4_8,prec) )        ! 1/4
    u           = gadd(u,x)                                     ! 1/4 + i*t/2
    u           = glngamma(u,prec)                              ! log(gamma( 1/4 + i*t/2 ))
    u           = gimag( u )                                    ! imag( above )

    x           = glog(y,prec)                                  ! log(pi)
    x           = gdiv(x, stor(-2_8,prec) )                     ! -log(pi)/2
    x           = gmul(x,t)                                     ! -t*log(pi)/2
    theta_pari  = gadd(u,x)                                     ! imag(log(gamma( 1/4 + i*t/2 ))) - t*log(pi)/2
    rsphase_pari = gmod(theta_pari, v)                          ! rsphase = (above) mod 2pi

    CALL rtoquad(rsphase_pari,rsphase)

    u           = gdiv(t, y)                                    ! t/pi
    v           = gmul(u, stor(8_8,prec))                       ! 8*t/pi
    a_pari      = gsqrt(v ,prec)                                ! sqrt(8*t/pi)

    CALL rtoquad(a_pari,a)
    
  CALL pari_close()

END subroutine pari_phases

subroutine rtoquad(pari_input,var)


  use ISO_C_BINDING, only : C_PTR, C_DOUBLE
  use PARI
  implicit none



  integer, parameter   :: dp = selected_real_kind(33)
  integer, parameter   :: dp1 = selected_int_kind(16)
  integer(kind=C_LONG) :: prec   = 25 ! 152 decimal digits
  integer(kind=C_LONG) :: m

  type(C_PTR)          :: n,x,y,y1,y2,pari_input
  type(C_PTR)          :: expo, scale_d, scale_u

  real(dp)             :: var, var1

  integer(dp1)         :: num 

  num = 10000000
 
!  CALL pari_init(10000000_8, 2_8)



     y           = pari_input


    m           = sizedigit(  y  ) -1

    if(m .le. 0) then
      m =  -sizedigit( gdiv( stoi(1_8), y  ))
    endif


    n           =  stoi( m )     !gsub(sizedigit(y), stoi(1_8) ) ! n-1    (for 1000, n = 4 as there are 4 digits)

 

    scale_d     = gpow( stoi(10_8), gmul( stor(-1_8,prec), n) , prec)
    scale_u     = gpow( stoi(10_8), n, prec)
    expo        = gsqr( stor(10000000_8, prec) ) ! 10^14




    y           = gmul(y,scale_d)

    y2          = gmul(y, expo )  !expo*pi
    x           = gfloor(y2)      ! integer part (expo * pi)
    y2          = gsub(y2, x)     ! decimal part (expo * pi)

    y1          = gdiv(x,  expo ) ! first part of pi
    y2          = gdiv(y2, expo ) ! second part of pi

    var   = 0.0
    var1  = 0.0

    var = rtodbl(y1)
    var = var * num
    var = var * num
    
    var = nint(var,dp1)
    var1= rtodbl(y2)
    var1= var1* num
    var1= var1* num
    var = var + var1
    var = var / num
    var = var / num



    if(m>0) then    
    var = var *  10**m
    else   
    var = var / 10**(abs(m))
    endif

 END subroutine rtoquad 
