!*********************************************************************************************
!
 MODULE AtomDiatomskvp
!
!*********************************************************************************************
!
! This module contains variable description for the SKVP Atom-Diatom program
! 
!=============================================================================================
!
!
!
 IMPLICIT NONE
!
 SAVE
!
!
 integer, parameter :: dp = kind(1.0d0)
 REAL(8)                                        :: alpha=5.0E-01, test_int, E, pi = dacos(-1d0)
 INTEGER                                        :: p=2, istatus, dim_x
 INTEGER                                        :: ngqp_x, ngqp_y, ncf, n_open, Jtot !=0
 INTEGER                                        :: n_xmb, k_xmb, jmxmb, k_ymb, n_openmb
!
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)        :: BAM_1, BAM_2, BAM_3, BAM_4, BAM_x1, BAM_x2, BAM_x3, BAM_x4, BAM_V, M10_V, M00_V
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)        :: BAM_xx1, BAM_xx2, BAM_xx3, BAM_xx4, mat_B, mat_C, Smat, BAM_r101, BAM_r102
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)        :: BAM_r103, BAM_xx1_test
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)        :: BAM_xx1b, BAM_xx2b, BAM_xx3b, BAM_xx4b, BAM_r001, BAM_r002, BAM_r003, M00_V2
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)        :: mat_M0, mat_M10, mat_M00, M0_V, M10_V2, M_V2, BAM_r01, BAM_r02, BAM_r03
 REAL(8), ALLOCATABLE, DIMENSION(:)             :: gq_weight_x, gq_root_x, norm, knots_x, norm_vec
 REAL(8), ALLOCATABLE, DIMENSION(:)             :: gq_weight_y, gq_root_y, eval1d, kvec, x, wx
 REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: M_V
 REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: BAM_r1, BAM_r2, BAM_r3
 REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: BAM_theta0, BAM_theta1, BAM_theta2
 REAL(8), ALLOCATABLE, DIMENSION(:,:,:)         :: BAM_vterm
 REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: BAM_theta00, BAM_theta01, BAM_theta02, BAM_theta000, BAM_theta001
 REAL(8), ALLOCATABLE, DIMENSION(:,:)           :: BAM_theta002, BAM_theta100, BAM_theta101, BAM_theta102, M_K, mat_M
 INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: quant_mat
 LOGICAL                                        :: bsp_optimized = .TRUE.

 REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: BAM_r
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:) :: BAM_r0
 REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: BAM_theta
 COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:) :: BAM_r00, BAM_r10
 INTEGER, allocatable, dimension(:) :: open_idx




!
!
 CONTAINS
!========
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!  LOGICAL FUNCTION iseven(m)
! !
! !*********************************************************************************************
! !
! !	This function determines whether m is even or odd
! !
! !=============================================================================================
! !
! 	INTEGER m
! !
! 	IF (m==0) THEN
! 		iseven = .TRUE.
! 	ELSE	IF (mod(m,2)==0) THEN
! 				iseven = .TRUE.
! 			ELSE
! 				iseven = .FALSE.
! 	ENDIF
! !
! !=============================================================================================
!  END FUNCTION iseven
!=============================================================================================





!***************************************************************


!*********************************************************************************************
!
 REAL(8) FUNCTION delta(i,j)
!
!*********************************************************************************************
!
!       This function computes the delta function 
!
!=============================================================================================
!
        INTEGER i,j
!
        delta=0d0
        IF (i==j) THEN
           delta=1d0
        ENDIF
!
!=============================================================================================
 END FUNCTION delta
!=============================================================================================
!
 
!*********************************************************************************************
!*********************************************************************************************
!
!  REAL(8) FUNCTION deltar(i,j)
! !
! !*********************************************************************************************
! !
! !       This function computes the delta function for reals
! !
! !=============================================================================================
! !
!         REAL(8) :: i,j
! !
!         deltar=0d0
!         IF (dabs(i-j)<1d-16) THEN
!            deltar=1d0
!         ENDIF
! !
! !=============================================================================================
!  END FUNCTION deltar
!=============================================================================================
!
!*********************************************************************************************
!*********************************************************************************************
!
 REAL(8) FUNCTION factorial(n)
!
!*********************************************************************************************
!
!       Factorial function of an integer
!
!=============================================================================================
!
        INTEGER i, n
!
        factorial=1d0
        do i=2,n
           factorial=factorial*i
        enddo
!
!=============================================================================================
 END FUNCTION factorial
!=============================================================================================
!
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!                                                                                             
!                                                               B-Splines functions          
!
!!*********************************************************************************************
!!*********************************************************************************************
!!
! SUBROUTINE searchleft_x(xx,lefte)
!!
!!               This function search the left value useful for bsp and dbsp
!!
!!=============================================================================================
!!
!        USE generateparameters
!!
!        IMPLICIT NONE
!!--------------------
!        INTEGER lefte
!        REAL(8) xx
!!
!        lefte=pbasst(1)%pb_pa1
!        DO WHILE (knots_x(lefte+1)<xx)
!                lefte=lefte+1
!        ENDDO
!!
!        IF (xx>knots_x(pbasst(1)%pb_nbr-1+pbasst(1)%pb_pa1)) THEN
!                lefte = pbasst(1)%pb_nbr-1+pbasst(1)%pb_pa1+1
!        ELSE IF (xx<knots_x(1)) THEN
!                lefte = 0
!        ENDIF
!!=============================================================================================
! END SUBROUTINE searchleft_x
!!=============================================================================================
!!
!!
!!*********************************************************************************************
!!*********************************************************************************************
!!
! SUBROUTINE searchleft_y(xx,lefte)
!!
!!               This function search the left value useful for bsp and dbsp
!!
!!=============================================================================================
!!
!        IMPLICIT NONE
!!--------------------
!        INTEGER lefte
!        REAL(8) xx
!!
!        lefte=k_y
!        DO WHILE (knots_y(lefte+1)<xx)
!                lefte=lefte+1
!        ENDDO
!!
!        IF (xx>knots_y(n_y-1+k_y)) THEN
!                lefte = n_y-1+k_y+1
!        ELSE IF (xx<knots_y(1)) THEN
!                lefte = 0
!        ENDIF
!!=============================================================================================
! END SUBROUTINE searchleft_y
!!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!
 REAL(8) FUNCTION bsp_x(j,xx)
!
!*********************************************************************************************
!
!                       This function computes B_j(x) with x knot sequence, order k  
!                        and emax values obtained from the module global_definition
!  
!=============================================================================================
!
        USE generateparameters
        USE bspline_mod
!
        IMPLICIT NONE
!--------------------
        INTEGER j
        REAL(8) xx
        !REAL(8) xx, biatx(k_x)
        !bsp_x=b_spline(j,k_x,xx,knots_x)
        bsp_x=b_spline(j,pbasst(1)%pb_pa1-1,xx,knots_x)
        return 
!!
!        call searchleft_x(xx,left)
!        bsp_x=0d0
!        IF ((left-k_x+1.le.j).and.(j.le.left)) THEN
!                call bsplvbx(knots_x,k_x,1,xx,left,biatx)
!                bsp_x=biatx(j-left+k_x)
!        ENDIF
!=============================================================================================
 END FUNCTION bsp_x
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!  REAL(8) FUNCTION bsp_y(j,x)
! !
! !*********************************************************************************************
! !
! !                       This function computes B_j(x) with x knot sequence, order k  
! !                        and emax values obtained from the module global_definition
! !  
! !=============================================================================================
! !
! !       USE NUMERICAL_LIBRARIES
! !
!         IMPLICIT NONE
! !--------------------
!         INTEGER left, j
!         REAL(8) x, pas, biaty(k_y)
! !
!         call searchleft_y(x,left)
!         bsp_y=0d0
!         IF ((left-k_y+1.le.j).and.(j.le.left)) THEN
!                 call bsplvby(knots_y,k_y,1,x,left,biaty)
!                 bsp_y=biaty(j-left+k_y)
!         ENDIF
! !=============================================================================================
!  END FUNCTION bsp_y
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!
 REAL(8) FUNCTION dbsp_x(j,xx)
!
!*********************************************************************************************
!
!                       This function computes B_j'(x) with  x knot sequence, order k  
!                        and emax values obtained from the module global_definition
!  
!=============================================================================================
!
        USE generateparameters    
        USE bspline_mod
!
        IMPLICIT NONE
!--------------------
        INTEGER j
        REAL(8) xx
        !REAL(8) xx, dbiatx(k_x,2), a(k_x,k_x)
        !dbsp_x=b_spline_deriv(j,k_x,xx,knots_x,1)
        dbsp_x=b_spline_deriv(j,pbasst(1)%pb_pa1-1,xx,knots_x,1)
        return 
!!
!        call searchleft_x(xx,left)
!        dbsp_x=0d0
!        IF ((left-k_x+1.le.j).and.(j.le.left)) THEN
!                call bsplvdx(knots_x,k_x,xx,left,a,dbiatx,2)
!                dbsp_x=dbiatx(j-left+k_x,2)
!        ENDIF
!=============================================================================================
 END FUNCTION dbsp_x
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!  REAL(8) FUNCTION dbsp_y(j,x)
! !
! !*********************************************************************************************
! !
! !                       This function computes B_j'(x) with  x knot sequence, order k  
! !                        and emax values obtained from the module global_definition
! !  
! !=============================================================================================
! !
! !       USE NUMERICAL_LIBRARIES
! !
!         IMPLICIT NONE
! !----------------
!         INTEGER left, j
!         REAL(8) x, pas, dbiaty(k_y,2), a(k_y,k_y)
! !
!         call searchleft_y(x,left)
!         dbsp_y=0d0
!         IF ((left-k_y+1.le.j).and.(j.le.left)) THEN
!                 call bsplvdy(knots_y,k_y,x,left,a,dbiaty,2)
!                 dbsp_y=dbiaty(j-left+k_y,2)
!         ENDIF
! !=============================================================================================
!  END FUNCTION dbsp_y
!=============================================================================================
!
!
!=============================================================================================
!
!							Integration Routines
!
!*********************************************************************************************
!*********************************************************************************************
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!   SUBROUTINE qgauss(func,a,b,s,w,x,n,i,k)
! ! !
! ! !*********************************************************************************************
! ! !
! ! !					Subroutine that performs gaussian quadrature
! ! !
! ! !*********************************************************************************************
! ! !	
!  	IMPLICIT NONE
! ! ----------------
!  	REAL(8) a, b, s
!  	INTEGER i, k, n, nn
!  	REAL(8) dx, xm, xr, w(n), x(n), func
!  	INTEGER j
! ! 
!  	nn = n/2
!  	xm = .5*(b+a)
!  	xr = .5*(b-a)
!  	s = 0.
! 	DO j = 1, nn
! 		dx = xr*x(j+nn)
! 		s = s + w(j+nn) * (func(i,k,xm+dx) + func(i,k,xm-dx))
! 	ENDDO
! 	s = xr * s
! 	RETURN
! !
! !=============================================================================================
!  END SUBROUTINE qgauss
!  (C) Copr. 1986-92 Numerical Recipes Software &.
!=============================================================================================
! 	
!
!*********************************************************************************************
!*********************************************************************************************
!
 SUBROUTINE gauleg(x1,x2,x,w,n)
!
!*********************************************************************************************
!
!	This routine gives the weights and roots for Gauss-Legendre integration
!
!=============================================================================================
!
    INTEGER n
    REAL(8) x1,x2,x(n),w(n)
    DOUBLE PRECISION EPS
    PARAMETER (EPS=3.d-14)
    INTEGER i,j,m
    DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO 12 i=1,m
        z=DCOS(pi*(i-.25d0)/(n+.5d0))
1       CONTINUE
          p1=1.d0
          p2=0.d0
          DO 11 j=1,n
                p3=p2
                p2=p1
                p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        CONTINUE
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        IF(DABS(z-z1).gt.EPS)GOTO 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12  CONTINUE
    RETURN
!
!=============================================================================================
 END SUBROUTINE gauleg
!  (C) Copr. 1986-92 Numerical Recipes Software &.
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!  SUBROUTINE piksr2(n,arr,brr)
! !
! !*********************************************************************************************
! !
! !       This routine sorts simultaneaously the arrays arr and brr
! !
! !=============================================================================================
! !
!     INTEGER n
!     REAL(8) arr(n),brr(n)
!     INTEGER i,j
!     REAL(8) a,b
! !
!     DO 12 j = 2, n
!        a=arr(j)
!        b=brr(j)
!        DO 11 i = j-1,1,-1
!              IF(arr(i).le.a)GOTO 10
!              arr(i+1)=arr(i)
!              brr(i+1)=brr(i)
! 11     ENDDO
!        i=0
! 10     arr(i+1)=a
!        brr(i+1)=b
! 12  ENDDO 
!     RETURN
! !
! !=============================================================================================
!  END SUBROUTINE piksr2
!  (C) Copr. 1986-92 Numerical Recipes Software &.
!=============================================================================================
!
!
!*********************************************************************************************
!
! Function for Legendre Polynomials
!
!*********************************************************************************************

REAL(8) FUNCTION plgndr(l,m,x)
        INTEGER l,m,mm
        REAL(8) x
        INTEGER i,ll
        REAL(8) fact,pll,pmm,pmmp1,somx2
        LOGICAL :: is_negative

        is_negative = .FALSE.
        !IF(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.)PAUSE

        mm = m
        IF (mm .lt. 0) THEN
           mm = -mm
           is_negative = .TRUE.
        ENDIF

        pmm=1.
        IF (mm .gt. 0) THEN
           somx2=dsqrt((1d0-x)*(1d0+x))
           fact=1d0
           DO i=1,m
              pmm=-pmm*fact*somx2
              fact=fact+2d0
           END DO
        END IF
        IF (l .eq. m) THEN
           plgndr=pmm
        ELSE
           pmmp1=x*(2*m+1)*pmm
           IF (l.eq.m+1) THEN
              plgndr=pmmp1
           ELSE
              DO ll=m+2,l
                 pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                 pmm=pmmp1
                 pmmp1=pll
              END DO
              plgndr=pll
           END IF
        END IF
        RETURN

END FUNCTION plgndr
!   (C) Copr. 1986-92 Numerical Recipes Software &.

!*********************************************************************************************
!
! Function for normalized u_l functions
!
!*********************************************************************************************

! REAL(8) FUNCTION u_function(alpha, l, r)
!         REAL(8) :: alpha, l, r
!         REAL(8) :: norm_const, integral
!         REAL(8), ALLOCATABLE, DIMENSION(:) :: x,w
!         INTEGER :: i, n
!         n = ngqp_x
!         ALLOCATE(x(n), w(n))

!         ! Normalizing functions using gauss quad
!         CALL gauleg(0d0, 10d0, x, w, n)
!         integral = 0d0
!         DO i = 1, n
!                 integral = integral + w(i) * x(i)**(2d0*l - 2d0) * DEXP(-2.0d0 * alpha * x(i))
!         END DO

!         norm_const = 1d0/(DSQRT(integral))

!         u_function = norm_const * r**(l-1) * DEXP(-1d0 * alpha * r)


! END FUNCTION u_function

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!                                                                                            
!                                                               Phi functions          
!
!*********************************************************************************************
!*********************************************************************************************
!
!  REAL(8) FUNCTION phi_y(i,x)
! !
! !               This function evaluates the vibrational function phi
! !
! !=============================================================================================
! !
!         IMPLICIT NONE
! !--------------------
!         INTEGER i, j
!         REAL(8) x
! !
! !
!         phi_y=0d0
!         DO j=1, dim_y
!                 phi_y=phi_y+evec1d(j,i)*bsp_y(j+1,x)
!         ENDDO
! !        
! !=============================================================================================
!  END FUNCTION phi_y
!=============================================================================================
!
!
!
!*********************************************************************************************
!*********************************************************************************************
!
!  REAL(8) FUNCTION dphi_y(i,x)
! !
! !               This function evaluates the derivative of the vibrational function phi
! !
! !=============================================================================================
! !
!         IMPLICIT NONE
! !--------------------
!         INTEGER i, j
!         REAL(8) x
! !
! !
!         dphi_y=0d0
!         DO j=1, dim_y
!                 dphi_y=dphi_y+evec1d(j,i)*dbsp_y(j+1,x)
!         ENDDO
! !
! !=============================================================================================
!  END FUNCTION dphi_y
!=============================================================================================
!
!
!

 

REAL(8) FUNCTION lambda_plus(jt,j,k,l)
        INTEGER :: j, k, jt, lm, l
        !IF (1d0*(j*(j+1) - k*(k+1)) < 0) THEN
        !        PRINT*, "Complex Args in lambda_plus: j,k", j,k, 1d0*(j*(j+1) - k*(k+1))
        !ENDIF
        lm=min(j,jt)
        IF (k>=lm) THEN
        lambda_plus = 0d0
        ENDIF
        IF ((k>=-lm).AND.(k<lm)) THEN
        IF (l==1) THEN  ! With Jt
        lambda_plus = DSQRT(1d0*(jt*(jt+1) - k*(k+1)))
        ELSE IF (l==2) THEN  ! With j
        lambda_plus = DSQRT(1d0*(j*(j+1) - k*(k+1)))
        ENDIF
        ENDIF
END FUNCTION lambda_plus

REAL(8) FUNCTION lambda_minus(jt,j,k,l)
        INTEGER :: j, k, jt, lm, l
        !IF (1d0*(j*(j+1) - k*(k-1)) < 0) THEN
        !        PRINT*, "Complex Args in lambda_minus: j,k", j,k, 1d0*(j*(j+1) - k*(k-1))
        !ENDIF
        lm=min(j,jt)
        IF (k<=-lm) THEN
        lambda_minus = 0d0
        ENDIF
        IF ((k>-lm).AND.(k<=lm)) THEN
        IF (l==1) THEN
        lambda_minus = DSQRT(1d0*(jt*(jt+1) - k*(k-1)))
        ELSE IF (l==2) THEN
        lambda_minus = DSQRT(1d0*(j*(j+1) - k*(k-1)))
        ENDIF
        ENDIF
END FUNCTION lambda_minus

REAL(8) FUNCTION W(k,k_prime,jt,j)
        INTEGER :: jt, k, k_prime, j
        !W = 1d0*big_j*(big_j+1d0) 

        W = dble(jt*(jt+1) + j*(j+1) - 2*k**2)*delta(k,k_prime) - &
            lambda_plus(jt,j,k,1)*lambda_plus(jt,j,k,2)*DSQRT(1d0+delta(k,0))*delta(k+1,k_prime) - &
            lambda_minus(jt,j,k,1)*lambda_minus(jt,j,k,2)*DSQRT(1d0+delta(k,1))*delta(k-1,k_prime)
END FUNCTION W

REAL(8) FUNCTION Wdd(Jtot, j1, k1, j2, k2, j1p, k1p, j2p, k2p)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Jtot, j1, k1, j2, k2
        INTEGER, INTENT(IN) :: j1p, k1p, j2p, k2p
        INTEGER :: Ktot
        REAL(8) :: channel_delta

                Ktot = k1 + k2

                channel_delta = delta(j1,j1p) * delta(k1,k1p) * &
                                delta(j2,j2p) * delta(k2,k2p)

                Wdd = 0d0

                Wdd = Wdd + (Jtot*(Jtot+1d0) - 2d0*k1**2 - 2d0*k2**2 - 2d0*k1*k2) * channel_delta

                Wdd = Wdd + lambda_plus(0,j1,k1,2) * lambda_minus(0,j2,k2,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1+1,k1p) * delta(k2-1,k2p)

                Wdd = Wdd + lambda_minus(0,j1,k1,2) * lambda_plus(0,j2,k2,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1-1,k1p) * delta(k2+1,k2p)

                Wdd = Wdd - lambda_plus(Jtot,0,Ktot,1) * lambda_plus(0,j1,k1,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1+1,k1p) * delta(k2,k2p)

                Wdd = Wdd - lambda_plus(Jtot,0,Ktot,1) * lambda_plus(0,j2,k2,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1,k1p) * delta(k2+1,k2p)

                Wdd = Wdd - lambda_minus(Jtot,0,Ktot,1) * lambda_minus(0,j1,k1,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1-1,k1p) * delta(k2,k2p)

                Wdd = Wdd - lambda_minus(Jtot,0,Ktot,1) * lambda_minus(0,j2,k2,2) * &
                        delta(j1,j1p) * delta(j2,j2p) * &
                        delta(k1,k1p) * delta(k2-1,k2p)

END FUNCTION Wdd

! ================================================================================================
! Cut off function for the u0 function
REAL(8) FUNCTION h(R)
        USE generateparameters
        REAL(8) :: R
        h = 0.5d0*(1d0 + DTANH(alpha*(R-r0)))
END FUNCTION h

! u0 function
COMPLEX(8) FUNCTION u0(n,R)
        USE generateparameters
        REAL(8) :: R
        INTEGER :: n, quant_j
        !quant_j = quant_mat(1, n)
        u0 = ZEXP(-COMPLEX(0d0,1d0) * kvec(n) * R)*h(R)*&
        DSQRT(mu_R/kvec(n))        
END FUNCTION u0
!
! ================================================================================================
COMPLEX(8) FUNCTION d2dR2u0(n, R)
        USE generateparameters
        REAL(8) :: R
        INTEGER :: n, quant_j
        d2dR2u0 = DSQRT(mu_R/kvec(n)) * (-(kvec(n))**(2d0) * ZEXP(-(0d0, 1d0) * kvec(n) * R) *h(R) - ZEXP(-(0d0,1d0)*kvec(n)*R)&
        *alpha*((1d0/dcosh(alpha*(R-r0)))**2d0)*((0d0,1d0)*kvec(n) + alpha*dtanh(alpha*(R-r0))))
END FUNCTION d2dR2u0

function wigner3j(l1,l2,l3,m1,m2,m3) result(res)
    integer, intent(in) :: l1,l2,l3,m1,m2,m3
    real(kind=8) :: res
    integer :: k
    integer :: kmin, kmax
    real(kind=8) :: sum, term, pref, sign
    real(kind=8), parameter :: zero=0.0d0

    ! Selection rules
    if (m1 + m2 + m3 /= 0) then
      res = zero; return
    endif
    if (abs(m1) > l1 .or. abs(m2) > l2 .or. abs(m3) > l3) then
      res = zero; return
    endif
    if (l3 < abs(l1 - l2) .or. l3 > l1 + l2) then
      res = zero; return
    endif

    ! Prefactor (with triangle coefficient and phase)
    pref = (-1.0d0)**(l1 - l2 - m3)
    pref = pref * sqrt( &
             dble(factorial(l1 + l2 - l3)) * &
             dble(factorial(l1 - l2 + l3)) * &
             dble(factorial(-l1 + l2 + l3)) / &
             dble(factorial(l1 + l2 + l3 + 1)) )

    pref = pref * sqrt( &
             dble(factorial(l1 + m1)) * dble(factorial(l1 - m1)) * &
             dble(factorial(l2 + m2)) * dble(factorial(l2 - m2)) * &
             dble(factorial(l3 + m3)) * dble(factorial(l3 - m3)) )

    ! Summation limits for summation index kk (not to be confused with kmin/kmax)
    kmin = max(0, l2 - l3 - m1, l1 - l3 + m2)
    kmax = min(l1 + l2 - l3, l1 - m1, l2 + m2)

    sum = 0.0d0
    do k = kmin, kmax
      sign = (-1.0d0)**k
      term = sign / &
             ( dble(factorial(k)) * &
               dble(factorial(l1 + l2 - l3 - k)) * &
               dble(factorial(l1 - m1 - k)) * &
               dble(factorial(l2 + m2 - k)) * &
               dble(factorial(l3 - l2 + m1 + k)) * &
               dble(factorial(l3 - l1 - m2 + k)) )
      sum = sum + term
    end do

    res = pref * sum
  end function wigner3j

  REAL(8) FUNCTION gaunt_coeff(j, k, lam, m, jp, kp)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j, k, lam, m, jp, kp
  REAL(8) :: phase

  gaunt_coeff = 0d0

  IF (-k + m + kp /= 0) RETURN

  !Because we have one conjugate spherical harmonics, there is a phase factor
  phase = MERGE(1d0, -1d0, MOD(ABS(k),2) == 0)

  gaunt_coeff = phase * DSQRT(((2d0*j+1d0)*(2d0*lam+1d0)*(2d0*jp+1d0))/(4d0*pi)) * &
                wigner3j(j,lam,jp,0,0,0) * &
                wigner3j(j,lam,jp,-k,m,kp)
  END FUNCTION gaunt_coeff

function assoc_legendre(l, m, x) result(plm)
        implicit none

        integer, intent(in) :: l, m
        real(dp), intent(in) :: x
        real(dp) :: plm, norm
        real(dp) :: pmm, pmmp1, pll
        real(dp) :: somx2, fact
        integer :: i, ll

        if (m < 0) then
            plm = 0.0_dp
            print *, "Error: m must be non-negative."
            return
        else if (m > l) then
            plm = 0.0_dp
            print *, "Error: m must be less than or equal to l."
            return
        else if (abs(x) > 1.0_dp) then
            plm = 0.0_dp
            print *, "Error: x must be in the range [-1, 1]."
            return
        end if

        norm = sqrt((2.0_dp * real(l, dp) + 1.0_dp) / 2.0_dp * &
                    factorial(l - m) / factorial(l + m))

        pmm = 1.0_dp
        if (m > 0) then
            somx2 = sqrt((1.0_dp - x) * (1.0_dp + x))
            fact = 1.0_dp
            do i = 1, m
                pmm = -1.0_dp * pmm * fact * somx2
                fact = fact + 2.0_dp
            end do
        end if

        if (l == m) then
            plm = pmm * norm
            return
        end if

        pmmp1 = x * real(2 * m + 1, dp) * pmm

        if (l == m + 1) then
            plm = pmmp1 * norm
            return
        end if

        do ll = m + 2, l
            pll = (x * real(2 * ll - 1, dp) * pmmp1 - real(ll + m - 1, dp) * pmm) / real(ll - m, dp)
            pmm = pmmp1
            pmmp1 = pll
        end do

        plm = pll * norm
    end function assoc_legendre

    function legendre_poly(l, x) result(p)
        implicit none

        integer, intent(in) :: l
        real(dp), intent(in) :: x
        real(dp) :: p
        real(dp) :: p_prev, p_curr
        integer :: i

        if (l == 0) then
            p = 1.0_dp
        else if (l == 1) then
            p = x
        else
            p_prev = 1.0_dp
            p_curr = x
            do i = 1, l - 1
                p = ((2.0_dp * i + 1.0_dp) * x * p_curr - real(i, dp) * p_prev) / (i + 1.0_dp)
                p_prev = p_curr
                p_curr = p
            end do
        end if
    end function legendre_poly




  End MODULE AtomDiatomskvp
!
!*********************************************************************************************
!*********************************************************************************************
