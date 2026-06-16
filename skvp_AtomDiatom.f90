!*********************************************************************************************
!*********************************************************************************************
!
PROGRAM skvpAtomDiatom
!
!*********************************************************************************************
!
!
!=============================================================================================
!
!	Date		:	04.04.2022
!	Author		:	Steve Ndengue 
!	Last Change	:	03.29.2026
!
!	The present program is aimed at applying the S-Matrix Kohn Variational Principle 
!	       on the 2D scattering problem using the traditional method from Miller.
! 
!=============================================================================================		
! 
 USE omp_lib
 USE generateparameters
 USE AtomDiatomskvp
 USE Potential_Interface
!


 IMPLICIT NONE
!-------------
!
 INTEGER :: i, j
 INTEGER :: t
 INTEGER :: j1, j2, j1low, j2low
 INTEGER :: k1, k2, k1range, k2range
 INTEGER                                  :: step, n_steps 

 REAL(8)                                  :: tm1, tm2, pas_x, yy=0.020d0
 REAL(8)                                  :: x1, x2              !, DLAMCH
 !COMPLEX(8)                               :: varc
 !REAL(8), ALLOCATABLE, DIMENSION(:)       :: var
 !COMPLEX(8), ALLOCATABLE, DIMENSION(:)    :: biatx
!
! Call for input paramaters reading routine
!------------------------------------------
!      
CALL read_input

!CALL set_potential_backend('BMKP')
!CALL set_bmkp_filename('/Users/yuan/Documents/skvp_diatomdiaton_firstdrfat/coefficients.dat')


CALL build_potential_index

CALL calculate_RMS_potential_expansion

                DO t = 1, n_pot
                WRITE(*,'(4I6)') t, pot_mat(1,t), pot_mat(2,t), pot_mat(3,t)
                ENDDO

        CALL CPU_TIME(tm1)
!
        OPEN (100,file ='proba.dat',status='unknown',action='write',position='append')
!
!
! Call for input paramaters reading routine
!------------------------------------------
!       
        dim_x    = pbasst(1)%pb_nbr-2
        ngqp_x   = p*pbasst(1)%pb_pa1  !2*int(rhomax-rhomin)*p
! 
! Determination of the knots sequences
!-------------------------------------
        ! Allocation of matrices
        ALLOCATE ( knots_x(1:(pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)), gq_root_x(1:ngqp_x), gq_weight_x(1:ngqp_x), STAT = istatus )
        ! Determination of the abcisses and weights for gaussian quadrature
        CALL gauleg(-1D0,1D0,gq_root_x,gq_weight_x,ngqp_x)
        WRITE (1, 590) ! Write gaussian quadrature points in global output file
        DO i = 1, ngqp_x
                WRITE (1, 595)  gq_root_x(i), gq_weight_x(i)
        ENDDO
        ! y = 1d0
        pas_x = 0d0
        knots_x = 0d0
        pas_x=(pbasst(1)%pb_max-pbasst(1)%pb_min)/(dexp(yy*dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1)/ &
               dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2))-dble(1))
        knots_x=(/(pbasst(1)%pb_min,i=1,pbasst(1)%pb_pa1), &
                  ((dexp(yy*dble(i-1)/dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2))-dble(1))*pas_x+ &
                  pbasst(1)%pb_min,i=2,pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1), &
                  (pbasst(1)%pb_max,i=pbasst(1)%pb_nbr+1,pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)/)
        ! pas_x=(pbasst(1)%pb_max-pbasst(1)%pb_min)/dble(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1)
        ! knots_x=(/(pbasst(1)%pb_min,i=1,pbasst(1)%pb_pa1),(dble(i)*pas_x+pbasst(1)%pb_min,i=1,pbasst(1)%pb_nbr-pbasst(1)%pb_pa1), &
        !           (pbasst(1)%pb_max,i=pbasst(1)%pb_nbr+1,pbasst(1)%pb_nbr+pbasst(1)%pb_pa1)/)
        PRINT*,'knots_x'
        DO i = 1, SIZE(knots_x)
        IF (i < SIZE(knots_x)) THEN
            WRITE(*,'(F6.2,A)', ADVANCE='NO') knots_x(i), ','
        ELSE
            WRITE(*,'(F6.2)') knots_x(i)
        END IF
        END DO
!
! Calculate the number of steps (using Nearest Integer to be safe) in the energy range and loop
!----------------------------------------------------------------------------------------------
        n_steps = nint((en_final - en_start) / en_step)
        ! Energy loop -- loopnmax = 1001
        DO step = 0, n_steps
        E = en_start + (real(step,8)*en_step)
        !E = dble(loopn-1)*(0.0025d0/10d0)
        !!!E = 0.000d0 + dble(loopn-1)*(0.040d0/1000d0)
        !E = 0.0d0 + dble(loopnmax-loopn+1)*(0.0367493d0/10d0)
!
! Computation of the targets wavefunctions and loop over Jtot
!------------------------------------------------------------
        ! Jtot loop
        DO Jtot = Jmin, Jmax
        ! (j,k) -> n, Compacting this set of quantum numbers into one quantum number 
        !!!ncf = (((2*Jtot + 1)*pbasst(2)%pb_nbr)/2)+1
        ! Compute ncf

        ncf=0
        do j1=0, pbasst(2)%pb_nbr, 2 !=> Revisit later
                j1low=min(j1,pbasst(2)%pb_pa1)
                k1range=j1low
                !k1range=min(j1low,Jtot)
                do k1=-k1range, k1range
                        do j2=0, pbasst(3)%pb_nbr, 2 !=> Revisit later
                                j2low=min(j2,pbasst(3)%pb_pa1)
                                k2range=j2low
                                !k2range=min(j2low,Jtot)
                                do k2=-k2range, k2range
                                        ncf=ncf+1
                                enddo    
                        enddo
                enddo
        enddo
        ! Compute quant_mat
        ALLOCATE(quant_mat(4,ncf))
        ncf = 0
        do j1=0, pbasst(2)%pb_nbr, 2 !=> Revisit later
                j1low=min(j1,pbasst(2)%pb_pa1)
                k1range=j1low
                !k1range=min(j1low,Jtot)
                do k1=-k1range, k1range
                        do j2=0, pbasst(3)%pb_nbr, 2 !=> Revisit later
                                j2low=min(j2,pbasst(3)%pb_pa1)
                                k2range=j2low
                                !k2range=min(j2low,Jtot)
                                do k2=-k2range, k2range
                                        ncf=ncf+1
                                        quant_mat(1,ncf)=j1
                                        quant_mat(2,ncf)=k1
                                        quant_mat(3,ncf)=j2
                                        quant_mat(4,ncf)=k2
                                enddo  
                        enddo
                enddo
        enddo

        PRINT*, 'channel index: j1 k1 j2 k2'
        DO i = 1, ncf
                WRITE(*,'(I5,4I6)') i, quant_mat(1,i), quant_mat(2,i), quant_mat(3,i), quant_mat(4,i)
        ENDDO


        ! Print parameters
        PRINT*, " "
        PRINT*, "ngqp_x:", ngqp_x
        PRINT*, "ncf", ncf
        PRINT*, "ENERGY: ", E
        PRINT*, "alpha: ", alpha
        PRINT*, "R0: ", r0
        PRINT*, "rhomin: ", pbasst(1)%pb_min
        PRINT*, "rhomax", pbasst(1)%pb_max
        PRINT*, "mu_R: ", mu_R
        PRINT*, "Brot", Brot
        PRINT*, "Jtot", Jtot
        PRINT*, "quant_mat"
!       PRINT*, quant_mat(1, :) !quant_mat(1, 1:10)
 !       PRINT*, quant_mat(2, :) !quant_mat(2, 1:10)
 !       PRINT*, quant_mat(3, :) !quant_mat(3, 1:10)
 !       PRINT*, quant_mat(4, :) !quant_mat(4, 1:10)
!
        !ALLOCATE(quant_mat(2,ncf))
        !nn = 1
        !DO j = 0, pbasst(2)%pb_nbr, 2
        !   IF (j == 0) THEN
        !      quant_mat(1,nn) = 0
        !      quant_mat(2,nn) = 0
        !      nn = nn + 1
        !   ELSE
        !      DO k = -Jtot,Jtot
        !         quant_mat(1,nn) = j
        !         quant_mat(2,nn) = k
        !         nn = nn + 1
        !      ENDDO
        !   ENDIF
        !ENDDO
        !PRINT*," " 
        !!PRINT*, "bsp(0, 2.1)", bsp_x(0, 2.0d0)
        !PRINT*, "quant_mat"
        !PRINT*, quant_mat(1, 1:10)
        !PRINT*, quant_mat(2, 1:10)
!
        CALL solve_target_levels
        x1 = pbasst(1)%pb_min   
        x2 = pbasst(1)%pb_max
        ALLOCATE(x(1:ngqp_x), wx(1:ngqp_x), STAT = istatus) 
        CALL gauleg(x1, x2, x, wx, ngqp_x)
!
! Calculate the Cross Section
!-----------------------------------------
        CALL CrossSection

!
        WRITE(6,*) 'PROBABILITIES at SOME ENERGY - ', step
        WRITE(6,'(10E15.7)') E*27.211399d0,((abs(Smat(1,j))**2.D0),j=1,n_open)
        WRITE(100,'(16E15.7)') E*27.211399d0, ((abs(Smat(1,j))**2.D0),j=1,n_open)

        !DEALLOCATE(knots_x, gq_root_x, gq_weight_x, knots_y, gq_root_y, gq_weight_y, STAT = istatus)

        DEALLOCATE(quant_mat, STAT = istatus)
        DEALLOCATE(BAM_1, BAM_3, BAM_4, norm, STAT = istatus)
        DEALLOCATE(BAM_x1, BAM_x3, BAM_x4, STAT = istatus)
        DEALLOCATE(BAM_xx1, BAM_xx3, BAM_xx4, BAM_xx1b, BAM_xx3b, BAM_xx4b, STAT = istatus)
        DEALLOCATE(M_V, M0_V, M00_V, M10_V, M_K, Smat, STAT = istatus)
        DEALLOCATE(mat_M, mat_M0, mat_M00, mat_M10, STAT = istatus)

        IF (ALLOCATED(kvec)) THEN
                DEALLOCATE(kvec, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r)) THEN
        DEALLOCATE(BAM_r, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r0)) THEN
                DEALLOCATE(BAM_r0, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r00)) THEN
                DEALLOCATE(BAM_r00, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r10)) THEN
                DEALLOCATE(BAM_r10, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_theta)) THEN
                DEALLOCATE(BAM_theta, STAT = istatus)
        ENDIF

        DEALLOCATE(M_V, M0_V, M00_V, M10_V, M_K, Smat, STAT = istatus)
        DEALLOCATE(mat_M, mat_M0, mat_M00, mat_M10, STAT = istatus)

        ENDDO ! Jtot
        ENDDO ! E


        !print*,'legendre pol - l=2,0.5', plgndr(2,0,0.5d0), plgndr(2,-1,0.5d0), plgndr(2,1,0.5d0) 
        !print*,'legendre pol - l=2,-0.5', plgndr(2,0,-0.5d0), plgndr(2,-1,-0.5d0), plgndr(2,1,-0.5d0)
        !print*,'legendre pol - l=2,1', plgndr(2,0,1.0d0), plgndr(2,-1,1.0d0), plgndr(2,1,1.0d0)
        !print*,'legendre pol - l=2,-1', plgndr(2,0,-1.0d0), plgndr(2,-1,-1.0d0), plgndr(2,1,-1.0d0)
   

!
! Writing format
!---------------
! 500    FORMAT (8E15.7)
 590    FORMAT (10X,'INITIAL QUADRATURE ROOTS AND WEIGHTS FOR QUADRATURE:')
 595    FORMAT (8X, D12.5, D12.5)
! 600    FORMAT (6E15.7)                                                  
!
        CALL CPU_TIME(tm2)
        write(6,'(1A,F8.3,1A)') 'execution time:',(tm2-tm1)/6d1,'min'
        PRINT*," "
!
!*********************************************************************************************
!*********************************************************************************************
!
 END PROGRAM skvpAtomDiatom
!
!*********************************************************************************************
!=============================================================================================
!*********************************************************************************************
!
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
  SUBROUTINE CrossSection
!
!*********************************************************************************************
!
! This subroutine creates the total cross section using the probability result of the S-matrix
! as described in Ruthie (2013)
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
!
        IMPLICIT NONE
!------------------------------
        INTEGER :: i, j
        REAL(8), ALLOCATABLE, DIMENSION(:)   :: Kj2, Coeff
        REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Smat_prob, Dsec, Xsec
        REAL(8) :: E_rot

        
        !inquire(file='matrices.bin', exist=file_exists, size=file_size)
!
!        IF (file_size == 0) THEN
!                CALL noc_independent_calc
!        ELSE 
!        ALLOCATE ( BAM_1(1:n_x,1:n_x), BAM_3(1:n_x,1:n_x), BAM_4(1:n_x,1:n_x), &
!        BAM_r1(n_x, n_x), BAM_r2(n_x,n_x), BAM_r3(n_x,n_x), norm(n_x))
!
!        BAM_1=complex(0D0,0D0)
!        BAM_3=complex(0D0,0D0)
!        BAM_4=complex(0D0,0D0)
!        BAM_r1 = 0d0
!        BAM_r2 = 0d0
!        BAM_r3 = 0d0
!        norm = 0d0
!        ENDIF

!!! DO Jtot = Jmin, Jmax
!--------------------
! Creation of the basic auxiliary matrices
!-----------------------------------------
        CALL basic_aux_mat_calcul
!--------------------
! Creation of the potential matrices
!-----------------------------------------
        CALL potential_mat_calcul
!               
!               
! Creation of the main scattering matrices
!-----------------------------------------
        CALL make_scatt_mat

!
! Sample test routines call
!--------------------------
        CALL PhaseShift
!
!-------------------------------
!       Allocation
!------------------------
        ALLOCATE(Kj2(1:n_open), Coeff(1:n_open), STAT = istatus)
        ALLOCATE(Smat_prob(1:n_open,1:n_open), Dsec(1:n_open,1:n_open), Xsec(1:n_open,1:n_open), STAT = istatus)
!----------------------------
!       Find Coefficient for Xsec
!-----------------------------     
    DO j = 1, min(n_open,7)
        DO i = 1, min(n_open,7)
            Smat_prob(i,j) = (abs(Smat(i,j))**2.D0) !defining probability matrix
        ENDDO
        write(6,'(7E18.8)') (Smat_prob(i,j), i=1, min(n_open,7)) !!问一下史蒂夫
                
        E_rot = Brot*quant_mat(1, open_idx(j))*(quant_mat(1, open_idx(j)) + 1d0) + &
               Brot*quant_mat(3, open_idx(j))*(quant_mat(3, open_idx(j)) + 1d0)


       Kj2(j) = mu_R * (E - E_rot) !why there is no 2
        ! make this for a vector E instead of single value
        Coeff(j) = pi/((2*quant_mat(1,j)+1d0)*Kj2(j)) !coefficient for cross section
        IF (Coeff(j) < 0) THEN
                PRINT*, 'Warning: Negative Xsec Coefficient'
        ENDIF
    ENDDO

    !PRINT*, 'Kj2', Kj2
    !PRINT*, 'Coeff', Coeff
    !PRINT*, ' '

!       Angular momentum component of Xsec
!--------------------------------------------
    !Dsec = 0.0  
    !DO Jtotal = Jmin, Jmax !summing over Jtot
        !Dsec = Dsec + (2*Jtot+1)*Smat_prob
        !PRINT*, 'Dsec', Jtotal
        !DO i = 1, 3
            !DO j = 1, 3
                !PRINT*, Dsec(i,j)
            !ENDDO
        !ENDDO
    !ENDDO

!       Combining Results for Xsec result
!--------------------------------------------
    !DO j=1, min(n_open,7)
    !    Xsec(j,:) = Dsec(j,:) * Coeff(j) !Coeff number determined by which state we are coming from
    !ENDDO

    !PRINT*, 'Cross Section'
    !DO i = 1, min(n_open,7)
    !    write(6,'(7E18.8)') (Xsec(i,j), j=1, min(n_open,7))
    !ENDDO
    !PRINT*, ' '

!       Deallocating Matrices
! ---------------------------------
        DEALLOCATE ( Kj2, Coeff, Smat_prob, Dsec, STAT = istatus )

 !!!ENDDO

!=============================================================================================
!
 END SUBROUTINE CrossSection
!
!=============================================================================================
!
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
 SUBROUTINE make_scatt_mat
!
!*********************************************************************************************
!
! This subroutine creates the relevant scattering matrices M, M_0, M_0,0, M_1,0
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
        USE Potential_Interface
!
        IMPLICIT NONE
!--------------------
        INTEGER i, j, k, i1, i2, j1, j2, N, quant_j, quant_j_prime, &
        quant_k, quant_k_prime
        INTEGER :: quant_j1, quant_k1, quant_j2, quant_k2
        INTEGER :: quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime
        REAL(8) :: E_rot
        REAL(8) :: channel_delta
        !REAL(8) modul, shift2delta
        !INTEGER, ALLOCATABLE, DIMENSION(:)      :: ipivx
!
        N = dim_x*ncf
!
!       Allocation
!       ----------
        ALLOCATE ( mat_M(1:N,1:N), mat_M0(1:N,1:n_open), STAT = istatus )
        ALLOCATE ( mat_M00(1:n_open,1:n_open), mat_M10(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE ( M_K(N,N))

        
!
        mat_M=0d0
        mat_M0=COMPLEX(0d0,0d0)
        mat_M00=COMPLEX(0d0,0d0)
        mat_M10=COMPLEX(0d0,0d0)
        !M_V = 0d0
        !M0_V = (0d0,0d0)
        !M00_V = (0d0,0d0)
        !M10_V = 0d0
        M_K = 0d0
!
!       Formation of the various M matrices
!       -----------------------------------
        DO i1=1, dim_x   !1 refers to initial state
        DO i2=1, ncf     !2 refers to final state

        i = (i1-1)*ncf+i2
        quant_j1 = quant_mat(1,i2)
        quant_k1 = quant_mat(2,i2)
        quant_j2 = quant_mat(3,i2)
        quant_k2 = quant_mat(4,i2)

        DO j1=1, dim_x
        DO j2=1, ncf
        j = (j1-1)*ncf+j2
        quant_j1_prime = quant_mat(1,j2)
        quant_k1_prime = quant_mat(2,j2)
        quant_j2_prime = quant_mat(3,j2)
        quant_k2_prime = quant_mat(4,j2)
        

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        M_K(i,j) = (1d0/(2d0*mu_R)) * BAM_1(i1+1,j1+1) * channel_delta
        
        E_rot = Brot*quant_j1*(quant_j1+1d0) + Brot*quant_j2*(quant_j2+1d0)

        mat_M(i,j)= M_K(i,j) + M_V(i,j) &
        - E*BAM_3(i1+1,j1+1)*channel_delta &
        + E_rot*BAM_3(i1+1,j1+1)*channel_delta & 
        + (1d0/(2d0*mu_R)) * BAM_4(i1+1,j1+1) * (quant_j1 * (quant_j1 + 1d0) + &
          quant_j2 * (quant_j2 + 1d0)) * channel_delta &
        + (1d0/(2d0*mu_R))*BAM_4(i1+1,j1+1)* Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
             quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime) ! 1/R**2 term
        ENDDO
        ENDDO
!
        DO k=1, n_open

        quant_j1_prime = quant_mat(1,open_idx(k))
        quant_k1_prime = quant_mat(2,open_idx(k))
        quant_j2_prime = quant_mat(3,open_idx(k))
        quant_k2_prime = quant_mat(4,open_idx(k))

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                    delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        E_rot = Brot * quant_j1 * (quant_j1 + 1d0) + Brot * quant_j2 * (quant_j2 + 1d0)

        mat_M0(i,k) = (-1d0/(2d0*mu_R))*BAM_x1(i1+1,k)*channel_delta + M0_V(i,k) &
                      - E*BAM_x3(i1+1,k)*channel_delta & 
                      + E_rot*channel_delta*BAM_x3(i1+1,k) & 
                      + (1d0/(2d0*mu_R)) * BAM_x4(i1+1,k) * &
                        (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                      + (1d0/(2d0*mu_R)) * BAM_x4(i1+1,k) * &
                        Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime) ! 1/R**2 term
        ENDDO

        ENDDO
        ENDDO
!
        DO i = 1, n_open !n
        DO j = 1, n_open !n_prime
        quant_j1 = quant_mat(1,open_idx(i))
        quant_k1 = quant_mat(2,open_idx(i))
        quant_j2 = quant_mat(3,open_idx(i))
        quant_k2 = quant_mat(4,open_idx(i))

        quant_j1_prime = quant_mat(1,open_idx(j))
        quant_k1_prime = quant_mat(2,open_idx(j))
        quant_j2_prime = quant_mat(3,open_idx(j))
        quant_k2_prime = quant_mat(4,open_idx(j))

        channel_delta = delta(quant_j1,quant_j1_prime) * delta(quant_k1,quant_k1_prime) * &
                    delta(quant_j2,quant_j2_prime) * delta(quant_k2,quant_k2_prime)

        E_rot = Brot * quant_j1 * (quant_j1 + 1d0) + &
                Brot * quant_j2 * (quant_j2 + 1d0)

        mat_M00(i,j)= -(1d0/(2d0*mu_R))*BAM_xx1(i,j)*channel_delta+M00_V(i,j) &
                      - E*BAM_xx3(i,j)*channel_delta + &
                      E_rot*channel_delta*BAM_xx3(i,j) &
                      + (1d0/(2d0*mu_R)) * BAM_xx4(i,j) * &
                        (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                      + (1d0/(2d0*mu_R))*BAM_xx4(i,j)* &
                        Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime)

        mat_M10(i,j)= -(1d0/(2d0*mu_R))*BAM_xx1b(i,j)*channel_delta+M10_V(i,j) &
                         - E*BAM_xx3b(i,j)*channel_delta &
                       + E_rot*channel_delta*BAM_xx3b(i,j) &
                       + (1d0/(2d0*mu_R)) * BAM_xx4b(i,j) * &
                         (quant_j1*(quant_j1 + 1d0)+quant_j2*(quant_j2 + 1d0))*channel_delta &
                       + (1d0/(2d0*mu_R))*BAM_xx4b(i,j)* &
                         Wdd(Jtot, quant_j1, quant_k1, quant_j2, quant_k2, &
                              quant_j1_prime, quant_k1_prime, quant_j2_prime, quant_k2_prime)
        ENDDO
        ENDDO
!
!        PRINT*,'M_V'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)') (dble(M_V(i,j)), j=1, min(N,6))
!        enddo
!        PRINT*,'M0_V'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)') (dble(M0_V(i,j)), j=1, min(n_open,6))
!        enddo
!        PRINT*,'mat_M_V2'
!        do i=1, min(N,6)
!        write(6,'(8E18.8)'),(dble(M_V2(i,j)), j=1, min(N,6))
!        enddo

!        PRINT*,'M00_V'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)') (dble(M00_V(i,:)))
!        enddo
!        PRINT*,'M10_V'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)') (dble(M10_V(i,:)))
!        enddo
!         PRINT*,'mat_M00_V2'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)'),(dble(M00_V2(i,:)))
!        enddo
!        PRINT*,'mat_M10_V2'
!        do i=1, min(N,n_open)
!        write(6,'(8E18.8)'),(dble(M10_V2(i,:)))
!       enddo


!       PRINT*,'mat_M'
!       do i=1, min(N,5)
!       write(6,'(8E18.8)') (dble(mat_M(i,j)), j=1, min(N,5))
!       enddo
!       PRINT*,'mat_M0'
!       do i= N-10, N
!       write(6,'(8E18.8)') (dble(mat_M0(i,j)), j=1, min(n_open,4))
!       enddo
!        PRINT*,'mat_M0_V'
!        do i=1, min(N,4)
!        write(6,'(8E18.8)'),(dble(M0_V(i,j)), j=1, min(n_open,4))
!        enddo
!       PRINT*,'mat_M_K'
!       do i=1, min(N,4)
!       write(6,'(8E18.8)') (dble(M_K(i,j)), j=1, min(n_open,4))
!       enddo
!       PRINT*,'mat_M00'
!       do i=1, min(n_open,5)
!       write(6,'(8E18.8)') (dble(mat_M00(i,:)))
!       enddo
!       PRINT*,'mat_M10'
!       do i=1, min(n_open,5)
!       write(6,'(8E18.8)') (dble(mat_M10(i,:)))
!       enddo


        ! OPEN(unit=10, file="M_Matrix.txt", status='unknown', position='append', action='write')
        ! WRITE(10,*) 'Matrix M:'
        ! DO i = 1, 5
        ! DO j = 1, 5
        !     WRITE(10, '(ES14.5)', advance='no') mat_M(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)

        !  WRITE(10,*) 'Matrix M_V:'
        !   DO i = 1, 15
        ! DO j = 1, 15
        !     WRITE(10,'(F15.3)', advance='no') M_V(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)

        !  WRITE(10,*) 'Matrix M_K:'
        !   DO i = 1, 15
        ! DO j = 1, 15
        !     WRITE(10,'(F15.3)', advance='no') M_K(i,j)
        ! END DO
        ! WRITE(10,*)
        ! END DO
        !  WRITE(10,*)
        ! CLOSE(10)        


!
! 600    FORMAT (6E15.7)
!
!PRINT*, '==== matrix NaN/debug check ===='
!PRINT*, 'max M_V     = ', MAXVAL(ABS(M_V))
!PRINT*, 'max M0_V    = ', MAXVAL(ABS(M0_V))
!PRINT*, 'max M00_V   = ', MAXVAL(ABS(M00_V))
!PRINT*, 'max M10_V   = ', MAXVAL(ABS(M10_V))
!PRINT*, 'max mat_M   = ', MAXVAL(ABS(mat_M))
!PRINT*, 'max mat_M0  = ', MAXVAL(ABS(mat_M0))
!PRINT*, 'max mat_M00 = ', MAXVAL(ABS(mat_M00))
!PRINT*, 'max mat_M10 = ', MAXVAL(ABS(mat_M10))

DO i = 1, n_open
        PRINT*, 'open channel ', i, &
                ' j1 k1 j2 k2 = ', quant_mat(1,open_idx(i)), &
                quant_mat(2,open_idx(i)), &
                quant_mat(3,open_idx(i)), &
                quant_mat(4,open_idx(i)), &
                ' kvec = ', kvec(i)
ENDDO

PRINT*, '==============================='
!=============================================================================================
!
 END SUBROUTINE make_scatt_mat
!
!=============================================================================================
!
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************!
 
 SUBROUTINE PhaseShift
!
!*********************************************************************************************
!
! This subroutine creates the relevant matrices and values and computes the S-matrix
!
!=============================================================================================
!
        USE AtomDiatomskvp
        USE generateparameters
        USE omp_lib
!
        IMPLICIT NONE
!--------------------
        INTEGER i, j, infox, lworkx, N
        
        !REAL(8) tm11, tm12, shift2delta
        INTEGER, ALLOCATABLE, DIMENSION(:)      :: ipivx
        REAL(8), ALLOCATABLE, DIMENSION(:)      :: sum_array, workx
        REAL(8), ALLOCATABLE, DIMENSION(:,:)    :: modul, mat_M_inv
        !COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: ttdum
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: dum1vx, Bsub, Csub
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: dum2vx, dum3vx, dum4vx, Smat_st
        COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: mat_B_st_inv, dum5vx, dum6vx

!
        lworkx = 50*(pbasst(1)%pb_nbr**2)
        N = dim_x*ncf
!
!       Allocation
!       ----------
        ALLOCATE( mat_M_inv(1:N,1:N), ipivx(1:N), mat_B_st_inv(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE( workx(1:lworkx), dum3vx(1:n_open,1:N), dum4vx(1:n_open,1:N), Smat(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE( dum1vx(1:N,1:n_open), dum2vx(1:N,1:n_open), Smat_st(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE( mat_B(1:n_open,1:n_open), mat_C(1:n_open,1:n_open), modul(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE( Bsub(1:n_open,1:n_open), Csub(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE( sum_array(n_open), STAT = istatus )
!
        mat_M_inv=mat_M
!        

!       Formation of the various final matrices
!       ---------------------------------------

        !print*,'here2'
        !print*,'dimension mat_M_inv:', size(mat_M_inv, dim=1), size(mat_M_inv, dim=2)

        call dsytrf('L', N, mat_M_inv, N, ipivx, workx, N, infox)
        if (infox /= 0) then
        print *, 'Error: DSYTRF failed with info = ', infox
        stop
        end if

        !print*,'dimension mat_M_inv:', size(mat_M_inv, dim=1), size(mat_M_inv, dim=2)
        !print*,'dimension ipivx:', size(ipivx, dim=1)
        !print*,'dimension workx:', size(workx, dim=1)

        ! Compute the inverse using DSYTRI
        call dsytri('L', N, mat_M_inv, N, ipivx, workx, infox)
        if (infox /= 0) then
        print *, 'Error: DSYTRI failed with info = ', infox
        stop
        end if

        print*,'dimension mat_M_inv:', size(mat_M_inv, dim=1), size(mat_M_inv, dim=2)

        !=> Why the need to enforce the symmetry (below)? (DSYTRF/I routines)
        do i = 1, N
           do j = i+1, N
           mat_M_inv(i,j) = mat_M_inv(j,i)
           end do
        end do

        !print*,'here3'

!       CALL ZGETRF(N, N, mat_M_inv, N, ipivx, infox)
!       PRINT*,'infox=',infox
!       CALL ZGETRI(N, mat_M_inv, N, ipivx, workx, lworkx, infox)
!       PRINT*,'infox=',infox
!       PRINT*, 'workx', workx
       ! ttdum=MATMUL(mat_M,mat_M_inv)
!        PRINT*,'ttdum'
!        !print*,'here4 - n_open=',n_open
!        do i=1, min(n_open,5)
!        write(6,'(8E18.8)') (REAL(ttdum(i,j)), j=1, min(n_open,5))
!        enddo    
!        PRINT*'mat_M_inv'
!        do i=1, 5
!        write(6,'(8E18.8)') (dble(mat_M_inv(i,j)), j=1, 5)
!        enddo     

!      
        dum1vx = MATMUL(mat_M_inv,mat_M0)
        dum2vx = CONJG(mat_M0)
        dum3vx = TRANSPOSE(mat_M0)
        dum4vx = TRANSPOSE(dum2vx)
        Bsub=MATMUL(dum3vx,dum1vx)
        mat_B = mat_M00 - Bsub
!        mat_B = MATMUL(dum3vx,dum1vx)
!!        print*,'submat_B'
!!        do i=1, min(n_open,4)
!!        write(6,'(8E18.8)') (Bsub(i,j), j=1, min(n_open,4))
!!        enddo
!        PRINT*,'mat_B'
!        do i=1, min(n_open,7)
!        write(6,'(8E18.8)') (mat_B(i,j), j=1, min(n_open,7))
!        enddo

        Csub=MATMUL(dum4vx,dum1vx)
        mat_C = mat_M10 - Csub
!        mat_C = MATMUL(dum4vx,dum1vx)
!        print*,'submat_C'
!        do i=1, min(n_open,4)
!        write(6,'(8E18.8)') (Csub(i,j), j=1, min(n_open,4))
!        enddo
!        PRINT*,'mat_C'
!        do i=1, min(n_open,4)
!        write(6,'(8E18.8)') (mat_C(i,j), j=1, min(n_open,7))
!        enddo
!
!       Scattering Matrix
!       -----------------
        DEALLOCATE ( workx, ipivx, STAT = istatus )
        ALLOCATE ( dum5vx(1:n_open,1:n_open), dum6vx(1:n_open,1:n_open), STAT = istatus )
        ALLOCATE ( workx(1:lworkx), ipivx(1:n_open), STAT = istatus )
        mat_B_st_inv = CONJG(mat_B)
        CALL ZGETRF(n_open, n_open, mat_B_st_inv, n_open, ipivx, infox)
!       PRINT*,'infox=',infox
        CALL ZGETRI(n_open, mat_B_st_inv, n_open, ipivx, workx, lworkx, infox)
!       PRINT*,'infox=',infox
        !ttdum=MATMUL(CONJG(mat_B),mat_B_st_inv)
        dum5vx=TRANSPOSE(mat_C)
        dum6vx=MATMUL(mat_B_st_inv,mat_C)
        Smat = (0d0,1d0)*(mat_B - MATMUL(dum5vx,dum6vx))  !(val_C**2)/(CONJG(val_B))


      

        Smat_st = TRANSPOSE(CONJG(Smat))

      


!       print Sdegger to check 
       ! DO i = 1, n_open
       !         dum5vx(i,i) = dum5vx(i,i) - (1d0,0d0)
       ! ENDDO

       ! PRINT*, 'max |Sdagger*S - I| = ', MAXVAL(ABS(dum5vx))


        modul=MATMUL(Smat_st,Smat)
        PRINT*,' '
!        PRINT*,'Module'
!        do i=1, min(n_open,7)
!        write(6,'(7E18.8)') ((modul(i,j)), j=1, min(n_open,7))
!        enddo
        PRINT*,' '
        PRINT*, 'Probabilities of first', min(7,n_open), 'channel'

        do i=1, min(n_open,7)
        write(6,'(7E18.8)') ((abs(Smat(i,j))**2D0), j=1, min(n_open,7))
        enddo

        PRINT*, ' '
!*!        shift2delta=dacos(-1d0)+datan2(dimag(Smat)/modul,dble(Smat)/modul)
!*!        PRINT*,'Phase Shift=',shift2delta/2d0,(shift2delta-dacos(-1d0))/2d0
!        PRINT*, 'Probabilities: j1,k1,j2,k2 -> j1p,k1p,j2p,k2p'
!
 !       DO i = 1, n_open
 !       DO j = 1, n_open
  !              write(6,'(E18.8)', advance='no') ABS(Smat(i,j))**2D0
   !             PRINT*, quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
   !                     quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
   !                     ' -> ', &
   !                     quant_mat(1,open_idx(j)), quant_mat(2,open_idx(j)), &
   !                     quant_mat(3,open_idx(j)), quant_mat(4,open_idx(j))
    !            write(*,*)
    !    ENDDO
    !    ENDDO

        sum_array = 0d0
        DO i = 1, n_open
        DO j = 1, n_open
             sum_array(i) = sum_array(i) + (abs(Smat(i,j))**2)
        ENDDO
        ENDDO
        PRINT*, "Checking probabilities sum to one: ", sum_array(:)
        PRINT*, ' '
        ! Adding data to output file
        ! IF (n_open >= 2) THEN
        !         open(unit=10, file="2_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, (abs(Smat(2,1))**2D0)
        !         close(10)
        ! ELSE 
        !         open(unit=10, file="2_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, 0d0
        !         close(10)
        ! ENDIF
        ! IF (n_open >= 3) THEN
        !         open(unit=10, file="4_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, (abs(Smat(3,1))**2D0)
        !         close(10)
        ! ELSE
        !         open(unit=10, file="4_0.txt", status='unknown', position='append', action='write')
        !         write(10, '(ES24.16, 2X, ES24.16)') E, 0d0
        !         close(10)
        ! ENDIF


        IF (n_open >= 1) THEN
                open(unit=10, file="0_0_all.txt", status='unknown', position='append', action='write')
                DO i = 1, n_open
                        write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') &
                                quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
                                quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
                                E, ABS(Smat(1,i))**2D0
                ENDDO
                close(10)
        ELSE
                open(unit=10, file="0_0_all.txt", status='unknown', position='append', action='write')
                write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') 0, 0, 0, 0, E, 0d0
                close(10)
        ENDIF


        IF (n_open >= 4) THEN
                open(unit=10, file="open4_all.txt", status='unknown', position='append', action='write')
                DO i = 1, n_open
                        write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') &
                                quant_mat(1,open_idx(i)), quant_mat(2,open_idx(i)), &
                                quant_mat(3,open_idx(i)), quant_mat(4,open_idx(i)), &
                                E, ABS(Smat(4,i))**2D0
                ENDDO
                close(10)
        ELSE
                open(unit=10, file="open4_all.txt", status='unknown', position='append', action='write')
                write(10, '(4I5, 2X, ES24.16, 2X, ES24.16)') 0, 0, 0, 0, E, 0d0
                close(10)
        ENDIF
!
!       Deallocating Matrices
!       ---------------------
        DEALLOCATE ( workx, ipivx, STAT = istatus )
        DEALLOCATE ( mat_M_inv, STAT = istatus )
        DEALLOCATE ( dum1vx, dum2vx, Smat_st, STAT = istatus )
        DEALLOCATE ( dum3vx, dum4vx, dum5vx, dum6vx, STAT = istatus )
        DEALLOCATE ( mat_B, mat_C, mat_B_st_inv, modul, STAT = istatus )
        DEALLOCATE ( Bsub, Csub, sum_array, STAT = istatus )
!   
!
! 600    FORMAT (6E15.7)
!
!=============================================================================================
!
 END SUBROUTINE PhaseShift
!
!=============================================================================================
!*********************************************************************************************
!*********************************************************************************************
