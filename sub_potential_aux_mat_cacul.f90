!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE potential_mat_calcul

        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib
        use Potential_Interface

        IMPLICIT NONE

        INTEGER :: i, i_prime
        INTEGER :: n, n_prime
        INTEGER :: row, col
        INTEGER :: ipot
        INTEGER :: Ntot

!
!       Total closed-channel matrix dimension
!       -------------------------------------
        Ntot = dim_x * ncf

!
!       Safety check
!       ------------
        IF (.NOT. ALLOCATED(pot_mat)) THEN
                PRINT*, 'Error: pot_mat is not allocated before potential_mat_calcul.'
                PRINT*, 'Call build_potential_index first.'
                STOP
        ENDIF

!
!       Build radial and angular potential BAMs
!       ---------------------------------------
        CALL build_BAM_rs
        CALL build_BAM_thetas

!
!       Allocate matrices
!       -----------------
        IF (ALLOCATED(M_V)) DEALLOCATE(M_V, STAT = istatus)
        ALLOCATE(M_V(1:Ntot,1:Ntot), STAT = istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for M_V in potential_mat_calcul.'
                STOP
        ENDIF
        M_V = 0d0

        IF (ALLOCATED(M0_V)) DEALLOCATE(M0_V, STAT = istatus)
        ALLOCATE(M0_V(1:Ntot,1:n_open), STAT = istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for M0_V in potential_mat_calcul.'
                STOP
        ENDIF
        M0_V = (0d0, 0d0)

        IF (ALLOCATED(M00_V)) DEALLOCATE(M00_V, STAT = istatus)
        ALLOCATE(M00_V(1:n_open,1:n_open), STAT = istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for M00_V in potential_mat_calcul.'
                STOP
        ENDIF
        M00_V = (0d0, 0d0)

        IF (ALLOCATED(M10_V)) DEALLOCATE(M10_V, STAT = istatus)
        ALLOCATE(M10_V(1:n_open,1:n_open), STAT = istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for M10_V in potential_mat_calcul.'
                STOP
        ENDIF
        M10_V = (0d0, 0d0)

!
!       Build M_V: closed-closed potential block
!       ----------------------------------------
        DO i = 1, dim_x
        DO n = 1, ncf

                row = (i-1)*ncf + n

                DO i_prime = 1, dim_x
                DO n_prime = 1, ncf

                        col = (i_prime-1)*ncf + n_prime

                        DO ipot = 1, n_pot
                                M_V(row,col) = M_V(row,col) + &
                                        BAM_r(ipot,i+1,i_prime+1) * &
                                        BAM_theta(ipot,n,n_prime)
                        ENDDO

                ENDDO
                ENDDO

        ENDDO
        ENDDO

!
!       Build M0_V: closed-open potential block
!       ---------------------------------------
        DO i = 1, dim_x
        DO n = 1, ncf

                row = (i-1)*ncf + n

                DO n_prime = 1, n_open

                        DO ipot = 1, n_pot
                                M0_V(row,n_prime) = M0_V(row,n_prime) + &
                                        BAM_r0(ipot,i+1,n_prime) * &
                                        BAM_theta(ipot,n,open_idx(n_prime))
                        ENDDO

                ENDDO

        ENDDO
        ENDDO

!
!       Build M00_V and M10_V: open-open potential blocks
!       --------------------------------------------------
        DO n = 1, n_open
        DO n_prime = 1, n_open

                DO ipot = 1, n_pot

                        M00_V(n,n_prime) = M00_V(n,n_prime) + &
                                BAM_r00(ipot,n,n_prime) * &
                                BAM_theta(ipot,open_idx(n),open_idx(n_prime))

                        M10_V(n,n_prime) = M10_V(n,n_prime) + &
                                BAM_r10(ipot,n,n_prime) * &
                                BAM_theta(ipot,open_idx(n),open_idx(n_prime))

                ENDDO

        ENDDO
        ENDDO

END SUBROUTINE potential_mat_calcul
!
!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE build_BAM_rs

        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib
        use Potential_Interface


        IMPLICIT NONE

        INTEGER :: i, i_prime
        INTEGER :: n, n_prime
        INTEGER :: l, ipot
        INTEGER :: q
        REAL(8) :: xm, xr, dx, sxv
        COMPLEX(8) :: sxv0

        CALL read_A_cache(knots_x, gq_root_x, ngqp_x)
!
!       Allocate BAM_r
!       --------------
        IF (ALLOCATED(BAM_r)) THEN
                DEALLOCATE(BAM_r, STAT = istatus)
        ENDIF

        ALLOCATE(BAM_r(1:n_pot,1:pbasst(1)%pb_nbr,1:pbasst(1)%pb_nbr), STAT = istatus)

        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for BAM_r in build_BAM_rs.'
                STOP
        ENDIF

        BAM_r = 0d0

        DO i = 1, pbasst(1)%pb_nbr
        DO i_prime = i, MIN(i+pbasst(1)%pb_pa1-1, pbasst(1)%pb_nbr)
        DO q = i_prime, i+pbasst(1)%pb_pa1-1

                xm = 0.5d0 * (knots_x(q+1) + knots_x(q))
                xr = 0.5d0 * (knots_x(q+1) - knots_x(q))

                DO ipot = 1, n_pot
                        sxv = 0d0

                        DO l = 1, ngqp_x/2
                                dx = xr * gq_root_x(l+ngqp_x/2)

                               sxv = sxv + gq_weight_x(l+ngqp_x/2) * &
                                        (A_cache(ipot,q,l,1) * bsp_x(i,xm+dx) * bsp_x(i_prime,xm+dx) + &
                                                A_cache(ipot,q,l,2) * bsp_x(i,xm-dx) * bsp_x(i_prime,xm-dx) )
                        ENDDO

                        BAM_r(ipot,i,i_prime) = BAM_r(ipot,i,i_prime) + &
                                sxv * xr / DSQRT(norm(i)*norm(i_prime))

                        BAM_r(ipot,i_prime,i) = BAM_r(ipot,i,i_prime)
                ENDDO

        ENDDO
        ENDDO
        ENDDO

!
!       Allocate BAM_r0
!       ---------------
!       Important:
!       u0(n,R) is COMPLEX(8), so BAM_r0 must be COMPLEX(8).
!
        IF (ALLOCATED(BAM_r0)) THEN
                DEALLOCATE(BAM_r0, STAT = istatus)
        ENDIF

        ALLOCATE(BAM_r0(1:n_pot,1:pbasst(1)%pb_nbr,1:n_open), STAT = istatus)

        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for BAM_r0 in build_BAM_rs.'
                STOP
        ENDIF

        BAM_r0 = (0d0, 0d0)

!
!       Calculate closed-open radial potential matrix:
!
!       BAM_r0(ipot,i,n_prime)
!       =
!       Int dR B_i(R) A_ipot(R) u0_n_prime(R)
!       / sqrt(norm(i))
!
!       u0 already contains:
!       exp(-i*k*R) * h(R) * sqrt(mu_R/k)
!
        BAM_r0 = (0d0, 0d0)

        DO i = 1, pbasst(1)%pb_nbr
        DO n_prime = 1, n_open

                DO q = i, i+pbasst(1)%pb_pa1-1
                        xm = 0.5d0 * (knots_x(q+1) + knots_x(q))
                        xr = 0.5d0 * (knots_x(q+1) - knots_x(q))

                        DO ipot = 1, n_pot
                                sxv0 = (0d0, 0d0)

                                DO l = 1, ngqp_x/2
                                        dx = xr * gq_root_x(l+ngqp_x/2)

                                        sxv0 = sxv0 + gq_weight_x(l+ngqp_x/2) * &
                                                ( bsp_x(i,xm+dx) * A_cache(ipot,q,l,1) * u0(n_prime,xm+dx) + &
                                                bsp_x(i,xm-dx) * A_cache(ipot,q,l,2) * u0(n_prime,xm-dx) )
                                ENDDO

                                BAM_r0(ipot,i,n_prime) = BAM_r0(ipot,i,n_prime) + &
                                        sxv0 * xr / DSQRT(norm(i))
                        ENDDO
                ENDDO

        ENDDO
        ENDDO


        !
!       Allocate BAM_r00 and BAM_r10
!       ----------------------------
        IF (ALLOCATED(BAM_r00)) THEN
                DEALLOCATE(BAM_r00, STAT = istatus)
        ENDIF

        IF (ALLOCATED(BAM_r10)) THEN
                DEALLOCATE(BAM_r10, STAT = istatus)
        ENDIF

        ALLOCATE(BAM_r00(1:n_pot,1:n_open,1:n_open), STAT = istatus)

        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for BAM_r00 in build_BAM_rs.'
                STOP
        ENDIF

        ALLOCATE(BAM_r10(1:n_pot,1:n_open,1:n_open), STAT = istatus)

        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for BAM_r10 in build_BAM_rs.'
                STOP
        ENDIF

        BAM_r00 = (0d0, 0d0)
        BAM_r10 = (0d0, 0d0)

        DO n = 1, n_open
        DO n_prime = 1, n_open

                DO q = pbasst(1)%pb_pa1, pbasst(1)%pb_nbr
                        xm = 0.5d0 * (knots_x(q+1) + knots_x(q))
                        xr = 0.5d0 * (knots_x(q+1) - knots_x(q))

                        DO ipot = 1, n_pot
                                sxv0 = (0d0, 0d0)

                                DO l = 1, ngqp_x/2
                                        dx = xr * gq_root_x(l+ngqp_x/2)

                                        sxv0 = sxv0 + gq_weight_x(l+ngqp_x/2) * &
                                                ( u0(n,xm+dx) * A_cache(ipot,q,l,1) * u0(n_prime,xm+dx) + &
                                                        u0(n,xm-dx) * A_cache(ipot,q,l,2) * u0(n_prime,xm-dx) )
                                ENDDO

                                BAM_r00(ipot,n,n_prime) = BAM_r00(ipot,n,n_prime) + sxv0 * xr

                                sxv0 = (0d0, 0d0)

                                DO l = 1, ngqp_x/2
                                        dx = xr * gq_root_x(l+ngqp_x/2)

                                        sxv0 = sxv0 + gq_weight_x(l+ngqp_x/2) * &
                                                        ( CONJG(u0(n,xm+dx)) * A_cache(ipot,q,l,1) * u0(n_prime,xm+dx) + &
                                                        CONJG(u0(n,xm-dx)) * A_cache(ipot,q,l,2) * u0(n_prime,xm-dx) )
                                ENDDO

                                BAM_r10(ipot,n,n_prime) = BAM_r10(ipot,n,n_prime) + sxv0 * xr
                        ENDDO
                ENDDO

        ENDDO
        ENDDO

END SUBROUTINE build_BAM_rs
!
!*********************************************************************************************
!=============================================================================================

!*********************************************************************************************
!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE build_BAM_thetas
!
!*********************************************************************************************
!
! This subroutine builds the angular potential matrix:
!
!   BAM_theta(ipot,n,n_prime)
!
! where
!
!   ipot <----> (lambda1, lambda2, m)
!   n    <----> (j1, k1, j2, k2)
!   n'   <----> (j1', k1', j2', k2')
!
!=============================================================================================
!
        USE generateparameters
        USE AtomDiatomskvp
        USE omp_lib
        use Potential_Interface

!
        IMPLICIT NONE
!
        INTEGER :: ipot, n, n_prime
        INTEGER :: lambda1, lambda2, m
        INTEGER :: j1, k1, j2, k2
        INTEGER :: j1_prime, k1_prime, j2_prime, k2_prime
        REAL(8) :: G1p, G2p, G1m, G2m, basis_factor
!

      
!       Allocate BAM_theta
!       ------------------
        IF (ALLOCATED(BAM_theta)) THEN
                DEALLOCATE(BAM_theta, STAT = istatus)
        ENDIF

        ALLOCATE(BAM_theta(1:n_pot,1:ncf,1:ncf), STAT = istatus)

        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for BAM_theta in build_BAM_theta.'
                STOP
        ENDIF

        BAM_theta = 0d0
!
!       Calculate BAM_theta(ipot,n,n_prime)
!       -----------------------------------
        DO ipot = 1, n_pot

                        lambda1 = pot_mat(1,ipot)
                        lambda2 = pot_mat(2,ipot)
                        m       = pot_mat(3,ipot)


                DO n = 1, ncf

                        j1 = quant_mat(1,n)
                        k1 = quant_mat(2,n)
                        j2 = quant_mat(3,n)
                        k2 = quant_mat(4,n)

                        DO n_prime = 1, ncf

                                j1_prime = quant_mat(1,n_prime)
                                k1_prime = quant_mat(2,n_prime)
                                j2_prime = quant_mat(3,n_prime)
                                k2_prime = quant_mat(4,n_prime)

                                G1p = gaunt_coeff(j1,k1,lambda1, m,j1_prime,k1_prime)
                                G2p = gaunt_coeff(j2,k2,lambda2,-m,j2_prime,k2_prime)

                                G1m = gaunt_coeff(j1,k1,lambda1,-m,j1_prime,k1_prime)
                                G2m = gaunt_coeff(j2,k2,lambda2, m,j2_prime,k2_prime)

                                IF (m == 0) THEN
                                basis_factor = DSQRT(2d0*pi)
                                ELSE
                                basis_factor = (-1d0)**m * 2d0 * DSQRT(pi)
                                ENDIF

                                BAM_theta(ipot,n,n_prime) = basis_factor * 0.5d0 * (G1p*G2p + G1m*G2m)

                        ENDDO
                ENDDO
        ENDDO

        !PRINT*, 'BAM_theta built successfully.'
        !PRINT*, 'Shape: ', n_pot, ncf, ncf
!
!=============================================================================================
!
END SUBROUTINE build_BAM_thetas
!
!*********************************************************************************************
!=============================================================================================







