!*********************************************************************************************
!*********************************************************************************************
!
PROGRAM test_MV_direct
!
!*********************************************************************************************
!
! Standalone test driver for the closed-channel potential matrix M_V.
!
! The production M_V path is:
!
!   build_potential_index  -> Potential_Interface.f90
!   solve_target_levels    -> sub_basic_aux_mat_calcul.f90
!   potential_mat_calcul   -> sub_potential_aux_mat_cacul.f90
!
! The direct integral below is only a check inside this test file.
!
!=============================================================================================
!
        USE generateparameters
        USE AtomDiatomskvp
        USE Potential_Interface

!
        IMPLICIT NONE
!
        EXTERNAL :: solve_target_levels
        EXTERNAL :: potential_mat_calcul
!
!       User-controlled test parameters.
!       --------------------------------
        INTEGER, PARAMETER :: lambda1_max_test = 2
        INTEGER, PARAMETER :: lambda2_max_test = 2
        INTEGER, PARAMETER :: m_max_test       = 2
!
!       Direct angular quadrature size.
!       Increase these if Integral_V_direct is not converged.
!       -----------------------------------------------------
        INTEGER, PARAMETER :: nang_test = 5
        INTEGER, PARAMETER :: nphi_test = 5
        INTEGER, PARAMETER :: matrix_block_test = 5
        INTEGER, PARAMETER :: radial_closed_closed = 1
        INTEGER, PARAMETER :: radial_closed_open   = 2
        INTEGER, PARAMETER :: radial_open_open_00  = 3
        INTEGER, PARAMETER :: radial_open_open_10  = 4
!
!       Test run parameters.
!       --------------------
        INTEGER, PARAMETER :: ca_nmd_test = 4
        REAL(8), PARAMETER :: en_start_test = 1.0d-1
        INTEGER, PARAMETER :: Jtot_test = 2
        REAL(8), PARAMETER :: r0_test = 25.0d0
        REAL(8), PARAMETER :: mu_R_test = 1.8371695E+03
        REAL(8), PARAMETER :: Brot_test = 2.7d-4
        LOGICAL, PARAMETER :: print_term_details = .FALSE.
        LOGICAL, PARAMETER :: run_MV_direct_test = .FALSE.
        LOGICAL, PARAMETER :: run_M0V_direct_test = .TRUE.
        LOGICAL, PARAMETER :: run_M00V_direct_test = .TRUE.
        LOGICAL, PARAMETER :: run_M10V_direct_test = .TRUE.
!
        INTEGER :: row, col, Ntot, nshow, nshow_open
        INTEGER :: max_row, max_col
        REAL(8) :: direct_val, diff, rel_diff
        REAL(8) :: max_abs_diff
!
!       Initialize the same global modules used by the main SKVP calculation.
!       --------------------------------------------------------------------
        !CALL set_potential_backend('BMKP')
        !CALL set_bmkp_filename('/Users/yuan/Documents/skvp_diatomdiaton_firstdrfat/coefficients.dat')
        CALL set_test_parameters
        CALL setup_radial_grid_and_knots
        CALL build_quant_mat_for_test
        CALL build_norm_for_test
        CALL build_potential_index
        CALL solve_target_levels
        CALL potential_mat_calcul

        
!       Compare M_V(row,col) with Integral_V_direct(row,col).
!       -----------------------------------------------------
        Ntot = dim_x * ncf
        nshow = MIN(matrix_block_test, Ntot)
!
        IF (run_MV_direct_test) THEN
!
        PRINT*, ' '
        PRINT*, '=============================================================='
        PRINT*, 'Testing M_V(row,col) against Integral_V_direct(row,col)'
        PRINT*, 'Ntot  = ', Ntot
        PRINT*, 'ncf   = ', ncf
        PRINT*, 'dim_x = ', dim_x
        PRINT*, 'n_pot = ', n_pot
        PRINT*, 'nang  = ', nang_test
        PRINT*, 'nphi  = ', nphi_test
        PRINT*, '=============================================================='
        PRINT*, ' '
!
        CALL print_potential_index
        CALL print_quant_channels
!
        IF (print_term_details) THEN
                CALL print_MV_term_details(nshow)
        ENDIF
!
        max_abs_diff = 0d0
        max_row = 1
        max_col = 1
!
        PRINT*, ' row  col            M_V                 Integral_direct          abs_diff          rel_diff'
        PRINT*, '--------------------------------------------------------------------------------------------'
!
        DO row = 1, nshow
        DO col = 1, nshow
!
                direct_val = Integral_V_direct(row, col, nang_test, nphi_test)
                diff = DABS(M_V(row,col) - direct_val)
!
                IF (DABS(M_V(row,col)) > 1d-14) THEN
                        rel_diff = diff / DABS(M_V(row,col))
                ELSE
                        rel_diff = diff
                ENDIF
!
                WRITE(*,'(I5,1X,I5,4(1X,ES22.12))') row, col, M_V(row,col), &
                        direct_val, diff, rel_diff
!
                IF (diff > max_abs_diff) THEN
                        max_abs_diff = diff
                        max_row = row
                        max_col = col
                ENDIF
!
        ENDDO
        ENDDO
!
        PRINT*, ' '
        PRINT*, 'Maximum absolute difference in printed block:'
        PRINT*, 'max_abs_diff = ', max_abs_diff
        PRINT*, 'at row, col  = ', max_row, max_col
        PRINT*, ' '
!
        ENDIF
!
!       Compare M0_V, M00_V, and M10_V against direct integrals.
!       ---------------------------------------------------------
        IF (n_open > 0) THEN
                nshow_open = MIN(matrix_block_test, n_open)
                IF (run_M0V_direct_test) CALL print_M0V_direct_test(nshow, nshow_open)
                IF (run_M00V_direct_test) CALL print_M00V_direct_test(nshow_open)
                IF (run_M10V_direct_test) CALL print_M10V_direct_test(nshow_open)
        ELSE
                PRINT*, 'No open channels: skip M0_V, M00_V, and M10_V direct tests.'
        ENDIF
!
!=============================================================================================
!
CONTAINS


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE set_test_parameters
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        ca_type  = 'collision'
        ca_spec  = 1102
        ca_sym   = 'para'
        ca_nmd   = ca_nmd_test
        ca_nest  = 1
        r0       = r0_test
        en_start = en_start_test
        en_final = en_start_test
        en_step  = en_start_test
        Jmin     = Jtot_test
        Jmax     = Jtot_test
        E        = en_start
        Jtot     = Jtot_test
        mu_R     = mu_R_test
        Arot     = Brot_test
        Brot     = Brot_test
        Crot     = Brot_test
!
        IF (ALLOCATED(pbasst)) DEALLOCATE(pbasst, STAT=istatus)
        ALLOCATE(pbasst(1:ca_nmd), STAT=istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for pbasst.'
                STOP
        ENDIF
!
        pbasst(1)%pb_nme = 'R'
        pbasst(1)%pb_typ = 'bspl_fbr'
        pbasst(1)%pb_nbr = 20
        pbasst(1)%pb_pa1 = 4
        pbasst(1)%pb_min = 2.0d0
        pbasst(1)%pb_max = 20.0d0
        pbasst(1)%pb_pa2 = 0.0d0
!
        pbasst(2)%pb_nme = 'th1'
        pbasst(2)%pb_typ = 'asleg_fbr'
        pbasst(2)%pb_nbr = 4
        pbasst(2)%pb_pa1 = 2
        pbasst(2)%pb_min = 0.0d0
        pbasst(2)%pb_max = 3.1415d0
        pbasst(2)%pb_pa2 = 0.0d0
!
        pbasst(3)%pb_nme = 'th2'
        pbasst(3)%pb_typ = 'asleg_fbr'
        pbasst(3)%pb_nbr = 4
        pbasst(3)%pb_pa1 = 2
        pbasst(3)%pb_min = 0.0d0
        pbasst(3)%pb_max = 3.1415d0
        pbasst(3)%pb_pa2 = 0.0d0
!
        pbasst(4)%pb_nme = 'phi'
        pbasst(4)%pb_typ = 'exp'
        pbasst(4)%pb_nbr = 2
        pbasst(4)%pb_pa1 = 1
        pbasst(4)%pb_min = 0.0d0
        pbasst(4)%pb_max = 6.2830d0
        pbasst(4)%pb_pa2 = 0.0d0
!
        dim_x = pbasst(1)%pb_nbr - 2
        ngqp_x = 6
!
        PRINT*, ' '
        PRINT*, 'Test parameters'
        PRINT*, 'lambda1_max, lambda2_max, m_max = ', &
                lambda1_max_test, lambda2_max_test, m_max_test
        PRINT*, 'E, Jtot = ', E, Jtot
        PRINT*, 'mu_R, Brot, r0 = ', mu_R, Brot, r0
        PRINT*, 'R basis: n, k, min, max = ', pbasst(1)%pb_nbr, pbasst(1)%pb_pa1, &
                pbasst(1)%pb_min, pbasst(1)%pb_max
        PRINT*, ' '
!
END SUBROUTINE set_test_parameters


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE setup_radial_grid_and_knots
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER :: i
        REAL(8) :: yy, pas_x
!
        yy = 0.020d0
!
        IF (ALLOCATED(knots_x)) DEALLOCATE(knots_x, STAT=istatus)
        IF (ALLOCATED(gq_root_x)) DEALLOCATE(gq_root_x, STAT=istatus)
        IF (ALLOCATED(gq_weight_x)) DEALLOCATE(gq_weight_x, STAT=istatus)
!
        ALLOCATE(knots_x(1:pbasst(1)%pb_nbr+pbasst(1)%pb_pa1), &
                 gq_root_x(1:ngqp_x), gq_weight_x(1:ngqp_x), STAT=istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for knots_x/gq arrays.'
                STOP
        ENDIF
!
        CALL gauleg(-1d0, 1d0, gq_root_x, gq_weight_x, ngqp_x)
!
        pas_x = (pbasst(1)%pb_max - pbasst(1)%pb_min) / &
                (DEXP(yy * DBLE(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1) / &
                DBLE(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2)) - 1d0)
!
        knots_x = (/ (pbasst(1)%pb_min, i=1,pbasst(1)%pb_pa1), &
                ((DEXP(yy*DBLE(i-1)/DBLE(pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+2)) - 1d0) * &
                pas_x + pbasst(1)%pb_min, i=2,pbasst(1)%pb_nbr-pbasst(1)%pb_pa1+1), &
                (pbasst(1)%pb_max, i=pbasst(1)%pb_nbr+1,pbasst(1)%pb_nbr+pbasst(1)%pb_pa1) /)
!
        PRINT*, 'Radial grid initialized.'
        PRINT*, 'dim_x, ngqp_x = ', dim_x, ngqp_x
!
END SUBROUTINE setup_radial_grid_and_knots


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE build_quant_mat_for_test
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER :: j1, j2, j1low, j2low
        INTEGER :: k1, k2, k1range, k2range
        INTEGER :: n
!
        IF (ALLOCATED(quant_mat)) DEALLOCATE(quant_mat, STAT=istatus)
!
        ncf = 0
        DO j1 = 0, pbasst(2)%pb_nbr, 2
                j1low = MIN(j1, pbasst(2)%pb_pa1)
                k1range = MIN(j1low, Jtot)
                DO k1 = -k1range, k1range
                        DO j2 = 0, pbasst(3)%pb_nbr, 2
                                j2low = MIN(j2, pbasst(3)%pb_pa1)
                                k2range = MIN(j2low, Jtot)
                                DO k2 = -k2range, k2range
                                        ncf = ncf + 1
                                ENDDO
                        ENDDO
                ENDDO
        ENDDO
!
        ALLOCATE(quant_mat(4,ncf), STAT=istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for quant_mat.'
                STOP
        ENDIF
!
        n = 0
        DO j1 = 0, pbasst(2)%pb_nbr, 2
                j1low = MIN(j1, pbasst(2)%pb_pa1)
                k1range = MIN(j1low, Jtot)
                DO k1 = -k1range, k1range
                        DO j2 = 0, pbasst(3)%pb_nbr, 2
                                j2low = MIN(j2, pbasst(3)%pb_pa1)
                                k2range = MIN(j2low, Jtot)
                                DO k2 = -k2range, k2range
                                        n = n + 1
                                        quant_mat(1,n) = j1
                                        quant_mat(2,n) = k1
                                        quant_mat(3,n) = j2
                                        quant_mat(4,n) = k2
                                ENDDO
                        ENDDO
                ENDDO
        ENDDO
!
        PRINT*, 'quant_mat built. ncf = ', ncf
!
END SUBROUTINE build_quant_mat_for_test


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE build_norm_for_test
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER :: i_x, q, l
        REAL(8) :: xm, xr, dx, sx
!
        IF (ALLOCATED(norm)) DEALLOCATE(norm, STAT=istatus)
        ALLOCATE(norm(1:pbasst(1)%pb_nbr), STAT=istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for norm.'
                STOP
        ENDIF
!
        norm = 0d0
!
        DO i_x = 1, pbasst(1)%pb_nbr
        DO q = i_x, i_x + pbasst(1)%pb_pa1 - 1
                xm = 0.5d0 * (knots_x(q+1) + knots_x(q))
                xr = 0.5d0 * (knots_x(q+1) - knots_x(q))
                sx = 0d0
!
                DO l = 1, ngqp_x / 2
                        dx = xr * gq_root_x(l + ngqp_x/2)
                        sx = sx + gq_weight_x(l + ngqp_x/2) * &
                                (bsp_x(i_x,xm+dx)**2 + bsp_x(i_x,xm-dx)**2)
                ENDDO
!
                norm(i_x) = norm(i_x) + sx * xr
        ENDDO
        ENDDO
!
        IF (MINVAL(norm) <= 0d0) THEN
                PRINT*, 'Error: non-positive B-spline norm.'
                PRINT*, 'min(norm) = ', MINVAL(norm)
                STOP
        ENDIF
!
        PRINT*, 'B-spline norms built.'
!
END SUBROUTINE build_norm_for_test


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_potential_index
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER :: ipot
!
        PRINT*, ' '
        PRINT*, 'Potential index map:'
        PRINT*, 'ipot        lambda1        lambda2        m'
        DO ipot = 1, n_pot
                WRITE(*,'(4I12)') ipot, pot_mat(1,ipot), pot_mat(2,ipot), pot_mat(3,ipot)
        ENDDO
        PRINT*, ' '
!
END SUBROUTINE print_potential_index


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_quant_channels
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER :: n, nprint
!
        nprint = MIN(ncf, matrix_block_test)
!
        PRINT*, 'First quantum channels used by printed block:'
        PRINT*, '   n          j1          k1          j2          k2'
        DO n = 1, nprint
                WRITE(*,'(5I12)') n, quant_mat(1,n), quant_mat(2,n), &
                        quant_mat(3,n), quant_mat(4,n)
        ENDDO
        PRINT*, ' '
!
END SUBROUTINE print_quant_channels


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_MV_term_details(nshow)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: nshow
!
        INTEGER :: row, col, ipot
        INTEGER :: i, i_prime, n, n_prime
        REAL(8) :: contribution
!
        PRINT*, ' '
        PRINT*, 'Detailed separated M_V contributions'
        PRINT*, 'ipot lambda1 lambda2 m            BAM_r              BAM_theta          Total'
        PRINT*, '--------------------------------------------------------------------------------------------'
!
        DO row = 1, nshow
                i = (row - 1) / ncf + 1
                n = MOD(row - 1, ncf) + 1
!
                DO col = 1, nshow
                        i_prime = (col - 1) / ncf + 1
                        n_prime = MOD(col - 1, ncf) + 1
!
                        WRITE(*,'("row col =",2I6,"   i n ip np =",4I6)') &
                                row, col, i, n, i_prime, n_prime
!
                        DO ipot = 1, n_pot
                                contribution = BAM_r(ipot,i+1,i_prime+1) * &
                                        BAM_theta(ipot,n,n_prime)
!
                                WRITE(*,'(4I6,3(1X,ES22.12))') &
                                        ipot, pot_mat(1,ipot), pot_mat(2,ipot), pot_mat(3,ipot), &
                                        BAM_r(ipot,i+1,i_prime+1), BAM_theta(ipot,n,n_prime), &
                                        contribution
                        ENDDO
!
                        WRITE(*,'("      sum(M_V) =",1X,ES22.12)') M_V(row,col)
                        PRINT*, ' '
                ENDDO
        ENDDO
!
END SUBROUTINE print_MV_term_details


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_M0V_direct_test(nrow_show, nopen_show)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: nrow_show, nopen_show
!
        INTEGER :: row, n_prime
        INTEGER :: max_row, max_col
        REAL(8) :: diff, rel_diff, max_abs_diff
        COMPLEX(8) :: direct_val
!
        PRINT*, ' '
        PRINT*, '=============================================================='
        PRINT*, 'Testing M0_V(row,open_col) against Integral_M0V_direct'
        PRINT*, 'closed rows printed = ', nrow_show
        PRINT*, 'open cols printed   = ', nopen_show
        PRINT*, '=============================================================='
        PRINT*, ' '
        PRINT*, ' row open       Re(M0_V)          Im(M0_V)        Re(direct)        Im(direct)       abs_diff        rel_diff'
        PRINT*, '--------------------------------------------------------------------------------------------------------------'
!
        max_abs_diff = 0d0
        max_row = 1
        max_col = 1
!
        DO row = 1, nrow_show
        DO n_prime = 1, nopen_show
                direct_val = Integral_M0V_direct(row, n_prime, nang_test, nphi_test)
                diff = ABS(M0_V(row,n_prime) - direct_val)
!
                IF (ABS(M0_V(row,n_prime)) > 1d-14) THEN
                        rel_diff = diff / ABS(M0_V(row,n_prime))
                ELSE
                        rel_diff = diff
                ENDIF
!
                WRITE(*,'(I5,1X,I5,6(1X,ES17.8))') row, n_prime, &
                        DBLE(M0_V(row,n_prime)), AIMAG(M0_V(row,n_prime)), &
                        DBLE(direct_val), AIMAG(direct_val), diff, rel_diff
!
                IF (diff > max_abs_diff) THEN
                        max_abs_diff = diff
                        max_row = row
                        max_col = n_prime
                ENDIF
        ENDDO
        ENDDO
!
        PRINT*, ' '
        PRINT*, 'Maximum absolute difference in M0_V printed block:'
        PRINT*, 'max_abs_diff = ', max_abs_diff
        PRINT*, 'at row, open = ', max_row, max_col
        PRINT*, ' '
!
END SUBROUTINE print_M0V_direct_test


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_M00V_direct_test(nopen_show)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: nopen_show
!
        INTEGER :: n, n_prime
        INTEGER :: max_row, max_col
        REAL(8) :: diff, rel_diff, max_abs_diff
        COMPLEX(8) :: direct_val
!
        PRINT*, ' '
        PRINT*, '=============================================================='
        PRINT*, 'Testing M00_V(open_row,open_col) against Integral_M00V_direct'
        PRINT*, 'open block printed = ', nopen_show
        PRINT*, '=============================================================='
        PRINT*, ' '
        PRINT*, 'open open      Re(M00_V)         Im(M00_V)       Re(direct)        Im(direct)       abs_diff        rel_diff'
        PRINT*, '--------------------------------------------------------------------------------------------------------------'
!
        max_abs_diff = 0d0
        max_row = 1
        max_col = 1
!
        DO n = 1, nopen_show
        DO n_prime = 1, nopen_show
                direct_val = Integral_M00V_direct(n, n_prime, nang_test, nphi_test)
                diff = ABS(M00_V(n,n_prime) - direct_val)
!
                IF (ABS(M00_V(n,n_prime)) > 1d-14) THEN
                        rel_diff = diff / ABS(M00_V(n,n_prime))
                ELSE
                        rel_diff = diff
                ENDIF
!
                WRITE(*,'(I5,1X,I5,6(1X,ES17.8))') n, n_prime, &
                        DBLE(M00_V(n,n_prime)), AIMAG(M00_V(n,n_prime)), &
                        DBLE(direct_val), AIMAG(direct_val), diff, rel_diff
!
                IF (diff > max_abs_diff) THEN
                        max_abs_diff = diff
                        max_row = n
                        max_col = n_prime
                ENDIF
        ENDDO
        ENDDO
!
        PRINT*, ' '
        PRINT*, 'Maximum absolute difference in M00_V printed block:'
        PRINT*, 'max_abs_diff = ', max_abs_diff
        PRINT*, 'at open, open = ', max_row, max_col
        PRINT*, ' '
!
END SUBROUTINE print_M00V_direct_test


!*********************************************************************************************
!*********************************************************************************************
!
SUBROUTINE print_M10V_direct_test(nopen_show)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: nopen_show
!
        INTEGER :: n, n_prime
        INTEGER :: max_row, max_col
        REAL(8) :: diff, rel_diff, max_abs_diff
        COMPLEX(8) :: direct_val
!
        PRINT*, ' '
        PRINT*, '=============================================================='
        PRINT*, 'Testing M10_V(open_row,open_col) against Integral_M10V_direct'
        PRINT*, 'open block printed = ', nopen_show
        PRINT*, '=============================================================='
        PRINT*, ' '
        PRINT*, 'open open      Re(M10_V)         Im(M10_V)       Re(direct)        Im(direct)       abs_diff        rel_diff'
        PRINT*, '--------------------------------------------------------------------------------------------------------------'
!
        max_abs_diff = 0d0
        max_row = 1
        max_col = 1
!
        DO n = 1, nopen_show
        DO n_prime = 1, nopen_show
                direct_val = Integral_M10V_direct(n, n_prime, nang_test, nphi_test)
                diff = ABS(M10_V(n,n_prime) - direct_val)
!
                IF (ABS(M10_V(n,n_prime)) > 1d-14) THEN
                        rel_diff = diff / ABS(M10_V(n,n_prime))
                ELSE
                        rel_diff = diff
                ENDIF
!
                WRITE(*,'(I5,1X,I5,6(1X,ES17.8))') n, n_prime, &
                        DBLE(M10_V(n,n_prime)), AIMAG(M10_V(n,n_prime)), &
                        DBLE(direct_val), AIMAG(direct_val), diff, rel_diff
!
                IF (diff > max_abs_diff) THEN
                        max_abs_diff = diff
                        max_row = n
                        max_col = n_prime
                ENDIF
        ENDDO
        ENDDO
!
        PRINT*, ' '
        PRINT*, 'Maximum absolute difference in M10_V printed block:'
        PRINT*, 'max_abs_diff = ', max_abs_diff
        PRINT*, 'at open, open = ', max_row, max_col
        PRINT*, ' '
!
END SUBROUTINE print_M10V_direct_test


!*********************************************************************************************
!*********************************************************************************************
!
REAL(8) FUNCTION Integral_V_direct(row, col, nang, nphi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: row, col
        INTEGER, INTENT(IN) :: nang, nphi
!
        Integral_V_direct = DBLE(Integral_potential_direct(row, col, 0, 0, &
                nang, nphi, radial_closed_closed))
!
END FUNCTION Integral_V_direct


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION Integral_M0V_direct(row, open_col, nang, nphi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: row, open_col
        INTEGER, INTENT(IN) :: nang, nphi
!
        Integral_M0V_direct = Integral_potential_direct(row, 0, 0, open_col, &
                nang, nphi, radial_closed_open)
!
END FUNCTION Integral_M0V_direct


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION Integral_M00V_direct(open_row, open_col, nang, nphi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: open_row, open_col
        INTEGER, INTENT(IN) :: nang, nphi
!
        Integral_M00V_direct = Integral_potential_direct(0, 0, open_row, open_col, &
                nang, nphi, radial_open_open_00)
!
END FUNCTION Integral_M00V_direct


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION Integral_M10V_direct(open_row, open_col, nang, nphi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: open_row, open_col
        INTEGER, INTENT(IN) :: nang, nphi
!
        Integral_M10V_direct = Integral_potential_direct(0, 0, open_row, open_col, &
                nang, nphi, radial_open_open_10)
!
END FUNCTION Integral_M10V_direct


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION Integral_potential_direct(row, col, open_row, open_col, &
        nang, nphi, radial_mode)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: row, col
        INTEGER, INTENT(IN) :: open_row, open_col
        INTEGER, INTENT(IN) :: nang, nphi
        INTEGER, INTENT(IN) :: radial_mode
!
        INTEGER :: i_basis, i_prime_basis
        INTEGER :: n, n_prime
        INTEGER :: j1, k1, j2, k2
        INTEGER :: j1_prime, k1_prime, j2_prime, k2_prime
        INTEGER :: q, q_left, q_right, l
        INTEGER :: ix1, ix2, iphi1, iphi2
        INTEGER :: Ntot
!
        REAL(8), ALLOCATABLE :: xang(:), wang(:)
        REAL(8) :: R, x1ang, x2ang
        REAL(8) :: xm, xr
        REAL(8) :: phi1, phi2, dphi, wphi
        COMPLEX(8) :: radial_part
        COMPLEX(8) :: sum_complex
        COMPLEX(8) :: Vang
        COMPLEX(8) :: y_left_1, y_right_1
        COMPLEX(8) :: y_left_2, y_right_2
!
        Ntot = dim_x * ncf
!
        SELECT CASE (radial_mode)
        CASE (radial_closed_closed)
                IF (row < 1 .OR. row > Ntot) THEN
                        PRINT*, 'Error: row is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                IF (col < 1 .OR. col > Ntot) THEN
                        PRINT*, 'Error: col is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                i_basis = (row - 1) / ncf + 2
                n       = MOD(row - 1, ncf) + 1
!
                i_prime_basis = (col - 1) / ncf + 2
                n_prime       = MOD(col - 1, ncf) + 1
!
                q_left  = MAX(i_basis, i_prime_basis)
                q_right = MIN(i_basis+pbasst(1)%pb_pa1-1, &
                        i_prime_basis+pbasst(1)%pb_pa1-1)
!
        CASE (radial_closed_open)
                IF (row < 1 .OR. row > Ntot) THEN
                        PRINT*, 'Error: row is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                IF (open_col < 1 .OR. open_col > n_open) THEN
                        PRINT*, 'Error: open_col is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                i_basis = (row - 1) / ncf + 2
                i_prime_basis = 0
                n       = MOD(row - 1, ncf) + 1
                n_prime = open_idx(open_col)
!
                q_left  = i_basis
                q_right = i_basis + pbasst(1)%pb_pa1 - 1
!
        CASE (radial_open_open_00, radial_open_open_10)
                IF (open_row < 1 .OR. open_row > n_open) THEN
                        PRINT*, 'Error: open_row is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                IF (open_col < 1 .OR. open_col > n_open) THEN
                        PRINT*, 'Error: open_col is outside range in Integral_potential_direct.'
                        STOP
                ENDIF
!
                i_basis = 0
                i_prime_basis = 0
                n       = open_idx(open_row)
                n_prime = open_idx(open_col)
!
                q_left  = pbasst(1)%pb_pa1
                q_right = pbasst(1)%pb_nbr
!
        CASE DEFAULT
                PRINT*, 'Error: unknown radial_mode in Integral_potential_direct.'
                STOP
        END SELECT
!
        j1 = quant_mat(1,n)
        k1 = quant_mat(2,n)
        j2 = quant_mat(3,n)
        k2 = quant_mat(4,n)
!
        j1_prime = quant_mat(1,n_prime)
        k1_prime = quant_mat(2,n_prime)
        j2_prime = quant_mat(3,n_prime)
        k2_prime = quant_mat(4,n_prime)
!
        ALLOCATE(xang(1:nang), wang(1:nang), STAT=istatus)
        IF (istatus /= 0) THEN
                PRINT*, 'Error: allocation failed for angular quadrature.'
                STOP
        ENDIF
!
        CALL gauleg(-1d0, 1d0, xang, wang, nang)
!
        dphi = 2d0 * pi / DBLE(nphi)
        wphi = dphi
        sum_complex = (0d0,0d0)
!
        DO q = q_left, q_right
                xm = 0.5d0 * (knots_x(q+1) + knots_x(q))
                xr = 0.5d0 * (knots_x(q+1) - knots_x(q))
!
                DO l = 1, ngqp_x
                        R = xm + xr * gq_root_x(l)
!
                        SELECT CASE (radial_mode)
                        CASE (radial_closed_closed)
                                radial_part = bsp_x(i_basis,R) * bsp_x(i_prime_basis,R) / &
                                        DSQRT(norm(i_basis) * norm(i_prime_basis))
                        CASE (radial_closed_open)
                                radial_part = bsp_x(i_basis,R) * u0(open_col,R) / &
                                        DSQRT(norm(i_basis))
                        CASE (radial_open_open_00)
                                radial_part = u0(open_row,R) * u0(open_col,R)
                        CASE (radial_open_open_10)
                                radial_part = CONJG(u0(open_row,R)) * u0(open_col,R)
                        END SELECT
!
                        DO ix1 = 1, nang
                                x1ang = xang(ix1)
!
                                DO iphi1 = 1, nphi
                                        phi1 = (DBLE(iphi1)-0.5d0) * dphi
!
                                        y_left_1 = Ylm_norm(j1, k1, x1ang, phi1)
                                        y_right_1 = Ylm_norm(j1_prime, k1_prime, x1ang, phi1)
!
                                        DO ix2 = 1, nang
                                                x2ang = xang(ix2)
!
                                                DO iphi2 = 1, nphi
                                                        phi2 = (DBLE(iphi2)-0.5d0) * dphi
!
                                                        y_left_2 = Ylm_norm(j2, k2, x2ang, phi2)
                                                        y_right_2 = Ylm_norm(j2_prime, k2_prime, x2ang, phi2)
!
                                                        Vang = potential_angular_expansion(R, x1ang, &
                                                                phi1, x2ang, phi2)
!
                                                        sum_complex = sum_complex + &
                                                                gq_weight_x(l) * xr * wang(ix1) * &
                                                                wang(ix2) * wphi * wphi * radial_part * &
                                                                CONJG(y_left_1 * y_left_2) * Vang * &
                                                                y_right_1 * y_right_2
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDDO
                ENDDO
        ENDDO
!
        Integral_potential_direct = sum_complex
!
        DEALLOCATE(xang, wang, STAT=istatus)
!
END FUNCTION Integral_potential_direct


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION potential_angular_expansion(R, x1ang, phi1, x2ang, phi2)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        REAL(8), INTENT(IN) :: R, x1ang, phi1, x2ang, phi2
!
        INTEGER :: ipot, lambda1, lambda2, m
        REAL(8) :: basis_factor
        COMPLEX(8) :: term_plus, term_minus
!
        potential_angular_expansion = (0d0,0d0)
!
        DO ipot = 1, n_pot
                lambda1 = pot_mat(1,ipot)
                lambda2 = pot_mat(2,ipot)
                m       = pot_mat(3,ipot)
!
                term_plus = C_lm_scaled(lambda1, m, x1ang, phi1) * &
                        C_lm_scaled(lambda2, -m, x2ang, phi2)
                term_minus = C_lm_scaled(lambda1, -m, x1ang, phi1) * &
                        C_lm_scaled(lambda2, m, x2ang, phi2)
!
                IF (m == 0) THEN
                        basis_factor = DSQRT(2d0*pi)
                ELSE
                        basis_factor = (-1d0)**m * 2d0 * DSQRT(pi)
                ENDIF
!
                potential_angular_expansion = potential_angular_expansion + &
                        A_pot(R, ipot) * basis_factor * 0.5d0 * (term_plus + term_minus)
        ENDDO
!
END FUNCTION potential_angular_expansion


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION Ylm_norm(l, m, xval, phi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: l, m
        REAL(8), INTENT(IN) :: xval, phi
!
        INTEGER :: mm
        REAL(8) :: norm_const, pval
        COMPLEX(8) :: im, ypos
!
        im = (0d0,1d0)
!
        IF (ABS(m) > l) THEN
                Ylm_norm = (0d0,0d0)
                RETURN
        ENDIF
!
        IF (m >= 0) THEN
                norm_const = DSQRT((2d0*DBLE(l)+1d0)/(4d0*pi) * &
                        factorial(l-m) / factorial(l+m))
                pval = plgndr(l,m,xval)
                Ylm_norm = norm_const * pval * EXP(im*DBLE(m)*phi)
        ELSE
                mm = -m
                norm_const = DSQRT((2d0*DBLE(l)+1d0)/(4d0*pi) * &
                        factorial(l-mm) / factorial(l+mm))
                pval = plgndr(l,mm,xval)
                ypos = norm_const * pval * EXP(im*DBLE(mm)*phi)
!
                IF (MOD(mm,2) == 0) THEN
                        Ylm_norm = CONJG(ypos)
                ELSE
                        Ylm_norm = -CONJG(ypos)
                ENDIF
        ENDIF
!
END FUNCTION Ylm_norm


!*********************************************************************************************
!*********************************************************************************************
!
COMPLEX(8) FUNCTION C_lm_scaled(l, m, xval, phi)
!
!*********************************************************************************************
!
        IMPLICIT NONE
!
        INTEGER, INTENT(IN) :: l, m
        REAL(8), INTENT(IN) :: xval, phi
!
        C_lm_scaled  =  Ylm_norm(l,m,xval,phi)
!
END FUNCTION C_lm_scaled


!=============================================================================================
!
END PROGRAM test_MV_direct
!
!*********************************************************************************************
!*********************************************************************************************
