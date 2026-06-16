MODULE Potential_Interface

        USE potential_BMKP
        IMPLICIT NONE
        SAVE

        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pot_mat
        INTEGER :: n_pot = 0
        REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:)       :: A_cache


        INTEGER, PARAMETER :: lambda1_max_default = 2
        INTEGER, PARAMETER :: lambda2_max_default = 2
        INTEGER, PARAMETER :: m_max_default       = 2


        CHARACTER(LEN=16)  :: potential_backend = 'DJ'
        CHARACTER(LEN=512) :: bmkp_filename = 'coefficients.dat'

CONTAINS

        SUBROUTINE set_potential_backend(name)
                IMPLICIT NONE
                CHARACTER(*), INTENT(IN) :: name

                potential_backend = ADJUSTL(name)
        END SUBROUTINE set_potential_backend


        SUBROUTINE set_bmkp_filename(filename)
                IMPLICIT NONE
                CHARACTER(*), INTENT(IN) :: filename

                bmkp_filename = ADJUSTL(filename)
        END SUBROUTINE set_bmkp_filename


        SUBROUTINE build_potential_index
                IMPLICIT NONE

                INTEGER :: lambda1_max_in, lambda2_max_in, m_max_in
                INTEGER :: lambda1, lambda2, m
                INTEGER :: ipot
                INTEGER :: alloc_status

                lambda1_max_in = lambda1_max_default
                lambda2_max_in = lambda2_max_default
                m_max_in       = m_max_default

                IF (TRIM(potential_backend) == 'BMKP') THEN
                        IF (.NOT. potential_ready) THEN
                                CALL init_potential_BMKP(TRIM(bmkp_filename))
                        ENDIF

                        n_pot = n_vterms

                        IF (ALLOCATED(pot_mat)) THEN
                                DEALLOCATE(pot_mat, STAT = alloc_status)
                        ENDIF

                        ALLOCATE(pot_mat(3,n_pot), STAT = alloc_status)

                        IF (alloc_status /= 0) THEN
                                PRINT*, 'Error: allocation failed for BMKP pot_mat.'
                                STOP
                        ENDIF

                        DO ipot = 1, n_pot
                                pot_mat(1,ipot) = pot_lam1(ipot)
                                pot_mat(2,ipot) = pot_lam2(ipot)
                                pot_mat(3,ipot) = pot_m(ipot)
                        ENDDO

                        RETURN
                ENDIF

                n_pot = 0

                DO lambda1 = 0, lambda1_max_in, 2
                DO lambda2 = 0, lambda2_max_in, 2
                DO m = 0, m_max_in
                        IF (m <= lambda1 .AND. m <= lambda2) THEN
                                n_pot = n_pot + 1
                        ENDIF
                ENDDO
                ENDDO
                ENDDO

                IF (ALLOCATED(pot_mat)) THEN
                        DEALLOCATE(pot_mat, STAT = alloc_status)
                ENDIF

                ALLOCATE(pot_mat(3,n_pot), STAT = alloc_status)

                IF (alloc_status /= 0) THEN
                        PRINT*, 'Error: allocation failed for pot_mat in build_potential_index.'
                        STOP
                ENDIF

                pot_mat = 0
                ipot = 0

                DO lambda1 = 0, lambda1_max_in, 2
                DO lambda2 = 0, lambda2_max_in, 2
                DO m = 0, m_max_in
                        IF (m <= lambda1 .AND. m <= lambda2) THEN
                                ipot = ipot + 1
                                pot_mat(1,ipot) = lambda1
                                pot_mat(2,ipot) = lambda2
                                pot_mat(3,ipot) = m
                        ENDIF
                ENDDO
                ENDDO
                ENDDO

        END SUBROUTINE build_potential_index


        REAL(8) FUNCTION A_pot(R, ipot)
                USE djpes_A_functions
                IMPLICIT NONE

                REAL(8), INTENT(IN) :: R
                INTEGER, INTENT(IN) :: ipot

                INTEGER :: lambda1, lambda2, m

                IF (.NOT. ALLOCATED(pot_mat)) THEN
                        PRINT*, 'Error: pot_mat is not allocated. Call build_potential_index first.'
                        STOP
                ENDIF

                IF (ipot < 1 .OR. ipot > n_pot) THEN
                        PRINT*, 'Error: ipot is outside the allowed range.'
                        PRINT*, 'ipot = ', ipot
                        PRINT*, 'n_pot = ', n_pot
                        STOP
                ENDIF

                lambda1 = pot_mat(1,ipot)
                lambda2 = pot_mat(2,ipot)
                m       = pot_mat(3,ipot)

                IF (TRIM(potential_backend) == 'BMKP') THEN
                        IF (.NOT. potential_ready) THEN
                                PRINT*, 'Error: BMKP potential has not been initialized.'
                                STOP
                        ENDIF

                        IF (R < pot_R(1) .OR. R > pot_R(n_rpoints)) THEN
                                A_pot = 0d0
                                RETURN
                        ENDIF

                        A_pot = VtermBMKP(ipot, R)
                        RETURN
                ENDIF

                A_pot = A_lambda1_lambda2_m(R, lambda1, lambda2, m)

        END FUNCTION A_pot

                     
  SUBROUTINE build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)

        IMPLICIT NONE

        REAL(8), INTENT(IN) :: knots_x_in(:)
        REAL(8), INTENT(IN) :: gq_root_x_in(:)
        INTEGER, INTENT(IN) :: ngqp_x_in

        INTEGER :: ipot, q, l
        INTEGER :: n_q
        INTEGER :: alloc_status
        REAL(8) :: xm, xr, dx

!
!       Safety checks
!       -------------
        IF (.NOT. ALLOCATED(pot_mat)) THEN
                PRINT*, 'Error: pot_mat is not allocated before build_A_cache.'
                PRINT*, 'Call build_potential_index first.'
                STOP
        ENDIF

        IF (n_pot <= 0) THEN
                PRINT*, 'Error: n_pot <= 0 in build_A_cache.'
                STOP
        ENDIF

!
!       Number of knot intervals
!       ------------------------
        n_q = SIZE(knots_x_in) - 1

!
!       Allocate A_cache
!       ----------------
        IF (ALLOCATED(A_cache)) THEN
                DEALLOCATE(A_cache, STAT = alloc_status)
        ENDIF

        ALLOCATE(A_cache(1:n_pot, 1:n_q, 1:ngqp_x_in/2, 1:2), STAT = alloc_status)

        IF (alloc_status /= 0) THEN
                PRINT*, 'Error: allocation failed for A_cache in build_A_cache.'
                STOP
        ENDIF

        A_cache = 0d0

!
!       Build A_cache(ipot,q,l,side)
!       ----------------------------
        DO ipot = 1, n_pot
        DO q = 1, n_q

                xm = 0.5d0 * (knots_x_in(q+1) + knots_x_in(q))
                xr = 0.5d0 * (knots_x_in(q+1) - knots_x_in(q))

                DO l = 1, ngqp_x_in/2

                        dx = xr * gq_root_x_in(l+ngqp_x_in/2)

                        A_cache(ipot,q,l,1) = A_pot(xm+dx, ipot)
                        A_cache(ipot,q,l,2) = A_pot(xm-dx, ipot)

                ENDDO

        ENDDO
        ENDDO

END SUBROUTINE build_A_cache

END MODULE Potential_Interface
