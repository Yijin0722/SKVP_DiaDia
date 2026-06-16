
 MODULE Potential_Interface

        use skvpAtomDiatom
!
!
 IMPLICIT NONE



SAVE

!
!

 INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: pot_mat
 INTEGER                                        :: n_pot
 REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:)       :: A_cache


INTEGER, PARAMETER :: lambda1_max_default = 2
INTEGER, PARAMETER :: lambda2_max_default = 2
INTEGER, PARAMETER :: m_max_default       = 2

!
 CONTAINS


        SUBROUTINE build_potential_index
        !
        !*********************************************************************************************
        !
        ! This subroutine builds an index map for the potential expansion terms.
        !
        ! Each potential term is labeled by one compact index ipot:
        !
        !       ipot  <---->  (lambda1, lambda2, m)
        !
        ! The map is stored as:
        !
        !       pot_mat(1,ipot) = lambda1
        !       pot_mat(2,ipot) = lambda2
        !       pot_mat(3,ipot) = m
        !
        !=============================================================================================
        !
                IMPLICIT NONE
        !
                INTEGER :: lambda1_max_in, lambda2_max_in, m_max_in
                INTEGER :: lambda1, lambda2, m
                INTEGER :: ipot
                INTEGER :: alloc_status
        !
                lambda1_max_in = lambda1_max_default
                lambda2_max_in = lambda2_max_default
                m_max_in       = m_max_default
        !
        
        !
        !       Count number of potential terms
        !       --------------------------------
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
        !
        !       Allocate potential index map
        !       ----------------------------
                IF (ALLOCATED(pot_mat)) THEN
                        DEALLOCATE(pot_mat, STAT = alloc_status)
                ENDIF

                ALLOCATE(pot_mat(3,n_pot), STAT = alloc_status)

                IF (alloc_status /= 0) THEN
                        PRINT*, 'Error: allocation failed for pot_mat in build_potential_index.'
                        STOP
                ENDIF

                pot_mat = 0
        !
        !       Fill potential index map
        !       ------------------------
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
        !
        !       Print potential index map
        !       -------------------------
            ! PRINT*, ' '
            ! PRINT*, 'Potential index map:'
            ! PRINT*, 'n_pot = ', n_pot
            ! PRINT*, 'ipot        lambda1        lambda2        m'

            ! DO ipot = 1, n_pot
            !         PRINT*, ipot, pot_mat(1,ipot), pot_mat(2,ipot), pot_mat(3,ipot)
            ! ENDDO

            ! PRINT*, ' '
        !
        !=============================================================================================
        !
        END SUBROUTINE build_potential_index
        !
        !=============================================================================================
        !*********************************************************************************************
        !*********************************************************************************************

        !*********************************************************************************************
        !*********************************************************************************************
        !!*********************************************************************************************
        !*********************************************************************************************
        !
        REAL(8) FUNCTION A_pot(R, ipot)
        !
        !*********************************************************************************************
        !
        ! This function returns the DJPES radial coefficient A(R,ipot).
        !
        !       ipot  <---->  (lambda1, lambda2, m)
        !
        !       A_pot(R,ipot) = A_lambda1_lambda2_m(R, lambda1, lambda2, m)
        !
        ! where
        !
        !       lambda1 = pot_mat(1,ipot)
        !       lambda2 = pot_mat(2,ipot)
        !       m       = pot_mat(3,ipot)
        !
        !=============================================================================================
        ! 
                USE djpes_A_functions

        !
                IMPLICIT NONE
        !
                REAL(8), INTENT(IN) :: R
                INTEGER, INTENT(IN) :: ipot
        !
                INTEGER :: lambda1, lambda2, m
        !
        !       Safety checks
        !       -------------
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

!       Read lambda1, lambda2, m from pot_mat
        !       -------------------------------------
                lambda1 = pot_mat(1,ipot)
                lambda2 = pot_mat(2,ipot)
                m       = pot_mat(3,ipot)
        !
        !       Return A(R,ipot)
        !       ----------------
                A_pot = A_lambda1_lambda2_m(R, lambda1, lambda2, m)
        !
        !=============================================================================================
        !
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


         

end MODULE Potential_Interface
