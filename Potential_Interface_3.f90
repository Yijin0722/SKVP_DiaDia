MODULE Potential_Interface

        USE generateparameters, ONLY: pbasst
        USE AtomDiatomskvp

        IMPLICIT NONE
        SAVE
        
!       Potential index map
!       -------------------
        INTEGER, ALLOCATABLE, DIMENSION(:,:)           :: pot_mat
        INTEGER                                        :: n_pot

!
!       A_cache(ipot,q,l,side)
!       ----------------------
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:)      :: A_cache

!
!       Default angular expansion limits
!       --------------------------------
        INTEGER, PARAMETER :: lambda1_max_default = 2
        INTEGER, PARAMETER :: lambda2_max_default = 2
        INTEGER, PARAMETER :: m_max_default       = 2

!
!       Gaussian quadrature parameters for angular projection
!       -----------------------------------------------------
        INTEGER, PARAMETER :: nx_A   = 24
        INTEGER, PARAMETER :: nphi_A = 24

        REAL(dp), PARAMETER :: pi_A = 3.1415926535897932384626433832795_dp

        INTEGER, PARAMETER :: A_cache_version = 1
        INTEGER, PARAMETER :: cache_id_len = 128
        CHARACTER(LEN=*), PARAMETER :: A_cache_filename = 'A_cache.dat'
        !CHARACTER(LEN=cache_id_len), PARAMETER :: potential_id = 'BMKP_cm-1_v1'
        !CHARACTER(LEN=cache_id_len), PARAMETER :: potential_id = 'BMKP_hartree_v1'
        !CHARACTER(LEN=cache_id_len), PARAMETER :: potential_id = 'DJ_K_v1'
        CHARACTER(LEN=cache_id_len), PARAMETER :: potential_id = 'DJ_hartree_v1'
        REAL(dp), PARAMETER :: cache_tol = 1.0d-12
!       Interface to external potential_v function
!       ------------------------------------------
!       potential_v must return V(R,theta1,theta2,phi).
!
        INTERFACE
                REAL(8) FUNCTION potential_v(R, theta1, theta2, phi) RESULT(V)
                        
                        IMPLICIT NONE

                        integer, parameter :: dp = kind(1.0d0)
                        REAL(dp), INTENT(IN) :: R
                        REAL(dp), INTENT(IN) :: theta1
                        REAL(dp), INTENT(IN) :: theta2
                        REAL(dp), INTENT(IN) :: phi
                END FUNCTION potential_v
        END INTERFACE


CONTAINS


        SUBROUTINE build_potential_index

                IMPLICIT NONE

                INTEGER :: lambda1_max_in, lambda2_max_in, m_max_in
                INTEGER :: lambda1, lambda2, m
                INTEGER :: ipot
                INTEGER :: alloc_status

                lambda1_max_in = lambda1_max_default
                lambda2_max_in = lambda2_max_default
                m_max_in       = m_max_default

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

        END SUBROUTINE build_potential_index


        REAL(dp) FUNCTION A_pot(R, ipot)

!
!       This function returns A(R,ipot).
!
!       ipot -> (lambda1,lambda2,m)
!
!       New version:
!       A_pot no longer calls djpes_A_functions.
!       It directly computes A_all from potential_v by Gaussian quadrature,
!       then extracts the requested component.
!
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: R
                INTEGER,  INTENT(IN) :: ipot

                INTEGER :: lambda1, lambda2, m
                REAL(dp) :: A_all(0:lambda1_max_default, &
                                  0:lambda2_max_default, &
                                  0:m_max_default)

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

!
!       Read lambda1, lambda2, m from pot_mat
!       -------------------------------------
                lambda1 = pot_mat(1,ipot)
                lambda2 = pot_mat(2,ipot)
                m       = pot_mat(3,ipot)

!
!       Compute all A(lambda1,lambda2,m) at this R
!       ------------------------------------------
                CALL calculate_A_from_potential_v(R, A_all)

!
!       Return selected A component
!       ---------------------------
                A_pot = A_all(lambda1, lambda2, m)

        END FUNCTION A_pot

        SUBROUTINE read_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)

!
!       Read A_cache from A_cache.dat if compatible.
!       If the file is missing or incompatible, build_A_cache is called.
!
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: knots_x_in(:)
                REAL(dp), INTENT(IN) :: gq_root_x_in(:)
                INTEGER,  INTENT(IN) :: ngqp_x_in

                LOGICAL :: file_exists
                LOGICAL :: cache_ok

                INTEGER :: n_q, n_knots, n_gq
                INTEGER :: unit_cache
                INTEGER :: io_status
                INTEGER :: alloc_status

                INTEGER :: file_version
                INTEGER :: file_lambda1_max
                INTEGER :: file_lambda2_max
                INTEGER :: file_m_max
                INTEGER :: file_nx_A
                INTEGER :: file_nphi_A
                INTEGER :: file_n_pot
                INTEGER :: file_n_q
                INTEGER :: file_ngqp_x
                INTEGER :: file_n_knots
                INTEGER :: file_n_gq

                CHARACTER(LEN=cache_id_len) :: file_potential_id

                INTEGER,  ALLOCATABLE :: file_pot_mat(:,:)
                REAL(dp), ALLOCATABLE :: file_knots_x(:)
                REAL(dp), ALLOCATABLE :: file_gq_root_x(:)

!
!       Basic checks
!       ------------
                IF (.NOT. ALLOCATED(pot_mat)) THEN
                        PRINT*, 'Error: pot_mat is not allocated before read_A_cache.'
                        PRINT*, 'Call build_potential_index first.'
                        STOP
                ENDIF

                IF (n_pot <= 0) THEN
                        PRINT*, 'Error: n_pot <= 0 in read_A_cache.'
                        STOP
                ENDIF

                IF (MOD(ngqp_x_in, 2) /= 0) THEN
                        PRINT*, 'Error: ngqp_x_in must be even in read_A_cache.'
                        PRINT*, 'ngqp_x_in = ', ngqp_x_in
                        STOP
                ENDIF

                n_q     = SIZE(knots_x_in) - 1
                n_knots = SIZE(knots_x_in)
                n_gq    = SIZE(gq_root_x_in)

!
!       If cache file does not exist, build it.
!       ---------------------------------------
                INQUIRE(file = A_cache_filename, exist = file_exists)

                IF (.NOT. file_exists) THEN
                        PRINT*, 'A_cache file not found. Building fresh A_cache.'
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

!
!       Open cache file
!       ---------------
                OPEN(newunit = unit_cache, file = A_cache_filename, &
                     form = 'unformatted', access = 'stream', &
                     status = 'old', action = 'read', iostat = io_status)

                IF (io_status /= 0) THEN
                        PRINT*, 'Could not open A_cache file. Building fresh A_cache.'
                        PRINT*, 'iostat = ', io_status
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

!
!       Read header
!       -----------
                READ(unit_cache, iostat = io_status) file_version
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) file_potential_id
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) &
                        file_lambda1_max, file_lambda2_max, file_m_max
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) file_nx_A, file_nphi_A
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) &
                        file_n_pot, file_n_q, file_ngqp_x, file_n_knots, file_n_gq
                IF (io_status /= 0) GOTO 900

!
!       Check header compatibility
!       --------------------------
                cache_ok = .TRUE.

                IF (file_version /= A_cache_version) cache_ok = .FALSE.

                IF (TRIM(file_potential_id) /= TRIM(potential_id)) cache_ok = .FALSE.

                IF (file_lambda1_max /= lambda1_max_default) cache_ok = .FALSE.
                IF (file_lambda2_max /= lambda2_max_default) cache_ok = .FALSE.
                IF (file_m_max       /= m_max_default)       cache_ok = .FALSE.

                IF (file_nx_A   /= nx_A)   cache_ok = .FALSE.
                IF (file_nphi_A /= nphi_A) cache_ok = .FALSE.

                IF (file_n_pot   /= n_pot)      cache_ok = .FALSE.
                IF (file_n_q     /= n_q)        cache_ok = .FALSE.
                IF (file_ngqp_x  /= ngqp_x_in)  cache_ok = .FALSE.
                IF (file_n_knots /= n_knots)    cache_ok = .FALSE.
                IF (file_n_gq    /= n_gq)       cache_ok = .FALSE.

                IF (.NOT. cache_ok) THEN
                        CLOSE(unit_cache)
                        PRINT*, 'A_cache header does not match current calculation.'
                        PRINT*, 'Building fresh A_cache.'
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

!
!       Read and compare pot_mat, knots_x, and gq_root_x
!       ------------------------------------------------
                ALLOCATE(file_pot_mat(3,file_n_pot), &
                         file_knots_x(file_n_knots), &
                         file_gq_root_x(file_n_gq), &
                         STAT = alloc_status)

                IF (alloc_status /= 0) THEN
                        CLOSE(unit_cache)
                        PRINT*, 'Could not allocate temporary A_cache header arrays.'
                        PRINT*, 'Building fresh A_cache.'
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

                READ(unit_cache, iostat = io_status) file_pot_mat
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) file_knots_x
                IF (io_status /= 0) GOTO 900

                READ(unit_cache, iostat = io_status) file_gq_root_x
                IF (io_status /= 0) GOTO 900

                cache_ok = .TRUE.

                IF (ANY(file_pot_mat /= pot_mat)) cache_ok = .FALSE.

                IF (MAXVAL(ABS(file_knots_x - knots_x_in)) > cache_tol) THEN
                        cache_ok = .FALSE.
                ENDIF

                IF (MAXVAL(ABS(file_gq_root_x - gq_root_x_in)) > cache_tol) THEN
                        cache_ok = .FALSE.
                ENDIF

                IF (.NOT. cache_ok) THEN
                        CLOSE(unit_cache)
                        DEALLOCATE(file_pot_mat, file_knots_x, file_gq_root_x, &
                                   STAT = alloc_status)
                        PRINT*, 'A_cache grid or pot_mat does not match.'
                        PRINT*, 'Building fresh A_cache.'
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

!
!       Read A_cache
!       ------------
                IF (ALLOCATED(A_cache)) THEN
                        DEALLOCATE(A_cache, STAT = alloc_status)
                ENDIF

                ALLOCATE(A_cache(1:n_pot, 1:n_q, 1:ngqp_x_in/2, 1:2), &
                         STAT = alloc_status)

                IF (alloc_status /= 0) THEN
                        CLOSE(unit_cache)
                        DEALLOCATE(file_pot_mat, file_knots_x, file_gq_root_x, &
                                   STAT = alloc_status)
                        PRINT*, 'Could not allocate A_cache in read_A_cache.'
                        STOP
                ENDIF

                READ(unit_cache, iostat = io_status) A_cache

                CLOSE(unit_cache)

                DEALLOCATE(file_pot_mat, file_knots_x, file_gq_root_x, &
                           STAT = alloc_status)

                IF (io_status /= 0) THEN
                        PRINT*, 'Error while reading A_cache array.'
                        PRINT*, 'Building fresh A_cache.'
                        CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)
                        RETURN
                ENDIF

                PRINT*, 'A_cache loaded from file: ', A_cache_filename
                PRINT*, 'potential_id = ', potential_id

                RETURN

!
!       Read failure fallback
!       ---------------------
900             CONTINUE

                CLOSE(unit_cache)

                IF (ALLOCATED(file_pot_mat)) THEN
                        DEALLOCATE(file_pot_mat, STAT = alloc_status)
                ENDIF

                IF (ALLOCATED(file_knots_x)) THEN
                        DEALLOCATE(file_knots_x, STAT = alloc_status)
                ENDIF

                IF (ALLOCATED(file_gq_root_x)) THEN
                        DEALLOCATE(file_gq_root_x, STAT = alloc_status)
                ENDIF

                PRINT*, 'Could not read A_cache file correctly.'
                PRINT*, 'Building fresh A_cache.'
                PRINT*, 'iostat = ', io_status

                CALL build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)

        END SUBROUTINE read_A_cache



        SUBROUTINE build_A_cache(knots_x_in, gq_root_x_in, ngqp_x_in)

!
!       Fresh build of A_cache, then write it to A_cache.dat.
!
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: knots_x_in(:)
                REAL(dp), INTENT(IN) :: gq_root_x_in(:)
                INTEGER,  INTENT(IN) :: ngqp_x_in

                INTEGER :: q, l
                INTEGER :: n_q
                INTEGER :: n_knots, n_gq
                INTEGER :: alloc_status
                INTEGER :: unit_cache
                INTEGER :: io_status

                REAL(dp) :: xm, xr, dx
                REAL(dp) :: R_plus, R_minus

                REAL(dp) :: A_all(0:lambda1_max_default, &
                                  0:lambda2_max_default, &
                                  0:m_max_default)

                REAL(dp) :: t_cache_start, t_cache_end

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

                IF (MOD(ngqp_x_in, 2) /= 0) THEN
                        PRINT*, 'Error: ngqp_x_in must be even in build_A_cache.'
                        PRINT*, 'ngqp_x_in = ', ngqp_x_in
                        STOP
                ENDIF

!
!       Sizes
!       -----
                n_q     = SIZE(knots_x_in) - 1
                n_knots = SIZE(knots_x_in)
                n_gq    = SIZE(gq_root_x_in)

!
!       Allocate A_cache
!       ----------------
                IF (ALLOCATED(A_cache)) THEN
                        DEALLOCATE(A_cache, STAT = alloc_status)
                ENDIF

                ALLOCATE(A_cache(1:n_pot, 1:n_q, 1:ngqp_x_in/2, 1:2), &
                         STAT = alloc_status)

                IF (alloc_status /= 0) THEN
                        PRINT*, 'Error: allocation failed for A_cache in build_A_cache.'
                        STOP
                ENDIF

                A_cache = 0.0_dp

!
!       Build A_cache(ipot,q,l,side)
!       ----------------------------
!       side = 1 : R = xm + dx
!       side = 2 : R = xm - dx
!
                CALL CPU_TIME(t_cache_start)

                PRINT*, 'Building fresh A_cache from potential_v.'

                DO q = 1, n_q

                        xm = 0.5_dp * (knots_x_in(q+1) + knots_x_in(q))
                        xr = 0.5_dp * (knots_x_in(q+1) - knots_x_in(q))

                        DO l = 1, ngqp_x_in/2

                                dx = xr * gq_root_x_in(l+ngqp_x_in/2)

                                R_plus = xm + dx
                                CALL calculate_A_from_potential_v(R_plus, A_all)
                                CALL store_A_all_in_cache(A_all, q, l, 1)

                                R_minus = xm - dx
                                CALL calculate_A_from_potential_v(R_minus, A_all)
                                CALL store_A_all_in_cache(A_all, q, l, 2)

                        ENDDO

                ENDDO

!
!       Write A_cache to file
!       ---------------------
                OPEN(newunit = unit_cache, file = A_cache_filename, &
                     form = 'unformatted', access = 'stream', &
                     status = 'replace', action = 'write', iostat = io_status)

                IF (io_status /= 0) THEN
                        PRINT*, 'Warning: could not open A_cache file for writing.'
                        PRINT*, 'A_cache_filename = ', A_cache_filename
                        PRINT*, 'iostat = ', io_status
                        RETURN
                ENDIF

                WRITE(unit_cache, iostat = io_status) A_cache_version
                WRITE(unit_cache, iostat = io_status) potential_id

                WRITE(unit_cache, iostat = io_status) &
                        lambda1_max_default, lambda2_max_default, m_max_default

                WRITE(unit_cache, iostat = io_status) nx_A, nphi_A

                WRITE(unit_cache, iostat = io_status) &
                        n_pot, n_q, ngqp_x_in, n_knots, n_gq

                WRITE(unit_cache, iostat = io_status) pot_mat
                WRITE(unit_cache, iostat = io_status) knots_x_in
                WRITE(unit_cache, iostat = io_status) gq_root_x_in
                WRITE(unit_cache, iostat = io_status) A_cache

                CLOSE(unit_cache)

                CALL CPU_TIME(t_cache_end)

                PRINT*, 'time in fresh A_cache build/write: ', &
                        t_cache_end - t_cache_start, ' sec'

                IF (io_status /= 0) THEN
                        PRINT*, 'Warning: error while writing A_cache file.'
                        PRINT*, 'iostat = ', io_status
                ELSE
                        PRINT*, 'A_cache written to file: ', A_cache_filename
                        PRINT*, 'potential_id = ', potential_id
                ENDIF

        END SUBROUTINE build_A_cache

        SUBROUTINE store_A_all_in_cache(A_all, q, l, side)

!
!       Store full A(lambda1,lambda2,m) array into A_cache(ipot,q,l,side)
!       according to pot_mat.
!
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: A_all(0:lambda1_max_default, &
                                              0:lambda2_max_default, &
                                              0:m_max_default)

                INTEGER, INTENT(IN) :: q, l, side

                INTEGER :: ipot
                INTEGER :: lambda1, lambda2, m

                DO ipot = 1, n_pot

                        lambda1 = pot_mat(1,ipot)
                        lambda2 = pot_mat(2,ipot)
                        m       = pot_mat(3,ipot)

                        A_cache(ipot,q,l,side) = A_all(lambda1, lambda2, m)

                ENDDO

        END SUBROUTINE store_A_all_in_cache


        SUBROUTINE calculate_A_from_potential_v(R, A)

!
!       Compute all A(lambda1,lambda2,m) at one fixed R.
!
!       No theta1-theta2 exchange symmetry is assumed.
!
!       A(lambda1,lambda2,m)
!       =
!       integral V(R,theta1,theta2,phi)
!       P_lambda1^m(cos theta1)
!       P_lambda2^m(cos theta2)
!       Phi_m(phi)
!       d(cos theta1) d(cos theta2) dphi
!
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: R

                REAL(dp), INTENT(OUT) :: A(0:lambda1_max_default, &
                                           0:lambda2_max_default, &
                                           0:m_max_default)

                INTEGER :: i, j, k
                INTEGER :: lambda1, lambda2, m
                INTEGER :: lambda_max

                REAL(dp) :: xgrid(nx_A), wx(nx_A)
                REAL(dp) :: phigrid(nphi_A), wphi(nphi_A)
                REAL(dp) :: theta_grid(nx_A)

                REAL(dp) :: weight
                REAL(dp) :: basis_value

                REAL(dp) :: Plm(0:MAX(lambda1_max_default,lambda2_max_default), &
                                0:m_max_default, nx_A)

                REAL(dp) :: phi_part(0:m_max_default, nphi_A)

                REAL(dp) :: Vgrid(nx_A, nx_A, nphi_A)

                lambda_max = MAX(lambda1_max_default, lambda2_max_default)

                A        = 0.0_dp
                Plm      = 0.0_dp
                phi_part = 0.0_dp
                Vgrid    = 0.0_dp

!
!       Build Gaussian quadrature grids
!       -------------------------------
                CALL gauleg(-1.0_dp, 1.0_dp, xgrid, wx, nx_A)
                CALL gauleg(0.0_dp, 2.0_dp * pi_A, phigrid, wphi, nphi_A)

!
!       Convert x-grid to theta-grid
!       ----------------------------
!       x = cos(theta), theta = acos(x)
!
                DO i = 1, nx_A
                        theta_grid(i) = ACOS(MAX(-1.0_dp, MIN(1.0_dp, xgrid(i))))
                ENDDO

!
!       Precompute associated Legendre functions
!       ----------------------------------------
!       Plm(lambda,m,i) = P_lambda^m(x_i)
!
                DO lambda1 = 0, lambda_max
                DO m = 0, MIN(lambda1, m_max_default)
                DO i = 1, nx_A

                        Plm(lambda1,m,i) = assoc_legendre(lambda1, m, xgrid(i))

                ENDDO
                ENDDO
                ENDDO

!
!       Precompute phi basis functions
!       ------------------------------
                DO k = 1, nphi_A

                        phi_part(0,k) = 1.0_dp / SQRT(2.0_dp * pi_A)

                        DO m = 1, m_max_default
                                phi_part(m,k) = COS(REAL(m,dp) * phigrid(k)) / SQRT(pi_A)
                        ENDDO

                ENDDO

!
!       Precompute potential values on the full angular grid
!       ----------------------------------------------------
!       No theta1-theta2 symmetry is used.
!
                DO i = 1, nx_A
                DO j = 1, nx_A
                DO k = 1, nphi_A

                        Vgrid(i,j,k) = potential_v(R, theta_grid(i), theta_grid(j), phigrid(k))

                ENDDO
                ENDDO
                ENDDO

!
!       Gaussian quadrature projection
!       ------------------------------
                DO lambda1 = 0, lambda1_max_default, 2
                DO lambda2 = 0, lambda2_max_default, 2
                DO m = 0, MIN(lambda1, lambda2, m_max_default)

                        A(lambda1,lambda2,m) = 0.0_dp

                        DO i = 1, nx_A
                        DO j = 1, nx_A
                        DO k = 1, nphi_A

                                weight = wx(i) * wx(j) * wphi(k)

                                basis_value = Plm(lambda1,m,i) * &
                                              Plm(lambda2,m,j) * &
                                              phi_part(m,k)

                                A(lambda1,lambda2,m) = A(lambda1,lambda2,m) + &
                                        weight * Vgrid(i,j,k) * basis_value

                        ENDDO
                        ENDDO
                        ENDDO

                ENDDO
                ENDDO
                ENDDO

        
        
            
        END SUBROUTINE calculate_A_from_potential_v



        SUBROUTINE calculate_RMS_potential_expansion

!
!       This subroutine tests the potential expansion used in
!       Potential_Interface.
!
!       It uses the default expansion limits:
!
!               lambda1_max_default
!               lambda2_max_default
!               m_max_default
!
!       It compares:
!
!               V_expansion(R,theta1,theta2,phi)
!
!       reconstructed from A(lambda1,lambda2,m), against:
!
!               potential_v(R,theta1,theta2,phi)
!
                
                IMPLICIT NONE

                REAL(dp):: sum_squared_error
                REAL(dp) :: rms_error

!
!       Sampling parameters
!       -------------------
                INTEGER, PARAMETER :: nR     = 10
                INTEGER, PARAMETER :: ntheta = 10
                INTEGER, PARAMETER :: nphi   = 10

                INTEGER :: iR, itheta1, itheta2, iphi
                INTEGER :: lambda1, lambda2, m
                INTEGER :: npoints

                REAL(dp) :: R
                REAL(dp) :: theta1, theta2, phi
                REAL(dp) :: x1, x2
                REAL(dp) :: phi_part

                REAL(dp) :: V_expansion
                REAL(dp) :: V_potential
                REAL(dp) :: error

                REAL(dp) :: max_abs_error
                REAL(dp) :: max_error
                REAL(dp) :: max_R
                REAL(dp) :: max_theta1
                REAL(dp) :: max_theta2
                REAL(dp) :: max_phi
                REAL(dp) :: max_V_expansion
                REAL(dp) :: max_V_potential

                REAL(dp) :: A_all(0:lambda1_max_default, &
                                  0:lambda2_max_default, &
                                  0:m_max_default)

                REAL(dp) :: Rmin, Rmax

!
!       Initialize error accumulators
!       -----------------------------
                sum_squared_error = 0.0_dp
                rms_error         = 0.0_dp
                npoints           = 0

                max_abs_error     = -1.0_dp
                max_error         = 0.0_dp

                max_R             = 0.0_dp
                max_theta1        = 0.0_dp
                max_theta2        = 0.0_dp
                max_phi           = 0.0_dp

                max_V_expansion   = 0.0_dp
                max_V_potential   = 0.0_dp


        Rmin = REAL(pbasst(1)%pb_min, dp)
        Rmax = REAL(pbasst(1)%pb_max, dp)
!
!       Loop over test grid
!       -------------------
                DO iR = 1, nR

!
!               R grid from about 1 to 10.
!               The +0.03_dp avoids exact endpoint singularities.
!
                        R = Rmin + 0.03_dp + REAL(iR - 1, dp) * (Rmax - Rmin - 0.06_dp) / REAL(nR - 1, dp)
!
!               Compute all A(lambda1,lambda2,m) at this R using
!               the current Potential_Interface convention.
!
                        CALL calculate_A_from_potential_v(R, A_all)

                        DO itheta1 = 1, ntheta

                                theta1 = REAL(itheta1 - 1, dp) * pi_A / REAL(ntheta - 1, dp)
                                x1     = COS(theta1)

                        DO itheta2 = 1, ntheta

                                theta2 = REAL(itheta2 - 1, dp) * pi_A / REAL(ntheta - 1, dp)
                                x2     = COS(theta2)

                        DO iphi = 1, nphi

                                phi = REAL(iphi - 1, dp) * 2.0_dp * pi_A / REAL(nphi - 1, dp)

!
!                               Reconstruct V from A coefficients
!                               ----------------------------------
                                V_expansion = 0.0_dp

                                DO lambda1 = 0, lambda1_max_default, 2
                                DO lambda2 = 0, lambda2_max_default, 2
                                DO m = 0, MIN(lambda1, lambda2, m_max_default)

                                        IF (m == 0) THEN
                                                phi_part = 1.0_dp / SQRT(2.0_dp * pi_A)
                                        ELSE
                                                phi_part = COS(REAL(m, dp) * phi) / SQRT(pi_A)
                                        ENDIF

                                        V_expansion = V_expansion + &
                                                A_all(lambda1,lambda2,m) * &
                                                assoc_legendre(lambda1, m, x1) * &
                                                assoc_legendre(lambda2, m, x2) * &
                                                phi_part

                                ENDDO
                                ENDDO
                                ENDDO

!
!                               Original potential
!                               ------------------
                                V_potential = REAL(potential_v(R, theta1, theta2, phi), dp)

!
!                               Error
!                               -----
                                error = V_expansion - V_potential

                                sum_squared_error = sum_squared_error + error * error
                                npoints = npoints + 1

                                IF (ABS(error) > max_abs_error) THEN

                                        max_abs_error   = ABS(error)
                                        max_error       = error

                                        max_R           = R
                                        max_theta1      = theta1
                                        max_theta2      = theta2
                                        max_phi         = phi

                                        max_V_expansion = V_expansion
                                        max_V_potential = V_potential

                                ENDIF

                        ENDDO
                        ENDDO
                        ENDDO

                ENDDO

!
!       RMS error
!       ---------
                IF (npoints > 0) THEN
                        rms_error = SQRT(sum_squared_error / REAL(npoints, dp))
                ELSE
                        PRINT*, 'Error: npoints <= 0 in calculate_RMS_potential_expansion.'
                        STOP
                ENDIF

!
!       Print diagnostics
!       -----------------
                PRINT*, '----------------------------------------------'
                PRINT*, 'RMS test for Potential_Interface expansion'
                PRINT*, 'lambda1_max_default = ', lambda1_max_default
                PRINT*, 'lambda2_max_default = ', lambda2_max_default
                PRINT*, 'm_max_default       = ', m_max_default
                PRINT*, 'nR                  = ', nR
                PRINT*, 'ntheta              = ', ntheta
                PRINT*, 'nphi                = ', nphi
                PRINT*, 'npoints             = ', npoints
                PRINT*, 'sum_squared_error   = ', sum_squared_error
                PRINT*, 'rms_error           = ', rms_error
                PRINT*, 'max_abs_error       = ', max_abs_error
                PRINT*, 'max_error           = ', max_error
                PRINT*, 'max_error point:'
                PRINT*, 'R      = ', max_R
                PRINT*, 'theta1 = ', max_theta1
                PRINT*, 'theta2 = ', max_theta2
                PRINT*, 'phi    = ', max_phi
                PRINT*, 'V_expansion = ', max_V_expansion
                PRINT*, 'V_potential = ', max_V_potential
                PRINT*, '----------------------------------------------'

        END SUBROUTINE calculate_RMS_potential_expansion



END MODULE Potential_Interface