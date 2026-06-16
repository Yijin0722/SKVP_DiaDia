MODULE potential_BMKP
  IMPLICIT NONE
  SAVE

  INTEGER :: n_rpoints = 0
  INTEGER :: n_vterms  = 0
  LOGICAL :: potential_ready = .FALSE.

  REAL(8), ALLOCATABLE :: pot_R(:)
  REAL(8), ALLOCATABLE :: pot_coef(:,:)   ! pot_coef(i,t)
  REAL(8), ALLOCATABLE :: pot_y2(:,:)     ! spline helper

  INTEGER, ALLOCATABLE :: pot_lam1(:)
  INTEGER, ALLOCATABLE :: pot_lam2(:)
  INTEGER, ALLOCATABLE :: pot_m(:)

CONTAINS

  SUBROUTINE init_potential_BMKP(filename)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER :: unit, ios
    INTEGER :: i, t

    OPEN(newunit=unit, file=filename, status='old', action='read', iostat=ios)

    IF (ios /= 0) THEN
      PRINT*, 'Error: could not open potential file: ', TRIM(filename)
      STOP
    ENDIF

    READ(unit,*) n_rpoints, n_vterms

    ALLOCATE(pot_R(n_rpoints))
    ALLOCATE(pot_coef(n_rpoints, n_vterms))
    ALLOCATE(pot_y2(n_rpoints, n_vterms))
    ALLOCATE(pot_lam1(n_vterms))
    ALLOCATE(pot_lam2(n_vterms))
    ALLOCATE(pot_m(n_vterms))

    DO i = 1, n_rpoints
      READ(unit,*) pot_R(i)
    ENDDO

    DO t = 1, n_vterms
      READ(unit,*) pot_lam1(t), pot_lam2(t), pot_m(t)

      DO i = 1, n_rpoints
        READ(unit,*) pot_coef(i,t)
      ENDDO

      CALL spline(pot_R, pot_coef(:,t), n_rpoints, 1d30, 1d30, pot_y2(:,t))
    ENDDO

    CLOSE(unit)

    potential_ready = .TRUE.

    PRINT*, 'Loaded diatom-diatom potential coefficients'
    PRINT*, 'Number of R points: ', n_rpoints
    PRINT*, 'Number of terms:    ', n_vterms
  END SUBROUTINE init_potential_BMKP


  REAL(8) FUNCTION VtermBMKP(t, R)
    INTEGER, INTENT(IN) :: t
    REAL(8), INTENT(IN) :: R
    REAL(8) :: value

    IF (.NOT. potential_ready) THEN
      PRINT*, 'Error: potential_dd has not been initialized.'
      STOP
    ENDIF

    IF (t < 1 .OR. t > n_vterms) THEN
      PRINT*, 'Error: invalid potential term index t = ', t
      STOP
    ENDIF

    CALL splint(pot_R, pot_coef(:,t), pot_y2(:,t), n_rpoints, R, value)

    VtermBMKP = value
  END FUNCTION VtermBMKP


  REAL(8) FUNCTION VlamBMKP(lambda1, lambda2, m, R)
    INTEGER, INTENT(IN) :: lambda1, lambda2, m
    REAL(8), INTENT(IN) :: R
    INTEGER :: t

    IF (.NOT. potential_ready) THEN
      PRINT*, 'Error: potential_dd has not been initialized.'
      STOP
    ENDIF

    DO t = 1, n_vterms
      IF (pot_lam1(t) == lambda1 .AND. &
          pot_lam2(t) == lambda2 .AND. &
          pot_m(t)    == m) THEN

        VlamBMKP = VtermBMKP(t, R)
        RETURN
      ENDIF
    ENDDO

    VlamBMKP = 0d0
  END FUNCTION VlamBMKP


  LOGICAL FUNCTION has_VlamBMKP(lambda1, lambda2, m)
    INTEGER, INTENT(IN) :: lambda1, lambda2, m
    INTEGER :: t

    has_VlamBMKP = .FALSE.

    IF (.NOT. potential_ready) RETURN

    DO t = 1, n_vterms
      IF (pot_lam1(t) == lambda1 .AND. &
          pot_lam2(t) == lambda2 .AND. &
          pot_m(t)    == m) THEN

        has_VlamBMKP = .TRUE.
        RETURN
      ENDIF
    ENDDO
  END FUNCTION has_VlamBMKP


  SUBROUTINE print_potential_terms()
    INTEGER :: t

    IF (.NOT. potential_ready) THEN
      PRINT*, 'Error: potential_dd has not been initialized.'
      STOP
    ENDIF

    PRINT*, 'Potential terms:'
    DO t = 1, n_vterms
      PRINT*, t, pot_lam1(t), pot_lam2(t), pot_m(t)
    ENDDO
  END SUBROUTINE print_potential_terms

END MODULE potential_BMKP

