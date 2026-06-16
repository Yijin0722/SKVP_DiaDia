
subroutine jacobi_to_cartesian(R, r1, r2, x1, x2, phi, cc)
    implicit none

    real(8), intent(in) :: R, r1, r2
    real(8), intent(in) :: x1, x2, phi
    real(8), intent(out) :: cc(4,3)

    real(8) :: s1, s2
    real(8) :: r1vec(3), r2vec(3), Rvec(3)


    s1 = sqrt(max(0.0d0, 1.0d0 - x1*x1))
    s2 = sqrt(max(0.0d0, 1.0d0 - x2*x2))

    ! Molecule 1 bond vector
    r1vec(1) = r1 * s1
    r1vec(2) = 0.0d0
    r1vec(3) = r1 * x1

    ! Molecule 2 bond vector
    r2vec(1) = r2 * s2 * cos(phi)
    r2vec(2) = r2 * s2 * sin(phi)
    r2vec(3) = r2 * x2

    ! Center-of-mass separation vector
    Rvec(1) = 0.0d0
    Rvec(2) = 0.0d0
    Rvec(3) = R

    ! Atom 1 and 2: first H2 molecule
    cc(1,:) = -0.5d0 * r1vec(:)
    cc(2,:) =  0.5d0 * r1vec(:)

    ! Atom 3 and 4: second H2 molecule
    cc(3,:) = Rvec(:) - 0.5d0 * r2vec(:)
    cc(4,:) = Rvec(:) + 0.5d0 * r2vec(:)

end subroutine jacobi_to_cartesian


real(8) function potential_v(R, theta1, theta2, phi)
    implicit none

    real(8), intent(in) :: R, theta1, theta2, phi

    real(8), parameter :: r1 = 1.449d0
    real(8), parameter :: r2 = 1.449d0
    REAL(8), PARAMETER :: hartree_to_cm = 219474.6313705d0

    REAL(8) :: R_B, r1_B, r2_B
    real(8) :: x1, x2
    real(8) :: cc(4,3)
    real(8) :: dVdcc(4,3)

    external :: jacobi_to_cartesian
    external :: h4bmkp_cc

    R_B  = R  * 1.889726125d0
    r1_B = r1* 1.889726125d0
    r2_B = r2 * 1.889726125d0

    x1 = cos(theta1)
    x2 = cos(theta2)

    call jacobi_to_cartesian(R_B, r1_B, r2_B, x1, x2, phi, cc)
    call h4bmkp_cc(cc, potential_v, dVdcc, 0)

    !potential_v = potential_v * hartree_to_cm

end function potential_v
