module legendre_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

contains

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

    function factorial(n) result(fact)
        implicit none

        integer, intent(in) :: n
        real(dp) :: fact
        integer :: i

        if (n < 0) then
            fact = 0.0_dp
            print *, "Error: n must be non-negative."
            return
        else if (n == 0 .or. n == 1) then
            fact = 1.0_dp
            return
        end if

        fact = 1.0_dp
        do i = 2, n
            fact = fact * real(i, dp)
        end do
    end function factorial

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

    subroutine export_assoc_legendre_csv(filename, l, m)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in) :: l, m
        real(dp) :: x, val
        integer :: iunit, i

        open(newunit=iunit, file=filename, status='replace', action='write', form='formatted')
        write(iunit, '(A)') 'l,m,x,P_l^m(x)'

        do i = -100, 100
            x = real(i, dp) / 100.0_dp
            val = assoc_legendre(l, m, x)
            write(iunit, '(I0,A,I0,A,F6.2,A,ES16.8)') l, ',', m, ',', x, ',', val
        end do

        close(iunit)
    end subroutine export_assoc_legendre_csv

end module legendre_mod


real function potential_v(R, theta1, theta2, phi)
    implicit none

    real, intent(in) :: R, theta1, theta2, phi
    real :: V000, V022, V224
    real :: g000, g022, g202, g224
    real :: spline_value

    real, dimension(27) :: r_table, v000_table, v022_table, v224_table
    real, dimension(27) :: y2_v000, y2_v022, y2_v224

    r_table = (/ 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, &
                 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 6.25, &
                 6.50, 6.75, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.00 /)

    v000_table = (/ 3765.99, 1474.07, 533.08, 158.79, 19.85, -24.30, -32.61, &
                    -29.03, -22.75, -16.97, -12.45, -9.10, -6.68, -4.94, &
                    -3.70, -2.80, -2.14, -1.66, -1.29, -1.02, -0.81, -0.53, &
                    -0.36, -0.24, -0.17, -0.12, -0.09 /)

    v022_table = (/ 118.71, 51.89, 20.66, 6.94, 1.42, -0.48, -0.93, -0.87, &
                    -0.68, -0.50, -0.35, -0.25, -0.18, -0.13, -9.61E-2, &
                    -7.18E-2, -5.45E-2, -4.19E-2, -3.25E-2, -2.56E-2, &
                    -2.03E-2, -1.31E-2, -8.77E-3, -6.01E-3, -4.22E-3, &
                    -3.02E-3, -2.27E-3 /)

    v224_table = (/ 55.44, 30.76, 18.17, 11.28, 7.30, 4.89, 3.38, 2.39, &
                    1.73, 1.28, 0.96, 0.73, 0.57, 0.44, 0.35, 0.28, 0.23, &
                    0.19, 0.15, 0.13, 0.11, 7.48E-2, 5.41E-2, 4.00E-2, &
                    3.00E-2, 2.29E-2, 1.82E-2 /)

    call natural_spline(r_table, v000_table, 27, y2_v000)
    call natural_spline(r_table, v022_table, 27, y2_v022)
    call natural_spline(r_table, v224_table, 27, y2_v224)

    g000 = 1.0
    g202 = 2.5 * (3.0 * cos(theta1)**2 - 1.0)
    g022 = 2.5 * (3.0 * cos(theta2)**2 - 1.0)
    g224 = 45.0 / (4.0 * sqrt(70.0)) * ( &
         2.0 * (3.0 * cos(theta1)**2 - 1.0) * (3.0 * cos(theta2)**2 - 1.0) - &
         16.0 * sin(theta1) * cos(theta1) * sin(theta2) * cos(theta2) * cos(phi) + &
         sin(theta1)**2 * sin(theta2)**2 * cos(2.0 * phi) )

    if (R > 10.0) then
        V000 = 0.0
        V022 = 0.0
        V224 = 0.0
    else if (R < 2.0) then
        V000 = v000_table(1)
        V022 = v022_table(1)
        V224 = v224_table(1)
    else
        V000 = spline_value(R, r_table, v000_table, y2_v000, 27)
        V022 = spline_value(R, r_table, v022_table, y2_v022, 27)
        V224 = spline_value(R, r_table, v224_table, y2_v224, 27)
    end if

    potential_v = V000 * g000 + V022 * g022 + V022 * g202 + V224 * g224

end function potential_v



subroutine natural_spline(x, y, n, y2)
    implicit none

    integer :: n, i, k
    real :: p, qn, sig, un
    real :: x(n), y(n), y2(n), u(n)

    y2(1) = 0.0
    u(1) = 0.0

    do i = 2, n - 1
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
        p = sig * y2(i - 1) + 2.0
        y2(i) = (sig - 1.0) / p
        u(i) = (6.0 * ( &
             (y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
             (y(i) - y(i - 1)) / (x(i) - x(i - 1)) ) / &
             (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do

    qn = 0.0
    un = 0.0
    y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1.0)

    do k = n - 1, 1, -1
        y2(k) = y2(k) * y2(k + 1) + u(k)
    end do

end subroutine natural_spline


real function spline_value(x, x_table, y_table, y2_table, n)
    implicit none

    integer :: n, klo, khi, k
    real :: x, h, a, b
    real :: x_table(n), y_table(n), y2_table(n)

    klo = 1
    khi = n

    do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (x_table(k) > x) then
            khi = k
        else
            klo = k
        end if
    end do

    h = x_table(khi) - x_table(klo)
    a = (x_table(khi) - x) / h
    b = (x - x_table(klo)) / h

    spline_value = a * y_table(klo) + b * y_table(khi) + &
                   ((a**3 - a) * y2_table(klo) + (b**3 - b) * y2_table(khi)) * &
                   h**2 / 6.0

end function spline_value


SUBROUTINE gauleg(x1,x2,x,w,n)
implicit none
INTEGER n
DOUBLE PRECISION x1,x2,x(n),w(n)
DOUBLE PRECISION EPS
PARAMETER (EPS=3.d-14)
INTEGER i,j,m
DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
m=(n+1)/2
xm=0.5d0*(x2+x1)
xl=0.5d0*(x2-x1)

do i=1,m
    z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1   continue
    p1=1.d0
    p2=0.d0
    do j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    end do

    pp=n*(z*p1-p2)/(z*z-1.d0)
    z1=z
    z=z1-p1/pp

    if(abs(z-z1).gt.EPS)goto 1

    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
    w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
    w(n+1-i)=w(i)
end do
return
END


subroutine triple_integral(Vval, lambda1, lambda2, m, integral_value)
    use legendre_mod
    implicit none

    integer, intent(in) :: lambda1, lambda2, m
    real(dp), intent(out) :: integral_value

    integer, parameter :: nx = 12
    integer, parameter :: nphi = 12
    real(dp), parameter :: pi = 3.14159265358979323846_dp

    real(dp) :: xgrid(nx), wx(nx)
    real(dp) :: phigrid(nphi), wphi(nphi)
    real(dp) :: theta_grid(nx)
    real(dp) :: leg1(nx), leg2(nx)
    real(dp) :: phi_part(nphi)
    real(dp) :: phi
    real(dp) :: weight
    real(dp) :: basis_pair
    integer :: i, j, k

    interface
        function Vval(theta1, theta2, phi) result(val)
            use legendre_mod
            implicit none
            real(dp), intent(in) :: theta1, theta2, phi
            real(dp) :: val
        end function Vval
    end interface

    external gauleg

    call gauleg(-1.0_dp, 1.0_dp, xgrid, wx, nx)
    call gauleg(0.0_dp, 2.0_dp * pi, phigrid, wphi, nphi)

    do i = 1, nx
        theta_grid(i) = acos(max(-1.0_dp, min(1.0_dp, xgrid(i))))
        leg1(i) = assoc_legendre(lambda1, m, xgrid(i))
        leg2(i) = assoc_legendre(lambda2, m, xgrid(i))
    end do

    do k = 1, nphi
        if (m == 0) then
            phi_part(k) = 1.0_dp / sqrt(2.0_dp * pi)
        else
            phi = phigrid(k)
            phi_part(k) = cos(real(m, dp) * phi) / sqrt(pi)
        end if
    end do

    integral_value = 0.0_dp

    do i = 1, nx
        do j = i, nx
            do k = 1, nphi
                phi = phigrid(k)
                weight = wx(i) * wx(j) * wphi(k)
                if (i == j) then
                    basis_pair = leg1(i) * leg2(j) * phi_part(k)
                else
                    basis_pair = (leg1(i) * leg2(j) + leg1(j) * leg2(i)) * phi_part(k)
                end if
                integral_value = integral_value + weight * &
                                 Vval(theta_grid(i), theta_grid(j), phi) * basis_pair
            end do
        end do
    end do

end subroutine triple_integral


subroutine Calculate_one_DJPES_Coeff(R, lambda1, lambda2, m, integral_value)
    use legendre_mod
    implicit none

    integer, intent(in) :: lambda1, lambda2, m
    real(dp), intent(out) :: integral_value
    real, intent(in) :: R

    interface
        subroutine triple_integral(Vval, lambda1, lambda2, m, integral_value)
            use legendre_mod
            implicit none

            integer, intent(in) :: lambda1, lambda2, m
            real(dp), intent(out) :: integral_value

            interface
                function Vval(theta1, theta2, phi) result(val)
                    use legendre_mod
                    implicit none

                    real(dp), intent(in) :: theta1, theta2, phi
                    real(dp) :: val
                end function Vval
            end interface
        end subroutine triple_integral

        function potential_v(R, theta1, theta2, phi) result(val)
            implicit none

            real, intent(in) :: R, theta1, theta2, phi
            real :: val
        end function potential_v
    end interface

    call triple_integral(DJPES_Vval, lambda1, lambda2, m, integral_value)

contains

    function DJPES_Vval(theta1, theta2, phi) result(val)
        use legendre_mod
        implicit none

        real(dp), intent(in) :: theta1, theta2, phi
        real(dp) :: val

        val = real(potential_v(real(R), real(theta1), real(theta2), real(phi)), dp)
    end function DJPES_Vval

end subroutine Calculate_one_DJPES_Coeff


module djpes_A_functions
    use legendre_mod
    implicit none

    integer, parameter :: nr_table = 27
    real(dp), parameter :: r_table(nr_table) = (/ &
        2.00_dp, 2.25_dp, 2.50_dp, 2.75_dp, 3.00_dp, 3.25_dp, 3.50_dp, &
        3.75_dp, 4.00_dp, 4.25_dp, 4.50_dp, 4.75_dp, 5.00_dp, 5.25_dp, &
        5.50_dp, 5.75_dp, 6.00_dp, 6.25_dp, 6.50_dp, 6.75_dp, 7.00_dp, &
        7.50_dp, 8.00_dp, 8.50_dp, 9.00_dp, 9.50_dp, 10.00_dp /)

contains

    subroutine calculate_DJPES_A_values(R, lambda1_max, lambda2_max, m_max, Aijk)
        implicit none

        real(dp), intent(in) :: R
        integer, intent(in) :: lambda1_max, lambda2_max, m_max
        real(dp), intent(out) :: Aijk(0:lambda1_max, 0:lambda2_max, 0:m_max)
        integer :: lambda1, lambda2, m

        Aijk = 0.0_dp

        do lambda1 = 0, lambda1_max
            do lambda2 = 0, lambda2_max
                do m = 0, min(lambda1, lambda2, m_max)
                    if (lambda2 < lambda1 .and. lambda1 <= lambda2_max) then
                        Aijk(lambda1, lambda2, m) = Aijk(lambda2, lambda1, m)
                    else
                        Aijk(lambda1, lambda2, m) = &
                            A_lambda1_lambda2_m(R, lambda1, lambda2, m)
                    end if
                end do
            end do
        end do
    end subroutine calculate_DJPES_A_values

    

    subroutine calculate_DJPES_A_value(lambda1, lambda2, m, R, A_value)
        implicit none

        integer, intent(in) :: lambda1, lambda2, m
        real(dp), intent(in) :: R
        real(dp), intent(out) :: A_value

        A_value = A_lambda1_lambda2_m(R, lambda1, lambda2, m)
    end subroutine calculate_DJPES_A_value


    subroutine calculate_array_A(R, lambda1_max, lambda2_max, m_max, A)
        implicit none

        real(dp), intent(in) :: R
        integer, intent(in) :: lambda1_max, lambda2_max, m_max
        real(dp), intent(out) :: A(0:lambda1_max, 0:lambda2_max, 0:m_max)
        integer :: lambda1, lambda2, m
        integer :: i
        real :: r_single
        real :: A_table(nr_table), A_y2(nr_table)

        interface
            subroutine natural_spline(x, y, n, y2)
                implicit none

                integer :: n
                real :: x(n), y(n), y2(n)
            end subroutine natural_spline

            real function spline_value(x, x_table, y_table, y2_table, n)
                implicit none

                integer :: n
                real :: x, x_table(n), y_table(n), y2_table(n)
            end function spline_value
        end interface

        A = 0.0_dp
        r_single = real(R)

        do lambda1 = 0, lambda1_max
            do lambda2 = 0, lambda2_max
                do m = 0, min(lambda1, lambda2, m_max)
                    if (lambda2 < lambda1 .and. lambda1 <= lambda2_max) then
                        A(lambda1, lambda2, m) = A(lambda2, lambda1, m)
                    else
                        do i = 1, nr_table
                            A_table(i) = real(A_lambda1_lambda2_m(r_table(i), lambda1, lambda2, m))
                        end do

                        call natural_spline(real(r_table), A_table, nr_table, A_y2)

                        if (R > r_table(nr_table)) then
                            A(lambda1, lambda2, m) = 0.0_dp
                        else if (R < r_table(1)) then
                            A(lambda1, lambda2, m) = real(A_table(1), dp)
                        else
                            A(lambda1, lambda2, m) = &
                                real(spline_value(r_single, real(r_table), A_table, A_y2, nr_table), dp)
                        end if
                    end if
                end do
            end do
        end do
    end subroutine calculate_array_A

    

    function A_lambda1_lambda2_m(R, lambda1, lambda2, m) result(A_value)
        implicit none

        real(dp), intent(in) :: R
        integer, intent(in) :: lambda1, lambda2, m
        real(dp) :: A_value

        call Calculate_one_DJPES_Coeff(real(R), lambda1, lambda2, m, A_value)
    end function A_lambda1_lambda2_m

end module djpes_A_functions



module dj_potential_expansion_mod
    use legendre_mod
    use djpes_A_functions
    implicit none

    real(dp), parameter :: kelvin_to_cm_inv = 1.0_dp / 1.438776877_dp

contains

    subroutine calculate_DJ_potential_expansion(R, theta1, theta2, phi, &
                                                lambda1_max, lambda2_max, m_max, V_expansion)
        implicit none

        real(dp), intent(in) :: R, theta1, theta2, phi
        integer, intent(in) :: lambda1_max, lambda2_max, m_max
        real(dp), intent(out) :: V_expansion

        real(dp) :: A(0:lambda1_max, 0:lambda2_max, 0:m_max)

        call calculate_array_A(R, lambda1_max, lambda2_max, m_max, A)

        call calculate_DJ_potential_expansion_from_A(theta1, theta2, phi, &
                                                     lambda1_max, lambda2_max, &
                                                     m_max, A, V_expansion)
    end subroutine calculate_DJ_potential_expansion

    subroutine calculate_DJ_potential_expansion_from_A(theta1, theta2, phi, &
                                                       lambda1_max, lambda2_max, &
                                                       m_max, A, V_expansion)
        implicit none

        real(dp), intent(in) :: theta1, theta2, phi
        integer, intent(in) :: lambda1_max, lambda2_max, m_max
        real(dp), intent(in) :: A(0:lambda1_max, 0:lambda2_max, 0:m_max)
        real(dp), intent(out) :: V_expansion

        integer :: lambda1, lambda2, m
        real(dp) :: phi_part
        real(dp), parameter :: pi = 3.14159265358979323846_dp

        V_expansion = 0.0_dp

        do lambda1 = 0, lambda1_max
            do lambda2 = 0, lambda2_max
                do m = 0, min(lambda1, lambda2, m_max)
                    if (m == 0) then
                        phi_part = 1.0_dp / sqrt(2.0_dp * pi)
                    else
                        phi_part = cos(real(m, dp) * phi) / sqrt(pi)
                    end if

                    V_expansion = V_expansion + A(lambda1, lambda2, m) * &
                                  assoc_legendre(lambda1, m, cos(theta1)) * &
                                  assoc_legendre(lambda2, m, cos(theta2)) * &
                                  phi_part
                end do
            end do
        end do

        V_expansion = V_expansion * kelvin_to_cm_inv
    end subroutine calculate_DJ_potential_expansion_from_A

    subroutine DJ_Potential_expansion(R, theta1, theta2, phi, V_expansion)
        implicit none

        real(dp), intent(in) :: R, theta1, theta2, phi
        integer, parameter :: lambda1_max = 2, lambda2_max = 2, m_max = 2
        real(dp), intent(out) :: V_expansion

        call calculate_DJ_potential_expansion(R, theta1, theta2, phi, &
                                              lambda1_max, lambda2_max, m_max, V_expansion)
    end subroutine DJ_Potential_expansion



end module dj_potential_expansion_mod
