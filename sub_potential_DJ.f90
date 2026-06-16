
REAL(8) FUNCTION spline_value(x, x_table, y_table, y2_table, n)
    implicit none

    INTEGER, INTENT(IN) :: n
    INTEGER :: klo, khi, k
    REAL(8), INTENT(IN) :: x
    REAL(8), INTENT(IN) :: x_table(n), y_table(n), y2_table(n)
    REAL(8) :: h, a, b

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
                   h**2 / 6.0d0

end function spline_value

subroutine natural_spline(x, y, n, y2)
    implicit none

    INTEGER, INTENT(IN) :: n
    INTEGER :: i, k
    REAL(8), INTENT(IN) :: x(n), y(n)
    REAL(8), INTENT(OUT) :: y2(n)
    REAL(8) :: p, qn, sig, un
    REAL(8) :: u(n)

    y2(1) = 0.0d0
    u(1) = 0.0d0

    do i = 2, n - 1
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
        p = sig * y2(i - 1) + 2.0d0
        y2(i) = (sig - 1.0d0) / p
        u(i) = (6.0d0 * ( &
             (y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
             (y(i) - y(i - 1)) / (x(i) - x(i - 1)) ) / &
             (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do

    qn = 0.0d0
    un = 0.0d0
    y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1.0d0)

    do k = n - 1, 1, -1
        y2(k) = y2(k) * y2(k + 1) + u(k)
    end do

end subroutine natural_spline

REAL(8) FUNCTION potential_v(R, theta1, theta2, phi)
    implicit none

    REAL(8), INTENT(IN) :: R, theta1, theta2, phi
    REAL(8) :: V000, V022, V224
    REAL(8) :: g000, g022, g202, g224
    REAL(8) :: spline_value

    REAL(8), DIMENSION(27) :: r_table, v000_table, v022_table, v224_table
    REAL(8), DIMENSION(27) :: y2_v000, y2_v022, y2_v224

    r_table = (/ 2.00d0, 2.25d0, 2.50d0, 2.75d0, 3.00d0, 3.25d0, 3.50d0, 3.75d0, 4.00d0, &
                 4.25d0, 4.50d0, 4.75d0, 5.00d0, 5.25d0, 5.50d0, 5.75d0, 6.00d0, 6.25d0, &
                 6.50d0, 6.75d0, 7.00d0, 7.50d0, 8.00d0, 8.50d0, 9.00d0, 9.50d0, 10.00d0 /)

    v000_table = (/ 3765.99d0, 1474.07d0, 533.08d0, 158.79d0, 19.85d0, -24.30d0, -32.61d0, &
                    -29.03d0, -22.75d0, -16.97d0, -12.45d0, -9.10d0, -6.68d0, -4.94d0, &
                    -3.70d0, -2.80d0, -2.14d0, -1.66d0, -1.29d0, -1.02d0, -0.81d0, -0.53d0, &
                    -0.36d0, -0.24d0, -0.17d0, -0.12d0, -0.09d0 /)

    v022_table = (/ 118.71d0, 51.89d0, 20.66d0, 6.94d0, 1.42d0, -0.48d0, -0.93d0, -0.87d0, &
                    -0.68d0, -0.50d0, -0.35d0, -0.25d0, -0.18d0, -0.13d0, -9.61d-2, &
                    -7.18d-2, -5.45d-2, -4.19d-2, -3.25d-2, -2.56d-2, &
                    -2.03d-2, -1.31d-2, -8.77d-3, -6.01d-3, -4.22d-3, &
                    -3.02d-3, -2.27d-3 /)

    v224_table = (/ 55.44d0, 30.76d0, 18.17d0, 11.28d0, 7.30d0, 4.89d0, 3.38d0, 2.39d0, &
                    1.73d0, 1.28d0, 0.96d0, 0.73d0, 0.57d0, 0.44d0, 0.35d0, 0.28d0, 0.23d0, &
                    0.19d0, 0.15d0, 0.13d0, 0.11d0, 7.48d-2, 5.41d-2, 4.00d-2, &
                    3.00d-2, 2.29d-2, 1.82d-2 /)

    call natural_spline(r_table, v000_table, 27, y2_v000)
    call natural_spline(r_table, v022_table, 27, y2_v022)
    call natural_spline(r_table, v224_table, 27, y2_v224)

    g000 = 1.0d0
    g202 = 2.5d0 * (3.0d0 * cos(theta1)**2 - 1.0d0)
    g022 = 2.5d0 * (3.0d0 * cos(theta2)**2 - 1.0d0)
    g224 = 45.0d0 / (4.0d0 * sqrt(70.0d0)) * ( &
         2.0d0 * (3.0d0 * cos(theta1)**2 - 1.0d0) * (3.0d0 * cos(theta2)**2 - 1.0d0) - &
         16.0d0 * sin(theta1) * cos(theta1) * sin(theta2) * cos(theta2) * cos(phi) + &
         sin(theta1)**2 * sin(theta2)**2 * cos(2.0d0 * phi) )

    if (R > 10.0d0) then
        V000 = 0.0d0
        V022 = 0.0d0
        V224 = 0.0d0
    else if (R < 2.0d0) then
        V000 = v000_table(1)
        V022 = v022_table(1)
        V224 = v224_table(1)
    else
        V000 = spline_value(R, r_table, v000_table, y2_v000, 27)
        V022 = spline_value(R, r_table, v022_table, y2_v022, 27)
        V224 = spline_value(R, r_table, v224_table, y2_v224, 27)
    end if

    potential_v = V000 * g000 + V022 * g022 + V022 * g202 + V224 * g224

    potential_v = (V000 * g000 + V022 * g022 + V022 * g202 + V224 * g224) / 320000d0

end function potential_v

