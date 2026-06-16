subroutine spline(x, y, n, yp1, ypn, y2)
    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n)
    double precision, intent(in) :: yp1, ypn
    double precision, intent(out) :: y2(n)

    double precision :: u(n)
    double precision :: sig, p, qn, un
    integer :: i, k

    if (yp1 > 0.99d30) then
        y2(1) = 0.d0
        u(1) = 0.d0
    else
        y2(1) = -0.5d0
        u(1) = (3.d0/(x(2)-x(1))) * &
               ((y(2)-y(1))/(x(2)-x(1)) - yp1)
    end if

    do i = 2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.d0
        y2(i) = (sig - 1.d0)/p

        u(i) = (6.d0 * ( &
               (y(i+1)-y(i))/(x(i+1)-x(i)) - &
               (y(i)-y(i-1))/(x(i)-x(i-1)) ) / &
               (x(i+1)-x(i-1)) - sig*u(i-1)) / p
    end do

    if (ypn > 0.99d30) then
        qn = 0.d0
        un = 0.d0
    else
        qn = 0.5d0
        un = (3.d0/(x(n)-x(n-1))) * &
             (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
    end if

    y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.d0)

    do k = n-1, 1, -1
        y2(k) = y2(k)*y2(k+1) + u(k)
    end do

end subroutine spline


subroutine splint(xa, ya, y2a, n, x, y)
    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: xa(n), ya(n), y2a(n)
    double precision, intent(in) :: x
    double precision, intent(out) :: y

    integer :: klo, khi, k
    double precision :: h, a, b

    klo = 1
    khi = n

    do while (khi-klo > 1)
        k = (khi+klo)/2

        if (xa(k) > x) then
            khi = k
        else
            klo = k
        end if
    end do

    h = xa(khi) - xa(klo)

    if (h == 0.d0) then
        print *, "Error: bad xa input in splint"
        stop
    end if

    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h

    y = a*ya(klo) + b*ya(khi) + &
        ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) * h*h / 6.d0

end subroutine splint