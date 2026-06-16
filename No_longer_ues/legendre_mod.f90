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
