module bspline_mod
    implicit none
    integer, parameter :: dp = kind(0.0d0)

    public :: b_spline
    public :: b_spline_deriv
    
    private :: find_span, get_basis_table

contains

    !===================================================================
    ! FUNCTION 1: b_spline
    ! Returns the value of B-spline 'i' of degree 'k' at 'x'
    !===================================================================
    function b_spline(i, k, x, knots) result(val)
        integer, intent(in) :: i, k
        real(dp), intent(in) :: x
        real(dp), intent(in) :: knots(:)
        real(dp) :: val
        integer :: span, n_knots
        real(dp), allocatable :: N_table(:,:)

        n_knots = size(knots)

        ! 1. Domain Check
        if (x < knots(i) .or. x > knots(i+k+1)) then
            val = 0.0_dp
            ! Handle Closed Interval at very end of domain
            if (x == knots(i+k+1) .and. x == knots(n_knots)) then
                ! Proceed to calc to handle boundary correctly
            else
                return 
            end if
        end if

        ! 2. Find Span
        span = find_span(n_knots, k, x, knots)
        
        ! 3. Check if 'i' is active in this span
        ! Active indices in span 's' are: s-k, s-k+1, ... s
        if (i < span - k .or. i > span) then
            val = 0.0_dp
            return
        end if

        ! 4. Compute Basis Table (Values only)
        allocate(N_table(0:k, 0:k))
        call get_basis_table(span, x, k, knots, N_table)

        ! 5. Retrieve Value
        ! N_table(d, j) stores N_{span-d+j, d}
        ! We want N_{i, k}. 
        ! Mapping: i = span - k + j  =>  j = i - span + k
        val = N_table(k, i - span + k)
        
        deallocate(N_table)
    end function b_spline


    !===================================================================
    ! FUNCTION 2: b_spline_deriv
    ! Returns the 'order'-th derivative of B-spline 'i'
    !===================================================================
    function b_spline_deriv(i, k, x, knots, order) result(val)
        integer, intent(in) :: i, k, order
        real(dp), intent(in) :: x
        real(dp), intent(in) :: knots(:)
        real(dp) :: val
        
        real(dp), allocatable :: N_table(:,:) ! Stores values of all degrees
        real(dp), allocatable :: D_table(:)   ! Stores derivatives
        real(dp), allocatable :: D_prev(:)    ! Temp for recursion
        integer :: span, n_knots
        integer :: d, j, local_idx
        real(dp) :: num1, num2, den1, den2, term1, term2
        integer :: p, idx_global

        ! Base checks
        if (order == 0) then
            val = b_spline(i, k, x, knots)
            return
        end if
        if (order > k) then
            val = 0.0_dp
            return
        end if

        n_knots = size(knots)

        ! 1. Domain Check
        if (x < knots(i) .or. x > knots(i+k+1)) then
            val = 0.0_dp
            if (x == knots(i+k+1) .and. x == knots(n_knots)) then
                 ! Proceed
            else
                 return 
            end if
        end if

        ! 2. Find Span
        span = find_span(n_knots, k, x, knots)

        ! 3. Check Active
        if (i < span - k .or. i > span) then
            val = 0.0_dp
            return
        end if

        ! 4. Compute Full Basis Table (Values of all degrees 0..k)
        ! This is essential because D^m N_{k} depends on N_{k-1}, then N_{k-2}...
        allocate(N_table(0:k, 0:k))
        call get_basis_table(span, x, k, knots, N_table)

        ! 5. Compute Derivatives Iteratively
        ! We transform the values in N_table into derivatives layer by layer.
        ! We only need to store the current column of derivatives to save memory,
        ! but we need the LOWER degree basis values from N_table.
        
        allocate(D_table(0:k)) ! Will hold derivatives of degree p
        allocate(D_prev(0:k))  ! Will hold derivatives of degree p-1
        
        ! Initialize D_prev with degree (k-order) values from N_table
        ! Why? because D^k of degree k reduces to D^0 of degree 0 eventually.
        ! Actually, easiest way: 
        ! Start with degree p = k - order. The 'order'-th derivative of this is NOT useful.
        ! Standard Formula: D(N_{i,p}) = p * ( N_{i,p-1}/d1 - N_{i+1,p-1}/d2 )
        
        ! Let's implement the recurrence exactly:
        ! We want D^order of N_{i,k}. 
        ! This requires D^{order-1} of N_{..., k-1}.
        ! ...
        ! This requires D^0 (Values) of N_{..., k-order}.
        
        ! Load the "Values" (0-th deriv) of the starting degree: p_start = k - order
        do j = 0, k - order
            D_prev(j) = N_table(k - order, j) 
            ! N_table(d, j) is N_{span-d+j, d}
        end do

        ! Now loop to increase degree from (k-order+1) up to k
        ! In each step, we compute the derivatives of the HIGHER degree 
        ! using the values of the LOWER degree.
        
        do p = k - order + 1, k
             ! p is the degree we are constructing.
             ! We are constructing D^{something} of degree p.
             ! Actually, simply applying the formula 'order' times is correct?
             ! No. The factor 'p' in the formula changes.
             
             ! Algorithm:
             ! 1. Start with values of degree (k-order). These are effectively 
             !    the "0-th derivatives of the basis functions of degree k-order".
             ! 2. Apply formula to get "1-st derivatives of basis functions of degree k-order+1".
             ! 3. ...
             ! 4. Get "order-th derivatives of basis functions of degree k".
             
             do j = 0, p ! The valid indices in the active span for degree p
                  ! We want to compute D_table(j) corresponding to N_{span-p+j, p}
                  ! Formula terms depend on N_{span-p+j, p-1} (index j in prev)
                  ! and N_{span-p+j+1, p-1} (index j+1 in prev??)
                  
                  ! Let's trace indices:
                  ! Global idx = span - p + j
                  ! Term 1: N_{idx, p-1}.   Global index in p-1: span - (p-1) + (j-1).
                  ! So Term 1 uses D_prev(j). (Wait, j ranges 0..p-1 in prev)
                  
                  ! Term 2: N_{idx+1, p-1}. Global index: span - (p-1) + j.
                  ! So Term 2 uses D_prev(j+1)? No, D_prev(j).
                  
                  ! Let's be precise:
                  ! D_prev has size 0..p-1.
                  ! We need D_table of size 0..p.
                  
                  idx_global = span - p + j
                  
                  ! Term 1: (p / (u_{i+p} - u_i)) * Prev_Val(i)
                  ! Denom 1 = knots(idx_global + p) - knots(idx_global)
                  den1 = knots(idx_global + p) - knots(idx_global)
                  
                  if (den1 /= 0.0_dp) then
                      ! Which index in D_prev corresponds to global 'idx_global'?
                      ! D_prev stores degree p-1. 
                      ! D_prev(m) corresponds to global: span - (p-1) + m
                      ! We want global = idx_global.
                      ! idx_global = span - p + j
                      ! So: span - p + 1 + m = span - p + j  =>  m = j - 1
                      if (j - 1 >= 0) then
                          num1 = D_prev(j - 1)
                          term1 = (real(p,dp) / den1) * num1
                      else
                          term1 = 0.0_dp
                      end if
                  else
                      term1 = 0.0_dp
                  end if

                  ! Term 2: (p / (u_{i+p+1} - u_{i+1})) * Prev_Val(i+1)
                  ! Denom 2 = knots(idx_global + p + 1) - knots(idx_global + 1)
                  den2 = knots(idx_global + p + 1) - knots(idx_global + 1)
                  
                  if (den2 /= 0.0_dp) then
                      ! We want global = idx_global + 1
                      ! span - p + 1 + m = span - p + j + 1 => m = j
                      if (j <= p - 1) then
                          num2 = D_prev(j)
                          term2 = (real(p,dp) / den2) * num2
                      else
                          term2 = 0.0_dp
                      end if
                  else
                      term2 = 0.0_dp
                  end if

                  D_table(j) = term1 - term2
             end do
             
             ! Copy D_table to D_prev for next iteration
             ! D_prev now holds derivatives for degree p
             do j = 0, p
                 D_prev(j) = D_table(j)
             end do
        end do

        ! 6. Extract Result
        ! D_prev now holds the 'order'-th derivatives of degree 'k'
        ! Mapping: i = span - k + j  => j = i - span + k
        val = D_prev(i - span + k)

        deallocate(N_table)
        deallocate(D_table)
        deallocate(D_prev)
    end function b_spline_deriv


    !===================================================================
    ! HELPER 1: Find Span (Same as before)
    !===================================================================
    function find_span(n_knots, k, x, knots) result(idx)
        integer, intent(in) :: n_knots, k
        real(dp), intent(in) :: x
        real(dp), intent(in) :: knots(:)
        integer :: idx, low, high, mid

        if (x >= knots(n_knots - k)) then
            idx = n_knots - k - 1
            return
        end if

        low = k + 1
        high = n_knots - k 
        idx = low

        do while (low <= high)
            mid = (low + high) / 2
            if (x >= knots(mid)) then
                idx = mid
                low = mid + 1
            else
                high = mid - 1
            end if
        end do
    end function find_span


    !===================================================================
    ! HELPER 2: Compute Full Basis Table (Lower triangle)
    ! N_table(d, j) stores the value of N_{span-d+j, d}(x)
    !===================================================================
    subroutine get_basis_table(i, x, k, knots, N_tab)
        integer, intent(in) :: i, k
        real(dp), intent(in) :: x
        real(dp), intent(in) :: knots(:)
        real(dp), intent(out) :: N_tab(0:k, 0:k)
        
        real(dp) :: left(k+1), right(k+1)
        real(dp) :: saved, temp
        integer :: j, r

        ! Degree 0
        N_tab(0, 0) = 1.0_dp
        
        do j = 1, k
            left(j)  = x - knots(i + 1 - j)
            right(j) = knots(i + j) - x
            saved    = 0.0_dp
            
            do r = 0, j-1
                ! The standard iterative update
                temp = N_tab(j-1, r) / (right(r+1) + left(j-r))
                
                ! Update the "current" degree values
                ! But we also need to SAVE them into our 2D table
                N_tab(j, r) = saved + right(r+1) * temp
                saved       = left(j-r) * temp
            end do
            N_tab(j, j) = saved
        end do
    end subroutine get_basis_table

end module bspline_mod
