module util
    use types 
    implicit none 

    contains 

function spins2binrep(spins) result(binrep_integer)
!
! Purpose:
! ========
! Convert an array of spin values \in [-1,1] to 
! an integer whose binary representation corresponds to 
! the spin array (with -1 values replaced by 0).
! 
! Arguments:
! ==========
    integer, intent(in) :: spins(:)
    integer :: binrep_integer

! ... Local variables ...
    integer :: n
    integer :: i

! ... Executable ...
    n = size(spins, dim=1)
    binrep_integer = 0
    do i=1,n   
        if (spins(i) == +1) then 
            binrep_integer = binrep_integer + 2**(i-1)
        endif 
    enddo

end function 
    
subroutine rotate(A, l)
    ! Purpose:
    ! --------
    !    Cyclically left-shift the elements of integer array 'A'
    !    by 'l' positions. 
    ! Arguments:
    ! ----------
        integer, intent(inout) :: A(:)  
        integer, intent(in)  :: l
    ! ... Local variables ...
        integer :: ii, n, s, temp(size(A,1))
    ! ... Executable ...
        n = size(A,1)
        s = mod(l, n)
        temp(1:s) = A(1:s)
        do ii = 1, n-s
            A(ii) = A(ii+s)
        enddo 
        A(n-s+1:n) = temp(1:s)
    end subroutine rotate

    subroutine random_permutation(A)
    ! Purpose:
    ! --------
    !   Randomly permute the elements in the integer array A(:)
    !   (in place).
    !
    ! Arguments:
    ! ----------
        integer, intent(inout) :: A(:)
    ! ... Local variables ...
        integer :: n, k, l
        real(dp) :: eta
    ! ... Executable ...
        n = size(A,1)
        do k = n, 2, -1
            call random_number(eta)
            l = ceiling(eta * (k-1))
            call rotate(A(1:k), l)
        enddo 
    end subroutine 

    subroutine init_RNG(MPI_rank, DETERMINISTIC)
    ! Purpose:
    ! --------
    !    Initialize the standard pseudo-random number generator 
    !    with a seed obtained from the system time (at the millisecond level)
    !    and the MPI rank.
    !    Call the random number generator:
    !
    !        integer :: eta
    !        call random_number(eta)
    ! Arguments:
    ! ----------
        integer, intent(in) :: MPI_rank
        logical, intent(in) :: DETERMINISTIC
    ! ... Local variables ...
        integer :: n, values(1:8)
        integer, allocatable :: seed(:)

        call random_seed(size=n)
        allocate(seed(n))

        if(DETERMINISTIC) then 
            seed(1:n) = 42
            call random_seed(put=seed)
        else
            call date_and_time(values=values)
            seed(1:n) = values(8) + MPI_rank
            call random_seed(put=seed)
        endif 

    end subroutine init_RNG


subroutine JackKnife(n, avg, err, x, y, sgn, sum_sgn)
! 
! Purpose:
! ========
! This subroutine implements the JackKnife method. 
!    
! Reference: 
! ==========
!    A.P. Young: "Everything you wanted to know about Data Analysis 
!                  and Fitting but were afraid to ask"
!    arXiv:1210.3781v3  (2014) 
!
! Arguments:
! ==========
integer, intent(in) :: n
real(dp), intent(out) :: avg
real(dp), intent(out) :: err
real(dp), intent(in) :: x(n)
real(dp), intent(out) :: y(n)
real(dp), intent(out) :: sgn(n)
real(dp), intent(out) :: sum_sgn

! ... Parameter ..
real(dp), parameter :: TOL = 1.0D-12

! ... Local variable ...
real(dp) :: sum_x 

! ... Executable ...
sum_x = sum(x)
avg = sum_x / n

sgn = (sum_x-x(1:n))
sum_sgn = sum_x

! compute y (Jackknife sample)
y = sgn/(n-1)

! avg_y = avg_x
y = y - avg 
y = y * y
err = (sum(y)*(n-1))/n
err = sqrt(err)

! If error is small enough, the regard it as 0.
if (err .lt. TOL*abs(avg)) then 
    err = ZERO
end if 
        
end subroutine JackKnife


subroutine GetError(n, err, avg, list)
! 
! Purpose:
! ========
!   This subroutine computes error of the measurements.
! 
! Arguments:
! ==========
integer, intent(in) :: n
real(dp), intent(out) :: err
real(dp), intent(in) :: avg 
real(dp), intent(inout) :: list(n)

! ... Parameter ...
real(dp), parameter :: TOL = 1.0D-12

! ... Local variable ...
real(dp) :: tmp 

! ... Executable ...

! compute average 
tmp = sum(list)/n

! standard deviation 
list = list - tmp 
list = list * list 
err = (sum(list)*(n-1))/n
err = sqrt(err)

! If error is small enough, then regard it as 0.
if (err .lt. TOL*abs(avg)) then 
    err = ZERO
end if 

end subroutine GetError  

end module 
