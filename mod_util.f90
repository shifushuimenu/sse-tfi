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
