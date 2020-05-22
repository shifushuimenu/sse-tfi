module types
    ! Here, the data type of all real variables can be set.
    implicit none 
    integer, parameter :: dp = kind(1.d0)

    real(dp), parameter :: ZERO = 0.0_dp  ! constant 0 
    real(dp), parameter :: ONE = 1.0_dp  ! constant 1 
    real(dp), parameter :: TWO = 2.0_dp  ! constant 2 
    real(dp), parameter :: HALF = 0.5_dp  ! constant 1/2

end module types 
