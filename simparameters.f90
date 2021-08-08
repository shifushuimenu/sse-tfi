    module simparameters 
    use types 
    implicit none 
    
    ! **************************************************
    ! simulation parameters:           
    ! ----------------------
    ! Their values are overwritten from the input file.
    ! **************************************************

type t_Simparams
    real(dp) :: J_1 = +1.0_dp    
    real(dp) :: hx = 0.60_dp       ! transverse field (in units of J_1)
    real(dp) :: temp = 0.1_dp      ! temperature (in units of J_1)
    real(dp) :: beta               ! inverse temperature 
    integer :: nx, ny, n_sites             

    integer :: nmeas_step = 10000
    integer :: ntherm_step = 10000 
    integer :: Nbin = 100
    character(len=10) :: lattice_type = "triangular"
    logical :: ignore_Jmatrix = .FALSE.
    character(len=30) :: Jmatrix_file = "Jmatrix.txt"
    logical :: translat_invar = .TRUE.
    character(len=12) :: paramscan
    real(dp) :: scan_min = 0.0, scan_max = 0.0
    logical :: heavy_use = .FALSE.
    logical :: deterministic = .FALSE.
    real(dp) :: hz = 0.1_dp  
    character(len=30) :: hz_fields_file = "hz_fields.txt"
    logical :: ignore_hz_fields = .FALSE.
    real(dp) :: Rb                 ! Rydberg blockade radius
    ! (tuneable) hyperparameter of the algorithm with longitudinal field    
    real(dp) :: C_par_hyperparam = 0.1_dp  ! allowed range [0, +\infty]    
end type


! *******************************************************************
! The same as above, CAREFUL !
! These variables are used in a NAMELIST statement. 
! *******************************************************************
real(dp) :: J_1 = +1.0_dp    
real(dp) :: hx = 0.60_dp       ! transverse field (in units of J_1)
real(dp) :: temp = 0.1_dp      ! temperature (in units of J_1)
real(dp) :: beta               ! inverse temperature 
integer :: nx, ny, n_sites             

integer :: nmeas_step = 10000
integer :: ntherm_step = 10000 
integer :: Nbin = 100
character(len=10) :: lattice_type = "triangular"
logical :: ignore_Jmatrix = .FALSE.
character(len=30) :: Jmatrix_file = "Jmatrix.txt"
logical :: translat_invar = .TRUE.
character(len=12) :: paramscan
real(dp) :: scan_min = 0.0, scan_max = 0.0
logical :: heavy_use = .FALSE.
logical :: deterministic = .FALSE.
real(dp) :: hz = 0.1_dp  
character(len=30) :: hz_fields_file = "hz_fields.txt"
logical :: ignore_hz_fields = .FALSE.
! Rydberg blockade radius 
real(dp) :: Rb = 1.0 
! (tuneable) hyperparameter of the algorithm with longitudinal field    
real(dp) :: C_par_hyperparam = 0.1_dp  ! allowed range [0, +\infty]   
! *******************************************************************

end module 