module SSE_configuration
    implicit none  
  
    type :: tConfig 
    ! Configuration of the simulation cell 
    ! for a given SSE configuration 
      integer :: n_sites   ! Number of lattice sites 
      integer :: LL        ! Total number of SSE propagation steps (including idientities)
      integer :: n_exp     ! Expansion order, i.e. the total number of non-trivial operators
      integer :: n_legs    ! Total number of vertex legs: n_legs=6*n6leg+4*n4leg+2*n2leg
      integer :: n6leg     ! Number of triangular plaquette (i.e. 6-leg) vertices
      integer :: n4leg     ! Number of Ising (i.e. 4-leg) vertices 
      integer :: n2leg     ! Number of constant or spin-flip (i.e. 2-leg) vertices 
    end type  
  
    type :: tBondOperator
    ! operator in the SSE string; includes all operator types 
    integer :: i
    integer :: j
    ! the following variables are only relevant for triangular plaquette operators
    ! Plaquette operators are signalled by i<0, j<0.
    integer :: k  ! for triangular plaquette operators abs(i) = ir_A, abs(j) = ir_B, abs(k) = ir_C 
                  ! denote the three linearly stored sites on which the plaquette sits
    logical :: PRIVILEGED_LEG_IS_MAJORITY_LEG   ! true if the privileged leg sits on a majority spin configuration               
    integer :: privleg ! privleg \ in [1,2,3]: leg number of the privileged leg. 
          ! 1 = A-site, 2 = B-site, 3 = C-site
          ! If privleg = 1 (2,3), then 4 (5,6) is also a privileged leg. 
    end type tBondOperator
  
    ! for plaquette-based cluster update
    type :: tPlaquette
        ! triangular plaquette operator 
        integer :: Asite
        integer :: Bsite
        integer :: Csite
    end type tPlaquette    

end module 
  