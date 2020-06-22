module SSE_configuration
    implicit none  
  
    ! labelling of legs 
    integer, parameter :: MAX_GHOSTLEGS = 6
    ! operator types 
    integer, parameter :: IDENTITY = 100
    integer, parameter :: ISING_BOND = 10
    integer, parameter :: TWO_LEG = 20, SPIN_FLIP = 21, CONSTANT = 22
    integer, parameter :: TRIANGULAR_PLAQUETTE = 30

    type :: t_Config 
    ! Configuration of the simulation cell 
    ! for a given SSE configuration 
      integer :: n_sites      ! Number of lattice sites 
      integer :: n_plaquettes ! Number of frustrated plaquettes which are tiling the lattice 
      integer :: LL           ! Total number of SSE propagation steps (including idientities)
      integer :: n_exp        ! Expansion order, i.e. the total number of non-trivial operators
      integer :: n_legs       ! Total number of physical (!) vertex legs: n_legs=6*n6leg+4*n4leg+2*n2leg
      integer :: n6leg        ! Number of triangular plaquette (i.e. 6-leg) vertices
      integer :: n4leg        ! Number of Ising (i.e. 4-leg) vertices 
      integer :: n2leg        ! Number of constant or spin-flip (i.e. 2-leg) vertices 
      integer :: n_ghostlegs  ! Total number of all legs, including 'ghostlegs': 
                              !     n_ghostlegs = n_exp * MAX_GHOSTLEGS
    end type  
  
    type :: t_BondOperator
    ! operator in the SSE string; includes all operator types (Ising, spin-flip, constant, triangular plaquette)
    integer :: i
    integer :: j
    ! the following variables are only relevant for triangular plaquette operators
    ! Plaquette operators are signalled by i<0, j<0.
    integer :: k  ! for triangular plaquette operators abs(i) = ir_A, abs(j) = ir_B, abs(k) = ir_C 
                  ! denote the three linearly stored sites on which the plaquette sits
    logical :: PRIVILEGED_LEG_IS_MAJORITY_LEG   ! true if the privileged leg sits on a majority spin configuration               
    ! The privileged site is always the 'A-site'. The lower leg on the 'A-site' always 
    ! has the lowest leg index around the vertex. To which sublattice an 'A-site' belongs
    ! is chosen during the diagonal update. 
    end type t_BondOperator
  
    ! for plaquette-based cluster update
    type :: t_Plaquette
        ! triangular plaquette operator 
        integer :: Asite
        integer :: Bsite
        integer :: Csite
    end type t_Plaquette    

end module 
  