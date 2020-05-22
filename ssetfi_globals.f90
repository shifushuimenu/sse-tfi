module ssetfi_globals
    use types 
    use SSE_configuration
    use diagonal_update 
    implicit none 
    
    integer, allocatable  :: sublattice(:)
    integer, allocatable  :: neigh(:,:)
    
    integer, parameter :: A_UPDATE=1, B_UPDATE=2, C_UPDATE=3
    type(t_Plaquette), allocatable :: plaquettes(:)
    
    real(dp), allocatable :: J_interaction_matrix(:,:)
    integer, allocatable  :: Jij_sign(:,:)
    type(t_ProbTable) :: probtable 
        
    type(t_Config) :: config 
    integer, allocatable :: spins(:)
    type(t_BondOperator), allocatable :: opstring(:)
    integer, allocatable :: vertexlink(:)
    logical, allocatable :: leg_visited(:)
    
    contains 
    
    subroutine deallocate_globals
    
        deallocate(sublattice)
        deallocate(neigh)
        deallocate(plaquettes)
        deallocate(J_interaction_matrix)
        deallocate(Jij_sign)
        deallocate(spins)
        deallocate(opstring)
        deallocate(vertexlink)
        deallocate(leg_visited)
        
    end subroutine deallocate_globals 
    
end module ssetfi_globals 