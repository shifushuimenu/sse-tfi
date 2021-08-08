module ssetfi_globals
    use types 
    use SSE_configuration
    use diagonal_update 
    implicit none 
    
    integer, allocatable  :: sublattice(:)
    integer, allocatable  :: neigh(:,:)
    real(dp), allocatable :: rvec(:,:)   ! rvec(1:rdim, 1:nsites), real space positions 
    
    ! for lattice with triangular plaquettes 
    type(t_Plaquette), allocatable :: plaquettes(:)
    
    real(dp), allocatable :: J_interaction_matrix(:,:)
    real(dp), allocatable :: J_translat_invar(:,:)
    integer, allocatable  :: Jij_sign(:,:)
    integer, allocatable  :: hz_fields_sign(:)
    real(dp), allocatable :: hz_fields(:)
    type(t_ProbTable) :: probtable 
        
    type(t_Config) :: config 
    integer, allocatable :: spins(:)
    type(t_BondOperator), allocatable :: opstring(:)
    integer, allocatable :: vertexlink(:)
    logical, allocatable :: leg_visited(:)
    
    contains 
    
    subroutine deallocate_globals
    
        if( allocated(sublattice) ) deallocate(sublattice)
        if( allocated(neigh) ) deallocate(neigh)
        if( allocated(rvec) ) deallocate(rvec) 
        if( allocated(plaquettes) ) deallocate(plaquettes)
        if( allocated(J_interaction_matrix) ) deallocate(J_interaction_matrix)
        if( allocated(Jij_sign) ) deallocate(Jij_sign)
        if( allocated(hz_fields) ) deallocate(hz_fields)
        if( allocated(spins) ) deallocate(spins)
        if( allocated(opstring) ) deallocate(opstring)
        if( allocated(vertexlink) ) deallocate(vertexlink)
        if( allocated(leg_visited) ) deallocate(leg_visited)
        
    end subroutine deallocate_globals 
    
end module ssetfi_globals 