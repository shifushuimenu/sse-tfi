module Rydberg
    use types
    implicit none

contains 
    subroutine init_Rydberg_interactions(S, Rb, Jmatrix)
        use lattice 
        implicit none 

        type(Struct), intent(in)      :: S                ! Lattice structure
        real(dp), intent(in)          :: Rb               ! Rydberg blockade radius (in units of the Bravais lattice constant)
        real(dp), intent(out)         :: Jmatrix(:,:)     ! Interaction matrix, not translationally invariant

        ! ... Local variables ...
        integer :: ir, jr
        real(dp) :: V0
    
        select case(trim(S%lattice_type))

        case("chain")
            V0 = Rb**6
            do ir = 1, S%nsites
            do jr = 1, S%nsites 
                if (ir .ne. jr) then 
                    ! translationally invariant 
                    Jmatrix(ir,jr) = 0.5*(V0 / (abs(ir -jr)**6) + V0 / ((S%nsites - abs(ir- jr))**6))
                endif 
            enddo 
            enddo
    

        case default
            print*, "init_Rydberg_interactions(): Unknown lattice type"
            stop 

        end select 

    end subroutine 

end module 