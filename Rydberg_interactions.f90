! This module rewrites the Rydberg Hamiltonian in terms of a transverse-
! field Ising model. 

module Rydberg
    use types
    implicit none

contains 
    subroutine init_Rydberg_interactions(S, Sim, hz_fields, Jmatrix)
        use lattice 
        use simparameters, only: t_Simparams
        implicit none 

        type(Struct), intent(in)      :: S                ! Lattice structure
        type(t_Simparams), intent(inout) :: Sim           ! Structure of simulation parameters 
        real(dp), intent(out)         :: hz_fields(:)     ! site-dependent longitudinal fields 
        real(dp), intent(out)         :: Jmatrix(:,:)     ! Interaction matrix, not translationally invariant

        ! ... Local variables ...
        integer :: ir, jr
        real(dp) :: V0
        
        real(dp) :: ss

        select case(trim(S%lattice_type))

        case("chain")

            ! Rydberg -> TFI dictionary             
            V0 =  Sim%Rb**6 / 4.0_dp !Sim%Omega * Sim%Rb**6 / 4.0_dp
            Jmatrix(:,:) = 0.0_dp
            do ir = 1, S%nsites
            do jr = 1, S%nsites 
                if (ir .ne. jr) then 
                    ! NOT translationally invariant 
                    Jmatrix(ir,jr) =  0.5*(V0 / (abs(ir -jr)**6) + V0 / ((S%nsites - abs(ir- jr))**6)) ! V0 / (abs(ir -jr)**6)
                endif 
            enddo 
            enddo

            Sim%hx = Sim%Omega / 2.0_dp

            do ir = 1, S%nsites 
                ss = 0.0_dp
                do jr = 1, S%nsites 
                    if (ir /= jr) then 
                        ss = ss + Jmatrix(ir,jr)
                    endif 
                enddo
                hz_fields(ir) = - Sim%delta / 2.0_dp + ss
            enddo

            ss = 0.0_dp
            do ir = 1, S%nsites
            do jr = ir+1, S%nsites 
                ss = ss + Jmatrix(ir,jr)
            enddo 
            enddo 
            ! add this off-set to the final energy per site 
            Sim%Rydberg_energy_offset = - Sim%delta / 2.0_dp + ss / float(S%nsites)

            print*, "Rydberg energy offset ", Sim%Rydberg_energy_offset

        case default
            print*, "init_Rydberg_interactions(): Unknown lattice type"
            stop 

        end select 

    end subroutine 

end module 