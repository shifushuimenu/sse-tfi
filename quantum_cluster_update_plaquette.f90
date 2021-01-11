module cluster_update

#ifdef DEBUG_CLUSTER_UPDATE
use test_helper, only: output_SSE_config
#endif
use types, only: dp, TWO
use SSE_configuration 
implicit none 
private

public quantum_cluster_update_plaquette

! convenient labels 
integer, parameter :: UP=+1, DOWN=-1
! for triangular-plaquette based update 
integer, parameter :: A_LEG=1, B_LEG=2, C_LEG=3

! Total weight ratio for all flipped clusters
! due to the exchange field as a result of a 
! longitudinal field. 
! (module-level variable)
real(dp) :: wratio_exchange_field

contains 

pure function gleg_to_ir(op, gleg) result(ir)
! ***************************************************************
! Purpose: 
! --------
! Implements the mapping from the global leg number 'gleg'
! of the operator 'op' to the linearly stored lattice site 'ir'
! on which the leg sits. 
! 
! This function is only called if a winding macrospin 
! is detected to update the initial spin configuration. 
! ***************************************************************

    implicit none 
    type(t_BondOperator), intent(in) :: op   ! operator at propagation step ip 
    integer, intent(in) :: gleg ! global leg number 
    integer :: ir               ! linearly stored index of a lattice site 

    integer :: vleg ! leg number around a vertex 
    integer :: ip   ! SSE propagation step
    ! Upper and lower legs of a vertex, which sit on the same site,
    ! have the same value of 'vleg_mod'        
    integer :: vleg_mod   
    
    ip = (gleg-1) / MAX_GHOSTLEGS + 1
    vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1
    vleg_mod = mod(vleg-1, MAX_GHOSTLEGS_HALF) + 1
    
    select case( op%optype )
        case(TRIANGULAR_PLAQUETTE)
            if( vleg_mod == A_LEG ) then
                ir = -op%i
            elseif( vleg_mod == B_LEG ) then 
                ir = -op%j
            elseif( vleg_mod == C_LEG ) then 
                ir = -op%k
            endif 
        case(ISING_BOND)
            if( vleg_mod == 1 ) then 
                ir = op%i
            elseif( vleg_mod == 2) then 
                ir = op%j
            endif 
        case(CONSTANT)       
            ir = op%i
        case(SPIN_FLIP)
            ir = op%i
        case(LONGITUDINAL)
            ir = op%i
        case default
            ! No checking for strange input, since stop and print statements 
            ! are not allowed in a pure function. 
    end select 
        

end function


pure function leg_direction(gleg) result(dir)
    implicit none 
    ! Arguments:
    ! ==========
    ! gleg is the global leg number
    integer, intent(in) :: gleg     
    integer :: dir

    ! ... Local variables ...
    integer :: ip 
    ! vleg is the leg number around a vertex
    integer :: vleg     

    ! ... Executable ...
    ip = (gleg-1) / MAX_GHOSTLEGS + 1
    vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1

    if( vleg > MAX_GHOSTLEGS_HALF ) then 
        dir = UP
    else 
        dir = DOWN
    endif 

end function leg_direction 
        

subroutine quantum_cluster_update_plaquette( &
    spins, opstring, vertexlink, leg_visited, config, &
    hz_fields, C_par_hyperparam )
! Purpose:
! --------
! Swendsen-Wang variant of the quantum cluster update.
! 
! For longitudinal field `hz` = 0:     
! Build all clusters, which is a deterministic process,
! and flip each cluster with probability 1/2. 
!
! For `hz` != 0:     
! Flip each cluster with probability 1/2 and accept the final
! configuration with Metropolis acceptance rate according to 
! the exchange fields of all flipped cluster. 

use class_Stack
implicit none
! Arguments:
! ----------
integer, intent(inout)              :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
integer, intent(in)                 :: vertexlink(:)
! Note: leg_visited(:) must have been initialized during the construction 
! of the linked list so as to mark 'ghostlegs' as visited. 
logical, intent(inout)              :: leg_visited(:) 
type(t_Config), intent(in)          :: config  
real(dp), intent(in)                :: hz_fields(:)      ! site-dependent longitudinal fields
real(dp), intent(in)                :: C_par_hyperparam  ! hyperparameter for algorithm with longitudinal field

    
! ... Local variables ...
type(t_BondOperator) :: old_opstring(size(opstring))
type(t_Stack) :: stack
integer :: smallest_unvisited_leg
logical :: LEGS_TO_BE_PROCESSED
! Has the spin at site ir been touched by a cluster ? 
! If not, flip it with probability 1/2.
logical :: touched( config%n_sites ) 

#ifdef DEBUG_CLUSTER_UPDATE
! Has the operator at propagation step ip been visited by the current Swendsen-Wang cluster ?
! (only for visualization purposes)
logical :: visited_ip( config%LL )
! Space-time lattice of spins 
integer :: spins_spacetime(config%n_sites, config%LL)
#endif 

! FLIPPING indicates whether a Swendsen-Wang cluster should be flipped 
! or not (with probability 1/2). 
logical :: FLIPPING
logical :: WINDING_MACROSPIN( config%n_sites )  

real(dp), parameter :: ONE_HALF = 0.5
real(dp) :: prob

integer :: leg_start, leg, leg_next
integer :: dir 
integer :: ip, ir, l  

#ifdef DEBUG_CLUSTER_UPDATE
integer :: spins2(size(spins,dim=1))
integer :: i1, i2, i3
    print*, "initial spin config=", spins(:)

    spins2(:) = spins(:)
    do ip=1,config%LL
        i1 = opstring(ip)%i
        i2 = opstring(ip)%j
        i3 = opstring(ip)%k
        spins_spacetime(:, ip) = spins2(:)
        ! Propagate spins as spin-flip operators are encountered.
        if ((i1.ne.0).and.(i2.eq.0)) then
            spins2(i1) = -spins2(i1)
        endif
    enddo

#endif 

call stack%init( config%n_legs )  
touched(:) = .FALSE.
WINDING_MACROSPIN(:) = .FALSE.
! first cluster is always built from the smallest unvisited leg
! (which need not be leg 1 if identities also carry ghostlegs)
smallest_unvisited_leg = 1 
LEGS_TO_BE_PROCESSED = .TRUE.

do while(leg_visited(smallest_unvisited_leg))
    smallest_unvisited_leg = smallest_unvisited_leg + 1
    if( smallest_unvisited_leg == config%n_ghostlegs) then 
        LEGS_TO_BE_PROCESSED = .FALSE.
        exit
    endif 
enddo 


deterministic_cluster_construction: do while( LEGS_TO_BE_PROCESSED )

    ! ! ! ! Swendsen-Wang: flip cluster or not 
    ! ! ! call random_number(prob)
    ! ! ! if (prob < ONE_HALF) then 
    ! ! !     FLIPPING = .TRUE.
    ! ! ! else
    ! ! !     FLIPPING = .FALSE.
    ! ! ! endif 

    ! Weight ratio for flipping the current cluster
    ! due to the exchange field as a result of a 
    ! longitudinal field.
    wratio_exchange_field = 1.0_dp
    ! Copy the operator string in order to restore it 
    ! if a cluster flip is rejected  
    old_opstring(:) = opstring(:)
    WINDING_MACROSPIN(:) = .false.
    FLIPPING = .true.
        
    leg_start = smallest_unvisited_leg    
    call stack%push( leg_start )
    leg_visited( leg_start ) = .TRUE.        
#ifdef DEBUG_CLUSTER_UPDATE
    print*, "============================="    
    print*, "process start_leg", leg_start
    print*, "FLIPPING=", FLIPPING
    visited_ip(:) = .FALSE.
    ip = (leg_start-1) / MAX_GHOSTLEGS + 1 
    visited_ip(ip) = .TRUE.
#endif    
    call process_leg( leg_start, FLIPPING, opstring, &
                      stack, leg_visited, touched, &
                      hz_fields, C_par_hyperparam )
    
    do while( .not.stack%is_empty() )
        leg = stack%pop()   
        leg_next = vertexlink(leg)

        ! Check for winding macrospin 
        dir = leg_direction(leg)     
        if( (leg_next - leg)*dir < 0 ) then
            ip = (leg_next-1) / MAX_GHOSTLEGS + 1 
            ir = gleg_to_ir( opstring(ip), leg_next )
            if( FLIPPING ) then
                ! Note: The way the code is written a site ir can be checked 
                ! twice for a winding macrospin. Nonetheless, each site is only 
                ! visited at most once by any cluster branch. 
                WINDING_MACROSPIN(ir) = .TRUE.
#ifdef DEBUG_CLUSTER_UPDATE 
                print*, "winding macrospin at ir=", ir
#endif
            endif 
        endif 
        
        if( .not.leg_visited(leg_next) ) then 
#ifdef DEBUG_CLUSTER_UPDATE            
            print*, "process next_leg", leg_next
            ip = (leg_next-1) / MAX_GHOSTLEGS + 1 
            visited_ip(ip) = .TRUE.
#endif             
            leg_visited(leg_next) = .TRUE.                        
            call process_leg( gleg=leg_next, FLIPPING=FLIPPING, opstring=opstring, &
                stack=stack, leg_visited=leg_visited, touched=touched, &
                hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam )                        
        endif 
    enddo

    ! Decide whether to flip the just constructed cluster
    ! with heat bath probability
    call random_number(prob)
    if( prob < wratio_exchange_field / (1.0_dp + wratio_exchange_field) ) then 
    ! if( prob < wratio_exchange_field ) then 
        ! Accept the flipped clusters and 
        ! update the initial spin configuration.
#ifdef DEBUG_CLUSTER_UPDATE
        print*, "flip the winding macrospins, weight_ratio=", wratio_exchange_field
#endif 
        do ir = 1, config%n_sites
            if (WINDING_MACROSPIN(ir)) then
                spins(ir) = -spins(ir)
            endif
        enddo
    else
        ! Reject the flipping of the cluster
        opstring(:) = old_opstring(:)
#ifdef DEBUG_CLUSTER_UPDATE        
        print*, "reject cluster flip, weight_ratio=", wratio_exchange_field
#endif         
    endif 

    ! Which legs have not been processed yet ?
    do l = smallest_unvisited_leg, config%n_ghostlegs + 1
        if( l == config%n_ghostlegs + 1) then 
            LEGS_TO_BE_PROCESSED = .FALSE.
            exit
        endif 
        ! Note: ghostlegs should have been marked as visited 
        !       during construction of the linked list.
        if( .NOT.leg_visited(l) ) then
           smallest_unvisited_leg = l
           exit
        endif 
    enddo
    
#ifdef DEBUG_CLUSTER_UPDATE
    ! Output the SSE configuration after each Swendsed-Wang cluster construction
    call output_SSE_config(config, opstring, spins, visited_ip, filename="SSEconfig.dat")
#endif 

enddo deterministic_cluster_construction

        ! Flip all spins that are not part of a cluster with probability 1/2.
do ir = 1, config%n_sites
    if (.not.touched(ir)) then
      call random_number(prob)
      if( prob.le. ONE_HALF ) then
        spins(ir) = -spins(ir)
      endif    
    endif
enddo    

#ifdef DEBUG_CLUSTER_UPDATE
! Output the final SSE configuration after all custers have been built
! and after the initial spin configuration has been updated.
    ! call  output_SSE_config(config, opstring, spins, visited_ip, filename="SSEconfig.dat")
#endif 

#ifdef DEBUG_CLUSTER_UPDATE    
        ! Check after all cluster have been built that 
        ! plaquette operators still sit on minimally frustrated plaquettes 
        spins2(:) = spins(:)
        do ip=1, config%LL
            i1 = opstring(ip)%i
            i2 = opstring(ip)%j
            i3 = opstring(ip)%k
            if (i1 < 0) then 
                if(abs(spins2(abs(i1)) + spins2(abs(i2)) + spins2(abs(i3))) /= 1) then 
                    print*, "ERROR: plaquette operator on maximally frustrated spin config."
                    print*, spins2(abs(i1)), spins2(abs(i2)), spins2(abs(i3))
                    print*, "i1=", abs(i1), "i2=", abs(i2), "i3=", abs(i3)
                    print*, "This should not happen. Exiting ..."
                    stop
                endif 
            endif        
            ! Propagate spins as spin-flip operators are encountered.
            if ((i1 > 0).and.(i2 == 0)) then
                spins2(i1) = -spins2(i1)
            endif
        enddo
#endif

#ifdef DEBUG_CLUSTER_UPDATE
    print*, "final spin config=", spins(:)
#endif 

end subroutine


subroutine process_leg( &
    gleg, FLIPPING, opstring, stack, leg_visited, touched, &
    hz_fields, C_par_hyperparam )
! Purpose:
! --------  
! Process a leg that has been popped from the stack:
! - Find the other legs connected to the same vertex
!   and put them on the stack if the vertex is an Ising vertex. 
! - Implement the cluster construction rules, i.e. the 
!   cluster stops if a spin-flip (constant) vertex is encountered
!   and the spin-flip (constant) vertex is transformed into
!   constant (spin-flip) vertex. If a longitudinal vertex is encountered,
!   record the exchange field and accumulated the weight ratio for 
!   the flipped vertices. 
  
use class_Stack
implicit none

! Arguments:
! ----------
integer, intent(in)                 :: gleg 
logical, intent(in)                 :: FLIPPING ! whether the Swendsen-Wang cluster is to be flipped or not 
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Stack), intent(inout)        :: stack
logical, intent(out)                :: leg_visited(:)
logical, intent(out)                :: touched(:)
real(dp), intent(in)                :: hz_fields(:)      ! site-dependent longitudinal fields
real(dp), intent(in)                :: C_par_hyperparam  ! hyperparameter for algorithm with longitudinal field


! local variables 
integer :: ip, i1, i2, i3
integer :: dir 
integer :: vleg, vleg_mod
integer :: leg1, leg2, leg3, leg4, leg5 
integer :: alignment 

! Find the operator to which the leg is connected
! and which leg it is in a numbering scheme 
! relative to the given vertex. 
ip = (gleg-1) / MAX_GHOSTLEGS + 1
vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1

! i1, i2 and i3 are just needed to mark the sites 
! as "touched" by a cluster. Spins at the initial propagation step
! that are not touched by any operator, form a cluster by themselves
! and can be flipped with probability 1/2. 
i1 = opstring(ip)%i
dir = leg_direction(gleg)

operator_types: select case( opstring(ip)%optype )
    case( TRIANGULAR_PLAQUETTE )
#ifdef DEBUG_CLUSTER_UPDATE
        print*, "leg connected to triangular plaquette"
#endif 
        ! Triangular plaquette operator encountered.
        i2 = opstring(ip)%j
        i3 = opstring(ip)%k
        ! For triangular plaquettes the cluster construction rules 
        ! are those described in Ref. [1] 
        touched((/abs(i1), abs(i2), abs(i3)/)) = .TRUE.

        vleg_mod = mod(vleg-1, MAX_GHOSTLEGS_HALF) + 1
        ! Whether an update is an A-update, B-update, or C-update
        ! is determined (in the current implementation) during the 
        ! diagonal update where the components of opstring(:) 
        ! are initialized depending on the chosen update type. 
        if (opstring(ip)%PRIVILEGED_LEG_IS_MAJORITY_LEG) then  
            ! A-site takes part in majority spin configuration  
            IF( vleg_mod == A_LEG ) THEN  
                ! entrance leg is privileged leg
                ! put only the other privileged leg onto the stack 
                ! (Which vleg number this is depends 
                ! on whether the entrance leg points down or up.)
                leg1 = gleg - dir*3
                call stack%push(leg1)
                leg_visited(leg1) = .TRUE.

            ELSE ! entrance leg is ordinary leg 
                ! put all ordinary legs except for the entrance leg onto the stack 
                if (dir.eq.DOWN) then 
                    IF( vleg_mod == B_LEG ) THEN
                        ! put legs 
                        !      gleg + 1, gleg + 3, gleg + 4   onto stack 
                        leg1 = gleg + 1
                        leg2 = gleg + 3
                        leg3 = gleg + 4		  		
                    ELSEIF( vleg_mod == C_LEG ) THEN
                        ! put legs 
                        !      gleg - 1, gleg + 2, gleg + 3   onto stack 
                        leg1 = gleg - 1
                        leg2 = gleg + 2
                        leg3 = gleg + 3
                    ENDIF 
                else ! dir.eq.UP
                    IF( vleg_mod == B_LEG ) THEN
                        ! put legs 
                        !      gleg - 3, gleg - 2, gleg + 1   onto stack 
                        leg1 = gleg - 3
                        leg2 = gleg - 2
                        leg3 = gleg + 1		
                    ELSEIF( vleg_mod == C_LEG ) THEN
                        ! put legs 
                        !      gleg - 1, gleg - 3, gleg - 4   onto stack
                        leg1 = gleg - 1		  
                        leg2 = gleg - 3
                        leg3 = gleg - 4      
                    ENDIF 
                endif ! dir = [UP | DOWN]
                call stack%push_many( (/ leg1, leg2, leg3 /) )
                leg_visited(leg1) = .TRUE.          
                leg_visited(leg2) = .TRUE.          
                leg_visited(leg3) = .TRUE.     
            ENDIF ! ENDIF: entrance leg is ordinary leg   

        else ! A-site on minority spin configuration
            ! => Put all legs except for the entrance leg onto the stack. 
            if (dir.eq.DOWN) then
                IF( vleg_mod == A_LEG ) THEN 
                    ! put gleg +1, +2, +3, +4, +5 onto stack          
                    leg1 = gleg + 1        
                    leg2 = gleg + 2
                    leg3 = gleg + 3
                    leg4 = gleg + 4
                    leg5 = gleg + 5
                ELSEIF( vleg_mod == B_LEG ) THEN
                    ! put gleg -1, +1, +2, +3, +4 onto stack          
                    leg1 = gleg - 1
                    leg2 = gleg + 1
                    leg3 = gleg + 2
                    leg4 = gleg + 3
                    leg5 = gleg + 4
                ELSEIF( vleg_mod == C_LEG ) THEN
                    ! put gleg -2, -1, +1, +2, +3 onto stack 
                    leg1 = gleg - 2
                    leg2 = gleg - 1
                    leg3 = gleg + 1
                    leg4 = gleg + 2
                    leg5 = gleg + 3      
                ENDIF
            else ! dir == UP
                IF( vleg_mod == A_LEG ) THEN 
                    ! put gleg -3, -2, -1, +1 ,+2 onto stack       
                    leg1 = gleg - 3	  
                    leg2 = gleg - 2
                    leg3 = gleg - 1         
                    leg4 = gleg + 1
                    leg5 = gleg + 2      
                ELSEIF( vleg_mod == B_LEG ) THEN 
                    ! put gleg -4, -3, -2, -1, +1 onto stack       
                    leg1 = gleg - 4	  
                    leg2 = gleg - 3
                    leg3 = gleg - 2
                    leg4 = gleg - 1
                    leg5 = gleg + 1      
                ELSEIF( vleg_mod == C_LEG ) THEN 
                    ! put gleg -5, -4, -3, -2, -1 onto stack     
                    leg1 = gleg - 5
                    leg2 = gleg - 4
                    leg3 = gleg - 3
                    leg4 = gleg - 2
                    leg5 = gleg - 1
                ENDIF   
            endif ! dir = [UP | DOWN]
            call stack%push_many( (/ leg1, leg2, leg3, leg4, leg5 /) )
            leg_visited(leg1) = .TRUE.
            leg_visited(leg2) = .TRUE.
            leg_visited(leg3) = .TRUE.
            leg_visited(leg4) = .TRUE.
            leg_visited(leg5) = .TRUE.   

#ifdef DEBUG_CLUSTER_UPDATE
            print*, "putting on the stack legs ", leg1, leg2, leg3, leg4, leg5
#endif

        endif ! A-site on majority spin configuration ?

    case( ISING_BOND )
    ! Ising operator encountered
        i2 = opstring(ip)%j
    ! Extract the three other legs, put them on the stack and mark them.
        if( vleg == 1 ) then 
            leg1 = gleg + 1
            leg2 = gleg + MAX_GHOSTLEGS_HALF
            leg3 = gleg + MAX_GHOSTLEGS_HALF + 1       
        elseif( vleg == 2 ) then 
            leg1 = gleg - 1
            leg2 = gleg + MAX_GHOSTLEGS_HALF - 1
            leg3 = gleg + MAX_GHOSTLEGS_HALF
        elseif( vleg == MAX_GHOSTLEGS_HALF + 1 ) then 
            leg1 = gleg - MAX_GHOSTLEGS_HALF
            leg2 = gleg - MAX_GHOSTLEGS_HALF + 1
            leg3 = gleg + 1        
        elseif( vleg == MAX_GHOSTLEGS_HALF + 2) then 
            leg1 = gleg - MAX_GHOSTLEGS_HALF - 1
            leg2 = gleg - MAX_GHOSTLEGS_HALF
            leg3 = gleg - 1
        endif 
        call stack%push_many( (/ leg1, leg2, leg3 /) )
        leg_visited(leg1) = .TRUE.
        leg_visited(leg2) = .TRUE.       
        leg_visited(leg3) = .TRUE.    
        touched((/i1, i2/)) = .TRUE.

#ifdef DEBUG_CLUSTER_UPDATE
        print*, "putting on the stack legs ", leg1, leg2, leg3
#endif

    case( SPIN_FLIP )
    ! The leg is connected to a spin flip operator.
    ! Update: spin flip -> const. The cluster branch STOPS,
    ! i.e. no new legs are put on the stack.
#ifdef DEBUG_CLUSTER_UPDATE
        print*, "hitting on a 2-leg vertex -> STOP"
#endif
        if (FLIPPING) then
            opstring(ip)%i = i1
            opstring(ip)%j = i1
            opstring(ip)%optype = CONSTANT 
        endif    
        touched(i1) = .TRUE.

    case( CONSTANT )
    ! The leg is connected to a constant operator.
    ! Update: const -> spin flip. The cluster branch STOPS.
#ifdef DEBUG_CLUSTER_UPDATE
        print*, "hitting on a 2-leg vertex -> STOP"
#endif
        if( FLIPPING ) then
            opstring(ip)%i = i1
            opstring(ip)%j = 0
            opstring(ip)%optype = SPIN_FLIP
        endif
        touched(i1) = .TRUE.

    case( LONGITUDINAL )
        ! The leg is connected to a longitudinal field operator. 
        ! The cluster does NOT stop, i.e. put the other leg onto
        ! the stack. Accumulate the product of the weight ratios
        ! due to the exchange field that arises if the cluster is flipped.
        ! dir = +1 => UP; dir = -1 => DOWN 
        leg2 = gleg - dir * MAX_GHOSTLEGS_HALF

        call stack%push(leg2)
        leg_visited(leg2) = .true.
        touched(i1) = .true.    

        if( FLIPPING ) then 
            ! For longitudinal field (`hz`) operators, the structure component
            ! opstring(ip)%k is used to store whether the hz operator 
            ! is aligned with the spin it is sitting on or not. opstring(ip)%k = +1(-1) means
            ! aligned (anti-aligned).
            alignment = opstring(ip)%k
            if( alignment > 0 ) then 
                ! spin aligned with the field:
                ! accumulate weight ratio weight_new / weight_old 
                wratio_exchange_field = wratio_exchange_field &
                    * C_par_hyperparam / (TWO * abs(hz_fields(i1)) + C_par_hyperparam)
            else 
                ! spin not aligned or zero longitudinal field 
                if( C_par_hyperparam > 0 ) then 
                    wratio_exchange_field = wratio_exchange_field &
                        * (TWO * abs(hz_fields(i1)) + C_par_hyperparam) / C_par_hyperparam
                else
                    wratio_exchange_field = 1.0_dp
                endif 
            endif 
            ! Convert (hz,aligned) into (hz,antialigned) vertex and vice versa.
            ! Note: (hz,aligned) and (hz,antialigned) are two different vertices with 
            !       different weights. 
            opstring(ip)%k = -opstring(ip)%k
        endif    
        
    case default
        print*, "cluster update: strange operator detected"
        print*, "ip=", ip, "i1=", opstring(ip)%j, "i2=",i2, "i3=", opstring(ip)%k
        print*, "optype=", opstring(ip)%optype
        stop

    end select operator_types

end subroutine process_leg

end module 
