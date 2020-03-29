module cluster_update

use types, only: dp
use SSE_configuration 
implicit none 
private 

public quantum_cluster_update_plaquette, unit_test

! convenient labels 
integer, parameter :: UP=+1, DOWN=-1
! for triangular-plaquette based update 
integer, parameter :: A_LEG=1, B_LEG=2, C_LEG=3

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
    
    ip = gleg / MAX_GHOSTLEGS + 1
    vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1
    vleg_mod = mod(vleg-1, MAX_GHOSTLEGS/2) + 1
    
    select case( operator_type(op) )
        case(TRIANGULAR_PLAQUETTE)
            if( vleg_mod == A_LEG ) then
                ir = op%i
            elseif( vleg_mod == B_LEG ) then 
                ir = op%j
            elseif( vleg_mod == C_LEG ) then 
                ir = op%k
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
    end select 
        
end function

! IMPROVE: This function should be removed in favour 
! of a new labelling scheme of legs around a vertex
! where 
!    vleg > MAX_GHOSTLEGS/2 points UP
!    and vleg <= MAX_GHOSTLEGS/2 points DOWN
! even for two-leg and four-leg vertices. 
pure function leg_direction(opstring, gleg) result(dir)
    implicit none 
    type(t_BondOperator), intent(in) :: opstring(:)
    ! gleg is the global leg number
    integer, intent(in) :: gleg       
    
    integer :: dir
    integer :: ip, site_i, site_j
    ! vleg is the leg number around a vertex
    integer :: vleg 

    ip = gleg / MAX_GHOSTLEGS + 1
    ! Site %i and %j are just needed to determine the operator type.
    site_i = opstring(ip)%i 
    site_j = opstring(ip)%j

    vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1
    if( site_i < 0 ) then 
    ! triangular plaquette
        if(vleg <= 3) then
            dir = DOWN
        else
            dir = UP
        endif 
    elseif( site_i > 0 ) then 
        if( site_j > site_i) then 
        ! Ising operator
            if( vleg <= 2 ) then 
                dir = DOWN
            else
                dir = UP
            endif 
        elseif( site_j <= site_i) then
        ! constant or spin-flip operator 
            if ( vleg == 1) then 
                dir = DOWN
            else
                dir = UP
            endif 
        endif 
    endif 

end function leg_direction 
        
pure function operator_type(op) result(t)
    implicit none 
    type(t_BondOperator), intent(in) :: op
    integer :: t
                    
    if( op%i < 0 ) then 
#ifdef DEBUG_CLUSTER_UPDATE
        if( (op%j > 0).OR.(op%k > 0) ) then
            STOP "operator_type(): ERROR: strange operator detected"
        endif 
#endif         
        t = TRIANGULAR_PLAQUETTE
    elseif( op%i > 0) then 
        if( op%j > op%i ) then 
            t = ISING_BOND
        elseif( op%j == op%i ) then 
            t = CONSTANT 
        elseif( op%j == 0) then
            t = SPIN_FLIP
        else
            STOP "operator_type(): ERROR: strange operator detected"
        endif 
    elseif( op%i == 0) then 
        t = IDENTITY
    endif 
            
end function operator_type 
        

subroutine quantum_cluster_update_plaquette( &
    spins, opstring, vertexlink, config )
! ***************************************************
! Swendsen-Wang variant of the quantum cluster update
! ***************************************************
use class_Stack
implicit none

integer, intent(inout) :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
integer, intent(in) :: vertexlink(:)
type(t_Config), intent(in) :: config  
    
type(t_Stack) :: stack
logical :: leg_visited( config%n_ghostlegs ) 
integer :: smallest_unvisited_leg
logical :: LEGS_TO_BE_PROCESSED
! Has the spin at site ir been touched by a cluster ? 
! If not, flip it with probability 1/2.
logical :: touched( config%n_sites ) 
! FLIPPING indicates whether a Swendsen-Wang cluster should be flipped 
! or not (with probability 1/2). 
logical :: FLIPPING
logical :: WINDING_MACROSPIN( config%n_sites )  

real(dp), parameter :: ONE_HALF = 0.5
real(dp) :: prob

integer :: leg_start, leg, leg_next
integer :: dir 
integer :: ip, ir, l  

call stack%init( config%n_legs )  
leg_visited(:) = .FALSE.
touched(:) = .FALSE.
WINDING_MACROSPIN(:) = .FALSE.
smallest_unvisited_leg = 1
LEGS_TO_BE_PROCESSED = .TRUE.

do while( LEGS_TO_BE_PROCESSED )

    ! Swendsen-Wang: flip cluster or not 
    prob = ran()
    if (prob < ONE_HALF) then 
        FLIPPING = .TRUE.
    else
        FLIPPING = .FALSE.
    endif 
    
    leg_start = smallest_unvisited_leg    
    call stack%push( leg_start )    
    leg_visited( leg_start ) = .TRUE.        
    call process_leg( leg_start )
    
    do while( .not.stack%is_empty() )
        leg = stack%pop()   
        leg_next = vertexlink(leg)
        
        ! Check for winding macrospin 
        dir = leg_direction(opstring, leg)     
        if( (leg_next - leg)*dir < 0 ) then
            ip = leg_next / MAX_GHOSTLEGS + 1   
            ir = gleg_to_ir( opstring(ip), leg_next )
            if( FLIPPING ) WINDING_MACROSPIN(ir) = .TRUE.
        endif 
        
        if( .not.leg_visited(leg_next) ) then 
            
            leg_visited(leg_next) = .TRUE.
            
            call process_leg( leg_next )
            
            
        endif 
        
    enddo

    ! Which legs have not been processed yet ?
    ! Loop over all legs only once. 
    do l = smallest_unvisited_leg, config%n_ghostlegs + 1
        if( l == config%n_ghostlegs + 1) then 
            LEGS_TO_BE_PROCESSED = .FALSE.
            exit
        endif 
        if( .NOT.leg_visited(l) ) then
           smallest_unvisited_leg = l
           exit
        endif 
    enddo

enddo

! AT THE VERY END, WHEN ALL CLUSTERS HAVE BEEN BUILT...
! Update the initial spin configuration 
do ir = 1, config%n_sites
    if (WINDING_MACROSPIN(ir)) then
      spins(ir) = -spins(ir)
    endif
enddo
! Flip all spins that are not part of a cluster with probability 1/2.
do ir = 1, config%n_sites
  if (.not.touched(ir)) then
    prob = ran()
    if( prob.le. ONE_HALF ) then
      spins(ir) = -spins(ir)
    endif    
  endif
enddo

end subroutine


subroutine process_leg( &
    gleg, FLIPPING, opstring, stack, leg_visited, touched )
! *******************************************************************    
! Process a leg that has been popped from the stack. 
! *******************************************************************    
use class_Stack
implicit none

integer, intent(in)                :: gleg 
logical, intent(in)                :: FLIPPING ! whether the Swendsen-Wang cluster is to be flipped or not 
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Stack), intent(out)         :: stack
logical, intent(out)               :: leg_visited(:)
logical, intent(out)               :: touched(:)

! local variables 
integer :: ip, i1, i2
integer :: dir 
integer :: vleg, vleg_mod
integer :: leg1, leg2, leg3, leg4, leg5 

! Find the operator to which the leg is connected
! and which leg it is in a numbering scheme 
! relative to the given vertex. 
ip = gleg / MAX_GHOSTLEGS + 1
vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1

! i1 and i2 are just needed to determine the operator type
i1 = opstring(ip)%i
i2 = opstring(ip)%j

dir = leg_direction(opstring, gleg)

if( i1.lt.0 ) then
    ! Triangular plaquette operator encountered.
    ! For plaquettes the cluster construction rules 
    ! are those described in Ref. [1] 

    vleg_mod = mod(vleg-1, MAX_GHOSTLEGS/2) + 1
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
            ! on whether in entrance leg points down or up.)
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

    endif ! A-site on majority spin configuration ?

elseif( (i1.gt.0).and.(i2.gt.0).and.(i1.ne.i2) ) then
    ! Ising operator encountered

    touched(i1) = .TRUE.
    touched(i2) = .TRUE.

    ! Extract the three other legs, put them on the stack and mark them.
    if( vleg == 1 ) then 
        leg1 = gleg + 1
        leg2 = gleg + 2
        leg3 = gleg + 3       
    elseif( vleg == 2 ) then 
        leg1 = gleg - 1
        leg2 = gleg + 1
        leg3 = gleg + 2
    elseif( vleg == 3 ) then 
        leg1 = gleg - 2
        leg2 = gleg - 1
        leg3 = gleg + 1        
    elseif( vleg == 4) then 
        leg1 = gleg - 1
        leg2 = gleg - 2
        leg3 = gleg - 3 
    endif 
    call stack%push_many( (/ leg1, leg2, leg3 /) )
    leg_visited(leg1) = .TRUE.
    leg_visited(leg2) = .TRUE.       
    leg_visited(leg3) = .TRUE.

#ifdef DEBUG_CLUSTER_UPDATE
    print*, "putting on the stack legs ", leg1, leg2, leg3
#endif

elseif( (i1.ne.0).and.(i2.eq.0) ) then
    ! The leg is connected to a spin flip operator.
    ! Update: spin flip -> const. The cluster branch STOPS,
    ! i.e. no new legs are put on the stack.

    touched(i1) = .TRUE.

#ifdef DEBUG_CLUSTER_UPDATE
    print*, "starting on a 2-leg vertex"
#endif
    if (FLIPPING) then
        opstring(ip)%i = i1
        opstring(ip)%j = i1
    endif    

elseif( (i1.ne.0).and.(i2.eq.i1) ) then
    ! The leg is connected to a constant operator.
    ! Update: const -> spin flip. The cluster branch STOPS.

    touched(i1) = .TRUE.

#ifdef DEBUG_CLUSTER_UPDATE
    print*, "starting on a 2-leg vertex"
#endif
    if (FLIPPING) then
        opstring(ip)%i = i1
        opstring(ip)%j = 0
    endif
else
    print*, "cluster update: strange operator detected"
    print*, "ip=", ip, "i1=", i1, "i2=",i2
    stop
endif

end subroutine process_leg


subroutine unit_test

end subroutine     

end module 