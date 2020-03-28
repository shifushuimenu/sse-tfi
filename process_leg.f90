subroutine process_leg( gleg, opstring, stack, leg_visited, touched )
implicit none

integer, intent(in)                :: gleg 
type(tBondOperator), intent(inout) :: opstring(:)
type(t_Stack), intent(out)         :: stack
integer, intent(out)               :: leg_visited(:)
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

! IMPROVE: change this function 
dir = leg_direction(i1, i2, gleg)

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