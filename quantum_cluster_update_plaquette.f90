! TODO:
! - 'gleg' is the "global" leg number 
!   whereas 'leg...' denotes the leg number on a given vertex. 
! - MAX_GHOSTLEGS = 6 "number of legs around a vertex, including ghostlegs"


module class_Stack
! Pseudo object oriented design 
    implicit none 
    private 
    public t_Stack, push, push_many, pop, init
    
    type t_Stack
        integer, allocatable :: vals(:)
        integer :: stack_position 
    end type 

    contains 
    
    function pop(this) result(v)
        type(t_Stack), intent(inout) :: this
        integer :: v
        v = this%vals(this%stack_position)        
        this%vals(this%stack_position) = 0
        this%stack_position = this%stack_position - 1         
    end function 
    
    subroutine push(this, v)
        type(t_Stack), intent(inout) :: this
        integer, intent(in) :: v
        this%stack_position = this%stack_position + 1 
        this%vals(this%stack_position) = v 
    end subroutine 
    
    subroutine push_many(this, array_of_vals)
        type(t_Stack), intent(inout) :: this
        integer, intent(in) :: array_of_vals(:)
        integer :: i
        do i=1,size(array_of_vals, 1)
            this%stack_position = this%stack_position + 1 
            this%vals(this%stack_position) = array_of_vals(i)            
        enddo        
    end subroutine 
    
    subroutine init(this, n)
        type(t_Stack), intent(inout) :: this
        integer, intent(in) :: n                
        allocate(this%vals(n))
        this%vals(:) = 0
        this%stack_position = 0        
    end subroutine 
    
! As soon as the suroutine in which the stack is allocated exits,
! the memory is freed. 
! !     subroutine delete(this)
! !         type(t_Stack), intent(inout) :: this 
! !         deallocate(this%vals)        
! !     end subroutine 
    
end module class_Stack


module cluster_update

    implicit none 

    ! convenient labels 
    integer, parameter :: UP=+1, DOWN=-1
    integer, parameter :: MAX_GHOSTLEGS = 6
    ! operator types 
    integer, parameter :: IDENTITY = 100
    integer, parameter :: ISING_BOND = 10
    integer, parameter :: TWO_LEG = 20, SPIN_FLIP = 21, CONSTANT = 22
    integer, parameter :: TRIANGULAR_PLAQUETTE = 30
    
    ! for triangular-plaquette based update 
    integer, parameter :: A_LEG=1, B_LEG=2, C_LEG=3
    
    contains 
    
    use class_Stack
    implicit none 
    
    contains         
    
    pure function gleg_to_ir(op, gleg) result(ir)
    ! ***************************************************************
    ! Purpose: 
    ! --------
    ! Iplements the mapping from the global leg number 'gleg'
    ! of the operator 'op' to the linealry stored lattice site 'ir'
    ! on which it sits. 
    ! 
    ! This function is only called if a winding macrospin 
    ! is detected. 
    ! ***************************************************************
    
        implicit none 
        type(tBondOperator) :: op   ! operator at propagation step ip 
        integer, intent(in) :: gleg ! global leg number 
        
        integer :: vleg ! leg number around a vertex 
        integer :: ir   ! linearly stored index of a lattice site 
         
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
    
    
    pure function leg_direction(site_i, site_j, vleg) result(dir)
        ! opstring(ip)%i == site_i
        ! opstring(ip)%j == site_j
        implicit none 
        ! vleg is the leg number around a vertex
        integer, intent(in) :: vleg 
        integer, intent(in) :: site_i, site_j
        integer :: dir
    
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
        type(tBondOperator), intent(in) :: op
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
        
        
!-----------------------------------------------!
subroutine quantum_cluster_update_plaquette( &
    spins, opstring, vertexlink, config )
!-----------------------------------------------!

! Swendsen-Wang variant of the cluster update
! ============================================

  use configuration
  use legposition
  use various ! needed for variable DEBUG_CLUSTER_UPDATE
  implicit none
  
  integer, intent(inout) :: spins(:)
  type(tBondOperator), intent(inout) :: opstring(:)
  type(tConfig), intent(in) :: config  
   
  REAL(dp), PARAMETER :: ONE_HALF = 0.5
   
  INTEGER :: stack( config%n_legs )
  INTEGER :: stack_pos 	! index of the last element on the stack
  
  ! The cluster to which the leg belongs has been visited and its legs
  ! have been put on the stack   	     
  LOGICAL :: leg_visited( config%n_legs ) 
  INTEGER :: smallest_unvisited_leg
  
  ! Has the spin at site ir been touched by a cluster ? 
  ! If not, flip it with probability 1/2.
  LOGICAL :: touched( config%n_sites ) 

  ! FLIPPING indicates whether a Swendsen-Wang cluster should be flipped 
  ! or not (with probability 1/2). 
  LOGICAL :: FLIPPING
  LOGICAL :: WINDING_MACROSPIN( config%n_sites )  
  
  integer :: l, i1, i2, ip, ir, leg, leg_next, leg1, leg2, leg3, leg4, leg5, &
	     ir_start, leg_start, other_site, current_site
    
  ! variables for Swendsen-Wang algorithm
  integer :: n_clusters
  integer :: dir
  
  double precision :: eta
  
  double precision :: ran2 !!!!! Is this necessary ????
  
  ! for plaquette-based quantum cluster update
  integer :: ir_A, ir_B, ir_C

! ******************************************************************* 
! ALGORITHM: 
! *******************************************************************
! Build all Swendsen-Wang clusters on the space-time lattice formed 
! by the linked list of vertex-legs. 
! 
! 1. Choose with probability 1/2 whether the cluster to be 
!    constructed should be flipped (FLIPPING==.TRUE.) or not.
!
! 2. Choose an unvisited starting leg at random.
!    Put it on the stack. Set the variable start_macrospin=[leg] and dir=[up|down].
!    The linked list will give the next entrance leg.
!    If the vertex to which the entrance leg belongs is 
!       -- an Ising operator
!       -- or a plaquette operator,
!    then put all remaining legs of the operator on the global stack and carry
!    on with the other leg of the current site.
!    If the new vertex is 
!       -- a constant 
!       -- or a spin flip operator, 
!    then stop and exchange the operator (const <-> spin flip). 
!    Set stop_macrospin=[leg].
! Invert the direction and carry on with the last leg on the onsite stack.
!
! Once a cluster has been completely traversed, start another cluster from 
! a leg which has not yet been visited, until all clusters are constructed.
!
! Note: Single constants or spin-flip operators acting on a site
! are visited twice by the growing cluster and are thus left unaffected.
!
! Criterion for a winding macrospin:
!   The leg number decreases when going up or
!   the leg number increases when going down.
! 
! **********
! GHOSTLEGS
! **********
! Thanks to 'ghostlegs' the mapping from the 'global' leg number 
! 'gleg' to the coordinates of a leg 
! 
!      gleg  <--->  (ip,ir)  
! 
! can be computed on the fly as 
! 
!      ip = gleg // MAX_GHOSTLEGS + 1          gleg = (ip-1)*MAX_GHOSTLEGS + vleg,
!
! where 'vleg' is the local leg number around a given vertex, and 
! 
!      opstring(ip) ---> (i, j, k)
!
!      vleg = mod(gleg-1, MAX_GHOSTLEGS) + 1
!
! Note: It is implicitly assumed that 'ghostlegs' 
!       can never appear as global leg numbers.    
!
! Depending on the operator type sitting at ip, the local leg 
! can be associated to a lattice site (one of i,j,or k).  
! 
! *******************************************************************

! Initialization 
leg_visited(:) = .FALSE.
! Note: Each spin belongs only to one cluster, 
! i.e. a spin cannot be flipped twice.
touched(:) = .FALSE.
! The first leg of the first cluster is chosen to be leg=1.
smallest_unvisited_leg = 1

#ifdef DEBUG_CLUSTER_UPDATE
! Do a cluster update if there is at least 
! one vertex in the operator list.
if (n_legs < 2) then
    STOP "cluster_update: ERROR: operator list is empty"
endif 
#endif 

#ifdef CLUSTER_STATS
 n_clusters = 0
#endif 

! *****************************************************************
! Build all clusters according to Swendsen-Wang
! *****************************************************************
! Note: The building of clusters is completely deterministic:
! The shape of the clusters is determined by the operator list.
! Each vertex leg belongs to exactly one cluster. 
! *****************************************************************

! Initialize WINDING_MACROSPIN(:) before constructing any cluster 
! Note: each site should be visited at most by one cluster 
!       branch. 
WINDING_MACROSPIN(:) = .FALSE.

do while( smallest_unvisited_leg <= n_legs ) 

#ifdef CLUSTER_STATS
 n_clusters = n_clusters + 1
#endif  

! Decide with probability 1/2 whether this cluster is to 
! be flipped or not.
eta = ran2(idum)
IF( eta <= ONE_HALF ) THEN
    ! flip the cluster
    FLIPPING = .TRUE.
#ifdef DEBUG_CLUSTER_UPDATE
      print*, "========================="
      print*, "flipping cluster"
      print*, "leg visited?",(leg_visited(l), l=1,n_legs)
#endif 
ELSE
    ! do not flip the cluster
    FLIPPING = .FALSE.
#ifdef DEBUG_CLUSTER_UPDATE
      print*, "========================="
      print*, "non-flipping cluster"
      print*, "leg visited?", (leg_visited(l), l=1,n_legs)
#endif 
ENDIF
	    	    
! After building a cluster completely the per-cluster stack should be empty.
! Initialize anyway, or check that the stack is really empty.
#ifdef DEBUG_CLUSTER_UPDATE
    if( (stack_pos /= 0).OR.(.NOT.ALL(stack == 0)) ) then
        STOP "Stack is not empty at the beginning of a new cluster."
    endif 
#endif 
stack_pos = 0
stack(:) = 0

! Choose as a STARTING LEG the smallest unvisited leg.
! The processing of the starting leg is slightly different (?????), therefore
! it is done separately from the legs which appear when the cluster 
! branches out.
! The difference is that, for the starting leg, we look to which operator 
! it is connected and put the other legs of this operator on the stack
! for processing while, for subsequent legs, which are popped from the stack,
! we follow the linked list to the next entrance leg.

leg_start = smallest_unvisited_leg

! add the first leg to the stack and mark it
stack_pos = stack_pos + 1
stack(stack_pos) = leg_start
leg_visited(leg_start) = .TRUE.

! Find the operator to which the leg is connected
! and which leg it is in a numbering scheme 
! relative to the given vertex. 
ip = leg_start / MAX_GHOSTLEGS + 1
vleg = mod(leg_start-1, MAX_GHOSTLEGS) + 1

! i1 and i2 are just needed to determine the operator type
i1 = opstring(ip)%i
i2 = opstring(ip)%j

! IMPROVE: change this function 
dir = leg_direction(i1, i2, vleg)

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
      leg1 = leg_start - dir*3
      stack_pos = stack_pos + 1
      stack(stack_pos) = leg1
      leg_visited(leg1) = .TRUE.
	              
	    else ! entrance leg is ordinary leg 
	      ! put all ordinary legs except for the entrance leg onto the stack 
	      if (dir.eq.DOWN) then 
	        IF( vleg_mod == B_LEG ) THEN
		! put legs 
		!      leg_start + 1, leg_start + 3, leg_start + 4   onto stack 
		  leg1 = leg_start + 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start + 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 4
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		
		 ELSEIF( vleg_mod == C_LEG ) THEN
		! put legs 
		!      leg_start - 1, leg_start + 2, leg_start + 3   onto stack 
		  leg1 = leg_start - 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start + 2
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		 
		 ENDIF 
	      else ! dir.eq.UP
	        IF( vleg_mod == B_LEG ) THEN
		! put legs 
		!      leg_start - 3, leg_start - 2, leg_start + 1   onto stack 
		  leg1 = leg_start - 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start - 2
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		
		ELSEIF( vleg_mod == C_LEG ) THEN
		! put legs 
		!      leg_start - 1, leg_start - 3, leg_start - 4   onto stack 
		  leg1 = leg_start - 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start - 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start - 4
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
				
		ENDIF 
	      
	      endif 
	      
	    endif 
        
 else ! A-site on minority spin configuration
   ! => Put all legs except for the entrance leg onto the stack. 
   if (dir.eq.DOWN) then
      IF( vleg_mod == A_LEG ) THEN 
         ! put leg_start +1, +2, +3, +4, +5 onto stack 
         
          leg1 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 5
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  

      ELSEIF( vleg_mod == B_LEG ) THEN
         ! put leg_start -1, +1, +2, +3, +4 onto stack 
         
          leg1 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ELSEIF( vleg_mode == C_LEG ) THEN
         ! put leg_start -2, -1, +1, +2, +3 onto stack 
                           
          leg1 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ENDIF
   else ! dir == UP
      IF( vleg_mod == A_LEG ) THEN 
      ! put leg_start -3, -2, -1, +1 ,+2 onto stack 
      
          leg1 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ELSEIF( vleg_mod == B_LEG ) THEN 
      ! put leg_start -4, -3, -2, -1, +1 onto stack 
      
	  leg1 = leg_start - 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  

      
      ELSEIF( vleg_mode == C_LEG ) THEN 
      ! put leg_start -5, -4, -3, -2, -1 onto stack 
      
          leg1 = leg_start - 5
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
            
      ENDIF

   
   endif ! dir = [UP | DOWN]
   
 endif ! A-site on majority spin configuration ?
 

elseif( (i1.gt.0).and.(i2.gt.0).and.(i1.ne.i2) ) then
! Ising operator encountered

   touched(i1) = .TRUE.
   touched(i2) = .TRUE.

   ! Extract the three other legs, put them on the stack and mark them.
   if( vleg == 1 ) then 
       leg1 = leg_start + 1
       leg2 = leg_start + 2
       leg3 = leg_start + 3
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
   elseif( vleg == 2 ) then 
       leg1 = leg_start - 1
       leg2 = leg_start + 1
       leg3 = leg_start + 2
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    elseif( vleg == 3 ) then 
       leg1 = leg_start - 2
       leg2 = leg_start - 1
       leg3 = leg_start + 1
        
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    elseif( vleg == 4) then 
       leg1 = leg_start - 1
       leg2 = leg_start - 2
       leg3 = leg_start - 3 
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    endif 
           
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

! *************************************************************************
! Let the cluster BRANCH OUT until all connected legs have been processed.
! While the stack is not empty ...
! *************************************************************************

!remove
#ifdef DEBUG_CLUSTER_UPDATE
  print*, "starting leg =", leg_start
#endif

do while ( stack_pos.ge.1 )

  ! Pop a leg from the stack and note the direction in which the cluster branches out
  leg = stack(stack_pos)
  stack(stack_pos) = 0
  stack_pos = stack_pos - 1
#ifdef DEBUG_CLUSTER_UPDATE
    print*, "popped from the stack: leg = ", leg
#endif

  ! branch out starting from this leg
  dir = leg_direction(leg)

  ! follow the linked list 
  leg_next = vertexlink(leg)
#ifdef DEBUG_CLUSTER_UPDATE
  print*, "following the linked list: leg_next = ", leg_next
#endif
  ! Check for winding macrospin 
  IF( (leg_next - leg)*dir < 0 ) THEN
        ir = gleg_to_ir(leg_next)
        IF( FLIPPING ) WINDING_MACROSPIN(ir) = .TRUE.
#ifdef DEBUG_CLUSTER_UPDATE
          print*, "winding macrospin at ir = ", ir
#endif
  ENDIF
	          
	          	          
  if (.not.leg_visited(leg_next)) then
        ! mark the next leg as belonging to the cluster
	leg_visited(leg_next) = .TRUE.
	
        ! Find the operator to which the leg is connected
        ! and which leg it is in a numbering scheme 
        ! relative to the given vertex. 
        ip = leg_start / MAX_GHOSTLEGS + 1
        vleg = mod(leg_start-1, MAX_GHOSTLEGS) + 1
        
	dir = leg_direction(leg_next)

	i1 = opstring(ip)%i
	i2 = opstring(ip)%j	

	if (i1.lt.0) then
	! Triangular plaquette operator encountered.
	! For plaquettes the cluster construction rules 
	! are those described in Ref. [1] 
	
	ir_A = abs(i1)
	ir_B = abs(i2)
	ir_C = abs(opstring(ip)%k)
	
	touched(ir_A) = .TRUE.
	touched(ir_B) = .TRUE.
	touched(ir_C) = .TRUE.
	
    
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
      leg1 = leg_start - dir*3
      stack_pos = stack_pos + 1
      stack(stack_pos) = leg1
      leg_visited(leg1) = .TRUE.
	              
	    else ! entrance leg is ordinary leg 
	      ! put all ordinary legs except for the entrance leg onto the stack 
	      if (dir.eq.DOWN) then 
	        IF( vleg_mod == B_LEG ) THEN
		! put legs 
		!      leg_start + 1, leg_start + 3, leg_start + 4   onto stack 
		  leg1 = leg_start + 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start + 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 4
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		
		 ELSEIF( vleg_mod == C_LEG ) THEN
		! put legs 
		!      leg_start - 1, leg_start + 2, leg_start + 3   onto stack 
		  leg1 = leg_start - 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start + 2
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		 
		 ENDIF 
	      else ! dir.eq.UP
	        IF( vleg_mod == B_LEG ) THEN
		! put legs 
		!      leg_start - 3, leg_start - 2, leg_start + 1   onto stack 
		  leg1 = leg_start - 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start - 2
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start + 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
		
		ELSEIF( vleg_mod == C_LEG ) THEN
		! put legs 
		!      leg_start - 1, leg_start - 3, leg_start - 4   onto stack 
		  leg1 = leg_start - 1
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg1
		  leg_visited(leg1) = .TRUE.          
		  
		  leg2 = leg_start - 3
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg2
		  leg_visited(leg2) = .TRUE.          

		  leg3 = leg_start - 4
		  stack_pos = stack_pos + 1
		  stack(stack_pos) = leg3
		  leg_visited(leg3) = .TRUE.          	  
				
		ENDIF 
	      
	      endif 
	      
	    endif 
        
 else ! A-site on minority spin configuration
   ! => Put all legs except for the entrance leg onto the stack. 
   if (dir.eq.DOWN) then
      IF( vleg_mod == A_LEG ) THEN 
         ! put leg_start +1, +2, +3, +4, +5 onto stack 
         
          leg1 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 5
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  

      ELSEIF( vleg_mod == B_LEG ) THEN
         ! put leg_start -1, +1, +2, +3, +4 onto stack 
         
          leg1 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ELSEIF( vleg_mode == C_LEG ) THEN
         ! put leg_start -2, -1, +1, +2, +3 onto stack 
                           
          leg1 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ENDIF
   else ! dir == UP
      IF( vleg_mod == A_LEG ) THEN 
      ! put leg_start -3, -2, -1, +1 ,+2 onto stack 
      
          leg1 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
      
      ELSEIF( vleg_mod == B_LEG ) THEN 
      ! put leg_start -4, -3, -2, -1, +1 onto stack 
      
	  leg1 = leg_start - 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start + 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  

      
      ELSEIF( vleg_mode == C_LEG ) THEN 
      ! put leg_start -5, -4, -3, -2, -1 onto stack 
      
          leg1 = leg_start - 5
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg1
	  leg_visited(leg1) = .TRUE.          
	  
	  leg2 = leg_start - 4
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg2
	  leg_visited(leg2) = .TRUE.          

	  leg3 = leg_start - 3
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg3
	  leg_visited(leg3) = .TRUE.          	  
         
	  leg4 = leg_start - 2
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg4
	  leg_visited(leg4) = .TRUE.          	  

	  leg5 = leg_start - 1
	  stack_pos = stack_pos + 1
	  stack(stack_pos) = leg5
	  leg_visited(leg5) = .TRUE.          	  
            
      ENDIF
   
   endif ! dir = [UP | DOWN]
   
 endif ! A-site on majority spin configuration ?


      	! The leg is connected to an Ising operator
	elseif ((i1.gt.0).and.(i2.gt.0).and.(i1.ne.i2)) then	
	  touched(i1) = .TRUE.
	  touched(i2) = .TRUE.
	
	  ! Extract the three other legs and put them on the stack. Don't mark them
   ! Extract the three other legs, put them on the stack and mark them.
   if( vleg == 1 ) then 
       leg1 = leg_start + 1
       leg2 = leg_start + 2
       leg3 = leg_start + 3
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
   elseif( vleg == 2 ) then 
       leg1 = leg_start - 1
       leg2 = leg_start + 1
       leg3 = leg_start + 2
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    elseif( vleg == 3 ) then 
       leg1 = leg_start - 2
       leg2 = leg_start - 1
       leg3 = leg_start + 1
        
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    elseif( vleg == 4) then 
       leg1 = leg_start - 1
       leg2 = leg_start - 2
       leg3 = leg_start - 3 
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg1
       leg_visited(leg1) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg2
       leg_visited(leg2) = .TRUE.
       
       stack_pos = stack_pos + 1
       stack(stack_pos) = leg3
       leg_visited(leg3) = .TRUE.
    endif 
    	  
#ifdef DEBUG_CLUSTER_UPDATE
	    print*, "putting on the stack legs ", leg1, leg2, leg3
#endif
	  
	! The leg is connected to a spin flip operator.
	elseif ((i1.ne.0).and.(i2.eq.0)) then	
	  touched(i1) = .TRUE.	  
#ifdef DEBUG_CLUSTER_UPDATE
	    print*, "branch STOPS"
#endif
          ! Update: spin flip -> const. The cluster branch STOPS,
          ! i.e. no new legs are put on the stack.
	  if (FLIPPING) then
	    opstring(ip)%i = i1
	    opstring(ip)%j = i1
	  endif
	  
	! The leg is connected to a constant operator.
	elseif ((i1.gt.0).and.(i2.eq.i1)) then
	  touched(i1) = .TRUE.
#ifdef DEBUG_CLUSTER_UPDATE 
	    print*, "branch STOPS"
#endif
          ! Update: const -> spin flip. The cluster branch STOPS.
	  if (FLIPPING) then
	    opstring(ip)%i = i1
	    opstring(ip)%j = 0
	  endif
	else
	    print*, "cluster update: strange operator detected"
	    stop
	endif
      
      else
#ifdef DEBUG_CLUSTER_UPDATE 
	print*, "leg already visited"
#endif
      endif ! if(.not.leg_visited(leg_next))
        
enddo ! while (stack_pos.ge.1)

! Which legs have not been processed yet ?
! The next cluster starts with the leg which has the highest index
! of all unvisited legs (convention).
legs_proc_pos = 0

DO l = smallest_unvisited_leg, n_legs
  IF( .not.leg_visited(l) ) THEN
     smallest_unvisited_leg = l
     BREAK
  ENDIF
ENDDO
smallest_unvisited_leg = n_legs + 1 ! invalid value so that the loop stops 

#ifdef DEBUG_CLUSTER_UPDATE
  print*, "winding ?", (WINDING_MACROSPIN(ir), ir=1,NN)
  print*, "spins", (spins(ir), ir=1,NN) 
#endif

enddo ! build all clusters

! Update the initial spin configuration 
! Note: Since each site at the initial propagation step can only be touched 
!    by one cluster, the update of the initial spin configuration can be postponed 
!    to the moment when all clusters have been built. 
do ir = 1, config%n_sites
    if (WINDING_MACROSPIN(ir)) then
      spins(ir) = -spins(ir)
    endif
enddo

! Flip all spins that are not part of a cluster with probability 1/2.
! (Very unlikely)
do ir = 1, config%n_sites
  if (.not.touched(ir)) then
    eta = ran2(idum)
    if( eta.le. ONE_HALF ) then
      spins(ir) = -spins(ir)
    endif    
  endif
enddo

end subroutine quantum_cluster_update_plaquette
!-----------------------------------------------!

end module cluster_update 