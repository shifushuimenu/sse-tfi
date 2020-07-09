! TODO:
! - relabel legs around two-leg and four-leg vertices so 
!   that the direction into which a leg is pointing can be 
!   determined easily from the vertex leg number.
! - Ultimately, replace the array leg_visited( 1:MAX_GHOSTLEGS*n_exp )
!   by a different data structure which requires less memory. 
! - check self-consistency of the linked list 

! - Implement the linked list via pointers and derived data types 
!   so that the number of nodes is equal to the number of physical legs.
!   If the linked list is implemented via an array, then the number of 
!   entries of the array must be equal to the number of 'ghostlegs' because
!   the leg number is used as an index into the array, with the array value 
!   giving the connecting leg. 

MODULE linked_list

    USE SSE_configuration 
    
    CONTAINS

subroutine build_linkedlist_plaquette( &
  opstring, config, vertexlink, leg_visited )
!---------------------------------------------------------------!
! convention for numbering of vertex legs (`vleg`):             !
!                                                               !
!   Ising operator:		Spin-flip operator or constant:           !
!  (4)           (5)  [3,6]           (4)   [2,3,5,6]           !
!   |_____________|                    |                        !
!   |             |                    |                        !
!  (1)           (2)                  (1)                       !
! site-i        site-j	                                        !
!                                                               !
!                                                               !
!    Triangular plaquette operator:                             !
!       (4)           (5)           (6)                         !
!        |_____________|_____________|                          !
!        |A-site       |B-site       |C-site                    !
!       (1)           (2)           (3)                         !
!                                                               !
! Note 1: "ghostlegs" in square brackets are not connected to   !
!       any other legs.                                         !
! Note 2: The labelling scheme around a vertex is chosen such   ! 
!       that                                                    !
!          vleg > MAX_GHOSTLEGS/2 points UP                     !
!          and vleg <= MAX_GHOSTLEGS/2 points DOWN              !
!       even for two-leg and four-leg vertices.                 !
!---------------------------------------------------------------!
implicit none 

! automatic array 
type(t_BondOperator), intent(in) :: opstring(:)
type(t_Config), intent(in) :: config  
integer, allocatable, intent(out) :: vertexlink(:)
logical, allocatable, intent(out) :: leg_visited(:)

! automatic arrays used for building the linked list 
integer :: firstleg( config%n_sites )
integer :: lastleg( config%n_sites ) 

integer :: i, ip, i1, i2, leg_counter
integer :: ir_A, ir_B, ir_C

if( allocated(vertexlink) ) deallocate(vertexlink)
allocate(vertexlink(config%n_ghostlegs))
if( allocated(leg_visited) ) deallocate(leg_visited)
allocate(leg_visited(config%n_ghostlegs))

leg_visited(:) = .FALSE.
firstleg(:) = 0
lastleg(:) = 0 

! initialize linked list with invalid leg numbers 
vertexlink(:) = -1

! The first leg has index 1.
leg_counter = 0

! **********************************************
! Traverse operator list leaving out identities
! **********************************************
do ip=1, config%LL

! Determine the type of vertex.
  i1 = opstring(ip)%i  
  i2 = opstring(ip)%j

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Identity operator encountered
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if((i1 == 0) .and. (i2 == 0)) then 
  ! necessary for the mapping leg_number -> ip
  leg_visited(leg_counter+1:leg_counter+MAX_GHOSTLEGS) = .TRUE.
  leg_counter = leg_counter + MAX_GHOSTLEGS
endif 

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Six-leg vertex, i.e.
! triangular plaquette operator encountered
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (i1.lt.0) then 
! By construction the opstring(ip) structure components %i, %j, and %k 
! contain the linearly stored site index of the A-site, B-site, and C-site, 
! respectively. 
   ir_A = abs(i1)
   ir_B = abs(i2)
   ir_C = abs(opstring(ip)%k)

   ! lower three legs 
   if (lastleg(ir_A).ne.0) then
   ! double-linked list
       ! A-site always has lowest leg number on a triangular plaquette 
       vertexlink(lastleg(ir_A)) = leg_counter + 1 
       vertexlink(leg_counter + 1) = lastleg(ir_A)
   else
       firstleg(ir_A) = leg_counter + 1
   endif
      
   if (lastleg(ir_B).ne.0) then
       vertexlink(lastleg(ir_B)) = leg_counter + 2
       vertexlink(leg_counter + 2) = lastleg(ir_B)
   else
       firstleg(ir_B) = leg_counter + 2
   endif
   
   if (lastleg(ir_C).ne.0) then
       vertexlink(lastleg(ir_C)) = leg_counter + 3
       vertexlink(leg_counter + 3) = lastleg(ir_C)
   else
       firstleg(ir_C) = leg_counter + 3
   endif   
    
   ! upper three legs 
   lastleg(ir_A) = leg_counter + 4
   lastleg(ir_B) = leg_counter + 5
   lastleg(ir_C) = leg_counter + 6
   leg_counter = leg_counter + MAX_GHOSTLEGS   
endif
  
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Four-leg vertex, i.e.  
! Ising operator between i1 and i2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( (i1.gt.0).and.(i2.gt.0).and.(i1.ne.i2) ) then 
! ! IMPROVE: IF( (i1>0) .and. (i2>i1) ) THEN 
#ifdef DEBUG
  ! Check whether i1 < i2:
  if (i1.gt.i2) then
    print*,'Error: i1 > i2 => exiting'
    stop
  endif
#endif    

  ! update linked list 
  ! lower two legs
  if (lastleg(i1).ne.0) then 
     ! double-linked list
     vertexlink(lastleg(i1)) = leg_counter + 1
     vertexlink(leg_counter + 1) = lastleg(i1)
  else
     firstleg(i1) = leg_counter + 1
  endif
  
  if (lastleg(i2).ne.0) then
     vertexlink(lastleg(i2)) = leg_counter + 2
     vertexlink(leg_counter + 2) = lastleg(i2)    
  else
     firstleg(i2) = leg_counter + 2
  endif  

  ! upper two legs
  ! Note the convention for the labelling of legs 
  ! around a 4-leg vertex
  lastleg(i2) = leg_counter + 5
  lastleg(i1) = leg_counter + 4  
  ! mark the 'ghostlegs' as 'visited' so that they are 
  ! not used as starting legs for constructing a cluster
  ! during the off-diagonal update 
  leg_visited(leg_counter + 3) = .TRUE.
  leg_visited(leg_counter + 6) = .TRUE.
  ! increment leg counter using 'ghostlegs' so that 
  ! ever vertex is associated with MAX_GHOSTLEGS legs.  
  leg_counter = leg_counter + MAX_GHOSTLEGS
endif

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Two-leg vertex, i.e.
! spin-flip operator or constant operator at i1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IF( (i1 /= 0).and.((i2 == 0).or.(i2 == i1)) ) THEN 
  ! lower leg
  IF(lastleg(i1).ne.0) THEN
    vertexlink(lastleg(i1)) = leg_counter + 1
    vertexlink(leg_counter + 1) = lastleg(i1)
  ELSE
    firstleg(i1) = leg_counter + 1
  ENDIF  
  ! upper leg
  lastleg(i1) = leg_counter + 4  
  ! mark 'ghostlegs' as 'visited'
  leg_visited(leg_counter + 2) = .TRUE.
  leg_visited(leg_counter + 3) = .TRUE.
  leg_visited(leg_counter + 5) = .TRUE.
  leg_visited(leg_counter + 6) = .TRUE.
  leg_counter = leg_counter + MAX_GHOSTLEGS
ENDIF

#ifdef DEBUG 
!------------------------------------
if ((i1.eq.0).and.(i2.ne.0)) then
  print*, "Error: i1.eq.0 and i2.ne.0"
  stop
endif
!------------------------------------
#endif

enddo

! ***************************************
! Implement periodic boundary conditions 
! in imaginary time for the linked lis. 
! ***************************************
do i=1, config%n_sites
 if (lastleg(i).ne.0) then
    vertexlink(lastleg(i)) = firstleg(i)
    vertexlink(firstleg(i)) = lastleg(i)
  endif
enddo

! TODO: Traverse linked list and check for missing links
! TODO: Also check that there are no identities among the operators 
! connected by the linked list. 
#if defined(DEBUG_CLUSTER_UPDATE)
  leg_counter = 0
  do i=1, config%n_ghostlegs
    if (vertexlink(i) /= -1) then 
      leg_counter = leg_counter + 1 
      ip = (i-1)/MAX_GHOSTLEGS + 1
      if (i /=  vertexlink(vertexlink(i))) then 
        print*, "ERROR: linked list is inconsistent"
        print*, "leg=",i, "vertexlink(vertexlink(i))=",vertexlink(vertexlink(i))
        stop
      endif 
    endif 
  enddo 
  if (leg_counter /= config%n_legs) then 
    print*, "ERROR: linked list has missing links"
    print*, "leg_counter /= config%n_legs"
    stop
  endif 
#endif 


END SUBROUTINE build_linkedlist_plaquette


SUBROUTINE unit_test

END SUBROUTINE unit_test

END MODULE linked_list 
