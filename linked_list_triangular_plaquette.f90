MODULE linked_list

!IMPROVE: move type definition to another module 
! 
! Configuration of the simulation cell 
! for a given SSE configuration 
TYPE :: tConfig 
   INTEGER :: n_sites   ! Number of lattice sites 
   INTEGER :: LL        ! Total number of SSE propagation steps (including idientities)
   INTEGER :: n_legs    ! Total number of vertex legs: n_legs=6*n6leg+4*n4leg+2*n2leg
   INTEGER :: n6leg     ! Number of triangular plaquette (i.e. 6-leg) vertices
   INTEGER :: n4leg     ! Number of Ising (i.e. 4-leg) vertices 
   INTEGER :: n2leg     ! Number of constant or spin-flip (i.e. 2-leg) vertices 
END TYPE 

    CONTAINS

SUBROUTINE build_linkedlist_plaquette( vertexlink, opstring, config )
!---------------------------------------------------------------!
! convention for numbering of vertex legs:                      !
!                                                               !
!   Ising operator:		Spin-flip operator or constant: !
!  (3)           (4)  [5,6]           (2)   [3,4,5,6]           !
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
! Note: "ghostlegs" in square brackets are not connected to     !
!       any other legs.                                         !
!---------------------------------------------------------------!
IMPLICIT NONE 

INTEGER, POINTER, INTENT(OUT) :: vertexlink(:)
TYPE(tBondOperator), POINTER, INTENT(IN) :: opstring(:)
TYPE(tConfig), INTENT(IN) :: config  

! GENERALIZE: This depends on the type of frustrated plaquettes.
INTEGER, PARAMETER :: MAX_GHOSTLEGS = 6

! automatic arrays used for building the linked list 
INTEGER :: firstleg( config%n_sites )
INTEGER :: lastleg( config%n_sites ) 

integer :: i, ip, i1, i2, leg_counter
integer :: ir_A, ir_B, ir_C

if(allocated(vertexlink)) deallocate(vertexlink)
allocate(vertexlink( config%n_legs )) 
firstleg(:) = 0
lastleg(:) = 0 

! The first leg has index 1.
leg_counter = 0

! **********************************************
! Traverse operator list leaving out identities
! **********************************************
do ip=1, config%LL

! The structure components i1 and i2 are needed to determine the type 
! of vertex.
  i1 = opstring(ip)%i  
  i2 = opstring(ip)%j

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
#ENDIF   

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
  lastleg(i2) = leg_counter + 4
  lastleg(i1) = leg_counter + 3  
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
  lastleg(i1) = leg_counter + 2  
  leg_counter = leg_counter + MAX_GHOSTLEGS
ENDIF

#ifdef DEBUG 
!------------------------------------
if ((i1.eq.0).and.(i2.ne.0)) then
  print*, "Error: i1.eq.0 and i2.ne.0"
  stop
endif
!------------------------------------
#endif DEBUG 

enddo

! ***************************************
! Implement periodic boundary conditions 
! in imaginary time for the linked lis. 
! ***************************************
do i=1,n_sites
 if (lastleg(i).ne.0) then
    vertexlink(lastleg(i)) = firstleg(i)
    vertexlink(firstleg(i)) = lastleg(i)
  endif
enddo

END SUBROUTINE build_linkedlist_plaquette


SUBROUTINE unit_test
! ! ! #ifdef DEBUG

! ! ! integer :: l, leg, ir, i3
! ! ! integer, allocatable :: spins_tmp(:)

! ! !     allocate(spins_tmp(n_sites))
! ! ! 
! ! !     do ir=1,n_sites
! ! !       spins_tmp(ir) = spins(ir)
! ! !     enddo
! ! !     
! ! !     print*, "OPERATOR LIST      ip          i1            i2"
! ! !     do ip=LL,1, -1
! ! ! 	i1 = opstring(ip)%i
! ! ! 	i2 = opstring(ip)%j
! ! !       if ( (i1.eq.0).and. (i2.eq.0) ) then
! ! ! 	print*, "identity", ip, i1, i2
! ! !       elseif ( (i1.gt.0).and.(i2.eq.0) ) then
! ! ! 	print*, "spinflip", ip, i1, i2
! ! ! 	spins_tmp(i1) = - spins_tmp(i1)
! ! !       elseif ( (i1.gt.0).and.(i1.eq.i2) ) then
! ! ! 	print*, "constant", ip, i1, i2
! ! !       elseif ( (i1.gt.0).and.(i2.gt.0).and.(i1.ne.i2) ) then
! ! ! 	print*, "Ising   ", ip, i1, i2
! ! !       elseif ( (i1.lt.0).and.(i2.lt.0) ) then
! ! !         i3 = opstring(ip)%k
! ! !         print*, "Plaquette", ip, abs(i1), abs(i2), abs(i3)
! ! !       else
! ! ! 	print*, "strange operator detected, ip=", ip
! ! !       endif
! ! !       !print current spin configuration
! ! !   !     print*, "spins", ip, (spins_tmp(ir), ir=1,n_sites)
! ! !     enddo
! ! ! 
! ! !     deallocate(spins_tmp)
! ! ! 
! ! ! ! display linked list
! ! !     print*, "LINKED LIST DONE - > display"
! ! !     do leg =1,n_legs
! ! !       print*, leg, " --> ", vertexlink(leg)
! ! !     enddo
! ! ! 
! ! ! #endif 
! ! ! 
! ! ! ! TO DO: Traverse linked list and check for missing links
! ! ! ! 	Make sure that the linked list is closed.
! ! ! ! 	Check that the arrays ipir_to_leg(:) and  leg_to_ipir(:) are consistent.
! ! ! ! REMOVE ALL THE CHECK ROUTINES ONCE THE PROGRAM IS USED FOR COMPUTATIONS

END SUBROUTINE unit_test

END MODULE linked_list 
