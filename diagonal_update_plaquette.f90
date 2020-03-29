! TODO:
!   - Implement translation invariance of probability tables 
!     and array indicating whether a bond is FM or AFM.

module diagonal_update
  use types
  implicit none 

type t_ProbTable 
  ! cumulative probability tables for 
  !  - the first index (=first site on which the Ising operator acts)
  !  - the second index (given selection of the first index)
  ! of Ising (4-leg) bond operators 
  real(dp), allocatable :: P_cumulfirst(:)
  real(dp), allocatable :: P_cumulsecond(:,:)
  real(dp) :: sum_all_diagmatrix_elements_2or4leg
  real(dp) :: sum_all_diagmatrix_elements
end type 

  contains 

SUBROUTINE diagonal_update_plaquette( beta, Jij_sign, &
     spins, opstring, config, probtable, plaquettes, update_type )
! **************************************************
! Attempt to carry out replacements of the kind
!     id <--> diag. operator
! at every propagation step.
! **************************************************

use SSE_configuration
use lattice

! use probtables 
implicit none

real(dp), intent(in) :: beta        ! inverse temperature 
integer, intent(in) :: Jij_sign(:,:)  ! sign of the interaction bond (i,j): FM (<0) or AFM (>0)
                                      ! Note: If J_ij(i,j) == 0, then the corresponding bond will never be sampled
                                      !       from the cumulative probability table. 
integer, intent(in) :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Config), intent(inout) :: config
type(t_ProbTable), intent(in) :: probtable
type(t_Plaquette), intent(in) :: plaquettes(:)
integer, intent(in) :: update_type   ! if update_type = 1 => "A-site update, =2 => "B-site update", =3 => "C-site update"

integer :: ip, i1, i2, index1, index2, k

integer :: plaq_idx

real(dp) :: P_plus, P_minus, P_add, P_remove
real(dp) :: prob

! propagated spin configuration at a given propagation step
integer, allocatable :: spins2(:) 

logical :: FOUND, OP_INSERTED

! for plaquette-based cluster update 
integer :: ir_A, ir_B, ir_C


ALLOCATE(spins2( SIZE(spins,1) ))

spins2(:) = spins(:)

! Initialize P_add and P_remove
P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )

do ip=1, config%LL

  i1 = opstring(ip)%i 
  i2 = opstring(ip)%j
  
!identity encountered
IF( (i1 == 0).AND.(i2 == 0) ) THEN

    ! heat bath solutions to the detailed balance condition equations
    P_plus = P_add    
    call random_number(prob)
        
  IF( prob <= P_plus ) THEN 
  ! try to insert an operator  
  OP_INSERTED = .FALSE.
 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Choice among different classes of diagonal operators 
   ! First, decide whether to insert 
   !     - a triangular plaquette operator
   !     - or a diagonal 2-leg or 4-leg vertex.
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call random_number(prob)
   
   IF( prob <= (probtable%sum_all_diagmatrix_elements_2or4leg / probtable%sum_all_diagmatrix_elements) ) THEN
  ! Try inserting an Ising operator or a constant at propagation step ip
  ! assuming that all insertions are allowed.
    
    ! Heat bath algorithm for first index (=first site on which the Ising operator acts)
    call random_number(prob)
    index1 = binary_search( probtable%P_cumulfirst(:), prob )
            
    ! Heat bath algorithm for the second index given the selection of the first index
    call random_number(prob)
    index2 = binary_search( probtable%P_cumulsecond(index1, :), prob )
    

! Check whether the insertion of the Ising operator is allowed at propagation step ip,
! i.e. whether the FM or AFM nature of the Ising operator is compatible with the
! spin configuration at ip. If it is forbidden, repeat the selection process until an operator 
! is inserted. There is no risk of the search never ending as a constant can always be inserted (for h unequal 0 !).

    if (index1.ne.index2) then	!Ising operator, index1 and index2 not equal 0 by construction
    if ( Jij_sign(index1, index2).gt.0 ) then !AFM
      if (spins2(index1).ne.spins2(index2)) then
        OP_INSERTED = .TRUE.
        if (index1.lt.index2) then
          opstring(ip)%i = index1
          opstring(ip)%j = index2
        else
          opstring(ip)%i = index2
          opstring(ip)%j = index1
        endif
        ! 	Expansion order changes as n_exp -> n_exp + 1
	      config%n_exp = config%n_exp + 1; config%n4leg = config%n4leg + 1
        ! update P_add and  P_remove
        P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
        P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
    endif
    elseif ( Jij_sign(index1,index2).lt.0 ) then !FM	
      if (spins2(index1).eq.spins2(index2)) then
        OP_INSERTED = .TRUE.
        if (index1.lt.index2) then
          opstring(ip)%i = index1
          opstring(ip)%j = index2
        else
          opstring(ip)%i = index2
          opstring(ip)%j = index1
        endif
        config%n_exp = config%n_exp + 1; config%n4leg = config%n4leg + 1
        ! update P_add and  P_remove
        P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
        P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
    endif
    endif
    else !index1.eq.index2 => constant operator to be inserted
    ! There is no constraint on inserting constants
        OP_INSERTED = .TRUE.
        opstring(ip)%i = index1
        opstring(ip)%j = index2

        config%n_exp = config%n_exp + 1; config%n2leg = config%n2leg + 1	
        ! update P_add and  P_remove
        P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
        P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
   endif !index1.ne.index2

   else
    ! Try to insert a triangular plaquette operator 
    ! Select the plaquette number 
    call random_number(prob)

    ! IMPROVE: There must be a simpler expression if all plaquettes have the 
    ! same weight. 
    FOUND = .FALSE.; k=1
    do while (.not.FOUND)
      if (k <= (config%n_plaquettes * prob) ) then
        k=k+1
      else
        plaq_idx = k
        FOUND=.TRUE.
      endif
    enddo  
  
  ! Check whether the spin configuration is minimally frustrated so that 
  ! the insertion of a plaquette operator is allowed.
  ! Plaquettes on maximally frustrated spin configurations have zero weight (see Ref. [1]).
  
  ! Choose between A-site, B-site, or C-site update 
  ! by identifying a particular sublattice as the "A-sites", which are 
  ! by definition the 'privileged' sites (see Ref. [1]).
  ! It is only here that the update type of the plaquette update enters. 
  if (update_type == A_UPDATE) then
    ir_A = plaquettes(plaq_idx)%Asite
    ir_B = plaquettes(plaq_idx)%Bsite
    ir_C = plaquettes(plaq_idx)%Csite
  elseif (update_type == B_UPDATE) then 
    ir_A = plaquettes(plaq_idx)%Csite
    ir_B = plaquettes(plaq_idx)%Asite
    ir_C = plaquettes(plaq_idx)%Bsite  
  elseif (update_type == C_UPDATE) then
    ir_A = plaquettes(plaq_idx)%Bsite
    ir_B = plaquettes(plaq_idx)%Csite
    ir_C = plaquettes(plaq_idx)%Asite  
  else
     print*, "Error in diagonal_update: Invalid update_type"
     print*, "Exiting ..."
     stop
  endif 
    
  if (abs(spins2(ir_A) + spins2(ir_B) + spins2(ir_C)) .eq. 1) then 
  ! The spin config is minimally frustrated so that the plaquette operator 
  ! can be inserted. 
     OP_INSERTED = .TRUE.
     ! %i and %j are assigned negative values in order to differentiate 
     ! between Ising bond operators and plaquette operators  
     opstring(ip)%i = -ir_A  
     opstring(ip)%j = -ir_B
     opstring(ip)%k = -ir_C
     
     ! Is the A-site a majority site ?
     if (spins2(ir_B).ne.spins2(ir_C)) then
       ! A-site takes part in majority spin configuration 
       opstring(ip)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .TRUE.   
     else
       opstring(ip)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .FALSE.
     endif     
     config%n_exp = config%n_exp + 1; config%n6leg = config%n6leg + 1 
     ! update P_add and  P_remove
     P_remove = float(( config%LL - config%n_exp + 1)) &
      / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
     P_add = beta*probtable%sum_all_diagmatrix_elements &
      / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )

  endif 
  
   endif ! plaquette or (2-leg or 4-leg) vertex ?

  endif !if(prob <= P_plus)
  
endif !identity encountered

  ! Ising bond operator, Ising triangular plaquette, or constant encountered: 
  ! Try replacement:  Ising plaquette / Ising bond / const => identity
  ! Expansion order changes as n_exp -> n_exp - 1
  if ((i1.ne.0).and.(i2.ne.0)) then

    P_minus = P_remove
    call random_number(prob)
     
    if (prob <= P_minus) then
      opstring(ip)%i = 0
      opstring(ip)%j = 0
      config%n_exp = config%n_exp - 1
      ! update P_add and  P_remove
      P_remove = float(( config%LL - config%n_exp + 1)) &
        / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
      P_add = beta*probtable%sum_all_diagmatrix_elements &
        / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
  
    if (i1 .lt. 0) then ! Triangular plaquette operators are identified by opstring(ip)%i < 0 and opatring(ip)%j < 0. 
       config%n6leg = config%n6leg - 1  ! plaquette operator removed
    else ! constant or plaquette operator 
      if ( i1.ne.i2 ) then ! Ising operator removed
        config%n4leg = config%n4leg - 1
      elseif (i1.eq.i2) then ! constant removed ! IMPROVE: This "elseif" can be replace directly by "else"
        config%n2leg = config%n2leg - 1
      else
          STOP "diagonal_update(): ERROR: trying to remove unknown operator type"
      endif 
    endif  
      
    endif
  endif

  ! Propagate spins as spin-flip operators are encountered.
  if ((i1.ne.0).and.(i2.eq.0)) then
    spins2(i1) = -spins2(i1)
  endif

enddo !do ip=1,LL

deallocate(spins2)

END SUBROUTINE diagonal_update_plaquette 


PURE FUNCTION binary_search(cumul_prob_table, prob) RESULT(idx)
! ***********************************************************************
! Purpose:
! --------
! Sample the normalized discrete probability distribution 
! given in the array `cumul_prob_table(:)` 
! via a binary search. 
! 
! Input:
! ------
!    cumul_prob_table(:): 1D array of normalized cumulative 
!          probabilities. The index into the array starts at 1. 
!    IMPROVE: use intrinsic function lbound() and ubound()
!    prob: probability  0 < prob < 1 
! 
! Return:
! ------- 
!    idx:  The sampled index.
! ***********************************************************************

    USE types
    IMPLICIT NONE    
    REAL(dp), INTENT(IN) :: cumul_prob_table(:)
    REAL(dp), INTENT(IN) :: prob    
    INTEGER :: idx
    
    LOGICAL :: FOUND
    INTEGER :: k1, k2, k
            
#ifdef DEBUG
    IF(SUM(cumul_prob_table(:)) /= 1.d0) STOP &
          "binary_search(): ERROR: cumul. prob. table is not normalized."
#endif 
    FOUND = .FALSE.
    k1=1
    k2=SIZE(cumul_prob_table, 1);
    DO WHILE( .NOT.FOUND )
      k = INT(REAL(k1 + k2)/2.0)
      IF( cumul_prob_table(k) <= prob ) THEN
        k1 = k ! k2 remains the same
        IF( (k2-k1) <= 1 ) THEN
          FOUND = .TRUE.
          idx = k2
        ENDIF
      ELSE
        k2=k ! k1 remains the same
        IF( (k2-k1) <= 1 ) THEN
          FOUND = .TRUE.
          ! The case prob < cumul_prob_table(k=1) is special 
          ! because k=1 cannot be bracketed by k1 and k2.
          IF( prob < cumul_prob_table(1) ) THEN
            idx = 1
          ELSE
            idx = k2
          ENDIF
        ENDIF
      ENDIF    
    ENDDO

END FUNCTION binary_search

end module 