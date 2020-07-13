! TODO:
!   - Implement translation invariance of probability tables 
!     and array indicating whether a bond is FM or AFM.

module diagonal_update
  use types
  implicit none 
  private 

  public t_ProbTable
  public diagonal_update_plaquette
  public init_probtables
  public extend_cutoff

type t_ProbTable 
  ! cumulative probability tables for 
  !  - the first index (=first site on which the Ising operator acts)
  !  - the second index (given selection of the first index)
  ! of Ising (4-leg) bond operators 
  real(dp), allocatable :: P_cumulfirst(:)
  real(dp), allocatable :: P_cumulsecond(:,:)
  real(dp) :: sum_all_diagmatrix_elements_2or4leg
  real(dp) :: sum_all_diagmatrix_elements
  real(dp) :: consts_added 
end type 

  contains 

SUBROUTINE diagonal_update_plaquette( S, beta, Jij_sign, &
     spins, opstring, config, probtable, plaquettes, update_type, &
     TRANSLAT_INVAR )
! **************************************************
! Attempt to carry out replacements of the kind
!     id <--> diag. operator
! at every propagation step.
! **************************************************

use SSE_configuration
use lattice

! use probtables 
implicit none

integer, parameter :: A_UPDATE=111, B_UPDATE=112, C_UPDATE=113

type(Struct), intent(in) :: S                     ! lattice structure object 
real(dp), intent(in) :: beta          ! inverse temperature 
integer, allocatable, intent(in) :: Jij_sign(:,:)  ! sign of the interaction bond (i,j): FM (<0) or AFM (>0)
                                      ! Note: If J_ij(i,j) == 0, then the corresponding bond will never be sampled
                                      !       from the cumulative probability table. 
integer, intent(in) :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Config), intent(inout) :: config
type(t_ProbTable), intent(in) :: probtable
type(t_Plaquette), intent(in) :: plaquettes(:)
! if update_type == A_UPDATE => "A-site update
integer, intent(in) :: update_type   
logical, intent(in) :: TRANSLAT_INVAR

! ... Local variables ...
integer :: ip, i1, i2, index1, index2
integer :: r(2)

integer :: plaq_idx

real(dp) :: P_plus, P_minus, P_add, P_remove
real(dp) :: prob

! propagated spin configuration at a given propagation step
integer, allocatable :: spins2(:) 

logical :: OP_INSERTED

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
    
    ! Heat bath algorithm for first index (=first site on which the Ising operator / const acts)
    call random_number(prob)
    if( TRANSLAT_INVAR ) then 
      index1 = ceiling( config%n_sites * prob )
    else
      index1 = binary_search( probtable%P_cumulfirst(:), prob )
    endif 
            
    ! Heat bath algorithm for the second index given the selection of the first index
    call random_number(prob)
    index2 = binary_search( probtable%P_cumulsecond(index1, :), prob )
    

! Check whether the insertion of the Ising operator is allowed at propagation step ip,
! i.e. whether the FM or AFM nature of the Ising operator is compatible with the
! spin configuration at ip. If it is forbidden, repeat the selection process until an operator 
! is inserted. There is no risk of the search never ending as a constant can always be inserted (for h unequal 0 !).

    if (index1.ne.index2) then	!Ising operator, index1 and index2 not equal 0 by construction
      ! for translationally invariant system 
      r = translate_kagome( S, index1, index2 ) ! IMPROVE: so far, only kagome 
      if ( Jij_sign(r(1), r(2)).gt.0 ) then !AFM
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
#ifdef DEBUG_DIAGONAL_UPDATE
    print*, "insert an AFM Ising operator"
#endif       
          config%n_exp = config%n_exp + 1; config%n4leg = config%n4leg + 1
          ! update P_add and  P_remove
          P_remove = float(( config%LL - config%n_exp + 1)) &
            / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
          P_add = beta*probtable%sum_all_diagmatrix_elements &
            / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
        endif
      else !FM; Note that index1 and index2 such that Jij_sign(index1,index2)==0 can never be draw from the prob. tables.
        if (spins2(index1).eq.spins2(index2)) then
          OP_INSERTED = .TRUE.
          if (index1.lt.index2) then
            opstring(ip)%i = index1
            opstring(ip)%j = index2
          else
            opstring(ip)%i = index2
            opstring(ip)%j = index1
          endif
#ifdef DEBUG_DIAGONAL_UPDATE
    print*, "insert an FM Ising operator"
#endif             
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
#ifdef DEBUG_DIAGONAL_UPDATE
  print*, "insert a const. operator"
#endif 
        config%n_exp = config%n_exp + 1; config%n2leg = config%n2leg + 1	
        ! update P_add and  P_remove
        P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
        P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
    endif !index1.ne.index2

   ELSE
    ! Try to insert a triangular plaquette operator 
    ! Select the plaquette number 
    call random_number(prob)
    plaq_idx = ceiling(config%n_plaquettes * prob)
  
  ! Check whether the spin configuration is minimally frustrated so that 
  ! the insertion of a plaquette operator is allowed.
  ! Plaquettes on maximally frustrated spin configurations have zero weight (see Ref. [1]).
  
  ! Choose between A-site, B-site, or C-site update 
  ! by identifying a particular sublattice as the "A-sites", which are 
  ! by definition the 'privileged' sites (see Ref. [1]).
  ! It is only here that the update type ("A-update", "B-update", or "C-update") 
  ! of the plaquette update enters. 
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
#ifdef DEBUG_DIAGONAL_UPDATE
  print*, "insert a minimally frustrated plaquette operator, plaq_idx=", plaq_idx, "ut=", update_type
#endif          
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

  ! Diagonal operator encountered: 
  ! Try replacement:  Ising plaquette / Ising bond / const => identity
  ! Expansion order changes as n_exp -> n_exp - 1
  if ((i1.ne.0).and.(i2.ne.0)) then

    P_minus = P_remove
    call random_number(prob)
     
    if (prob <= P_minus) then
      opstring(ip)%i = 0
      opstring(ip)%j = 0
#ifdef DEBUG_DIAGONAL_UPDATE
  print*, "remove a diagonal operator"
#endif                 
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
      elseif (i1.eq.i2) then ! constant removed 
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

! Update important variables that have not been explicitly 
! updated
    config%n_legs = 2*config%n2leg+4*config%n4leg+6*config%n6leg
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL

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


subroutine init_probtables( S, J_interaction_matrix, hx, &
    probtable, J_1, n_plaquettes, TRANSLAT_INV )
! *******************************************************************
! Purpose:
! --------
!    Precompute cumulative probability tables needed for sampling 
!    SSE operators during the diagonal update.

! Input:
! ------
!    S: lattice structure object
!    J_interaction_matrix: 2D array of shape (n_sites, n_sitess)
!       The element [i,j] contains the 
!       interaction between linearly stored sites i and j.
!    hx: value of the transverse field 
!    n_plaquettes: number of triangular plaquettes
!    TRANSLAT_INV: logical
!       Indicates whether the interactions are translationally 
!       invariant. In this case the interaction matrix can be 
!       replaced by a vector. 
!       !!! Not implemented yet !!!
!
! Output:
! -------
!    probtable: structure of type(t_ProbTable)
!       Precomputed cumulative probability tables
!    Jij_sign: integer array
!       Sign of the interaction bonds. 
! 
! *******************************************************************
  use SSE_configuration 
  use lattice
  implicit none 

  type(Struct), intent(in)       :: S
  real(dp), intent(in)           :: J_interaction_matrix(:,:)
  real(dp), intent(in)           :: hx
  type(t_ProbTable), intent(out) :: probtable 
  real(dp), intent(in)           :: J_1  ! nearest neighbour interactions inside a plaquette, subroutine expects: J_1 = +1
  integer,intent(in)             :: n_plaquettes 
  logical, intent(in)            :: TRANSLAT_INV

  ! automatic helper arrays
  ! array of the matrix elements of the bond operators, i.e. h for i=j and 2|J_ij| else
  real(dp) :: M_ij( size(J_interaction_matrix, 1), size(J_interaction_matrix, 1) )
  ! relative probabilities for choosing the first index of an Ising operator
  real(dp) :: P_first( size(J_interaction_matrix, 1) )
  integer :: n
  real(dp) :: cc, ss

  integer :: ir, jr, k
  real(dp) :: norm

  if( J_1 /= +1 ) then 
    print*, "For the triangular plaquette update to be meaningful we need J_1 = +1 (i.e. AFM)."
    print*, "Furthermore, J_1 sets the energy scale."
    if (J_1 == 0) then 
      print*, "It is assumed that you set J_1 = 0 in order to suppress the plaquette update" 
      print*, "and only use the bond-based update. The interaction matrix for pairwise interacionts"
      print*, "should be read from some file."
      print*, "Continuing ..."
    else
      print*, "Exiting ..."
      stop
    endif 
  endif 


  if( TRANSLAT_INV ) then 
    print*, "Translational invariance of probability tables not implemented yet."
    print*, "Exiting ..."
    stop

  endif 

  ! number of lattice sites 
  n = size(J_interaction_matrix, 1)
  if( n /= S%Nsites ) then 
    print*, "ERROR: init_probtables(): n /= S%Nsites "
    stop
  endif 

  ! Hamiltonian matrix elements on the computational, i.e. the Sz - basis
  ! Matrix elements have only two values, TWO*abs(J_ij) and hx:
  ! The factor TWO is not necessary if a bond can be represented both as (i1, j1) and (j1, i1),
  ! which is the choice taken in this code.  
  M_ij(:,:) = dabs(J_interaction_matrix(:,:)) ! TWO*dabs(J_interaction_matrix(:,:)) 

  do ir=1,n
    M_ij(ir,ir) = hx
  enddo

! 1. The constants which have been added to the Hamiltonian artificially have to be
! subtracted from the energy in the end.
! 2. Sum over all matrix elements used in the transition probabilities 
! when deciding whether to insert or remove an operator

cc = n*hx
ss = n*hx
do ir=1, n
  do jr=1, ir-1 ! avoid double counting and count only ir.ne.jr
    cc = cc + abs(J_interaction_matrix(ir,jr)) 
    ss = ss + TWO*abs(J_interaction_matrix(ir,jr))
  enddo
enddo

! see Ref. [1]
if (trim(S%lattice_type) == "triangular") then 
  ! Triangular lattice 
  cc = cc + (3.0_dp/2.0_dp) * abs(J_1) * n_plaquettes
elseif (trim(S%lattice_type) == "kagome") then 
  ! Kagome Llattice 
  cc = cc + 3.0_dp * abs(J_1) * n_plaquettes
else
  stop "init_probtables(): Unknown lattice type."
endif 

probtable%consts_added = cc
probtable%sum_all_diagmatrix_elements_2or4leg = ss 
! include matrix elements from triangular plaquettes (of the Hamiltonian shifted by the constants)
if (trim(S%lattice_type) == "triangular") then 
  probtable%sum_all_diagmatrix_elements = ss + 2.0_dp * abs(J_1) * n_plaquettes
elseif (trim(S%lattice_type) == "kagome") then 
  probtable%sum_all_diagmatrix_elements = ss + 4.0_dp * abs(J_1) * n_plaquettes
else
  stop "init_probtables(): Unknown lattice type."
endif 

! Calculate the cumulative probability tables used in the diagonal update
! for the insertion of diagonal 2leg and 4leg vertices. 
allocate(probtable%P_cumulfirst(n))
allocate(probtable%P_cumulsecond(n,n))

do ir = 1, n
  P_first(ir) = 0.0_dp; probtable%P_cumulfirst(ir) = 0.0_dp
  do k = 1, n
    probtable%P_cumulsecond(ir,k) = 0.0_dp
  enddo
enddo

! probabilities for selecting the first index of an Ising/constant operator
! relative probability for selecting index i:
do ir = 1, n
  P_first(ir) = sum(abs(M_ij(ir,:)))
enddo

! cumulative probabilities used in the heat bath algorithm
norm = sum(P_first(:))
do k = 1, n 
  do ir = 1, k
    probtable%P_cumulfirst(k) = probtable%P_cumulfirst(k) + P_first(ir)
  enddo
  probtable%P_cumulfirst(k) = probtable%P_cumulfirst(k)/norm 
enddo

! The relative probability for choosing j as the second index given i as the first index
! is M_ij (conditional probability). 
do ir = 1, n  
  do k = 1, n
    do jr = 1, k
      probtable%P_cumulsecond(ir,k) = probtable%P_cumulsecond(ir,k) + abs(M_ij(ir,jr))
    enddo
  enddo
enddo

!normalize
do ir = 1, n
  norm = 0.0_dp
  do jr = 1, n
    norm = norm + abs(M_ij(ir,jr))
  enddo
  do k = 1, n
    probtable%P_cumulsecond(ir,k) = probtable%P_cumulsecond(ir,k) / norm
  enddo
enddo

! !REMOVE
! print*, "P_cumulfirst", probtable%P_cumulfirst(:)
! print*, "P_cumulsecond"
! do ir = 1,n
!   print*, ir, probtable%P_cumulsecond(ir,:)
! enddo
! stop
! !REMOVE

end subroutine 


subroutine extend_cutoff(opstring, config)
!****************************************************************
! Extend the cut-off of the fixed length operator string to
!   LL_new = n_exp + n_exp / 2. 
!****************************************************************
  use SSE_configuration 
  use types 
  implicit none 

  type(t_BondOperator), allocatable, intent(inout) :: opstring(:)
  type(t_Config), intent(inout) :: config 
  
  logical, allocatable :: p_taken(:)
  type(t_BondOperator), allocatable :: opstring_new(:)
  integer :: LL, LL_new, n_exp_new
  integer :: i, ip, pos, n    
  real(dp) :: eta

  LL = config%LL
  LL_new = config%n_exp + int(float(config%n_exp)/2.0_dp)
  
  allocate(p_taken(LL_new))
  allocate(opstring_new(LL_new))
  do ip = 1, LL_new
      p_taken(ip) = .FALSE.
  enddo
  
  ! Extract a random sequence of (LL_new - LL) positions
  ! in the interval [0,LL_new] where
  ! the identity operators will be inserted. 
  i=1
  do while (i <= (LL_new-LL))
    call random_number(eta)
    pos = (int(eta*LL_new) + 1) 
    if (.not.p_taken(pos)) then
      p_taken(pos) = .TRUE.
      i=i+1
    endif
  enddo
  
  n=1   
  do ip = 1, LL_new
    if (p_taken(ip)) then
      ! insert an additional identity
      opstring_new(ip)%i = 0
      opstring_new(ip)%j = 0
    else  ! carry over the old operator string
      opstring_new(ip) = opstring(n)
      n = n + 1
    endif
  enddo
  
  deallocate(opstring); allocate(opstring(LL_new))
  
  do ip = 1, LL_new
    opstring(ip) = opstring_new(ip)
  enddo
  
  ! Check that the expansion order has been maintained
  n_exp_new = 0
  do ip = 1, LL_new
    if( .not.((opstring(ip)%i == 0).and.(opstring(ip)%j == 0)) ) then
      ! non-trivial operator encountered
      n_exp_new = n_exp_new + 1
    endif
  enddo
  
  if (n_exp_new /= config%n_exp) then
    print*, "subroutine increase_cutoff: n_exp_new = ", n_exp_new, "  n_exp = ", config%n_exp
    print*, "Exiting ..."
    stop
  endif
  
  ! Update the only entries of the structure 'config'
  ! which have changed.
  config%LL = LL_new
  config%n_ghostlegs = MAX_GHOSTLEGS*config%LL
  
  deallocate(p_taken); deallocate(opstring_new)
  
end subroutine extend_cutoff

end module 