! TODO:
!   - Implement translation invariance of probability tables 
!     and array indicating whether a bond is FM or AFM.
!   - if ( Jij_sign(r(1), r(2)).gt.0 ) then !AFM (line 142) 
!     => write such that it works also for non-translationally invariant systems
!   - inhomongeneous longitudinal field, sample in diagonal update w


module diagonal_update
  use types
  use util

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
  real(dp) :: consts_added 
  real(dp) :: sum_all_diagmatrix_elements 
  integer :: n_opclass                           ! number of classes of diagonal operators to be chosen from in the diagonal update
  real(dp), allocatable :: cumprob_opclass(:)    ! cumulative sum of the probabilities to insert an operator from a certain class
  integer, allocatable  :: idx_opclass(:)        ! integer indices for operator classes (used in the diagonal update)
  real(dp), allocatable :: hz_matrix_element(:)  ! ( 2*abs(hz(i)) + 2*C_par )
  real(dp), allocatable :: hz_cumul(:)           ! cumulative probability table, hz_matrix_element(ir) are relative probabilities 

  ! needed for updating the probability tables 
  ! whenever the instantaneous spin configuration changes 
  real(dp), allocatable :: sum_all_diagmatrix_elements_peropclass(:)
  real(dp), allocatable :: consts_added_per_opclass(:)
  real(dp), allocatable :: prob_opclass(:)
end type 

  contains 

pure function choose_diagopclass(eta, probtable) result(optype)
  ! Purpose:
  ! --------
  !
  ! Arguments:
  ! ----------
  real(dp), intent(in) :: eta   ! random number, uniformly distributed on [0,1]
  type(t_ProbTable), intent(in) :: probtable   
  integer :: optype 
  integer :: c

  do c = 1, probtable%n_opclass
    if( eta < probtable%cumprob_opclass(c) ) then 
      optype = probtable%idx_opclass(c)
      exit
    endif 
  enddo 
end function 


SUBROUTINE diagonal_update_plaquette( S, beta, Jij_sign, hz_fields_sign, &
     spins, opstring, config, probtable, plaquettes, update_type, &
     TRANSLAT_INVAR, hz_fields, C_par_hyperparam )
! **************************************************
! Attempt to carry out replacements of the kind
!     id <--> diag. operator
! at every propagation step.
! **************************************************

use SSE_configuration
use lattice
use util, only: assert

! use probtables 
implicit none
integer, parameter :: A_UPDATE=111, B_UPDATE=112, C_UPDATE=113

type(Struct), intent(in) :: S                     ! lattice structure object 
real(dp), intent(in) :: beta          ! inverse temperature 
integer, allocatable, intent(in) :: Jij_sign(:,:)  ! sign of the interaction bond (i,j): FM (<0) or AFM (>0)
                                      ! Note: If J_ij(i,j) == 0, then the corresponding bond will never be sampled
                                      !       from the cumulative probability table. 
integer, intent(in) :: hz_fields_sign(:)  ! sign of the longitudinal fields: FM (<0) or AFM(>0).
integer, intent(in) :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Config), intent(inout) :: config
type(t_ProbTable), intent(inout) :: probtable  ! intent(inout) because probtables can be changed if they depend on instantaneous spin config
type(t_Plaquette), intent(in) :: plaquettes(:)
! if update_type == A_UPDATE => "A-site update
integer, intent(in) :: update_type   
logical, intent(in) :: TRANSLAT_INVAR
real(dp), intent(in) :: C_par_hyperparam
real(dp), intent(in) :: hz_fields(:) 

! ... Local variables ...
integer :: ip, i1, i2, index1, index2, index2_from_origin
integer :: ii1
integer :: r(2)
integer :: sub1
integer :: optype, optype_old

integer :: plaq_idx

real(dp) :: P_add, P_remove
real(dp) :: prob, eta
integer :: alignment 

! REMOVE
real(dp) :: norm 
real(dp) :: cumprob_hz(1:S%Nsites)
integer :: ir
! REMOVE

! propagated spin configuration at a given propagation step
integer, allocatable :: spins2(:) 

! for plaquette-based cluster update 
integer :: ir_A, ir_B, ir_C

! "instantaneous" spin configuration at a given propagation step 
allocate(spins2( 1:size(spins,1) ))
spins2(:) = spins(:)

! Initialize P_add and P_remove given the current expansion order.
P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )

do ip=1, config%LL

  i1 = opstring(ip)%i 
  i2 = opstring(ip)%j
  optype_old = opstring(ip)%optype
  
!identity encountered
if( (i1 == 0).and.(i2 == 0) ) then
    
  call random_number(eta)
        
  if( eta <= P_add ) then
  ! try to insert an operator  

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Choice among different classes of diagonal operators 
   ! First, decide whether to insert 
   !     - a triangular plaquette operator
   !     - or a diagonal 2-leg CONSTANT or 4-leg ISING vertex
   !     - or a diagonal 2-leg LONGITUDINAL operator 
   ! IMPROVE: maybe make 2-leg CONSTANT a separate class ? 
   
   ! choose with heat bath probability 
   call random_number(eta)
   optype = choose_diagopclass(eta, probtable)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   choose_insert_optype: select case( optype )

  case( TWOLEGCONST_OR_FOURLEG )
  ! Try inserting an Ising operator or a constant at propagation step `ip`
  ! assuming that all insertions are allowed.
    
    ! Heat bath algorithm for first index (=first site on which the Ising operator / const acts)
    call random_number(prob)
    if( TRANSLAT_INVAR ) then 
      index1 = ceiling( config%n_sites * prob )
      call random_number(prob)
      sub1 = mod(index1-1, S%Nbasis) + 1
      index2_from_origin = binary_search( probtable%P_cumulsecond(sub1, :), prob )
      index2 = translate_from_origin( S=S, i=index1, j=index2_from_origin )
    else
      index1 = binary_search( probtable%P_cumulfirst(:), prob )
      ! Heat bath algorithm for the second index given the selection of the first index
      call random_number(prob)
      index2 = binary_search( probtable%P_cumulsecond(index1, :), prob )
    endif           
  
! Check whether the insertion of the Ising operator is allowed at propagation step ip,
! i.e. whether the FM or AFM nature of the Ising operator is compatible with the
! spin configuration at ip. If it is forbidden, the replacement
!     identity => Ising
! is not carried out. 
! There is no risk of the search never ending as a constant can always be inserted (for hx unequal 0 !).

    if (index1.ne.index2) then	! Ising operator, index1 and index2 not equal 0 by construction
      if( TRANSLAT_INVAR ) then 
        r = translate_to_origin( S, index1, index2 ) 
      else
        r = (/ index1, index2 /)
      endif 
      if ( Jij_sign(r(1), r(2)).gt.0 ) then !AFM
        if (spins2(index1).ne.spins2(index2)) then
          if (index1.lt.index2) then
            opstring(ip)%i = index1
            opstring(ip)%j = index2
          else
            opstring(ip)%i = index2
            opstring(ip)%j = index1
          endif
          opstring(ip)%optype = ISING_BOND
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
          if (index1.lt.index2) then
            opstring(ip)%i = index1
            opstring(ip)%j = index2
          else
            opstring(ip)%i = index2
            opstring(ip)%j = index1
          endif
          opstring(ip)%optype = ISING_BOND
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
    else  ! index1.eq.index2 => constant operator to be inserted
    ! There is no constraint on inserting constants
        opstring(ip)%i = index1
        opstring(ip)%j = index2
        opstring(ip)%optype = CONSTANT
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

  case( TRIANGULAR_PLAQUETTE )
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
     ! %i and %j are assigned negative values in order to differentiate 
     ! between Ising bond operators and plaquette operators  
     opstring(ip)%i = -ir_A  
     opstring(ip)%j = -ir_B
     opstring(ip)%k = -ir_C
     opstring(ip)%optype = TRIANGULAR_PLAQUETTE
     
     ! Is the A-site a majority site ?
     if (spins2(ir_B).ne.spins2(ir_C)) then
       ! A-site takes part in majority spin configuration 
       opstring(ip)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .TRUE.   
     else
       opstring(ip)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .FALSE.
     endif     
#ifdef DEBUG_DIAGONAL_UPDATE
  print*, "insert a minimally frustrated plaquette operator, plaq_idx=", plaq_idx, "ut=", update_type
  print*, "spins2(ir_A/B/C)=", spins2(ir_A), spins2(ir_B), spins2(ir_C), ir_A, ir_B, ir_C
  print*, "spins2(:)=", spins2(:)
#endif          
     config%n_exp = config%n_exp + 1; config%n6leg = config%n6leg + 1 
     ! update P_add and  P_remove
     P_remove = float(( config%LL - config%n_exp + 1)) &
      / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
     P_add = beta*probtable%sum_all_diagmatrix_elements &
      / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
  endif 
  
  case( LONGITUDINAL )

!     ! -----------------------------------------
!     ! For inhomogeneous longitudinal fields (variant 1)
!     ! Probability tables need to be recomputed after each spin flip. 
!     ! -----------------------------------------
! #ifdef DEBUG_DIAGONAL_UPDATE
!     print*, "insert longitudinal operator"
! #endif     
!     norm = 0.0_dp
!     cumprob_hz(:) = 0.0_dp
!     do ir = 1, S%Nsites 
!          norm = norm + probtable%hz_matrix_element(ir)
!          cumprob_hz(ir) = norm 
!     enddo 
!     cumprob_hz(:) = cumprob_hz(:) / norm

!     call assert(all(cumprob_hz>=0), 'all(cumprob_hz>0) failed')
    
!     call random_number(eta)
!     ii1 = binary_search(cumprob_hz, eta)

!     if ( spins2(ii1) == hz_fields_sign(ii1) ) then           
!         alignment = +1
!     else 
!         alignment = -1
!     endif 

!     opstring(ip)%optype = LONGITUDINAL 
!     opstring(ip)%i = ii1
!     opstring(ip)%j = -ii1         ! Actually, not necessary to set this entry, it is not used. 
!     opstring(ip)%k = alignment   ! The spin state is needed in the cluster update to accumulate the exchange fields.

!     config%n_exp = config%n_exp + 1; config%n2leg_hz = config%n2leg_hz + 1 

!     P_remove = float(( config%LL - config%n_exp + 1)) &
!         / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
!     P_add = beta*probtable%sum_all_diagmatrix_elements &
!         / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )

    ! --------------------------------------------------
    ! For inhomogeneous longitudinal fields (variant 2)
    ! Probability tables need not be recomputed => better 
    ! --------------------------------------------------        
    call random_number(eta)

    ! ! homogeneous (in absolute value) field: pick a site at random 
    ! ii1 = int(eta * config%n_sites) + 1
    ! inhomogeneous (in magnitude and sign) field:

    ii1 = binary_search(probtable%hz_cumul(:), eta)

    ! choose one of the matrix elements of hz at random 
    call random_number(eta)
    if ( eta < ( 2*abs(hz_fields(ii1)) + C_par_hyperparam ) &
        / ( 2*abs(hz_fields(ii1)) + 2*C_par_hyperparam ) ) then 
      alignment = + 1
    else
      alignment = -1
    endif 

    ! check whether the selected operator fits with the spin configuration 
    if ( alignment == spins2(ii1) * hz_fields_sign(ii1) ) then 

#ifdef DEBUG_DIAGONAL_UPDATE
      print*, "insert longitudinal operator"
#endif 
  
      opstring(ip)%optype = LONGITUDINAL 
      opstring(ip)%i = ii1
      opstring(ip)%j = -ii1         ! Actually, not necessary to set this entry, it is not used. 
      opstring(ip)%k = alignment   ! The spin state is needed in the cluster update to accumulate the exchange fields.

      config%n_exp = config%n_exp + 1; config%n2leg_hz = config%n2leg_hz + 1 

      P_remove = float(( config%LL - config%n_exp + 1)) &
          / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
      P_add = beta*probtable%sum_all_diagmatrix_elements &
          / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )

    else
      ! do nothing and go to next propagation step  
    endif     
    ! -------------------------------------     

  case default
      print*, "Diagonal update: Trying to insert unknown diagonal operator type"
      print*, "optype = ", optype
      print*, "Exiting ..."
      stop 

  end select choose_insert_optype

  endif !if(prob <= P_add)
  
endif !identity encountered

  ! Diagonal operator encountered: 
  ! Try replacement:  Ising plaquette / Ising bond / const (hx) / longitudinal (hz) => identity
  ! Expansion order changes as n_exp -> n_exp - 1
  if ((i1.ne.0).and.(i2.ne.0)) then

    call random_number(prob)
     
    if (prob <= P_remove) then
#ifdef DEBUG_DIAGONAL_UPDATE
      print*, "remove a diagonal operator"
#endif      
      optype = opstring(ip)%optype
           
      config%n_exp = config%n_exp - 1
      ! update P_add and  P_remove because n_exp has changed 
      P_remove = float(( config%LL - config%n_exp + 1)) &
        / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
      P_add = beta*probtable%sum_all_diagmatrix_elements &
        / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
  
      select case( optype )
        case( TRIANGULAR_PLAQUETTE )
          config%n6leg = config%n6leg - 1  ! plaquette operator removed
        case( ISING_BOND )
          config%n4leg = config%n4leg - 1
        case( CONSTANT )
          config%n2leg = config%n2leg - 1
        case( LONGITUDINAL )
          config%n2leg_hz = config%n2leg_hz - 1
        case default 
          print*, "optype=", optype
          stop "diagonal_update(): ERROR: trying to remove unknown operator type"
      end select 
      
      opstring(ip)%i = 0
      opstring(ip)%j = 0
      opstring(ip)%k = 0 ! to be sure nothing bad happens
      opstring(ip)%optype = IDENTITY  

    endif
  endif

  ! Off-diagonal operators are not exchanged in the diagonal update.
  ! Propagate spins as spin-flip operators are encountered.
  if ((i1.ne.0).and.(i2.eq.0)) then
    spins2(i1) = -spins2(i1)
!    call update_probtables(S=S, probtable=probtable, spins=spins2, &
!        C_par_hyperparam=C_par_hyperparam, hz_fields=hz_fields)
    ! `sum_all_diagmatrix_elements` has changed => update 
    P_remove = float(( config%LL - config%n_exp + 1)) &
        / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
    P_add = beta*probtable%sum_all_diagmatrix_elements &
        / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )       
  endif

enddo !do ip=1,LL

! Update important variables that have not been explicitly 
! updated
    config%n_legs = 2*config%n2leg+2*config%n2leg_hz+4*config%n4leg+6*config%n6leg
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL

deallocate(spins2)

#ifdef DEBUG_CLUSTER_UPDATE
  print*, "Diagonal update completed"
#endif 

END SUBROUTINE diagonal_update_plaquette 


PURE FUNCTION binary_search(cumul_prob_table, prob) RESULT(idx)
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
!
    USE types
    IMPLICIT NONE    
    REAL(dp), INTENT(IN) :: cumul_prob_table(:)
    REAL(dp), INTENT(IN) :: prob    

! ... Local variables ...
    INTEGER :: idx
    LOGICAL :: FOUND
    INTEGER :: k1, k2, k
            
#ifdef DEBUG
    IF(SUM(cumul_prob_table(:)) /= 1.0_dp) STOP &
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
  probtable, J_1, n_plaquettes, TRANSLAT_INV, &
  hz_fields, C_par_hyperparam, spins )
! Purpose:
! --------
!    Precompute cumulative probability tables needed for sampling 
!    SSE operators during the diagonal update.
!
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
!       replaced by a vector for a Bravais lattice (or by a matrix 
!        of shape (S%Nbasis, S%Nsites) for a lattice with a unit cell).
!
! Note: Longitudinal field operators (of strength hz_fields(:)) have, depending on the spin configuration,
!       two different nonzero matrix elements. Therefore the probability of inserting an hz operator 
!       at a given propagation step depends on the instantaneous spin configuration and the 
!       corresponding part of the probability table needs to be recomputed after each propagation 
!       step with a spin-flip operator (which changes the spin configuration).
!
! Output:
! -------
!    probtable: structure of type(t_ProbTable)
!       Precomputed cumulative probability tables
!    Jij_sign: integer array
!       Sign of the interaction bonds. 
use SSE_configuration 
use lattice
implicit none 

! Arguments:
! ----------
type(Struct), intent(in)       :: S
real(dp), intent(in)           :: J_interaction_matrix(:,:)
real(dp), intent(in)           :: hx
type(t_ProbTable), intent(out) :: probtable 
real(dp), intent(in)           :: J_1  ! nearest neighbour interactions inside a plaquette, subroutine expects: J_1 = +1
integer,intent(in)             :: n_plaquettes 
logical, intent(in)            :: TRANSLAT_INV
real(dp), intent(in)           :: hz_fields(:)
real(dp), intent(in)           :: C_par_hyperparam
integer, intent(in)            :: spins(:)  ! instantaneous spin configuration at a given propagation step (needed for hz_fiels(:)!=0)

! ... Local variables ...
! array of the matrix elements of the bond operators, i.e. h for i=j and 2|J_ij| else
real(dp) :: M_ij( size(J_interaction_matrix, 1), size(J_interaction_matrix, 1) )
real(dp) :: M_ij_sub( 1:S%Nbasis, 1:S%Nsites )
! relative probabilities for choosing the first index of an Ising operator
real(dp) :: P_first( size(J_interaction_matrix, 1) )

integer, parameter :: n_opclass = 3 
real(dp) :: sum_all_diagmatrix_elements_peropclass(1:n_opclass)
real(dp) :: consts_added_per_opclass(1:n_opclass)
real(dp) :: prob_opclass(1:n_opclass)

real(dp) :: cc, ss, norm
integer :: ir, jr, k, n, subi

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
! number of lattice sites 
n = size(J_interaction_matrix, 1)
if( n /= S%Nsites ) then 
  print*, "ERROR: init_probtables(): n /= S%Nsites "
  stop
endif 
call assert( all(shape(J_interaction_matrix) == (/S%Nsites, S%Nsites/)) ) 


M_ij(:,:) = abs(J_interaction_matrix(:,:)) 
do ir=1,n
! Overwrite Ising self-interactions that might arise from Ewald summation.
! (They would only contribute a constant energy offset.)
  M_ij(ir,ir) = abs(hx)
enddo

! In our decomposition of the Hamiltonian there are 3 classes 
! of diagonal operators (`n_opclass`)
! (1) 2-leg CONSTANT or 4-leg ISING operators 
! (2) 6-leg TRIANGULAR_PLAQUETTE operators 
! (3) 2-leg LONGITUDINAL (i.e. `hz`) operators
probtable%n_opclass = n_opclass
allocate(probtable%prob_opclass(1:n_opclass))
probtable%prob_opclass(:) = 0.0_dp
allocate(probtable%sum_all_diagmatrix_elements_peropclass(1:n_opclass))
probtable%sum_all_diagmatrix_elements_peropclass(:) = 0.0_dp
allocate(probtable%consts_added_per_opclass(1:n_opclass))
probtable%consts_added_per_opclass(:) = 0.0_dp
allocate(probtable%idx_opclass(1:n_opclass))
probtable%idx_opclass(1) = TWOLEGCONST_OR_FOURLEG
probtable%idx_opclass(2) = TRIANGULAR_PLAQUETTE
probtable%idx_opclass(3) = LONGITUDINAL 
allocate(probtable%cumprob_opclass(1:n_opclass))
probtable%cumprob_opclass(:) = 0.0_dp
allocate(probtable%hz_matrix_element(1:S%Nsites))
probtable%hz_matrix_element(:) = 0.0_dp
allocate(probtable%hz_cumul(1:S%nsites))
probtable%hz_cumul(:) = 0.0_dp

! 1. The constants which have been added to the Hamiltonian artificially have to be
! subtracted from the energy in the end.
! 2. Sum over all matrix elements used in the transition probabilities 
! when deciding whether to insert or remove an operator

! (1) 2-leg CONSTANT or 4-leg ISING operators
! -------------------------------------------
cc = n*hx  ! constant operators added 
ss = n*hx  ! sum of all matrix elements 
do ir=1, n
do jr=1, ir-1 ! avoid double counting and count only ir.ne.jr
  cc = cc + M_ij(ir,jr)      ! abs(J_interaction_matrix(ir,jr)) 
  ss = ss + TWO*M_ij(ir,jr)  ! TWO*abs(J_interaction_matrix(ir,jr))
enddo
enddo
probtable%sum_all_diagmatrix_elements_peropclass(1) = ss 
probtable%consts_added_per_opclass(1) = cc

! (2) TRIANGULAR PLAQUETTE operators
! ----------------------------------
! see Ref. [1]
! matrix elements from triangular plaquettes 
! (of the Hamiltonian shifted by the constants)
if (trim(S%lattice_type) == "triangular") then 
    ! Triangular lattice 
    probtable%consts_added_per_opclass(2) = (3.0_dp/2.0_dp) * abs(J_1) * n_plaquettes
    probtable%sum_all_diagmatrix_elements_peropclass(2) = 2.0_dp * abs(J_1) * n_plaquettes
elseif (trim(S%lattice_type) == "kagome") then 
    ! Kagome lattice 
    probtable%consts_added_per_opclass(2) = 3.0_dp * abs(J_1) * n_plaquettes
    probtable%sum_all_diagmatrix_elements_peropclass(2) = 4.0_dp * abs(J_1) * n_plaquettes
elseif (trim(S%lattice_type) == "chain") then    
    ! No plaquettes at all. 
    probtable%consts_added_per_opclass(2) = 0.0_dp
    probtable%sum_all_diagmatrix_elements_peropclass(2) = 0.0_dp
else
    stop "init_probtables(): Unknown lattice type."
endif 
 
! (3) LONGITUDINAL field operators:
! ---------------------------------
ss = 0.0_dp
do ir = 1, S%Nsites 
  if ( abs(hz_fields(ir))>0 ) then 
    ss = ss + C_par_hyperparam + abs(hz_fields(ir))
  endif
enddo 
probtable%consts_added_per_opclass(3) = ss

! after the consts for all operator classes have been initialized: 
probtable%consts_added = sum(probtable%consts_added_per_opclass(:))

! For hz (longitudinal field) operators, the insertion probability depends on the 
! instantaneous spin configuration. 
! `update_probtables()` needs to be called whenever the 
! instantaneous spin configuration has changed.
call update_probtables(S=S, probtable=probtable, hz_fields=hz_fields, &
    C_par_hyperparam=C_par_hyperparam, spins=spins)

! Calculate the cumulative probability tables used in the diagonal update
! for the insertion of diagonal 2leg CONSTANT and 4leg ISING vertices. 
if( TRANSLAT_INV) then 

  allocate( probtable%P_cumulfirst(1) ) ! P_cumulfirst is not needed for translationally invariant systems 
  allocate( probtable%P_cumulsecond(1:S%Nbasis, 1:S%Nsites) )
  probtable%P_cumulfirst(1) = 1.0_dp / S%Nsites ! actually, not needed 
  probtable%P_cumulsecond(:,:) = 0.0_dp

  ! assuming that the interaction matrix is indeed translationally 
  ! invariant => IMPROVE
  do subi = 1, S%Nbasis
    M_ij_sub(subi, :) = M_ij(subi, :)
  enddo 

  ! The relative probability for choosing j as the second index 
  ! given i=1,...,S%Nbasis as the first index in a translationally 
  ! invariant system.
  do subi = 1, S%Nbasis
      do k = 1, S%Nsites
          do jr = 1, k
              probtable%P_cumulsecond(subi, k) = probtable%P_cumulsecond(subi, k) &
                  + abs(M_ij_sub(subi,jr))
          enddo
      enddo
  enddo

  ! normalize
  do subi = 1, S%Nbasis
      norm = sum(M_ij_sub(subi,1:S%Nsites))
      probtable%P_cumulsecond(subi, :) = probtable%P_cumulsecond(subi, :) / norm      
  enddo

else 
      allocate(probtable%P_cumulfirst(n))
      allocate(probtable%P_cumulsecond(n,n))

      P_first(:) = 0.0_dp
      probtable%P_cumulfirst(:) = 0.0_dp
      probtable%P_cumulsecond(:,:) = 0.0_dp

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

      ! normalize
      do ir = 1, n
          norm = 0.0_dp
          do jr = 1, n
              norm = norm + abs(M_ij(ir,jr))
          enddo
          do k = 1, n
              probtable%P_cumulsecond(ir,k) = probtable%P_cumulsecond(ir,k) / norm
          enddo
      enddo
endif

end subroutine     

subroutine update_probtables(S, probtable, hz_fields, C_par_hyperparam, spins)
  use lattice
  ! Purpose:
  ! --------
  ! Update those parts of the probability tables which depend on the instantaneous spin configuration.
  ! 
  ! When the matrix element of a diagonal operator thus the insertion probability depends
  ! on the instantaneous spin configuration, the corresponding part of the probability 
  ! table needs to be recomputed at each propagation step with a spin-flip operator.   
  implicit none 
  type(Struct), intent(in)         :: S
  type(t_ProbTable), intent(inout) :: probtable 
  real(dp)                         :: hz_fields(:)
  real(dp)                         :: C_par_hyperparam     
  integer, intent(in)              :: spins(:)
  
  ! ... Local variables ...
  integer :: ir, k

  probtable%hz_matrix_element(:) = 0.0_dp
  do ir = 1, S%Nsites
    ! if (abs(hz_fields(ir)) > 0.0_dp) then 
    !   ! The hz operators have two different matrix elements depending on the spin
    !   ! configuration. Add only the matrix elements of the hz operator for the given 
    !   ! instantaneous spin configuration.
    !   ! Only add constants for non-zero fields. 
    !   if (spins(ir) * hz_fields(ir) < 0.0_dp) then ! anti-aligned with the field 
    !     probtable%hz_matrix_element(ir) = C_par_hyperparam
    !   else if (spins(ir) * hz_fields(ir) > 0.0_dp) then ! aligned with the field
    !     probtable%hz_matrix_element(ir) = 2*abs(hz_fields(ir)) + C_par_hyperparam
    !   endif
    ! endif 
    ! REMOVE
    ! Include both possible matrix elements on a site => Then, actually, the probability tables needn't be recomputed !!!
    probtable%hz_matrix_element(ir) = 2*abs(hz_fields(ir)) + 2*C_par_hyperparam
    ! REMOVE 
  enddo
  probtable%sum_all_diagmatrix_elements_peropclass(3) = sum(probtable%hz_matrix_element(:))

  ! Cumulative probability table 
  do ir = 1, S%nsites
    probtable%hz_cumul(ir) = sum(probtable%hz_matrix_element(1:ir))
  enddo 
  probtable%hz_cumul(:) = probtable%hz_cumul(:) / sum(probtable%hz_matrix_element(:))

  
  probtable%sum_all_diagmatrix_elements = sum(probtable%sum_all_diagmatrix_elements_peropclass(:))
    
    ! Probability for inserting an operator from one of the
    ! `n_opclass` classes. 
    probtable%prob_opclass(:) = probtable%sum_all_diagmatrix_elements_peropclass(:) &
        / probtable%sum_all_diagmatrix_elements
    ! call assert( all(probtable%prob_opclass >= 0.0_dp), "neg. probs. in prob_opclass" )
     
    ! Cumulative probabilities for sampling
    do k = 1, probtable%n_opclass
        probtable%cumprob_opclass(k) = sum(probtable%prob_opclass(1:k))
    enddo
  
end subroutine 

subroutine extend_cutoff(opstring, config)
! Purpose:
! ---------
! Extend the cut-off of the fixed length operator string to
!   LL_new = n_exp + n_exp / 2. 
! 
  use SSE_configuration 
  use types 
  implicit none 
! Arguments:
! ---------_
  type(t_BondOperator), allocatable, intent(inout) :: opstring(:)
  type(t_Config), intent(inout) :: config 
  
! ... Local variables ...
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
      opstring_new(ip)%optype = IDENTITY
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