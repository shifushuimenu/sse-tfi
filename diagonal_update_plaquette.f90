! TODO:
!   - Implement translation invariance of probability tables 
!     and array indicating whether a bond is FM or AFM.
!   - if ( Jij_sign(r(1), r(2)).gt.0 ) then !AFM (line 142) 
!     => write such that it works also for non-translationally invariant systems


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
  integer, allocatable  :: idx_opclass(:)      ! integer indices for operator classes (used in the diagonal update)
  real(dp), allocatable :: hz_insert_aligned(:)  ! ( 2*abs(hz(i)) + C_par ) / ( 2*abs(hz(i)) + 2*C_par )
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
integer, allocatable, intent(in) :: hz_fields_sign(:)  ! sign of the longitudinal fields: FM (<0) or AFM(>0).
integer, intent(in) :: spins(:)
type(t_BondOperator), intent(inout) :: opstring(:)
type(t_Config), intent(inout) :: config
type(t_ProbTable), intent(in) :: probtable
type(t_Plaquette), intent(in) :: plaquettes(:)
! if update_type == A_UPDATE => "A-site update
integer, intent(in) :: update_type   
logical, intent(in) :: TRANSLAT_INVAR

! ... Local variables ...
integer :: ip, i1, i2, index1, index2, index2_from_origin
integer :: r(2)
integer :: sub1
integer :: optype

integer :: plaq_idx

real(dp) :: P_plus, P_minus, P_add, P_remove
real(dp) :: prob, eta
integer :: alignment 

! propagated spin configuration at a given propagation step
integer, allocatable :: spins2(:) 

! for plaquette-based cluster update 
integer :: ir_A, ir_B, ir_C

! "instantaneous" spin configuration
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
  
!identity encountered
if( (i1 == 0).and.(i2 == 0) ) then

    P_plus = P_add    
    call random_number(prob)
        
  if( prob <= P_plus ) then
  ! try to insert an operator  

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Choice among different classes of diagonal operators 
   ! First, decide whether to insert 
   !     - a triangular plaquette operator
   !     - or a diagonal 2-leg CONSTANT or 4-leg ISING vertex
   !     - or a diagonal 2-leg LONGITUDINAL operator 
   ! IMPROVE: maybe make 2-leg CONSTANT a separate class ? 
   
   ! choose with heat bath probability 
   call random_number(prob)
   optype = choose_diagopclass(prob, probtable)
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
! There is no risk of the search never ending as a constant can always be inserted (for h unequal 0 !).

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
    else !index1.eq.index2 => constant operator to be inserted
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
#endif          
     config%n_exp = config%n_exp + 1; config%n6leg = config%n6leg + 1 
     ! update P_add and  P_remove
     P_remove = float(( config%LL - config%n_exp + 1)) &
      / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
     P_add = beta*probtable%sum_all_diagmatrix_elements &
      / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
  endif 
  
  case( LONGITUDINAL )
     ! for homogeneous longitudinal fields
     call random_number(eta)
     i1 = int(eta*S%Nsites) + 1
     ! MISSING: sample positions for inhomogeneous fields
     ! based on hz_fields(:) 

     ! Choose between an aligned or antialigned 
     ! `hz` vertex according to the strength and orientation of the field,
     ! but independently of the spin alignment. 
     call random_number(eta)
     if ( eta <= probtable%hz_insert_aligned(i1) ) then
        alignment = +1
     else
        alignment = -1
     endif 
     !  if( spins2(i1) == hz_fields_sign(i1) ) then
     !     prob_insert = probtable%hz_insert_aligned(i1)
     !  else ! spin at the propagation step is not aligned with the longitudinal field 
     !     prob_insert = 1.0_dp - probtable%hz_insert_aligned(i1)
     !  endif      

     ! (Of course, it is possible that no hz vertex is inserted.)     
     if( spins2(i1) == alignment * hz_fields_sign(i1) ) then 
        opstring(ip)%optype = LONGITUDINAL 
        opstring(ip)%i = i1
        opstring(ip)%j = -i1 ! actually, not necessary to set this entry 
        opstring(ip)%k = spins2(i1)   ! The spin state is needed in the cluster update to accumulate the exchange fields 

        config%n_exp = config%n_exp + 1; config%n2leg_hz = config%n2leg_hz + 1 
        ! update P_add and  P_remove
        P_remove = float(( config%LL - config%n_exp + 1)) &
            / ( float(config%LL- config%n_exp + 1) + beta*probtable%sum_all_diagmatrix_elements )
        P_add = beta*probtable%sum_all_diagmatrix_elements &
            / ( float(config%LL-config%n_exp) + beta*probtable%sum_all_diagmatrix_elements )
     else
        ! Don't insert any operator. 
     endif 

  case default
      print*, "Diagonal update: Trying to insert unknown diagonal operator type"
      print*, "optype = ", optype
      print*, "Exiting ..."
      stop 

  end select choose_insert_optype

  endif !if(prob <= P_plus)
  
endif !identity encountered

  ! Diagonal operator encountered: 
  ! Try replacement:  Ising plaquette / Ising bond / const (hx) / longitudinal (hz) => identity
  ! Expansion order changes as n_exp -> n_exp - 1
  if ((i1.ne.0).and.(i2.ne.0)) then

    P_minus = P_remove
    call random_number(prob)
     
    if (prob <= P_minus) then
#ifdef DEBUG_DIAGONAL_UPDATE
      print*, "remove a diagonal operator"
#endif      
      optype = opstring(ip)%optype
           
      config%n_exp = config%n_exp - 1
      ! update P_add and  P_remove
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
          stop "diagonal_update(): ERROR: trying to remove unknown operator type"
      end select 
      
      opstring(ip)%i = 0
      opstring(ip)%j = 0
      opstring(ip)%k = 0 ! to be sure nothing bad happens
      opstring(ip)%optype = IDENTITY  

    endif
  endif

  ! Propagate spins as spin-flip operators are encountered.
  if ((i1.ne.0).and.(i2.eq.0)) then
    spins2(i1) = -spins2(i1)
  endif

enddo !do ip=1,LL

! Update important variables that have not been explicitly 
! updated
    config%n_legs = 2*config%n2leg+2*config%n2leg_hz+4*config%n4leg+6*config%n6leg
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL

deallocate(spins2)

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
  hz_fields, C_par_hyperparam )
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
!       replaced by a vector for a Bravais lattice (or by a matrix 
!        of shape (S%Nbasis, S%Nsites) for a lattice with a unit cell).
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
sum_all_diagmatrix_elements_peropclass(1) = ss 
consts_added_per_opclass(1) = cc

! (2) TRIANGULAR PLAQUETTE operators
! ----------------------------------
! see Ref. [1]
! matrix elements from triangular plaquettes 
! (of the Hamiltonian shifted by the constants)
if (trim(S%lattice_type) == "triangular") then 
    ! Triangular lattice 
    consts_added_per_opclass(2) = (3.0_dp/2.0_dp) * abs(J_1) * n_plaquettes
    sum_all_diagmatrix_elements_peropclass(2) = 2.0_dp * abs(J_1) * n_plaquettes
elseif (trim(S%lattice_type) == "kagome") then 
    ! Kagome lattice 
    consts_added_per_opclass(2) = 3.0_dp * abs(J_1) * n_plaquettes
    sum_all_diagmatrix_elements_peropclass(2) = 4.0_dp * abs(J_1) * n_plaquettes
else
    stop "init_probtables(): Unknown lattice type."
endif 
 
! (3) LONGITUDINAL field operators 
! --------------------------------
cc = 0
ss = 0
do ir = 1, n
  cc = cc + abs(hz_fields(ir)) + C_par_hyperparam
  ! `hz` vertex on an up-spin and `hz` vertex on a down-spin are regarded as 
  ! two different types of vertices. One has the matrix element `C_par_hyperparam`
  ! and the other has the matrix element 2*abs(hz_fields(i)) + `C_par_hyperparam`. 
  ! This is why `C_par_hyperparam` appears twice in the sum of all matrix elements. 
  ss = ss + 2*abs(hz_fields(ir)) + 2*C_par_hyperparam
enddo
sum_all_diagmatrix_elements_peropclass(3) = ss 
consts_added_per_opclass(3) = cc


! Probability for inserting an operator from one of the
! `n_opclass` classes. 
prob_opclass(:) = sum_all_diagmatrix_elements_peropclass(:) / sum(sum_all_diagmatrix_elements_peropclass(:))
call assert( all(prob_opclass >= 0.0_dp), "neg. probs. in prob_opclass" )

! Cumulative probabilities for sampling
allocate(probtable%cumprob_opclass(1:n_opclass))
do k = 1, n_opclass
    probtable%cumprob_opclass(k) = sum(prob_opclass(1:k))
enddo

allocate(probtable%idx_opclass(1:n_opclass))
probtable%idx_opclass(1) = TWOLEGCONST_OR_FOURLEG
probtable%idx_opclass(2) = TRIANGULAR_PLAQUETTE
probtable%idx_opclass(3) = LONGITUDINAL 

probtable%sum_all_diagmatrix_elements = sum(sum_all_diagmatrix_elements_peropclass(:))
probtable%consts_added = sum(consts_added_per_opclass(:))
probtable%n_opclass = n_opclass

! Probabilitiy for inserting an hz operator on an aligned spin state,
! (having already decided to insert and hz operator either on an aligned or 
! on an anti-aligned spin state).
! Conversely: probtable%hz_insert_antialigned = 1.0 - probtable%hz_insert_aligned
allocate(probtable%hz_insert_aligned(1:S%Nsites))
probtable%hz_insert_aligned(:) = (2*abs(hz_fields(:)) + C_par_hyperparam) / (2*abs(hz_fields(:)) + 2*C_par_hyperparam)


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