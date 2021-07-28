! TODO:
!  - Do not allocate arrays of length LL, find a more appropriate way since LL can be very large.
!  - Longitudinal field operators.
!  - Code verification for a single spin in a longitudinal and transverse field. 

module local_update

 contains

subroutine local_offdiagonal_update(spins, opstring, config, &
                      hz_fields, C_par_hyperparam)

use SSE_configuration
use util, only: assert 
use types
implicit none
    
! Arguments:
! ==========
integer, intent(inout)               :: spins(:)
type(t_BondOperator), intent(inout)  :: opstring(:) 
type(t_Config), intent(in)           :: config
real(dp), intent(in)                 :: hz_fields(:)
real(dp), intent(in)                 :: C_par_hyperparam

! ... Local definitions and variables ...

! ipsite associates for each 2-leg vertex its position ip in the operator string
! with the linearly stored site ir on which it acts. Ising delimiters are defined
! to sit at both positions ir=i1 and ir=i2. This information is needed to rebuild
! the operator string frmo the sub-strings. 
type :: t_IpSite
  integer :: ip
  integer :: ir
end type t_IpSite

type (t_IpSite), allocatable :: ipsite(:)

type t_weights
  real(dp) :: weight_new 
  real(dp) :: weight_old
end type t_weights 

type(t_weights), allocatable :: opsubstring_weights(:,:)

integer, allocatable  :: opsubstring_site(:, :)	
! opsubstring_site(ir, :) contains for each site ir a sequence of 2-leg operators 
! that act on this site. Subsequences of unconstrained 2-leg vertices are
! separated by one or several Ising operators (which are represented 
! by a single "Ising delimiter").
! The notation for the elements of opsubstring_site(ir,:) is (see SSE_configuration.f90)
!	-1: Ising delimiter (Ising bond or triangular plaquette)
!	 22: constant operator
!	 21: spin-flip operator
! NOTE: The arrays opsubstring_site(ir,:) have by default length LL. Most of the entries are not used and therefore
! not initialized. Here it is assumed that they are initialized by default (i.e. when allocating) to a value
! different from -1, 1 or 2.

integer :: optype
integer, allocatable :: n_opsubstring_site(:) ! number of elements in opsubstring_site(ir,:) (2-leg vertices and Ising delimiters)
integer, allocatable :: counter(:) ! allows to jump over Ising operators in the operator 
				   ! sequences on each site when rebuilding the operator string
logical, allocatable :: ISING_LAST(:)   ! .TRUE. if last operator in the substring was an Ising delimiter. ISING_LAST(ir) is only set 
           ! to .FALSE. by detection of spin-flip or constant operators at site ir, but not by longitudinal field operators, which are simply ignored. 
logical, allocatable :: TWOLEG_LAST(:)  ! .TRUE. if last operator in the substring was a constant or spin-flip operator 
logical, allocatable :: LONGITUDINAL_LAST(:)   ! .TRUE. if last operator in the substring was a longitudinal operator  => combined weights 
logical, allocatable :: touched(:)  ! touched(ir) indicates whether the spin ir has any operator acting on it 
				    ! If not, it can be flipped with probability 1/2 after the local off-diagonal update in order to 
				    ! ensure ergodicity, which might be violated, especially for small transverse fields

logical :: WRAPPING ! WRAPPING=.TRUE. means that two operators that are to be replaced in the off-diagonal update
		    ! are connected through imaginary time

integer :: n_tot, n_test ! = n_2leg + n_Ising_delimiter (not equal to n_exp !)
integer :: i, ii, l, site, ip, subip, subip_next, ir, i1,i2
integer :: op1, op2, op1_new, op2_new
integer :: k1, k2, l1, l2
integer :: n_exp_new
integer :: ix, i3
integer :: sites(1:3) ! a plaquette has three sites, a bond operator only two sites  

integer :: alignment 
real(dp) :: weight_new, weight_old, prob
real(dp) :: eta

if( config%n2leg < 2 ) then 
    return 
endif 

allocate(opsubstring_site(1:config%n_sites, 1:5*config%LL))  ! IMPROVE
allocate(opsubstring_weights(0:config%n_sites, 0:5*config%LL))
opsubstring_site(:,:) = NOT_USED ! initialize to invalid values 

do ip=1,config%LL 
  do ir=1,config%n_sites 
    opsubstring_weights(ir, ip)%weight_new = -13 ! invalid value
    opsubstring_weights(ir, ip)%weight_old = -13 
  enddo 
enddo

allocate(n_opsubstring_site(config%n_sites)) ! Mustn't the length of opsubstring_site and ipsite be 2*LL ??????
allocate(counter(config%n_sites))
counter(:) = 0
allocate(ISING_LAST(config%n_sites))
allocate(TWOLEG_LAST(config%n_sites))
allocate(LONGITUDINAL_LAST(config%n_sites))
ISING_LAST(:) = .FALSE.  
TWOLEG_LAST(:) = .FALSE.
LONGITUDINAL_LAST(:) = .FALSE.
allocate(ipsite(5*config%LL))
allocate(touched(config%n_sites))
touched(:) = .FALSE.
  
n_tot = 0 ! n_tot is used also as a counter 

! ****************************************************************
! 1. Compile the operator substrings for each site and mark 
! whether the site has any operator acting on it or not.
! ****************************************************************
do ip = 1, config%LL
  ! extract sites where the operator acts
  i1 = opstring(ip)%i
  i2 = opstring(ip)%j
  optype = opstring(ip)%optype
  
  select case(optype)
      case(CONSTANT)
        touched(i1) = .TRUE.
        counter(i1) = counter(i1) + 1
        opsubstring_site(i1, counter(i1)) = CONSTANT
        ISING_LAST(i1) = .FALSE.
        LONGITUDINAL_LAST(i1) = .FALSE.
        n_tot = n_tot + 1
        ipsite(n_tot)%ip = ip
        ipsite(n_tot)%ir = i1
      case(SPIN_FLIP)
        touched(i1) = .TRUE.
        counter(i1) = counter(i1) + 1
        opsubstring_site(i1, counter(i1)) = SPIN_FLIP
        ISING_LAST(i1) = .FALSE.
        LONGITUDINAL_LAST(i1) = .FALSE.
        n_tot = n_tot + 1
        ipsite(n_tot)%ip = ip
        ipsite(n_tot)%ir = i1
      case(ISING_BOND)
        sites(1:2) = (/i1, i2/)
        do ix = 1,2
          touched(sites(ix)) = .TRUE.
          if (.not.ISING_LAST(sites(ix))) then 
              counter(sites(ix)) = counter(sites(ix)) + 1
              opsubstring_site(sites(ix), counter(sites(ix))) = DELIMITER
              ISING_LAST(sites(ix)) = .TRUE.
              n_tot = n_tot + 1
              ipsite(n_tot)%ip = ip
              ipsite(n_tot)%ir = sites(ix)
          endif 
        enddo       
      case(TRIANGULAR_PLAQUETTE)
        i3 = opstring(ip)%k
        sites(1:3) = (/abs(i1), abs(i2), abs(i3)/)
        do ix = 1,3
          touched(sites(ix)) = .TRUE.
          if (.not.ISING_LAST(sites(ix))) then 
              counter(sites(ix)) = counter(sites(ix)) + 1
              opsubstring_site(sites(ix), counter(sites(ix))) = DELIMITER
              ISING_LAST(sites(ix)) = .TRUE.
              n_tot = n_tot + 1
              ipsite(n_tot)%ip = ip
              ipsite(n_tot)%ir = sites(ix)
          endif 
        enddo 
      
      case(LONGITUDINAL)
        ! Do not reset ISING_LAST(i1).
        if ( .not.ISING_LAST(i1) ) then  ! excludes DELIMITER
          ! Store the weight changes due to longitudinal field operators between *unconstrained* 
          ! spin-flip operators. Longitudinal field operators on constrained segments do not 
          ! matter since no spin-flip will occur there, and so they are not recorded. The new/old
          ! weights of consecutive longitudinal field operators are combined and the combined weights 
          ! are stored. 

          if( .not.LONGITUDINAL_LAST(i1) ) then ! could have been CONSTANT, SPIN_FLIP
            ! reset weights 
            weight_new = 1.0_dp
            weight_old = 1.0_dp
          endif 
            
          ! For longitudinal field (`hz`) operators, the structure component
          ! opstring(ip)%k is used to store whether the hz operator 
          ! is aligned with the spin it is sitting on or not. opstring(ip)%k = +1(-1) means
          ! aligned (anti-aligned).
          alignment = opstring(ip)%k
          if (alignment > 0) then 
              ! spin aligned with the field 
              weight_new = weight_new * C_par_hyperparam
              weight_old = weight_old * (TWO * abs(hz_fields(i1)) + C_par_hyperparam)
          else
              weight_old = weight_old * C_par_hyperparam
              weight_new = weight_new * (TWO * abs(hz_fields(i1)) + C_par_hyperparam)
          endif 
          LONGITUDINAL_LAST(i1) = .TRUE.
          
          ! opsubstring_weights(:,:) only has valid values at positions where a LONGITUDINAL_WEIGHT 
          ! sits. Consecutive longitudinal weights are combined into one and their weights 
          ! are combined. 
          ! counter(i1) is the same as for the last non-LONGITUDINAL operator
          opsubstring_weights(i1, counter(i1))%weight_new = weight_new
          opsubstring_weights(i1, counter(i1))%weight_old = weight_old 

          ! if( (opsubstring_site(i1, counter(i1)) /= SPIN_FLIP) .and. &
          !     (opsubstring_site(i1, counter(i1)) /= CONSTANT) ) then 
          !   print*, "Error:  opsubstring_weights"
          !   stop
          ! endif 

          ! ! Convert (hz,aligned) into (hz,antialigned) vertex and vice versa.
          ! ! Note: (hz,aligned) and (hz,antialigned) are two different vertices with 
          ! !       different weights. 
          ! opstring(ip)%k = -opstring(ip)%k

        endif 

      case default
    !     print*, "jumping over identities"
        if ((i1.ne.0).or.(i2.ne.0)) then
          print*,"You probably forgot an if-case. i1=, i2= ", i1, i2
          stop
        endif
  end select 

enddo
  
n_opsubstring_site(1:config%n_sites) = counter(1:config%n_sites)
n_test = sum(counter(1:config%n_sites))
  
call assert( (n_test == n_tot), "n_tot /= n_test" )
  

do ir=1, config%n_sites
  ! ****************************************************************
! 2. Perform the local off-diagonal update for each operator substring
  ! ****************************************************************
  do i=1, n_opsubstring_site(ir)/2 ! divide number of update steps by two, since in each step potentially two operators are affected !!!!!! REMOVE
! Attempt off-diagonal update for a number of propagation steps proportional to the 
! length of the operator substring    
    if (n_opsubstring_site(ir) > 1) then

    WRAPPING = .FALSE.
    call random_number(eta)
    subip = int( eta*n_opsubstring_site(ir) ) + 1 ! subip is an element of {1,2,..., n_opsubstring_site(ir) }

! select two operators on consecutive positions at random
! NOTE: subip is only defined modulo the length of the operator substrings
    if (subip == n_opsubstring_site(ir)) then
    ! periodic boundary conditions in imaginary time
      subip_next = 1
      WRAPPING = .TRUE.
    else
      subip_next = subip + 1
    endif

    op1 = opsubstring_site(ir, subip)
    op2 = opsubstring_site(ir, subip_next)

    if (op1 == op2) then 

       if( opsubstring_weights(ir, subip)%weight_new /= -13 ) then 
          weight_new = opsubstring_weights(ir, subip)%weight_new
          weight_old = opsubstring_weights(ir, subip)%weight_old
          prob = weight_new / (weight_new + weight_old)
       else 
          prob = 0.5_dp
       endif 
       call random_number(eta)

! Two adjacent spin flip operators or two adjacent constants
! encountered => Perform update: (const)_p1,(const)_p2 <=> (spin flip)_p1,(spin flip)_p2
! and update the spin at site ir if WRAPPING = .TRUE.
! NOTE that - by construction - there should never be two delimiter operators on consecutive positions
! and likewise by construction there should not be two longitudinal weight operators on consecutive positions. 
      if (op1 == CONSTANT) then
        if (eta < prob) then 
          op1_new = SPIN_FLIP
          op2_new = SPIN_FLIP
          if (WRAPPING) then
            spins(ir) = - spins(ir)
          endif
        else
          op1_new = op1
          op2_new = op2 
        endif 
      elseif (op1 == SPIN_FLIP) then
        if (eta < prob) then 
          op1_new = CONSTANT
          op2_new = CONSTANT
          if (WRAPPING) then
            spins(ir) = - spins(ir)
          endif
        else
          op1_new = op1
          op2_new = op2 
        endif 
      else
!!!!!?????? Two Ising delimiters in a row across the boundary in imaginary time => don't exit ?????
        if ( .not.(subip == n_opsubstring_site(ir)) ) then 
          print*, "Error: Two strange operators in a row opsubstring_site(ir,:) ", "ir =", ir
          print*, "ir=", ir, " subip=", subip, " op(subip)=", op1, " op(subip+1)=", op2
          print*, "Here is the operator substring on site ir=", ir
          do ii=1, n_opsubstring_site(ir)
              print*, ii, opsubstring_site(ir, ii)
          enddo 
          ! print*, "And here is the full operator string filtered to show only Ising operators"
          ! do ip = 1, config%LL
          !     if (opstring(ip)%optype == ISING_BOND) then 
          !         i1 = opstring(ip)%i
          !         i2 = opstring(ip)%j
          !         print*, "i1, i2=", i1, i2
          !         if ((i1 == 0) .or. (i2 == 0)) then 
          !             stop 
          !         endif 
          !     endif 
          ! enddo
	        stop
        endif
! Do nothing, i.e. keep the adjacent Ising delimiters 
! across the boundary in imaginary time in their positions
        op1_new = op1
        op2_new = op2  
      endif

   elseif (op1 /= op2) then   
! If the two operators are not the same and there is no Ising delimiter
! among them, permute the operators.
! If one of the operators is a delimiter, no update procedure can be performed.

    if( opsubstring_weights(ir, subip)%weight_new /= -13 ) then 
      weight_new = opsubstring_weights(ir, subip)%weight_new
      weight_old = opsubstring_weights(ir, subip)%weight_old
      prob = weight_new / (weight_new + weight_old)
    else 
      prob = 0.5_dp
    endif 
    call random_number(eta)

! When swapping operators that are connected through the boundary in imaginary time,
! the initial spin configuration has to be updated. Only the (unconstrained) spin sequence between
! those two operators accross imaginary time is affected by this update.
      if ((op1 == CONSTANT).and.(op2 == SPIN_FLIP)) then
        if (eta < prob) then 
          op1_new = SPIN_FLIP
          op2_new = CONSTANT
          if (WRAPPING) then
            spins(ir) = - spins(ir)
          endif
        else
          op1_new = op1
          op2_new = op2 
        endif
      elseif ((op1 == SPIN_FLIP).and.(op2 == CONSTANT)) then
        if (eta < prob) then 
          op1_new = CONSTANT
          op2_new = SPIN_FLIP
          if (WRAPPING) then
            spins(ir) = - spins(ir)
          endif
        else
          op1_new = op1
          op2_new = op2 
        endif
      elseif ((op1 == DELIMITER).and.(op2 /= DELIMITER)) then
	! Do nothing
        op1_new = op1
        op2_new = op2
      elseif ((op2 == DELIMITER).and.(op1 /= DELIMITER)) then
	! Do nothing
        op1_new = op1
        op2_new = op2
      endif         
   endif

! re-assign the operators to the operator sequence at positions ip and ip+1 (at site ir)
    opsubstring_site(ir, subip) = op1_new
    opsubstring_site(ir, subip_next) = op2_new 

! Check
! Do not raise an error message if there are two consecutive Ising delimiters
! across the boundary in imaginary time.
  if ( .not.(subip == n_opsubstring_site(ir)) ) then
  if ((op1 == DELIMITER).and.(op2 == DELIMITER)) then
    print*, "Error: Two delimiter operators in a row"
    !stop
  endif
  endif

  endif ! only if n_opsubstring_site(ir) > 1

  enddo ! do n_opsubstring_site(ir) times ...
enddo ! Local off-diagonal update for each site ir

! ****************************************************************
! 3. Rebuild the operators string after the off-diagonal update:
! ****************************************************************
! Insert again all the 2-leg vertices from each array opsubstring_site(ir,:)
! into the total operator string thereby jumping over Ising delimiters
! encountered in each array opsubstring_site(ir,:).

 counter(1:config%n_sites) = 1

! n_tot is the total number of 2-leg operators (CONSTANT, SPIN_FLIP, DELIMITER) in all operator
! substrings {opsubstring_site(ir,:), ir=1,n_sites}. It does not correspond to the total number of operators 
! in the original operator string since consecutive Ising operators acting on one site are not recorded
! in the operator substrings.

do l = 1, n_tot
    ip = ipsite(l)%ip
    ir = ipsite(l)%ir
    if (opsubstring_site(ir, counter(ir)) == CONSTANT) then ! constant operator
      opstring(ip)%i = ir
      opstring(ip)%j = ir
      opstring(ip)%optype = CONSTANT
    elseif (opsubstring_site(ir, counter(ir)) == SPIN_FLIP ) then ! spin flip operator 
      opstring(ip)%i = ir
      opstring(ip)%j = 0
      opstring(ip)%optype = SPIN_FLIP
    endif
! ignore Ising delimiters and longitudinal weights at each site ir
    counter(ir) = counter(ir) + 1
enddo


! The off-diagonal update does not change the expansion order (=energy).
! It only changes the magnetization through the exchange of constants and spin flip operators,
! which (due to the periodicity in imaginary time) takes effect in the initial spin configuration.

! Check that the total expansion order has been maintained
n_exp_new = 0
do ip = 1, config%LL
  if ( .not.( opstring(ip)%optype .eq. IDENTITY ) ) then
    ! non-trivial operator encountered
    n_exp_new = n_exp_new + 1
  endif
enddo

if (n_exp_new /= config%n_exp) then
  print*, "subroutine local_offdiagonal_update: ERROR: n_exp_new = ", n_exp_new, "  n_exp = ", config%n_exp
  stop
endif


! Flip spins that do not have any operator acting on them with probability 1/2
do ir=1,config%n_sites
  if (.not.touched(ir)) then
    call random_number(eta)
    if (eta <= 0.5) spins(ir) = -spins(ir)
  endif
enddo


! Deallocate arrays
deallocate(opsubstring_site); deallocate(n_opsubstring_site); deallocate(counter)
deallocate(ISING_LAST); deallocate(ipsite); deallocate(touched)

end subroutine local_offdiagonal_update

end module local_update