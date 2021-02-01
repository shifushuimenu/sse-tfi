! IMPROVE: do not allocate arrays of length LL, find a more appropriate way since LL can be very large.

module local_update

 contains

subroutine local_offdiagonal_update(spins, opstring, config)

use SSE_configuration
use util, only: assert 
use types, only: dp
implicit none
    
! Arguments:
! ==========
integer, intent(inout)               :: spins(:)
type(t_BondOperator), intent(inout)  :: opstring(:) 
type(t_Config), intent(in)           :: config

! ... Local definitions and variables ...

! ipsite associates for each 2-leg vertex its position ip in the operator string
! with the linearly stored site ir on which it acts. Ising delimiters are defined
! to sit at both positions ir=i1 and ir=i2.
type :: tIpSite
  integer :: ip
  integer :: ir
end type tIpSite

real(dp) :: eta

integer, allocatable  :: opstring_site(:, :)	
! opstring_site(ir, :) contains for each site ir a sequence of 2-leg operators 
! that act on this site. Subsequences of unconstrained 2-leg vertices are
! separated by one or several Ising operators (which are represented 
! by a single "Ising delimiter").
! The notation for the elements of opstring_site(ir,:) is
!	-1: Ising delimiter (Ising bond or triangular plaquette)
!	 1: constant operator
!	 2: spin-flip operator
! NOTE: The arrays opstring_site(ir,:) have by default length LL. Most of the entries are not used and therefore
! not initialized. Here it is assumed that they are initialized by default (i.e. when allocating) to a value
! different from -1, 1 or 2.
integer :: DELIMITER, CONST_OP, FLIP_OP, NOT_USED
parameter (DELIMITER = -1)
parameter (CONST_OP = 1)
parameter (FLIP_OP = 2)
parameter (NOT_USED = 0)

integer :: optype
integer, allocatable :: n_opstring_site(:) ! number of elements in opstring_site(ir,:) (2-leg vertices and Ising delimiters)
type (tIpSite), allocatable :: ipsite(:)
integer, allocatable :: counter(:) ! allows to jump over Ising operators in the operator 
				   ! sequences on each site when rebuilding the operator string
logical, allocatable :: ISING_LAST(:) ! .TRUE. if last operator in the substring was an Ising delimiter

logical, allocatable :: touched(:)  ! touched(ir) indicates whether the spin ir has any operator acting on it 
				    ! If not, it can be flipped with probability 1/2 after the local off-diagonal update in order to 
				    ! ensure ergodicity, which might be violated, especially for small transverse fields

logical :: WRAPPING ! WRAPPING=.TRUE. means that two operators that are to be replaced in the off-diagonal update
		    ! are connected through imaginary time

integer :: n_tot, n_test ! = n_2leg + n_Ising_delimiter (not equal to n_exp !)
integer :: i, ii, l, site, ip, ip_site, ip_site_next, ir, i1,i2
integer :: op1, op2, op1_new, op2_new
integer :: k1, k2, l1, l2
integer :: n_exp_new
integer :: ix, i3

integer :: sites(1:3) ! a plaquette has three sites, a bond operator only two sites  

if( config%n2leg < 2 ) then 
    return 
endif 

allocate(opstring_site(config%n_sites,5*config%LL))
opstring_site(:,:) = NOT_USED ! initialize to invalid values 
allocate(n_opstring_site(config%n_sites)) ! Mustn't the length of opstring_site and ipsite be 2*LL ??????
allocate(counter(config%n_sites))
 counter(:) = 0
allocate(ISING_LAST(config%n_sites))
ISING_LAST(:) = .FALSE.  
allocate(ipsite(5*config%LL))
allocate(touched(config%n_sites))
touched(:) = .FALSE.
  
n_tot = 0

! 1. Compile the operator substrings for each site and mark whether the site
! has any operator acting on it or not
do ip = 1, config%LL
  ! extract sites where the operator acts
  i1 = opstring(ip)%i
  i2 = opstring(ip)%j
  optype = opstring(ip)%optype
  
  select case(optype)
      case(CONSTANT)
        touched(i1) = .TRUE.
        counter(i1) = counter(i1) + 1
        opstring_site(i1, counter(i1)) = CONST_OP
        ISING_LAST(i1) = .FALSE.
        n_tot = n_tot + 1
    ! Improve: confusingly n_tot is used as a counter here
        ipsite(n_tot)%ip = ip
        ipsite(n_tot)%ir = i1
      case(SPIN_FLIP)
        touched(i1) = .TRUE.
        counter(i1) = counter(i1) + 1
        opstring_site(i1, counter(i1)) = FLIP_OP
        ISING_LAST(i1) = .FALSE.
        n_tot = n_tot + 1
        ipsite(n_tot)%ip = ip
        ipsite(n_tot)%ir = i1
      case(ISING_BOND)
        sites(1:2) = (/i1, i2/)
        do ix = 1,2
          touched(sites(ix)) = .TRUE.
          if (.not.ISING_LAST(sites(ix))) then 
              counter(sites(ix)) = counter(sites(ix)) + 1
              opstring_site(sites(ix), counter(sites(ix))) = DELIMITER
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
              opstring_site(sites(ix), counter(sites(ix))) = DELIMITER
              ISING_LAST(sites(ix)) = .TRUE.
              n_tot = n_tot + 1
              ipsite(n_tot)%ip = ip
              ipsite(n_tot)%ir = sites(ix)
          endif 
        enddo 
      
      case default
    !     print*, "jumping over identities"
        if ((i1.ne.0).or.(i2.ne.0)) then
          print*,"You probably forgot an if-case. i1=, i2= ", i1, i2
          stop
        endif
  end select 

enddo
  
n_opstring_site(1:config%n_sites) = counter(1:config%n_sites)
n_test = sum(counter(1:config%n_sites))
  
call assert( (n_test == n_tot), "n_tot /= n_test" )
  

do ir=1, config%n_sites
! 2. Perform the local off-diagonal update for each operator substring
  do i=1, n_opstring_site(ir)/2 ! divide number of update steps by two, since in each step potentially two operators are affected !!!!!! REMOVE
! Attempt off-diagonal update for a number of propagation steps proportional to the 
! length of the operator substring    
    if (n_opstring_site(ir) > 1) then

    WRAPPING = .FALSE.
    call random_number(eta)
    ip_site = int( eta*n_opstring_site(ir) ) + 1 ! ip_site is an element of {1,2,..., n_opstring_site(ir) }

! select two operators on consecutive positions at random
! NOTE: ip_site is only defined modulo the length of the operator substrings
    if (ip_site == n_opstring_site(ir)) then
    ! periodic boundary conditions in imaginary time
      ip_site_next = 1
      WRAPPING = .TRUE.
    else
      ip_site_next = ip_site + 1
    endif

    op1 = opstring_site(ir, ip_site)
    op2 = opstring_site(ir, ip_site_next)
    
    if (op1 == op2) then 
! Two adjacent spin flip operators or two adjacent constants
! encountered => Perform update: (const)_p1,(const)_p2 <=> (spin flip)_p1,(spin flip)_p2
! and update the spin at site ir if WRAPPING = .TRUE.
! NOTE that - by construction - there should never be two delimiter operators on consecutive positions
      if (op1 == CONST_OP) then
	op1_new = FLIP_OP
	op2_new = FLIP_OP
	if (WRAPPING) then
	  spins(ir) = - spins(ir)
	endif
      elseif (op1 == FLIP_OP) then
	op1_new = CONST_OP
	op2_new = CONST_OP
	if (WRAPPING) then
	  spins(ir) = - spins(ir)
	endif
      else
!!!!!?????? Two Ising delimiters in a row across the boundary in imaginary time => don't exit ?????
        if ( .not.(ip_site == n_opstring_site(ir))) then 
	  print*, "Error: Two strange operators in a row opstring_site(ir,:) ", "ir =", ir
	  print*, "ir=", ir, " ip_site=", ip_site, " op(ip_site)=", op1, " op(ip_site+1)=", op2
	  print*, "Here is the operator substring on site ir=", ir
          do ii=1, n_opstring_site(ir)
              print*, ii, opstring_site(ir, ii)
          enddo 
!           print*, "And here is the full operator string filtered to show only Ising operators"
!           do ip = 1, config%LL
!               if (opstring(ip)%optype == ISING_BOND) then 
!                   i1 = opstring(ip)%i
!                   i2 = opstring(ip)%j
!                   print*, "i1, i2=", i1, i2
!                   if ((i1 == 0) .or. (i2 == 0)) then 
!                       stop 
!                   endif 
!               endif 
!           enddo
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

! When swapping operators that are connected through the boundary in imaginary time,
! the initial spin configuration has to be updated. Only the (unconstrained) spin sequence between
! those two operators accross imaginary time is affected by this update.
      if ((op1 == CONST_OP).and.(op2 == FLIP_OP)) then
	op1_new = FLIP_OP
	op2_new = CONST_OP
	if (WRAPPING) then
	  spins(ir) = - spins(ir)
	endif
      elseif ((op1 == FLIP_OP).and.(op2 == CONST_OP)) then
	op1_new = CONST_OP
	op2_new = FLIP_OP
	if (WRAPPING) then
	  spins(ir) = - spins(ir)
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
    opstring_site(ir, ip_site) = op1_new
    opstring_site(ir, ip_site_next) = op2_new 

! Check
! Do not raise an error message if there are two consecutive Ising delimiters
! across the boundary in imaginary time.
  if ( .not.(ip_site == n_opstring_site(ir)) ) then
  if ((op1 == DELIMITER).and.(op2 == DELIMITER)) then
    print*, "Error: Two delimiter operators in a row"
    stop
  endif
  endif

  endif ! only if n_opstring_site(ir) > 1

  enddo ! do n_opstring_site(ir) times ...
enddo ! Local off-diagonal update for each site ir

! 3. Rebuild the operators string after the off-diagonal update:
! Insert again all the 2-leg vertices from each array opstring_site(ir,:)
! into the total operator string thereby jumping over Ising delimiters
! encountered in each array opstring_site(ir,:).

 counter(1:config%n_sites) = 1

! n_tot is the total number of 2-leg operators (CONST_OP, FLIP_OP, DELIMITER) in all operator
! substrings {opstring_site(ir,:), ir=1,n_sites}. It does not correspond to the total number of operators 
! in the original operator string since consecutive Ising operators acting on one site are not recorded
! in the operator substrings.

do l = 1, n_tot
    ip = ipsite(l)%ip
    ir = ipsite(l)%ir
    if (opstring_site(ir, counter(ir)) == CONST_OP) then ! constant operator
      opstring(ip)%i = ir
      opstring(ip)%j = ir
    elseif (opstring_site(ir, counter(ir)) == FLIP_OP ) then ! spin flip operator 
      opstring(ip)%i = ir
      opstring(ip)%j = 0
    endif
! ignore Ising delimiters at each site ir
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
  print*, "subroutine local_offdiagonal_update: Caramba! n_exp_new = ", n_exp_new, "  n_exp = ", config%n_exp
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
deallocate(opstring_site); deallocate(n_opstring_site); deallocate(counter)
deallocate(ISING_LAST); deallocate(ipsite); deallocate(touched)

end subroutine local_offdiagonal_update

end module local_update