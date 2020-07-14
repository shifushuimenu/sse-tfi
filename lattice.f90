! IMPROVE: Use the same subroutine for initializing the triangular
!   Bravais lattice of the kagome lattice as for the ordinary kagome lattice.
!   (use pbc() function)
! TODO: Initialize enlarged structure S for the triangular lattice. 
! Check that everywhere this version is used: 
!         subj = mod(j-1, S%Nbasis) + 1
! Add assert for translate_from_origin() and translate_to_origin()

module lattice 
use SSE_configuration, only: t_Plaquette
use types 
use util, only: PI
implicit none 

public t_Plaquette

public Struct 
public t_Kgrid
public init_lattice_triangular
public init_lattice_kagome
public momentum_grid_triangular_Bravais

! Information about lattice structure 
type Struct   
    character(len=30) :: lattice_type     
    integer :: coord ! coordination number
    integer :: rdim = 2  ! Number of spatial dimensions 
    integer :: Nx_bravais
    integer :: Ny_bravais 
    integer :: Nbravais_sites 
    integer :: Nbasis  ! number of sites in the motif
    integer :: Nsites  ! Nsites = Nbravais_sites * Nbasis

    real(dp), allocatable :: a1_p(:), a2_p(:) ! basis vectors of the real space Bravais lattice 
    real(dp), allocatable :: b1_p(:), b2_p(:) ! basis vectors of the reciprocal lattice
    real(dp), allocatable :: r_p(:,:) ! coordinates of the sites in the motif (=basis) with which the Bravais lattice is decorated
end type Struct 

! Momentum grid 
type t_Kgrid
  integer :: Nq      ! number of relevant momentum points
  integer, allocatable :: listk(:,:)
  real(dp), allocatable :: cosqr(:,:)
  real(dp), allocatable :: sinqr(:,:)
end type t_Kgrid 

contains 

! Comparison of arrays of plaquette structures   
!!! not tested yet !!!  
FUNCTION compare(p1, p2) RESULT(b)        
    TYPE(t_Plaquette), INTENT(IN) :: p1(:), p2(:)
    LOGICAL, allocatable :: b(:)
    integer :: i, n
    n = size(p1,1)
    allocate(b(n))
    i=0
    DO WHILE( i<n )
        i=i+1
        b(i) = (   (p1(i)%Asite == p2(i)%Asite) &
             .AND. (p1(i)%Bsite == p2(i)%Bsite) &
             .AND. (p1(i)%Csite == p2(i)%Csite) )         
    ENDDO
END FUNCTION      

   
SUBROUTINE init_lattice_triangular( &
            nx, ny, &
            S, neigh, sublattice, & 
            plaquettes )    
            
! *******************************************************************
! Purpose:
! --------
!    Initialize variables specifying a triangular lattice with
!    periodic boundary conditions in a way which is suitable 
!    for plaquette-based updates, see Ref. [1]
! 
!    With periodic boundary conditions the number of plaquettes
!    is equal to twice the number of sites. 
!
! Input: 
! ------
!    nx, ny: dimensions of the lattice 

! Output:
! -------
!    S: type(Struct), lattice structure object 
!    neigh(0:6, 1:nx*ny): lineraly stored indices of the 
!       nearest neighbours of a given site index
!    sublattice(1:nx*ny) \\in [1,2,3]: sublattice index of a 
!       linearly stored site index 
!    plaquettes(1:2*nx*ny): Array of plaquette structs.
!       Contains linearly stored indices of the `A-sites`, 
!       `B-sites`, and `C-sites` of each plaquette. 
!       By construction the `A-site`(`B-site`, `C-site`) 
!       belongs to sublattice 1 (2, 3). See Ref. [1] for a 
!       justification of why this is important. 
!
! NOTE: With periodic boundary conditions a triangular lattice
! has 2*n_sites plaquettes (including up and down-triangles).
! *******************************************************************          
  
 IMPLICIT NONE
   
 INTEGER, INTENT(IN)  :: nx, ny
 TYPE(Struct), INTENT(OUT) :: S
 INTEGER, ALLOCATABLE, INTENT(OUT) :: neigh(:,:)
 INTEGER, ALLOCATABLE, INTENT(OUT) :: sublattice(:)
 TYPE (t_Plaquette), ALLOCATABLE, INTENT(OUT) :: plaquettes(:)    
 
 ! missing: information about momentum grid of the Bravais lattice 

 INTEGER :: n
 integer :: ir

 ! Variable for plaquette-based cluster update
 integer :: irA, irB, irC
 integer :: plaq_idx
 integer :: plaq_type
 integer :: arr(3)     

 n = nx*ny   

 S%lattice_type = "triangular"
 S%rdim = 2
 S%coord = 6
 S%Nsites = n
 S%Nx_bravais = nx
 S%Ny_bravais = ny
 S%Nbravais_sites = n
 S%Nbasis = 1
  
 allocate(S%a1_p(1:S%rdim)); allocate(S%a2_p(1:S%rdim))
 allocate(S%b1_p(1:S%rdim)); allocate(S%b2_p(1:S%rdim))
 allocate(S%r_p(1:S%Nbasis, 1:S%rdim))

 ! Basis vectors of the real-space triangular Bravais lattice.
 S%a1_p(1) = 1.d0;  S%a1_p(2) = 0.d0
 S%a2_p(1) = 0.5d0; S%a2_p(2) = sqrt(3.d0) / 2.d0
  
 ! Reciprocal lattice of the Bravais lattice is also triangular.                    
 ! delta_k in direction of the unit vectors of the reciprocal lattice
 S%b1_p(1) = 1.d0*PI; S%b1_p(2) = -1.d0*PI / sqrt(3.d0)
 S%b2_p(1) = 0.0;     S%b2_p(2) = 2.d0*PI / sqrt(3.d0)

 ! Basis is just one site for triangular lattice 
 S%r_p(1,:) = 0.d0


 CALL neighbours_triangular( ny, ny, S, neigh ) 
 CALL sublattice_triangular( nx, ny, sublattice ) 
  
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
  ! ASSIGN NEAREST NEIGHBOUR BONDS TO PLAQUETTE OPERATORS 
  ! (for triangular lattice)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IF( .NOT.(ALLOCATED(plaquettes)) ) ALLOCATE(plaquettes( 2*n )) 
  IF (SIZE(plaquettes) /= 2*n) THEN 
      PRINT*, "init_lattice_triangular(): ERROR, inconsistent input"
      PRINT*, "Exiting ..." 
     STOP
  ENDIF 
  
  plaq_idx = 0
  DO ir = 1, n
     
     ! Sites belonging to the "upside" triangular plaquette who's 
     ! lower left corner is site `ir`.
     ! (see depiction in the subroutine  'neighbours_triangular')
     plaq_idx = plaq_idx + 1

     irA = neigh(0, ir)
     irB = neigh(1, ir)
     irC = neigh(2, ir)
                 
     ! Make sure the `A-site`of a plaquette is always assigned 
     ! the lattice site from the A sublattice and similarly 
     ! for the `B-site`and `C-site`. 
     SELECT CASE(sublattice(ir))
         CASE(1)
             plaq_type = 0
         CASE(2) 
             plaq_type = 1
         CASE(3)
             plaq_type = 2
     END SELECT
     ! can be simplified since plaq_type = sublattice(ir) - 1
     ! Note the minus sign for anticyclic shift.
     arr = CSHIFT(array=(/irA, irB, irC/), shift=-plaq_type)
     
     plaquettes(plaq_idx)%Asite = arr(1)
     plaquettes(plaq_idx)%Bsite = arr(2)
     plaquettes(plaq_idx)%Csite = arr(3)

     ! Sites belonging to the "downside" triangular plaquette who's 
     ! upper right corner is site `ir`.
     plaq_idx = plaq_idx + 1

     irA = neigh(0, ir)
     irB = neigh(5, ir)
     irC = neigh(4, ir)
     
     ! Note the minus sign for anticyclic shift.
     arr = CSHIFT(array=(/irA, irB, irC/), shift=-plaq_type)
     
     plaquettes(plaq_idx)%Asite = arr(1)
     plaquettes(plaq_idx)%Bsite = arr(2)
     plaquettes(plaq_idx)%Csite = arr(3)

  ENDDO
 
#ifdef DEBUG_COMPLETED
  ! Write a testfile with valid values for 
  ! parameters and corresponding data structures.
  ! See subroutine 'unit_test()' for details. 
  OPEN(UNIT=4, FILE="ut.test", ACTION="write")
  WRITE(4,*) "&PARAMS"
  WRITE(4,*) "nx =" , nx
  WRITE(4,*) "ny =" , ny
  WRITE(4,*) "/"
  WRITE(4,*) ""
  WRITE(4,*) "&DATASTRUCTURES"   
  WRITE(4,*) "neigh_exp = ", neigh
  WRITE(4,*) "sublattice_exp = ", sublattice
  WRITE(4,*) "plaquettes_exp = ", plaquettes 
  WRITE(4,*) "/"
  WRITE(4,*) ""
  CLOSE(4)
#endif   

END SUBROUTINE init_lattice_triangular    
   
          
SUBROUTINE cyclic_shift(A, B, C, n)

! *************************************************************
! Purpose:
! --------
! Cyclically re-assign the values of `A`,`B`, and `C` in place:
!     A B C  ----->   C A B
! Perform `n` cyclic shift operations. 
! `n` == 0 corresponds to no shift.   
! *************************************************************
 
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: A, B, C
  INTEGER, INTENT(IN) :: n
  
  INTEGER :: i, tmp
  
  IF ( n < 0 .OR. n > 2) THEN 
      PRINT*, "cyclic_shift(): WARNING: meaningless input n=", n
      PRINT*, "Exiting ..."
      STOP
  ENDIF 
  
  DO i = 1, n
      tmp = A
      A = C 
      C = B
      B = tmp
  ENDDO  
       
END SUBROUTINE


SUBROUTINE neighbours_triangular( nx, ny, S, neigh )

! *****************************************************************
! Purpose:
! --------
! Allocate and assign the nearest neighbours of each lattice site 
! for a triangular lattice of dimensions (nx, ny).
! Fix the lattice geometry with periodic boundary conditions,
! by fixing the nearest-neighbour (and next-nearest neighbour) matrix
!
!                 8
!           | 3       | 2
!      9      |     |       7
!              |  |
!      ---4----- 0 ------1--- 
!              |  |
!     10     |      |       12
!          | 5        | 6
!               11
! *****************************************************************
 
 IMPLICIT NONE 

 INTEGER, INTENT(IN) :: nx, ny
 TYPE(Struct), INTENT(IN) :: S
 INTEGER, ALLOCATABLE, INTENT(OUT) :: neigh(:,:)
 
 INTEGER :: n
 INTEGER :: ix, iy, ir 
 
 ! local automatic arrays for obtaining the nearest neighbour matrix 
 INTEGER :: x_forw( nx*ny )
 INTEGER :: x_back( nx*ny )
 INTEGER :: diag_rightup( nx*ny )
 INTEGER :: diag_leftdown( nx*ny )
 INTEGER :: diag_leftup( nx*ny )
 INTEGER :: diag_rightdown( nx*ny )   
   
 n = nx*ny
 
 IF( .NOT.(ALLOCATED(neigh)) ) ALLOCATE(neigh(0:2*S%coord, n))
 ! In case the array has been allocated incorrectly before ...
 IF( SIZE(neigh,1) /= (2*S%coord+1) .OR. SIZE(neigh,2) /= n ) THEN
     PRINT*, "neighbours_triangular(): ERROR, inconsistent input"
     PRINT*, "Exiting ..."      
     STOP
 ENDIF  
                                                             
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! PERIODIC BOUNDARY CONDITIONS
 ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ! bulk of the lattice
  DO ix = 2,nx-1
     DO iy = 2,ny-1        
        ir = ix + (iy-1)*nx
        x_forw(ir) = ir + 1 
        x_back(ir) = ir - 1
        diag_rightup(ir) = ir + nx
        diag_leftdown(ir) = ir - nx 
        diag_leftup(ir) = ir + nx - 1 
        diag_rightdown(ir) = ir - nx + 1 
     ENDDO
  ENDDO    
  ! first row
  DO ir=2,nx-1
    x_forw(ir) = ir+1; x_back(ir) = ir-1
     diag_rightup(ir) = ir + nx
     diag_leftup(ir) = ir + nx - 1
     diag_rightdown(ir) = (ny-1)*nx + ir + 1
     diag_leftdown(ir) = (ny-1)*nx + ir
  ENDDO
  ! last row 
  DO ix = 2,nx-1
    ir = (ny-1)*nx + ix
    x_forw(ir) = ir+1; x_back(ir) = ir-1
    diag_rightup(ir) = ix
    diag_leftup(ir) = ix - 1
    diag_rightdown(ir) = ir - nx + 1
    diag_leftdown(ir) = ir - nx
  ENDDO
  ! left boundary
  DO iy = 2, ny-1
    ir = (iy-1)*nx + 1
    x_forw(ir) = ir + 1; x_back(ir) = ir + nx - 1
    diag_rightup(ir) = ir + nx 
    diag_leftup(ir) =  ir + 2*nx - 1
    diag_rightdown(ir) = ir - nx + 1
    diag_leftdown(ir) = ir - nx 
  ENDDO
  ! right boundary
  DO iy = 2, ny-1
    ir = (iy-1)*nx + nx
    x_forw(ir) = ir - nx +1; x_back(ir) = ir -1 
    diag_rightup(ir) = ir + nx
    diag_leftup(ir) = ir + nx -1
    diag_rightdown(ir) = ir - 2*nx + 1
    diag_leftdown(ir) = ir - nx
  ENDDO  
  ! corners	
  ! ir=1
  x_forw(1) = 2; x_back(1) = nx 
  diag_rightup(1) = nx+1
  diag_leftup(1) = 2*nx
  diag_rightdown(1) = nx*(ny-1) +2
  diag_leftdown(1) = nx*(ny-1) + 1
 ! ir=nx  
  x_forw(nx) = 1; x_back(nx) = nx-1
  diag_rightup(nx) = 2*nx 
  diag_leftup(nx) =  2*nx-1
  diag_rightdown(nx) = nx*(ny-1) + 1
  diag_leftdown(nx) = nx*ny  
 ! 
  ir=nx*(ny-1) + 1
  x_forw(ir) = ir + 1 ; x_back(ir) = nx*ny
  diag_rightup(ir) = 1
  diag_leftup(ir) = nx
  diag_rightdown(ir) = ir - nx +1
  diag_leftdown(ir) = ir - nx  
 !
  ir = nx*ny
  x_forw(ir) = nx*(ny-1) + 1; x_back(ir) = ir - 1
  diag_rightup(ir) = nx
  diag_leftup(ir) = nx - 1
  diag_rightdown(ir) = nx*(ny-2) + 1
  diag_leftdown(ir) = ir - nx
    
  DO ir = 1, n
     neigh(0,ir) = ir
     neigh(1,ir) = x_forw(ir)
     neigh(2,ir) = diag_rightup(ir)
     neigh(3,ir) = diag_leftup(ir)
     neigh(4,ir) = x_back(ir)
     neigh(5,ir) = diag_leftdown(ir)
     neigh(6,ir) = diag_rightdown(ir)
     neigh(7,ir) = diag_rightup(x_forw(ir))
     neigh(8,ir) = diag_leftup(diag_rightup(ir))
     neigh(9,ir) = diag_leftup(x_back(ir))
     neigh(10,ir) = diag_leftdown(x_back(ir))
     neigh(11,ir) = diag_leftdown(diag_rightdown(ir))
     neigh(12,ir) = diag_rightdown(x_forw(ir))
  ENDDO 

END SUBROUTINE     


SUBROUTINE sublattice_triangular( nx, ny, sublattice )        

! *************************************************************
! Purpose:
! --------
! Assign to each linearly stored site its sublattice index.  
! 
! Input:
! ------
!    nx, ny: dimensions of the lattice 
! Output:
! -------
!    sublattice(1:nx*ny) \\in [1,2,3]: array assigning each 
!       lattice site its sublattice index. 
!       sublattice(:) is an assumed shape dummy argument which 
!       must have been allocated by the calling function. 
! *************************************************************
   
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nx, ny
   INTEGER, ALLOCATABLE, INTENT(OUT) :: sublattice(:)
   
   INTEGER :: xg, yg, ir, ir_row, sub

   IF( .NOT.(ALLOCATED(sublattice)) )  ALLOCATE(sublattice( nx*ny ))
   IF (nx*ny /= size(sublattice, 1)) THEN  
       PRINT*, "sublattice_triangular(): ERROR, inconsistent input"
       STOP
   ENDIF 

DO xg = 1, nx
 DO yg = 1, ny
 
   ir = xg + (yg-1)*nx
   ir_row = mod(ir, nx)
   
   if( mod(yg,3).eq.1 ) then                                         
     if( mod(ir_row,3).eq.1 ) then                                   
       sub = 1
     elseif( mod(ir_row,3).eq.2 ) then
       sub = 2
     else
       sub = 3
     endif
   elseif( mod(yg,3).eq.2 ) then
     if( mod(ir_row,3).eq.1 ) then
       sub = 3
     elseif( mod(ir_row,3).eq.2 ) then
       sub = 1
     else
       sub= 2
     endif
   elseif( mod(yg,3).eq.0 ) then
     if( mod(ir_row,3).eq.1 ) then
       sub = 2
     elseif( mod(ir_row,3).eq.2 ) then
       sub = 3
     else
       sub = 1
     endif      
   endif
   
   sublattice(ir) = sub

 ENDDO
ENDDO    

END SUBROUTINE 
       

SUBROUTINE unit_test_triangular(testfile)
   IMPLICIT NONE 
   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: testfile
   
   CHARACTER(LEN=101) :: default_testfile
   INTEGER :: ioerr
   
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Copy-paste the interface of the routine to be tested 
   ! and create variables for 'expected' results, which
   ! are read from a testfile. 
   INTEGER :: nx, ny
   INTEGER, ALLOCATABLE :: neigh(:,:), neigh_exp(:,:)
   INTEGER, ALLOCATABLE :: sublattice(:), sublattice_exp(:)
   TYPE (t_Plaquette), ALLOCATABLE :: plaquettes(:), plaquettes_exp(:)      
   TYPE (Struct) :: S
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   NAMELIST /PARAMS/ nx, ny
   NAMELIST /DATASTRUCTURES/ neigh_exp, sublattice_exp, plaquettes_exp

               
   IF( PRESENT(testfile) ) THEN
       default_testfile = testfile
   ELSE
       default_testfile="../unit_tests/triangular_lattice.test"
   ENDIF 
         
   print*, TRIM(default_testfile)
   OPEN( UNIT=5, FILE=TRIM(default_testfile), ACTION="read", STATUS="old", IOSTAT=ioerr )
   IF( ioerr /= 0) STOP "File <testfile> could not be opened."
  
   print*, "reading params"
   READ(5, NML=PARAMS)

   ! generated values           ; ! expected values read from file
   ALLOCATE(neigh( 0:6, nx*ny )); ALLOCATE(neigh_exp( 0:6, nx*ny ))
   ALLOCATE(sublattice( nx*ny )); ALLOCATE(sublattice_exp( nx*ny ))
   ALLOCATE(plaquettes( nx*ny )); ALLOCATE(plaquettes_exp( nx*ny ))
           
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
   ! routine to be checked: Generate data structures and compare with 
   ! expected values read from file. 
   CALL init_lattice_triangular( nx, ny, &
         S, neigh, sublattice, plaquettes )
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
   
   ! expected values read from file 
   print*, "reading data structures"
   READ(5, NML=DATASTRUCTURES)
   CLOSE(5)
   
   ! compare
   IF( ALL(neigh == neigh_exp) .AND. &
       ALL(sublattice == sublattice_exp) .AND. &
       ALL(compare(plaquettes, plaquettes_exp)) ) PRINT*, "Test passed. Done."
          
END SUBROUTINE unit_test_triangular        


! IMPROVE: make this a pure function 
integer function pbc(nr,l)
! ==============================================
! periodic boundary conditions in one dimension 
! ==============================================
! nr: site index [0, ..., l-1]
! l: length of the chain 
! pbc returns a site index [0, ... ,l-1]
  implicit none
  integer, intent(in) :: nr
  integer, intent(in) :: l
  pbc = nr
  if (nr.gt.l-1) pbc = nr - l
  if (nr.lt.0) pbc = nr + l
  if ((nr.ge.2*l).or.(nr.lt.-l)) then 
     print*, "Error in function pbc(nr,l)"
     print*, "invalid arguments: nr=", nr, "l=", l
     print*, "Exiting..."
     stop
  endif 
end function pbc    


subroutine init_lattice_kagome( &          
    nx, ny, &
    S, neigh, sublattice, & 
    plaquettes )  
! Kagome lattice with nearest neighbour interactions and periodic boundary conditions. 
! Nx_brav and Ny_brav are the dimensions of the triangular Bravais lattice
! so that N_unit_cells = Nx_brav * Ny_brav is the number of unit cells. 

! Parameters:
! -----------
implicit none

integer, intent(in) :: nx, ny
type(Struct), intent(out) :: S 
integer, allocatable, intent(out) :: neigh(:,:)
integer, allocatable, intent(out) :: sublattice(:)
type(t_Plaquette), allocatable, intent(out) :: plaquettes(:)

! ... Local variables ...
integer :: ir,jr,ix,iy,imj,sub
integer, allocatable :: neigh_tmp(:,:) ! helper array, which is zero-indexed and only used inside this function 
    
integer :: irB, jrB
integer :: imj_ix, imj_iy,  subi, subj, i, j, nk, iq
integer, parameter :: rdim = 2! number of spatial dimensions 
double precision :: rimj_p(1:rdim)
integer, allocatable :: leftB(:), rightB(:), left_upB(:), right_upB(:), left_downB(:), right_downB(:)

double precision ::  qvec(1:rdim), ri(1:rdim)

integer :: plaq_idx, ir1, ir2, ir3

double precision, allocatable :: rvec(:,:)
double precision, allocatable :: local_quant_axis(:,:)

! CLEAN UP !
! Bravais lattice reciprocal lattice for kagome lattice 
integer, parameter :: coordB=6 ! coordination number of the Bravais lattice 
integer, allocatable :: listB(:,:) ! second index of listB(1:2, 0:N_unit_cells-1) starts at "0" !
integer, allocatable :: invlistB(:,:) ! returns an index from 0 to N_unit_cells-1
integer, allocatable :: neighB(:,:)
integer, allocatable :: lattB_imj(:,:) ! map distance between site i and j to distance from origin.
                                       ! (This data structure refers to the underlying Bravais lattice).
integer, allocatable :: imjBdeg(:) ! How often does the distance imj occur ?                                           
double precision, allocatable :: x_pos(:), y_pos(:)
! CLEAN UP !


S%lattice_type = "kagome"
S%coord = 4                 ! coordination number of the kagome lattice
S%Nbasis = 3                ! number of sites in the basis of the kagome lattice (=motif)
S%Nbravais_sites = nx*ny    ! number of points in the (triangular) Bravais lattice 
S%Nx_bravais = nx
S%Ny_bravais = ny
S%Nsites = S%Nbravais_sites * S%Nbasis   !! total number of lattice sites

allocate(neigh(S%coord, 1:S%Nsites))
allocate(neigh_tmp(S%coord,0:S%Nsites-1))
allocate(sublattice(1:S%Nsites))

! ===========================================================
! TRIANGULAR BRAVAIS LATTICE     
! build a triangular Bravais lattice (B) with periodic boundary conditions 
! arrays and variables labelled by "B" are just helper variables 
allocate(lattB_imj(1:S%Nbravais_sites, 1:S%Nbravais_sites))
allocate(imjBdeg(1:S%Nbravais_sites))       
allocate(neighB(coordB, 0:S%Nbravais_sites-1))
allocate(listB(1:rdim, 0:S%Nbravais_sites-1))
allocate(invlistB(0:S%Nx_bravais-1, 0:S%Ny_bravais-1))
allocate(leftB(0:S%Nbravais_sites-1)); allocate(rightB(0:S%Nbravais_sites-1))
allocate(left_upB(0:S%Nbravais_sites-1)); allocate(right_upB(0:S%Nbravais_sites-1)); 
allocate(left_downB(0:S%Nbravais_sites-1)); allocate(right_downB(0:S%Nbravais_sites-1))
do ix=0,S%Nx_bravais-1
  do iy=0,S%Ny_bravais-1
    irB = ix + iy*S%Nx_bravais
    listB(1,irB) = ix
    listB(2,irB) = iy
    invlistB(ix,iy) = irB   
  enddo
enddo
! nearest neighbours on the Bravais lattice 
do irB = 0,S%Nbravais_sites-1 !!!!! ZERO-INDEXING NECESSARY 
    ix = listB(1, irB)
    iy = listB(2, irB) 
    leftB(irB)       = invlistB( pbc(ix-1, S%Nx_bravais), iy )
    rightB(irB)      = invlistB( pbc(ix+1, S%Nx_bravais), iy )
    left_upB(irB)    = invlistB( pbc(ix-1, S%Nx_bravais), pbc(iy+1, S%Ny_bravais) )
    right_upB(irB)   = invlistB( pbc(ix, S%Nx_bravais), pbc(iy+1, S%Ny_bravais) )
    left_downB(irB)  = invlistB( pbc(ix, S%Nx_bravais), pbc(iy-1, S%Ny_bravais) )
    right_downB(irB) = invlistB( pbc(ix+1, S%Nx_bravais), pbc(iy-1, S%Ny_bravais) )        
    neighB(4,irB) = leftB(irB)       
    neighB(3,irB) = rightB(irB)      
    neighB(5,irB) = left_upB(irB)    
    neighB(6,irB) = right_upB(irB)   
    neighB(1,irB) = left_downB(irB)  
    neighB(2,irB) = right_downB(irB) 
enddo    

! ! lattB_imj
! do jrB = 0, S%Nbravais_sites-1
!     do irB = 0, S%Nbravais_sites-1
!         imj_ix = pbc( listB(1,irB) - listB(1,jrB), S%Nx_bravais )
!         imj_iy = pbc( listB(2,irB) - listB(2,jrB), S%Ny_bravais )
!         lattB_imj(irB+1,jrB+1) = invlistB(imj_ix, imj_iy) + 1
!     enddo
! enddo
! ! the degeneracy of LattB_imj
!  imjBdeg(:) = 0
!  do jrB = 0, S%Nbravais_sites-1
!     do irB = 0, S%Nbravais_sites-1
!        imj = lattB_imj(irB+1,jrB+1)
!        imjBdeg(imj) = imjBdeg(imj) + 1
!     end do
!  end do
 
 allocate(S%a1_p(1:S%rdim)); allocate(S%a2_p(1:S%rdim))
 allocate(S%b1_p(1:S%rdim)); allocate(S%b2_p(1:S%rdim))
 allocate(S%r_p(1:S%Nbasis, 1:S%rdim))

 ! Basis vectors of the real-space triangular Bravais lattice.
 ! The length unit is the lattice constant "a" of the kagome lattice
 ! such that the basis vectors of the Bravais lattice have length 2*a.
 S%a1_p(1) = 1.d0;  S%a1_p(2) = 0.d0
 S%a2_p(1) = 0.5d0; S%a2_p(2) = sqrt(3.d0) / 2.d0
  
 ! Reciprocal lattice of the Bravais lattice is also triangular.                    
 ! delta_k in direction of the unit vectors of the reciprocal lattice
 S%b1_p(1) = 2.d0*PI; S%b1_p(2) = -2.d0*PI / sqrt(3.d0)
 S%b2_p(1) = 0.0;     S%b2_p(2) = 4.d0*PI / sqrt(3.d0)

 ! Basis vectors of the motif
 S%r_p(1,1) = 0.d0
 S%r_p(1,2) = 0.d0
 S%r_p(2,:) = 0.5d0 * S%a1_p(:)
 S%r_p(3,:) = 0.5d0 * S%a2_p(:)
 
 ! real space position of the sites 
 allocate(rvec(rdim,S%Nsites))
 allocate(local_quant_axis(rdim,S%Nsites))
 do ir = 1, S%Nsites
     sub = mod(ir-1,3) + 1
     ! map site index to index in the Bravais lattice 
     irB = (ir-1)/3
     ! real space position of the site with index ir
     rvec(:,ir) = listB(1,irB)*S%a1_p + listB(2,irB)*S%a2_p + S%r_p(sub,:)     
     
     ! vector of the local quantization axis => only used for visualization of kagome ice loops
     select case (sub)
         case(1)
             local_quant_axis(1,ir) = cos(PI/6.d0)
             local_quant_axis(2,ir) = sin(PI/6.d0)
         case(2)
             local_quant_axis(1,ir) = -cos(PI/6.d0)
             local_quant_axis(2,ir) = sin(PI/6.d0)
         case(3)
             local_quant_axis(1,ir) = 0.d0
             local_quant_axis(2,ir) = 1.d0
     end select
 enddo

! =======================================
! Nearest neighbours of each kagome site.
! =======================================    
! CONVENTION: The numbering of neighbours of ir is such that neighbours 1 and 2 belong
! to the same up-plaquette as ir. Neighbours (2 and 3) belong to one sublattice, neighbours (1 and 4)
! belong to the other sublattice, ir itself belongs to the third sublattice. 
! CONVENTION NECESSARY FOR GEOMETRIC CLUSTER UPDATE: 
! Take ir as the inversion center. Then the numbering of neighbours is such that 
!  inv_neigh(ir,neigh(ir)) = 5 - neigh(ir). 
! In general: inv_neigh(ir, neigh(ir)) = (coord+1) - neigh(ir)
do ir=0, S%Nsites-1
   irB = int(ir/3) ! point on the Bravais lattice 
   sub = mod(ir,3) + 1 ! sublattice(ir) can have values 1,2,3
   sublattice(ir+1) = sub !!!! INDEX STARTS AT ONE
   if (sub.eq.1) then
      neigh_tmp(2,ir) = ir + 2 ! right-up
      neigh_tmp(4,ir) = leftB(irB)*3 + 1 ! left
      neigh_tmp(1,ir) = ir + 1 ! right
      neigh_tmp(3,ir) = left_downB(irB)*3 + 2 ! left-down
   elseif (sub.eq.2) then 
      neigh_tmp(2,ir) = ir + 1 ! left-up
      neigh_tmp(1,ir) = ir - 1 ! left
      neigh_tmp(4,ir) = rightB(irB)*3 ! right 
      neigh_tmp(3,ir) = right_downB(irB)*3 + 2 ! right-down 
   elseif (sub.eq.3) then 
      neigh_tmp(3,ir) = left_upB(irB)*3 + 1  ! left-up
      neigh_tmp(4,ir) = right_upB(irB)*3     ! right-up 
      neigh_tmp(1,ir) = ir - 2  ! left-down
      neigh_tmp(2,ir) = ir - 1  ! right-down
   endif 
enddo

! shift all site indices by 1 to conform with Fortran indexing style 
do ir = 0, S%Nsites-1
   neigh(:, ir+1) = neigh_tmp(:, ir) + 1
enddo
deallocate(neigh_tmp)

deallocate(leftB); deallocate(rightB); deallocate(left_upB)
deallocate(right_upB); deallocate(left_downB); deallocate(right_downB)   

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
! ASSIGN NEAREST NEIGHBOUR BONDS TO PLAQUETTE OPERATORS 
! (for kagome lattice)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if( .not.(allocated(plaquettes)) ) allocate( plaquettes( 2*S%NBravais_sites ) ) 
if (size(plaquettes) /= 2*S%Nbravais_sites) then  
    print*, "init_lattice_kagome(): ERROR, inconsistent input:"
    print*, "size(plaquettes) /= 2*S%Nbravais_sites => Exiting ..." 
    stop
endif

plaq_idx = 1
do irB = 0, S%NBravais_sites-1
    ! Loop over sites on A-sublattice. 
    ir1 = irB*3 + 1
    ! "UP"-plaquettes have an odd plaquette index 
    ir2 = neigh(1, ir1)
    ir3 = neigh(2, ir1)
    plaquettes(plaq_idx) = t_Plaquette(ir1, ir2, ir3)
    plaq_idx = plaq_idx + 1

    ! "DOWN"-plaquettes have an even plaquette index 
    ir2 = neigh(4, ir1)
    ir3 = neigh(3, ir1)
    plaquettes(plaq_idx) = t_Plaquette(ir1, ir2, ir3)
    plaq_idx = plaq_idx + 1    
enddo 

! REMOVE
! Check the subroutine translate_kagome 
! print*, "Check translate_kagome"
! do ir1 = 1, 27
!   do ir2 = 1, 27 
!     print*, ir1, ir2, " => ", translate_kagome(S, ir1 ,ir2)
!   enddo 
! enddo 
! stop

end subroutine init_lattice_kagome


subroutine make_translat_invar( S, J_interaction_matrix, J_translat_invar )
  implicit none 
! Purpose:
! ========
!   Given a translationally invariant interaction matrix of shape (S%Nsites, S%Nsites),
!   generate a more compact interaction matrix of dimension 
!         (0:S%Nbravais_sites-1, 1:S%Nbasis*S%Nbasis),
!   which takes into account translational invariance. 
!
!   This routine should be generic for any lattice structure S.
!   If the input interaction matrix is not translationally invariant, an 
!   error is thrown. 
! 
! Arguments:
! ==========  
  type(Struct) :: S
  real(dp), allocatable, intent(in) :: J_interaction_matrix(:,:)
  real(dp), allocatable, intent(out) :: J_translat_invar(:,:)

  ! ... Local variables ...
  integer :: ir, jr 
  integer :: r(2)
  real(dp), parameter :: INIT = -1000000

  ! ... Executable ...
  allocate( J_translat_invar(0:S%Nbravais_sites-1, 1:S%Nbasis*S%Nbasis) )
  J_translat_invar(:,:) = INIT 

  do jr = 1, S%Nsites 
    do ir = 1, S%Nsites ! self-interaction is possible as a result of Ewald summation 
      r = translate_to_origin(S, ir, jr)
      if( J_translat_invar( r(1), r(2) ) == INIT) then 
        J_translat_invar( r(1), r(2) ) = J_interaction_matrix( ir, jr )
      ! Check for translational invariance 
      elseif( abs(J_translat_invar( r(1), r(2) ) - J_interaction_matrix( ir, jr )) &
             > epsilon(1.0_4) ) then 
        print*, "Error: make_translat_invar_kagome(): interaction matrix is not"
        print*, "       translationally invariant."
        print*, J_translat_invar( r(1), r(2) ), J_interaction_matrix( ir, jr )
        stop
      endif 
    enddo 
  enddo

  if( any(J_translat_invar == INIT) ) then 
    print*, "J_translat_invar(:,:) not fully initialized"
  endif 
 
end subroutine make_translat_invar

! IMPROVE: make this a pure function 
function translate_to_origin(S, i, j) result(r)
  implicit none 
! Purpose:
! ========
!   With `i` and `j` the linearly stored sites on the kagome 
!   lattice, return the linearly stored site index `d` of the 
!   distance vector between the two sites `i` and `j`
!   and the necessary additional sublattice information. 
!   ("translation to the origin"). Thus, this function provides 
!   the two indices (distance, sublattice info) into the translationally 
!   invariant interaction matrix. 
! 
!   site labelling for e.g. kagome lattice:
!
!        3           6           9
! 
!     1 --- 2 --- 4 --- 5 --- 7 --- 8 
! (origin)
! 
!   This subroutine should be generic for any lattice structure S. 
!
! Usage: Indexing into a translationally invariant interaction matrix. 
! ====== to get the sign of the interactions. 
!
!      Jsign_transinvar = make_translat_invar_kagome( J_sign )
!      sign = Jsign_transinvar( translate_kagome(i,j) )
! 
! Arguments:
! ==========
  type(Struct), intent(in) ::  S
  integer, intent(in) :: i, j
! r(1): Bravais distance (ZERO-INDEXED !)
! r(2): sublattice "difference"
  integer :: r(2)

! ... Local variables ...
  integer :: iB, jB
  integer :: ixB, iyB, jxB, jyB

! ... Executable ...
! Linearly stored Bravais site indices start at zero. 
  iB = (i-1) / S%Nbasis
  jB = (j-1) / S%Nbasis
! Bravais coordinates between sites i and j   
  ixB = mod(iB, S%Nx_bravais)
  iyB = iB / S%Nx_bravais
  jxB = mod(jB, S%Nx_bravais)
  jyB = jB / S%Nx_bravais

! Bravais distance ( linearly stored index )
  r(1) = pbc(jxB-ixB, S%Nx_bravais) + pbc(jyB-iyB, S%Nx_bravais) * S%Nx_Bravais
! sublattice "distance" ( for a Bravais lattice, r(2) = 1 )
  r(2) = sublattice_distance(S, i, j)

end function translate_to_origin   

! IMPROVE: make this a pure function 
function translate_from_origin( S, i, j) result(k)
! Input: index2, which is relative to the origin 
! Return index2 such that it is relative to index 1. 
  type(Struct), intent(in) :: S
  integer, intent(in) :: i, j 
  integer :: k

! ... Local variables ...
  integer :: iB, jB
  integer :: ixB, iyB, jxB, jyB  
  integer :: d, subj

! Linearly stored Bravais site indices start at zero. 
  iB = (i-1) / S%Nbasis
  jB = (j-1) / S%Nbasis
! Bravais coordinates of sites i and j 
  ixB = mod(iB, S%Nx_bravais)
  iyB = iB / S%Nx_bravais
  jxB = mod(jB, S%Nx_bravais)
  jyB = jB / S%Nx_bravais  

  d = pbc(jxB+ixB, S%Nx_bravais) + pbc(jyB+iyB, S%Nx_bravais) * S%Nx_Bravais  
  subj = mod(j-1, S%Nbasis) + 1

  ! linearly stored index after translation from `i`
  k = d * S%Nbasis + subj

end function translate_from_origin


pure function sublattice_distance(S, i, j) result(d)
  type(Struct), intent(in) :: S 
  integer, intent(in) :: i, j
  integer :: d
  integer :: subi, subj 

  subi = mod(i-1, S%NBasis) + 1
  subj = mod(j-1, S%Nbasis) + 1
  d = subi + S%Nbasis * (subj -1)

end function sublattice_distance

subroutine momentum_grid_triangular_Bravais(S, Kgrid)
  ! Purpose:
  ! --------
  ! Initialize data structure needed for Fourier transformation 
  ! on a triangular Bravais lattice (i.e. triangular or kagome lattice).
  ! 
  ! Arguments:
  ! ----------
  implicit none 

  type(Struct), intent(in) :: S
  type(t_Kgrid), intent(out) :: Kgrid
  
  ! ... Local variables ...
  integer :: Nq
  integer :: ik, i, j
  integer :: iq, ir
  real(dp) :: qvec(1:S%rdim), ri(1:S%rdim)
  
  ! IMPROVE !!!!!:  This part can be a source of inconsistency 
  ! since listB(:,:) is also defined locally in init_lattice_kagome()
  integer :: ix, iy, irB, sub
  integer :: listB(2,0:S%Nbravais_sites-1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Here, we could select to compute just some high-symmetry momentum points.
  Nq = S%Nbravais_sites
  Kgrid%Nq = Nq

  if (.not.allocated(Kgrid%listk) ) allocate( Kgrid%listk(1:S%rdim, 1:Nq) )
  if (.not.allocated(Kgrid%cosqr) ) allocate( Kgrid%cosqr(1:S%Nsites, 1:Nq) )
  if (.not.allocated(Kgrid%sinqr) ) allocate( Kgrid%sinqr(1:S%Nsites, 1:Nq) )

  ! listk: grid coordinates of points in reciprocal space 

  ik = 0 
  do j = 0, S%Ny_bravais-1 
     do i = 0, S%Nx_bravais-1 
         ik = ik+1
         Kgrid%listk(1,ik) = i-S%Nx_bravais/2
         Kgrid%listk(2,ik) = j-S%Ny_bravais/2
     end do
  end do
  if( ik /= Nq ) then
      stop " Error: nk /= Nq"
  end if


  ! IMPROVE: listB(:,:) appears also as local variable in the subroutine 
  ! init_lattice_kagome. THIS CAN BE A SOURCE OF INCONSISTENCY !
  do ix=0,S%Nx_bravais-1
    do iy=0,S%Ny_bravais-1
      irB = ix + iy*S%Nx_bravais
      listB(1,irB) = ix
      listB(2,irB) = iy
    enddo
  enddo

  ! Matrix elements of Fourier transform (Basis = Motif)
  do iq = 1, Nq
    qvec = Kgrid%listk(1,iq)*S%b1_p + Kgrid%listk(2,iq)*S%b2_p
    do irB = 0, S%Nbravais_sites-1 ! sum over all sites: Bravais lattice ... 
      do sub=1,S%Nbasis,1   ! ... and basis
        ri = listB(1,irB)*S%a1_p + listB(2,irB)*S%a2_p + S%r_p(sub,:)
        ir = irB*S%Nbasis + sub  ! IMPROVE: This is only valid for my number inf scheme of the kagome lattice ! MAKE THIS MORE STABLE !
        Kgrid%sinqr(ir, iq) = sin( qvec(1)*ri(1) + qvec(2)*ri(2) )
        Kgrid%cosqr(ir, iq) = cos( qvec(1)*ri(1) + qvec(2)*ri(2) )
      enddo 
    enddo
  enddo
  
end subroutine momentum_grid_triangular_Bravais

end module lattice 

