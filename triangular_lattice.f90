! TODO:
! ~~~~~
! 1. Overload comparison operator for plaquette structures: 
! 
!      ! Operator function for comparison of plaquette structures
!      INTERFACE OPERATOR (==)
!          FUNCTION comparison(p1, p2) RESULT(b)        
!              TYPE(tPlaquette), INTENT(IN) :: p1, p2
!              LOGICAL :: b           
!      END INTERFACE 
! 
! FUNCTION compare(p1, p2) RESULT(b)        
!      TYPE(tPlaquette), INTENT(IN) :: p1, p2
!      LOGICAL :: b
!      b = (      (p1%Asite == p2%Asite) &
!           .AND. (p1%Bsite == p2%Bsite) &
!           .AND. (p1%Csite == p2%Csite) )         
! END FUNCTION   
!    
! 3. Object-oriented design 
! 4. Bravais lattice as parent class

MODULE lattice
    USE SSE_configuration, only: t_Plaquette
    USE types
    IMPLICIT NONE
    PRIVATE

    PUBLIC t_Plaquette
    PUBLIC init_lattice_triangular
    PUBLIC unit_test
    PUBLIC Struct
        
    ! Information about lattice structure 
    TYPE Struct       
        integer :: coord ! coordination number
        integer :: Nbravais_sites 
        integer :: Nbasis
        integer :: Nsites ! Nsites = Nbravais_sites * Nbasis
    END TYPE Struct 

    CONTAINS

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
!    for plaquette-based updates, see Ref. [1]-
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

  S%coord = 6
  S%Nsites = n
  S%Nbravais_sites = n
  S%Nbasis = 1
   
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
       PRINT*, "cyclic_shift(): WARNING: meaninless input n=", n
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
        

SUBROUTINE unit_test(testfile)
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
           
END SUBROUTINE unit_test        


END MODULE                   
