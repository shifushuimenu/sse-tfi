PROGRAM driver
    USE lattice
    IMPLICIT NONE 
    
    INTEGER, ALLOCATABLE :: neigh(:,:)
    INTEGER, ALLOCATABLE :: sublattice(:)
    TYPE (t_Plaquette), ALLOCATABLE :: plaquettes(:)
        
    integer :: nx=3, ny=3
    
    !allocate(neigh(0:6, nx*ny))
    !allocate(sublattice(nx*ny))
    !allocate(plaquettes(nx*ny))
        
    CALL init_lattice_triangular( &
        nx=nx, ny=ny,&
        neigh=neigh, sublattice=sublattice,&
        plaquettes=plaquettes )

    print*, "testing..."
    CALL unit_test('../unit_tests/triangular_lattice.test')
            
END PROGRAM 
