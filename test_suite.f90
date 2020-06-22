! TODO:
!   - make sure no data structures contain unitialized values
!     (even if these values are never needed)
!     => opstring(ip) for two-leg and four-leg operators has no 
!        entry for %PRIVILEGED_LEG_IS_MAJORITY_LEG
!   - unit tests should fully cover all branches of the code

module test_suite
    use SSE_configuration
    implicit none 

    ! global variables 
    type(t_Config) :: config
    type(t_BondOperator), allocatable :: opstring(:)
    integer, allocatable :: spins(:)

    integer, allocatable :: vertexlink(:)
    logical, allocatable  :: leg_visited(:)

    logical, allocatable :: visited_ip(:)

    contains 

    subroutine init_test_SSEconfig

        ! prepare a test SSE configuration
        ! (see page 56 of my PhD thesis for a graphical depiction
        !  of the specific configuration)         
        config%n_sites = 4
        config%n_exp = 8
        config%LL = config%n_exp + 1 
        config%n2leg = 4
        config%n4leg = 2 
        config%n6leg = 2
        config%n_ghostlegs = MAX_GHOSTLEGS*config%n_exp
        config%n_legs = 2*config%n2leg+4*config%n4leg+6*config%n6leg
        config%n_plaquettes = 1

        allocate( opstring(config%LL) )
        allocate( spins(config%n_sites) )

        allocate( visited_ip(config%LL) )
        visited_ip(:) = .true.

        opstring(1)%i=1; opstring(1)%j=2
        opstring(2)%i=4; opstring(2)%j=0
        opstring(3)%i=1; opstring(3)%j=1
        opstring(4)%i=2; opstring(4)%j=4
        opstring(5)%i=1; opstring(5)%j=1
        opstring(6)%i=4; opstring(6)%j=0
        ! site indices of triangular plaquette operators are negative 
        opstring(7)%i=-2; opstring(7)%j=-3; opstring(7)%k=-4
            opstring(7)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .FALSE.
        opstring(8)%i=-1; opstring(8)%j=-2; opstring(8)%k=-3
            opstring(8)%PRIVILEGED_LEG_IS_MAJORITY_LEG = .TRUE.
        opstring(9)%i=0; opstring(9)%j=0


        spins = (/-1,-1, 1, 1/)

        allocate(vertexlink( config%n_ghostlegs ))
        allocate(leg_visited( config%n_ghostlegs) )
        vertexlink(:) = 0
        leg_visited(:) = .FALSE.

    end subroutine 

    subroutine test_cluster_update
        use linked_list
        use cluster_update

        call init_test_SSEconfig
        call build_linkedlist_plaquette(&
          opstring, config, vertexlink, leg_visited )
        call quantum_cluster_update_plaquette( &
            spins, opstring, vertexlink, leg_visited, config )

    end subroutine 

    subroutine test_diagonal_update 
        use diagonal_update 
        use lattice         

    end subroutine 

end module     

program unit_test
    use test_suite
    use test_helper
    implicit none
    ! call test_cluster_update
    call init_test_SSEconfig
    call output_SSE_config(config, opstring, spins, visited_ip, 'SSEconfig.dat')
end program 