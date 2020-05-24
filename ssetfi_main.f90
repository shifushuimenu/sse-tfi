module ssetfi_main
    use types
    use SSE_configuration
    use diagonal_update 
    use lattice 
    use cluster_update
    use linked_list

    implicit none
    private 

    public one_MCS_plaquette
    public init_SSEconfig_hostart

    contains

subroutine one_MCS_plaquette(beta, Jij_sign, spins, opstring, &
                config, probtable, plaquettes, vertexlink, &
                leg_visited)
    ! *****************************************************
    ! Perform one Monte-Carlo step (MCS), consisting of 
    ! diagonal and off-diagonal update. 
    ! For the case of a plaquette update with triangular 
    ! plaquettes, we cycle over A-, B-, and C-update for 
    ! completing one MCS. 
    ! *****************************************************   

    real(dp), intent(in) :: beta
    integer, allocatable, intent(in) :: Jij_sign(:,:)
    integer, allocatable, intent(inout) :: spins(:)
    type(t_BondOperator), allocatable, intent(inout) :: opstring(:)    
    type(t_Config), intent(inout) :: config
    type(t_ProbTable), intent(in) :: probtable      
    type(t_Plaquette), allocatable, intent(in) :: plaquettes(:)
    integer, allocatable, intent(inout) :: vertexlink(:)
    logical, allocatable, intent(inout) :: leg_visited(:)  
        
    integer :: ut ! update type

    ! IMPROVE: no obscure number codes 
    do ut = 111, 113, 1
        ! loop over A-update (=111), B-update (=112) and C-update (=113)
        call diagonal_update_plaquette(beta=beta, &
            Jij_sign=Jij_sign, spins=spins, opstring=opstring, &
            config=config, probtable=probtable, plaquettes=plaquettes,&
            update_type=ut)
        call build_linkedlist_plaquette( &
            opstring=opstring, config=config, &
            vertexlink=vertexlink, leg_visited=leg_visited )        
        call quantum_cluster_update_plaquette( &
            spins=spins, opstring=opstring, vertexlink=vertexlink, &
            leg_visited=leg_visited, config=config )
    enddo

end subroutine 


subroutine init_SSEconfig_hostart( n_sites, LL, config, &
    opstring, spins, vertexlink, leg_visited )
    ! ***********************************************************
    ! Initialize the operator string with identities.
    ! Set the spin configuration at propagation step 1 randomly. 
    ! ***********************************************************

    integer, intent(in) :: n_sites
    integer, intent(in) :: LL
    type(t_Config), intent(out) :: config 
    type(t_BondOperator), allocatable, intent(out) :: opstring(:)
    integer, allocatable, intent(out) :: spins(:)
    integer, allocatable, intent(out) :: vertexlink(:)
    logical, allocatable, intent(out) :: leg_visited(:)

    integer :: ip, ir 
    real(dp) :: prob

    config%n_sites = n_sites
    config%n_exp = 0
    config%LL = LL
    config%n2leg = 0
    config%n4leg = 0 
    config%n6leg = 0
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL
    config%n_legs = 2*config%n2leg+4*config%n4leg+6*config%n6leg

    ! IMPROVE: For triangular lattice with PBC only:
    ! Assign n_plaquettes from the lattice structure S
    config%n_plaquettes=2*n_sites

    allocate(opstring(LL))
    do ip=1, LL
        opstring(ip)%i = 0; opstring(ip)%j = 0; opstring(ip)%k = 0
    enddo 
    allocate( spins(n_sites) )
    
    ! hot start 
    do ir = 1, n_sites
        call random_number(prob)
        if( prob < 0.5 ) then 
            spins(ir) = -1
        else
            spins(ir) = +1
        endif 
    enddo
    
    allocate(vertexlink( config%n_ghostlegs ))
    allocate(leg_visited( config%n_ghostlegs) )

end subroutine 

end module ssetfi_main

program ssetfi 
    use types
    use SSE_configuration
    use diagonal_update 
    use lattice 
    use cluster_update
    use linked_list
    use measurements
    use ssetfi_globals
    use ssetfi_main
    use MPI_parallel
    use util, only: init_RNG
    implicit none 

    ! ****************************************
    ! simulation parameters                  !
    ! ****************************************
    REAL(dp) :: J_1 = +1.0_dp    
    real(dp) :: hx = 0.60_dp   
    real(dp) :: beta = 5.0_dp    
    integer :: nx, ny, n_sites             

    integer :: nmeas_step = 10000
    integer :: ntherm_step = 10000 
    integer :: Nbin = 100
    ! ****************************************

    type(Phys) :: P0
    type(Struct) :: S

    integer ::  ir, jr, k, i
    integer :: it, im 
    integer :: obs
    integer :: ioerr 
    NAMELIST /SIMPARAMS/ J_1, hx, beta, nx, ny, n_sites, nmeas_step, ntherm_step, Nbin

#if defined (USE_MPI)
    include "mpif.h"
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_size, ierr)
#else
    MPI_rank = 0
    MPI_size = 1
#endif 
    write(chr_rank, "(i5.5)") MPI_rank

    OPEN( UNIT=5, FILE='simparams.in', ACTION="read", STATUS="old", IOSTAT=ioerr)
    IF( ioerr /=0 ) STOP "File simparams.in could not be opened."
    print*, "reading params"
    READ(5, NML=SIMPARAMS)

    ! TODO: Check input parameters ...

    ! seed random number generator with the system time (at the millisecond level)
    call init_RNG(MPI_rank) 

    hx = 0.05 + MPI_rank * 0.1
    ! beta = 1.0 + MPI_rank * 0.2

    if (nmeas_step < Nbin) then 
        print*, "Need nmeas_step >= Nbin."
        stop
    endif 

    if (mod(nmeas_step, Nbin) /= 0) then 
        print*, "`nmeas_step` should be an integer multiple of `Nbin`."
        stop
    endif 

    S%Nsites = n_sites

    call init_lattice_triangular(nx=nx, ny=ny, &
        neigh=neigh, sublattice=sublattice, plaquettes=plaquettes) 
    ! Specify interactions beyond nearest neighbours 
    ! Nearest neighbour interactions are already taken care 
    ! of by the plaquette operators. 
    allocate( J_interaction_matrix(n_sites,n_sites) )
    allocate( Jij_sign(n_sites, n_sites) )
    J_interaction_matrix(:,:) = 0.0_dp
    do ir = 1, n_sites
        do jr = 1, n_sites
            do k = 1, coord
                if (neigh(k, ir) == jr) then 
                  ! nearest neighbour interactions are already taken 
                  ! care of by the plaquette operators 
                  J_interaction_matrix(ir, jr) = -0.0_dp
                endif 
            enddo
        enddo         
    enddo   

    ! REMOVE
    open(100, file="Jmatrix.dat", status="unknown", action="write")
    do ir = 1, n_sites
        write(100, *) ( J_interaction_matrix(ir, jr), jr = 1, n_sites )
    enddo 
    close(100)
    open(101, file="sublattice.dat", status="unknown", action="write")
    do ir = 1, n_sites 
        write(101, *) ir, sublattice(ir)
    enddo 
    close(101)
    ! REMOVE

    ! Precompute the probability tables from which diagonal operators 
    ! will be sampled. 
    call init_probtables( J_interaction_matrix=J_interaction_matrix, &
        hx=hx, probtable=probtable, Jij_sign=Jij_sign, J_1=J_1, &
        n_plaquettes=size(plaquettes,dim=1), TRANSLAT_INV=.FALSE.)

    call init_SSEconfig_hostart( n_sites=nx*ny, LL=10, config=config, &
        opstring=opstring, spins=spins, vertexlink=vertexlink, leg_visited=leg_visited )

    do it = 1, ntherm_step
    
        call one_MCS_plaquette( beta=beta, Jij_sign=Jij_sign, spins=spins, &
            opstring=opstring, config=config, probtable=probtable, &
            plaquettes=plaquettes, vertexlink=vertexlink, &
            leg_visited=leg_visited )

        if ( (float(config%n_exp) / float(config%LL)) > 2.0_dp / 3.0_dp ) then 
            print*, "Extending cutoff LL_old=",  config%LL
            call extend_cutoff( opstring=opstring, config=config )
            print*, "LL_new=", config%LL
        endif 
    enddo

    call Phys_Init(P0, S, beta, Nbin)

    do im = 1, nmeas_step
        call one_MCS_plaquette( beta=beta, Jij_sign=Jij_sign, spins=spins, &
            opstring=opstring, config=config, probtable=probtable, &
            plaquettes=plaquettes, vertexlink=vertexlink, &
            leg_visited=leg_visited )
        call Phys_Measure(P0, config, spins, opstring, beta, probtable%consts_added)

        if (mod(im, nmeas_step / Nbin) == 0) then 
            call Phys_Avg(P0)
        endif 
    enddo

    call Phys_GetErr(P0)

    ! Output (clean up !)
    print*, hx, beta, &
            P0%meas(P0_ENERGY, P0%avg), P0%meas(P0_ENERGY, P0%err), &
            P0%meas(P0_MAGNETIZATION, P0%avg), P0%meas(P0_MAGNETIZATION, P0%err), &
            P0%meas(P0_COPARAM, P0%avg), P0%meas(P0_COPARAM, P0%err)
    open(500, file='averages'//chr_rank//'.dat', position='append', status='unknown')
    write(500, *) hx, beta, &
                !   P0%meas(P0_ENERGY, P0%avg), P0%meas(P0_ENERGY, P0%err), &
                !   P0%meas(P0_MAGNETIZATION, P0%avg), P0%meas(P0_MAGNETIZATION, P0%err), &
                !   P0%meas(P0_COPARAM, P0%avg), P0%meas(P0_COPARAM, P0%err)
                ( P0%meas(obs, P0%avg), P0%meas(obs, P0%err), obs = 1, P0%Nscalar_prop )
    close(500)

    if( MPI_rank == root_rank) then 
        open(700, file="output.txt", status="unknown", action="write")
        write(700, *) "The columns in the file averageXXXXX.dat have the following meaning:"
        write(700, *) "Simulation parameters:"
        write(700, *) "transverse field        : 1"
        write(700, *) "inverse temperatuer beta: 2"
        write(700, *) "Scalar observables: avg, err"
        do obs = 1, P0%Nscalar_prop
            write(700, *) P0_STR(obs), 2 +(obs-1)*2 + 1, 2 +(obs-1)*2 + 2
        enddo
        close(700)
    endif 


    call deallocate_globals

#if defined(USE_MPI)
    call MPI_FINALIZE(ierr)
#endif 

end program 
