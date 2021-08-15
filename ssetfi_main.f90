module ssetfi_main
    use types
    use SSE_configuration
    use diagonal_update 
    use lattice 
    use cluster_update
    use local_update
    use linked_list
    use simparameters

    implicit none
    private 

    public one_MCS_plaquette
    public init_SSEconfig_hotstart

    contains

subroutine one_MCS_plaquette(S, beta, Jij_sign, hz_fields_sign, TRANSLAT_INVAR, &
                spins, opstring, config, probtable, plaquettes, vertexlink, &
                leg_visited, hz_fields, C_par_hyperparam)
    ! Purpose:
    ! --------
    ! Perform one Monte-Carlo step (MCS), consisting of 
    ! diagonal and off-diagonal update. 
    ! For the case of a plaquette update with triangular 
    ! plaquettes, we cycle over A-, B-, and C-update for 
    ! completing one MCS. 
    !
    ! Arguments:
    ! ----------
    type(Struct), intent(in) :: S
    real(dp), intent(in) :: beta
    integer, allocatable, intent(in) :: Jij_sign(:,:)
    integer, allocatable, intent(in) :: hz_fields_sign(:)
    logical, intent(in) :: translat_invar
    integer, allocatable, intent(inout) :: spins(:)
    type(t_BondOperator), allocatable, intent(inout) :: opstring(:)    
    type(t_Config), intent(inout) :: config
    type(t_ProbTable), intent(inout) :: probtable      
    type(t_Plaquette), allocatable, intent(in) :: plaquettes(:)
    integer, allocatable, intent(inout) :: vertexlink(:)
    logical, allocatable, intent(inout) :: leg_visited(:)  
    real(dp), intent(in) :: hz_fields(:)
    real(dp), intent(in) :: C_par_hyperparam
        
    integer :: ut ! update type

    ! IMPROVE: no obscure number codes 
    do ut = 111, 113, 1
        ! loop over A-update (=111), B-update (=112) and C-update (=113)
        call diagonal_update_plaquette( S=S, beta=beta, &
            Jij_sign=Jij_sign, hz_fields_sign=hz_fields_sign, spins=spins, opstring=opstring, &
            config=config, probtable=probtable, plaquettes=plaquettes,&
            update_type=ut, TRANSLAT_INVAR=translat_invar, hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam )
        call build_linkedlist_plaquette( &
            opstring=opstring, config=config, &
            vertexlink=vertexlink, leg_visited=leg_visited )        
    call quantum_cluster_update_plaquette( &
            spins=spins, opstring=opstring, vertexlink=vertexlink, &
            leg_visited=leg_visited, config=config, &
            hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam )
! ! ! Local off-diagonal update is not working correctly.
        ! call local_offdiagonal_update(spins=spins, opstring=opstring, config=config, &
        !                 hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam)
    enddo

end subroutine 


subroutine init_SSEconfig_hotstart( S, LL, config, &
    opstring, spins, vertexlink, leg_visited )
    ! ***********************************************************
    ! Initialize the operator string with identities.
    ! Set the spin configuration at propagation step 1 randomly. 
    ! ***********************************************************

    type(Struct), intent(in) :: S
    integer, intent(in) :: LL
    type(t_Config), intent(out) :: config 
    type(t_BondOperator), allocatable, intent(out) :: opstring(:)
    integer, allocatable, intent(out) :: spins(:)
    integer, allocatable, intent(out) :: vertexlink(:)
    logical, allocatable, intent(out) :: leg_visited(:)

    integer :: ip, ir 
    real(dp) :: prob

    config%n_sites = S%Nsites  ! IMPROVE: avoid duplication of information in different structs 
    config%n_exp = 0
    config%LL = LL
    config%n2leg = 0
    config%n2leg_hz = 0
    config%n4leg = 0 
    config%n6leg = 0
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL
    config%n_legs = 2*config%n2leg+2*config%n2leg_hz+4*config%n4leg+6*config%n6leg

    allocate(opstring(LL))
    do ip=1, LL
        opstring(ip)%i = 0; opstring(ip)%j = 0; opstring(ip)%k = 0
        opstring(ip)%optype = IDENTITY
    enddo 
    allocate( spins(S%Nsites) )
    
    ! hot start 
    do ir = 1, S%Nsites
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
    use util
    use tau_embedding 
    use Rydberg
    implicit none 


    type(t_Simparams) :: Sim

    ! Imaginary time correlations 
    type(t_MatsuGrid)  :: MatsuGrid
    
    type(Phys) :: P0
    type(Struct) :: S
    type(t_Kgrid) :: Kgrid

    integer ::  ir, jr, k, i
    integer :: iit, iim
    integer :: ioerr 
    integer :: t1, t2, rate 

    !REMOVE
    character(len=3) :: chr_L
    !REMOVE

    real(dp), allocatable :: J_matrix_out(:,:)

    NAMELIST /SIMPARAMS/ J_1, hx, temp, nx, ny, n_sites, nmeas_step, ntherm_step, Nbin, &
        lattice_type, ignore_Jmatrix, Jmatrix_file, translat_invar, paramscan, &
        scan_min, scan_max, heavy_use, deterministic, &
        hz, hz_fields_file, ignore_hz_fields, Rb, delta, Omega, C_par_hyperparam

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
    ! Check that MPI works correctly
    print*, "My MPI_rank is rank=", MPI_rank, "of size=", MPI_size

    OPEN( UNIT=5, FILE='simparams.in', ACTION="read", STATUS="old", IOSTAT=ioerr)
    IF( ioerr /=0 ) STOP "File simparams.in could not be opened."
    print*, "reading params"
    READ(5, NML=SIMPARAMS)

    ! TODO: Check input parameters ...
    if( J_1==0 .and. ignore_Jmatrix ) then 
        print*, "INCONSISTENT INPUT:"
        print*, "You set `J_1=0` and `ignore_Jmatrix=.true.`. This means that the interactions"
        print*, "are all set to zero, which is probably not what you want."
        print*, "Exiting ..."
        ! stop
    endif 
    call assert( C_par_hyperparam >=0.0_dp )

    if(trim(paramscan) == "parampoint") then 
        print*, "Simulating hx=", hx, "temp=", temp
    elseif(trim(paramscan) == "paramscan_hx")  then 
        hx = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)
    elseif(trim(paramscan) == "paramscan_hz")  then 
        hz = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)        
    elseif(trim(paramscan) == "paramscan_T") then 
        temp = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)
    elseif(trim(paramscan) == "paramscan_delta") then 
        delta = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)
    elseif(trim(paramscan) == "paramscan_L") then 
        nx = 3 + 3 * MPI_rank
        ny = 3 + 3 * MPI_rank
        n_sites = nx * ny * 3  ! only for kagome ! REMOVE
        if (nx < 10 ) then 
            write(chr_L, '(i1)') nx
        else
            write(chr_L, '(i2)') nx
        endif 
        Jmatrix_file = "J"//trim(chr_L)//"x"//trim(chr_L)//"_inplane.txt"
        print*, Jmatrix_file
    else 
        print*, "ERROR: Unknown value of input parameter `paramscan`"
        stop
    endif 

    beta = 1.0_dp / temp
    Sim%J_1=J_1 
    Sim%hx=hx
    Sim%temp=temp
    Sim%beta=1.0_dp/temp
    Sim%nx=nx
    Sim%ny=ny
    Sim%n_sites=n_sites
    Sim%nmeas_step=nmeas_step
    Sim%ntherm_step=ntherm_step
    Sim%Nbin=Nbin
    Sim%lattice_type=lattice_type
    Sim%ignore_Jmatrix=ignore_Jmatrix
    Sim%Jmatrix_file=Jmatrix_file
    Sim%translat_invar=translat_invar
    Sim%paramscan=paramscan
    Sim%scan_min=scan_min
    Sim%scan_max=scan_max
    Sim%heavy_use=heavy_use
    Sim%deterministic=deterministic
    Sim%hz=hz
    Sim%hz_fields_file=hz_fields_file
    Sim%ignore_hz_fields=ignore_hz_fields
    Sim%Rb=Rb
    Sim%delta = delta 
    Sim%Omega = Omega 
    Sim%C_par_hyperparam=C_par_hyperparam

    if (nmeas_step < Nbin) then 
        print*, "Need nmeas_step >= Nbin."
        stop
    endif 

    if (mod(nmeas_step, Nbin) /= 0) then 
        print*, "`nmeas_step` should be an integer multiple of `Nbin`."
        stop
    endif 

    if (trim(Sim%lattice_type) == "triangular") then         
        call init_lattice_triangular(nx=nx, ny=ny, &
            S=S, neigh=neigh, sublattice=sublattice, plaquettes=plaquettes) 
    elseif (trim(Sim%lattice_type) == "kagome") then 
        call init_lattice_kagome(nx=nx, ny=ny, &
            S=S, neigh=neigh, sublattice=sublattice, plaquettes=plaquettes, &
            rvec=rvec) 
    elseif (trim(Sim%lattice_type) == "chain") then 
        call init_lattice_chain(nsites=Sim%n_sites, S=S, neigh=neigh)    
        ! We need a dummy for plaquettes array so that the code does not crash    
        allocate(plaquettes(0))
    else 
        print*, "Error: unknown lattice type"
        stop
    endif 
    if (S%Nsites /= n_sites) then 
        print*, "Error: S%Nsites /= n_sites. Please check the lattice type."
        stop
    endif
    ! IMPROVE: Assign n_plaquettes from the lattice structure S
    ! or better replace config%n_plaquettes everywhere by a component of S
    config%n_plaquettes=size(plaquettes, dim=1)

    ! Specify interactions beyond nearest neighbours 
    ! Nearest neighbour interactions are already taken care 
    ! of by the plaquette operators (for triangular lattice)
    allocate( J_interaction_matrix(S%Nsites,S%Nsites) )
    allocate( hz_fields(1:S%Nsites) )
    if(ignore_Jmatrix) then 
        J_interaction_matrix(:,:) = ZERO
        call init_Rydberg_interactions(S=S, Sim=Sim, &
                hz_fields=hz_fields, Jmatrix=J_interaction_matrix)

        if (MPI_rank == root_rank) then 
            ! Output Jmatrix as input for exact diagonalization code 
            open(101, file="Jmatrix_Rydberg.txt", status="unknown", action="write")
            do ir = 1, n_sites
                write(101, *) ( J_interaction_matrix(ir, jr), jr = 1, n_sites )
            enddo 
            close(101)
        endif 
    else
        if (MPI_rank == root_rank) then 
            open(100, file=trim(Jmatrix_file), action="read", status="old")
            print*, "reading Jmatrix_file", Jmatrix_file
            do i=1,S%Nsites
                read(100, *) J_interaction_matrix(i,1:S%Nsites)
            enddo
            close(100)
        endif
#if defined(USE_MPI)        
        ! IMPROVE: replace MPI_DOUBLE in some way by real(dp) so that 
        ! the code does not break if dp is set to single precision 
        call MPI_BCAST( J_interaction_matrix, size(J_interaction_matrix), &
            MPI_DOUBLE, root_rank, MPI_COMM_WORLD, ierr )
#endif             
    endif 

    ! ! Output the interaction matrix, which is the combination of the 
    ! ! input interaction matrix and the nearest neighbour interactions `J_1`.
    ! allocate( J_matrix_out(S%Nsites,S%Nsites) )   
    ! ! ! REMOVE
    ! ! Interaction matrix with FM next-nearest neighbour interactions
    ! J_matrix_out(:,:) = 0.0_dp
    ! do ir = 1, S%Nsites
    !     do jr = 1, S%Nsites
    !         do k = 1, S%coord ! S%coord + 1, 2*S%coords  => next-nearest neighbours 
    !             if (neigh(k, ir) == jr) then ! nearest neighbours 
    !               J_matrix_out(ir, jr) = -1.0_dp !-0.1_dp
    !             endif 
    !         enddo
    !     enddo         
    ! enddo       
    ! open(101, file="Jmatrix_nnAFM.txt", status="unknown", action="write")
    ! do ir = 1, n_sites
    !     write(101, *) ( J_matrix_out(ir, jr), jr = 1, n_sites )
    ! enddo 
    ! close(101)
    ! J_matrix_out(:,:) = 0.0_dp
    ! ! REMOVE

    ! J_matrix_out(:,:) = J_interaction_matrix(:,:)
    ! do ir = 1, S%Nsites
    !     do jr = 1, S%Nsites
    !         do k = 1, S%coord
    !             if (neigh(k, ir) == jr) then 
    !               ! Nearest neighbour interactions are already taken 
    !               ! care of by the plaquette operators. Include them 
    !               ! here so that J_matrix_out(:,:) can be used as input 
    !               ! for exact diagonalization.
    !               J_matrix_out(ir, jr) = J_1 + J_interaction_matrix(ir, jr)
    !             endif 
    !         enddo
    !     enddo         
    ! enddo   
    ! if( MPI_rank == root_rank) then 
    !     print*, "writing Jmatrix.dat"
    !     open(100, file="Jmatrix.dat", status="unknown", action="write")
    !     do ir = 1, S%Nsites
    !         write(100, *) ( J_matrix_out(ir, jr), jr = 1, S%Nsites )
    !     enddo 
    !     close(100)
    !     deallocate( J_matrix_out )
    !     if( .not. (trim(Sim%lattice_type) == "chain") ) then 
    !         open(201, file="sublattice.dat", status="unknown", action="write")
    !         do ir = 1, S%Nsites 
    !             write(201, *) ir, sublattice(ir)
    !         enddo 
    !         close(201)
    !     endif 
    ! endif 
    ! ! REMOVE

    if (translat_invar) then 
        call make_translat_invar( S, J_interaction_matrix, J_translat_invar )
        allocate( Jij_sign(0:S%Nbravais_sites-1, S%Nbasis*S%Nbasis) )
        where( J_translat_invar > 0 )
            Jij_sign = +1
        elsewhere( J_translat_invar < 0 )
            Jij_sign = -1
        elsewhere 
            Jij_sign = 0
        endwhere        
    else ! not translationally invariant system (e.g. for disorder realization or open BC)
        allocate( Jij_sign(1:S%Nsites, 1:S%Nsites) )
        where( J_interaction_matrix > 0 )
            Jij_sign = +1
        elsewhere( J_interaction_matrix < 0 )
            Jij_sign = -1
        elsewhere 
            Jij_sign = 0
        endwhere
        ! IMPROVE: deallocate( J_interaction_matrix )
    endif 

    ! read the site-dependent longitudinal fields from file 
    ! unless this field is to be ignored 
    if( ignore_hz_fields ) then 
        !!! hz_fields(:) = Sim%hz
        ! hz_fields(:) is written in init_Rydberg_interactions()
    else
        if( MPI_rank == root_rank ) then 
            open(200, file=trim(hz_fields_file), action="read", status="old")
            read(200, *) hz_fields(1:S%Nsites)
            close(200)
        endif 
#if defined(USE_MPI)
        call MPI_BCAST( hz_fields, size(hz_fields), &
            MPI_DOUBLE, root_rank, MPI_COMM_WORLD, ierr )
#endif 
    endif

    allocate(hz_fields_sign(1:size(hz_fields)))
    where( hz_fields > 0 ) 
        hz_fields_sign = +1
    elsewhere( hz_fields < 0 )
        hz_fields_sign = -1
    elsewhere
        hz_fields_sign = 0
    endwhere

    ! seed random number generator:
    !  - with a fixed seed or 
    !  - with the system time (at the millisecond level)
    call init_RNG(MPI_rank, DETERMINISTIC=deterministic) 

    call init_SSEconfig_hotstart( S=S, LL=10, config=config, &
        opstring=opstring, spins=spins, vertexlink=vertexlink, leg_visited=leg_visited )

    ! Precompute the probability tables from which diagonal operators 
    ! will be sampled. 
        call init_probtables( S=S, J_interaction_matrix=J_interaction_matrix, &
        hx=Sim%hx, probtable=probtable, J_1=J_1, &
        n_plaquettes=config%n_plaquettes, TRANSLAT_INV=translat_invar, &
        hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam, spins=spins)
    ! J_interaction_matrix is not needed anymore.
    if( allocated(J_interaction_matrix) ) deallocate( J_interaction_matrix )

    do iit = 1, ntherm_step    
        call one_MCS_plaquette( S=S, beta=beta, Jij_sign=Jij_sign, hz_fields_sign=hz_fields_sign, &
            TRANSLAT_INVAR=translat_invar, spins=spins, opstring=opstring, config=config, & 
            probtable=probtable, plaquettes=plaquettes, vertexlink=vertexlink, leg_visited=leg_visited, &
            hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam )

        if ( (float(config%n_exp) / float(config%LL)) > 2.0_dp / 3.0_dp ) then 
            print*, "Extending cutoff LL_old=",  config%LL
            call extend_cutoff( opstring=opstring, config=config )
            print*, "LL_new=", config%LL
        endif 
    enddo

    call Phys_Init(P0=P0, S=S, Kgrid=Kgrid, MatsuGrid=MatsuGrid, Sim=Sim,&
         Nbin=Nbin, nmeas=nmeas_step)

    call system_clock(count=t1)

    do iim = 1, nmeas_step
        call one_MCS_plaquette( S=S, beta=beta, Jij_sign=Jij_sign, hz_fields_sign=hz_fields_sign, &
            TRANSLAT_INVAR=translat_invar, spins=spins, opstring=opstring, config=config, & 
            probtable=probtable, plaquettes=plaquettes, vertexlink=vertexlink, leg_visited=leg_visited, &
            hz_fields=hz_fields, C_par_hyperparam=C_par_hyperparam )       
        call Phys_Measure(P0, S, Kgrid, MatsuGrid, config, spins, opstring, &
            Sim, probtable%consts_added, heavy_use=heavy_use)

        if ( mod(iim, nmeas_step / Nbin) == 0 ) then 
            call Phys_Avg(P0)
        endif 
    enddo

    ! IMPROVE: Add cluster statistics here (distribution of cluster sizes)
    call system_clock(count=t2, count_rate=rate)
    open(unit=777, file="log_ncpu"//chr_rank//".log", status="Unknown", position="rewind")
    write(777,*) "Running time per measurement MCS: ", (t2 - t1) / real(rate) / real(nmeas_step), "(seconds)"
    close(777)

    call Phys_GetErr(P0)

    call Phys_Print(P0=P0, Kgrid=Kgrid, S=S, MatsuGrid=MatsuGrid, Sim=Sim)

    call deallocate_globals

#if defined(USE_MPI)
    call MPI_FINALIZE(ierr)
#endif 

end program 
