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


subroutine init_SSEconfig_hostart( S, LL, config, &
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
    config%n4leg = 0 
    config%n6leg = 0
    config%n_ghostlegs = MAX_GHOSTLEGS*config%LL
    config%n_legs = 2*config%n2leg+4*config%n4leg+6*config%n6leg

    allocate(opstring(LL))
    do ip=1, LL
        opstring(ip)%i = 0; opstring(ip)%j = 0; opstring(ip)%k = 0
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
    use util, only: init_RNG
    use tau_embedding 
    implicit none 

    ! ****************************************
    ! simulation parameters                  !
    ! ***********************************************
    REAL(dp) :: J_1 = +1.0_dp    
    real(dp) :: hx = 0.60_dp   
    real(dp) :: temp = 0.1_dp
    real(dp) :: beta 
    integer :: nx, ny, n_sites             

    integer :: nmeas_step = 10000
    integer :: ntherm_step = 10000 
    integer :: Nbin = 100
    character(len=10) :: lattice_type = "triangular"
    logical :: ignore_Jmatrix = .FALSE.
    character(len=30) :: Jmatrix_file = "Jmatrix.txt"
    character(len=12) :: paramscan
    real(dp) :: scan_min = 0.0, scan_max = 0.0
    logical :: heavy_use = .FALSE.
    logical :: deterministic = .FALSE.

    ! Imaginary time correlations 
    type(t_MatsuGrid)  :: MatsuGrid
    ! ***********************************************

    type(Phys) :: P0
    type(Struct) :: S
    type(t_Kgrid) :: Kgrid

    integer ::  ir, jr, k, i
    integer :: iit, iim
    integer :: ioerr 
    integer :: t1, t2, rate 

    real(dp), allocatable :: J_matrix_out(:,:)

    NAMELIST /SIMPARAMS/ J_1, hx, temp, nx, ny, n_sites, nmeas_step, ntherm_step, Nbin, &
        lattice_type, ignore_Jmatrix, Jmatrix_file, paramscan, &
        scan_min, scan_max, heavy_use, deterministic

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

    if(trim(paramscan) == "parampoint") then 
        print*, "Simulating hx=", hx, "temp=", temp
    elseif(trim(paramscan) == "paramscan_hx")  then 
        hx = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)
    elseif(trim(paramscan) == "paramscan_T") then 
        temp = scan_min + MPI_rank * (scan_max - scan_min)/ float(MPI_size)
    elseif(trim(paramscan) == "paramscan_L") then 
        nx = 3 + 3 * MPI_rank
        ny = 3 + 3 * MPI_rank
        n_sites = nx * ny * 3  ! only for kagome ! REMOVE
    else 
        print*, "ERROR: Unknown value of input parameter `paramscan`"
        stop
    endif 
    beta = 1.d0 / temp


    if (nmeas_step < Nbin) then 
        print*, "Need nmeas_step >= Nbin."
        stop
    endif 

    if (mod(nmeas_step, Nbin) /= 0) then 
        print*, "`nmeas_step` should be an integer multiple of `Nbin`."
        stop
    endif 

    if (trim(lattice_type) == "triangular") then         
        call init_lattice_triangular(nx=nx, ny=ny, &
            S=S, neigh=neigh, sublattice=sublattice, plaquettes=plaquettes) 
    elseif (trim(lattice_type) == "kagome") then 
        call init_lattice_kagome(nx=nx, ny=ny, &
            S=S, neigh=neigh, sublattice=sublattice, plaquettes=plaquettes) 
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
    ! of by the plaquette operators. 
    allocate( J_interaction_matrix(S%Nsites,S%Nsites) )
    allocate( Jij_sign(S%Nsites, S%Nsites) )
    if(ignore_Jmatrix) then 
        J_interaction_matrix(:,:) = ZERO
    else
        if (MPI_rank == root_rank) then 
            open(100, file=trim(Jmatrix_file), action="read", status="old")
            print*, "reading ", Jmatrix_file
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
    ! J_interaction_matrix(:,:) = 0.0_dp
    ! do ir = 1, S%Nsites
    !     do jr = 1, S%Nsites
    !         do k = 1, S%coord
    !             if (neigh(k, ir) == jr) then 
    !               ! nearest neighbour interactions are already taken 
    !               ! care of by the plaquette operators 
    !               J_interaction_matrix(ir, jr) = -0.0_dp
    !             endif 
    !         enddo
    !     enddo         
    ! enddo   


    ! Output the interaction matrix, which is the combination of the 
    ! input interaction matrix and the nearest neighbour interactions `J_1`.
    allocate( J_matrix_out(S%Nsites,S%Nsites) )   

    ! ! ! REMOVE
    ! ! Interaction matrix with FM next-nearest neighbour interactions
    ! J_matrix_out(:,:) = 0.0_dp
    ! do ir = 1, S%Nsites
    !     do jr = 1, S%Nsites
    !         do k = S%coord+1, 2*S%coord 
    !             if (neigh(k, ir) == jr) then ! next-nearest neighbours 
    !               J_matrix_out(ir, jr) = -0.1_dp
    !             endif 
    !         enddo
    !     enddo         
    ! enddo       
    ! open(101, file="Jmatrix_nnFM.dat", status="unknown", action="write")
    ! do ir = 1, n_sites
    !     write(101, *) ( J_matrix_out(ir, jr), jr = 1, n_sites )
    ! enddo 
    ! close(101)
    J_matrix_out(:,:) = 0.0_dp
    ! REMOVE
    J_matrix_out(:,:) = J_interaction_matrix(:,:)
    do ir = 1, S%Nsites
        do jr = 1, S%Nsites
            do k = 1, S%coord
                if (neigh(k, ir) == jr) then 
                  ! Nearest neighbour interactions are already taken 
                  ! care of by the plaquette operators. Include them 
                  ! here so that J_matrix_out(:,:) can be used as input 
                  ! for exact diagonalization.
                  J_matrix_out(ir, jr) = J_1 + J_interaction_matrix(ir, jr)
                endif 
            enddo
        enddo         
    enddo   
    if( MPI_rank == root_rank) then 
        print*, "writing Jmatrix.dat"
        open(100, file="Jmatrix.dat", status="unknown", action="write")
        do ir = 1, S%Nsites
            write(100, *) ( J_matrix_out(ir, jr), jr = 1, S%Nsites )
        enddo 
        close(100)
        deallocate( J_matrix_out )
        open(201, file="sublattice.dat", status="unknown", action="write")
        do ir = 1, n_sites 
            write(201, *) ir, sublattice(ir)
        enddo 
        close(201)
    endif 
    ! REMOVE

    ! seed random number generator with the system time (at the millisecond level)
    call init_RNG(MPI_rank, DETERMINISTIC=deterministic) 

    ! Precompute the probability tables from which diagonal operators 
    ! will be sampled. 
    call init_probtables( S=S, J_interaction_matrix=J_interaction_matrix, &
        hx=hx, probtable=probtable, Jij_sign=Jij_sign, J_1=J_1, &
        n_plaquettes=size(plaquettes,dim=1), TRANSLAT_INV=.FALSE.)

    call init_SSEconfig_hostart( S=S, LL=10, config=config, &
        opstring=opstring, spins=spins, vertexlink=vertexlink, leg_visited=leg_visited )

    do iit = 1, ntherm_step    
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

    call Phys_Init(P0=P0, S=S, Kgrid=Kgrid, MatsuGrid=MatsuGrid, beta=beta, Nbin=Nbin)

    call system_clock(count=t1)

    do iim = 1, nmeas_step
        call one_MCS_plaquette( beta=beta, Jij_sign=Jij_sign, spins=spins, &
            opstring=opstring, config=config, probtable=probtable, &
            plaquettes=plaquettes, vertexlink=vertexlink, &
            leg_visited=leg_visited )
        call Phys_Measure(P0, S, Kgrid, MatsuGrid, config, spins, opstring, &
            beta, probtable%consts_added, heavy_use=heavy_use)

        if (mod(iim, nmeas_step / Nbin) == 0) then 
            call Phys_Avg(P0)
        endif 
    enddo

    call system_clock(count=t2, count_rate=rate)
    open(unit=777, file="log_ncpu"//chr_rank//".log", status="Unknown", position="rewind")
    write(777,*) "Running time per measurement MCS: ", (t2 - t1) / real(rate) / real(nmeas_step), "(seconds)"
    close(777)

    call Phys_GetErr(P0)

    call Phys_Print(P0=P0, Kgrid=Kgrid, S=S, MatsuGrid=MatsuGrid, hx=hx, &
            temp=temp, heavy_use=heavy_use)

    call deallocate_globals

#if defined(USE_MPI)
    call MPI_FINALIZE(ierr)
#endif 

end program 
