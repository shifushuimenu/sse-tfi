! TODO:
!  - Error routines are taken from QUEST code.

module measurements
    use types 
    use SSE_configuration 
    use lattice, only: Struct, t_Kgrid, momentum_grid_triangular_Bravais
    use tau_embedding
    implicit none 

    ! Number of array-valued properties 
    integer, parameter :: narrays = 2

    ! Index of the array variables 
    integer, parameter :: IMEAS = 0  ! Scalar measurements 
    integer, parameter :: ISFA0 = 1  ! Equal-time structure factor 
    integer, parameter :: ISFAM = 2  ! Corr. func. <Az(\tau)Bz(0)> at selected Matsubara frequencies and selected momenta

    integer :: array_sizes(0:narrays)

    ! Parameter for the indexing of scalar variables 
    integer, parameter :: P0_ENERGY = 1
    integer, parameter :: P0_MAGNETIZATION = 2   ! uniform magnetization per site 
    integer, parameter :: P0_COPARAM = 3         ! Clock-order parameter (triangular lattice)

    ! Total number of scalar variables 
    integer, parameter :: P0_N = 3

    ! Name of scalar variables 
    character(len=*), parameter :: P0_STR(P0_N) = (/&
        "          Energy :               ", &
        "          Magnetization :        ", &
        "          Clock order parameter: "/)

    type Phys
        ! Monte Carlo measurements 
        integer :: Nbin         ! number of bins 
        integer :: Nscalar_prop ! total number of scalar measured properties
        integer :: Narray_prop  ! total number of elements of all array-like properties 
        integer :: avg        ! avg = Nbin + 1
        integer :: err        ! err = Nbin + 2 

        integer :: cnt     ! current number of MC measurements *inside* current bin 
        integer :: idx     ! current bin index 

        integer :: Nsites  ! number of sites 
        real(dp) :: beta   ! inverse temperature 

        integer :: N_Matsubara   ! number of selected Matsubara indices for calculating dyn. corr. func.
        integer :: Nq            ! number of selected momentum points  

        ! Scalar array
        real(dp), pointer :: meas(:,:)      ! scalar variables such as energy, etc.

        ! Indices 
        integer :: IARR(0: narrays + 1)

        real(dp), pointer :: AllProp(:,:)   ! array of all real properties 

        ! Equal-time structure factor 
        real(dp), pointer :: Szq(:,:)          
        ! Corr. func. <Az(\tau)Bz(0)> at selected Matsubara frequencies and selected momenta
        ! IMPROVE: Theoretically, this correlation function is complex, but ultimately only 
        ! the real part is needed. 
        real(dp), pointer :: AzBzq_Matsu(:,:)


        logical :: init 
        logical :: compSFA  ! whether to compute structure factor
        logical :: compTAU  ! whether to compute imaginary-time dependent quantities

    end type Phys

    contains 

    subroutine Phys_Init(P0, S, Kgrid, MatsuGrid, beta, Nbin)
    ! 
    ! Purpose:
    ! ========
    ! Initialize the measurement structure Phys.
    !
    ! Precondition:
    ! =============
    !   init_MatsuGrid() and momentum_grid_triangular_Bravais() has been 
    !   called. 
    !        
    ! Arguments:
    ! ==========
    type(Phys), intent(inout) :: P0     ! Phys to be initialized 
    type(Struct), intent(in) :: S
    integer, intent(in) :: Nbin
    real(dp), intent(in) :: beta
    type(t_Kgrid), intent(inout) :: Kgrid 
    type(t_MatsuGrid), intent(inout) :: MatsuGrid

    ! ... Local vars ...
    integer :: i, n

    ! ... Executable ...

    P0%Nbin = Nbin 
    P0%Nsites = S%Nsites 
    P0%beta = beta 

    P0%avg = Nbin + 1
    P0%err = Nbin + 2 
    P0%cnt = 0
    P0%idx = 1 

    P0%Nscalar_prop = P0_N
    
    ! heavy use
    call momentum_grid_triangular_Bravais(S=S, Kgrid=Kgrid)

    call init_MatsuGrid(beta=beta, MatsuGrid=MatsuGrid)
    P0%Nq = Kgrid%Nq
    P0%N_Matsubara = MatsuGrid%N_Matsubara

    ! Allocate storage for properties 
    array_sizes(IMEAS) = 0               ! Scalar variables, do not modify entry 0
    array_sizes(ISFA0) = P0%Nq
    array_sizes(ISFAM) = P0%N_Matsubara * P0%Nq  

    P0%Narray_prop = sum(array_sizes(1:narrays))

    n = P0%Nscalar_prop + P0%Narray_prop
    allocate(P0%AllProp(n, P0%err))

    ! Pointer to beginning of each array 
    P0%IARR(IMEAS) = 1
    ! do i = 1, narrays + 1
        ! P0%IARR(i) = P0%Nscalar_prop + 1 + (i - 1) * array_sizes(i - 1)
    ! enddo
    do i = 1, narrays + 1
        P0%IARR(i) = P0%Nscalar_prop + 1 + array_sizes(i - 1)
    enddo

    P0%meas => P0%AllProp(P0%IARR(IMEAS):P0%IARR(IMEAS + 1) - 1, :)
    P0%Szq  => P0%AllProp(P0%IARR(ISFA0):P0%IARR(ISFA0 + 1) - 1, :)
    P0%AzBzq_Matsu => P0%AllProp(P0%IARR(ISFAM):P0%IARR(ISFAM + 1) - 1, :)

    ! Initialize
    P0%meas = ZERO
    P0%Szq  = ZERO
    P0%AzBzq_Matsu = ZERO

    P0%init = .true.

    end subroutine Phys_Init 


    subroutine Phys_Avg(P0)
    !
    ! Purpose:
    ! ========
    !   Average data within a bin and reset the counter in the bin.
    !
    ! Arguments:
    ! ==========
    type(Phys), intent(inout) :: P0

    ! ... local scalar ...
    real(dp) :: factor 
    integer :: idx

    ! ... Executable ...
    idx = P0%idx 

    ! Compute the normalization factor = 1/cnt 
    if (P0%cnt == 0) then 
        print*, "Phys_Avg: Error: cnt = 0"
        stop 
    endif 
    factor = ONE / P0%cnt 

    ! Average scalar and array-like quantities
    P0%meas(:, idx) = P0%meas(:, idx) * factor 
    P0%AzBzq_Matsu(:, idx) = P0%AzBzq_Matsu(:, idx) * factor

    ! Advance bin
    P0%idx = P0%idx + 1

    ! Reset the counter 
    P0%cnt = 0

    end subroutine Phys_Avg 


    subroutine Phys_GetErr(P0)
        use util, only: JackKnife
    !
    ! Purpose:
    ! ========
    !   Compute averages and errors of the measurements.
    !
    ! Arguments:
    ! ==========
    type(Phys), intent(inout) :: P0

    ! .. Local Scalar ...
    integer :: i, n
    integer :: avg, err
    integer :: nproc 

    ! ... Local Array ...
    real(dp) :: sum_sgn, sgn(P0%Nbin), y(P0%Nbin), data(P0%Nbin)

    ! ... Executable ...
    n = P0%Nbin
    avg = P0%avg
    err = P0%err 
    nproc = 1 

    if (nproc == 1) then 
        
        ! Error of  scalar quantities 
        do i = 1, P0%Nscalar_prop
            data = P0%meas(i, 1:n)
            call JackKnife(n, P0%meas(i, avg), P0%meas(i, err), data, &
                y, sgn, sum_sgn)
        enddo

        ! Error of array-like quantitites 
        do i = 1, P0%Narray_prop
            data = P0%AllProp(i, 1:n)
            call JackKnife(n, P0%AllProp(i, avg), P0%AllProp(i, err), data, &
             y, sgn, sum_sgn )
        enddo 

    else

    endif 

    end subroutine Phys_GetErr

    ! subroutine Phys_Print()

    ! end subroutine 


    subroutine Phys_Measure(P0, S, Kgrid, MatsuGrid, config, spins, opstring, &
                beta, consts_added)
        use util, only: spins2binrep
        use ssetfi_globals, only: sublattice 
        ! Arguments:
        ! ==========
        type(Phys), intent(inout) :: P0
        type(Struct), intent(in) :: S
        type(t_Kgrid), intent(in) :: Kgrid
        type(t_MatsuGrid), intent(in) :: MatsuGrid
        type(t_Config) :: config 
        integer, intent(in) :: spins(:)
        type(t_BondOperator), intent(in) :: opstring(:)
        real(dp), intent(in) :: beta
        real(dp), intent(in) :: consts_added

        ! ... Local variables ...
        real(dp) :: energy, magnz, magnz2
        complex(dp) :: COparam_                  ! complex clock order parameter 
        real(dp) :: COparam                      ! absolute value of the clock order parameter
        complex(dp) :: COphase(3)                ! sublattice phase factors for clock order parameter
        real(dp) :: magnz_tmp
        integer :: tmp_idx, idx 
        integer :: ip, ir, m, q, i1, i2, LL, Nsites
        integer :: spins_tmp(config%N_sites)
        integer :: l_nochange 
        real(dp) :: factor 

        ! temporary help variable => REMOVE later 
        complex(dp), allocatable  :: AzBzq_temp(:,:)
        if (.not.allocated(AzBzq_temp)) allocate( AzBzq_temp(MatsuGrid%N_Matsubara, Kgrid%Nq) ) 

        ! ... Executable ...

        idx = P0%idx
        tmp_idx = P0%avg    ! use the bin which is designed for averages as a temporary storage
                            ! to avoid introducing new variables, e.g. energy, magnetization etc. 

        LL = config%LL
        Nsites = config%N_sites

        energy = (- config%n_exp / beta + consts_added) / float(Nsites)

        spins_tmp(:) = spins(:)
        l_nochange = 0
        factor = ONE / float(LL)

        ! sublattice phases for clock order parameter 
        COphase = (/cmplx(1.0, 0.0), cmplx(-0.5, sqrt(3.0)/2.0), cmplx(-0.5, -sqrt(3.0)/2.0) /)

        magnz = 0.0_dp
        magnz2 = 0.0_dp
        COparam = 0.0_dp
        do ip = 1, LL
            l_nochange = l_nochange + 1
            i1 = opstring(ip)%i
            i2 = opstring(ip)%j
            ! Propagate spins as spin-flip operators are encountered.
            if ((i1.ne.0).and.(i2.eq.0)) then
                ! At propagation steps between spin-flip operators the spin configuration
                ! does not change. Count for how many propagation steps the spin configuration 
                ! stays constant (=l_nochange) and weight the spin config before the next 
                ! spin flip operator by that number. 
                magnz_tmp = sum(spins_tmp(:)) / float(Nsites)
                magnz = magnz + magnz_tmp * l_nochange * factor 
                magnz2 = magnz2 + magnz_tmp**2 * l_nochange * factor 
                COparam_ = complex(0.0_dp, 0.0_dp)
                do ir = 1, Nsites 
                    COparam_ = COparam_ + spins_tmp(ir)*COphase(sublattice(ir))
                enddo 
                COparam = COparam + abs(COparam_) * l_nochange * factor 

                spins_tmp(i1) = -spins_tmp(i1)
                ! reset 
                l_nochange = 0
            endif
        enddo 
        ! Take care of the segment of propagation steps from the last spin 
        ! flip operator up to ip=LL.
        magnz_tmp = sum(spins_tmp(:)) / float(Nsites)
        magnz = magnz + magnz_tmp * l_nochange * factor 
        magnz2 = magnz2 + magnz_tmp**2 * l_nochange * factor 
        COparam_ = complex(0.0_dp, 0.0_dp)
        do ir = 1, Nsites 
            COparam_ = COparam_ + spins_tmp(ir)*COphase(sublattice(ir))
        enddo 
        COparam = COparam + abs(COparam_) * l_nochange * factor 

        COparam = COparam * 3.0_dp / float(Nsites)

        P0%meas(P0_ENERGY, tmp_idx) = energy
        P0%meas(P0_MAGNETIZATION, tmp_idx) = magnz2
        P0%meas(P0_COPARAM, tmp_idx) = COparam

        ! Accumulate result to P0(:, idx)
        P0%meas(:, idx) = P0%meas(:, idx) + P0%meas(:, tmp_idx)


        ! open(100, file='TS.dat', position='append', status='unknown')
        ! write(100, *) energy, magnz, magnz2, spins2binrep(spins), COparam
        ! close(100)   

        ! ===========
        ! heavy use
        ! ===========
        call measure_SzSzTimeCorr(Matsu=MatsuGrid, Kgrid=Kgrid, config=config, S=S, &
            opstring=opstring, spins=spins, beta=beta, chiqAzBz=AzBzq_temp)
        
        ! Flatten the array (and transpose). The index of the flattened array runs for each momentum point
        ! over all Matsubara indices. This convention must be remembered when outputting the array. 
        P0%AzBzq_Matsu(:, tmp_idx) = reshape(source=real(transpose(AzBzq_temp)), shape=(/ size(AzBzq_temp) /))
        ! accumulate result 
        P0%AzBzq_Matsu(:, idx) = P0%AzBzq_Matsu(:, idx) + P0%AzBzq_Matsu(:, tmp_idx)

        P0%cnt = P0%cnt + 1 

    end subroutine Phys_Measure

end module 