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
    integer, parameter :: P0_COPARAM = 3         ! clock-order parameter (triangular lattice)
    integer, parameter :: P0_AV_NEXP = 4         ! average expansion order 
    integer, parameter :: P0_AV_NEXP2 = 5        ! average expansion order squared 
    integer, parameter :: P0_SPECIFIC_HEAT = 6   ! specific heat per site 

    ! Total number of scalar variables 
    integer, parameter :: P0_N = 6

    ! Name of scalar variables 
    character(len=*), parameter :: P0_STR(P0_N) = (/&
        "          Energy (per site):     ", &
        "          Magnetization :        ", &
        "          Clock order parameter: ", &
        "          av. expansion order:   ", &
        "   av. expansion order  squared: ", &
        "       specific heat (per site): "/)

    type Phys
        ! Monte Carlo measurements 
        integer :: Nbin         ! number of bins 
        integer :: nmeas        ! total number of measurement steps 
        integer :: Nscalar_prop ! total number of scalar measured properties
        integer :: Narray_prop  ! total number of elements of all array-like properties 
        integer :: avg        ! avg = Nbin + 1
        integer :: err        ! err = Nbin + 2 
        integer :: ac_time    ! ac_time = Nbin + 3, P0%meas(:, ac_time) stores the running <x^2> and,finally, the AC time. 

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

    subroutine Phys_Init(P0, S, Kgrid, MatsuGrid, beta, Nbin, nmeas)
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
    type(t_Kgrid), intent(inout) :: Kgrid 
    type(t_MatsuGrid), intent(inout) :: MatsuGrid
    real(dp), intent(in) :: beta
    integer, intent(in) :: Nbin         ! Number of bins 
    integer, intent(in) :: nmeas        ! total number of measurement steps 


    ! ... Local vars ...
    integer :: i, n

    ! ... Executable ...

    if( mod(nmeas, Nbin) /= 0 ) then 
         stop "Error, Phys_Init(): `nmeas_step` must be an integer multiple of the number of bins" 
    endif 
    P0%Nbin = Nbin 
    P0%nmeas = nmeas 
    P0%Nsites = S%Nsites 
    P0%beta = beta 

    P0%avg = Nbin + 1
    P0%err = Nbin + 2 
    P0%ac_time = Nbin + 3 
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
    array_sizes(ISFA0) = P0%Nq           ! equal-time structure factor 
    array_sizes(ISFAM) = P0%N_Matsubara * P0%Nq   ! dynamical correlaction function 

    P0%Narray_prop = sum(array_sizes(1:narrays))

    n = P0%Nscalar_prop + P0%Narray_prop
    allocate(P0%AllProp(n, P0%ac_time))

    ! Pointer to beginning of each array 
    P0%IARR(IMEAS) = 1
    P0%IARR(ISFA0) = P0%Nscalar_prop + 1
    do i = 2, narrays + 1
        P0%IARR(i) = P0%IARR(i - 1) + array_sizes(i - 1)
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

    ! Compute fluctuation-dissipation quantity per bin
    P0%meas(P0_SPECIFIC_HEAT, idx) = ( P0%meas(P0_AV_NEXP2, idx) - P0%meas(P0_AV_NEXP, idx)**2 &
                                     - P0%meas(P0_AV_NEXP, idx) ) / P0%Nsites

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
    !   Also compute autocorrelation times with the binning method.
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
        
        ! Error of scalar quantities 
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

        ! Compute autocorrelation times for scalar quantities 
        do i = 1, P0%Nscalar_prop
            data = P0%meas(i, 1:n)
            P0%meas(i, P0%ac_time) = ACT_BinningMethod(  &
                P0%nmeas / P0%Nbin, data, P0%meas(i, P0%ac_time) / real(P0%nmeas, kind=dp) )
        enddo 

    else 
        ! MPI parallelization over measurement steps 
        ! In this case the autocorrelation time is not computed. 

    endif 

    end subroutine Phys_GetErr


    pure function ACT_BinningMethod(Nsamples_per_bin, bin_vals, x2_av_tot) result(AC_time)
    ! Purpose:
    ! ========
    !  Compute the autocorrelation time of quantity x 
    !  with the binning method. 
    ! 
    ! Arguments:
    ! ==========
        integer, intent(in) :: Nsamples_per_bin   ! number of values in each bin 
        real(dp), intent(in) :: bin_vals(:)       ! array of bin averages <x>_{bin_i}
        real(dp), intent(in) :: x2_av_tot         ! <x^2> over the entire Markov chain 
        real(dp) :: AC_time 

    ! ... Local variables ...
        real(dp) :: bin_av_variance    ! variance of the bin averages 
        real(dp) :: total_variance     ! variance over the entire Markov chain 
        real(dp) :: x_av, x2_av        ! average and variance of bin averages 
        real(dp) :: factor  
        integer :: i, n

        n = size(bin_vals, dim=1)
        factor = 1.0 / real(n, kind=dp)
        x_av = 0.0_dp; x2_av = 0.0_dp
        do i = 1, n
            x_av = x_av + bin_vals(i)
            x2_av = x2_av + bin_vals(i)**2
        enddo 
        x_av = x_av * factor 
        bin_av_variance = x2_av * factor - x_av**2
        total_variance = x2_av_tot - x_av**2

        AC_time = Nsamples_per_bin * (bin_av_variance / total_variance)

    end function ACT_BinningMethod

    
    subroutine Phys_Measure(P0, S, Kgrid, MatsuGrid, config, spins, opstring, &
                beta, consts_added, heavy_use)
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
        logical, intent(in) :: heavy_use

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
        P0%meas(P0_AV_NEXP, tmp_idx) = config%n_exp
        P0%meas(P0_AV_NEXP2, tmp_idx) = config%n_exp**2
        ! fluctuation-dissipation quantities are only defined *over* a bin,
        ! not in a single configuration 
        P0%meas(P0_SPECIFIC_HEAT, tmp_idx) = 0.0_dp   

        ! Accumulate result to P0(:, idx)
        P0%meas(:, idx) = P0%meas(:, idx) + P0%meas(:, tmp_idx)
        P0%meas(:, P0%ac_time) = P0%meas(:, P0%ac_time) + P0%meas(:, tmp_idx)**2


        ! open(100, file='TS.dat', position='append', status='unknown')
        ! write(100, *) energy, magnz, magnz2, spins2binrep(spins), COparam
        ! close(100)   

        if( heavy_use ) then 
            call measure_SzSzTimeCorr(Matsu=MatsuGrid, Kgrid=Kgrid, config=config, S=S, &
                opstring=opstring, spins=spins, beta=beta, chiqAzBz=AzBzq_temp)
            
            ! Flatten the array (and transpose). The index of the flattened array runs for each momentum point
            ! over all Matsubara indices. This convention must be remembered when outputting the array. 
            P0%AzBzq_Matsu(:, tmp_idx) = reshape(source=real(transpose(AzBzq_temp)), shape=(/ size(AzBzq_temp) /))
            ! accumulate result 
            P0%AzBzq_Matsu(:, idx) = P0%AzBzq_Matsu(:, idx) + P0%AzBzq_Matsu(:, tmp_idx)
        endif 

        P0%cnt = P0%cnt + 1 

    end subroutine Phys_Measure


    subroutine Phys_Print(P0, Kgrid, S, MatsuGrid, hx, temp, hz, heavy_use)
        use MPI_parallel, only: MPI_rank, chr_rank, root_rank     ! global variables, IMPROVE
        implicit none 

        ! Arguments:
        ! ==========
        type(Phys), intent(in)        :: P0 
        type(t_Kgrid), intent(in)     :: Kgrid 
        type(Struct), intent(in)      :: S
        type(t_MatsuGrid), intent(in) :: MatsuGrid 
        logical, intent(in)           :: heavy_use
        ! IMPROVE: combine parameters into a struct
        real(dp), intent(in) :: hx, temp, hz

        ! ... Local variables ...
        integer :: obs, m, q, idx

        real(dp) :: qvec(2)

        print*, hx, temp, hz, &
        P0%meas(P0_ENERGY, P0%avg), P0%meas(P0_ENERGY, P0%err), &
        P0%meas(P0_MAGNETIZATION, P0%avg), P0%meas(P0_MAGNETIZATION, P0%err), &
        P0%meas(P0_COPARAM, P0%avg), P0%meas(P0_COPARAM, P0%err), &
        P0%meas(P0_SPECIFIC_HEAT, P0%avg), P0%meas(P0_SPECIFIC_HEAT, P0%err)

        ! Averages, error bars and autocorrelation times for all scalar quantities 
        open(500, file='averages'//chr_rank//'.dat', position='append', status='unknown')
        write(500, *) hx, temp, hz,  &
            ( P0%meas(obs, P0%avg), P0%meas(obs, P0%err), P0%meas(obs, P0%ac_time), obs = 1, P0%Nscalar_prop )
        close(500)

        if( heavy_use ) then 
            ! Write out average of imaginary time correlation function 
            open(100, file="Sqz_matsu"//chr_rank//".dat", position="append", status="unknown")
            do m = 1, MatsuGrid%N_Matsubara
                write(100, *) MatsuGrid%im_Matsubara(m), ( P0%AzBzq_Matsu((m-1)*Kgrid%Nq + q, P0%avg), &
                            P0%AzBzq_Matsu((m-1)*Kgrid%Nq + q, P0%err), q=1, Kgrid%Nq )
            enddo 
            write(100, *)
            write(100, *)
            close(100)

            ! Write out all bins of the imaginary time correlation function 
            ! in a format which is useful for postprocessing by a code for analytical continuation.
            ! For each momentum point write out all bins, separated by empty lines. 
            open(200, file="corr_matsu"//chr_rank//".dat", position="append", status="unknown")
            do q=1, Kgrid%Nq
                qvec = Kgrid%listk(1,q)*S%b1_p + Kgrid%listk(2,q)*S%b2_p
                do idx = 1, size(P0%AzBzq_Matsu, dim=2) - 2   ! loop over bins, last two bins are avg and err      
                    write(200, *) q, qvec(1), qvec(2)
                    do m = 1, MatsuGrid%N_Matsubara
                        write(200, *) MatsuGrid%im_Matsubara(m),  P0%AzBzq_Matsu((m-1)*Kgrid%Nq + q, idx)
                    enddo 
                    write(200, *)
                enddo
            enddo 
            write(200, *)
            write(200, *)
            close(200)
        endif 

        if( MPI_rank == root_rank) then 
            open(700, file="output.txt", status="unknown", action="write")
            write(700, *) "The columns in the file averageXXXXX.dat have the following meaning:"
            write(700, *) "======================"
            write(700, *) "Simulation parameters:"
            write(700, *) "======================"
            write(700, *) "transverse field        : 1"
            write(700, *) "temperature             : 2"
            write(700, *) "longitudinal field      : 3"
            write(700, *) "==============================================="
            write(700, *) "Scalar observables: avg, err, (binning) AC time"
            write(700, *) "==============================================="        
            do obs = 1, P0%Nscalar_prop
                write(700, *) P0_STR(obs), 3 +(obs-1)*3 + 1, 3 +(obs-1)*3 + 2, 3 +(obs-1)*3 + 3
            enddo
            close(700)
        endif 


    end subroutine     

end module 