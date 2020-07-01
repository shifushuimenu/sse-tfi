module tau_embedding
    use types 
implicit none 
private 

public t_MatsuGrid 
public measure_SzSzTimeCorr
public init_MatsuGrid

type t_MatsuGrid
    integer :: N_Matsubara 
    integer, allocatable :: im_Matsubara(:)
    integer :: L_tau   ! IMPROVE: move this vairable to another, logically more suitable place 
    real(dp) :: dtau   ! IMPROVE: move this vairable to another, logically more suitable place 
end type t_MatsuGrid

contains

subroutine measure_SzSzTimeCorr(Matsu, Kgrid, config, S, &
    opstring, spins, beta, chiqAzBz)
! *********************************************************************
! Purpose:
! --------
! Compute dynamical correlation functions in Matsubara frequency space.
! (via "embedding of SSE configuration in imaginary tme")

! Method:
! -------
! From a uniform distribution over [0,beta] assign n random
! imaginary times to the n SSE operators in an operator string.
! Since the SSE ensembles and the ensembles in imaginary time are
! equivalent, Green's functions in imaginary time can be computed 
! via Fourier transform.
! 
! Reference: 
! ~~~~~~~~~~
!    Michel and Evertz, arXiv:0705.0799v2.
!
! Furthermore, compute the auxiliary correlation function <Sx(l*dtau)Sx(0)> .
! For an explanation of this method see:
!
! Reference:    
! ~~~~~~~~~~    
!    Appendix A in Phys. Rev. B 92, 245137 (2015).
!
! Precondition:
! -------------
! The lattice initialization routine has been called and the Matsubara
! frequency grid has been initialized.
!
! Arguments:
! ----------
!   Matsu:  Matsubara struct object containing information of
!           the (potentially inhomogeneous) grid of Matsubara frequencies.
!   Kgrid:  Kgrid object containing the momentum points for which 
!           the correlation functions are to be computed. 
!   config: SSE configuration struct object
!   S:      Lattice struture object 
!   chiqAzBz: The imaginary-time correlation function for selected Matsubara
!             frequencies and selected momentum points.
use SSE_configuration
use lattice, only: Struct, t_Kgrid
use mod_sort, only: quicksort
use util, only: PI
implicit none 

type(t_MatsuGrid), intent(in) ::  Matsu
type(t_Kgrid), intent(in) :: Kgrid
type(t_Config), intent(in) :: config 
type(Struct), intent(in) :: S
type(t_Bondoperator), intent(in) :: opstring(:)
integer, intent(in) :: spins(:)
real(dp), intent(in) :: beta
complex(dp), intent(out) :: chiqAzBz(1:Matsu%N_Matsubara,1:Kgrid%Nq) 

! IMPROVE: This should acutally be a dummy argument 
real(dp), allocatable :: SxSxTimeCorr(:)

! ... Local variables ...
integer :: itau, im, ir, ip, i1, i2, q
real(dp) :: tau_p1, tau_p2

integer :: spins2(size(spins,1))
complex(dp) :: Asum( Matsu%N_Matsubara, S%Nsites )
complex(dp) :: Bsum( Matsu%N_Matsubara, S%Nsites )
! Spatial cosine and sine transforms of Asum(im,1:NN) and Bsum(im,1:NN)
complex(dp) :: CA( Matsu%N_Matsubara, S%Nsites )
complex(dp) :: CB( Matsu%N_Matsubara, S%Nsites )
complex(dp) :: SA( Matsu%N_Matsubara, S%Nsites )
complex(dp) :: SB( Matsu%N_Matsubara, S%Nsites )  
real(dp) :: tau_values( config%n_exp )
real(dp) :: harvest( config%n_exp )

real(dp) :: AsumRe_magnz, BsumRe_magnz

! Instead of the imaginary time correlation function <Sx(tau)Sx(0)> 
! we measure the auxiliary correlation function
! that is defined in terms of a discrete time grid laid over the tau-interval [0,beta].
integer :: itau1, itau2, lPos
real(dp) :: tau_dist
integer :: counter( S%Nsites )
real(dp), allocatable :: tau_values_of_Sx(:,:)
integer, allocatable :: N_SxSx(:,:)

real(dp)  :: m_ip, magnz_
real(dp) :: chi_magnz
integer :: itau_Sx

integer :: n, Nq, LL, n_exp, n2leg, N_Matsubara, L_tau
real(dp) :: dtau 

! ... Executable ...

N_Matsubara = Matsu%N_Matsubara
L_tau = Matsu%L_tau
dtau = Matsu%dtau

n = S%Nsites
Nq = Kgrid%Nq
LL = config%LL
n_exp = config%N_exp
n2leg = config%n2leg

if (.not.allocated(SxSxTimeCorr)) allocate(SxSxTimeCorr(L_tau))

allocate(tau_values_of_Sx(1:max(10,4*int(n2leg/n)), n)) ! It is assumed that there are not more than 4*int(n2leg/NN) spin-flip operators sitting on one site.
! N_SxSx(:,:) is number of pairs (Sx(tau1),Sx(tau2)) of Sx operators in the operator string at site ir where one operator acts
! lPos time intervals later than the other one
allocate(N_SxSx(n, L_tau))  

 counter(:) = 0
 chiqAzBz(:,:) = 0.d0
 ! Kubo susceptibility:  int_0^{beta} magnz(tau) magnz(0) dtau 	where magnz = (S_1 + ... S_NN) / NN 
 chi_magnz = 0.d0

! 1. Choose n_exp random real numbers between 0 and beta
call random_number(harvest)
tau_values(:) = beta*harvest(:)

! 2. Order them in plcase and assign them to all 
! operators (including identities) as their position in continuous imaginary time
call quicksort(tau_values)

! 3. Compute spin-spin correlation function in terms of Matsubara frequencies
Asum(:,:)=0.d0; Bsum(:,:)=0.d0
AsumRe_magnz=0.d0; BsumRe_magnz=0.d0

spins2(:) = spins(:)

! To compute the average magnetization, do not sum over all sites
! at each propagation step since there is only one bond active at each step
! and the other spins do not change.
! It is more efficient to determine the magnetization in the initial spin configuration
! and adapt it only whenever a spin flip operator is encountered.
m_ip = sum(spins2(1:n))

itau=0
do ip=1,LL

  ! imaginary times are only assigned to non-trivial operators
  i1 = opstring(ip)%i
  i2 = opstring(ip)%j
  if ( ((i1.ne.0).or.(i2.ne.0)).or.(ip.eq.LL)) then ! non-trivial operator detected or end of imaginary time is reached

    if (ip.ne.LL) then
      if (itau.eq.0) then
        tau_p1 = 0.d0
      else
        tau_p1 = tau_values(itau)
      endif
      tau_p2 = tau_values(itau+1)
      itau = itau + 1
    else ! tau == beta is reached. The last step of Riemann integration over the interval [tau_p1, beta] also needs to be performed.
      if (itau == 0) then ! not a single spin-flip operator detected (possible at high temp.)
        tau_p1 = 0.d0
        tau_p2 = beta
      else 
        tau_p1 = tau_values(itau)
        tau_p2 = beta
      endif 
    endif

    ! IMPROVE: Call a function to compute A_z(ip) on the spin configuration at propagation step ip.
    ! 1. The function (or even several functions) should be passed to the routine for calculating imaginary 
    ! time correlations as a parameter. Thus, the inner workings of this routine need not be touched 
    ! by the user, who only needs to modify the function. Fortran routines can return arrays, i.e. 
    ! there is no problem to compute e.g. the structure factor for all values of momenta. 
    ! 2. The spin configuration between propagation step ip1 and ip2 onlz changes if there 
    ! is the spin flip operator between ip1 and ip2. The structure factor should not be recomputed 
    ! at every non-trivial propagation step. Rather the user should provide a function how to update 
    ! the structure factor (or any other diagonal quantitiy signified by A_z(ip) ) if the spin a single position changes. 
    ! Actually, the user needs to supply two functions, A_z({S_i}^{z}) and B_z({S_i}^{z}), which both are complex,
    ! as well as functions update_A_z(ir) and update_B_z(ir) for updating the values of A and B after a spin flip
    ! at position ir: 
    !       A^{\prime} = A + update_A(ir) 
    !       B^{\prime} = B + update_B(ir)
    ! In the case of the dynamical structure factor, A- and B-operators are identical. However, in general they could 
    ! be different. 

    ! average magnetization at propagation step ip
    magnz_= m_ip / dble(n)

    AsumRe_magnz = AsumRe_magnz +  magnz_*(tau_p2 - tau_p1)
    BsumRe_magnz = BsumRe_magnz +  magnz_*(tau_p2 - tau_p1)

    ! propagate spins if spin-flip operators are encountered
    if ( (i1.ne.0).and.(i2.eq.0) ) then
      spins2(i1) = -spins2(i1)
      ! update the instantaneous magnetization at propagation step ip
      m_ip = m_ip + 2*spins2(i1)
      ! Record the imaginary time that is assigned to the spin flip operator => for computing <Sx(tau)Sx(0)>
      counter(i1) = counter(i1) + 1
      tau_values_of_Sx(counter(i1),i1) = tau_values(itau)
    endif

   endif  ! non-trivial operator detected
 enddo 

 spins2(:) = spins(:)

 do ir=1,n
 do itau_Sx=0,counter(ir)
    ! This sum is only over non-trivial propagation steps where a 
    ! spin-flip operator sits. 

    if (itau_Sx.eq.0) then
        ! first interval of the Riemann sum [0, tau_first] where tau_first
        ! is the first spin-flip operator acting on site ir
        tau_p1 = 0.d0
    else
        tau_p1 = tau_values_of_Sx(itau_Sx,ir)
        spins2(ir) = -spins2(ir) ! get the right spin state of the imaginary time segment
    endif

    if (itau_Sx.ne.counter(ir)) then
        tau_p2 = tau_values_of_Sx(itau_Sx+1,ir)
    else
	      ! last interval of the Riemann sum [tau_last, beta] where tau_last
        ! is the last spin-flip operator acting on site ir
        tau_p2 = beta
    endif

    ! Matsubara frequency omega_m = 0
    Asum(1,ir) = Asum(1,ir) +  spins2(ir)*(tau_p2 - tau_p1)
    Bsum(1,ir) = Bsum(1,ir) +  spins2(ir)*(tau_p2 - tau_p1)
    ! The other Matsubara frequencies ...
    do im=1,N_Matsubara-1 ! Fortran indices start at 1
      Asum(im+1,ir) = Asum(im+1,ir) + spins2(ir)*( exp(cmplx(0.0_dp, -2*PI*Matsu%im_Matsubara(im+1)*tau_p2/beta, kind=dp)) - &
                        exp(cmplx(0.0_dp, -2*PI*Matsu%im_Matsubara(im+1)*tau_p1/beta, kind=dp)) )
      Bsum(im+1,ir) = Bsum(im+1,ir) + spins2(ir)*( exp(cmplx(0.0_dp, 2*PI*Matsu%im_Matsubara(im+1)*tau_p2/beta, kind=dp)) - &
                        exp(cmplx(0.0_dp, 2*PI*Matsu%im_Matsubara(im+1)*tau_p1/beta, kind=dp)) )
    enddo
     
 enddo
 enddo
 
! Spatial cosine (C) and sine (S) Fourier transformation
! IMPROVE: Replace by discrete Fourier transform in real space
CA(:,:)=cmplx(0.d0, 0.d0, kind=dp); CB(:,:)=cmplx(0.d0, 0.d0, kind=dp)
SA(:,:)=cmplx(0.d0, 0.d0, kind=dp); SB(:,:)=cmplx(0.d0, 0.d0, kind=dp)
do q=1,Nq
  do ir=1,n   
    do im=1, N_Matsubara
      CA(im,q)=CA(im,q)+Kgrid%cosqr(ir,q)*Asum(im,ir)
      CB(im,q)=CB(im,q)+Kgrid%cosqr(ir,q)*Bsum(im,ir)
      SA(im,q)=SA(im,q)+Kgrid%sinqr(ir,q)*Asum(im,ir)
      SB(im,q)=SB(im,q)+Kgrid%sinqr(ir,q)*Bsum(im,ir)      
    enddo
  enddo
enddo  

! Spin-spin correlation function in Matsubara and momentum space. 
do q=1,Nq
  ! Matsubara frequency omega_m = 0
  im=1
  chiqAzBz(1,q) = real( CA(im,q)*CB(im,q) + SA(im,q)*SB(im,q) ) / beta 
  do im=2,N_Matsubara
    chiqAzBz(im,q) = real( CA(im,q)*CB(im,q) + SA(im,q)*SB(im,q) ) &
                    * beta / (2*PI*Matsu%im_Matsubara(im)*2*PI*Matsu%im_Matsubara(im)) 
  enddo 
enddo
chiqAzBz(:,:) = chiqAzBz(:,:) / float(n)

! Kubo static susceptibility:  int_0^{beta} magnz(tau) magnz(0) dtau 	where magnz = (S_1 + ... S_L) / L
  chi_magnz = AsumRe_magnz*BsumRe_magnz / beta
  
! Compute the auxiliary correlation function <Sx(l*dtau)Sx(0)> .
! For an explanation of the method see Appendix A in Phys. Rev. B 92, 245137 (2015).

  N_SxSx(:,:) = 0
  do ir=1,n ! for each site
    do itau1=1,counter(ir)
      do itau2=itau1+1,counter(ir)
        tau_dist = abs( tau_values_of_Sx(itau2,ir) - tau_values_of_Sx(itau1,ir) )
        lPos = int(tau_dist/dtau)+1
        N_SxSx(ir,lPos) = N_SxSx(ir,lPos) + 1 ! number of pairs (Sx(tau1),Sx(tau2)) of Sx operators in the operator string at site ir where one operator acts
            ! lPos time intervals later than the other one
        lPos = L_tau - int(tau_dist/dtau) ! count the same pair of operators again, but this time the distance between them wraps around imaginary time.
        N_SxSx(ir,lPos) = N_SxSx(ir,lPos) + 1
      enddo
    enddo
  enddo

  ! average over all sites
  SxSxTimeCorr(:) = 0.d0
  do lPos=1,L_tau
    SxSxTimeCorr(lPos) = sum(N_SxSx(1:n,lPos)) / (float(S%Nsites)*beta*dtau)
  enddo

deallocate(tau_values_of_Sx)
deallocate(N_SxSx)

end subroutine measure_SzSzTimeCorr


subroutine init_MatsuGrid(beta, MatsuGrid)
    implicit none 
    ! Purpose:
    ! --------
    ! Generate an inhomogeneous grid of Matsubara indices. 
    ! 
    ! Arguments:
    ! ----------
    real(dp), intent(in) :: beta
    type(t_MatsuGrid), intent(out) :: MatsuGrid

    ! ... Local variables ...
    integer :: im, r

    MatsuGrid%N_Matsubara = 123 
    ! indices of Matsubara frequency which are to be computed
    allocate(MatsuGrid%im_Matsubara(MatsuGrid%N_Matsubara)) 
    do im=0,80
        MatsuGrid%im_Matsubara(im+1) = im
    enddo
    r=0
    do im=81,102
    r=r+1
        MatsuGrid%im_Matsubara(im+1) = 80+r*10
    enddo
    r=0
    do im=103,122
    r=r+1
        MatsuGrid%im_Matsubara(im+1) = 300+r*50
    enddo

    ! Evenly-spaced imaginary time grid for computing <Sx(l*dtau)Sx(0)>
    MatsuGrid%L_tau = 200
    MatsuGrid%dtau = beta / float(MatsuGrid%L_tau) ! spacing of the imaginary time grid

end subroutine init_MatsuGrid    

end module tau_embedding 
