module lattice_kagome
    use SSE_configuration, only: t_Plaquette
    use types
    implicit none
    private 
    save
    
    public t_Plaquette
    public init_lattice_kagome
    public Struct 

    ! Information about lattice structure 
    TYPE Struct   
        character(len=30) :: lattice_type     
        integer :: coord ! coordination number
        integer :: Nbravais_sites 
        integer :: Nbasis
        integer :: Nsites ! Nsites = Nbravais_sites * Nbasis
    END TYPE Struct 

    ! ! for plaquette-based cluster update
    ! type :: t_Plaquette
    !     ! triangular plaquette operator 
    !     integer :: Asite
    !     integer :: Bsite
    !     integer :: Csite
    ! end type t_Plaquette  

    public :: NN, N_unit_cells, N_basis, coord, PI 
    public :: pbc 
    public :: rvec, local_quant_axis
    public :: listB, a1_p, a2_p, r_p, b1_p, b2_p
    public :: coordB, invlistB, neighB, Nx_brav, Ny_brav, lattB_imj

    integer :: NN  ! total number of lattice sites: NN = N_basis * Nx_brav * Ny_brav
    integer :: N_unit_cells ! number of points in the Bravais lattice 
    integer, parameter :: N_basis = 3 ! number of sites in the basis of the kagome lattice (=motif)
    integer, parameter :: coord = 4 ! coordination number of the lattice
    integer, parameter :: rdim = 2! number of spatial dimensions <

    double precision, allocatable :: rvec(:,:)
    double precision, allocatable :: local_quant_axis(:,:)
    
    double precision, parameter :: PI = dacos(-1.d0)
        
    ! Bravais lattice reciprocal lattice for kagome lattice 
    integer :: Nx_brav, Ny_brav
    integer, parameter :: coordB=6 ! coordination number of the Bravais lattice 
    integer, allocatable :: listB(:,:) ! second index of listB(1:2, 0:N_unit_cells-1) starts at "0" !
    integer, allocatable :: invlistB(:,:) ! returns an index from 0 to N_unit_cells-1
    integer, allocatable :: neighB(:,:)
    integer, allocatable :: lattB_imj(:,:) ! map distance between site i and j to distance from origin.
                                           ! (This data structure refers to the underlying Bravais lattice).
    integer, allocatable :: imjBdeg(:) ! How often does the distance imj occur ?                                           
    double precision, allocatable :: x_pos(:), y_pos(:)
    double precision :: a1_p(2), a2_p(2) ! basis vectors of the real space Bravais lattice 
    double precision :: xk_p(2)
    double precision :: r_p(1:N_basis,2) ! coordinates of the sites in the motif (=basis) with which the Bravais lattice is decorated 
    double precision :: b1_p(2), b2_p(2) ! basis vectors of the reciprocal lattice
    double precision, allocatable :: zexpiqrB(:,:)
    double precision, allocatable :: zexpiqrM(:,:,:)
 
    contains

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
    double precision :: rimj_p(1:rdim)
    integer, allocatable :: leftB(:), rightB(:), left_upB(:), right_upB(:), left_downB(:), right_downB(:)

    double precision ::  qvec(1:rdim), ri(1:rdim)

    integer :: plaq_idx, ir1, ir2, ir3
    
    Nx_brav=Nx; Ny_brav=Ny
    N_unit_cells = Nx_brav * Ny_brav
    NN = N_basis * Nx_brav * Ny_brav

    S%lattice_type = "kagome"
    S%coord = 4
    S%Nbasis = 3
    S%Nbravais_sites = nx*ny
    S%Nsites = S%Nbravais_sites * S%Nbasis 

    allocate(neigh(coord, 1:NN))
    allocate(neigh_tmp(coord,0:NN-1))
    allocate(sublattice(1:NN))

    ! ===========================================================
    ! TRIANGULAR BRAVAIS LATTICE     
    ! build a triangular Bravais lattice (B) with periodic boundary conditions 
    ! arrays and variables labelled by "B" are just helper variables 
    allocate(lattB_imj(1:N_unit_cells, 1:N_unit_cells))
    allocate(imjBdeg(1:N_unit_cells))
    allocate(zexpiqrB(1:N_unit_cells,1:N_unit_cells))
    allocate(zexpiqrM(1:N_basis,1:N_basis,1:N_unit_cells))        
    allocate(neighB(coordB, 0:N_unit_cells-1))
    allocate(listB(1:rdim, 0:N_unit_cells-1))
    allocate(invlistB(0:Nx_brav-1, 0:Ny_brav-1))
    allocate(leftB(0:N_unit_cells-1)); allocate(rightB(0:N_unit_cells-1))
    allocate(left_upB(0:N_unit_cells-1)); allocate(right_upB(0:N_unit_cells-1)); 
    allocate(left_downB(0:N_unit_cells-1)); allocate(right_downB(0:N_unit_cells-1))
    do ix=0,Nx_brav-1
      do iy=0,Ny_brav-1
        irB = ix + iy*Nx_brav
        listB(1,irB) = ix
        listB(2,irB) = iy
        invlistB(ix,iy) = irB   
      enddo
    enddo
    ! nearest neighbours on the Bravais lattice 
    do irB = 0,N_unit_cells-1 !!!!! ZERO-INDEXING NECESSARY 
        ix = listB(1, irB)
        iy = listB(2, irB) 
        leftB(irB)       = invlistB( pbc(ix-1, Nx_brav), iy )
        rightB(irB)      = invlistB( pbc(ix+1, Nx_brav), iy )
        left_upB(irB)    = invlistB( pbc(ix-1, Nx_brav), pbc(iy+1, Ny_brav) )
        right_upB(irB)   = invlistB( pbc(ix, Nx_brav), pbc(iy+1, Ny_brav) )
        left_downB(irB)  = invlistB( pbc(ix, Nx_brav), pbc(iy-1, Ny_brav) )
        right_downB(irB) = invlistB( pbc(ix+1, Nx_brav), pbc(iy-1, Ny_brav) )        
        neighB(4,irB) = leftB(irB)       
        neighB(3,irB) = rightB(irB)      
        neighB(5,irB) = left_upB(irB)    
        neighB(6,irB) = right_upB(irB)   
        neighB(1,irB) = left_downB(irB)  
        neighB(2,irB) = right_downB(irB) 
    enddo    
    ! lattB_imj
    do jrB = 0, N_unit_cells-1
        do irB = 0, N_unit_cells-1
            imj_ix = pbc( listB(1,irB) - listB(1,jrB), Nx_brav )
            imj_iy = pbc( listB(2,irB) - listB(2,jrB), Ny_brav )
            lattB_imj(irB+1,jrB+1) = invlistB(imj_ix, imj_iy) + 1
        enddo
    enddo
    ! the degeneracy of LattB_imj
     imjBdeg(:) = 0
     do jrB = 0, N_unit_cells-1
        do irB = 0, N_unit_cells-1
           imj = lattB_imj(irB+1,jrB+1)
           imjBdeg(imj) = imjBdeg(imj) + 1
        end do
     end do
     
     ! Basis vectors of the real-space triangular Bravais lattice.
     ! The length unit is the lattice constant "a" of the kagome lattice
     ! such that the basis vectors of the Bravais lattice have length 2*a.
     a1_p(1) = 1.d0; a1_p(2) = 0.d0
     a2_p(1) = 0.5d0; a2_p(2) = sqrt(3.d0) / 2.d0
      
     ! Reciprocal lattice of the Bravais lattice is also triangular.                    
     ! delta_k in direction of the unit vectors of the reciprocal lattice
     b1_p(1) = 2.d0*PI; b1_p(2) = -2.d0*PI / sqrt(3.d0)
     b2_p(1) = 0.0; b2_p(2) = 4.d0*PI / sqrt(3.d0)

     ! Basis vectors of the motif
     r_p(1,1) = 0.d0; r_p(1,2) = 0.d0
     r_p(2,:) = 0.5d0 * a1_p(:)
     r_p(3,:) = 0.5d0 * a2_p(:)
     
     ! real space position of the sites 
     allocate(rvec(rdim,NN))
     allocate(local_quant_axis(rdim,NN))
     do ir = 1, NN
         sub = mod(ir-1,3) + 1
         ! map site index to index in the Bravais lattice 
         irB = (ir-1)/3
         ! real space position of the site with index ir
         rvec(:,ir) = dble(listB(1,irB))*a1_p + dble(listB(2,irB))*a2_p + r_p(sub,:)     
         
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

     ! listk: grid coordinates of points in reciprocal space 
                         
! ! !      ! zexpiqrB
! ! !      ! Matrix elements of Fourier transform (Bravais lattice)
! ! !      do iq = 1, N_unit_cells
! ! !         qvec = dble(listk(1,iq))*b1_p + dble(listk(2,iq))*b2_p
! ! !         do i = 1, N_unit_cells
! ! !             ri = dble(listB(1,i))*a1_p + dble(listB(2,i))*a2_p
! ! !             zexpiqrB(i, iq) = exp( dcmplx( 0.d0, qvec(1)*ri(1) + qvec(2)*ri(2) ) )
! ! !         enddo
! ! !      enddo
! ! ! 
! ! !                    
! ! !      
! ! !      ! zexpiqrM
! ! !      ! Matrix elements of Fourier transform (Basis = Motif)
! ! !      do iq = 1, N_unit_cells
! ! !         qvec = dble(listk(1,iq))*b1_p + dble(listk(2,iq))*b2_p
! ! !         do subi = 1, N_basis
! ! !            do subj = 1, N_basis 
! ! !               rimj_p = r_p(subi,:) - r_p(subj,:)
! ! !               zexpiqrM(subi, subj, nk) = exp( dcmplx( 0.d0, qvec(1)*rimj_p(1) + qvec(2)*rimj_p(2) ) ) 
! ! !            enddo
! ! !         enddo
! ! !      enddo

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
    do ir=0, NN-1
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
    do ir = 0, NN-1
       neigh(:, ir+1) = neigh_tmp(:, ir) + 1
    enddo
    deallocate(neigh_tmp)
    
    deallocate(leftB); deallocate(rightB); deallocate(left_upB)
    deallocate(right_upB); deallocate(left_downB); deallocate(right_downB)   
    
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
   ! ASSIGN NEAREST NEIGHBOUR BONDS TO PLAQUETTE OPERATORS 
   ! (for triangular lattice)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( .not.(allocated(plaquettes)) ) allocate( plaquettes( S%NBravais_sites ) ) 
    if (size(plaquettes) /= S%Nbravais_sites) then  
        print*, "init_lattice_triangular(): ERROR, inconsistent input"
        print*, "Exiting ..." 
        stop
    endif
    
    plaq_idx = 1
    do irB = 0, S%NBravais_sites-1
        ! Loop over sites on A-sublattice. 
        ir1 = irB*3 + 1
        ir2 = neigh(1, ir1)
        ir3 = neigh(2, ir1)

        plaquettes(plaq_idx) = t_Plaquette(ir1, ir2, ir3)

        plaq_idx = plaq_idx + 1
    enddo 

    end subroutine init_lattice_kagome
    
end module lattice_kagome
    
