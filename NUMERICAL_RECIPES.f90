! ***********************************
! Taken from: 
!       Numerical Recipes
!       in Fortran 90
!       Second Edition
! ***********************************
module nrutil
    use types
    implicit none  

   INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
   INTEGER, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
   INTEGER, PARAMETER :: NPAR_CUMSUM=16
   INTEGER, PARAMETER :: NPAR_CUMPROD=8
   INTEGER, PARAMETER :: NPAR_POLY=8
   INTEGER, PARAMETER :: NPAR_POLYTERM=8

   REAL(dp), PARAMETER :: PI = dacos(-1.d0)
   REAL(dp), PARAMETER :: TWOPI = 2.0_dp * dacos(-1.d0)

    contains 

SUBROUTINE assert(n1,string)
    ! Report and die if any logical is false (used for arg range checking).
    CHARACTER(LEN=*), INTENT(IN) :: string
    LOGICAL, INTENT(IN) :: n1
    if (.not. n1) then
        write (*,*) "nrerror: an assertion failed with this tag:", &
        string
        STOP "program terminated by assert1"
    end if
END SUBROUTINE assert
        
    
SUBROUTINE swap(a,b)
    COMPLEX(dp), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(dp), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
END SUBROUTINE swap
    
    

FUNCTION arth(first,increment,n)
    INTEGER, INTENT(IN) :: first,increment,n
    INTEGER, DIMENSION(n) :: arth
    INTEGER :: k,k2,temp
    if (n > 0) arth(1)=first
    if (n <= NPAR_ARTH) then
        do k=2,n
            arth(k)=arth(k-1)+increment
        end do
    else
        do k=2,NPAR2_ARTH
            arth(k)=arth(k-1)+increment
        end do
        temp=increment*NPAR2_ARTH
        k=NPAR2_ARTH
        do
            if (k >= n) exit
            k2=k+k
            arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
END FUNCTION arth
    
        
FUNCTION assert_eq(n1,n2,string)
    ! Report and die if integers not all equal (used for size checking).
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq
    if (n1 == n2) then
        assert_eq=n1
    else
        write (*,*) "nrerror: an assert_eq failed with this tag:", &
        string
        STOP "program terminated by assert_eq2"
    end if
END FUNCTION assert_eq    
    

FUNCTION zroots_unity(n,nn)
    ! Complex function returning nn powers of the nth root of unity.
    INTEGER, INTENT(IN) :: n,nn
    COMPLEX(dp), DIMENSION(nn) :: zroots_unity
    INTEGER :: k
    REAL(dp) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
        if (k >= nn) exit
        zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),dp)
        zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
            zroots_unity(2:min(k,nn-k))
        k=2*k
    end do
    END FUNCTION zroots_unity
    

RECURSIVE FUNCTION cumsum(arr,seed) RESULT(ans)
    ! Cumulative sum on an array, with optional additive seed.
    REAL(dp), DIMENSION(:), INTENT(IN) :: arr
    REAL(dp), OPTIONAL, INTENT(IN) :: seed
    REAL(dp), DIMENSION(size(arr)) :: ans
    INTEGER :: n,j
    REAL(dp) :: sd
    n=size(arr)
    if (n == 0) RETURN
    sd=0.0_dp
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
        do j=2,n
            ans(j)=ans(j-1)+arr(j)
        end do
    else
        ans(2:n:2)=cumsum(arr(2:n:2)+arr(1:n-1:2),sd)
        ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
END FUNCTION cumsum

end module nrutil

module numerical_recipes
    implicit none 

    public cosft1, sinft 

    contains 
    
SUBROUTINE fourrow(data,isign)
    USE types 
    USE nrutil, ONLY : assert,swap, PI
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    INTEGER :: n,i,istep,j,m,mmax,n2
    REAL(dp) :: theta
    COMPLEX(dp), DIMENSION(size(data,1)) :: temp
    COMPLEX(dp) :: w,wp
    COMPLEX(dp) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, "n must be a power of 2 in fourrow_dp")
    n2=n/2
    j=n2
    do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
            if (m < 2 .or. j < m) exit
            j=j-m
            m=m/2
        end do
        j=j+m
    end do
    mmax=1
    do
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI/(isign*mmax)
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dp)
        w=cmplx(1.0_dp,0.0_dp,kind=dp)
        do m=1,mmax
            ws=w
            do i=m,n,istep
                j=i+mmax
                temp=ws*data(:,j)
                data(:,j)=data(:,i)-temp
                data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
        end do
        mmax=istep
    end do
END SUBROUTINE fourrow
        

SUBROUTINE four1(data,isign)
    USE types
    USE nrutil, ONLY : arth,assert, TWOPI
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(dp), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, "n must be a power of 2 in four1_dp")
    m1=2**ceiling(0.5_dp*log(real(n,dp))/0.693147_dp)
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dp)
    w=cmplx(1.0_dp,0.0_dp,kind=dp)
    do j=2,m2
        w=w*wp+w
        dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
END SUBROUTINE four1
        
SUBROUTINE realft(data,isign,zdata)
    use types
    USE nrutil, ONLY : assert,assert_eq,zroots_unity
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: isign
    COMPLEX(dp), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER :: n,ndum,nh,nq
    COMPLEX(dp), DIMENSION(size(data)/4) :: w
    COMPLEX(dp), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(dp), DIMENSION(:), POINTER :: cdata
    COMPLEX(dp) :: z
    REAL(dp) :: c1=0.5_dp,c2
    n=size(data)
    call assert(iand(n,n-1)==0, "n must be a power of 2 in realft_dp")
    nh=n/2
    nq=n/4
    if (present(zdata)) then
        ndum=assert_eq(n/2,size(zdata),"realft_dp")
        cdata=>zdata
        if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=dp)
    else
        allocate(cdata(n/2))
        cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=dp)
    end if
    if (isign == 1) then
        c2=-0.5_dp
        call four1(cdata,+1)
    else
        c2=0.5_dp
    end if
    w=zroots_unity(sign(n,isign),n/4)
    w=cmplx(-aimag(w),real(w),kind=dp)
    h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)=h1+w(2:nq)*h2
    cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
    z=cdata(1)
    if (isign == 1) then
        cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dp)
    else
        cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dp)
        call four1(cdata,-1)
    end if
    if (present(zdata)) then
        if (isign /= 1) then
            data(1:n-1:2)=real(cdata)
            data(2:n:2)=aimag(cdata)
        end if
    else
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
        deallocate(cdata)
    end if
END SUBROUTINE realft
        


SUBROUTINE cosft1(y)
    USE types
    USE nrutil, ONLY : assert,cumsum,zroots_unity
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
    ! Calculates the cosine transform of a set of N +1 real-valued data points y. The transformed
    ! data replace the original data in array y. N must be a power of 2. This program, without
    ! changes, also calculates the inverse cosine transform, but in this case the output array
    ! should be multiplied by 2/N .
    COMPLEX(dp), DIMENSION((size(y)-1)/2) :: w
    REAL(dp), DIMENSION((size(y)-1)/2-1) :: y1,y2
    REAL(dp) :: summ
    INTEGER :: n,nh
    n=size(y)-1
    call assert(iand(n,n-1)==0, "n must be a power of 2 in cosft1")
    nh=n/2
    w=zroots_unity(n+n,nh)
    summ=0.5_dp*(y(1)-y(n+1))
    y(1)=0.5_dp*(y(1)+y(n+1))
    !  Construct the two pieces of the auxiliary array.
    y1=0.5_dp*(y(2:nh)+y(n:nh+2:-1))
    y2=y(2:nh)-y(n:nh+2:-1)
    !  Carry along this sum for later use in unfolding
    ! the transform.
    summ=summ+sum(real(w(2:nh))*y2)
    y2=y2*aimag(w(2:nh))
    !  Calculate the auxiliary function.
    y(2:nh)=y1-y2
    y(n:nh+2:-1)=y1+y2
    !  Calculate the transform of the auxiliary function.        
    call realft(y(1:n),1)
    y(n+1)=y(2)
    !  summ is the value of F1 in equation (12.3.21).
    y(2)=summ
    ! Equation (12.3.20).
    y(2:n:2)=cumsum(y(2:n:2))
END SUBROUTINE cosft1
        

SUBROUTINE sinft(y)
    USE types
    USE nrutil, ONLY : assert,cumsum,zroots_unity
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
    ! Calculates the sine transform of a set of N real-valued data points stored in array y. The
    ! number N must be a power of 2. On exit y is replaced by its transform. This program,
    ! without changes, also calculates the inverse sine transform, but in this case the output array
    ! should be multiplied by 2/N .
    REAL(dp), DIMENSION(size(y)/2+1) :: wi
    REAL(dp), DIMENSION(size(y)/2) :: y1,y2
    INTEGER :: n,nh
    n=size(y)
    call assert(iand(n,n-1)==0, "n must be a power of 2 in sinft")
    nh=n/2
    ! Calculate the sine for the auxiliary array.
    wi=aimag(zroots_unity(n+n,nh+1))
    y(1)=0.0
    ! Construct the two pieces of the auxiliary array.
    y1=wi(2:nh+1)*(y(2:nh+1)+y(n:nh+1:-1))
    ! Put them together to make the auxiliary array.
    y2=0.5_dp*(y(2:nh+1)-y(n:nh+1:-1))
    y(2:nh+1)=y1+y2
    y(n:nh+1:-1)=y1-y2
    ! Transform the auxiliary array.
    call realft(y,+1)
    ! Initialize the sum used for odd terms.
    y(1)=0.5_dp*y(1)
    y(2)=0.0
    ! Odd terms are determined by this running sum.
    y1=cumsum(y(1:n-1:2))
    ! Even terms in the transform are determined directly.
    y(1:n-1:2)=y(2:n:2)
    y(2:n:2)=y1

END SUBROUTINE sinft
    
end module numerical_recipes 