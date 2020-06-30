!-------------------------!
module mod_sort
!---------------------!
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

implicit none
public :: quicksort
private :: partition

contains

recursive subroutine quicksort(A)
  double precision, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call partition(A, iq)
     call quicksort(A(:iq-1))
     call quicksort(A(iq:))
  endif
end subroutine quicksort

subroutine partition(A, marker)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i = 0
  j = size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine partition

end module mod_sort

