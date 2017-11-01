module LinearAlgebra
contains
!___________________________________________________________
function MatProd(mat1,mat2,r1,c2)
!Need something magical
!Dot product of 2D matrices
implicit none
!input
integer                             :: r1 !number of rows of matrix 1
integer                             :: c2 !number of columns of matrix 2
real*8, dimension(:,:)              :: mat1
real*8, dimension(:,:)              :: mat2
!output
real*8, dimension(r1,c2)            :: MatProd
!auxillaries
integer                             :: i,j,k
integer                             :: c1r2 !number of columns of matrix 1

c1r2 = size(mat1,2)
do i=1,r1
    do j=1,c2
        MatProd(i,j)=0
        do k=1,c1r2
            MatProd(i,j)=MatProd(i,j)+mat1(i,k)*mat2(k,j)
        end do
    end do
end do
end function


!______________Copy from other people_______________________
subroutine inverse(n,a,c)
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
!===========================================================
implicit none 
integer n
real*8 a(n,n), c(n,n), a_copy(n,n)
real*8 L(n,n), U(n,n), b(n), d(n), x(n)
real*8 coeff
integer i, j, k


a_copy=a
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
a=a_copy
end subroutine inverse

end module
