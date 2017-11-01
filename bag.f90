
!program bag
!implicit none
!real*8,dimension(3,2) :: bag1
!real*8,dimension(4,2) :: bag2
!integer             :: n = 2, num1 = 3, num2 = 0
!integer             :: i,j,bt

!call random_seed()
!bag1(1,:) = (/1.0,2.0/)
!bag1(2,:) = (/3.0,4.0/)
!bag1(3,:) = (/5.0,6.0/)

!do i = 1,4
!    do j = 1,2
!        bag2(i,j) = -1.0
!    end do
!@end do

!call BagTransfer(bag1,bag2,n,num1,num2)

!contains

module bag
contains

subroutine BagTransfer(bag1,bag2,n,num1,num2)
!randomly transfers n elements from bag1 to mag2
implicit none
!input
real*8, dimension(:,:),intent(inout)  :: bag1
real*8, dimension(:,:),intent(inout)  :: bag2 !2 bags. Empty slot marked by (-1,0).
integer, intent(in)                 :: n
integer, intent(inout)                :: num1, num2 !number of elements in bag1 and bag2
!auxillaries
integer               :: i, j, k
real*8                :: rdm
integer               :: rdmInt

do i = 1,n 
    call random_number(rdm)
    rdmInt = ceiling(rdm*dfloat(num1-i+1))
    bag2(num2+i,:) = bag1(rdmInt,:)
    do j = rdmInt,num1-i
        bag1(j,:) = bag1(j+1,:)
    end do
    do j = num1-i+1,size(bag1,1)
        do k = 1,size(bag1,2)
            bag1(j,k) = -1
        end do
    end do
end do

num1 = num1 - n
num2 = num2 + n
end subroutine

!end program
end module
