!This is for linear regression that determines A1_3, B1_3, with alphas
!and betas fixed.

module linearReg 
contains
!_______________________________________
function IfOneCentralOnly_linear(central, pointX, alphas, betas, AB)
!The contribution of one central data to point X, if unweighed...
implicit none
!input
real*8, dimension(:)           :: central !dimension 4
real*8, dimension(:)           :: pointX  !dimension 3, the unknown point
real*8, dimension(:)           :: alphas  !dimension 3
real*8, dimension(:)           :: betas   !dimension 3
real*8, dimension(:)           :: AB      !dimension 6
!output
real*8                         :: IfOneCentralOnly_linear
!auxillaries
real*8, dimension(3)           :: d

d(:) = pointX(:)-central(:)

IfOneCentralOnly_linear = AB(1)*d(1)*exp(-alphas(1)*d(1)**2) &
                        +AB(2)*d(2)*exp(-alphas(2)*d(2)**2) &
                        +AB(3)*d(3)*exp(-alphas(3)*d(3)**2) &
                        +AB(4)*(exp(-betas(1)*d(1)**2)-1)   &
                        +AB(5)*(exp(-betas(2)*d(2)**2)-1)   &
                        +AB(6)*(exp(-betas(3)*d(3)**2)-1)   &
                        +central(4)

end function
!_______________________________________________________
function CentralRMSE(centralData,nClosest,alphas,betas,n,AB)
!RMSE at a sample point
implicit none
!input
real*8, dimension(:)   :: centralData !dimension(4)
real*8, dimension(:,:) :: nClosest !dimension(n,4)
real*8, dimension(:)   :: alphas, betas
integer                :: n
real*8, dimension(:)   :: AB      !dimension 6
!output
real*8                 :: CentralRMSE
!zuxillaries
real*8, dimension(n)   :: err    !the error vector Y-Yhat

err = error(centralData,nClosest,alphas,betas,n,AB)
CentralRMSE = CalcRMSE(err)

end function
!_______________________________________________________
function CalcRMSE(err)
!the standard RMSE function
implicit none
!input
real*8, dimension(:) :: err ! the error vector, dimension = n
!output
real*8               :: CalcRMSE
!auxillaries
integer              :: n   ! dimension of err
integer              :: i

n = size(err)
CalcRMSE = 0
do i = 1,n
    CalcRMSE = CalcRMSE + err(i)**2
end do
CalcRMSE = dsqrt(CalcRMSE/n)

end function
!_______________________________________________________
function error(centralData,nClosest,alphas,betas,n,AB)
!tge error vector
implicit none
!input
real*8, dimension(:)   :: centralData !dimension(4)
real*8, dimension(:,:) :: nClosest !dimension(n,4)
real*8, dimension(:)   :: alphas, betas
integer                :: n
real*8, dimension(:)   :: AB      !dimension 6
!output
real*8, dimension(n)   :: error    !the error vector Y-Yhat
!auxillaries
real*8, dimension(n)   :: esti
integer                :: i

do i = 1,n
    esti(i) = IfOneCentralOnly_linear(centralData, nClosest(i,1:3), alphas, betas, AB)
    error(i) = nClosest(i,4) - esti(i)
end do

end function
!_______________________________________________________
function AB_calc(centralData,nClosest,alphas,betas,n,p)
implicit none
!input
real*8, dimension(:)   :: centralData !dimension(4)
real*8, dimension(:,:) :: nClosest !dimension(n,4)
integer                :: n        !number of surrounding data
integer                :: p        !number of parameters (currently = 6)
real*8, dimension(3)   :: alphas, betas
!output
real*8, dimension(p)   :: AB_calc
!auxillaries
real*8, dimension(n,p) :: X
real*8, dimension(p)   :: Y

call XY(centralData,nClosest,n,p,alphas,betas,X,Y)
AB_calc = PEstimate(X,Y,n,p)

end function
!=======================================================
function PEstimate(X,Y,n,p)
use LinearAlgebra
!Use linear regression to estimate parameters A1_3, B1_3
implicit none
!input
integer                :: n        !number of surrounding data (lphas,betas=10)
integer                :: p        !number of parameters (=6)
real*8, dimension(:,:) :: X        !The X matrix, definition see next subroutine or my notebook, dimension (n,p)
real*8, dimension(:)   :: Y        !The Y vector, definition see next subroutine or my notebook, dimension (n)

!output
real*8, dimension(p)   :: PEstimate !estimated A1_3, B1_3 (dimension = p = 6)
!auxillaries
real*8, dimension(p,n) :: XT   !transpose of derimat
real*8, dimension(p,p) :: XTX,XTX_Inv
real*8, dimension(p,n) :: XTX_Inv_XT
real*8, dimension(n,1) :: Yvec
real*8, dimension(p,1) :: resultVec

XT  = transpose(X)
XTX = MatProd(XT,X,p,p)
call inverse(p,XTX,XTX_Inv)
XTX_Inv_XT   = MatProd(XTX_Inv,XT,p,n)
Yvec(:,1)    = Y(:)
resultVec    = MatProd(XTX_Inv_XT,Yvec,p,1)
PEstimate(:) = resultVec(:,1)

end function
!________________________________________________________
subroutine XY(centralData,nClosest,n,p,alphas,betas,X,Y)
implicit none
!computes the X matrix and the Y vector, alphas and betas given and fixed
!input
real*8, dimension(:)   :: centralData !dimension(4)
real*8, dimension(:,:) :: nClosest !dimension(n,4)
integer                :: n        !number of surrounding data (currently =
integer                :: p        !number of parameters (currently = 6)
real*8, dimension(:)   :: alphas, betas
!output
real*8, dimension(:,:) :: X        !matrix for independent variables
real*8, dimension(:)   :: Y        !vector of the dependent variables minus the energy of the central data
!auxillaries
integer                :: i,j      !iterator
real*8                 :: E0       !energy of the central data
real*8                 :: d        !distance in 


E0 = centralData(4)

!Construct Y
do i = 1,n
    Y(i) = nClosest(i,4) - E0
end do

do i = 1,n
    do j = 1,3
        d = nClosest(i,j)-centralData(j)
        X(i,j)   = d*exp(-alphas(j)*d**2)
        X(i,j+3) = exp(-betas(j)*d**2)-1
    end do
end do

end subroutine
!________________________________________________________

end module

