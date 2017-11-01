!This is the Gauss-Newton nonlinear regression for Ar-rigid H2o the model specified
!in the manuscript
module GaussNewton
use ChangeCoord
use Pot_h2oar

contains
!_______________________________________
function IfOneCentralOnly(central, pointX, paras, a123b)
!The contribution of one central data to point X, if unweighed...
implicit none
!input
real*8, dimension(:)           :: central !dimension 4
real*8, dimension(:)           :: pointX  !dimension 3
real*8, dimension(:)           :: paras   !dimension 6
real*8, dimension(:)           :: a123b   !dimension 4
!output
real*8                         :: IfOneCentralOnly
!auxillaries
real*8, dimension(3)           :: X

X(:) = pointX(:)-central(:)

IfOneCentralOnly = a123b(1)*X(1)*exp(-paras(1)*X(1)**2) &
                  +a123b(2)*X(2)*exp(-paras(2)*X(2)**2) &
                  +a123b(3)*X(3)*exp(-paras(3)*X(3)**2) &
                  +a123b(4)*exp(-paras(4)*X(1)**2)        &
                  +a123b(4)*exp(-paras(5)*X(2)**2)        &
                  +a123b(4)*exp(-paras(6)*X(3)**2)        &
                  -2*a123b(4)

end function
!________________________________________________________
function RandomClosestMirror(centralData,n)
!find the n close data to the central data that are added to symmetry
implicit none
!input
real*8, dimension(:)   :: centralData !in internal coords
integer                :: n
!output
real*8, dimension(n,4) :: RandomClosestMirror
!auxillaries
real*8, dimension(4)   :: centralMirror
real*8, dimension(n,4) :: resultMirror
integer                :: i

centralMirror(1) = centralData(1)
centralMirror(2) = centralData(3)
centralMirror(3) = centralData(2)
centralMirror(4) = centralData(4)

resultMirror = RandomClosest(centralMirror,n)

do i=1,n
    RandomClosestMirror(i,1) = resultMirror(i,1)
    RandomClosestMirror(i,2) = resultMirror(i,3)
    RandomClosestMirror(i,3) = resultMirror(i,2)
    RandomClosestMirror(i,4) = resultMirror(i,4)
end do

end function
!________________________________________________________
function RandomClosest(centralData,n)
!find random neighbors of the central data
implicit none
!input
real*8, dimension(:)   :: centralData !in internal coords
integer                :: n
!output
real*8, dimension(n,4) :: RandomClosest
!auxillaries
real*8, dimension(3)   :: deviation, temp
integer                :: i,j
real*8                 :: delta = 0.3

do i = 1,n
    do j = 1,3
        call random_number(deviation(j))
        RandomClosest(i,j) = centralData(j) + (deviation(j)-0.5)*delta
    end do
    temp = RandomClosest(i,1:3)
    RandomClosest(i,4) = PotHartree(IntToSph(temp))
end do
end function
!________________________________________________________
function nClosestUltra(centralData,internalData,n)
implicit none
!input
real*8, dimension(:)        :: centralData
real*8, dimension(:,:)      :: internalData
integer                     :: n
!output
real*8, dimension(n,4)      :: nClosestUltra
!auxillaries
real*8, dimension(2968,4)   :: internalDataDuplicated
integer                     :: i

!duplicate internal data
do i = 1,1484
    internalDataDuplicated(i,:) = internalData(i,:)
end do

do i = 1485,2968
    internalDataDuplicated(i,1) = internalData(i,1)
    internalDataDuplicated(i,2) = internalData(i,3)
    internalDataDuplicated(i,3) = internalData(i,2)
    internalDataDuplicated(i,4) = internalData(i,4)
end do

nClosestUltra = nClosest(centralData,internalDataDuplicated,n)

end function
!________________________________________________________
function nClosestMirror(centralData,internalData,n)
implicit none
!Find the n closest data for the central data added due to symmetry
!input
real*8, dimension(:)   :: centralData
real*8, dimension(:,:) :: internalData
integer                :: n
!output
real*8, dimension(n,4) :: nClosestMirror
!auxillaries
real*8, dimension(4)   :: centralMirror
real*8, dimension(n,4) :: resultMirror
integer                :: i

centralMirror(1) = centralData(1)
centralMirror(2) = centralData(3)
centralMirror(3) = centralData(2)
centralMirror(4) = centralData(4)

resultMirror = nClosest(centralMirror,internalData,n)

do i=1,n
    nClosestMirror(i,1) = resultMirror(i,1)
    nClosestMirror(i,2) = resultMirror(i,3)
    nClosestMirror(i,3) = resultMirror(i,2)
    nClosestMirror(i,4) = resultMirror(i,4)
end do

end function
!________________________________________________________
function nClosest(centralData,internalData,n)
implicit none
!input
real*8, dimension(:)   :: centralData
real*8, dimension(:,:) :: internalData !the big data base to search from
integer                :: n
!output
real*8, dimension(n,4) :: nClosest !The closest n configurations to the central datum
!auxillary
integer                :: i, whereLarge, canonOrNot
real*8, dimension(n)   :: acceptedDist
real*8                 :: thisDist
real*8, dimension(4)   :: thisData 

!sphericals = IntToSph(centralData(1:3))


do i = 1,n+1
    acceptedDist(i) = 1000 !any large number lah
end do

do i = 1,size(internalData,1)
    thisData(:) = internalData(i,:)
    thisDist = dsqrt((thisData(1)-centralData(1))**2+(thisData(2)-centralData(2))**2+(thisData(3)-centralData(3))**2)
    if (thisDist < maxval(acceptedDist) .and. thisDist > 0.00005) then
        whereLarge = maxloc(acceptedDist,1)
        nClosest(whereLarge,:) = thisData(:)
        acceptedDist(whereLarge) = thisDist
    end if
end do
end function

!_______________________________________________
function MultiStart(centralData,nClosest,p)
!try multiple initial guesses of the parameters
implicit none
!input
real*8, dimension(:)   :: centralData !dimension 4
real*8, dimension(:,:) :: nClosest    !dimension (n,4)
integer                :: p           !number of parameters
!output
real*8, dimension(p)   :: MultiStart
!auxillaries
integer                :: i,j,n
real*8                 :: SSE,min_SSE,RMSE
real*8, dimension(p)   :: initialGuess 
real*8, dimension(p+1) :: currentParas
integer                :: numGuess = 200    !number of guesse

n = size(nClosest,1)

min_SSE = 1000
do i = 1,numGuess
   !randomly choose an intial guess
   do j = 1,p
       call random_number(initialGuess(j))
       initialGuess(j) = initialGuess(j)*10
   end do
   !iterate to find the best fit starting from the current intial guess
   currentParas = Iterate(centralData,n,nClosest,p,initialGuess)
   !compare with the result of other intial guesses
   SSE = currentParas(p+1)
   if (SSE<min_SSE) then
       min_SSE = SSE
       MultiStart(:) = currentParas(1:p)
   end if
end do

RMSE = dsqrt(min_SSE/n)
!print *,RMSE
!print *,'The best fit parameters are',MultiStart
end function
!_______________________________________________________
function Iterate(centralData,n,nClosest,p,paras)
!Iterates once only when user say so, just to see how good the method works
implicit none
!input
real*8, dimension(:)   :: centralData !dimension 4
integer                :: n           !number of data
real*8, dimension(:,:) :: nClosest    !dimension (n,4)
integer                :: p           !number of parameters
real*8, dimension(:)   :: paras       !dimension p = 6, it's the initial guess of the parameters and will be renewed after each iteration.
!output
real*8, dimension(p+1)   :: Iterate   !the last slot for SSE
!auxillaries
real*8, dimension(p)   :: correction    !correction to the previous  guess of parameters
real*8, dimension(n,3) :: X             !data configurations
real*8, dimension(n)   :: dY            !energy difference between true and predicted at the n data points
real*8, dimension(n,p) :: DD            !the derivative matrix
real*8, dimension(4)   :: a123b         !(A1,A2,A3,B)
real*8                 :: SSE,min_SSE,RMSE      !sum square of error, root mean square error
integer                :: i, j, k
!control parameters 
integer                :: rounds = 100       !number of rounds
!real*8                 :: tol = 0.0000001 !tolerance on SSE


a123b = A1A2A3B(centralData) !this remains unchanged during parameter optiization


min_SSE = 1000
k=0
do while (k <= rounds .and. isNaN(SSE)==.false.)
    call XY(centralData,nClosest,n,X,dY,a123b,paras)
    SSE = 0
    do i = 1,n
        SSE = SSE + dY(i)**2
    end do
    if (isNaN(SSE)==.false. .and. SSE < min_SSE) then
        min_SSE = SSE
        Iterate(p+1) = min_SSE
        Iterate(1:p) = paras(:)
    end if
!rint *,'. Enter 1 to continue. Enter 0 to exit.'
!rint *,'Current parameters: ',paras 
!ead *, conti
    DD = DeriMat(paras,X,n,p,a123b)
    correction = StepwiseCorrection(n,p,DD,dY)
    do j = 1,p
        paras(j) = paras(j) + correction(j)
        if (paras(j)<0) then
            paras(j)=(paras(j)-correction(j))/2
        end if
        if (paras(j)>5000) then
           paras(j)=5000-(5000-(paras(j)-correction(j)))/2
        end if
    end do
    k = k+1
end do

RMSE = dsqrt(min_SSE/n)
end function
!_______________________________________________________
function StepwiseCorrection(n,p,derimat,dY)
use LinearAlgebra
!This is essentially the amount the the parameter estimate should change in each
!refinement step in the Gauss-Newton method. (the vector b in the textbook by
!Kutner, Nachtsheim, Neter)
implicit none
!input
integer                :: n        !number of data
integer                :: p        !number of parameters
real*8, dimension(:,:) :: derimat  !the derivative matrix (dimension (n,p))
real*8, dimension(:)   :: dY       !the difference between true energies of the data and those predicted by the model (dimension n)
!output
real*8, dimension(p)   :: StepwiseCorrection
!auxillaries
real*8, dimension(p,n) :: deriT   !transpose of derimat
real*8, dimension(p,p) :: dTd,dTdInv
real*8, dimension(p,n) :: dTdInvdT
real*8, dimension(n,1) :: dYvec
real*8, dimension(p,1) :: resultVec
INTEGER                :: I

deriT  = transpose(derimat)
dTd    = MatProd(deriT,derimat,p,p)
call inverse(p,dTd,dTdInv)
dTdInvdT = MatProd(dTdInv,deriT,p,n)
dYvec(:,1) = dY(:)
resultVec = MatProd(dTdInvdT,dYvec,p,1)
StepwiseCorrection(:) = resultVec(:,1)

end function
!________________________________________________________
subroutine XY(centralData,nClosest,n,X,dY,a123b,paras)
!get the data ready for regression. X is an (n,3) matrix giving the geometry
!at each of the n data surrounding the central data. The components
!in each row is (r_x_i - r_central_i), i = 1,2,3.
!dY is the difference between true and predicted energies at each of the n surrounding data.
implicit none
!input
real*8, dimension(:)   :: centralData !dimension(4)
real*8, dimension(:,:) :: nClosest    !dimension(10,4)
integer                :: n           !number of surrounding data 
real*8, dimension(:)   :: a123b       !A1, A2, A3, B
real*8, dimension(:)   :: paras       !alpha1, alpha2 alpha3, beta1, beta2, beta3

!output
real*8, dimension(n,3) :: X
real*8, dimension(n)   :: dY
!auxillaries
integer                :: i,j
real*8, dimension(n)   :: Y ! true energy
real*8, dimension(n)   :: predict !energy predicted by model

do i = 1,n
    do j = 1,3
        X(i,j) = nClosest(i,j)-centralData(j)
    end do
    Y(i) = nClosest(i,4)
    predict(i) = a123b(1)*X(i,1)*exp(-paras(1)*X(i,1)**2) &
                +a123b(2)*X(i,2)*exp(-paras(2)*X(i,2)**2) &
                +a123b(3)*X(i,3)*exp(-paras(3)*X(i,3)**2) &
                +a123b(4)*exp(-paras(4)*X(i,1)**2)        &
                +a123b(4)*exp(-paras(5)*X(i,2)**2)        &
                +a123b(4)*exp(-paras(6)*X(i,3)**2)        &
                -2*a123b(4)
    dY(i) = Y(i)-predict(i)
end do

end subroutine
!_______________________________________________________
function DeriMat(paras,dataX,nData,pPara,a123b)
!The derivative matrix, specific to our model!!! Don't use it anyhow!!!
implicit none
!input
real*8, dimension(:)   :: paras !The guessed parameters of last step (dimension is pPara)
real*8, dimension(:,:) :: dataX !The molecular geometries of the data to be fitted (dimension is (nData, 3).)
integer                :: nData !Number of data
integer                :: pPara !Number of parameters (6 for our model)
real*8, dimension(:)   :: a123b !(A1,A2,A3,B)
!output
real*8, dimension(nData,pPara) :: DeriMat !The derivative matrix as predicted by the currently guessed values of parameters
                                          !Its dimension is (nData, 6) in our model.
!auxillaries
integer                        :: i, j
real*8                         :: B

B  = a123b(4)
!The first 3 columns. corresponding to Del E/Del alpha_1,2,3
do i =1,nData
    do j = 1,3
        DeriMat(i,j) = -a123b(j)*dataX(i,j)**3*exp(-paras(j)*dataX(i,j)**2)
    end do
end do
!The right 3 columns, corresponding to Del E/Del beta_1,2,3 
do j = 1,3
    do i = 1,nData
        DeriMat(i,j+3) = -B*dataX(i,j)**2*exp(paras(j+3)*dataX(i,j)**2)
    end do
end do

end function
!_______________________________________________________
function A1A2A3B(centralData)
!compute A1,A2,A3 and B in the model wrt the central data
implicit none
!input
real*8, dimension(:) :: centralData !Central data arranged as (ArO, ArH1,ArH2,E)
!output
real*8, dimension(4) :: A1A2A3B !The parameter set (A1,A2,A3,B)
!auxillary
real*8               :: delta = 0.00001
real*8,dimension(3,2):: GradSeconds 
integer              :: i

!A1,A2,A3 are just the derivatives wrt each dimension of the molecular geometry
gradSeconds = GradientAndDiagSeconds(centralData,delta)
do i = 1,3
    A1A2A3B(i) = gradSeconds(i,1)
end do
!B is the energy at central data
A1A2A3B(4)   = centralData(4)
end function
!________________________________________________________
function GradientAndDiagSeconds(centralData,delta)
!3D gradient and diagonal of the Hessian metrix.
implicit none
!computes the gradient of the central data in internal coordinates
!input
real*8, dimension(:) :: centralData !Central data arranged as (ArO, ArH1,ArH2,E)
real*8               :: delta       !size of step up and down, must be >0.
!output
real*8, dimension(3,2) :: GradientAndDiagSeconds    !(First column) The partial derivatives on each coordinate, the unit is au for energy/length
                                                  !(Second column)  the diagonal second derivatives
!auxillaries
real*8, dimension(6,4) :: upDowns
real*8, dimension(3,2) :: upDownDerivatives
integer                :: num_dimensions = 3, i, j
real*8                 :: centralEnergy
real*8, dimension(3)   :: temp

!construct the upDowns
do i = 1,num_dimensions
    do j = 1,num_dimensions
        if (j==i) then
            upDowns(2*i-1,j) = centralData(j)-delta
            upDowns(2*i,j)   = centralData(j)+delta
        else
            upDowns(2*i-1,j) = centralData(j)
            upDowns(2*i,j)   = centralData(j)
        end if
    end do
end do

!convert internals to sphericals and calculate the energy of the upDowns geometries (the last column)
do i = 1,num_dimensions*2
    temp = IntToSph(upDowns(i,1:num_dimensions))
    upDowns(i,num_dimensions+1) = potHartree(temp)
end do

!calculate the derivative along each coordinate
do i = 1,num_dimensions
    GradientAndDiagSeconds(i,1) = (upDowns(2*i,num_dimensions+1)-upDowns(2*i-1,num_dimensions+1))/(2*delta)
end do

!calculate the up and down derivatives.
centralEnergy = potHartree(IntToSph(centralData(1:num_dimensions)))

do i = 1,num_dimensions
    upDownDerivatives(i,1) = (centralEnergy-upDowns(2*i-1,num_dimensions+1))/delta
    upDownDerivatives(i,2) = (upDowns(2*i,num_dimensions+1)-centralEnergy)/delta
end do

!calculate the diagonal second derivatives.
do i = 1,num_dimensions
    GradientAndDiagSeconds(i,2) = (upDownDerivatives(i,2)-upDownDerivatives(i,1))/delta
end do
end function

end module

