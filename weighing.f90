module weighing
use GaussNewton
use LinearReg
contains
!_____________________________________________________________
function Predict_linear(centrals,InternalData,pointX,numCentrals,expo,alphas,betas)
implicit none
!input
real*8, dimension(:,:)        :: centrals      !dimension (s,4) where s is the number numCentrals
real*8, dimension(:,:)        :: InternalData  !dimension (15..,4)
real*8, dimension(:)          :: pointX        !dimension 3
real*8                        :: expo
integer                       :: numCentrals
real*8, dimension(:)          :: alphas, betas !dimension 3, they are essentials centralsParas in the Predict function
!output
real*8                        :: Predict_linear  !Predicted energy at pointX
!auxillaries
integer                                 :: n=12,p=6
real*8, dimension(numCentrals,6)        :: centralsAB !dimension(s,6)
real*8, dimension(numCentrals)          :: centralsRMSE
real*8, dimension(numCentrals)          :: weights
integer                                 :: i
real*8, dimension(12,4)                 :: closest
real*8, dimension(5,2)                  :: shareholders
real*8                                  :: averageShareholderRMSE

do i = 1,numCentrals
    if (mod(i,2)==1) then
        closest = nClosest(centrals(i,:),InternalData,n)
        !closest  = RandomClosest(centrals(i,:),n)
    else
        closest = nClosestMirror(centrals(i,:),InternalData,n)
        !closest  = RandomClosestMirror(centrals(i,:),n)
    end if
    centralsAB(i,:) = AB_calc(centrals(i,:),closest,alphas,betas,n,p) !n=10,p=6
    centralsRMSE(i) = CentralRMSE(centrals(i,:),closest,alphas,betas,n,centralsAB(i,:))
end do

weights = InverseDistanceWeight(numCentrals,pointX,centrals,expo)
shareholders = FindShareholders(weights)
write (20,*)           'average shareholder RMSE:'
averageShareholderRMSE = 0
do i = 1,5
    averageShareholderRMSE = averageShareholderRMSE + centralsRMSE(shareholders(i,1))
end do
averageShareholderRMSE = averageShareholderRMSE/5*1000000
write (20,"(10f10.2)") averageShareholderRMSE

Predict_linear = 0
do i = 1,numCentrals
    Predict_linear = Predict_linear &
                   + IfOneCentralOnly_linear(centrals(i,:),pointX,alphas,betas,centralsAB(i,:)) * weights(i)
end do

end function
!_______________________________________
function FindShareholders(weights)
!all the large shareholders of a test point
implicit none
!input
integer                                 :: numCentrals
real*8, dimension(:)                    :: weights
!output
real*8, dimension(5,2)                 :: FindShareholders
!auxillaries
integer                                 :: i, wheresmall
real*8                                  :: howsmall

numCentrals = size(weights,1)

do i = 1,5
    FindShareholders(i,2) = -1.0
end do

wheresmall = 1
howsmall   = -1

do i = 1,numCentrals
    if (weights(i)>howsmall) then
        Findshareholders(wheresmall,1) = dfloat(i)
        FindShareholders(wheresmall,2) = weights(i)
        howsmall = minval(FindShareholders(:,2))
        wheresmall = minloc(FindShareholders(:,2),1)
    end if    
end do

end function
!_______________________________________
function Predict(centrals,InternalData,pointX,numCentrals,expo)
!this is nonlinear regression
use GaussNewton
!seriously predicting the energy at any pointX.
implicit none
!input
real*8, dimension(:,:)        :: centrals !dimension (s,4) where s is the number numCentrals
real*8, dimension(:,:)        :: InternalData
real*8, dimension(:)          :: pointX   !dimension 3
real*8                        :: expo
integer                       :: numCentrals
!output
real*8                        :: Predict  !Predicted energy lolololol
!auxillaries
real*8, dimension(numCentrals,4)        :: centralsA123b !dimension(s,4)
real*8, dimension(numCentrals,6)        :: centralsParas 
real*8, dimension(numCentrals)          :: weights
integer                                 :: i
real*8, dimension(10,4)                 :: closest

do i = 1,numCentrals
    if (mod(i,2)==1) then
        closest = nClosest(centrals(i,:),InternalData,10)
    else
        closest = nClosestMirror(centrals(i,:),InternalData,10)
    end if 
    centralsParas(i,:) = MultiStart(centrals(i,:),closest,6)
    centralsA123B(i,:) = A1A2A3B(centrals(i,:))
end do

weights = InverseDistanceWeight(numCentrals,pointX,centrals,expo)

Predict = 0
do i = 1,numCentrals
    Predict = Predict + IfOneCentralOnly(centrals(i,:),pointX,centralsParas(i,:),centralsA123b(i,:)) * weights(i)
end do
end function

!_______________________________________
function InverseDistanceWeight(numCentrals, pointX, centrals, expo)
implicit none
!input
integer                :: numCentrals
real*8, dimension(:,:) :: centrals !(numCentrals,4)
real*8, dimension(:)   :: pointX   !dimension 3
real*8                 :: expo
!output
real*8, dimension(numCentrals) :: InverseDistanceWeight
!auxillary
real*8, dimension(numCentrals) :: unnormalizedWeight
integer                        :: i, j
real*8                         :: normalizingFactor

do i = 1,numCentrals
    if (Distance(pointX,centrals(i,:))==0) then
        do j = 1,numCentrals
            if (i==j) then
                unnormalizedWeight(j) = 1
            else
                unnormalizedWeight(j) = 0
            end if
        end do
        exit
    else
        unnormalizedWeight(i) = (Distance(pointX,centrals(i,:)))**(-expo)
    end if
end do

normalizingFactor = sum(unnormalizedWeight)
InverseDistanceWeight(:) = unnormalizedWeight(:)/normalizingFactor

end function
!____________________________________
function Distance(pointX,central)
implicit none
!input
real*8, dimension(:) :: pointX
real*8, dimension(:) :: central
!output
real*8               :: Distance
!auxillaries
integer              :: i, numDim

distance = 0
numDim = size(pointX)
do i = 1,numDim
    Distance = Distance + (pointX(i)-central(i))**2
end do
Distance = dsqrt(Distance)

end function
!_____________________________________
end module
