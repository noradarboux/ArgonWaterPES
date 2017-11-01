program FitForAr
use ChangeCoord
use Pot_h2oar
use GaussNewton 
use weighing
use bag
implicit none
!Original ab-intio calcs data list
real*8, dimension(1584,4) :: InternalData, InternalDataExceptTestSet
!iterators
integer :: i,j,k,ii,ij,ik
!experimental field
real*8, dimension(3) :: sph,cart
real*8, dimension(3,2) :: gradNSeconds
real*8, dimension(3,3) :: WHEREISWATER !PARAMETER
!========================================
!       playground (for declaration)
!========================================
real*8, dimension(3)                :: pointX, weight
real*8                              :: expo, predictE, predictEohne, trueE, RMSE, RMSEohne
real*8, dimension(2000,4)           :: centrals,CentralsOhne
integer                             :: numCentrals,num1,num2
real*8, dimension(3)                :: alphas, betas
real*8, dimension(200,4)            :: testSet    
real*8, dimension(200)              :: errorSigned, errorSignedOhne
real*8, dimension(20,4)             :: tempbag

!========================================
!             end playground
!========================================

call random_seed()
!following is the cartesian coordinate of the rigid water. The origin being the
!CM of water. DO NOT MODIFY.
WHEREISWATER(1,:) =(/0.06512,0.0,0.0/) !POSITION OF O
WHEREISWATER(2,:) =(/0.65131,0.75722,0.0/) !position of H1
WHEREISWATER(3,:) =(/0.65131,-0.75722,0.0/) !position of H2
InternalData = GetInternalData(WHEREISWATER)

!****************************************************
!           PLAYGROUND. CAN MODIFY NOW.
!++++++++++++++++++++++++++++++++++++++++++++++++++++

!arbitrarily specify alphas and betas
alphas = (/0.5,0.5,0.5/)
betas  = (/0.5,0.5,0.5/)
!end specify alphas and betas

open(20,file="Ar-H2O predict linear_fixed test set and sample set_exploring error prediction")

!-------------------------------=============--------------=============
!construct the test set
num1 = size(InternalData,1)
num2 = 0
call RandomCentralDataByBagTransfer(InternalData,testSet,200,num1,num2)
InternalDataExceptTestSet = InternalData !just a copy...
!end construct test set
!----------------------================---------------==================
expo = 15
numCentrals = 0

do j = 1,5
    !sample centrals in the config space with test set excluded
    !!!!numCentrals = 500
    !!!!expo = 2+(j-1)*2
    call RandomCentralDataByBagTransfer(InternalData,centrals,100,num1,numCentrals)
    do k = numCentrals-99,numCentrals
        centralsOhne(k,:) = centrals(k,:)
    end do

    do k = 1,10
        write(20,*) "sample size = ",numCentrals/2,'*2'
        !!!!write(20,*) "expo = ",expo
        !estimate energies of the configs in the test set
        do i = 1,size(testSet,1)
            pointX = testSet(i,1:3)
            trueE = testSet(i,4)
            predictE = Predict_linear(centrals(1:numCentrals,:),InternalDataExceptTestSet,pointX,numCentrals,expo,alphas,betas)
            predictEohne = Predict_linear(centralsOhne(1:numCentrals,:),InternalDataExceptTestSet,pointX,numCentrals,expo,alphas,betas)
            errorSigned(i) = predictE - trueE
            errorSignedOhne(i) = predictEohne - trueE
        end do
        !and calculate RMSE
        RMSE = CalcRMSE(errorSigned)
        RMSEohne = CalcRMSE(errorSignedOhne)
        print *, numCentrals/2,RMSE,RMSEohne
        !!call ErrorDistribution(errorSigned)
        !pick the point with the largest predicted error and add it to the sample set
        call ErrorProbing(InternalData,centrals,num1,numCentrals,InternalDataExceptTestSet,alphas,betas,100)
        numCentrals = numCentrals - 10
        do i = 1,10    
            InternalData(num1+i,:) = centrals(numCentrals+i,:)
            tempbag(i,:) = centrals(numCentrals+i,:)
        end do
        
        call RandomCentralDataByBagTransfer(InternalData,centralsOhne,10,num1,numCentrals)
        ii = num1 - 9
        do while (ii <= num1)
            do ij = 1,10
                if (VectorEq(InternalData(ii,:),tempbag(ij,:))) then
                    num1 = num1 - 1
                    do ik = 1, num1
                        InternalData(ii+ik-1,:) = InternalData(ii+ik,:) !delete this one cause it's in tempbag
                    end do
                    ii = ii -1
                    exit
                end if
            end do
            ii = ii + 1
        end do
        
    end do

    write(21,*) '**********************************************'
    print *,    '**********************************************'
end do

close(20)
!======================================================================
!#####END OF MAIN PROGRAM.####################
!======================================================================
contains
function VectorEq(vec1,vec2)
implicit none
real*8, dimension(:) :: vec1,vec2
logical                :: VectorEq      
integer                :: j

VectorEq = .true.

    do j = 1 , size(vec1)
        if (vec1(j)/=vec2(j)) then
            VectorEq = .false.
            exit
        end if 
    end do

end function
!________________________________________________________________________
function GetInternalData(WHEREISWATER)
implicit none
real*8, dimension(1584,4) :: OriginalData
real*8, dimension(1584,4) :: GetInternalData
character                 :: trash
real*8, dimension(3*3)    :: WHEREISWATER

!read the original data
open(10,file="PESofAr-H2O.txt")
read(10,*) trash
do i=1,1584
    read(10,*) OriginalData(i,1), OriginalData(i,2), OriginalData(i,3), OriginalData(i,4)
end do

!convert the spherical coordinates in the data to internal coordinates. Needed
!for all operations. DO NOT MODIFY.
do i=1,1584
    GetInternalData(i,1:3) = SphToInt(OriginalData(i,1:3),WHEREISWATER)
    GetInternalData(i,4)   = potHartree(IntToSph(GetInternalData(i,1:3)))
end do
end function
!________________________________________________________________________
function RandomTestset(WHEREISWATER,testsize)
implicit none
!input
real*8, dimension(3,3)        :: WHEREISWATER
integer                       :: testsize
!output
real*8, dimension(testsize,4) :: RandomTestset
!auxi
real*8, dimension(3)          :: randomreal
integer                       :: i, j

do j=1,testsize
    do i=1,3
        call random_number(randomreal(i))
    end do
    RandomTestset(j,1) = randomreal(1)*90 !theta
    RandomTestset(j,2) = randomreal(2)*180 !phi
    RandomTestset(j,3) = randomreal(3)*20 !r
    RandomTestset(j,4) = potHartree(RandomTestset(j,1:3))
    RandomTestset(j,1:3) = SphToInt(RandomTestset(j,1:3),WHEREISWATER)
end do

end function
!________________________________________________________________________
subroutine ErrorProbing(InternalData,centrals,IntDataLeftHowmany,numCentrals,IntDataExTest,alphas,betas,sizeprob)
implicit none
!Pick the point with the largest predicted error and add it to the sample set
!input and output
real*8, dimension(:,:)               :: InternalData, centrals, IntDataExTest
integer                              :: IntDataLeftHowmany, numCentrals
real*8, dimension(:)                 :: alphas, betas
!auxillaries
integer                              :: sizeProb !size of the probing set
real*8, dimension(sizeProb,4)        :: probingSet
integer                              :: i, j, numMaCon, n, p, num1, num2, chosenID, chooseHowMany
real*8, dimension(5,4)               :: majorContributors !major 5 contributors
real*8, dimension(5,6)               :: maConAB !dimension(s,6)
real*8, dimension(5)                 :: maConRMSE
real*8, dimension(12,4)              :: closest
real*8, dimension(sizeProb)          :: avrMaConRMSE
real*8, dimension(4)                 :: chosenProb
real*8, dimension(1,4)               :: unchosen

numMaCon=5
n=12
p=6
num2 = 0
chooseHowMany = 5

!100 random probe points from the internal data set that are not already in the sample set
call BagTransfer(InternalData,probingSet,sizeProb,IntDataLeftHowmany,num2)

!For each probe point
do i = 1,sizeProb
    !find its major 5 contributors. They're the 5 sample points closest to the probe point
    majorContributors = nClosest(probingSet(i,:),centrals,numMaCon)
    !for each major contributor
    do j = 1,numMaCon
        closest = nClosestUltra(majorContributors(j,:),IntDataExTest,n)
        maConAB(j,:) = AB_calc(majorContributors(j,:),closest,alphas,betas,n,p)
        maConRMSE(j) = CentralRMSE(majorContributors(j,:),closest,alphas,betas,n,maConAB(j,:))
    end do
    !calculate the average major contributors' RMSE
    avrMaConRMSE(i) = sum(maConRMSE)/dfloat(numMaCon)
end do

do i = 1,chooseHowMany
    !find the probe point with the largest avrMaConRMSE
    chosenID      = maxloc(avrMaConRMSE,1)
    if (avrMaConRMSE(chosenID)<=10.0**(-4.0)) then
    else
        chosenProb(:) = probingSet(chosenID,:)
        !transfer it to centrals
        call AddOnePtToCentrals(centrals,numCentrals,chosenProb)
        !mark the ID as chosen
        avrMaConRMSE(chosenID) = 0
    end if
end do

!transfer the rest back to InternalData
do i=1,sizeProb
    if (avrMaConRMSE(i)/=0) then
        unchosen(1,:) = probingSet(i,:)
        num1 = 1
        call BagTransfer(unchosen,InternalData,1,num1,IntDataLeftHowmany)
    end if
end do
end subroutine
!________________________________________________________________________
subroutine AddOnePtToCentrals(centrals,numCentrals,pt)
implicit none
!input
real*8, dimension(:,:) :: centrals
integer                :: numCentrals
real*8, dimension(:)   :: pt

centrals(numCentrals+1,:) = pt(:)
centrals(numCentrals+2,1) = pt(1)
centrals(numCentrals+2,3) = pt(2)
centrals(numCentrals+2,2) = pt(3)
centrals(numCentrals+2,4) = pt(4)
numCentrals = numCentrals + 2

end subroutine
!________________________________________________________________________
subroutine ErrorDistribution(errorSigned)
implicit none
!distribution of error in log. scale
!input
real*8, dimension(:) :: errorSigned ! error vector of the prediction on the test set, dimension(n)
!output
real*8, dimension(12,2) :: ErrorDistri
!auxillaries
integer                :: n, logError, i

open(21,file = 'H2O-Ar error distribution2')

do i = 1,12
    ErrorDistri(i,1) = 10.0**(dfloat(i)-12.0) 
    ErrorDistri(i,2) = 0
end do

n = size(errorSigned)

do i = 1,n
    logError = ceiling(log10(abs(errorSigned(i))))
    if (logError<-11) then
        ErrorDistri(1,2) = ErrorDistri(1,2) + 1
    else if (logError>-1) then
        ErrorDistri(12,2) = ErrorDistri(12,2) + 1
    else
        ErrorDistri(logError+12,2) = ErrorDistri(logError+12,2) + 1
    end if 
end do

do i = 1,12
    ErrorDistri(i,2) = (ErrorDistri(i,2)/n)*100.0
    write (21,"(es10.1,f10.3)") ErrorDistri(i,:)
end do

write (21,*) ''
end subroutine

!________________________________________________________
function randomCentralData(InternalData,numCentrals)
implicit none
!to pick any numCentrals central data
!input
real*8, dimension(:,:), intent(in) :: InternalData !The ab-initio data in internal coords
integer                            :: numCentrals  !The sample size we want
!output
!the central datas where the basis function is situated, in internal coords
!(ArO,ArH1,ArH2)
real*8, dimension(numCentrals*2,4) :: randomCentralData !times 2 because of perm symmetry
!auxillaries
integer                         :: dataSize, i, randomInd, j, duplicate
integer, dimension(numCentrals) :: indices
real*8                          :: randomReal

dataSize = size(InternalData,1)

i=0
do while (i<numCentrals)
    call random_number(randomReal)
    randomInd = ceiling(randomReal * dataSize)
    !check duplicate
    duplicate = 0
    do j = 1,i
        if (indices(j)==randomInd) then
            duplicate = 1
            exit
        end if
    end do
    !if no duplicate, add to the indices list
    if (duplicate==0) then
        i = i+1
        indices(i) = randomInd
    end if
end do

!writing the coordinates and energies of the central data in to the output
!and incorporating the per. inv. symmetry
do i = 1,numCentrals
    randomCentralData(2*i-1,:) = InternalData(indices(i),:)
    randomCentralData(2*i,1)   = randomCentralData(2*i-1,1)
    randomCentralData(2*i,2)   = randomCentralData(2*i-1,3)
    randomCentralData(2*i,3)   = randomCentralData(2*i-1,2)
    randomCentralData(2*i,4)   = randomCentralData(2*i-1,4)
end do
end function randomCentralData
!____________________________________________________
subroutine RandomCentralDataByBagTransfer(InternalData,centrals,n,num1,currentNumCentrals)
use bag
implicit none
real*8, dimension(:,:) :: InternalData
real*8, dimension(:,:) :: centrals
integer                :: n, num1, currentNumCentrals, num2
!auxillary
real*8, dimension(n,4) :: bag2
integer                :: i

num2 = 0
call BagTransfer(InternalData,bag2,n/2,num1,num2)

bag2 = CentralSym(bag2(1:n/2,:),n/2)

do i = 1,n
    centrals(currentNumCentrals+i,:) = bag2(i,:)
end do

currentNumCentrals = currentNumCentrals + n
end subroutine
!____________________________________________________
function CentralSym(centralOri,numCentrals)
implicit none
!input
integer                            :: numCentrals
real*8, dimension(numCentrals,4)   :: centralOri
!output
real*8, dimension(2*numCentrals,4) :: CentralSym
!auxillaries
integer                            :: i

do i = 1,numCentrals
    CentralSym(2*i-1,:) = centralOri(i,:)
    CentralSym(2*i,1)   = CentralSym(2*i-1,1)
    CentralSym(2*i,2)   = CentralSym(2*i-1,3)
    CentralSym(2*i,3)   = CentralSym(2*i-1,2)
    CentralSym(2*i,4)   = CentralSym(2*i-1,4)
end do

end function

end program
