program Asymfit
implicit none

real*8 :: r,r0,A,alpha,B,beta,mean,gamma1,gamma2
real*8, dimension(24,2) :: table
integer, dimension(2) :: table_size
real*8, dimension(1) :: key
real*8, dimension(5,2) :: closest_1
real*8, dimension(24) :: graphMean,graphVar
INTEGER :: I
gamma1 = 0.5
gamma2 = 0.5
print *, 'enter A, alpha, beta'
read(*,*) A,alpha,beta
!call asym_mean(r,r0,A,alpha,B,beta,mean)
!print *, 'the estimated energy at r is ',mean 

call read_abinitio_H2(table)
table_size = (/24,2/)

print *,'enter key between 0 and 10'
read(*,*) key(1)

!call find_closest(key,1,table,table_size,5,closest_1)
!do i=1,5
!    print *,closest_1(i,:)
!end do

call single_data_graph(key,1,A,alpha,beta,table,table_size,gamma1,gamma2,graphMean,graphVar)
end program

!__________________________________________________________________________
subroutine single_data_graph(in_r0,in_dimension,in_A,in_alpha,in_beta,in_table,in_tableSize,in_gamma1,in_gamma2,out_graphMean,out_graphVar)
implicit none
!input output
real*8                             :: in_A, in_alpha, in_beta, in_gamma1, in_gamma2
integer                            :: in_dimension
real*8, dimension(in_dimension)    :: in_r0
integer, dimension(2)              :: in_tableSize
real*8, dimension(in_tableSize(1),in_tableSize(2))    :: in_table
real*8, dimension(in_tableSize(1)) :: out_graphMean, out_graphVar
!auxis
real*8, dimension(in_dimension)       :: current_geometry, effective_r0
real*8                                :: B
real*8, dimension(1, in_tableSize(2)) :: closest
integer                               :: i, printma

!find the closest geometry to r0 in the table
call find_closest(in_r0,1,in_table,in_tableSize,1,closest)
effective_r0 = closest(1,1:in_dimension)
B            = closest(1,in_tableSize(2))
!calculate mean and variance at each geometry
do i=1,in_tableSize(1)
    current_geometry = in_table(i,1:in_dimension)
    call asym_mean(current_geometry,effective_r0,in_dimension,in_A,in_alpha,B,in_beta,out_graphMean(i))
    call asym_var(current_geometry,effective_r0,in_dimension,in_gamma1,in_gamma2,out_graphVar(i))
end do
!print
print *,'Do you want to print the result? 1 for yes, other for no'
read *,printma
if (printma==1) then
    do i=1,in_tableSize(1)
        write(*,"(4f10.6)") out_graphMean(i),out_graphVar(i)
    end do
end if

end subroutine
!__________________________________________________________________________
subroutine asym_mean(in_r,in_r0,in_dimension,in_A,in_alpha,in_B,in_beta,out_mean)
implicit none

!Input
!in_r     the arbitrary molecular geometry where we want to estimate the energy
!in_r0    the data geometry where the energy is know by ab-initio calc
!in_A     the asymmetric contributor
!in_alpha the asymmetric decay speed
!in_B     the symmetric contributor, essentially the ab-intio energy at r0
!in_beta  the symmetric decay speed
!output
!out_mean the estimated energy at geometry r.
integer                        :: in_dimension
real*8,dimension(in_dimension) :: in_r,in_r0
real*8                         :: in_A,in_alpha,in_B,in_beta,out_mean
real*8                         :: d, distance

d = in_r0(1)-in_r(1)
out_mean = in_A*d*exp(-in_alpha*d*d) + in_B*exp(-in_beta*d*d)

end subroutine
!___________________________________________________________________________
subroutine asym_var(in_r,in_r0,in_dimension,in_gamma1,in_gamma2,out_var)
implicit none

!input
!in_r     the arbitrary molecular geometry where we want to estimate the energy
!in_r0    the data geometry where the energy is know by ab-initio calc
!in_B     the symmetric contributor, essentially the ab-intio energy at r0
!in_gamma1 maximum variance allowed
!in_gamma2 increase speed with variance with geometry deviation of r from r0
!output
!out_var  the estimated variance at the geometry r
integer :: in_dimension
real*8,dimension(in_dimension) :: in_r, in_r0
real*8 :: in_gamma1, in_gamma2, out_var 
real*8 :: d, distance

d = distance(in_r,in_r0,in_dimension)
out_var = in_gamma1*(1.0-exp(-in_gamma2*d*d))

end subroutine
!__________________________________________________________________________
subroutine read_abinitio_H2(out_table)
implicit none
!output
!out_table Array 24*2. The first column is the internucleic distance r and the
!second column is the ab-intio energy...
real*8,dimension(24,2) :: out_table
integer                :: i

open(20,file='PES_of_H2.txt')
do i=1,24
    read(20,*) out_table(i,1),out_table(i,2)
    out_table(i,2) = out_table(i,2) + 1
end do
end subroutine
!__________________________________________________________________________
subroutine find_closest(in_key,in_dimension,in_table,in_table_size,in_num_slots,out_closest)
implicit none
!find the closest geometries in the data table to the geometry of interest
!input
!in_key       the geometry of interest
!in_table     table of data points
!in_num_slots how many closest geometries do we want to find?
!out_closest  the result. Including the found geometries and their energies
integer                                              :: in_dimension
real*8, dimension(in_dimension)                      :: in_key
integer, dimension(2)                                :: in_table_size
real*8, dimension(in_table_size(1),in_table_size(2)) :: in_table
integer                                              :: in_num_slots
real*8, dimension(in_num_slots,in_table_size(2))     :: out_closest

integer                          :: i, dimens, threshold_ind
real*8, dimension(in_num_slots)  :: smallest_distances
integer, dimension(in_num_slots) :: candidate_indices
real*8                           :: threshold_distance, current_distance
real*8                           :: distance

dimens = in_dimension
!initialize
do i=1,in_num_slots
    candidate_indices(i)=i
    smallest_distances(i)=distance(in_key,in_table(i,1:dimens),dimens)
end do
threshold_distance = maxval(smallest_distances)
threshold_ind      = maxloc(smallest_distances,1)
!iterate
do i=1+in_num_slots,in_table_size(1)
    current_distance = distance(in_key,in_table(i,1:dimens),dimens)
    if (current_distance<threshold_distance) then
        candidate_indices(threshold_ind) = i
        smallest_distances(threshold_ind) = current_distance
        threshold_distance = maxval(smallest_distances)
        threshold_ind = maxloc(smallest_distances,1)
    end if
end do
!clean up
do i=1,in_num_slots
    out_closest(i,:) = in_table(candidate_indices(i),:)
end do

end subroutine  
!_______________________________________________________________________
real*8 function distance(in_a,in_b,in_dimension)
implicit none
!"cartesian distances between two geometries"
!input
!in_a,in_b     two molecular geometries
!in_dimension  number of dimensions
!output
!distance  Cartesian distance
integer                         :: in_dimension
real*8, dimension(in_dimension) :: in_a, in_b

integer  :: i

distance = 0
do i=1,in_dimension
    distance = distance + (in_a(i)-in_b(i))**2
end do
distance = sqrt(distance)
return
end function












         
