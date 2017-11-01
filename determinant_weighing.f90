program dete_weigh
implicit none
real*8, dimension(24)   :: axis
real*8, dimension(5)    :: data_geom, k_xi
real*8, dimension(24,5) :: weight
real*8, dimension(5,5) :: K_mat
real*8                 :: gamma3
real*8, dimension(5,5) :: substitution_mat
real*8                 :: this_geom, k_dete, this_dete, this_k_xj, FindDet, M55DET
integer                :: i, j

print *,'enter gamma3'
read *,gamma3

!establish the x axis
do i = 1,15
    axis(i) = 0.2 * i
end do 
do i = 1,4
    axis(15+i) = 0.5 * i + 3.0
end do
do i = 1,5
    axis(19+i) = i + 5.0
end do
!end establishing the axis

!the data geometry
data_geom = (/2.0, 1.4, 5.0, 1.0, 4.0/)
!end

!the K_matrix
do i = 1,5
    do j = 1,5
        call similarity(gamma3,data_geom(i),data_geom(j),K_mat(i,j))
    end do
end do
!end

!determinant of the K_matrix
k_dete = M55DET(K_mat)

!for each geometry on the axis
do i = 1,24
    this_geom = axis(i)
    !for the k_xi vector
    call form_k_xi(gamma3,this_geom,data_geom,k_xi)
    !for each row
    !form the substitution matrix
    do j = 1,5
        !copy the K_mat first
        substitution_mat(:,:) = K_mat
        !substitute the jth row with k_xi
        substitution_mat(j,:) = k_xi
        !find the determinant
        this_dete = M55DET(substitution_mat)
        !the similarity between this and the jth geometry
        call similarity(gamma3,this_geom,data_geom(j),this_k_xj) 
        weight(i,j) = this_dete/(k_dete)
    end do
end do
!end

!print
do i = 1,24
    write(*,"(5e20.10)") weight(i,:)
end do

end program
!_______________________________________________________________
subroutine similarity(in_gamma3,in_x1,in_x2,out_similarity)
implicit none
real*8, intent(in)  :: in_gamma3, in_x1, in_x2
real*8, intent(out) :: out_similarity

out_similarity = exp(-in_gamma3*(in_x1-in_x2)*(in_x1-in_x2))

end subroutine
!_______________________________________________________________
subroutine form_k_xi(in_gamma3,in_this_geom, in_data_geom, out_k_xi)
implicit none
real*8, intent(in) :: in_gamma3,in_this_geom
real*8, dimension(5), intent(in) :: in_data_geom
real*8, dimension(5), intent(out) :: out_k_xi
integer                           :: i

do i = 1,5
    call similarity(in_gamma3,in_this_geom,in_data_geom(i),out_k_xi(i))
end do

end subroutine
!_______________________________________________________________
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using
!this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of
!the diagonal elements
!
REAL FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
    
END FUNCTION FindDet
!____________________________________________________________________________
FUNCTION M55DET (A) RESULT (DET)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(5,5), INTENT(IN)  :: A

      DOUBLE PRECISION :: DET, A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55


      A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
      A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
      A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
      A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
      A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

      DET = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+    &
         A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
         A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
         A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
         A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
         A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
         A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
         A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
         A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
         A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
         A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
         A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
         A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
         A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
         A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
         A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
         A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
         A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
         A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
         A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
         A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
         A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
         A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
         A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
         A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
         A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
         A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
         A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
         A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
         A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
         A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
         A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
         A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
         A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
         A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
         A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
         A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
         A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
         A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
         A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

      RETURN

      END FUNCTION M55DET

