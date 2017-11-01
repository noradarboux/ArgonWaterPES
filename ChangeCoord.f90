module ChangeCoord
contains

!__________________________________________________________

function IntToSph(internal)
implicit none
!***Assuming our standard water position already.*** This function is solved by Mathematica
!input
!internal coordinate (ArO,ArH1,ArH2)
real*8, dimension(3), intent(in) :: internal
!output
!spherical coordinate
real*8, dimension(3)             :: IntToSph !(theta,phi,r)
!Auxillaries
real*8, dimension(3)             :: cart !cartesian coordinates
real*8                           :: O, dH1, dH2, x, y, z

O = internal(1)
dH1 = internal(2)
dH2 = internal(3)

cart(1) = 0.84728 - 0.426481*dH1**2 - 0.426481*dH2**2 + 0.852961*O**2
cart(2) = -0.330157*dH1**2 + 0.330157*dH2**2
cart(3) = dsqrt(dsqrt((-0.611779 + 0.667155*dH1** 2 - 0.29089*dH1**4 + &
0.667155*dH2**2 - 0.145764*dH1**2*dH2**2 - 0.29089*dH2**4 - 0.334309*O**2 + &
0.727543*dH1**2*O**2 + 0.727543*dH2**2*O**2 - 0.727543*O**4)**2))

x = cart(1)
y = cart(2)
z = cart(3)

IntToSph(3) = dsqrt(x**2 + y**2 + z**2) !r
IntToSph(1) = atand(dsqrt(x**2 + y**2)/z) !theta
IntToSph(2) = atand(y/x) !phi
if (x<0) then
    intToSph(2) = 180 + intToSph(2)
else if (y<0 .and. x>=0) then
    intToSph(2) = 360 + intToSph(2)
end if


end function 
!______________________________________________
function SphToCart(sph)
implicit none
!input 
!sph: spherical coordinates in the order (theta,phi,r), (deg, deg, au)
real*8, dimension(3), intent(in)  :: sph
!output
!The cartesian products in x,y,z, (au, au, au)
!The origin is the CM of water
real*8, dimension(3)              :: SphToCart
!auxillariess
real*8                            :: r,theta,phi

r = sph(3)
theta = sph(1)
phi = sph(2)

SphToCart(1) = r * sind(theta) * cosd(phi)
SphToCart(2) = r * sind(theta) * sind(phi)
SphToCart(3) = r * cosd(theta)

end function
!___________________________________________
function CartToInt(cart,whereIsWater)
implicit none
!input
!cartesian coordinates of Ar
real*8, dimension(3), intent(in)   :: cart
!cartesian coordinates of water. Arranged like {{Ox,Oy,Oz},{H1x,H1y,H1z},{H2x,H2y,H2z}}.
real*8, dimension(3,3), intent(in) :: whereIsWater
!output
!internal coordinates of Ar, arranged like {Ar-O,Ar-H1,Ar-H2}
real*8, dimension(3) :: CartToInt
!auxillaries
integer              :: i

do i = 1,3
CartToInt(i)=dsqrt((cart(1)-whereIsWater(i,1))**2+(cart(2)-whereIsWater(i,2))**2+(cart(3)-whereIsWater(i,3))**2)
end do
end function
!_________________________________________
function SphToInt(sph,whereIsWater)
implicit none
!input
!sph: spherical coordinates in the order (theta,phi,r), (deg, deg, au)
real*8, dimension(3), intent(in)   :: sph
!cartesian coordinates of water. Arranged like {{Ox,Oy,Oz},{H1x,H1y,H1z},{H2x,H2y,H2z}}.
real*8, dimension(3,3), intent(in) :: whereIsWater
!Output
!internal coordinates of Ar, arranged like {Ar-O,Ar-H1,Ar-H2}
real*8, dimension(3)               :: SphToInt

SphToInt = CartToInt(SphToCart(sph),whereIsWater)

end function

end module ChangeCoord 
