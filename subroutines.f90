module subroutines
implicit none 
contains
	
	 
subroutine rij_and_rr(rij, rsqr, i, j)
	use global_variables
	!-----------------------------------------!
	!compute displacement vector and displacement of two particles
	!input : post(pos or pos1), i, j(particle number) 
	!output: rij(displacement vecter), rr(square of displacement)
	!External Variant: Lz(used in period condition)
	!note: including period condition
	!-----------------------------------------!
	implicit none

	real*8, dimension(3), intent(out) :: rij
	real*8, intent(out) :: rsqr
	integer, intent(in) :: i
	integer, intent(in) :: j

	rij=pos(i,1:3)-pos(j,1:3)
	if (rij(1)>Lx/2) then
		rij(1)=rij(1)-Lx
	elseif(rij(1)<=-Lx/2) then
		rij(1)=rij(1)+Lx
	end if
	if (rij(2)>Ly/2) then
		rij(2)=rij(2)-Ly
	elseif(rij(2)<=-Ly/2) then
		rij(2)=rij(2)+Ly
	end if
! 	rij(1)=rij(1)-floor(rij(1)/Lx+0.5)*Lx
! 	rij(2)=rij(2)-floor(rij(2)/Ly+0.5)*Ly
	rsqr=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
end subroutine rij_and_rr










end module subroutines
