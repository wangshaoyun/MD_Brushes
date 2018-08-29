module initialize_update
implicit none
contains

subroutine Initialize_position
	!-----------------------------------------!
	!input: pos
	!output: pos
	!External variants: R_bond, R_box, Lz, Npe, qq
	!note: the box's radius is (R_rod,R_box) and height is (-Lz/2, Lz/2)
	!-----------------------------------------!
	implicit none
	integer :: i, j, k, l, m, n, x, y, p
	real*8 :: theta, rnd1, rnd2, rnd3, rsqr
	real*8, dimension(3) :: rij

	!position of PE, random grafted and the chains are also random
	!对每个粒子赋值时，都要用while循环，判断其赋值的合理性
	l=0
	!star brushes
	write(*,*) 'star brushes'
	do i=1, Nga
		m=1
		l=l+1
		do while (m==1)			!anchor point
			m=0
			call random_number(rnd1)
			call random_number(rnd2)
			pos(l,1)=rnd1*Lx-Lx/2
			pos(l,2)=rnd2*Ly-Ly/2
			pos(l,3)=0
			do n=1,l-1
				call rij_and_rr(rij,rsqr,n,l)
				if (rsqr<0.8) then
					m=1
					cycle
				end if
			end do
		end do
		do j=1, arm
			do k=1, Nma
				l=l+1
				m=1
				p=0
				do while (m/=0)
					m=0
					call random_number(rnd1)
					call random_number(rnd2)
					if (k==1 .and. j/=1) then		!branching point
						pos(l,1)=pos(l-(j-2)*Nma-1,1)+R_bond*cos(2*pi*rnd2)*sin(pi*rnd1)
						pos(l,2)=pos(l-(j-2)*Nma-1,2)+R_bond*sin(2*pi*rnd2)*sin(pi*rnd1)
						pos(l,3)=pos(l-(j-2)*Nma-1,3)+R_bond*cos(pi*rnd1)
					else
						if(p<10) then
							rnd1=rnd1/1000
						end if
						pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
						pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
						pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
					end if
						!periodic condition
					if (pos(l,1)>Lx/2) then
						pos(l,1)=pos(l,1)-Lx
					elseif(pos(l,1)<=-Lx/2) then
							pos(l,1)=pos(l,1)+Lx
					end if
					if (pos(l,2)>Ly/2) then
						pos(l,2)=pos(l,2)-Ly
					elseif(pos(l,2)<=-Ly/2) then
						pos(l,2)=pos(l,2)+Ly
					end if
					do n=1,l-1
						call rij_and_rr(rij,rsqr,n,l)
						if (rsqr<0.7 .or. pos(l,3)<0.9) then
							m=1
							p=p+1
							cycle
						end if
					end do
				end do
!						write(*,*) l,j,k,pos(l,:)
    	end do
		end do
	end do

	!linear brushes
	write(*,*) 'linear brushes'
	do i=Nga+1, N_anchor
		m=1
		l=l+1
		do while (m==1)
			m=0
			call random_number(rnd1)
			call random_number(rnd2)
			pos(l,1)=rnd1*Lx+Lx/2
			pos(l,2)=rnd2*Ly+Ly/2
			pos(l,3)=0
			do n=1,l-1
				call rij_and_rr(rij,rsqr,n,l)
				if (rsqr<0.8) then
					m=1
					cycle
				end if
			end do
		end do
		do k=2, Nml
			l=l+1
			m=1
			p=0
			do while (m==1)
				m=0
				call random_number(rnd1)
! 					rnd1=rnd1**2
				if (p<10) then
					rnd1=rnd1/2
				else
					rnd1=rnd1**2
				end if
				call random_number(rnd2)
				pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
				pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
				pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
				if (pos(l,1)>Lx/2) then
					pos(l,1)=pos(l,1)-Lx
				elseif(pos(l,1)<=-Lx/2) then
					pos(l,1)=pos(l,1)+Lx
				end if
				if (pos(l,2)>Ly/2) then
					pos(l,2)=pos(l,2)-Ly
				elseif(pos(l,2)<=-Ly/2) then
					pos(l,2)=pos(l,2)+Ly
				end if
				do n=1,l-1
					call rij_and_rr(rij,rsqr,n,l)
					if (rsqr<0.7 .or. pos(l,3)<0.9) then
						m=1
						p=p+1
						cycle
					end if
				end do
			end do
! 				write(*,*) l,k,pos(l,:)
		end do
	end do


! 	!position of PE, random grafted and the chains are also random
! 	!对每个粒子赋值时，都要用while循环，判断其赋值的合理性
! 	l=0
! 	!star brushes
! 	write(*,*) 'star brushes'
! 	do i=1, N_anchor
! 		x=(i-1)/nint(sqrt(N_anchor*1.))+1						!the root of N_anchor must be integer
! 		y=mod(i-1,nint(sqrt(N_anchor*1.)))+1
! 		if (mod(x+y,2)==0) then
! 			m=1
! 			l=l+1
! 			do while (m==1)			!anchor point
! 				m=0
! 				pos(l,1)=Lx/nint(sqrt(N_anchor*1.))*(x-0.5)-Lx/2
! 				pos(l,2)=Ly/nint(sqrt(N_anchor*1.))*(y-0.5)-Ly/2
! 				pos(l,3)=0
! ! 				call random_number(rnd1)
! ! 				call random_number(rnd2)
! ! 				pos(l,1)=rnd1*Lx-Lx/2
! ! 				pos(l,2)=rnd2*Ly-Ly/2
! ! 				pos(l,3)=0
! ! 				do n=1,l-1
! ! 					call rij_and_rr(rij,rsqr,n,l)
! ! 					if (rsqr<0.8) then
! ! 						m=1
! ! 						cycle
! ! 					end if
! ! 				end do
! 			end do
! 			do j=1, arm
! 				do k=1, Nma
! 					l=l+1
! 					m=1
! 					p=0
! 					do while (m/=0)
! 						m=0
! 						call random_number(rnd1)
! 						call random_number(rnd2)
! 						if (k==1 .and. j/=1) then		!branching point
! 							pos(l,1)=pos(l-(j-2)*Nma-1,1)+R_bond*cos(2*pi*rnd2)*sin(pi*rnd1)
! 							pos(l,2)=pos(l-(j-2)*Nma-1,2)+R_bond*sin(2*pi*rnd2)*sin(pi*rnd1)
! 							pos(l,3)=pos(l-(j-2)*Nma-1,3)+R_bond*cos(pi*rnd1)
! 						else
! 							if(p<10) then
! 								rnd1=rnd1/1000
! 							end if
! 							pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
! 							pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
! 							pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
! 						end if
! 						!periodic condition
! 						if (pos(l,1)>Lx/2) then
! 							pos(l,1)=pos(l,1)-Lx
! 						elseif(pos(l,1)<=-Lx/2) then
! 							pos(l,1)=pos(l,1)+Lx
! 						end if
! 						if (pos(l,2)>Ly/2) then
! 							pos(l,2)=pos(l,2)-Ly
! 						elseif(pos(l,2)<=-Ly/2) then
! 							pos(l,2)=pos(l,2)+Ly
! 						end if
! 						do n=1,l-1
! 							call rij_and_rr(rij,rsqr,n,l)
! 							if (rsqr<0.7 .or. pos(l,3)<0.9) then
! 								m=1
! 								p=p+1
! 								cycle
! 							end if
! 						end do
! 					end do
! ! 					write(*,*) l,j,k,pos(l,:)
! 				end do
! 			end do
! 		end if
! 	end do
	
! 	!linear brushes
! 	write(*,*) 'linear brushes'
! 	do i=1, N_anchor
! 		x=(i-1)/nint(sqrt(N_anchor*1.))+1
! 		y=mod(i-1,nint(sqrt(N_anchor*1.)))+1
! 		if (mod(x+y,2)==1) then
! 			m=1
! 			l=l+1
! 			do while (m==1)
! 				m=0
! 				pos(l,1)=Lx/nint(sqrt(N_anchor*1.))*(x-0.5)-Lx/2
! 				pos(l,2)=Ly/nint(sqrt(N_anchor*1.))*(y-0.5)-Ly/2
! 				pos(l,3)=0
! ! 				call random_number(rnd1)
! ! 				call random_number(rnd2)
! ! 				pos(l,1)=rnd1*Lx-Lx/2
! ! 				pos(l,2)=rnd2*Ly-Ly/2
! ! 				pos(l,3)=0
! ! 				do n=1,l-1
! ! 					call rij_and_rr(rij,rsqr,n,l)
! ! 					if (rsqr<0.7) then
! ! 						m=1
! ! 						cycle
! ! 					end if
! ! 				end do
! 			end do
! 			do k=2, Nml
! 				l=l+1
! 				m=1
! 				p=0
! 				do while (m==1)
! 					m=0
! 					call random_number(rnd1)
! ! 					rnd1=rnd1**2
! 					if (p<10) then
! 						rnd1=rnd1/5
! 					else
! 						rnd1=rnd1**2
! 					end if
! 					call random_number(rnd2)
! 					pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
! 					pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
! 					pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
! 					if (pos(l,1)>Lx/2) then
! 						pos(l,1)=pos(l,1)-Lx
! 					elseif(pos(l,1)<=-Lx/2) then
! 						pos(l,1)=pos(l,1)+Lx
! 					end if
! 					if (pos(l,2)>Ly/2) then
! 						pos(l,2)=pos(l,2)-Ly
! 					elseif(pos(l,2)<=-Ly/2) then
! 						pos(l,2)=pos(l,2)+Ly
! 					end if
! 					do n=1,l-1
! 						call rij_and_rr(rij,rsqr,n,l)
! 						if (rsqr<0.7 .or. pos(l,3)<0.9) then
! 							m=1
! 							p=p+1
! 							cycle
! 						end if
! 					end do
! 				end do
! ! 				write(*,*) l,k,pos(l,:)
! 			end do
! 		end if
! 	end do
	
	write(*,*) 'ions'
	!position of anions of PE
	do i=Npe+1, NN
		m=1
		do while (m==1)
			m=0
			call random_number(rnd1)
			call random_number(rnd2)
			call random_number(rnd3)
			pos(i,1)=rnd1*Lx-Lx/2
			pos(i,2)=rnd2*Ly-Ly/2
			pos(i,3)=rnd3*(Lz-1.8)+0.9
			do j=1,i-1
				call rij_and_rr(rij,rsqr,i,j)
				if (rsqr<0.8) then
					m=1
					cycle
				end if
			end do
		end do
		pos(i,4)=-qq/abs(qq)
	end do
	
	do i=1,Nq/(nint(abs(qq))+1)
		pos(charge(i),4)=qq									!end charged
	end do
end subroutine Initialize_Position


subroutine initialize_velocity
!-----------------------------------------!
!input: pos
!output: pos
!External variants: R_bond, R_box, Lz, Npe, qq
!note: the box's radius is (R_rod,R_box) and height is (-Lz/2, Lz/2)
!-----------------------------------------!
implicit none
integer i,j,k
real*8 rnd

! generate gauss velocity distribution
do i=1, NN
  do j=1, 3
    call gauss_dist(0.D0,sqrt(1./Beta),rnd)
    vel(i,j)=rnd
  end do
end do
do i=1,N_anchor
	vel(anchor_list(i),:)=0
end do
call rescale_velocity
end subroutine initialize_velocity


subroutine rescale_velocity
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8 v2_sum     !sum of velocity square
	integer i,j
	real*8, dimension(3) :: v_sum
	
	v_sum=0.
	v2_sum=0.
	do i=1, NN
	  do j=1, 3
	    v_sum(j)=v_sum(j)+vel(i,j)
	  end do
	end do
	v_sum=v_sum/(NN-N_anchor);
	do i=1,NN
	  do j=1,3
	    vel(i,j)=vel(i,j)-v_sum(j)
	  end do
		v2_sum=v2_sum+vel(i,1)*vel(i,1)+vel(i,2)*vel(i,2)+vel(i,3)*vel(i,3)
	end do
	do i=1,N_anchor
		v2_sum=v2_sum-dot_product(vel(anchor_list(i),:),vel(anchor_list(i),:))
		vel(anchor_list(i),:)=0
	end do
  ! compute mean square velocity
  v2_sum=v2_sum/(NN-N_anchor-1)
  ! rescale_velocity
	vel=vel*sqrt(3./Beta/v2_sum)
end subroutine rescale_velocity


subroutine new_position
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, dimension(3) :: dri
	real*8 dr1,dr2,dmax1,dmax2,eta1,eta2,eta3
	integer i

	dmax1=0.
	dmax2=0.
	vel=acc*dt/2+(1-xi*dt/2)*vel
	do i=1,NN
		dri=vel(i,:)*dt
  	pos(i,1:3)=pos(i,1:3)+dri
		dr1=sqrt(dot_product(dri,dri))
		if (dmax1<=dr1) then
			dmax1=dr1
		end if
		if (pos(i,4) /= 0) then
			dr2=dr1
			if (dmax2<=dr2) then
				dmax2=dr2
			end if
		end if
	end do
	dr_max1=dr_max1+dmax1
	dr_max2=dr_max2+dmax2
  call period_condition
  call compute_force
	vel=(vel+acc*dt/2.)/(1+xi*dt/2.)
end subroutine new_position


end module update_position


