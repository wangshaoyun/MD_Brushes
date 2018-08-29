module global_variables
implicit none
save
!########################constants#########################!
	real*8, parameter:: pi=3.141592653589793D0			
	real*8, parameter:: gamma=.5772156649015329D0		!Euler Gamma
!########################constants#########################!

!####################systems coefficient###################!
	integer :: Npe			!Total monomers in Polyelectrolytes(PE)
	integer :: arm			!Arms of star brushes
											!including the chains anchored to the plate
	integer :: Nma			!Monomers of each arm
	integer :: Nml			!Monomers of each linear chain
	integer :: Nga			!Number of star chains grafted on plate
	integer :: Ngl			!Number of linear chains grafted on plate
	integer :: Nta			!Total monomers of star brushes
	integer :: Ntl			!Total monomers of linear brushes
	integer :: Nq				!Total charge in the system
	integer :: NN				!Total particles in the system
	integer :: N_anchor !Anchored chains
	integer :: man		  !Manning effect: every man pariticle have one charge
	real*8  :: Lx			  !Length of cell in x direction
	real*8  :: Ly       !Length of cell in y direction
	real*8  :: Lz			  !Distance of two plate
	real*8  :: Z_empty  !Empty space ratio of height and length in slab geometry
	real*8  :: ratio_xy !Rotio of length x and width y of the box
	real*8  :: sigmag	  !Grafting density of brushes on the plate
	real*8  :: Beta		  !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
	real*8  :: qq			  !Charge of charged monomers
!##################end systems coefficient#################!

!##################running and Histogram###################!
	integer :: restart_or_continue  !Restart or continue after breaking off
	integer :: StepNum0							!Steps of preheating
	integer :: StepNum							!Steps of running
	integer :: DeltaStep				    !steps of each calculation of physical 
																	!quantities
	integer :: step									!Ordering number
	integer :: dstep								!Ordering number
	integer :: multistep						!each multistep recalculate the coulomb force
	real*8  :: dt										!time of each move
	!
	!timing
	real*8  :: started							!time at starting
	real*8  :: finished							!time at finishing
	real*8  :: total_time=0					!total time of the simulation
	!
	!histogram
	integer :: SizeHist=500					!number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
	real*8, allocatable, dimension(:,:) :: pos		!array of position
	real*8, allocatable, dimension(:,:) :: vel		!array of velocity
	real*8, allocatable, dimension(:,:) :: acc		!array of accelaration
!########################end arrays########################!


contains
	


!##############charge_function#############!
subroutine charge_function
	!----------------------------------------!
	!note: mod(Nma,man) and mod(Nml-1,man) must be zero
	!
	!----------------------------------------!
	implicit none
	integer i,j,k,l

	l=0
	!star brushes
	if (qq/=0) then
		do i=1,Nga
			do k=2,Nma+1     !the chain anchored to the plate
				if (mod(k-1,man)==0) then
					l=l+1
					charge(l)=(i-1)*(arm*Nma+1)+k
				end if
			end do
			do j=2,arm	 
				do k=1,Nma
					if (mod(k,man)==0) then
						l=l+1
						charge(l)=(i-1)*(arm*Nma+1)+((j-1)*Nma+1)+k
					end if
				end do
			end do
		end do
		!linear brushes
		do i=1,Ngl
			do j=2,Nml
				if (mod(j-1,man)==0) then
					l=l+1
					charge(l)=Nta+(i-1)*Nml+j
				end if
			end do
		end do
		do i=Npe+1,NN
			l=l+1
			charge(l)=i
		end do
	end if
! 	write(*,*) charge(1:100)
! 	stop

	!anchors of the polymers 
	do i=1,Nga
		anchor_list(i)=(i-1)*(arm*Nma+1)+1
	end do
	do i=1,Ngl
		anchor_list(i+Nga)=Nta+(i-1)*Nml+1
	end do

	!pairs of monomers of the FENE force
	l=0
	do i=1,Nga
		do k=1,Nma  		!the anchor to the branching point
			l=l+1
			fene_list(l,1)=(i-1)*(arm*Nma+1)+k
			fene_list(l,2)=(i-1)*(arm*Nma+1)+k+1
		end do
		do j=2,arm		  !the arms
			l=l+1				  !the branching point to the first point of the arms
			fene_list(l,1)=(i-1)*(arm*Nma+1)+Nma+1
			fene_list(l,2)=(i-1)*(arm*Nma+1)+(Nma*(j-1)+1)+1
			do k=1,Nma-1	!the first point of the arms to the last point of the arms
				l=l+1
				fene_list(l,1)=(i-1)*(arm*Nma+1)+(Nma*(j-1)+1)+k
				fene_list(l,2)=(i-1)*(arm*Nma+1)+(Nma*(j-1)+1)+k+1
			end do
		end do
	end do
		
	do i=1,Ngl
		do j=1,Nml-1
			l=l+1
			fene_list(l,1)=Nta+(i-1)*Nml+j
			fene_list(l,2)=Nta+(i-1)*Nml+j+1
		end do
	end do
end subroutine charge_function
!##########################################!


!#############fourier_function#############!
subroutine fourier_function
	!----------------------------------------!
	!
	!----------------------------------------!
	implicit none
	integer i,j,k,l
	real*8 ksqr,k1,k2,k3,factor,kcut
	
	K_total=0
	do k=0,Kmax3
		do i=-Kmax1,Kmax1
			do j=-Kmax2,Kmax2
				kcut=(1.*i/Kmax1)*(1.*i/Kmax1)+(1.*j/Kmax2)*(1.*j/Kmax2)+(1.*k/Kmax3)*(1.*k/Kmax3)
				if ( kcut>1 .or. kcut==0) cycle
				K_total=K_total+1
			end do
		end do
	end do

	allocate(exp_ksqr(K_total))
	l=0
	do k=0,Kmax3
		if (k==0) then
			factor=1
		else
			factor=2
		end if
		do i=-Kmax1,Kmax1
			do j=-Kmax2,Kmax2
				kcut=(1.*i/Kmax1)*(1.*i/Kmax1)+(1.*j/Kmax2)*(1.*j/Kmax2)+(1.*k/Kmax3)*(1.*k/Kmax3)
				if ( kcut>1 .or. kcut==0) cycle
				k1=2*pi*i/Lx
				k2=2*pi*j/Ly
				k3=2*pi*k/Lz/Z_empty
				ksqr=k1*k1+k2*k2+k3*k3
				l=l+1
				exp_ksqr(l)=factor*4*pi/(Lx*Ly*Lz*Z_empty)*exp(-ksqr/4/alpha2)/ksqr*lb/Beta
			end do
		end do
	end do
end subroutine fourier_function
!##########################################!

!###############real_function##############!
subroutine real_function
	!----------------------------------------!
	!
	!----------------------------------------!
	implicit none
	real*8 :: rr,rsqr,c1,c2,rf
	integer :: N=100000,i,j
	real_fun=0.
	c1=lb/Beta
	c2=2*alpha/sqrt(pi)
	real_itv=(tol+0.1)**2/N/alpha2
	do i=1,N
		rsqr=i*real_itv
		rr=sqrt(rsqr)
		rsqr=rr*rr
		real_fun(i)=c1/rsqr*( erfc(alpha*rr)/rr + c2*exp(-alpha2*rsqr) )
	end do
end subroutine real_function
!##########################################!


!###############bspln_coeffs###############!
subroutine bspln_coeffs
	!----------------------------------------!
	!2010, Gradimir V. Milovanovic, Applied Mathematics Letters
	!----------------------------------------!
	implicit none
	integer :: h,i,j,k,l,m,g,ii
	real*8, allocatable, dimension(:,:,:) :: a
! 	real*8 :: factorial
	bspln_cof=0
	dspln_cof=0
	
	do h=1,3
		allocate(a(ordr(h),ordr(h),ordr(h)+1))
		a=0
		a(1,1,2)=1
		do m=2,ordr(h)
			a(m-1+1,m,-1+2)=0
			a(m-2+1,m,-1+2)=0
			g=floor(m/2.D0)-1
			do k=0,g
				a(m-1+1,m,k+2) = 1.D0/(m-1)/(m-k)* ( 1.D0*m*a(m-2+1,m-1,k+2) - 1.D0*k*(m-1)*a(m-1+1,m,k-1+2) &
												 - 1.D0*a(m-2+1,m,k-1+2) )
				a(m-1+1,m,m-k-1+2) = (-1.D0)**(m-1) * a(m-1+1,m,k+2)
				do i=m-2,0,-1
					a(i+1,m,k+2) = 1.D0*m/(i+1-m) * ( 1.D0*(i+1)*a(i+1+1,m,k+2)-1.D0*a(i+1,m-1,k+2) )
					do ii=0,m-i-1
						a(i+1,m,m-k-1+2) = 1.D0*a(i+1,m,m-k-1+2) + factorial(ii+i)/factorial(ii)  &
															 * a(i+ii+1,m,k+2)*m**ii
					end do
					a(i+1,m,m-k-1+2) = 1.D0*a(i+1,m,m-k-1+2) * (-1)**i/factorial(i)
				end do
			end do
			if (mod(m,2)/=0) then
				a(m-1+1,m,g+1+2) = 1.D0/(m-1)/(m-(g+1)) * ( 1.D0*m*a(m-2+1,m-1,g+1+2) - 1.D0*(g+1)*(m-1) &
													 * a(m-1+1,m,(g+1)-1+2) - 1.D0*a(m-2+1,m,(g+1)-1+2) )
				do i=m-2,0,-1
					a(i+1,m,g+1+2) = 1.D0*m/(i+1-m) * ( 1.D0*(i+1)*a(i+1+1,m,g+1+2)-1.D0*a(i+1,m-1,g+1+2) )
				end do
			end if
		end do

		do i=0,ordr(h)-1
			do k=0,ordr(h)-1
				bspln_cof(k+1,i+1,h)=a(ordr(h)-i,ordr(h),k+2) 
			end do
		end do
		do i=1,ordr(h)
			do j=1,ordr(h)
				dspln_cof(i,j,h)=bspln_cof(i,j,h)*(ordr(h)-j)
			end do
		end do
		deallocate(a)
	end do
! 	write(*,*) bspln_cof(:,:,1)
! 	stop
end subroutine bspln_coeffs
!##########################################!


!##################PME_BC##################!
subroutine PME_BC
	!----------------------------------------!
	!
	!----------------------------------------!
	implicit none
	complex (kind=8) :: bx(ng(1)), by(ng(2)), bz(ng(3))
	real*8 :: mx(ng(1)), my(ng(2)), mz(ng(3))
	real*8 :: ex(ng(1)), ey(ng(2)), ez(ng(3))
	real*8 dx, dxy, dxyz
	integer :: i,j,k
	
	BC_PME=0
	
	call pmeOrthoTabBC(bx, mx, ex, 1)
	call pmeOrthoTabBC(by, my, ey, 2)
	call pmeOrthoTabBC(bz, mz, ez, 3)
	
  do i=1,ng(1)
    dx = mx(i)*mx(i)
    do j=1,ng(2)
      dxy = dx + my(j)*my(j)
      do k = 1,ng(3)
        dxyz = dxy + mz(k)*mz(k)
				if (i+j+k>3) then													!避免除以0
        	BC_PME(i,j,k) = bx(i)*conjg(bx(i))*by(j)*conjg(by(j))*bz(k)*conjg(bz(k))*ex(i)*ey(j)*ez(k)/dxyz
				end if
      end do
    end do
  end do
end subroutine PME_BC
!##########################################!


!###############pmeOrthoTabBC##############!
subroutine pmeOrthoTabBC(bv, mv, ev, xyz)
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	integer, intent(in) :: xyz
	complex (kind=8), dimension(ng(xyz)), intent(out) :: bv
	real*8, dimension(ng(xyz)), intent(out) :: mv
	real*8, dimension(ng(xyz)), intent(out) :: ev
	real*8 :: tpS,pivol
	integer :: i,j,ord
	complex (kind=8) :: bsp
	ord=ordr(xyz)
	bv=0
	mv=0
	ev=0
	
  if (xyz == 1) then
    pivol = ng(1)*ng(2)*ng(3)/(pi*gdim(1)*gdim(2)*gdim(3))
  end if
	tpS=pi/alpha
  do i = 1,ng(xyz)
    if (i-1 <= ng(xyz)/2) then						!这个m的取值范围有些讲究
      mv(i) = (i-1)/gdim(xyz)
    else
      mv(i) = (i-1-ng(xyz))/gdim(xyz)
    end if
    if (xyz == 1) then
      ev(i) = pivol*exp(-tpS*tpS*mv(i)*mv(i))
    else
      ev(i) = exp(-tpS*tpS*mv(i)*mv(i))
    end if
  end do
	
	do i=0, ng(xyz)-1
		bsp=0
		do j=0, ord-2
			bsp=bsp+bspln(j+1.D0,ord)*cmplx(cos(2*pi*i*j/ng(xyz)),sin(2*pi*i*j/ng(xyz)))
		end do
		if (abs(bsp)<1e-10) then
			bsp=1e15
		end if
		bv(i+1)=cmplx(cos(2*pi*i*(ord-1)/ng(xyz)),sin(2*pi*i*(ord-1)/ng(xyz)))/bsp !complex number can be divided
	end do
	
end subroutine pmeOrthoTabBC
!###########################################!


!################factorial#################!
Function factorial(n)
	integer, intent(in) :: n
	real*8 :: factorial
	integer :: i
	
	factorial=1.D0
	if (n/=0) then
		do i=n,1,-1
			factorial=factorial*i
		end do
	end if
	
end function factorial
!##########################################!


!###################bspln###################!
recursive real*8 Function bspln(x,ord) result(y)
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, intent(in)  :: x
	integer, intent(in) :: ord

  if (ord > 2) then
    y = (x/(ord-1))*bspln(x,ord-1) + ((ord-x)/(ord-1))*bspln(x-1.0,ord-1)
  elseif (ord == 2) then
    if (x >= 0.0 .and. x <= 2.0) then
      y = 1.0 - abs(x-1.0)
    else
      y = 0.0
    end if
  else
    write(*,*) 'bspln error!  Order = ', ord
  end if

end Function bspln
!###########################################!


end module data_module
