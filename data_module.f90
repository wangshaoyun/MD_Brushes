module data_module
implicit none
save

!########################constants#########################!
real*8, parameter:: pi=3.141592653589793D0			
real*8, parameter:: gamma=.5772156649015329D0								!Euler gamma
!########################constants#########################!


!####################systems######################!
real*8 :: Lx					!radius of box
real*8 :: Ly
real*8 :: Lz					!length of box
real*8 :: ratio_xy  
real*8 :: sigmag			!density of brushes
real*8 :: Beta				!
real*8 :: qq					!charge of PE
real*8 :: R_bond			!length of chemical band
real*8 :: Z_empty

integer :: Npe				!number of molecular in Polyelectrolytes
integer :: arm				!number of arm in star brushes, including the chain anchored to the plate
integer :: Nma				!number of monomers of each arm
integer :: Nml				!number of monomers of each linear chain, Nml-1 must be the multiple of man
integer :: Nga				!number of chains
integer :: Ngl				!grafting number of linear chains
integer :: Nta				!total monomers of star brushes
integer :: Ntl				!total monomers of linear brushes
integer :: Nq					!number of molecular with charge
integer :: NN					!number of total moleculars
integer :: N_bond
integer :: N_anchor
integer :: man				!Manning effect: every man pariticle have one charge

integer, dimension(3) :: ordr				!Order of spline function
integer, dimension(3) :: ng					!wave number
real*8,  dimension(3) :: gdim				!length, width and height of the box
!####################systems######################!


!###################potential#####################!
real*8 :: rcl					!cut-off length of lj potential
real*8 :: rvl					!cut-off length of verlet list of lj potnetial
real*8 :: rcc					!cut-off length of coulomb force in real space
real*8 :: rcc2				!rcc2=rcc*rcc
real*8 :: rvc					!rvc=rcc+rsk
real*8 :: rsk					!skin thickness of of the cut-off sphere in real space
real*8 :: tol					!tolerance
real*8 :: tau_rf			!ratio of time in real space and fourier space
integer :: Kmax1			!max wave number of x direction in fourier space 
integer :: Kmax2			!max wave number of y direction in fourier space 
integer :: Kmax3			!max wave number of z direction in fourier space 
integer :: K_total
integer :: real_verlet!call real verlet list or not
real*8 :: dr_max1     !max displacement of lj verlet list
real*8 :: dr_max2			!max displacemnet of real verlet list
real*8 :: R0_2				!max length of FENE potential
real*8 :: kFENE				!coefficient of FENE potential
real*8 :: lb					!Bjerrum length
real*8 :: xi					!friction coefficient
real*8 :: EF					!electric field
real*8 :: alpha				!screening parameter of Gauss screening function 
real*8 :: alpha2			!alpha2=alpha*alpha
real*8 :: real_itv
!###################potential#####################!


!#############running and Histogram###############!
integer :: restart_or_continue  !restart or continue after breaking off
integer :: StepNum0			!steps of preheating
integer :: StepNum			!steps of running
integer :: DeltaStep		!steps of every calculation of physical quantities
integer :: step					!steps of calculate the physical quantities
integer :: dstep				!interval of each calculation of the physical quantities
integer :: longstep			!each longstep recalculate the coulomb force
real*8  :: dt						!time of each move
real*8 :: started				!time at starting
real*8 :: finished			!time at finishing
real*8 :: total_time=0	!total time of the simulation

integer :: npair1				!number of pairs in the sphere of radius of rvl
integer :: npair2				!number of pairs in the sphere of radius of rvc
integer :: SizeHist=500	!number of histogram which is equally divided
!#############running and Histogram###############!


!#################variables#######################!
real*8, allocatable, dimension(:,:) :: pos					!array of position
real*8, allocatable, dimension(:,:) :: posq					!array of position
real*8, allocatable, dimension(:,:) :: vel					!array of velocity
real*8, allocatable, dimension(:,:) :: acc					!array of accelaration
real*8, allocatable, dimension(:,:) :: acc_c				!array of accelaration
integer, allocatable, dimension(:,:):: lj_pair_list	!LJ potential verlet list
integer, allocatable, dimension(:,:):: real_pair_list!verlet list of coulomb force in real space
integer, allocatable, dimension(:) :: charge				!charge number to monomer number
integer, allocatable, dimension(:,:) :: fene_list		!pairs of two adjacent monomers
integer, allocatable, dimension(:) :: anchor_list		!number of the anchored monomer
real*8, allocatable, dimension(:) :: exp_ksqr	
real*8, allocatable, dimension(:) :: real_fun			
real*8, allocatable, dimension(:,:,:) :: bspln_cof
real*8, allocatable, dimension(:,:,:) :: dspln_cof
complex (kind=8), allocatable, dimension(:,:,:) :: Q_PME
complex (kind=8), allocatable, dimension(:,:,:) :: U_PME
complex (kind=8), allocatable, dimension(:,:,:) :: BC_PME

real*8, allocatable, dimension(:,:) :: phi_tot
real*8, allocatable, dimension(:,:) :: phi_l
real*8, allocatable, dimension(:,:) :: phi_le
real*8, allocatable, dimension(:,:) :: phi_s
real*8, allocatable, dimension(:,:) :: phi_sb
real*8, allocatable, dimension(:,:) :: phi_se
real*8, allocatable, dimension(:,:) :: phi_a
real*8, allocatable, dimension(:,:) :: phi_i
real*8, allocatable, dimension(:,:) :: phi_q
real*8, allocatable, dimension(:,:) :: delta_angle1
real*8, allocatable, dimension(:,:) :: delta_angle2
real*8, allocatable, dimension(:,:) :: delta_angle3
real*8, allocatable, dimension(:,:) :: theta_l
real*8, allocatable, dimension(:,:) :: theta_lz
real*8, allocatable, dimension(:,:) :: theta_ssl
real*8, allocatable, dimension(:,:) :: theta_sslz
real*8, allocatable, dimension(:,:) :: theta_sbl
real*8, allocatable, dimension(:,:) :: theta_sblz
real*8, allocatable, dimension(:,:) :: theta_bez
real*8, allocatable, dimension(:,:) :: force_l
real*8, allocatable, dimension(:,:) :: force_sy
real*8, allocatable, dimension(:,:) :: force_sn
real*8, allocatable, dimension(:,:) :: force_so
real*8, allocatable, dimension(:,:) :: force_l1
real*8, allocatable, dimension(:,:) :: force_sy1
real*8, allocatable, dimension(:,:) :: force_sn1
real*8, allocatable, dimension(:,:) :: force_so1
integer, allocatable, dimension(:,:) :: phi_zx
integer, allocatable, dimension(:,:) :: phi_xy
integer, allocatable, dimension(:,:) :: phi_yz
integer, allocatable, dimension(:,:) :: phi_lzx
integer, allocatable, dimension(:,:) :: phi_lxy
integer, allocatable, dimension(:,:) :: phi_lyz
integer, allocatable, dimension(:,:) :: phi_lezx
integer, allocatable, dimension(:,:) :: phi_lexy
integer, allocatable, dimension(:,:) :: phi_leyz
integer, allocatable, dimension(:,:) :: phi_szx
integer, allocatable, dimension(:,:) :: phi_sxy
integer, allocatable, dimension(:,:) :: phi_syz
integer, allocatable, dimension(:,:) :: phi_sbzx
integer, allocatable, dimension(:,:) :: phi_sbxy
integer, allocatable, dimension(:,:) :: phi_sbyz
integer, allocatable, dimension(:,:) :: phi_sezx
integer, allocatable, dimension(:,:) :: phi_sexy
integer, allocatable, dimension(:,:) :: phi_seyz
integer, allocatable, dimension(:,:) :: phi_azx
integer, allocatable, dimension(:,:) :: phi_axy
integer, allocatable, dimension(:,:) :: phi_ayz
integer, allocatable, dimension(:,:) :: phi_izx
integer, allocatable, dimension(:,:) :: phi_ixy
integer, allocatable, dimension(:,:) :: phi_iyz
integer, allocatable, dimension(:,:) :: phi_qzx
integer, allocatable, dimension(:,:) :: phi_qxy
integer, allocatable, dimension(:,:) :: phi_qyz

real*8 :: h_avg
real*8 :: hl_avg
real*8 :: hl_max
real*8 :: hl_end
real*8 :: hs_avg
real*8 :: hs_max
real*8 :: hs_end
real*8 :: hs_branch

real*8 :: Re_l
real*8 :: Re_lz
real*8 :: Re_s
real*8 :: Re_sz
real*8 :: Re_ss
real*8 :: Re_ssz
real*8 :: Re_sb
real*8 :: Re_sbz

real*8 :: Rg_l
real*8 :: Rg_lz
real*8 :: Rg_s
real*8 :: Rg_sz
real*8 :: Rg_ss
real*8 :: Rg_ssz
real*8 :: Rg_sb
real*8 :: Rg_sbz

real*8 :: kenetic_energy
real*8 :: ratio_stretch
real*8 :: ratio_collapse
real*8 :: ratio_other
!#################variables#######################!


contains
	
!#################read_data################!
subroutine read_data	
	implicit none
	logical alive
	open(10,file='./data.txt')
		read(10,*) Lz	
		read(10,*) sigmag						
		read(10,*) qq					
		read(10,*) arm
		read(10,*) Nma
		read(10,*) Nga
		read(10,*) Nml
		read(10,*) Ngl
		read(10,*) man
		read(10,*) tol
		read(10,*) tau_rf						
		read(10,*) EF
		read(10,*) StepNum0 					
		read(10,*) StepNum						
		read(10,*) DeltaStep 		
		read(10,*) longstep	
		read(10,*) dt
	close(10)
	
	Beta=1/1.2
	Z_empty=2
	ordr=(/6,6,6/)
	rcl=1.1224620
	rvl=1.8
	rsk=1
	R_bond=0.97
	R0_2=2.25
	kFENE=30
	lb=3
	xi=1
	
	Nta=(Nma*arm+1)*Nga
	Ntl=Nml*Ngl
	Npe=Nta+Ntl
	N_anchor=Nga+Ngl
	N_bond=arm*Nma*Nga+(Nml-1)*Ngl
	if ( man/=0 )	then
		!the anchor monomer is uncharged and the branching point is charged
		NN=Npe+(Npe-N_anchor)/man*nint(abs(qq))	
	else
		NN=Npe
	end if 
	if (abs(qq)==0) then
		Nq=0
	else
		Nq=(Npe-N_anchor)/man*(nint(abs(qq))+1)
	end if
	write(*,*) 'arm:',arm,'Nma:',Nma,'Nga:',Nga,'Nml:',Nml,'Ngl',Ngl
	write(*,*) 'total particles, NN:', NN
	write(*,*) 'total charged particles, Nq:', Nq
	write(*,*) 'total anchored particles, N_anchor:', N_anchor
	write(*,*) 'total brushes particles, Npe:', Npe
	write(*,*) 'total linear brushes particles, Ntl:', Ntl
	write(*,*) 'total star brushes particles, Nta:', Nta
	write(*,*) 'total chemical bonds, N_bond:', N_bond
	
	Lx=sqrt(N_anchor/sigmag)
	Ly=Lx
	Z_empty=(Lz+Lx*Z_empty)/Lz
	dr_max1=0
	dr_max2=0
	Inquire(file='start_time.txt',exist=alive)
	if (alive) then
		open(11,file='./start_time.txt')
			read(11,*) restart_or_continue
		close(11)
	else
		restart_or_continue=0
	end if
end subroutine read_data
!##########################################!

!##############data_operation##############!
subroutine data_operation
	alpha=(tau_rf*pi**3*Nq/(Lx*Ly*Lz)**2)**(1./6)
	alpha2=alpha*alpha
	rcc=tol/alpha
	rcc2=rcc*rcc
	rvc=rcc+rsk
	if ((int(Lx/rvc)*int(Ly/rvc)*int(Lz/rvc))>27 .and. qq/=0) then 
		Kmax1=ceiling(tol*Lx*alpha/pi)
		Kmax2=ceiling(tol*Ly*alpha/pi)
		Kmax3=ceiling(tol*Lz*Z_empty*alpha/pi)
		real_verlet=1
	else
		if (Lx>Ly) then
			rcc=Ly/2
		else
			rcc=Lx/2
		end if
 		rcc2=rcc*rcc
		Kmax1=ceiling(tol*tol/pi*Lx/rcc)
		Kmax2=ceiling(tol*tol/pi*Ly/rcc)
		Kmax3=ceiling(tol*tol/pi*Lz*Z_empty/rcc)
		real_verlet=0
	end if
	Kmax1=ceiling(Kmax1*3.0)
	Kmax2=ceiling(Kmax2*3.0)
	Kmax3=ceiling(Kmax3*3.0)
	ng=(/Kmax1,Kmax2,Kmax3/)
	gdim=(/Lx,Ly,Lz*Z_empty/)
end subroutine data_operation
!##########################################!

!##############data_allocate###############!
subroutine data_allocate
	!----------------------------------------!
	!allocate array and initialize them
	!----------------------------------------!
	implicit none

	allocate(pos(NN,4))
	allocate(vel(NN,3))
	allocate(acc(NN,3))
	allocate(acc_c(NN,3))
	allocate(lj_pair_list(NN*nint(8./3*pi*rvl**3),2))
	allocate(fene_list(N_bond,2))
	allocate(anchor_list(N_anchor))
	allocate(bspln_cof(maxval(ordr),maxval(ordr),3))
	allocate(dspln_cof(maxval(ordr),maxval(ordr),3))
	allocate(Q_PME(Kmax1,Kmax2,Kmax3))
	allocate(U_PME(Kmax1,Kmax2,Kmax3))
	allocate(BC_PME(Kmax1,Kmax2,Kmax3))
	allocate(real_fun(100000))
	if (qq/=0) then
		allocate(charge(Nq))
		allocate(posq(Nq,4))								!position of the charges
		allocate(real_pair_list(Nq*Nq,2))
		posq=0
		charge=0
		real_pair_list=0
	end if
	pos=0
	vel=0
	acc=0
	acc_c=0
	lj_pair_list=0
	fene_list=0
	anchor_list=0
	bspln_cof=0
	dspln_cof=0
	
	allocate(phi_tot(SizeHist,2))
	allocate(phi_l(SizeHist,2))
	allocate(phi_le(SizeHist,2))
	allocate(phi_s(SizeHist,2))
	allocate(phi_sb(SizeHist,2))
	allocate(phi_se(SizeHist,2))
	allocate(phi_a(SizeHist,2))
	allocate(phi_i(SizeHist,2))
	allocate(phi_q(SizeHist,2))
	allocate(delta_angle1(SizeHist,2))
	allocate(delta_angle2(SizeHist,2))
	allocate(delta_angle3(SizeHist,2))
	allocate(theta_l(SizeHist,2))
	allocate(theta_lz(SizeHist,2))
	allocate(theta_ssl(SizeHist,2))
	allocate(theta_sslz(SizeHist,2))
	allocate(theta_sbl(SizeHist,2))
	allocate(theta_sblz(SizeHist,2))
	allocate(theta_bez(SizeHist,2))
	allocate(force_l(Nml-1,2))
	allocate(force_sy(Nma*2,2))
	allocate(force_sn(Nma*2,2))
	allocate(force_so(Nma*2,2))
	allocate(force_l1(SizeHist,2))
	allocate(force_sy1(SizeHist,2))
	allocate(force_sn1(SizeHist,2))
	allocate(force_so1(SizeHist,2))
	allocate(phi_zx(SizeHist,SizeHist))
	allocate(phi_xy(SizeHist,SizeHist))
	allocate(phi_yz(SizeHist,SizeHist))
	allocate(phi_lzx(SizeHist,SizeHist))
	allocate(phi_lxy(SizeHist,SizeHist))
	allocate(phi_lyz(SizeHist,SizeHist))
	allocate(phi_lezx(SizeHist,SizeHist))
	allocate(phi_lexy(SizeHist,SizeHist))
	allocate(phi_leyz(SizeHist,SizeHist))
	allocate(phi_szx(SizeHist,SizeHist))
	allocate(phi_sxy(SizeHist,SizeHist))
	allocate(phi_syz(SizeHist,SizeHist))
	allocate(phi_sbzx(SizeHist,SizeHist))
	allocate(phi_sbxy(SizeHist,SizeHist))
	allocate(phi_sbyz(SizeHist,SizeHist))
	allocate(phi_sezx(SizeHist,SizeHist))
	allocate(phi_sexy(SizeHist,SizeHist))
	allocate(phi_seyz(SizeHist,SizeHist))
	allocate(phi_azx(SizeHist,SizeHist))
	allocate(phi_axy(SizeHist,SizeHist))
	allocate(phi_ayz(SizeHist,SizeHist))
	allocate(phi_izx(SizeHist,SizeHist))
	allocate(phi_ixy(SizeHist,SizeHist))
	allocate(phi_iyz(SizeHist,SizeHist))
	allocate(phi_qzx(SizeHist,SizeHist))
	allocate(phi_qxy(SizeHist,SizeHist))
	allocate(phi_qyz(SizeHist,SizeHist))
	
	phi_tot=0
	phi_l=0
	phi_le=0
	phi_s=0
	phi_sb=0
	phi_se=0
	phi_a=0
	phi_i=0
	phi_q=0
	delta_angle1=0
	delta_angle2=0
	delta_angle3=0
	theta_l=0
	theta_lz=0
	theta_ssl=0
	theta_sslz=0
	theta_sbl=0
	theta_sblz=0
	theta_bez=0
	force_l=0
	force_sy=0
	force_sn=0
	force_so=0
	force_l1=0
	force_sy1=0
	force_sn1=0
	force_so1=0
	phi_zx=0
	phi_xy=0
	phi_yz=0
	phi_lzx=0
	phi_lxy=0
	phi_lyz=0
	phi_lezx=0
	phi_lexy=0
	phi_leyz=0
	phi_szx=0
	phi_sxy=0
	phi_syz=0
	phi_sbzx=0
	phi_sbxy=0
	phi_sbyz=0
	phi_sezx=0
	phi_sexy=0
	phi_seyz=0
	phi_azx=0
	phi_axy=0
	phi_ayz=0
	phi_izx=0
	phi_ixy=0
	phi_iyz=0
	phi_qzx=0
	phi_qxy=0
	phi_qyz=0
	
	h_avg=0
end subroutine data_allocate
!##########################################!


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
