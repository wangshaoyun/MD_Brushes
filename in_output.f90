module input_output
implicit none 

save

	real*8,  allocatable, dimension(:,:) :: phi_tot
	real*8,  allocatable, dimension(:,:) :: phi_l
	real*8,  allocatable, dimension(:,:) :: phi_le
	real*8,  allocatable, dimension(:,:) :: phi_s
	real*8,  allocatable, dimension(:,:) :: phi_sb
	real*8,  allocatable, dimension(:,:) :: phi_se
	real*8,  allocatable, dimension(:,:) :: phi_a
	real*8,  allocatable, dimension(:,:) :: phi_i
	real*8,  allocatable, dimension(:,:) :: phi_q
	real*8,  allocatable, dimension(:,:) :: delta_angle1
	real*8,  allocatable, dimension(:,:) :: delta_angle2
	real*8,  allocatable, dimension(:,:) :: delta_angle3
	real*8,  allocatable, dimension(:,:) :: theta_l
	real*8,  allocatable, dimension(:,:) :: theta_lz
	real*8,  allocatable, dimension(:,:) :: theta_ssl
	real*8,  allocatable, dimension(:,:) :: theta_sslz
	real*8,  allocatable, dimension(:,:) :: theta_sbl
	real*8,  allocatable, dimension(:,:) :: theta_sblz
	real*8,  allocatable, dimension(:,:) :: theta_bez
	real*8,  allocatable, dimension(:,:) :: force_l
	real*8,  allocatable, dimension(:,:) :: force_sy
	real*8,  allocatable, dimension(:,:) :: force_sn
	real*8,  allocatable, dimension(:,:) :: force_so
	real*8,  allocatable, dimension(:,:) :: force_l1
	real*8,  allocatable, dimension(:,:) :: force_sy1
	real*8,  allocatable, dimension(:,:) :: force_sn1
	real*8,  allocatable, dimension(:,:) :: force_so1
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

contains

	subroutine initialize_parameters

	end subroutine initialize_parameters
	

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


end module input_output
