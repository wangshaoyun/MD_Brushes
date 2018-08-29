module compute_acc
	use data_module
	use subroutines
	implicit none
contains

	!###############error_analysis##############!
	subroutine error_analysis
		!---------------------------------------!
		!
		!---------------------------------------!
		implicit none
		real*8 f_r, f_k, q_tot, rmsf, tol1, tau_rf1, sumf1, sumf2,st,fn,tm1,tm2
		integer i,j,m
		real*8, dimension(NN,3):: acc_se
		
		q_tot=0
		do m=1,Nq
			i=charge(m)
			q_tot=q_tot+pos(i,4)*pos(i,4)
		end do 
		f_r=2*abs(qq)*sqrt(q_tot/rcc/Lx/Ly/Lz)*exp(-tol*tol)
		f_k=abs(qq)*alpha/(Lx*Ly*Lz)**(1./3)/pi*sqrt(8*q_tot/(Kmax1*Kmax2*Kmax3)**(1./3))*exp(-tol*tol)
		write(*,*) 'cut off in real space rcc',rcc
		write(*,*) 'length of box Lx', Lx
		write(*,*) 'length of box Ly', Ly
		write(*,*) 'length of box Lz', Lz*Z_empty
		write(*,*) 'max wave number in x direction Kmax1', Kmax1
		write(*,*) 'max wave number in y direction Kmax2', Kmax2
		write(*,*) 'max wave number in z direction Kmax3', Kmax3
		write(*,*) 'Total wave number', 2*Kmax1*Kmax2*Kmax3
		write(*,*) 'Ewald parameters alpha', alpha
		write(*,*) 'Total particle numbers NN', NN
		write(*,*) 'Total charges Nq', Nq
		write(*,*) 'whether the real space verlet list is used,0 or 1', real_verlet
		write(*,*) 'Standard Ewald error in real space:', f_r
		write(*,*) 'Standard Ewald error in fourier spqce', f_k
		
		tol1=tol
		tau_rf1=tau_rf
		tol=5										!the error is about 1e-6
		tau_rf=10
		call data_operation
		Kmax1=ceiling(Kmax1/3.)
		Kmax2=ceiling(Kmax2/3.)
		Kmax3=ceiling(Kmax3/3.)
		call fourier_function
		call real_function
		if (real_verlet /=0 ) then
			call real_verlet_list
		end if
		acc_c=0
		call Standard_Ewald
		call real_space
		acc_se=acc_c
		
		tol=tol1
		tau_rf=tau_rf1
		call data_operation
		call bspln_coeffs				!coefficients of bspline
		call PME_BC							!compute B*C
		if (real_verlet /=0 ) then
			call real_verlet_list
		end if
		acc_c=0
		posq(:,:)=pos(charge,:)
		call cpu_time(st)
		call SPME_Ewald
		call cpu_time(fn)
		tm1=fn-st
		write(*,*) 'time in SPME:', tm1
		call cpu_time(st)
		call real_space
		call cpu_time(fn)
		tm2=fn-st
		write(*,*) 'time in real space:', tm2
		
		sumf1=0
		sumf2=0
		do i=1,Nq
			do j=1,3
				sumf1=sumf1+(acc_se(i,j)-acc_c(i,j))**2
				sumf2=sumf2+acc_se(i,j)**2
			end do
		end do
		rmsf=sqrt(sumf1/sumf2)
		write(*,*) 'Rms force of coulomb force is:',rmsf
		acc_c=0

		open(60,position='append',file='./data/rms_force.txt')
			write(60,600) step*1., rmsf, npair1*1., npair2*1., tm1, tm2
			close(60)
		600 format(6F17.6) 
		
		deallocate(exp_ksqr)
	end subroutine
	!##########################################!


	!###############compute_force###############!
	subroutine Compute_Force
		!-----------------------------------------!
		!input: pos
		!output: force
		!including:
		!-----------------------------------------!
		implicit none
		real*8, dimension(3) :: rij
		real*8 :: ff,r2,eta1,eta2,eta3,st,fn
		integer :: i,j,k,n
		
		n=size(anchor_list,1)
		acc=0
	 	if (qq/=0 .and. mod(dstep,longstep)==0) then	  !each longstep renew the coulomb force
	 		acc_c=0
	 		posq(:,:)=pos(charge,:)
	 		call SPME_Ewald
	 		call real_space
	 	end if
	 	acc=acc_c
	 	call lj_force
	 	call fene_force
	 	do i=1, NN
	 		call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta1)
	 		call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta2)
	 		call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta3)
	 		acc(i,1)=acc(i,1)+eta1
	 		acc(i,2)=acc(i,2)+eta2
			acc(i,3)=acc(i,3)+eta3+pos(i,4)*EF 						!EF points from bottom to top
		end do
		do i=1,N_anchor
			j=anchor_list(i)
			acc(j,1:3)=0
		end do
	end subroutine Compute_Force
	!###########################################!


!#################lj_force##################!
subroutine lj_force
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	integer i,j,k,n
	real*8, dimension(3) :: rij
	real*8 :: rsqr,inv_r2,inv_r6,dl,cc,hLx,nhLx,hLy,nhLy,rcl2
	real*8, dimension(NN,3) :: acc_lj

	hLx=Lx/2
	nhLx=-Lx/2
	hLy=Ly/2
	nhLy=-Ly/2
	rcl2=rcl*rcl
	acc_lj=0
	do k=1,npair1
		i=lj_pair_list(k,1)
		j=lj_pair_list(k,2)
		rij=pos(i,1:3)-pos(j,1:3)
		if (rij(1)>hLx) then
			rij(1)=rij(1)-Lx
			elseif(rij(1)<=nhLx) then
				rij(1)=rij(1)+Lx
			end if
			if (rij(2)>hLy) then
				rij(2)=rij(2)-Ly
				elseif(rij(2)<=nhLy) then
					rij(2)=rij(2)+Ly
				end if
				rsqr=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
				if (rsqr<rcl2) then
					inv_r2=1/rsqr
					inv_r6=inv_r2*inv_r2*inv_r2
					acc_lj(i,:)=acc_lj(i,:)+inv_r2*inv_r6*(inv_r6-0.5)*rij
				end if
			end do
			acc_lj=acc_lj*48
	!wall force
	do i=1,NN
		if (abs(pos(i,3))<0.01 .or. abs(pos(i,3)-Lz)<0.01) 	then  !exclude the end monomer on the plate
			cycle
			elseif ( pos(i,3)<rcl ) then
				dl=pos(i,3)
				acc_lj(i,3)=acc_lj(i,3)+0.4*(3/(dl**11)-1/(dl**5))*dl
				elseif ( (Lz-pos(i,3))<rcl ) then
					dl=Lz-pos(i,3)
					acc_lj(i,3)=acc_lj(i,3)-0.4*(3/(dl**11)-1/(dl**5))*dl
				end if
			end do
			acc=acc+acc_lj
		end subroutine lj_force
!###########################################!


!################real_space#################!
subroutine real_space
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	integer i,j,k,m,n,x
	real*8 rsqr,Mz,rr,c1,c2,hLx,nhLx,hLy,nhLy
	real*8, dimension(3) :: rij,ff

	Mz=0
	c1=lb/Beta
	c2=2*alpha/sqrt(pi)
	hLx=Lx/2
	nhLx=-Lx/2
	hLy=Ly/2
	nhLy=-Ly/2
	if (real_verlet==0) then
		do m=1,Nq-1
			i=charge(m)
			do n=m+1,Nq
				j=charge(n)
				call rij_and_rr(rij,rsqr,i,j)
				if (rsqr<rcc2) then
					rr=sqrt(rsqr)
					ff=pos(i,4)*pos(j,4)/rsqr*( erfc(alpha*rr)/rr + c2*exp(-alpha2*rsqr) )*rij*c1
					acc_c(i,:)=acc_c(i,:)+ff
					acc_c(j,:)=acc_c(j,:)-ff
				end if
			end do
		end do
	else
		do k=1,npair2
			i=real_pair_list(k,1)
			j=real_pair_list(k,2)
			rij=pos(i,1:3)-pos(j,1:3)
			if (rij(1)>hLx) then
				rij(1)=rij(1)-Lx
				elseif(rij(1)<=nhLx) then
					rij(1)=rij(1)+Lx
				end if
				if (rij(2)>hLy) then
					rij(2)=rij(2)-Ly
					elseif(rij(2)<=nhLy) then
						rij(2)=rij(2)+Ly
					end if
					rsqr=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
					if (rsqr<rcc2) then
! 						x=nint(rsqr/real_itv)
! 						acc_c(i,:)=acc_c(i,:)+pos(i,4)*pos(j,4)*real_fun(x)*rij
						rr=sqrt(rsqr)
						ff=pos(i,4)*pos(j,4)/rsqr*( erfc(alpha*rr)/rr + c2*exp(-alpha2*rsqr) )*rij*c1
						acc_c(i,:)=acc_c(i,:)+ff
					end if
				end do
			end if
	!compute the modified force
	do m=1,Nq
		i=charge(m)
		Mz=Mz+pos(i,4)*pos(i,3)*lb/Beta
	end do
	do m=1,Nq
		i=charge(m)
		acc_c(i,3)=acc_c(i,3)-4*pi/Lx/Ly/Lz/Z_empty*pos(i,4)*Mz
	end do
end subroutine real_space
!###########################################!


!###############Standard_Ewald##############!
subroutine Standard_Ewald
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	complex eikx(1:Nq, -Kmax1:Kmax1)
	complex eiky(1:Nq, -Kmax2:Kmax2)
	complex eikz(1:Nq, 0:Kmax3)
	complex eikr(Nq),sum
	integer i,j,l,m,n,p,q,r,tot_k
	real*8 c1,c2,c3,ksqr,kcut
	real*8 zq(Nq),acc_f(1:Nq,3)
	real*8, dimension(3) :: kk
	complex, dimension(3) :: kk1 

	acc_f=0
	do m=1,Nq
		i=charge(m)
		zq(m)=pos(i,4)
	end do
	c1=2*pi/Lx
	c2=2*pi/Ly
	c3=2*pi/Lz/Z_empty
	do m=1,Nq
		i=charge(m)
		eikx(m,0)=(1,0)
		eiky(m,0)=(1,0)
		eikz(m,0)=(1,0)

		eikx(m,1)=cmplx( cos(c1*pos(i,1)), sin(c1*pos(i,1)) )
		eiky(m,1)=cmplx( cos(c2*pos(i,2)), sin(c2*pos(i,2)) )
		eikz(m,1)=cmplx( cos(c3*pos(i,3)), sin(c3*pos(i,3)) )

		eikx(m,-1)=conjg(eikx(m,1))
		eiky(m,-1)=conjg(eiky(m,1))
	end do

	do p=2, Kmax1
		do m=1, Nq
			eikx(m,p)=eikx(m,p-1)*eikx(m,1)
			eikx(m,-p)=conjg(eikx(m,p))
		end do
	end do
	do q=2, Kmax2
		do m=1, Nq
			eiky(m,q)=eiky(m,q-1)*eiky(m,1)
			eiky(m,-q)=conjg(eiky(m,q))
		end do
	end do
	do r=2, Kmax3
		do m=1, Nq
			eikz(m,r)=eikz(m,r-1)*eikz(m,1)
		end do
	end do

	tot_k=0
	do r=0,Kmax3
		do p=-Kmax1,Kmax1
			do q=-Kmax2,Kmax2
				kcut=(1.*p/Kmax1)*(1.*p/Kmax1)+(1.*q/Kmax2)*(1.*q/Kmax2)+(1.*r/Kmax3)*(1.*r/Kmax3)
				if ( kcut>1 .or. kcut==0) cycle
				kk(1)=p*c1
				kk(2)=q*c2
				kk(3)=r*c3
				ksqr=kk(1)*kk(1)+kk(2)*kk(2)+kk(3)*kk(3)
				sum=(0.0,0.0)
				tot_k=tot_k+1
				do m=1,Nq
					eikr(m)=zq(m)*eikx(m,p)*eiky(m,q)*eikz(m,r)
					sum=sum+eikr(m)
				end do
				kk1=conjg(sum)*kk
				do m=1,Nq
					acc_f(m,:)=acc_f(m,:)+aimag( eikr(m)*kk1 )*exp_ksqr(tot_k)
				end do
			end do
		end do
	end do
	do m=1,Nq
		i=charge(m)
		acc_c(i,:)=acc_c(i,:)+acc_f(m,:)
	end do
end subroutine Standard_Ewald
!###########################################!


!#################SPME_Ewald################!
subroutine SPME_Ewald
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, dimension(Nq,ordr(1)) :: SPx
	real*8, dimension(Nq,ordr(1)) :: dSPx
	real*8, dimension(Nq,ordr(2)) :: SPy
	real*8, dimension(Nq,ordr(2)) :: dSPy
	real*8, dimension(Nq,ordr(3)) :: SPz
	real*8, dimension(Nq,ordr(3)) :: dSPz
	integer, dimension(Nq,ordr(1)) :: iqmap
	integer, dimension(Nq,ordr(2)) :: jqmap
	integer, dimension(Nq,ordr(3)) :: kqmap
	real*8 :: st,fn,tm
	real*8, dimension(Nq,3) :: acc_f
	integer :: i,m
	
! 	call cpu_time(st)
call splcof(SPx, dSPx, iqmap, ordr(1), 1)
call splcof(SPy, dSPy, jqmap, ordr(2), 2)
call splcof(SPz, dSPz, kqmap, ordr(3), 3)
! 	call cpu_time(fn)
! 	write(*,*) fn-st

! 	call cpu_time(st)
call MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)
! 	call cpu_time(fn)
! 	write(*,*) fn-st

! 	call cpu_time(st)
call pmeOrthoConvBC
! 	call cpu_time(fn)
! 	write(*,*) fn-st

! 	call cpu_time(st)
call IFrc(acc_f, SPx, SPy, SPz, dSPx, dSPy, dSPz, iqmap, jqmap, kqmap)
! 	call cpu_time(fn)
! 	write(*,*) fn-st

call ZeroForce(acc_f)

! 	write(*,*) charge(1:100)
! 	write(*,*) acc_f(1:10,1)
! 	stop
	acc_c(charge,:)=acc_c(charge,:)+acc_f*lb/Beta   !Nq和NN不一样的时候这是不对的
! 	do m=1,Nq
! 		i=charge(m)
! 		acc_c(i,:)=acc_c(i,:)+acc_f(m,:)*lb/Beta
! 	end do
! 	write(*,*) acc(1:10,1)
! 	stop
end subroutine SPME_Ewald
!###########################################!


!################fene_force#################!
subroutine fene_force
	!-----------------------------------------!
	!Note: rij=ri-rj
	!-----------------------------------------!
	implicit none
	integer :: i,j,k,n
	real*8 :: rsqr
	real*8, dimension(3) :: rjk,ff
	do i=1, N_bond
		j=fene_list(i,1)
		k=fene_list(i,2)
		call rij_and_rr(rjk, rsqr, j, k)
		if(rsqr>R0_2) then
			write(*,*) 'The bond is break off!'
			write(*,*) k*1.,j*1.,sqrt(rsqr)
			call write_pos
			call write_vel
			call write_acc
			stop
		end if
		ff=-kfene*R0_2/(R0_2-rsqr)*rjk
		acc(j,:)=acc(j,:)+ff
		acc(k,:)=acc(k,:)-ff
	end do
end subroutine fene_force
!###########################################!

!###################splcof##################!
subroutine splcof(SP, dSP, qmap, ord, xyz)
	!-----------------------------------------!
	!ord指定数组大小可以吗，回去查查, 可以的
	!不能在定义时候初始化
	!-----------------------------------------!
	implicit none
	integer, intent(in) :: ord
	real*8, dimension(Nq,ord), intent(out) :: SP
	real*8, dimension(Nq,ord), intent(out) :: dSP
	integer, dimension(Nq,ord), intent(out) :: qmap
	integer, intent(in) :: xyz
	real*8, dimension(Nq) :: uu,ww
	integer, dimension(Nq) :: uuc
	integer :: i,j,k
	SP=0
	dSP=0
	qmap=0
	
  uu = ng(xyz)/gdim(xyz)*posq(:,xyz)						!化到标定后的分数坐标
  uuc = ceiling(uu)														!ceiling can also be used for array
  
  if (ord==4) then
  	ww = uu-uuc+1
  	
  	SP(:,3) = 0.5*ww*ww
  	SP(:,1) = 0.5*(1.0-ww)*(1.0-ww)
  	SP(:,2) = 1.0-SP(:,1)-SP(:,3)

  	dSP(:,1) = -SP(:,1)
  	dSP(:,2) = SP(:,1) - SP(:,2)
  	dSP(:,3) = SP(:,2) - SP(:,3)
  	dSP(:,4) = SP(:,3)
  	
  	SP(:,4) = ww*SP(:,3)/3
  	SP(:,3) = ((ww+1.0)*SP(:,2) + (3.0-ww)*SP(:,3))/3
  	SP(:,1) = (1.0-ww)*SP(:,1)/3
  	SP(:,2) = 1.0-SP(:,1)-SP(:,3)-SP(:,4)
  else
  	ww = uu-uuc+1
  	do k=1,ord
  		SP(:,ord-k+1)=bspln_cof(k,1,xyz)*ww
  		dSP(:,ord-k+1)=dspln_cof(k,1,xyz)*ww	
  		do j=2,ord-2
  			SP(:,ord-k+1)=ww*(SP(:,ord-k+1)+bspln_cof(k,j,xyz))
  			dSP(:,ord-k+1)=ww*(dSP(:,ord-k+1)+dspln_cof(k,j,xyz))
  		end do
  		SP(:,ord-k+1)=ww*(SP(:,ord-k+1)+bspln_cof(k,ord-1,xyz))
  		SP(:,ord-k+1)=SP(:,ord-k+1)+bspln_cof(k,ord,xyz)
  		dSP(:,ord-k+1)=dSP(:,ord-k+1)+dspln_cof(k,ord-1,xyz)
  		ww=ww+1		
  	end do
  end if

  do i = 0,ord-1
  	j = i + 1
  	qmap(:,j) = uuc(:) - ord + i
  	qmap(:,j) = qmap(:,j) - ng(xyz)*floor(1.*qmap(:,j)/ng(xyz))
  end do
  
  qmap = qmap + 1
  
  if (xyz == 1) then
  	do i = 1,ord
  		do j = 1,Nq
  			SP(j,i) = SP(j,i)*posq(j,4)
  			dSP(j,i) = dSP(j,i)*posq(j,4)
  		end do
  	end do
  end if
  
end subroutine splcof
!###########################################!


!#################MapCharges################!
subroutine MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, dimension(Nq,ordr(1)), intent(in) :: SPx
	real*8, dimension(Nq,ordr(2)), intent(in) :: SPy
	real*8, dimension(Nq,ordr(3)), intent(in) :: SPz
	integer, dimension(Nq,ordr(1)), intent(in) :: iqmap
	integer, dimension(Nq,ordr(2)), intent(in) :: jqmap
	integer, dimension(Nq,ordr(3)), intent(in) :: kqmap
	integer :: h,i,j,k,ii,jj,kk,t,ng1o,ng2o,ng3o,ih,jh,kh
	real*8 :: dqx,dqxy,dd,qmij(ordr(1),ordr(2)),SPxh(ordr(1),1),SPyh(1,ordr(2)),SPzh(ordr(3))
	real (kind=8), dimension(ordr(1),ordr(2),ordr(3)) :: qmijk
	integer :: imp(ordr(1)),jmp(ordr(2)),kmp(ordr(3))
	Q_PME=0
	
	ng1o=ng(1)-ordr(1)+2
	ng2o=ng(2)-ordr(2)+2
	ng3o=ng(3)-ordr(3)+2
	qmijk=0
	do h=1,Nq
		imp=iqmap(h,:)
		jmp=jqmap(h,:)
		kmp=kqmap(h,:)
		ih=imp(1)
		jh=jmp(1)
		kh=kmp(1)
		
! 		if (ih<ng1o .and. jh<ng2o .and. kh<ng3o) then
SPxh(:,1)=SPx(h,:)
SPyh(1,:)=SPy(h,:)
SPzh=SPz(h,:)
qmij=matmul(SPxh,SPyh)
do i=1,ordr(3)
	qmijk(:,:,i)=SPzh(i)*qmij
end do
Q_PME(imp,jmp,kmp)=Q_PME(imp,jmp,kmp)+qmijk
! 		else
! 			do i=1,ordr(1)
! 				ii=iqmap(h,i)
! 				dqx=SPx(h,i)
! 				do j=1,ordr(2)
! 					jj=jqmap(h,j)
! 					dqxy=dqx*SPy(h,j)
! 					do k=1,ordr(3)
! 						kk=kqmap(h,k)
! 						Q_PME(ii,jj,kk)=Q_PME(ii,jj,kk)+dqxy*SPz(h,k)
! 					end do
! 				end do
! 			end do
! 		end if
end do

end subroutine Mapcharges
!###########################################!

!###############pmeOrthoConvBC##############!
subroutine pmeOrthoConvBC
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	include "fftw3.f90"
	complex ( kind = 8 ), dimension(ng(1),ng(2),ng(3)) :: fQ_PME
	complex ( kind = 8 ), dimension(ng(1),ng(2),ng(3)) :: fQBC
	integer i,j,k
	integer ( kind = 8 ) plan_backward
	integer ( kind = 8 ) plan_forward
	integer ( kind = 8 ) plan
	
  call dfftw_plan_dft_3d_ ( plan_forward, ng(1), ng(2), ng(3), Q_PME, fQ_PME, FFTW_FORWARD, FFTW_Estimate ) !MEASURE !Estimate

  call dfftw_execute_ ( plan_forward )
  
  call dfftw_destroy_plan_ ( plan_forward )
  
	fQ_PME(1,1,1)=0				  !fft能变换出无穷的数吗
	
! 	do i=1, ng(1)
! 	  do j=1,ng(2)
! 	    do k = 1,ng(3)
! 				Enrg=Enrg+0.5*abs(fQ_PME(i,j,k)*fQ_PME(i,j,k))*B_PME(i,j,k)*C_PME(i,j,k)
! 	    end do
! 	  end do
! 	end do
! 	Enrg=Enrg/(ng(1)*ng(2)*ng(3))

fQBC=fQ_PME*BC_PME

call dfftw_plan_dft_3d_ ( plan_backward, ng(1), ng(2), ng(3), fQBC, U_PME, FFTW_BACKWARD, FFTW_Estimate )

call dfftw_execute_ ( plan_backward )	

call dfftw_destroy_plan_ ( plan_backward )

! 	U_PME=U_PME/(ng(1)*ng(2)*ng(3))

end subroutine pmeOrthoConvBC
!###########################################!


!####################IFrc###################!
subroutine IFrc(acc_f, Bx, By, Bz, dBx, dBy, dBz, iqm, jqm, kqm)
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, dimension(Nq,3), intent(out) :: acc_f
	real*8, dimension(Nq,ordr(1)), intent(in) :: Bx
	real*8, dimension(Nq,ordr(2)), intent(in) :: By
	real*8, dimension(Nq,ordr(3)), intent(in) :: Bz
	real*8, dimension(Nq,ordr(1)), intent(in) :: dBx
	real*8, dimension(Nq,ordr(2)), intent(in) :: dBy
	real*8, dimension(Nq,ordr(3)), intent(in) :: dBz
	integer, dimension(Nq,ordr(1)), intent(in) :: iqm
	integer, dimension(Nq,ordr(2)), intent(in) :: jqm
	integer, dimension(Nq,ordr(3)), intent(in) :: kqm
	real*8 :: dqx,tqx,tdqxy,dtqxy,ddqxy
	integer :: h,i,j,k,ii,jj,kk,ng1o,ng2o,ng3o,ih,jh,kh
	complex (kind=8) :: f1,f2,f3,ulm
	real*8 :: Bxh(ordr(1),1),dBxh(ordr(1),1),Byh(1,ordr(2)),dByh(1,ordr(2)),Bzh(ordr(3)),dBzh(ordr(3))
	real (kind=8), dimension(ordr(1),ordr(2),ordr(3)) :: dqmijkx,dqmijky,dqmijkz
	complex (kind=8), dimension(ordr(1),ordr(2),ordr(3)) :: Ureg
	real (kind=8), dimension(ordr(1),ordr(2)) :: dqmijx,dqmijy,dqmijz
	acc_f=0
	
	ng1o=ng(1)-ordr(1)+2;
	ng2o=ng(2)-ordr(2)+2;
	ng3o=ng(3)-ordr(3)+2;
	dqmijkx=0
	dqmijky=0
	dqmijkz=0
	do h=1,Nq
		f1=0
		f2=0
		f3=0
		ih=iqm(h,1)
		jh=jqm(h,2)
		kh=kqm(h,3)
! 		if (ih<ng1o .and. jh<ng2o .and. kh<ng3o) then
Bxh(:,1)=Bx(h,:)
Byh(1,:)=By(h,:)
Bzh=Bz(h,:)
dBxh(:,1)=dBx(h,:)
dByh(1,:)=dBy(h,:)
dBzh=dBz(h,:)
dqmijx=matmul(dBxh,Byh)
dqmijy=matmul(Bxh,dByh)
dqmijz=matmul(Bxh,Byh)
do i=1,ordr(3)
	dqmijkx(:,:,i) = Bzh(i)*dqmijx
	dqmijky(:,:,i) = Bzh(i)*dqmijy
	dqmijkz(:,:,i) = dBzh(i)*dqmijz
end do
Ureg = U_PME(iqm(h,:),jqm(h,:),kqm(h,:))
      f1 = -sum(dqmijkx*Ureg)									!the coefficient of the inverse fft in fftw
      f2 = -sum(dqmijky*Ureg)
      f3 = -sum(dqmijkz*Ureg)
! 		else
! 			ii=1
! 			do i=1,ordr(1)
! 				jj=1
! 				dqx=Bx(h,ii)
! 				tqx=dBx(h,ii)
! 				do j=1,ordr(2)
! 					kk=1
! 					tdqxy=tqx*By(h,jj)
! 					dtqxy=dqx*dBy(h,jj)
! 					ddqxy=dqx*By(h,jj)
! 					do k=1,ordr(3)
! 						ulm=U_PME(iqm(h,i),jqm(h,j),kqm(h,k))
! 						f1=f1-tdqxy*Bz(h,kk)*ulm            !the coefficient of the inverse fft in fftw
! 						f2=f2-dtqxy*Bz(h,kk)*ulm
! 						f3=f3-ddqxy*dBz(h,kk)*ulm
! 						kk=kk+1
! 					end do
! 					jj=jj+1
! 				end do
! 				ii=ii+1
! 			end do
! 		end if
acc_f(h,:)=(/real(f1)*ng(1)/gdim(1),real(f2)*ng(2)/gdim(2),real(f3)*ng(3)/gdim(3)/)/(ng(1)*ng(2)*ng(3))
end do

end subroutine IFrc
!###########################################!

!##################ZeroForce################!
subroutine ZeroForce(acc_f)
	!-----------------------------------------!
	!
	!-----------------------------------------!
	implicit none
	real*8, dimension(Nq,3),intent(inout) :: acc_f
	real*8, dimension(3) :: mean_acc
	integer :: i
	mean_acc=0
	
	do i=1,Nq
		mean_acc(1)=mean_acc(1)+acc_f(i,1)
		mean_acc(2)=mean_acc(2)+acc_f(i,2)
		mean_acc(3)=mean_acc(3)+acc_f(i,3)
	end do
	mean_acc=mean_acc/Nq
	
	acc_f(:,1)=acc_f(:,1)-mean_acc(1)
	acc_f(:,2)=acc_f(:,2)-mean_acc(2)
	acc_f(:,3)=acc_f(:,3)-mean_acc(3)
	
end subroutine ZeroForce
!###########################################!


end module compute_acc



