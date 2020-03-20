module compute_acceleration
  implicit none

  save
  !
  !lj potential
  real*8,  private :: epsilon     !Energy unit epsilon in lj potential
  real*8,  private :: sigma       !Distance sigma in lj potential
  real*8           :: sigma_s     !diameter of salt
  real*8,  private :: sigma_sp    !Distance sigma in lj potential
  real*8           :: rcl         !Cut off radius of LJ potential
  real*8           :: rcl_s       !Cut off radius of LJ potential
  real*8           :: rcl_sp      !Cut off radius of LJ potential
  real*8           :: rvl         !Verlet list radius of LJ potential
  integer, private :: npair1      !number of pairs in the lj verlet sphere
  real*8           :: dr_max1     !max displacement of lj verlet list
  !
  !fene potential
  real*8           :: R0_2        !max length of FENE potential
  real*8           :: kFENE       !coefficient of FENE potential
  integer          :: N_bond      !Number of Chemical bond of polymers
  !
  !coulomb potential
  real*8           :: lb          !Bjerrum length
  real*8           :: xi          !friction coefficient
  real*8,  private :: EF          !electric field
  real*8,  private :: tol         !tolerance
  real*8,  private :: tau_rf      !ratio of time between real and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  !
  !real space
  real*8,  private :: rcc         !Cut off radius of real space
  real*8,  private :: rcc2        !rcc2=rcc*rcc
  real*8,  private :: rvc         !rvc=rcc+rsk
  real*8           :: rsk         !Skin between cut off sphere and verlet list 
                                  !sphere in real space
  integer          :: real_verlet !call real verlet list or not
  real*8,  private :: real_itv    !Inteval of interpolation of Error Function
  integer, private :: npair2      !number of pairs in the real verlet sphere
  real*8           :: dr_max2     !max displacement of lj verlet list
  !
  !reciprocal space
  integer, private :: Kmax1       !max wave number of x direction
  integer, private :: Kmax2       !max wave number of y direction 
  integer, private :: Kmax3       !max wave number of z direction
  integer, private :: K_total     !Total wave number in reciprocal space
  integer, dimension(3) :: ordr   !Order of spline function
  integer, dimension(3) :: ng     !wave number
  real*8,  dimension(3) :: gdim   !length, width and height of the box
  !
  !arrays
  real*8,  allocatable, dimension(:,:), private :: posq        
                        !array of position of charged particle
  real*8,  allocatable, dimension(:,:), private :: acc_c
                        !array of accelaration
  integer, allocatable, dimension(:,:), private :: lj_pair_list
                        !LJ potential verlet list
  integer, allocatable, dimension(:,:), private :: real_pair_list
                        !real space verlet list
  integer, allocatable, dimension( : )          :: charge
                        !charge number to monomer number
  integer, allocatable, dimension( : )          :: charge1
                        !charge number to monomer number
  integer, allocatable, dimension(:,:), private :: fene_list
                        !pairs of two adjacent monomers
  integer, allocatable, dimension( : )          :: anchor_list
                        !number of the anchored monomer
  real*8,  allocatable, dimension( : ), private :: exp_ksqr 
                        !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: real_fun 
                        !Function list of erfc and coefficients in real space
  real*8,  allocatable, dimension(:,:,:), private :: bspln_cof
                        !Coefficients of b-spline basis
  real*8,  allocatable, dimension(:,:,:), private :: dspln_cof
                        !Coefficients of derivatives of b-spline basis
  complex (kind=8), allocatable, dimension(:,:,:), private :: Q_PME
                        !Q in SPME
  complex (kind=8), allocatable, dimension(:,:,:), private :: U_PME
                        !U in SPME
  complex (kind=8), allocatable, dimension(:,:,:), private :: BC_PME
                        !B*C in SPME

  contains


subroutine initialize_force_parameters
  !--------------------------------------!
  !
  !
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read force parameters from file
  call read_force_parameters
  !
  !Initialize lj parameters and array allocate.
  call initialize_lj_parameters
  !
  !Initialize fene parameters and array allocate.
  call build_fene_list
  !
  !
  if ( qq /= 0 ) then
    !
    !Initialize ewald parameters and array allocate.
    call Initialize_ewald_parameters
    !
    !Construct the coefficients vector in Fourier space
    call fourier_function
    !
    !build charge array
    call build_charge
    !
    !Construct the coefficients vector in real space
    call real_function
    !
    !build coefficient array of b-spline
    call bspln_coeffs
    !
    !build B*C in SPME
    call PME_BC
  end if

  allocate( acc_c(NN,3) )
  allocate( posq(Nq,4)  )
end subroutine initialize_force_parameters


subroutine read_force_parameters

  open(10,file='./force_data.txt')
    read(10,*) epsilon
    read(10,*) sigma  
    read(10,*) sigma_s
    read(10,*) rcl    
    read(10,*) rvl        
    read(10,*) rsk
    read(10,*) R0_2
    read(10,*) kFENE
    read(10,*) ordr(1)
    read(10,*) ordr(2)
    read(10,*) ordr(3)
    read(10,*) lb
    read(10,*) xi
    read(10,*) EF           
    read(10,*) tol
    read(10,*) tau_rf
  close(10)

  sigma_sp = (sigma+sigma_s)/2
  rcl_s = 2**(1./6)*sigma_s
  rcl_sp = 2**(1./6)*sigma_sp
  rvl = rcl_s + 1

end subroutine read_force_parameters


subroutine initialize_lj_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet
  !
  !allocate verlet list of LJ potential
  rho = NN / (Lx * Ly * Lz)
  v_verlet = 8.D0/3 * pi * rvl**3
!   allocate(  lj_pair_list(25*NN*ceiling(rho*v_verlet),2)  )
!   lj_pair_list = 0
  dr_max1 = 0
  npair1  = 0
end subroutine initialize_lj_parameters


subroutine lj_verlet_list
  !--------------------------------------!
  !Construct lj_pair_list and lj_point by the link list
  !method with the complexity of O(N)
  !   
  !Input
  !   pos
  !Output
  !   lj_pair_list
  !External Variables
  !   NN, Lx, Ly, Lz, rv_lj, lj_pair_list, lj_point, pos
  !Routine Referenced:
  !1. rij_and_rr(rij, rr, i, j)
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.550-552. 
  !--------------------------------------!
  use global_variables
  implicit none
  integer i,j,k,l,m,n,p,q,r
  integer icel,jcel,kcel,ncel1,ncel2,ncel3
  real*8, dimension(3) :: rij
  real*8 :: rsqr,rcel1,rcel2,rcel3
  integer, dimension(NN) :: cell_list
  integer, allocatable, dimension(:,:,:) :: hoc
  ncel1=int(Lx/rvl)
  ncel2=int(Ly/rvl)
  ncel3=int(Lz/rvl)
  allocate(hoc(0:ncel1-1,0:ncel2-1,0:ncel3-1))

  hoc=0
  rcel1=Lx/ncel1
  rcel2=Ly/ncel2
  rcel3=Lz/ncel3
  do i=1,NN
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    cell_list(i)=hoc(icel,jcel,kcel)
    hoc(icel,jcel,kcel)=i
  end do

  k=0
  do i=1,NN
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do l=-1,1
      if (icel+l .ge. ncel1) then
        p=icel+l-ncel1
      elseif(icel+l<0) then
        p=icel+l+ncel1
      else
        p=icel+l
      end if
      do m=-1,1
        if (jcel+m .ge. ncel2) then
          q=jcel+m-ncel2
        elseif(jcel+m<0) then
          q=jcel+m+ncel2
        else
          q=jcel+m
        end if
        do n=-1,1
          if (kcel+n .ge. ncel3) then
            cycle
          elseif(kcel+n<0) then
            cycle
          else
            r=kcel+n
          end if
          j=hoc(p,q,r)
          do while (j /= 0)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rvl*rvl) ) then
              k=k+1
            end if
            j=cell_list(j)
          end do
        end do
      end do
    end do
  end do
  npair1=k
  if ( allocated(lj_pair_list) ) deallocate(lj_pair_list)
  allocate(lj_pair_list(npair1,2))
  k=0
  do i=1,NN
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do l=-1,1
      if (icel+l .ge. ncel1) then
        p=icel+l-ncel1
      elseif(icel+l<0) then
        p=icel+l+ncel1
      else
        p=icel+l
      end if
      do m=-1,1
        if (jcel+m .ge. ncel2) then
          q=jcel+m-ncel2
        elseif(jcel+m<0) then
          q=jcel+m+ncel2
        else
          q=jcel+m
        end if
        do n=-1,1
          if (kcel+n .ge. ncel3) then
            cycle
          elseif(kcel+n<0) then
            cycle
          else
            r=kcel+n
          end if
          j=hoc(p,q,r)
          do while (j /= 0)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rvl*rvl) ) then
              k=k+1
              lj_pair_list(k,1)=i
              lj_pair_list(k,2)=j
            end if
            j=cell_list(j)
          end do
        end do
      end do
    end do
  end do
end subroutine lj_verlet_list


subroutine build_fene_list
  !--------------------------------------!
  !Construct the fene_list array.
  !   
  !Input
  !   
  !Output
  !   fene_list, anchor_list
  !External Variables
  !   N_bond, Npe, Nml, Ngl
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l

  N_bond=arm*Nma*Nga+(Nml-1)*Ngl
  allocate( fene_list(N_bond,2) )
  allocate( anchor_list(N_anchor) )
  !
  !anchors of the polymers 
  do i = 1, Nga
    anchor_list(i) = (i-1)*(arm*Nma+1) + 1
  end do
  do i = 1, Ngl
    anchor_list(i+Nga) = Nta + (i-1)*Nml + 1
  end do
  !
  !star brushes fene list
  l=0
  do i = 1, Nga
    do k = 1, Nma     !the anchor to the branching point
      l = l + 1
      fene_list(l,1) = (i-1)*(arm*Nma+1) + k
      fene_list(l,2) = (i-1)*(arm*Nma+1) + k+1
    end do
    do j = 2, arm     !the arms
      l = l + 1       !the branching point to the first point of the arms
      fene_list(l,1) = (i-1)*(arm*Nma+1) + Nma+1
      fene_list(l,2) = (i-1)*(arm*Nma+1) + (Nma*(j-1)+1) + 1
      do k=1,Nma-1  !the first point of the arms to the last point of the arms
        l = l + 1
        fene_list(l,1) = (i-1)*(arm*Nma+1) + (Nma*(j-1)+1) + k
        fene_list(l,2) = (i-1)*(arm*Nma+1) + (Nma*(j-1)+1) + k+1
      end do
    end do
  end do
  !
  !star brushes fene list   
  do i = 1, Ngl
    do j = 1, Nml-1
      l = l + 1
      fene_list(l,1) = Nta + (i-1)*Nml + j
      fene_list(l,2) = Nta + (i-1)*Nml + j+1
    end do
  end do
end subroutine build_fene_list


subroutine Initialize_ewald_parameters
  !--------------------------------------!
  !Calculate alpha, Kmax1, Kmax2, Kmax3, rcc
  !
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306. 
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rho, v_verlet

  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly*Lz)**2 ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rcc  = tol / alpha
  rcc2 = rcc * rcc
  rvc  = rcc + rsk
  !
  !use verlet list in real space
  if ( ( int(Lx/rvc) * int(Ly/rvc) * int(Lz/rvc) ) > 27 ) then 
    Kmax1 = ceiling(tol*Lx*alpha/pi)
    Kmax2 = ceiling(tol*Ly*alpha/pi)
    Kmax3 = ceiling(tol*Lz*Z_empty*alpha/pi)
    real_verlet = 1
  !
  !don't use verlet list in real space
  else
    if ( Lx > Ly ) then
      rcc = Ly/2
    else
      rcc = Lx/2
    end if
    rcc2 = rcc * rcc
    Kmax1    = ceiling(tol*tol/pi*Lx/rcc)
    Kmax2    = ceiling(tol*tol/pi*Ly/rcc)
    Kmax3    = ceiling(tol*tol/pi*Lz*Z_empty/rcc)
    alpha    = tol / rcc
    alpha2   = alpha * alpha
    rvc  = rcc + rsk
    real_verlet = 0
  end if
  Kmax1 = ceiling(Kmax1*3.0)
  Kmax2 = ceiling(Kmax2*3.0)
  Kmax3 = ceiling(Kmax3*3.0)
  ng    = (/Kmax1,Kmax2,Kmax3/)
  gdim  = (/Lx,Ly,Lz*Z_empty/)
  !
  !allocate verlet list of real space
  rho = Nq / (Lx * Ly * Lz)
  v_verlet = 8.D0/3 * pi * rvc**3
!   if ( allocated(real_pair_list) ) deallocate(real_pair_list)
!   allocate( real_pair_list(25*Nq*ceiling(rho*v_verlet),2) )
!   write(*,*) v_verlet,25*Nq*ceiling(rho*v_verlet)
!   real_pair_list = 0
  dr_max2 = 0
end subroutine Initialize_ewald_parameters


subroutine build_charge
  !--------------------------------------!
  !choose the charged particle from pos array to form
  !a new array.
  !   
  !Input
  !   pos
  !Output
  !   charge
  !External Variables
  !   pos, charge, NN, Nq
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l

  if ( allocated(charge) ) deallocate(charge)
  allocate( charge(Nq) )
  if ( allocated(charge1) ) deallocate(charge1)
  allocate( charge1(Nq) )
  
  l=0
  !star brushes
  do i=1,Nga
    if ( charged_arm(1,1) /= 0 ) then
      do k = 2, Nma+1     !the chain anchored to the plate
        if ( mod( k-1, man_s ) == 0 ) then
          l = l + 1
          charge(l) = (i-1) * (arm*Nma+1) + k
          charge1(l) = charged_arm(1,2)
        end if
      end do
    end if
    do j = 2, arm  
      if (charged_arm(j,1) /=0 ) then 
        do k = 1, Nma
          if ( mod( k, man_s ) == 0) then
            l = l + 1
            charge(l) = (i-1)*(arm*Nma+1) + ((j-1)*Nma+1) + k
            charge1(l) = charged_arm(j,2)
          end if
        end do
      end if
    end do
  end do
  !linear brushes
  do i=1,Ngl
    do j=2,Nml
      if (mod( j-1, man_l )==0) then
        l=l+1
        charge(l)=Nta+(i-1)*Nml+j
      end if
    end do
  end do
  !ions
  do i=Npe+1,NN
    l=l+1
    charge(l)=i
  end do

end subroutine build_charge


subroutine real_verlet_list
  !--------------------------------------!
  !Construct real_pair_list and real_point by the link list
  !method with the complexity of O(N)
  !   
  !Input
  !   pos
  !Output
  !   real_pair_list
  !External Variables
  !   Nq, Lx, Ly, Lz, rv_real,
  !   real_pair_list, real_point, pos, charge
  !Routine Referenced:
  !1. rij_and_rr(rij, rr, i, j)
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.550-552. 
  !--------------------------------------!
  use global_variables
  implicit none
  integer i,j,k,m,n,p,q,r,u,v,w
  integer icel,jcel,kcel,ncel1,ncel2,ncel3
  real*8, dimension(3) :: rij
  real*8 :: rsqr,rcel1,rcel2,rcel3
  integer, dimension(Nq) :: cell_list
  integer,allocatable,dimension(:,:,:)::hoc
  ncel1=int(Lx/rvc)
  ncel2=int(Ly/rvc)
  ncel3=int(Lz/rvc)
  allocate(hoc(0:ncel1-1,0:ncel2-1,0:ncel3-1))

  hoc=0
  rcel1=Lx/ncel1
  rcel2=Ly/ncel2
  rcel3=Lz/ncel3
  do m=1,Nq
    i=charge(m)
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    cell_list(m)=hoc(icel,jcel,kcel)
    hoc(icel,jcel,kcel)=m
  end do
  k=0
  do m=1,Nq
    i=charge(m)
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do u=-1,1
      if (icel+u>=ncel1) then
        p=icel+u-ncel1
      elseif(icel+u<0) then
        p=icel+u+ncel1
      else
        p=icel+u
      end if
      do v=-1,1
        if (jcel+v>=ncel2) then
          q=jcel+v-ncel2
        elseif(jcel+v<0) then
          q=jcel+v+ncel2
        else
          q=jcel+v
        end if
        do w=-1,1
          if (kcel+w>=ncel3) then
            cycle
          elseif(kcel+w<0) then
            cycle
          else
            r=kcel+w
          end if
          n=hoc(p,q,r)
          do while (n /= 0)
            j=charge(n)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rvc*rvc) ) then
              k=k+1
            end if
            n=cell_list(n)
          end do
        end do
      end do
    end do
  end do
  npair2=k
  if ( allocated(real_pair_list) ) deallocate(real_pair_list)
  allocate(real_pair_list(npair2,2))
  k=0
  do m=1,Nq
    i=charge(m)
    icel=int((pos(i,1)+Lx/2)/rcel1)
    jcel=int((pos(i,2)+Ly/2)/rcel2)
    kcel=int(pos(i,3)/rcel3)
    do u=-1,1
      if (icel+u>=ncel1) then
        p=icel+u-ncel1
      elseif(icel+u<0) then
        p=icel+u+ncel1
      else
        p=icel+u
      end if
      do v=-1,1
        if (jcel+v>=ncel2) then
          q=jcel+v-ncel2
        elseif(jcel+v<0) then
          q=jcel+v+ncel2
        else
          q=jcel+v
        end if
        do w=-1,1
          if (kcel+w>=ncel3) then
            cycle
          elseif(kcel+w<0) then
            cycle
          else
            r=kcel+w
          end if
          n=hoc(p,q,r)
          do while (n /= 0)
            j=charge(n)
            call rij_and_rr(rij,rsqr,i,j)
            if ( i/=j .and. rsqr<(rvc*rvc) ) then
              k=k+1
              real_pair_list(k,1)=i
              real_pair_list(k,2)=j
            end if
            n=cell_list(n)
          end do
        end do
      end do
    end do
  end do
  npair2=k
end subroutine real_verlet_list


subroutine fourier_function
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
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

  if ( allocated(exp_ksqr) ) deallocate( exp_ksqr )
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


subroutine real_function
  !----------------------------------------!
  !make a erfc funciton table to fast calculate the function.
  !however, the accuracy is very low, so we havn't use it in this
  !program.
  !----------------------------------------!
  use global_variables
  implicit none
  real*8 :: rr,rsqr,c1,c2,rf
  integer :: N=100000,i,j
  if ( allocated(real_fun) ) deallocate( real_fun )
  allocate( real_fun(N) )
  real_fun=0.
  c1=lb/Beta
  c2=2*alpha/sqrt(pi)
  real_itv=(tol+0.1)**2/N/alpha2
  do i=1, N
    rsqr=i*real_itv
    rr=sqrt(rsqr)
    rsqr=rr*rr
    real_fun(i)=c1/rsqr*( erfc(alpha*rr)/rr + c2*exp(-alpha2*rsqr) )
  end do
end subroutine real_function


subroutine update_verlet_list
  !--------------------------------------!
  !When the max displacement of movement is greater than
  !skin between cut off sphere and verlet list sphere,
  !we need to update verlet list.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !  lj_verlet_list
  !  real_verlet_list
  !--------------------------------------!
  implicit none

  if ( dr_max1>(rvl-rcl_s)/2 ) then
    call lj_verlet_list
    dr_max1=0
  end if
  if ( dr_max2>rsk/2 .and. real_verlet /=0 ) then
    call real_verlet_list
    dr_max2=0
  end if

end subroutine update_verlet_list


subroutine error_analysis
  !--------------------------------------!
  !
  !Reference:
  !1.The error is estimated by Kolafa1992,
  !  JIRI KOLAFA, JOHN W. PERRAM, 'CUTOFF ERRORS IN THE 
  !  EWALD SUMMATION FORMULAE FOR POINT CHARGE SYSTEMS',
  !  Molecular Simulation, Vol. 9(5), pp.351-368, (1992).
  !2.The true error e.g. RMS force is defined at:
  !  Ulrich Essmann, Lalith Perera et al, 'A smooth particle
  !  mesh Ewald method', Journal of Chemical Physics, Vol. 103(9),
  !  (1995).
  !--------------------------------------!
  use global_variables
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
  write(*,*) 'Standard Ewald error in real space:', f_r
  write(*,*) 'Standard Ewald error in fourier spqce', f_k
  tol1=tol
  tau_rf1=tau_rf
  tol=5                   !the error is about 1e-6
  tau_rf=10
  call Initialize_ewald_parameters
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
  call Initialize_ewald_parameters
!   call bspln_coeffs       !coefficients of bspline
!   call PME_BC             !compute B*C
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
  
end subroutine error_analysis


subroutine Compute_Force
  !--------------------------------------!
  !Each step update lj force and fene force
  !however, the coulomb is updated at each multisteps.
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !
  !Reference:
  !1.The calculation of random force is referred to:
  !  Tamar Schlick, 'Molecular Modeling and Simulation:
  !  An Interdisciplinary Guide', 2nd edition, pp.482,
  !  Springer, 2010.
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, dimension(3) :: rij
  real*8 :: ff,r2,eta1,eta2,eta3,st,fn
  integer :: i,j,k
  
  acc=0
  if ( qq /= 0 .and. mod(step,multistep) == 0 ) then    
    acc_c = 0
    posq(:,:) = pos(charge,:)
    call SPME_Ewald
    call real_space
  end if
  acc = acc_c
  call lj_force
  call fene_force
  do i=1, NN
    call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta1)
    call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta2)
    call gauss_dist(0.D0,sqrt(2./Beta/dt*xi),eta3)
    acc(i,1) = acc(i,1) + eta1
    acc(i,2) = acc(i,2) + eta2
    acc(i,3) = acc(i,3) + eta3 + pos(i,4)*EF          
  end do
  do i = 1, N_anchor
    j = anchor_list(i)
    acc(j,1:3) = 0
  end do
end subroutine Compute_Force


subroutine fene_force
  !-----------------------------------------!
  !Note: rij=ri-rj
  !-----------------------------------------!
 !   use input_output, only : write_pos, write_vel, write_acc
  use global_variables
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
      write(*,*) 'step:', step
!       call write_pos
!       call write_vel
!       call write_acc
      stop
    end if
    ff=-kfene*R0_2/(R0_2-rsqr)*rjk
    acc(j,:)=acc(j,:)+ff
    acc(k,:)=acc(k,:)-ff
  end do
end subroutine fene_force


subroutine lj_force
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  integer i,j,k,n,m
  real*8, dimension(3) :: rij
  real*8 :: rsqr,inv_r2,inv_r6,dl,cc,hLx,nhLx,hLy,nhLy
  real*8 :: rcl2,rcl3,rcl4,sigma2,sigma3
  real*8, dimension(NN,3) :: acc_lj

  hLx=Lx/2
  nhLx=-Lx/2
  hLy=Ly/2
  nhLy=-Ly/2
  rcl2=rcl*rcl
  rcl3=rcl_s*rcl_s
  rcl4=rcl_sp*rcl_sp
  sigma2=sigma_s*sigma_s
  sigma3=sigma_sp*sigma_sp
  acc_lj=0

  do k=1,npair1
    i=lj_pair_list(k,1)
    j=lj_pair_list(k,2)
    rij=pos(i,1:3)-pos(j,1:3)
    m = abs(pos(i,4))+abs(pos(j,4))
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
    if (m==6) then
      if (rsqr<rcl3) then
        inv_r2=sigma2/rsqr
        inv_r6=inv_r2*inv_r2*inv_r2
        acc_lj(i,:)=acc_lj(i,:)+inv_r6*(inv_r6-0.5)*rij/rsqr
      end if
    elseif (m<3) then
      if (rsqr<rcl2) then
        inv_r2=1/rsqr
        inv_r6=inv_r2*inv_r2*inv_r2
        acc_lj(i,:)=acc_lj(i,:)+inv_r6*(inv_r6-0.5)*rij/rsqr
      end if
    else
      if (rsqr<rcl4) then
        inv_r2=sigma3/rsqr
        inv_r6=inv_r2*inv_r2*inv_r2
        acc_lj(i,:)=acc_lj(i,:)+inv_r6*(inv_r6-0.5)*rij/rsqr
      end if          
    end if
  end do
  acc_lj=acc_lj*48
  !wall force
  do i=1,NN
    if (abs(pos(i,3))<0.01 .or. abs(pos(i,3)-Lz)<0.01) cycle
    if (pos(i,4)/=3) then
      if ( pos(i,3)<(rcl/(5**(1./6))) ) then
        dl=pos(i,3)
        acc_lj(i,3)=acc_lj(i,3)+0.4*(3/(dl**11)-1/(dl**5))*dl
      elseif ( (Lz-pos(i,3))<(rcl/(5**(1./6))) ) then
        dl=Lz-pos(i,3)
        acc_lj(i,3)=acc_lj(i,3)-0.4*(3/(dl**11)-1/(dl**5))*dl
      end if
    else
      if ( pos(i,3)<(rcl_s/(5**(1./6))) ) then
        dl=pos(i,3)
        acc_lj(i,3)=acc_lj(i,3)+0.4*(3*((sigma_s/dl)**11)-1*((sigma_s/dl)**5))*dl/sigma_s/sigma_s
      elseif ( (Lz-pos(i,3))<(rcl_s/(5**(1./6))) ) then
        dl=Lz-pos(i,3)
        acc_lj(i,3)=acc_lj(i,3)-0.4*(3*((sigma_s/dl)**11)-1*((sigma_s/dl)**5))*dl/sigma_s/sigma_s
      end if
    end if
  end do
  acc=acc+acc_lj
end subroutine lj_force


subroutine real_space
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
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
!         !
!         !Calculate acc by look up table.
!         x=nint(rsqr/real_itv)
!         acc_c(i,:)=acc_c(i,:)+pos(i,4)*pos(j,4)*real_fun(x)*rij
        !
        !Calculate acc by fortran function.
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


subroutine Standard_Ewald
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  complex(kind=8) :: eikx(1:Nq, -Kmax1:Kmax1)
  complex(kind=8) :: eiky(1:Nq, -Kmax2:Kmax2)
  complex(kind=8) :: eikz(1:Nq, 0:Kmax3)
  complex(kind=8) :: eikr(Nq), sum
  integer i,j,l,m,n,p,q,r,tot_k
  real*8 c1,c2,c3,ksqr,kcut
  real*8 zq(Nq),acc_f(1:Nq,3)
  real*8, dimension(3) :: kk
  complex(kind=8), dimension(3) :: kk1 

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

    eikx(m,1)=cmplx( cos(c1*pos(i,1)), sin(c1*pos(i,1)), 8 )
    eiky(m,1)=cmplx( cos(c2*pos(i,2)), sin(c2*pos(i,2)), 8 )
    eikz(m,1)=cmplx( cos(c3*pos(i,3)), sin(c3*pos(i,3)), 8 )

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

subroutine force_distribution(fene_f,lj_force_PE,lj_force_ions, &
                              coulomb_f,Bond_dist)
  use global_variables
  implicit none
  integer, dimension(500,500), intent(inout) :: fene_f
  integer, dimension(500,500), intent(inout) :: lj_force_PE
  integer, dimension(500,500), intent(inout) :: lj_force_ions
  integer, dimension(500,500), intent(inout) :: coulomb_f
  integer, dimension(500,500), intent(inout) :: Bond_dist
  integer :: i,j,k,l,m,n
  real*8 :: rsqr, fff, hz, inv_r2, inv_r6
  real*8, dimension(3) :: rjk,ff,rij

  do i=1, N_bond
    j=fene_list(i,1)
    k=fene_list(i,2)
    call rij_and_rr(rjk, rsqr, j, k)
    ff = -kfene*R0_2/(R0_2-rsqr)*rjk
    fff = abs(ff(3))
    hz = ( pos(j,3) + pos(k,3) ) / 2
    l = ceiling(hz/(Lz/500))
    m = ceiling(fff/(50./500))
    if (l>0 .and. l<=500 .and. m>0 .and. m<500) then
      fene_f(l,m) = fene_f(l,m) + 1
    end if
    n = ceiling(sqrt(rsqr)/(sqrt(R0_2)/500))
    if (l>0 .and. l<=500 .and. n>0 .and. n<500) then
      Bond_dist(l,n) = Bond_dist(l,n) + 1
    end if
  end do

  do k=1,npair1
    i=lj_pair_list(k,1)
    j=lj_pair_list(k,2)
    call rij_and_rr(rij,rsqr,i,j)
    if (rsqr<rcl*rcl) then
      inv_r2=1/rsqr
      inv_r6=inv_r2*inv_r2*inv_r2
      ff = inv_r2*inv_r6*(inv_r6-0.5)*rij
      fff = abs(ff(3))
      hz = pos(i,3)
      l = ceiling(hz/(Lz/500))
      m = ceiling(fff/(50./500))
      if (l>0 .and. l<=500 .and. m>0 .and. m<500) then
        if (i<=Npe .and. j<=Npe) then
          lj_force_PE(l,m) = lj_force_PE(l,m) + 1
        elseif(i<=Npe .and. j>Npe) then
          lj_force_ions(l,m) = lj_force_ions(l,m) + 1
        end if
      end if
    end if
  end do

  if ( qq /= 0 .and. mod(step,multistep) == 0 ) then    
    acc_c = 0
    posq(:,:) = pos(charge,:)
    call SPME_Ewald
    call real_space
    do i = 1, Npe
      ff = acc_c(i,:)
      fff = ff(3)
      hz = pos(i,3)
      l = ceiling(hz/(Lz/500))
      m = nint(fff/(100./500)) + 250
      if (l>0 .and. l<=500 .and. m>0 .and. m<500) then
        coulomb_f(l,m) = coulomb_f(l,m) + 1
      end if
    end do
  end if

end subroutine force_distribution

!
!The following program are the SPME program.
!

subroutine SPME_Ewald
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
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
  
  call splcof(SPx, dSPx, iqmap, ordr(1), 1)
  call splcof(SPy, dSPy, jqmap, ordr(2), 2)
  call splcof(SPz, dSPz, kqmap, ordr(3), 3)

  call MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)

  call pmeOrthoConvBC

  call IFrc(acc_f, SPx, SPy, SPz, dSPx, dSPy, dSPz, iqmap, jqmap, kqmap)

  call ZeroForce(acc_f)

  acc_c(charge,:)=acc_c(charge,:)+acc_f*lb/Beta   !Nq和NN不一样的时候这是不对的

end subroutine SPME_Ewald



subroutine splcof(SP, dSP, qmap, ord, xyz)
  !-----------------------------------------!
  !ord指定数组大小可以吗，回去查查, 可以的
  !不能在定义时候初始化
  !-----------------------------------------!
  use global_variables
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
  
  uu = ng(xyz)/gdim(xyz)*posq(:,xyz)            !化到标定后的分数坐标
  uuc = ceiling(uu)                           !ceiling can also be used for array
  
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


subroutine MapCharges(SPx, SPy, SPz, iqmap, jqmap, kqmap)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  real*8,  dimension(Nq,ordr(1)), intent(in) :: SPx
  real*8,  dimension(Nq,ordr(2)), intent(in) :: SPy
  real*8,  dimension(Nq,ordr(3)), intent(in) :: SPz
  integer, dimension(Nq,ordr(1)), intent(in) :: iqmap
  integer, dimension(Nq,ordr(2)), intent(in) :: jqmap
  integer, dimension(Nq,ordr(3)), intent(in) :: kqmap
  integer :: h,i,j,k,ii,jj,kk,t,ng1o,ng2o,ng3o,ih,jh,kh
  real*8  :: dqx,dqxy,dd,qmij(ordr(1),ordr(2))
  real*8  :: SPxh(ordr(1),1),SPyh(1,ordr(2)),SPzh(ordr(3))
  integer :: imp(ordr(1)),jmp(ordr(2)),kmp(ordr(3))
  real*8, dimension(ordr(1),ordr(2),ordr(3)) :: qmijk

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
    SPxh(:,1)=SPx(h,:)
    SPyh(1,:)=SPy(h,:)
    SPzh=SPz(h,:)
    qmij=matmul(SPxh,SPyh)
    do i=1,ordr(3)
      qmijk(:,:,i)=SPzh(i)*qmij
    end do
    Q_PME(imp,jmp,kmp)=Q_PME(imp,jmp,kmp)+qmijk
  end do

end subroutine Mapcharges


subroutine pmeOrthoConvBC
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
  implicit none
  include "fftw3.f90"
  complex ( kind = 8 ), dimension(ng(1),ng(2),ng(3)) :: fQ_PME
  complex ( kind = 8 ), dimension(ng(1),ng(2),ng(3)) :: fQBC
  integer i,j,k
  integer ( kind = 8 ) plan_backward
  integer ( kind = 8 ) plan_forward
  integer ( kind = 8 ) plan
  
  call dfftw_plan_dft_3d_ ( plan_forward, ng(1), ng(2), ng(3), Q_PME, fQ_PME, FFTW_FORWARD, FFTW_Estimate ) 

  call dfftw_execute_ ( plan_forward )
  
  call dfftw_destroy_plan_ ( plan_forward )
  
  fQ_PME(1,1,1)=0         

  fQBC=fQ_PME*BC_PME

  call dfftw_plan_dft_3d_ ( plan_backward, ng(1), ng(2), ng(3), fQBC, U_PME, FFTW_BACKWARD, FFTW_Estimate )

  call dfftw_execute_ ( plan_backward ) 

  call dfftw_destroy_plan_ ( plan_backward )

end subroutine pmeOrthoConvBC



subroutine IFrc(acc_f, Bx, By, Bz, dBx, dBy, dBz, iqm, jqm, kqm)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
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
  integer :: h,i,ng1o,ng2o,ng3o,ih,jh,kh
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
    f1 = -sum(dqmijkx*Ureg)     !the coefficient of the inverse fft in fftw
    f2 = -sum(dqmijky*Ureg)
    f3 = -sum(dqmijkz*Ureg)

    acc_f(h,1) = real(f1) * ng(1)/gdim(1) / ( ng(1)*ng(2)*ng(3) )
    acc_f(h,2) = real(f2) * ng(2)/gdim(2) / ( ng(1)*ng(2)*ng(3) )
    acc_f(h,3) = real(f3) * ng(3)/gdim(3) / ( ng(1)*ng(2)*ng(3) )

  end do

end subroutine IFrc


subroutine ZeroForce(acc_f)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
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


subroutine bspln_coeffs
  !----------------------------------------!
  !Calculate the coefficients of the b-spline basis function.
  !2010, Gradimir V. Milovanovic, Applied Mathematics Letters
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: h,i,j,k,m,g,ii
  real*8, allocatable, dimension(:,:,:) :: a
  if ( allocated(bspln_cof) ) deallocate(bspln_cof)
  if ( allocated(dspln_cof) ) deallocate(dspln_cof)
  allocate(bspln_cof(maxval(ordr),maxval(ordr),3))
  allocate(dspln_cof(maxval(ordr),maxval(ordr),3))
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

end subroutine bspln_coeffs


subroutine PME_BC
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  complex (kind=8) :: bx(ng(1)), by(ng(2)), bz(ng(3))
  real*8 :: mx(ng(1)), my(ng(2)), mz(ng(3))
  real*8 :: ex(ng(1)), ey(ng(2)), ez(ng(3))
  real*8 dx, dxy, dxyz
  integer :: i,j,k
  
  if ( allocated(BC_PME) ) deallocate(BC_PME)
  if ( allocated(Q_PME)  ) deallocate(Q_PME)
  if ( allocated(U_PME)  ) deallocate(U_PME)

  allocate( BC_PME(Kmax1, Kmax2, Kmax3) )
  allocate( Q_PME(Kmax1, Kmax2, Kmax3)  )
  allocate( U_PME(Kmax1, Kmax2, Kmax3)  )

  BC_PME = 0
  Q_PME  = 0
  U_PME  = 0
  
  call pmeOrthoTabBC(bx, mx, ex, 1)
  call pmeOrthoTabBC(by, my, ey, 2)
  call pmeOrthoTabBC(bz, mz, ez, 3)
  
  do i=1,ng(1)
    dx = mx(i)*mx(i)
    do j=1,ng(2)
      dxy = dx + my(j)*my(j)
      do k = 1,ng(3)
        dxyz = dxy + mz(k)*mz(k)
        if (i+j+k>3) then                         !避免除以0
          BC_PME(i,j,k) = bx(i)*conjg(bx(i)) * by(j)*conjg(by(j)) *        &
                          bz(k)*conjg(bz(k)) * ex(i)*ey(j)*ez(k) / dxyz
        end if
      end do
    end do
  end do
end subroutine PME_BC


subroutine pmeOrthoTabBC(bv, mv, ev, xyz)
  !-----------------------------------------!
  !
  !-----------------------------------------!
  use global_variables
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
    if (i-1 <= ng(xyz)/2) then            !这个m的取值范围有些讲究
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
      bsp = bsp + bspln(j+1.D0,ord) * &
            cmplx( cos(2*pi*i*j/ng(xyz)), sin(2*pi*i*j/ng(xyz)), 8 )
    end do
    if (abs(bsp)<1e-10) then
      bsp=1e15
    end if
    bv(i+1)=cmplx( cos(2*pi*i*(ord-1)/ng(xyz)), sin(2*pi*i*(ord-1)/ng(xyz)), &
                   8 ) / bsp !complex number can be divided
  end do
  
end subroutine pmeOrthoTabBC


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


end module compute_acceleration



