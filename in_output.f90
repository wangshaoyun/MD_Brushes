module input_output
implicit none 

save
  !ZhangFen
  real*8,  allocatable, dimension(:,:), private :: phi_branch
  real*8,  allocatable, dimension(:,:), private :: alpha_stem
  real*8,  allocatable, dimension(:,:), private :: alpha_branch
  real*8,  allocatable, dimension(:,:), private :: alpha_end
  !Histogram
  real*8,  allocatable, dimension(:,:), private :: phi_tot
  real*8,  allocatable, dimension(:,:), private :: phi_l
  real*8,  allocatable, dimension(:,:), private :: phi_le
  real*8,  allocatable, dimension(:,:), private :: phi_s
  real*8,  allocatable, dimension(:,:), private :: phi_sb
  real*8,  allocatable, dimension(:,:), private :: phi_se
  real*8,  allocatable, dimension(:,:), private :: phi_a
  real*8,  allocatable, dimension(:,:), private :: phi_i
  real*8,  allocatable, dimension(:,:), private :: phi_q
  real*8,  allocatable, dimension(:,:), private :: delta_angle1
  real*8,  allocatable, dimension(:,:), private :: delta_angle2
  real*8,  allocatable, dimension(:,:), private :: delta_angle3
  real*8,  allocatable, dimension(:,:), private :: theta_l
  real*8,  allocatable, dimension(:,:), private :: theta_lz
  real*8,  allocatable, dimension(:,:), private :: theta_ssl
  real*8,  allocatable, dimension(:,:), private :: theta_sslz
  real*8,  allocatable, dimension(:,:), private :: theta_sbl
  real*8,  allocatable, dimension(:,:), private :: theta_sblz
  real*8,  allocatable, dimension(:,:), private :: theta_bez
  real*8,  allocatable, dimension(:,:), private :: force_l
  real*8,  allocatable, dimension(:,:), private :: force_sy
  real*8,  allocatable, dimension(:,:), private :: force_sn
  real*8,  allocatable, dimension(:,:), private :: force_so
  real*8,  allocatable, dimension(:,:), private :: force_l1
  real*8,  allocatable, dimension(:,:), private :: force_sy1
  real*8,  allocatable, dimension(:,:), private :: force_sn1
  real*8,  allocatable, dimension(:,:), private :: force_so1
  real*8,  allocatable, dimension(:,:), private :: linear_h_dist
  real*8,  allocatable, dimension(:,:), private :: linear_Rg2_dist
  real*8,  allocatable, dimension(:,:), private :: linear_Rgz2_dist
  real*8,  allocatable, dimension(:,:), private :: linear_Rgxy2_dist
  real*8,  allocatable, dimension(:,:), private :: star_h_dist
  real*8,  allocatable, dimension(:,:), private :: star_Rg2_dist
  real*8,  allocatable, dimension(:,:), private :: star_Rgz2_dist
  real*8,  allocatable, dimension(:,:), private :: star_Rgxy2_dist
  real*8,  allocatable, dimension(:,:), private :: star_end_h_dist
  integer,  allocatable, dimension(:,:), private :: fene_f
  integer,  allocatable, dimension(:,:), private :: lj_force_PE
  integer,  allocatable, dimension(:,:), private :: lj_force_ions
  integer,  allocatable, dimension(:,:), private :: coulomb_f
  integer,  allocatable, dimension(:,:), private :: Bond_dist
  integer, allocatable, dimension(:,:), private :: phi_zx
  integer, allocatable, dimension(:,:), private :: phi_xy
  integer, allocatable, dimension(:,:), private :: phi_yz
  integer, allocatable, dimension(:,:), private :: phi_lzx
  integer, allocatable, dimension(:,:), private :: phi_lxy
  integer, allocatable, dimension(:,:), private :: phi_lyz
  integer, allocatable, dimension(:,:), private :: phi_lezx
  integer, allocatable, dimension(:,:), private :: phi_lexy
  integer, allocatable, dimension(:,:), private :: phi_leyz
  integer, allocatable, dimension(:,:), private :: phi_szx
  integer, allocatable, dimension(:,:), private :: phi_sxy
  integer, allocatable, dimension(:,:), private :: phi_syz
  integer, allocatable, dimension(:,:), private :: phi_sbzx
  integer, allocatable, dimension(:,:), private :: phi_sbxy
  integer, allocatable, dimension(:,:), private :: phi_sbyz
  integer, allocatable, dimension(:,:), private :: phi_sezx
  integer, allocatable, dimension(:,:), private :: phi_sexy
  integer, allocatable, dimension(:,:), private :: phi_seyz
  integer, allocatable, dimension(:,:), private :: phi_azx
  integer, allocatable, dimension(:,:), private :: phi_axy
  integer, allocatable, dimension(:,:), private :: phi_ayz
  integer, allocatable, dimension(:,:), private :: phi_izx
  integer, allocatable, dimension(:,:), private :: phi_ixy
  integer, allocatable, dimension(:,:), private :: phi_iyz
  integer, allocatable, dimension(:,:), private :: phi_qzx
  integer, allocatable, dimension(:,:), private :: phi_qxy
  integer, allocatable, dimension(:,:), private :: phi_qyz
  !Physical quantities of height
  real*8, private :: h_avg
  real*8, private :: hl_avg
  real*8, private :: hl_max
  real*8, private :: hl_end
  real*8, private :: hs_avg
  real*8, private :: hs_max
  real*8, private :: hs_end
  real*8, private :: hs_branch
  !Physical quantities of Re
  real*8, private :: Re_l
  real*8, private :: Re_lz
  real*8, private :: Re_s
  real*8, private :: Re_sz
  real*8, private :: Re_ss
  real*8, private :: Re_ssz
  real*8, private :: Re_sb
  real*8, private :: Re_sbz
  !Physical quantities of Rg
  real*8, private :: Rg_l
  real*8, private :: Rg_lz
  real*8, private :: Rg_s
  real*8, private :: Rg_sz
  real*8, private :: Rg_ss
  real*8, private :: Rg_ssz
  real*8, private :: Rg_sb
  real*8, private :: Rg_sbz
  !Physical quantities of energy
  real*8, private :: kenetic_energy
  real*8, private :: R_up1
  real*8, private :: R_down1
  real*8, private :: R_I1
  real*8, private :: R_up2
  real*8, private :: R_down2
  real*8, private :: R_I2
  real*8, private :: R_up3
  real*8, private :: R_down3
  real*8, private :: R_I3
contains

subroutine initialize_parameters
  !--------------------------------------!
  !Initialize system parameters
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
  implicit none

  call read_data

  call data_operation

  call data_allocate

  call write_system_parameters

end subroutine initialize_parameters


subroutine read_data
  !--------------------------------------!
  !read system parameters
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
  logical alive
  !
  !read system data
  open(10,file='./system_data.txt')
    read(10,*) uniform_or_random
    read(10,*) Lz
    read(10,*) Z_empty  
    read(10,*) sigmag   
    read(10,*) Beta       
    read(10,*) qq
    read(10,*) qqi
    read(10,*) ion_ratio       
    read(10,*) arm
    read(10,*) Nma
    read(10,*) Nga
    read(10,*) Nml
    read(10,*) Ngl
    read(10,*) man
    read(10,*) R_bond
    read(10,*) StepNum0           
    read(10,*) StepNum            
    read(10,*) DeltaStep1
    read(10,*) DeltaStep2     
    read(10,*) DeltaStep3         
    read(10,*) multistep  
    read(10,*) dt
  close(10)
end subroutine read_data


subroutine data_operation
  !--------------------------------------!
  !Initialize system parameters
  !and judge whether restarted or continue
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
  logical alive
  integer :: i
  integer :: arm_q          !Charged arms of star brushes
  integer :: Charge_ions    !电量而不是粒子数
  !
  !the situation of the arm whether is charge or not
  allocate( charged_arm(arm) )
  charged_arm = 0
  if ( Nga /= 0 ) then
    open( 100, file='q_arm.txt' )
      do i = 1, arm
        read( 100, * ) charged_arm(i)
      end do
    close( 100 )
  end if
  !
  !The total arms that are charged
  arm_q = sum( charged_arm )
  !
  !The total monomers of star brushes
  Nta = ( Nma*arm + 1 ) * Nga
  !
  !The total monomers of linear brushes
  Ntl = Nml * Ngl
  !
  !The total monomers of star and linear brushes
  Npe = Nta + Ntl
  !
  !The total chains that are anchored on the plate
  N_anchor = Nga + Ngl
  !
  !The total particles NN and total particles that are
  !charged Nq in the system
  if ( abs(qq) == 0 ) then
    Nq_PE = 0
  else
    if ( man /= 0 ) then
      Nq_PE = (Ntl - Ngl + arm_q * Nma * Nga) / man
    else
      Nq_PE = 0
    end if
  end if
  if ( abs(qqi) == 0 ) then
    Nq_salt_ions = 0
  else
    Charge_ions  = nint( ion_ratio * Nq_PE * nint(abs(qq)) )
    Nq_salt_ions = Charge_ions / nint(abs(qqi))
  end if
  Nq = Nq_PE * ( nint(abs(qq))+1 ) + Nq_salt_ions * ( nint(abs(qqi)) + 1 )
  NN = Npe + Nq_PE * nint(abs(qq)) + Nq_salt_ions * ( nint(abs(qqi)) + 1 )
  !
  !System size
  Lx = sqrt( N_anchor / sigmag )
  Ly = Lx
  Z_empty = ( Lz + Lx * Z_empty ) / Lz
  !
  !whether continue or restart
  Inquire( file='start_time.txt', exist=alive )
  if (alive) then
    open(11,file='./start_time.txt')
      read(11,*) restart_or_continue
    close(11)
  else
    restart_or_continue = 0
  end if
end subroutine data_operation


subroutine data_allocate
  !--------------------------------------!
  !Allocate pos, vel, acc and histogram arrays
  !Initialize histogram arrays
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
  !position, velocity, and acceleration
  allocate( pos(NN,4) )
  allocate( vel(NN,3) )
  allocate( acc(NN,3) )
  pos = 0
  vel = 0
  acc = 0
  !
  !ZhangFen
  allocate( phi_branch(SizeHist,2)      )
  allocate( alpha_stem(SizeHist,2)      )
  allocate( alpha_branch(SizeHist,2)    )
  allocate( alpha_end(SizeHist,2)    )
  !
  !Histogram
  allocate( phi_tot(SizeHist,2)         )
  allocate( phi_l(SizeHist,2)           )
  allocate( phi_le(SizeHist,2)          )
  allocate( phi_s(SizeHist,2)           )
  allocate( phi_sb(SizeHist,2)          )
  allocate( phi_se(SizeHist,2)          )
  allocate( phi_a(SizeHist,2)           )
  allocate( phi_i(SizeHist,2)           )
  allocate( phi_q(SizeHist,2)           )
  allocate( delta_angle1(SizeHist,2)    )
  allocate( delta_angle2(SizeHist,2)    )
  allocate( delta_angle3(SizeHist,2)    )
  allocate( theta_l(SizeHist,2)         )
  allocate( theta_lz(SizeHist,2)        )
  allocate( theta_ssl(SizeHist,2)       )
  allocate( theta_sslz(SizeHist,2)      )
  allocate( theta_sbl(SizeHist,2)       )
  allocate( theta_sblz(SizeHist,2)      )
  allocate( theta_bez(SizeHist,2)       )
  allocate( force_l(Nml-1,2)            )
  allocate( force_sy(Nma*2,2)           )
  allocate( force_sn(Nma*2,2)           )
  allocate( force_so(Nma*2,2)           )
  allocate( force_l1(SizeHist,2)        )
  allocate( force_sy1(SizeHist,2)       )
  allocate( force_sn1(SizeHist,2)       )
  allocate( force_so1(SizeHist,2)       )
  allocate( linear_h_dist(SizeHist,2)   )
  allocate( linear_Rg2_dist((floor(Lz/2))**2,2) )
  allocate( linear_Rgz2_dist((floor(Lz/2))**2,2))
  allocate(linear_Rgxy2_dist((floor(Lz/2))**2,2))
  allocate( star_h_dist(SizeHist,2)     )
  allocate( star_Rg2_dist((floor(Lz/2))**2,2)   )
  allocate( star_Rgz2_dist((floor(Lz/2))**2,2)  )
  allocate( star_Rgxy2_dist((floor(Lz/2))**2,2) )
  allocate( star_end_h_dist(SizeHist,2) )
  allocate( fene_f(SizeHist,SizeHist) )
  allocate( lj_force_PE(SizeHist,SizeHist) )
  allocate( lj_force_ions(SizeHist,SizeHist) )
  allocate( coulomb_f(SizeHist,SizeHist) )
  allocate( Bond_dist(SizeHist,SizeHist) )
  allocate( phi_zx(SizeHist,SizeHist)   )
  allocate( phi_xy(SizeHist,SizeHist)   )
  allocate( phi_yz(SizeHist,SizeHist)   )
  allocate( phi_lzx(SizeHist,SizeHist)  )
  allocate( phi_lxy(SizeHist,SizeHist)  )
  allocate( phi_lyz(SizeHist,SizeHist)  )
  allocate( phi_lezx(SizeHist,SizeHist) )
  allocate( phi_lexy(SizeHist,SizeHist) )
  allocate( phi_leyz(SizeHist,SizeHist) )
  allocate( phi_szx(SizeHist,SizeHist)  )
  allocate( phi_sxy(SizeHist,SizeHist)  )
  allocate( phi_syz(SizeHist,SizeHist)  )
  allocate( phi_sbzx(SizeHist,SizeHist) )
  allocate( phi_sbxy(SizeHist,SizeHist) )
  allocate( phi_sbyz(SizeHist,SizeHist) )
  allocate( phi_sezx(SizeHist,SizeHist) )
  allocate( phi_sexy(SizeHist,SizeHist) )
  allocate( phi_seyz(SizeHist,SizeHist) )
  allocate( phi_azx(SizeHist,SizeHist)  )
  allocate( phi_axy(SizeHist,SizeHist)  )
  allocate( phi_ayz(SizeHist,SizeHist)  )
  allocate( phi_izx(SizeHist,SizeHist)  )
  allocate( phi_ixy(SizeHist,SizeHist)  )
  allocate( phi_iyz(SizeHist,SizeHist)  )
  allocate( phi_qzx(SizeHist,SizeHist)  )
  allocate( phi_qxy(SizeHist,SizeHist)  )
  allocate( phi_qyz(SizeHist,SizeHist)  )
  !initialize histogram
  phi_tot      = 0
  phi_l        = 0
  phi_le       = 0
  phi_s        = 0
  phi_sb       = 0
  phi_se       = 0
  phi_a        = 0
  phi_i        = 0
  phi_q        = 0
  delta_angle1 = 0
  delta_angle2 = 0
  delta_angle3 = 0
  theta_l      = 0
  theta_lz     = 0
  theta_ssl    = 0
  theta_sslz   = 0
  theta_sbl    = 0
  theta_sblz   = 0
  theta_bez    = 0
  force_l      = 0
  force_sy     = 0
  force_sn     = 0
  force_so     = 0
  force_l1     = 0
  force_sy1    = 0
  force_sn1    = 0
  force_so1    = 0
  linear_h_dist = 0 
  linear_Rg2_dist = 0  
  linear_Rgz2_dist = 0
  linear_Rgxy2_dist = 0
  star_h_dist = 0
  star_Rg2_dist = 0
  star_Rgz2_dist = 0
  star_Rgxy2_dist = 0
  star_end_h_dist = 0
  fene_f = 0
  lj_force_PE = 0
  lj_force_ions = 0
  coulomb_f = 0
  Bond_dist = 0
  phi_zx       = 0
  phi_xy       = 0
  phi_yz       = 0
  phi_lzx      = 0
  phi_lxy      = 0
  phi_lyz      = 0
  phi_lezx     = 0
  phi_lexy     = 0
  phi_leyz     = 0
  phi_szx      = 0
  phi_sxy      = 0
  phi_syz      = 0
  phi_sbzx     = 0
  phi_sbxy     = 0
  phi_sbyz     = 0
  phi_sezx     = 0
  phi_sexy     = 0
  phi_seyz     = 0
  phi_azx      = 0
  phi_axy      = 0
  phi_ayz      = 0
  phi_izx      = 0
  phi_ixy      = 0
  phi_iyz      = 0
  phi_qzx      = 0
  phi_qxy      = 0
  phi_qyz      = 0
  
  h_avg        = 0
end subroutine data_allocate


subroutine write_system_parameters
  !--------------------------------------!
  !Write system parameters to screen.
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
  !write data
  write(*,*) '******************system_data***********************'
  write(*,*) 'Arms of star brushes             Arm:', arm
  write(*,*) 'Monomers of each arm             Nma:', Nma
  write(*,*) 'Number of grafted star chains    Nga:', Nga
  write(*,*) 'Monomers of each linear chain    Nml:', Nml
  write(*,*) 'Number of grafted linear chains  Ngl:', Ngl
  write(*,*) 'total particles,                  NN:', NN
  write(*,*) 'total charged particles,          Nq:', Nq
  write(*,*) 'total anchored particles,   N_anchor:', N_anchor
  write(*,*) 'total brushes particles,         Npe:', Npe
  write(*,*) 'total linear brushes particles,  Ntl:', Ntl
  write(*,*) 'total star brushes particles,    Nta:', Nta
  write(*,*) 'Length of cell in x direction,    Lx:', Lx
  write(*,*) 'Length of cell in y direction,    Ly:', Ly
  write(*,*) 'Length of cell in z direction,    Lz:', Lz
  write(*,*) 'Grafting density,             sigmag:', sigmag
  write(*,*) 'Beta=1/kT,                      Beta:', Beta
  write(*,*) '****************************************************'
  
  write(*,*) '******************running_steps*********************'
  write(*,*) 'Preheating steps            StepNum0:', StepNum0
  write(*,*) 'Running steps                StepNum:', StepNum
  write(*,*) 'Total steps         StepNum0+StepNum:', (StepNum0+StepNum)
  write(*,*) 'Step inteval              DeltaStep1:', DeltaStep1
  write(*,*) 'Step inteval              DeltaStep2:', DeltaStep2
  write(*,*) 'Step inteval              DeltaStep3:', DeltaStep3
  write(*,*) 'Multisteps of coulomb      MultiStep:', MultiStep
  write(*,*) 'Distance of each move             dt:', dt
  write(*,*) '****************************************************'

end subroutine write_system_parameters
  
  
subroutine continue_read_data(l)
  !--------------------------------------!
  !When the program is break off due to power off,
  !we need to continue the program from the former data.
  !This subroutine is used to read data from file to continue
  !calculating.
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
  integer, intent(out) :: l
  integer :: i,j
  real*8, dimension(SizeHist,8)  :: alpha_phi
  real*8, dimension(SizeHist,10) :: phi
  real*8, dimension(SizeHist,8)  :: theta
  real*8, dimension(2*Nma,4)     :: force
  real*8, dimension(SizeHist,4)  :: force1
  real*8, dimension(SizeHist,4)  :: delta_angle
  real*8, dimension((floor(Lz/2))**2,7) :: Rg_dist
  real*8, dimension(SizeHist,4)  :: h_dist
  !
  !read pos and vel
  open(20,file='./data/pos1.txt')
  open(21,file='./data/vel1.txt')
    read(20,*) ((pos(i,j),j=1,4),i=1,NN)
    read(21,*) ((vel(i,j),j=1,3),i=1,NN)
  close(20)
  close(21)
  !
  !read steps and total_time of the former calculation.
  open(19,file='./start_time.txt')
    read(19,*) 
    read(19,*) l
    read(19,*) total_time
    l = l + 1
  close(19)
  !
  !read physical quantities
  phi   = 0
  theta = 0
  force = 0
  delta_angle=0
  open(19,file='./data/alpha_phi.txt')
  open(20,file='./data/phi.txt')
  open(21,file='./data/theta.txt')
  open(22,file='./data/force_liear.txt')
  open(23,file='./data/force_star.txt')
  open(24,file='./data/delta_angle.txt')
  open(25,file='./data/force_liear1.txt')
  open(26,file='./data/force_star1.txt')
  open(27,file='./data/Rg_dist.txt')
  open(28,file='./data/h_dist.txt')
    read(19,*) ((alpha_phi(i,j),j=1,8),i=1,SizeHist)
      phi_branch(:,2)   = alpha_phi(:,2)
      alpha_stem(:,2)   = alpha_phi(:,4)
      alpha_branch(:,2) = alpha_phi(:,6)
      alpha_end(:,2)    = alpha_phi(:,8)
    read(20,*) ((phi(i,j),j=1,10),i=1,SizeHist)
      phi_tot(:,2)= phi(:,2)
      phi_l(:,2)  = phi(:,3)
      phi_le(:,2) = phi(:,4)
      phi_s(:,2)  = phi(:,5)
      phi_sb(:,2) = phi(:,6)
      phi_se(:,2) = phi(:,7)
      phi_a(:,2)  = phi(:,8)
      phi_i(:,2)  = phi(:,9)
      phi_q(:,2)  = phi(:,10)
    read(21,*) ((theta(i,j),j=1,8),i=1,SizeHist)
      theta_l(:,2)=theta(:,2)
      theta_lz(:,2)=theta(:,3)
      theta_ssl(:,2)=theta(:,4)
      theta_sslz(:,2)=theta(:,5)
      theta_sbl(:,2)=theta(:,6)
      theta_sblz(:,2)=theta(:,7)
      theta_bez(:,2)=theta(:,8)
    read(22,*) ((force_l(i,j),j=1,2),i=1,Nml-1)
    read(23,*) ((force(i,j),j=1,4),i=1,2*Nma)
      force_sy(:,2)=force(:,2)
      force_sn(:,2)=force(:,3)
      force_so(:,2)=force(:,4)
    read(24,*) ((delta_angle(i,j),j=1,4),i=1,SizeHist)
      delta_angle1(:,2)=delta_angle(:,2)
      delta_angle2(:,2)=delta_angle(:,3)
      delta_angle3(:,2)=delta_angle(:,4)
    read(25,*) ((force_l1(i,j),j=1,2),i=1,SizeHist)
    read(26,*) ((force1(i,j),j=1,4),i=1,SizeHist)
      force_sy1(:,2)=force1(:,2)
      force_sn1(:,2)=force1(:,3)
      force_so1(:,2)=force1(:,4)
    read(27,*) ((Rg_dist(i,j),j=1,7),i=1,(floor(Lz/2))**2)
      linear_Rg2_dist(:,2)=Rg_dist(:,2)
      linear_Rgz2_dist(:,2)=Rg_dist(:,3)
      linear_Rgxy2_dist(:,2)=Rg_dist(:,4)
      star_Rg2_dist(:,2)=Rg_dist(:,5)
      star_Rgz2_dist(:,2)=Rg_dist(:,6)
      star_Rgxy2_dist(:,2)=Rg_dist(:,7)
    read(28,*) ((h_dist(i,j),j=1,4),i=1,SizeHist)
      linear_h_dist(:,2)=h_dist(:,2)
      star_h_dist(:,2)=h_dist(:,3)
      star_end_h_dist(:,2)=h_dist(:,4)
  close(28)
  close(27)
  close(26)
  close(25)
  close(24)
  close(23)
  close(22)
  close(21)
  close(20)
  close(19)
  !
  !read histogram
  open(21,file='./data/phi_2d11.txt')
  open(22,file='./data/phi_2d12.txt')
  open(23,file='./data/phi_2d13.txt')
  open(24,file='./data/phi_2d21.txt')
  open(25,file='./data/phi_2d22.txt')
  open(26,file='./data/phi_2d23.txt')
  open(27,file='./data/phi_2d31.txt')
  open(28,file='./data/phi_2d32.txt')
  open(29,file='./data/phi_2d33.txt')
  open(30,file='./data/phi_2d41.txt')
  open(31,file='./data/phi_2d42.txt')
  open(32,file='./data/phi_2d43.txt')
  open(33,file='./data/phi_2d51.txt')
  open(34,file='./data/phi_2d52.txt')
  open(35,file='./data/phi_2d53.txt')
  open(36,file='./data/phi_2d61.txt')
  open(37,file='./data/phi_2d62.txt')
  open(38,file='./data/phi_2d63.txt')
  open(39,file='./data/phi_2d71.txt')
  open(40,file='./data/phi_2d72.txt')
  open(41,file='./data/phi_2d73.txt')
  open(42,file='./data/phi_2d81.txt')
  open(43,file='./data/phi_2d82.txt')
  open(44,file='./data/phi_2d83.txt')
  open(45,file='./data/phi_2d91.txt')
  open(46,file='./data/phi_2d92.txt')
  open(47,file='./data/phi_2d93.txt')
  open(48,file='./data/fene_f.txt')
  open(49,file='./data/lj_force_PE.txt')
  open(50,file='./data/lj_force_ions.txt')
  open(51,file='./data/coulomb_f.txt')
  open(52,file='./data/Bond_dist.txt')
    read(21,*) ((phi_zx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(22,*) ((phi_xy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(23,*) ((phi_yz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(24,*) ((phi_lzx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(25,*) ((phi_lxy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(26,*) ((phi_lyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(27,*) ((phi_lezx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(28,*) ((phi_lexy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(29,*) ((phi_leyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(30,*) ((phi_szx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(31,*) ((phi_sxy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(32,*) ((phi_syz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(33,*) ((phi_sbzx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(34,*) ((phi_sbxy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(35,*) ((phi_sbyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(36,*) ((phi_sezx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(37,*) ((phi_sexy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(38,*) ((phi_seyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(39,*) ((phi_azx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(40,*) ((phi_axy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(41,*) ((phi_ayz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(42,*) ((phi_izx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(43,*) ((phi_ixy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(44,*) ((phi_iyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(45,*) ((phi_qzx(i,j),j=1,SizeHist),i=1,SizeHist)
    read(46,*) ((phi_qxy(i,j),j=1,SizeHist),i=1,SizeHist)
    read(47,*) ((phi_qyz(i,j),j=1,SizeHist),i=1,SizeHist)
    read(48,*) ((fene_f(i,j),j=1,SizeHist),i=1,SizeHist)
    read(49,*) ((lj_force_PE(i,j),j=1,SizeHist),i=1,SizeHist)
    read(50,*) ((lj_force_ions(i,j),j=1,SizeHist),i=1,SizeHist)
    read(51,*) ((coulomb_f(i,j),j=1,SizeHist),i=1,SizeHist)
    read(52,*) ((Bond_dist(i,j),j=1,SizeHist),i=1,SizeHist)
  close(52)
  close(51)
  close(50)
  close(49)
  close(48)
  close(47)
  close(46)
  close(45)
  close(44)
  close(43)
  close(42)
  close(41)
  close(40)
  close(39)
  close(38)
  close(37)
  close(36)
  close(35)
  close(34)
  close(33)
  close(32)
  close(31)
  close(30)
  close(29)
  close(28)
  close(27)
  close(26)
  close(25)
  close(24)
  close(23)
  close(22)
  close(21)
end subroutine continue_read_data


subroutine height
  !--------------------------------------!
  !Initialize system parameters
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
  use compute_acceleration
  implicit none
  integer i,j,k,l,m,n,p,q,r,s,num_stretch,num_collapse
  real*8 :: rr,rsqr,rr1,rr2,rr3,maxh,Rg1,Rg1z,Rg2,Rg2z
  real*8 :: max_hs_end,min_hs_end,hs_avg_arm,hb_avg,he_avg,max_hs,min_hs
  real*8, dimension(3) :: rij
  real*8, dimension(3) :: rij1,rij2,rij3,f1
  real*8 :: h1

  !------------up and down-------------!
  R_up1   = 0
  R_down1 = 0
  R_I1    = 0
  hb_avg = 0
  n      = arm * Nma + 1 
  do i = 1, Nga
    hb_avg = hb_avg + pos( (i-1)*n+Nma+1, 3 ) 
  end do
  hb_avg = hb_avg / Nga
  do i = 1, Nga
    max_hs = 0
    min_hs = Lz
    do j = 2, arm
      do k = 1, Nma
        m = (i-1)*n + (j-1)*Nma + 1 + k
        h1 = pos(m, 3 )
        if ( max_hs < h1 ) then
          max_hs = h1
        end if
        if ( min_hs > h1 ) then
          min_hs = h1
        end if
      end do
    end do
    if (min_hs > hb_avg) then
      R_up1 = R_up1 + 1
    elseif (max_hs < hb_avg) then
      R_down1 = R_down1 + 1
    else
      R_I1 = R_I1 + 1
    end if
  end do
  R_up1   = R_up1 / Nga
  R_down1 = R_down1 / Nga
  R_I1    = R_I1 / Nga

  R_up1   = 0
  R_down1 = 0
  R_I1    = 0
  he_avg = 0
  n      = arm * Nma + 1 
  do i = 1, Nga
    do j = 2, arm
      he_avg = he_avg + pos( (i-1)*n+Nma*j+1, 3 ) 
    end do
  end do
  he_avg = he_avg / Nga / (arm-1)
  do i = 1, Nga
    max_hs = 0
    min_hs = Lz
    do j = 2, arm
      m = (i-1)*n + j*Nma + 1
      h1 = pos(m, 3)
      if ( max_hs < h1 ) then
        max_hs = h1
      end if
      if ( min_hs > h1 ) then
        min_hs = h1
      end if
    end do
    if (min_hs > he_avg) then
      R_up2 = R_up2 + 1
    elseif (max_hs < he_avg) then
      R_down2 = R_down2 + 1
    else
      R_I2 = R_I2 + 1
    end if
  end do
  R_up2   = R_up2 / Nga
  R_down2 = R_down2 / Nga
  R_I2    = R_I2 / Nga
  !------------------------------------!

  h_avg=0
  hl_avg=0
  hl_max=0
  hl_end=0
  hs_avg=0
  hs_max=0
  hs_end=0
  hs_branch=0

  Re_l=0
  Re_lz=0
  Re_s=0
  Re_sz=0
  Re_ss=0
  Re_ssz=0
  Re_sb=0
  Re_sbz=0

  Rg_l=0
  Rg_lz=0
  Rg_s=0
  Rg_sz=0
  Rg_ss=0
  Rg_ssz=0
  Rg_sb=0
  Rg_sbz=0
  
  kenetic_energy=0
  R_up3=0
  R_down3=0
  R_I3=0

!!----------------------h_avg-------------------------!
  do i=1,Npe
    h_avg=h_avg+pos(i,3)
  end do
  h_avg=h_avg/Npe                                   
!!---------------hl_max,hl_end,hl_avg-----------------!
  do i=1,Ngl                                        
    maxh=0
    do j=1,Nml
      hl_avg=hl_avg+pos(Nta+(i-1)*Nml+j,3)
      if (maxh<pos(Nta+(i-1)*Nml+j,3)) then
        maxh=pos(Nta+(i-1)*Nml+j,3)
      end if
    end do
    hl_max=hl_max+maxh
    hl_end=hl_end+pos(Nta+i*Nml,3)
  end do
  hl_avg=hl_avg/Ngl/Nml
  hl_max=hl_max/Ngl
  hl_end=hl_end/Ngl
!!------------hs_max,hs_end,hs_branch,hs_avg-----------!
!!-------------force_sy,force_sn,force_sn--------------!
!!------------force_sy1,force_sn1,force_sn1------------!
!!--------delta_angle1,delta_angle2,delta_angle3-------!
!!------------------R_up3,R_down3,R_I3-----------------!
  do i=1,Nga
    maxh=0
    max_hs_end=0
    min_hs_end=Lz
    hs_avg_arm=0
    do k=1,Nma+1
      hs_avg_arm=hs_avg_arm+pos((i-1)*(arm*Nma+1)+k,3)
      if (maxh<pos((i-1)*(arm*Nma+1)+k,3)) then
        maxh=pos((i-1)*(arm*Nma+1)+k,3)                         !max height
      end if
    end do
    hs_branch=hs_branch+pos((i-1)*(arm*Nma+1)+Nma+1,3)          !branch point
    do j=2,arm
      do k=1,Nma
        hs_avg_arm=hs_avg_arm+pos((i-1)*(arm*Nma+1)+((j-1)*Nma+1)+k,3)  !average height
        if (maxh<pos((i-1)*(arm*Nma+1)+((j-1)*Nma+1)+k,3)) then
          maxh=pos((i-1)*(arm*Nma+1)+((j-1)*Nma+1)+k,3)         !max height
        end if
      end do
      hs_end=hs_end+pos((i-1)*(arm*Nma+1)+1+j*Nma,3)            !end of arms
      if (max_hs_end<pos((i-1)*(arm*Nma+1)+1+j*Nma,3)) then
        max_hs_end=pos((i-1)*(arm*Nma+1)+1+j*Nma,3)
      end if
      if (min_hs_end>pos((i-1)*(arm*Nma+1)+1+j*Nma,3)) then
        min_hs_end=pos((i-1)*(arm*Nma+1)+1+j*Nma,3) 
      end if
    end do
    hs_avg=hs_avg+hs_avg_arm
    hs_avg=hs_avg/(Nga*(arm*Nma+1))
    hs_max=hs_max+maxh
    hs_avg_arm=hs_avg_arm/(arm*Nma+1)
    if (max_hs_end<hs_avg) then
      R_down3=R_down3+1
      do k=2,Nma+1
        m=(i-1)*(arm*Nma+1)+k
        n=(i-1)*(arm*Nma+1)+k-1
        call rij_and_rr(rij1,rr1,m,n)
        f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
        f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
        f1=f1-lb*Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
        force_sn(k-1,2)=force_sn(k-1,2)+f1(3)/Nga
        force_sn1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=  &
              force_sn1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga
      end do
      do j=2, arm
        do k=1, Nma
          m=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k
          if (k==1) then
            n=(i-1)*(arm*Nma+1)+Nma+1
          else
            n=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k-1
          end if
          call rij_and_rr(rij1,rr1,m,n)
          f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
          f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
          f1=f1-lb/Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
          force_sn(k+Nma,2)=force_sn(k+Nma,2)+f1(3)/Nga/(arm-1)
          force_sn1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=  &
                force_sn1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga/(arm-1)
        end do
        m=(i-1)*(arm*Nma+1)+Nma+1
        n=(i-1)*(arm*Nma+1)+j*Nma+1
        call rij_and_rr(rij,rr,m,n)
        p=ceiling(acos(rij(3)/sqrt(rr))/(pi/SizeHist))
        if (p<=0 .or. p>SizeHist) then
          write(*,*)'error in delta_angle, exceeding (0,SizeHist)'
          cycle
        end if
        delta_angle2(p,2)=delta_angle2(p,2)+1
      end do
    elseif (min_hs_end>hs_avg) then
      R_up3=R_up3+1
      do k=2,Nma+1
        m=(i-1)*(arm*Nma+1)+k
        n=(i-1)*(arm*Nma+1)+k-1
        call rij_and_rr(rij1,rr1,m,n)
        f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
        f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
        f1=f1-lb/Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
        force_sy(k-1,2)=force_sy(k-1,2)+f1(3)/Nga
        force_sy1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=&
              force_sy1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga
      end do
      do j=2, arm
        do k=1, Nma
          m=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k
          if (k==1) then
            n=(i-1)*(arm*Nma+1)+Nma+1
          else
            n=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k-1
          end if
          call rij_and_rr(rij1,rr1,m,n)
          f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
          f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
          f1=f1-lb/Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
          force_sy(k+Nma,2)=force_sy(k+Nma,2)+f1(3)/Nga/(arm-1)
          force_sy1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=  &
                force_sy1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga/(arm-1)
        end do
        m=(i-1)*(arm*Nma+1)+Nma+1
        n=(i-1)*(arm*Nma+1)+j*Nma+1
        call rij_and_rr(rij,rr,m,n)
        p=ceiling(acos(rij(3)/sqrt(rr))/(pi/SizeHist))
        if (p<=0 .or. p>SizeHist) then
          write(*,*)'error in delta_angle, exceeding (0,SizeHist)'
          cycle
        end if
        delta_angle1(p,2)=delta_angle1(p,2)+1
      end do
    else
      R_I3=R_I3+1
      do k=2,Nma+1
        m=(i-1)*(arm*Nma+1)+k
        n=(i-1)*(arm*Nma+1)+k-1
        call rij_and_rr(rij1,rr1,m,n)
        f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
        f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
        f1=f1-lb/Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
        force_so(k-1,2)=force_so(k-1,2)+f1(3)/Nga
        force_so1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=  &
              force_so1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga
      end do
      do j=2, arm
        do k=1, Nma
          m=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k
          if (k==1) then
            n=(i-1)*(arm*Nma+1)+Nma+1
          else
            n=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k-1
          end if
          call rij_and_rr(rij1,rr1,m,n)
          f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
          f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
          f1=f1-lb/Beta*(1/rr1**1.5)*rij1*pos(m,4)*pos(n,4)
          force_so(k+Nma,2)=force_so(k+Nma,2)+f1(3)/Nga/(arm-1)
          force_so1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=&
                force_so1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Nga/(arm-1)
        end do
        m=(i-1)*(arm*Nma+1)+Nma+1
        n=(i-1)*(arm*Nma+1)+1+j*Nma
        call rij_and_rr(rij,rr,m,n)
        p=ceiling(acos(rij(3)/sqrt(rr))/(pi/SizeHist))
        if (p<=0 .or. p>SizeHist) then
          write(*,*)'error in delta_angle, exceeding (0,SizeHist)'
          cycle
        end if
        delta_angle3(p,2)=delta_angle3(p,2)+1
      end do
    end if    
  end do
  hs_max=hs_max/Nga
  hs_end=hs_end/(Nga*(arm-1))
  hs_branch=hs_branch/Nga
  R_up3=R_up3/Nga
  R_down3=R_down3/Nga
  R_I3=R_I3/Nga
!!---------------Re_l,Re_lz,Rg_l,Rg_lz------------------! 
  do i=1,Ngl
    m=Nta+(i-1)*Nml+1
    n=Nta+i*Nml
    call rij_and_rr(rij,rr,m,n)
    Re_l=Re_l+rr
    Re_lz=Re_lz+rij(3)*rij(3)
    Rg1=0
    Rg1z=0
    do j=1,Nml-1
      do k=j+1,Nml
        m=Nta+(i-1)*Nml+j
        n=Nta+(i-1)*Nml+k
        call rij_and_rr(rij,rr,m,n)
        Rg1=Rg1+rr
        Rg1z=Rg1z+rij(3)*rij(3)
      end do
      m=Nta+(i-1)*Nml+j+1
      n=Nta+(i-1)*Nml+j
      call rij_and_rr(rij1,rr1,m,n)
      f1=48*rr1**(-4)*(rr1**(-3)-0.5)*rij1
      f1=f1-kfene*R0_2/(R0_2-rr1)*rij1
      f1=f1-lb*Beta*(1/rr1)*rij1*pos(m,4)*pos(n,4)
      force_l(j,2)=force_l(j,2)+f1(3)/Ngl
      force_l1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)=&
            force_l1(ceiling((pos(m,3)+pos(n,3))/2/(Lz/SizeHist)),2)+f1(3)/Ngl
    end do
    Rg_l=Rg_l+Rg1/(Nml+1)/(Nml+1)
    Rg_lz=Rg_lz+Rg1z/(Nml+1)/(Nml+1)
  end do
  Re_l=Re_l/Ngl
  Re_lz=Re_lz/Ngl
  Rg_l=Rg_l/Ngl
  Rg_lz=Rg_lz/Ngl
  !!----Re_ss,Re_ssz,Rg_ss,Rg_ssz,Re_sb,Resbz,Rg_sb,Rg_sbz,Re_s,Re_sz----!
  do i=1,Nga
    m=(i-1)*(arm*Nma+1)+1
    n=(i-1)*(arm*Nma+1)+Nma+1
    call rij_and_rr(rij,rr,m,n)
    Re_ss=Re_ss+rr
    Re_ssz=Re_ssz+rij(3)*rij(3)
    Rg1=0
    Rg1z=0
    do p=1,Nma
      do q=p+1,Nma+1
        m=(i-1)*(arm*Nma+1)+p
        n=(i-1)*(arm*Nma+1)+q
        call rij_and_rr(rij,rr,m,n)
        Rg1=Rg1+rr
        Rg1z=Rg1z+rij(3)*rij(3)
      end do
    end do
    Rg_ss=Rg_ss+Rg1/(Nma+2)/(Nma+2)
    Rg_ssz=Rg_ssz+Rg1z/(Nma+2)/(Nma+2)

    do j=2,arm
      m=(i-1)*(arm*Nma+1)+1
      n=(i-1)*(arm*Nma+1)+1+j*Nma
      call rij_and_rr(rij,rr,m,n)
      Re_s=Re_s+rr
      Re_sz=Re_sz+rij(3)*rij(3)
      m=(i-1)*(arm*Nma+1)+Nma+1
      n=(i-1)*(arm*Nma+1)+1+j*Nma
      call rij_and_rr(rij,rr,m,n)
      Re_sb=Re_sb+rr
      Re_sbz=Re_sbz+rij(3)*rij(3)
      Rg2=0
      Rg2z=0
      do p=1,Nma-1
        do q=p+1,Nma
          m=(i-1)*(arm*Nma+1)+1+((j-1)*Nma+p)
          n=(i-1)*(arm*Nma+1)+1+((j-1)*Nma+q)
          call rij_and_rr(rij,rr,m,n)
          Rg2=Rg2+rr
          Rg2z=Rg2z+rij(3)*rij(3)
        end do
      end do
      Rg_sb=Rg_sb+Rg2/(Nma+1)/(Nma+1)
      Rg_sbz=Rg_sbz+Rg2z/(Nma+1)/(Nma+1)
    end do
  end do
  Re_ss=Re_ss/Nga
  Re_ssz=Re_ssz/Nga
  Rg_ss=Rg_ss/Nga
  Rg_ssz=Rg_ssz/Nga
  Re_sb=Re_sb/Nga/(arm-1)
  Re_sbz=Re_sbz/Nga/(arm-1)
  Rg_sb=Rg_sb/Nga/(arm-1)
  Rg_sbz=Rg_sbz/Nga/(arm-1)
  Re_s=Re_s/Nga/(arm-1)
  Re_sz=Re_sz/Nga/(arm-1)
!!-------------Rg_s,Rg_sz-------------!
  do i=1,Nga
    Rg2=0
    Rg2z=0
    do j=1,arm*Nma
      do k=j+1,arm*Nma+1
        m=(i-1)*(arm*Nma+1)+j
        n=(i-1)*(arm*Nma+1)+k
        call rij_and_rr(rij,rr,m,n)
        Rg2=Rg2+rr
        Rg2z=Rg2z+rij(3)*rij(3)
      end do
    end do
    Rg_s=Rg_s+Rg2/(arm*Nma+1+1)/(arm*Nma+1+1)
    Rg_sz=Rg_sz+Rg2z/(arm*Nma+1+1)/(arm*Nma+1+1)
  end do
  Rg_s=Rg_s/Nga
  Rg_sz=Rg_sz/Nga
!!----------------kenentic_energy--------------!
  kenetic_energy=0.
  do i=1,NN
    kenetic_energy=kenetic_energy+dot_product(vel(i,:),vel(i,:))/2.
  end do
  kenetic_energy=kenetic_energy/(NN-N_anchor)
end subroutine height


subroutine write_height(j)
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
  integer, intent(in) :: j

  open(35,position='append', file='./data/up_and_down.txt')
    write(35,350) 1.*j, R_up1, R_down1, R_I1, R_up2, R_down2, R_I2,  &
                        R_up3, R_down3, R_I3
  close(35)
  350 format(10F15.6)
!     open(36,position='append', file='./data/height.txt')
!       write(36,361) 1.*j, h_avg
!     close(36)
!     361 format(2F17.6)
  open(36,position='append', file='./data/height.txt')
    write(36,360) 1.*j, h_avg, hl_max, hl_end, hl_avg, hs_max, hs_end, hs_branch, hs_avg
  close(36)
  open(36,position='append', file='./data/Re.txt')
    write(36,360) 1.*j, Re_l, Re_lz, Re_ss, Re_ssz, Re_sb, Re_sbz, Re_s, Re_sz
  close(36)
  open(36,position='append', file='./data/Rg.txt')
    write(36,360) 1.*j, Rg_l, Rg_lz, Rg_ss, Rg_ssz, Rg_sb, Rg_sbz, Rg_s, Rg_sz
  close(36)
  open(36,position='append', file='./data/KE.txt')
    write(36,370) 1.*j, kenetic_energy, R_up3, R_down3, R_I3
  close(36)
  360 format(9F17.6)
  370 format(4F17.6)
end subroutine write_height


subroutine histogram
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
  use compute_acceleration
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  real*8, dimension(3) :: rij1,rij2,rij3,f1
  real*8 :: rsqr,theta,rr1,rr2,rr3,he_min,he_max,hb
  real*8 :: h_avg, Rg2_avg, Rgz2_avg, Rgxy2_avg

  call force_distribution(fene_f,lj_force_PE,lj_force_ions,coulomb_f,Bond_dist)
  
  !----------------h_dist--------------!
  n = arm*Nma + 1
  do i = 1, Nga
    h_avg = 0
    do j = 1, n
      h_avg = h_avg + pos( (i-1)*n+j, 3 )
    end do
    h_avg = h_avg/n
    k = ceiling( h_avg / (Lz/SizeHist) )
    if ( k == 0 ) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in star_h_dist: k<0 or k>SizeHist!'
      cycle
    end if
    star_h_dist(k,2) = star_h_dist(k,2)+1
  end do

  n = arm*Nma + 1
  do i = 1, Nga
    h_avg = 0
    do j = 2, arm
      h_avg = h_avg + pos( (i-1)*n+j*Nma+1, 3 )
    end do
    h_avg = h_avg / (arm-1)
    k = ceiling( h_avg / (Lz/SizeHist) )
    if ( k == 0 ) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in star_end_h_dist: k<0 or k>SizeHist!'
      cycle
    end if
    star_end_h_dist(k,2) = star_end_h_dist(k,2)+1
  end do

  n = Nga * (arm*Nma + 1)
  do i = 1, Ngl
    h_avg = 0
    do j = 1, Nml
      h_avg = h_avg + pos( n+(i-1)*Nml+j, 3 )
    end do
    h_avg = h_avg/Nml
    k = ceiling( h_avg / (Lz/SizeHist) )
    if ( k == 0 ) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in linear_h_dist: k<0 or k>SizeHist!'
      cycle
    end if
    linear_h_dist(k,2) = linear_h_dist(k,2)+1
  end do
  !------------------------------------!

  !---------------Rg_dist--------------!
  do i = 1, Nga
    call star_Rg(i,Rg2_avg,Rgz2_avg,Rgxy2_avg)
    k = ceiling( Rg2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in star_Rg2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    star_Rg2_dist(k,2) = star_Rg2_dist(k,2) + 1
    k = ceiling( Rgz2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in star_Rgz2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    star_Rgz2_dist(k,2) = star_Rgz2_dist(k,2) + 1
    k = ceiling( Rgxy2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in star_Rgxy2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    star_Rgxy2_dist(k,2) = star_Rgxy2_dist(k,2) + 1
  end do
  do i = 1,Ngl
    call linear_Rg(i,Rg2_avg,Rgz2_avg,Rgxy2_avg)
    k = ceiling( Rg2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in linear_Rg2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    linear_Rg2_dist(k,2) = linear_Rg2_dist(k,2) + 1
    k = ceiling( Rgz2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in linear_Rgz2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    linear_Rgz2_dist(k,2) = linear_Rgz2_dist(k,2) + 1
    k = ceiling( Rgxy2_avg )
    if ( k == 0 ) cycle
    if (k<0 .or. k>(floor(Lz/2)**2)) then
      write(*,*) 'Wrong in linear_Rgxy2_dist: k<0 or k>(floor(Lz/2)**2)!'
      cycle
    end if
    linear_Rgxy2_dist(k,2) = linear_Rgxy2_dist(k,2) + 1
  end do
  !------------------------------------!

  !---------------ZhangFen-------------!
  do i = 1, Npe
    if ( pos(i,4) /= 0 ) then
      k = ceiling( pos(i,3) / (Lz/SizeHist) )
      if ( k == 0 ) cycle
      if (k<0 .or. k>SizeHist) then
        write(*,*) 'Wrong in phi_tot: k<0 or k>SizeHist!'
        cycle
      end if
      phi_branch(k,2) = phi_branch(k,2) + 1 
    end if
  end do
  n = arm * Nma + 1
  do i = 1, Nga
    j = (i-1) * n + 1
    k = (i-1) * n + Nma + 1
    l = (i-1) * n + Nma*2 + 1
    m = i * n

    call rij_and_rr(rij1, rr1, j, k)
    x = ceiling( asin( abs(rij1(3)) / sqrt(rr1) ) / (pi/2/SizeHist) )
    if ( x == 0 ) cycle
    if (x<0 .or. x>SizeHist) then
      write(*,*) 'Wrong in alpha_stem: x<0 or x>SizeHist!'
      cycle
    end if
    alpha_stem(x,2) = alpha_stem(x,2) + 1 

    call rij_and_rr(rij1, rr1, l, m)
    x = ceiling( asin( abs(rij1(3)) / sqrt(rr1) ) / (pi/2/SizeHist) )
    if ( x == 0 ) cycle
    if (x<0 .or. x>SizeHist) then
      write(*,*) 'Wrong in alpha_end: x<0 or x>SizeHist!'
      cycle
    end if
    alpha_end(x,2) = alpha_end(x,2) + 1 

    call rij_and_rr(rij1, rr1, l, m)
    call rij_and_rr(rij2, rr2, k, l)
    call rij_and_rr(rij3, rr3, k, m)
    x = ceiling( acos( (rr2+rr3-rr1)/sqrt(rr2*rr3)/2 ) / (pi/SizeHist) )
    if ( x == 0 ) cycle
    if (x<0 .or. x>SizeHist) then
      write(*,*) 'Wrong in alpha_branch: x<0 or x>SizeHist!'
      cycle
    end if
    alpha_branch(x,2) = alpha_branch(x,2) + 1 

  end do
  !------------------------------------!
  
  !!!!!!!!!!!!!!!!!!!!!!!!1d_height_distribution!!!!!!!!!!!!!!!!!!!!!
  !!-------------------phi_tot-----------------!
  do i=1,Npe
    k=ceiling(pos(i,3)/(Lz/SizeHist))
    if (k==0) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in phi_tot: k<0 or k>SizeHist!'
      cycle
    end if
    phi_tot(k,2)=phi_tot(k,2)+1
  end do
  !!-----------------phi_l,phi_le---------------!
  do i=Nta+1,Npe
    !total monomers
    k=ceiling(pos(i,3)/(Lz/SizeHist))
    if (k==0) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in phi_l: k<0 or k>SizeHist!'
      cycle
    end if
    phi_l(k,2)=phi_l(k,2)+1
    !end  monomers
    if (mod(i-Nta,Nml)==0) then
      phi_le(k,2)=phi_le(k,2)+1
    end if
  end do
  !!-----------------phi_s,phi_se---------------!
  do i=1, Nta
    !total monomers
    k=ceiling(pos(i,3)/(Lz/SizeHist))
    if(k==0) cycle
    if (k<0 .or. k>SizeHist) then
      write(*,*) 'Wrong in phi_s: k<0 or k>SizeHist!'
      cycle
    end if
    phi_s(k,2)=phi_s(k,2)+1
    !branching monomers
    if (mod(i,arm*Nma+1)==Nma+1) then
      phi_sb(k,2)=phi_sb(k,2)+1
    !end monomers
    elseif (mod(mod(i,arm*Nma+1)-Nma-1,Nma)==0) then
      phi_se(k,2)=phi_se(k,2)+1
    end if
  end do
  !!--------------phi_i,phi_q,phi_q--------------!
  do i=1,Nq
    !total charged monomers
    k=ceiling(pos(charge(i),3)/(Lz/SizeHist))
    if (k<=0 .or. k>SizeHist) then
      write(*,*) 'Wrong in phi_q: k<0 or k>SizeHist!'
      cycle
    end if
    phi_q(k,2)=phi_q(k,2)+pos(charge(i),4)
    if ( i<= Nq/(nint(abs(qq))+1) ) then    !ions
      phi_i(k,2)=phi_i(k,2)+1.
    else                                    !aions
      phi_a(k,2)=phi_a(k,2)+1.
    end if
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!1d_theta_distribution!!!!!!!!!!!!!!!!!!!!!
  !!------------------theta_l,theta_lz------------------!
  do i=1,Ngl
    do j=2,Nml-1
      l=Nta+(i-1)*Nml+j
      m=Nta+(i-1)*Nml+j-1
      n=Nta+(i-1)*Nml+j+1
      call rij_and_rr(rij1,rr1,l,m)
      call rij_and_rr(rij2,rr2,l,n)
      call rij_and_rr(rij3,rr3,m,n)
      theta=acos((rr1+rr2-rr3)/2/sqrt(rr1)/sqrt(rr2))
      k=ceiling(theta/(pi/SizeHist))
      if (k<=0 .or. k>SizeHist) then
        write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
        cycle
      end if
      theta_l(k,2)=theta_l(k,2)+1
      theta=acos(rij1(3)/sqrt(rr1))
      k=ceiling(theta/(pi/SizeHist))
      if (k<=0 .or. k>SizeHist) then
        write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
        cycle
      end if
      theta_lz(k,2)=theta_lz(k,2)+1
    end do
  end do
  !!----theta_ssl,theta_sslz,theta_sbl,theta_sblz,theta_bez----!
  do i=1,Nga
    do j=2,Nma
      l=(i-1)*(arm*Nma+1)+j
      m=(i-1)*(arm*Nma+1)+j-1
      n=(i-1)*(arm*Nma+1)+j+1
      call rij_and_rr(rij1,rr1,l,m)
      call rij_and_rr(rij2,rr2,l,n)
      call rij_and_rr(rij3,rr3,m,n)
      theta=acos((rr1+rr2-rr3)/2/sqrt(rr1)/sqrt(rr2))
      k=ceiling(theta/(pi/SizeHist))
      if (k<=0 .or. k>SizeHist) then
        write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
        cycle
      end if
      theta_ssl(k,2)=theta_ssl(k,2)+1
      theta=acos(rij1(3)/sqrt(rr1))
      k=ceiling(theta/(pi/SizeHist))
      if (k<=0 .or. k>SizeHist) then
        write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
        cycle
      end if
      theta_sslz(k,2)=theta_sslz(k,2)+1
    end do
    
    do j=2,arm
      l=(i-1)*(arm*Nma+1)+j*Nma+1
      m=(i-1)*(arm*Nma+1)+Nma+1
      call rij_and_rr(rij1,rr1,l,m)
      theta=acos(rij1(3)/sqrt(rr1))
      p=ceiling(theta/(pi/SizeHist))
      theta_bez(p,2)=theta_bez(p,2)+1
      do k=2,Nma-1
        l=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k
        m=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k-1
        n=(i-1)*(arm*Nma+1)+(j-1)*Nma+1+k+1
        call rij_and_rr(rij1,rr1,l,m)
        call rij_and_rr(rij2,rr2,l,n)
        call rij_and_rr(rij3,rr3,m,n)
        theta=acos((rr1+rr2-rr3)/2/sqrt(rr1)/sqrt(rr2))
        p=ceiling(theta/(pi/SizeHist))
        if (p<=0 .or. p>SizeHist) then
          write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
          cycle
        end if
        theta=acos(rij1(3)/sqrt(rr1))
        theta_sbl(p,2)=theta_sbl(p,2)+1
        p=ceiling(theta/(pi/SizeHist))
        if (p<=0 .or. p>SizeHist) then
          write(*,*) 'Wrong in phi_se: k<0 or k>SizeHist!'
          cycle
        end if
        theta_sblz(p,2)=theta_sblz(p,2)+1
      end do
    end do
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!2d_distribution!!!!!!!!!!!!!!!!!!!!!
  !!------------------phi_zx,phi_xy,phi_szx,phi_sxy------------------!
  do i=1,Nta
    x=ceiling((pos(i,1)+Lx/2)/(Lx/SizeHist))
    y=ceiling((pos(i,2)+Ly/2)/(Ly/SizeHist))
    z=ceiling(pos(i,3)/(Lz/SizeHist))
    if (x==0 .or. y==0 .or. z==0) cycle
    phi_szx(x,z)=phi_szx(x,z)+1           !total distribution of star brushes
    phi_sxy(x,y)=phi_sxy(x,y)+1           !total distribution of star brushes
    phi_syz(y,z)=phi_syz(y,z)+1           !total distribution of star brushes
    phi_zx(x,z)=phi_zx(x,z)+1             !total distribution of all brushes
    phi_xy(x,y)=phi_xy(x,y)+1             !total distribution of all brushes
    phi_yz(y,z)=phi_yz(y,z)+1             !total distribution of all brushes
    if(mod(i,arm*Nma+1)==Nma+1) then
      phi_sbzx(x,z)=phi_sbzx(x,z)+1       !distribution of branching monomers of the star brushes
      phi_sbxy(x,y)=phi_sbxy(x,y)+1       !distribution of barnching monomers of the star brushes
      phi_sbyz(y,z)=phi_sbyz(y,z)+1       !distribution of barnching monomers of the star brushes
    elseif (mod(mod(i,arm*Nma+1)-Nma-1,Nma)==0) then
      phi_sezx(x,z)=phi_sezx(x,z)+1       !distribution of end monomers of the star brushes
      phi_sexy(x,y)=phi_sexy(x,y)+1       !distribution of end monomers of the star brushes
      phi_seyz(y,z)=phi_seyz(y,z)+1       !distribution of end monomers of the star brushes
    end if
  end do
  !!------------------phi_zx,phi_xy,phi_lzx,phi_lxy------------------!
  do i=Nta+1,Npe
    x=ceiling((pos(i,1)+Lx/2)/(Lx/SizeHist))
    y=ceiling((pos(i,2)+Ly/2)/(Ly/SizeHist))
    z=ceiling(pos(i,3)/(Lz/SizeHist))
    if (x==0 .or. y==0 .or. z==0) cycle
    phi_lzx(x,z)=phi_lzx(x,z)+1           !total distribution of linear brushes
    phi_lxy(x,y)=phi_lxy(x,y)+1           !total distribution of linear brushes   
    phi_lyz(y,z)=phi_lyz(y,z)+1           !total distribution of linear brushes   
    phi_zx(x,z)=phi_zx(x,z)+1             !total distribution of all brushes
    phi_xy(x,y)=phi_xy(x,y)+1             !total distribution of all brushes
    phi_yz(y,z)=phi_yz(y,z)+1             !total distribution of all brushes
    if (mod(i-Nta,Nml)==0) then
      phi_lezx(x,z)=phi_lezx(x,z)+1       !distribution of end monomers of the linear brushes
      phi_lexy(x,y)=phi_lexy(x,y)+1       !distribution of end monomers of the linear brushes
      phi_leyz(y,z)=phi_leyz(y,z)+1       !distribution of end monomers of the linear brushes
    end if
  end do
  !!------------------phi_qzx,phi_qxy,phi_izx,phi_ixy,phi_azx,phi_axy------------------!
  do i=1,Nq
    x=ceiling((pos(charge(i),1)+Lx/2)/(Lx/SizeHist))
    y=ceiling((pos(charge(i),2)+Ly/2)/(Ly/SizeHist))
    z=ceiling(pos(charge(i),3)/(Lz/SizeHist))
    if (x==0 .or. y==0 .or. z==0) cycle
    phi_qzx(x,z)=phi_qzx(x,z)+pos(i,4)
    phi_qxy(x,y)=phi_qxy(x,y)+pos(i,4)
    phi_qyz(y,z)=phi_qyz(y,z)+pos(i,4)
    if (i<=Nq/(nint(abs(qq))+1)) then
      phi_izx(x,z)=phi_izx(x,z)+1
      phi_ixy(x,y)=phi_ixy(x,y)+1
      phi_iyz(y,z)=phi_iyz(y,z)+1
    else
      phi_azx(x,z)=phi_azx(x,z)+1
      phi_axy(x,y)=phi_axy(x,y)+1
      phi_ayz(y,z)=phi_ayz(y,z)+1
    end if
  end do
   
end subroutine histogram

subroutine star_Rg( i, Rg2, Rgz2, Rgxy2 )
  use global_variables
  implicit none
  integer, intent(in) :: i
  real*8, intent(out) :: Rg2
  real*8, intent(out) :: Rgz2
  real*8, intent(out) :: Rgxy2
  integer :: j,k,l,m,n
  real*8 :: rij(3),rr

  l = arm*Nma+1
  do j=1,l-1
    do k=j+1,l
      m=(i-1)*l+j
      n=(i-1)*l+k
      call rij_and_rr(rij,rr,m,n)
      Rg2=Rg2+rr
      Rgz2=Rgz2+rij(3)*rij(3)
      Rgxy2=Rgxy2+rij(1)*rij(1)+rij(2)*rij(2)
    end do
  end do
  Rg2 = Rg2 / (l+1) / (l+1)
  Rgz2 = Rgz2 / (l+1) / (l+1)
  Rgxy2 = Rgxy2 / (l+1) / (l+1)
end subroutine star_Rg

subroutine linear_Rg( i, Rg2, Rgz2, Rgxy2 )
  use global_variables
  implicit none
  integer, intent(in) :: i
  real*8, intent(out) :: Rg2
  real*8, intent(out) :: Rgz2
  real*8, intent(out) :: Rgxy2
  integer :: j,k,m,n
  real*8 :: rij(3),rr

  do j=1,Nml-1
    do k=j+1,Nml
      m=Nta+(i-1)*Nml+j
      n=Nta+(i-1)*Nml+k
      call rij_and_rr(rij,rr,m,n)
      Rg2=Rg2+rr
      Rgz2=Rgz2+rij(3)*rij(3)
      Rgxy2=Rgxy2+rij(1)*rij(1)+rij(2)*rij(2)
    end do
  end do
  Rg2 = Rg2 / (Nml+1) / (Nml+1)
  Rgz2 = Rgz2 / (Nml+1) / (Nml+1)
  Rgxy2 = Rgxy2 / (Nml+1) / (Nml+1)
end subroutine linear_Rg

subroutine write_hist
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
  !------------------------------------!
  use global_variables
  implicit none
  integer i,j

  !----------------h_hist--------------!
  open(28,file='./data/h_dist.txt')
    do i = 1, SizeHist
      write(28,280) i*Lz/SizeHist, linear_h_dist(i,2), star_h_dist(i,2), &
                    star_end_h_dist(i,2)
    end do
  close(28)
  280 format(4F17.8)
  !------------------------------------!

  !---------------Rg_hist--------------!
  open(29,file='./data/Rg_dist.txt')
    do i = 1, floor(Lz/2)**2
      write(29,290) i*1., linear_Rg2_dist(i,2), linear_Rgz2_dist(i,2), &
                    linear_Rgxy2_dist(i,2), star_Rg2_dist(i,2), &
                    star_Rgz2_dist(i,2), star_Rgxy2_dist(i,2)
    end do
  close(29)
  290 format(7F17.8)
  !------------------------------------!

  !--------------ZhangFen--------------!
  open(30,file='./data/alpha_phi.txt')
    do i = 1, SizeHist
      phi_branch(i,1) = i*Lz/SizeHist
      alpha_stem(i,1) = i*90./SizeHist
      alpha_branch(i,1) = i*180./SizeHist
      alpha_end(i,1) = i*90./SizeHist
      write(30,300) phi_branch(i,1), phi_branch(i,2),   alpha_stem(i,1),  &
                    alpha_stem(i,2), alpha_branch(i,1), alpha_branch(i,2), &
                    alpha_end(i,1), alpha_end(i,2) 
    end do
  close(30)
  300 format(8F17.6)
  !------------------------------------!

  open(31,file='./data/phi.txt')
    do i=1,SizeHist
      phi_tot(i,1)=i*Lz/SizeHist
      write(31,340) phi_tot(i,1),phi_tot(i,2),phi_l(i,2),phi_le(i,2),  &
                    phi_s(i,2),phi_sb(i,2),phi_se(i,2), phi_a(i,2),    &
                    phi_i(i,2), phi_q(i,2)      
    end do
    340 format(10F17.6)
  close(31)
  open(32,file='./data/theta.txt')
    do i=1,SizeHist
      theta_l(i,1)=i*pi/SizeHist
      write(32,350) theta_l(i,1),    theta_l(i,2),    theta_lz(i,2),    &
                    theta_ssl(i,2),  theta_sslz(i,2), theta_sbl(i,2),   &
                    theta_sblz(i,2), theta_bez(i,2)
    end do
    350 format(8F17.6)
  close(32)
  
  open(33,file='./data/force_liear.txt')
    do i=1,Nml-1
      write(33,330) i*1., force_l(i,2)
    end do
    330 format(2F17.6)
  close(33)

  open(34,file='./data/force_star.txt')
    do i=1,Nma*2
      write(34,320) i*1., force_sy(i,2), force_sn(i,2), force_so(i,2)
    end do
    320 format(4F17.6)
  close(34)

    open(331,file='./data/force_liear1.txt')
    do i=1,SizeHist
      write(331,3301) i*1., force_l1(i,2)
    end do
    3301 format(2F17.6)
  close(331)

  open(341,file='./data/force_star1.txt')
    do i=1,SizeHist
      write(341,3201) i*1., force_sy1(i,2), force_sn1(i,2), force_so1(i,2)
    end do
    3201 format(4F17.6)
  close(341)
  
  open(35,file='./data/delta_angle.txt')
    do i=1,SizeHist
      write(35,310) i*pi/SizeHist, delta_angle1(i,2), delta_angle2(i,2), delta_angle3(i,2)
    end do
    310 format(4F17.6)
  close(35)
  
  open(51,file='./data/phi_2d11.txt')
  open(52,file='./data/phi_2d12.txt')
  open(53,file='./data/phi_2d13.txt')
  open(54,file='./data/phi_2d21.txt')
  open(55,file='./data/phi_2d22.txt')
  open(56,file='./data/phi_2d23.txt')
  open(57,file='./data/phi_2d31.txt')
  open(58,file='./data/phi_2d32.txt')
  open(59,file='./data/phi_2d33.txt')
  open(60,file='./data/phi_2d41.txt')
  open(61,file='./data/phi_2d42.txt')
  open(62,file='./data/phi_2d43.txt')
  open(63,file='./data/phi_2d51.txt')
  open(64,file='./data/phi_2d52.txt')
  open(65,file='./data/phi_2d53.txt')
  open(66,file='./data/phi_2d61.txt')
  open(67,file='./data/phi_2d62.txt')
  open(68,file='./data/phi_2d63.txt')
  open(69,file='./data/phi_2d71.txt')
  open(70,file='./data/phi_2d72.txt')
  open(71,file='./data/phi_2d73.txt')
  open(72,file='./data/phi_2d81.txt')
  open(73,file='./data/phi_2d82.txt')
  open(74,file='./data/phi_2d83.txt')
  open(75,file='./data/phi_2d91.txt')
  open(76,file='./data/phi_2d92.txt')
  open(77,file='./data/phi_2d93.txt')
  open(78,file='./data/fene_f.txt')
  open(79,file='./data/lj_force_PE.txt')
  open(80,file='./data/lj_force_ions.txt')
  open(81,file='./data/coulomb_f.txt')
  open(82,file='./data/Bond_dist.txt')
    do i=1,SizeHist
      write(51,'(500I10)') (phi_zx(i,j),j=1,SizeHist) 
      write(52,'(500I10)') (phi_xy(i,j),j=1,SizeHist)  
      write(53,'(500I10)') (phi_yz(i,j),j=1,SizeHist)  
      write(54,'(500I10)') (phi_lzx(i,j),j=1,SizeHist)  
      write(55,'(500I10)') (phi_lxy(i,j),j=1,SizeHist) 
      write(56,'(500I10)') (phi_lyz(i,j),j=1,SizeHist) 
      write(57,'(500I10)') (phi_lezx(i,j),j=1,SizeHist)  
      write(58,'(500I10)') (phi_lexy(i,j),j=1,SizeHist)  
      write(59,'(500I10)') (phi_leyz(i,j),j=1,SizeHist)  
      write(60,'(500I10)') (phi_szx(i,j),j=1,SizeHist) 
      write(61,'(500I10)') (phi_sxy(i,j),j=1,SizeHist)  
      write(62,'(500I10)') (phi_syz(i,j),j=1,SizeHist)  
      write(63,'(500I10)') (phi_sbzx(i,j),j=1,SizeHist)  
      write(64,'(500I10)') (phi_sbxy(i,j),j=1,SizeHist) 
      write(65,'(500I10)') (phi_sbyz(i,j),j=1,SizeHist) 
      write(66,'(500I10)') (phi_sezx(i,j),j=1,SizeHist)  
      write(67,'(500I10)') (phi_sexy(i,j),j=1,SizeHist)  
      write(68,'(500I10)') (phi_seyz(i,j),j=1,SizeHist)  
      write(69,'(500I10)') (phi_azx(i,j),j=1,SizeHist) 
      write(70,'(500I10)') (phi_axy(i,j),j=1,SizeHist)  
      write(71,'(500I10)') (phi_ayz(i,j),j=1,SizeHist)  
      write(72,'(500I10)') (phi_izx(i,j),j=1,SizeHist)  
      write(73,'(500I10)') (phi_ixy(i,j),j=1,SizeHist) 
      write(74,'(500I10)') (phi_iyz(i,j),j=1,SizeHist) 
      write(75,'(500I10)') (phi_qzx(i,j),j=1,SizeHist)  
      write(76,'(500I10)') (phi_qxy(i,j),j=1,SizeHist)   
      write(77,'(500I10)') (phi_qyz(i,j),j=1,SizeHist)   
      write(78,'(500I10)') (fene_f(i,j),j=1,SizeHist)   
      write(79,'(500I10)') (lj_force_PE(i,j),j=1,SizeHist)   
      write(80,'(500I10)') (lj_force_ions(i,j),j=1,SizeHist)   
      write(81,'(500I10)') (coulomb_f(i,j),j=1,SizeHist)   
      write(82,'(500I10)') (Bond_dist(i,j),j=1,SizeHist)   
    end do
  close(51)
  close(52)
  close(53)
  close(54)
  close(55)
  close(56)
  close(57)
  close(58)
  close(59)
  close(60)
  close(61)
  close(62)
  close(63)
  close(64)
  close(65)
  close(66)
  close(67)
  close(68)
  close(69)
  close(70)
  close(71)
  close(72)
  close(73)
  close(74)
  close(75)
  close(76)
  close(77)
  close(78)
  close(79)
  close(80)
  close(81)
  close(82)
end subroutine write_hist


subroutine write_pos
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
  integer :: i

  open(30,file='./data/pos.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3), pos(i,4)
      300 format(4F17.6)
    end do
  close(30)
end subroutine write_pos


subroutine write_pos1
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
  integer :: i

  open(31,file='./data/pos1.txt')
    do i=1, NN
      write(31,310) pos(i,1), pos(i,2), pos(i,3), pos(i,4)
      310 format(4F17.6)
    end do
  close(31)
end subroutine write_pos1 


subroutine write_vel
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
  integer :: i

  open(32,file='./data/vel.txt')
    do i=1, NN
      write(32,320) vel(i,1), vel(i,2), vel(i,3)
      320 format(3F17.6)
    end do
  close(32)

end subroutine write_vel


subroutine write_vel1(j)
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
  integer :: i
  integer, intent(in) :: j
  
  open(33,file='./data/vel1.txt')
    do i=1, NN
      write(33,330) vel(i,1), vel(i,2), vel(i,3)
      330 format(3F17.6)
    end do
  close(33)

  open(32,file='./start_time.txt')
    write(32,*) 1
    write(32,*) j
    call cpu_time(finished)
    total_time=total_time+finished-started
    call cpu_time(started)
    write(32,*) total_time
    write(32,*) 'time:(minutes)', real(total_time/60)
    write(32,*) 'time:(hours)', real(total_time/3600)
    write(32,*) 'time:(days)', real(total_time/86400)
  close(32)

end subroutine write_vel1


subroutine write_acc
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
  integer :: i
  
  open(37,file='./data/acc.txt')
    do i=1, NN
      write(37,370) acc(i,1), acc(i,2), acc(i,3)
      370 format(3F17.6)
    end do
  close(37)

end subroutine write_acc


subroutine write_time(time)
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
  real*8, intent(in) :: time

  open(10,file='./data/time.txt')
    write(10,*) 'time:(seconds)', real(total_time)
    write(10,*) 'time:(minutes)', real(total_time/60)
    write(10,*) 'time:(hours)  ', real(total_time/3600)
    write(10,*) 'time:(days)   ', real(total_time/86400)
    write(10,*) 'Lx:           ', real(Lx)
    write(10,*) 'Ly:           ', real(Ly)
    write(10,*) 'Lz:           ', real(Lz)
    write(10,*) 'Nq:           ', Nq
    write(10,*) 'NN:           ', NN
    write(10,*) 'sigmag:       ', real(sigmag)
    write(10,*) 'qq:           ',nint(qq)
  close(10)

end subroutine write_time


end module input_output
