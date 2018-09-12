module initialize_update
implicit none
contains


subroutine Initialize_position
  !--------------------------------------!
  !Initialize position array
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  uniform_or_random, Nq, qq, charge
  !Routine Referenced:
  !1.random_star_brushes(l)
  !2.random_linear_brushes(l)
  !3.uniform_star_brushes(l)
  !4.uniform_linear_brushes(l)
  !5.initialize_ions
  !note:
  !1.Because initialize can be unreal, there are thousand initial
  !  method in thousand people's mind.
  !2.The monomer distributed in the space by random walk , so it may go to 
  !  cul-de-sac when the density is high.
  !--------------------------------------!
  use compute_acceleration, only : charge, lj_verlet_list, real_verlet_list, &
                                   real_verlet
  use global_variables
  implicit none
  integer :: l, i

  l=0
  !
  !Graft star and linear brushes
  if ( uniform_or_random == 0 ) then
    call uniform_linear_brushes(l)
    call uniform_star_brushes(l)
  else
    call random_star_brushes(l)
    call random_linear_brushes(l)
  end if
  call period_condition_pos
  !
  !random distribution in the box
  if ( qq /= 0 ) then
    call initialize_ions
  end if
  !
  !initialize chargen on PE
  do i=1,Nq/(nint(abs(qq))+1)
    pos(charge(i),4) = qq                 
  end do
  !
  !initialize lj_verlet_list
  call lj_verlet_list
  !
  !Construct the real verlet list and real_point vector
  if ( real_verlet == 1 ) then
    call real_verlet_list
  end if

  write(*,*) 'Inititalize position is finished!'

end subroutine Initialize_Position


subroutine random_star_brushes(l)
  !--------------------------------------!
  !Initialize star brushes.
  !
  !Input
  !  l, pos
  !Output
  !  l, pos
  !External Variables
  !  
  !Routine Referenced:
  !   rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, m, n, x, y, p
  integer, intent(inout) :: l
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij

  !star brushes
  write(*,*) 'initilize random star brushes'
  do i=1, Nga
    m=1
    l=l+1
    do while (m==1)     !anchor point
      m=0
      call random_number(rnd1)
      call random_number(rnd2)
      pos(l,1)=rnd1*Lx-Lx/2
      pos(l,2)=rnd2*Ly-Ly/2
      pos(l,3)=0
      !
      !Judge whether the particle is close the former paritcle
      !too much.
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
          if (k==1 .and. j/=1) then   !branching point
            !
            !uniform distribution on the sphere surface
            pos(l,1)=pos(l-(j-2)*Nma-1,1)+R_bond*cos(2*pi*rnd2)*sin(pi*rnd1)
            pos(l,2)=pos(l-(j-2)*Nma-1,2)+R_bond*sin(2*pi*rnd2)*sin(pi*rnd1)
            pos(l,3)=pos(l-(j-2)*Nma-1,3)+R_bond*cos(pi*rnd1)
          else
            if(p<10) then
              rnd1=rnd1/1000 !let the monomer go to the up direction
            end if
            pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
            pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
            pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
          end if
          !periodic condition
          call period_condition_rij(pos(l,1:3))
          !
          !Jugde whether the particle is close the former paritcle
          !too much.
          do n=1,l-1
            call rij_and_rr(rij,rsqr,n,l)
            if (rsqr<0.7 .or. pos(l,3)<0.9) then
              m=1
              p=p+1
              cycle
            end if
          end do
        end do
      end do
    end do
  end do
end subroutine random_star_brushes


subroutine random_linear_brushes(l)
  !--------------------------------------!
  !Initialize linear brushes.
  !
  !Input
  !  l, pos
  !Output
  !  l, pos
  !External Variables
  !  
  !Routine Referenced:
  !  rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, m, n, x, y, p
  integer, intent(inout) :: l
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij

  !linear brushes
  write(*,*) 'initialize random linear brushes'
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
      !
      !Jugde whether the particle is close the former paritcle
      !too much.
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
!           rnd1=rnd1**2
        if (p<10) then
          rnd1=rnd1/2
        else
          rnd1=rnd1**2
        end if
        call random_number(rnd2)
        pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
        pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
        pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
        !periodic condition
        call period_condition_rij(pos(l,1:3))
        !
        !Jugde whether the particle is close the former paritcle
        !too much.
        do n=1,l-1
          call rij_and_rr(rij,rsqr,n,l)
          if (rsqr<0.7 .or. pos(l,3)<0.9) then
            m=1
            p=p+1
            cycle
          end if
        end do
      end do
    end do
  end do

end subroutine random_linear_brushes


subroutine uniform_star_brushes(l)
  !--------------------------------------!
  !Initialize star brushes.
  !
  !Input
  !  l, pos
  !Output
  !  l, pos
  !External Variables
  !  
  !Routine Referenced:
  !   rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, m, n, x, y, p
  integer, intent(inout) :: l
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij
  l=0
  !
  !star brushes
  write(*,*) 'initialize uniform star brushes'
  do i=1, N_anchor
    if ( Nga == 0 ) cycle
    x=(i-1)/nint(sqrt(N_anchor*1.))+1    !the root of N_anchor must be integer
    y=mod(i-1,nint(sqrt(N_anchor*1.)))+1
    if (Nga == Ngl .and. mod(x+y,2)==0) cycle
    m=1
    l=l+1
    do while (m==1)     !anchor point
      m=0
      pos(l,1)=Lx/nint(sqrt(N_anchor*1.))*(x-0.5)-Lx/2
      pos(l,2)=Ly/nint(sqrt(N_anchor*1.))*(y-0.5)-Ly/2
      pos(l,3)=0
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
          if (k==1 .and. j/=1) then   !branching point
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
          call period_condition_rij(pos(l,1:3))
          !
          !Judge whether the particle is close the former paritcle
          !too much.
          do n=1,l-1
            call rij_and_rr(rij,rsqr,n,l)
            if (rsqr<0.7 .or. pos(l,3)<0.9) then
              m=1
              p=p+1
              cycle
            end if
          end do
        end do
      end do
    end do
  end do
  
end subroutine uniform_star_brushes


subroutine uniform_linear_brushes(l)
  !--------------------------------------!
  !Initialize star brushes.
  !
  !Input
  !  l, pos
  !Output
  !  l, pos
  !External Variables
  !  
  !Routine Referenced:
  !   rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, m, n, x, y, p
  integer, intent(inout) :: l
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij

  !linear brushes
  write(*,*) 'uniform linear brushes'
  do i=1, N_anchor
    if ( Ngl == 0 ) cycle
    x=(i-1)/nint(sqrt(N_anchor*1.))+1
    y=mod(i-1,nint(sqrt(N_anchor*1.)))+1
    if (Nga == Ngl .and. mod(x+y,2)==1) cycle
    m=1
    l=l+1
    do while (m==1)
      m=0
      pos(l,1)=Lx/nint(sqrt(N_anchor*1.))*(x-0.5)-Lx/2
      pos(l,2)=Ly/nint(sqrt(N_anchor*1.))*(y-0.5)-Ly/2
      pos(l,3)=0
    end do
    do k=2, Nml
      l=l+1
      m=1
      p=0
      do while (m==1)
        m=0
        call random_number(rnd1)
!           rnd1=rnd1**2
        if (p<10) then
          rnd1=rnd1/5
        else
          rnd1=rnd1**2
        end if
        call random_number(rnd2)
        pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
        pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
        pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
        !periodic condition
        call period_condition_rij(pos(l,1:3))
        !
        !Judge whether the particle is close the former paritcle
        !too much.
        do n=1,l-1
          call rij_and_rr(rij,rsqr,n,l)
          if (rsqr<0.7 .or. pos(l,3)<0.9) then
            m=1
            p=p+1
            cycle
          end if
        end do
      end do
    end do
  end do

end subroutine uniform_linear_brushes


subroutine initialize_ions
  !--------------------------------------!
  !initialize ions
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  
  !Routine Referenced:
  !  rij_and_rr
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, x, y, p
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij

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
    if ( i <= ( NN - Nq_salt_ions * nint( abs(qqi)+1 ) ) ) then
      pos(i,4) = - qq / abs(qq)
    elseif ( i <= ( NN - Nq_salt_ions ) ) then
      pos(i,4) = - qqi / abs(qqi)
    else
      pos(i,4) = qqi
    end if
  end do

end subroutine initialize_ions



subroutine initialize_velocity
  !-----------------------------------------!
  !Initialize velocity by Gauss distribution.
  !
  !input:
  !  vel
  !output:
  !  vel
  !External variants:
  !  
  !Routine Referenced:
  !  rescale_velocity
  !-----------------------------------------!
  use compute_acceleration, only : anchor_list
  use global_variables
  implicit none
  integer i,j,k
  real*8 rnd
  !
  !generate gauss velocity distribution
  do i = 1, NN
    do j = 1, 3
      call gauss_dist(0.D0,sqrt(1./Beta),rnd)
      vel(i,j) = rnd
    end do
  end do
  !
  !the anchored monomer is constrained on the plate.
  do i = 1, N_anchor
    vel(anchor_list(i),:) = 0
  end do
  call rescale_velocity
end subroutine initialize_velocity


subroutine rescale_velocity
  !--------------------------------------!
  !In Langevin dynamics, we needn't to rescale velocity to
  !achieve equal temperatrue condition.
  !Here, we rescale velocity to avoid the break of chemical bonds
  !
  !Input
  !  
  !Output
  !   
  !External Variables
  !  
  !Routine Referenced:
  !  Frenkel, Smit, 'Understanding molecular simulation: from
  !  algorithm to applications', Elsevier, 2002, pp.66
  !  (Algorithm 4).
  !--------------------------------------!
  use compute_acceleration, only : anchor_list
  use global_variables
  implicit none
  real*8 :: v2_sum     
  integer :: i, j
  real*8, dimension(3) :: v_sum
  
  v_sum  = 0.
  v2_sum = 0.
  !
  !average veloctiy
  do i = 1, NN
    do j = 1, 3
      v_sum(j) = v_sum(j) + vel(i,j)
    end do
  end do
  v_sum = v_sum / (NN-N_anchor);
  !
  !average magnitude of velocity
  do i = 1,NN
    do j = 1,3
      vel(i,j) = vel(i,j) - v_sum(j)
    end do
    v2_sum = v2_sum + vel(i,1)*vel(i,1) + vel(i,2)*vel(i,2) + vel(i,3)*vel(i,3)
  end do
  !
  !the anchored monomer is constrained on the plate.
  do i = 1, N_anchor
    v2_sum = v2_sum - dot_product(vel(anchor_list(i),:),vel(anchor_list(i),:))
    vel(anchor_list(i),:) = 0
  end do
  !
  !compute mean square velocity
  v2_sum = v2_sum / (NN-N_anchor-1)
  !
  !rescale_velocity
  vel = vel * sqrt(3./Beta/v2_sum)
end subroutine rescale_velocity


subroutine new_position
  !--------------------------------------!
  !Velocity Verlet algorithm of Langevin dynamics
  !
  !
  !Input
  !  pos, vel, acc
  !Output
  !  pos, vel, acc
  !External Variables
  !  dr_max1, dr_max2, xi, NN, dt
  !Routine Referenced:
  !1.Tamar Schlick, 'Molecular Modeling and Simulation:
  !  An Interdisciplinary Guide', 2nd edition, pp.482(14.30),
  !  Springer, 2010.
  !2.Gromacs, 'manul 5.0.7', pp.50(3.115-3.119), 2015.
  !--------------------------------------!
  use global_variables
  use compute_acceleration, only : dr_max1, dr_max2, xi, &
                                   compute_force
  implicit none
  real*8, dimension(3) :: dri
  real*8 dr1,dr2,dmax1,dmax2
  integer i
  !
  !move
  vel=acc*dt/2+(1-xi*dt/2)*vel
  !
  !calculate max displacement
  dmax1=0.
  dmax2=0.
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
  !
  !compute force and move
  call period_condition_pos
  call compute_force
  vel=(vel+acc*dt/2.)/(1+xi*dt/2.)
end subroutine new_position


subroutine period_condition_pos
  !--------------------------------------!
  !Period_condition of the postion array.
  !
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  Lx, Ly, NN
  !Routine Referenced:
  !
  !note:
  !  the range of the system is
  !  (-Lx/2, Lx/2)*(-Ly/2, Ly/2)
  !--------------------------------------!
  use global_variables
  implicit none
  integer k
  
  do k=1, NN
    if (pos(k,1)>Lx/2) then
      pos(k,1)=pos(k,1)-Lx
    elseif(pos(k,1)<=-Lx/2) then
      pos(k,1)=pos(k,1)+Lx
    end if
    if (pos(k,2)>Ly/2) then
      pos(k,2)=pos(k,2)-Ly
    elseif(pos(k,2)<=-Ly/2) then
      pos(k,2)=pos(k,2)+Ly
    end if
  end do
end subroutine period_condition_pos


subroutine period_condition_rij(rij)
  !--------------------------------------!
  !Period_condition of the postion array.
  !
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  Lx, Ly, NN
  !Routine Referenced:
  !
  !note:
  !  the range of the system is
  !  (-Lx/2, Lx/2)*(-Ly/2, Ly/2)
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, dimension(3), intent(inout) :: rij
  
  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if
end subroutine period_condition_rij


! subroutine avoid_particle_too_close(l, judge)
!   !--------------------------------------!
!   !Avoid this particle is so close to the all former particles
!   !that the force between them are very large.
!   !
!   !
!   !Input
!   !  l: 
!   !Output
!   !  judge: 0 or 1
!   !External Variables
!   !  pos
!   !Routine Referenced:
!   !  rij_and_rr
!   !--------------------------------------!
!   use global_variables
!   implicit none
!   integer, intent(in)  :: l
!   integer, intent(out) :: judge
!   integer :: n
!   real*8 :: rij(3), rsqr

!   do n=1,l-1
!     call rij_and_rr(rij,rsqr,n,l)
!     if (rsqr<0.7 .or. pos(l,3)<0.9) then
!       judge = 1
!       cycle
!     end if
!   end do

! end subroutine

! subroutine next_position(post, l, rnd1, rnd2)
!   !--------------------------------------!
!   !Generate the post, which is random distributed
!   !on the sphere of particle l with radius R_bond
!   !
!   !
!   !Input
!   !  l, rnd1, rnd2 
!   !Output
!   !  pos(l,1:3)
!   !External Variables
!   !  pos
!   !Routine Referenced:
!   !  rij_and_rr
!   !--------------------------------------!
!   use global_variables
!   implicit none
!   integer, intent(in)  :: l
!   real*8,  intent(in)  :: rnd1
!   real*8,  intent(in)  :: rnd2
!   real*8,  intent(out) :: post(3)

!   post(1)=pos(l,1) + R_bond * sin(pi*rnd1) * cos(2*pi*rnd2)
!   post(2)=pos(l,2) + R_bond * sin(pi*rnd1) * sin(2*pi*rnd2)
!   post(3)=pos(l,3) + R_bond * cos(pi*rnd1)

! end subroutine next_position

end module initialize_update


