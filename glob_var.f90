module global_variables
implicit none
save
!########################constants#########################!
  real*8, parameter:: pi=3.141592653589793D0      
  real*8, parameter:: gamma=.5772156649015329D0   !Euler Gamma
!########################constants#########################!

!####################systems coefficient###################!
  integer :: Npe      !Total monomers in Polyelectrolytes(PE)
  integer :: arm      !Arms of star brushes
                      !including the chains anchored to the plate
  integer :: Nma      !Monomers of each arm
  integer :: Nml      !Monomers of each linear chain
  integer :: Nga      !Number of star chains grafted on plate
  integer :: Ngl      !Number of linear chains grafted on plate
  integer :: Nta      !Total monomers of star brushes
  integer :: Ntl      !Total monomers of linear brushes
  integer :: Nq       !Total charge in the system
  integer :: NN       !Total particles in the system
  integer :: N_anchor !Anchored chains
  integer :: man_l    !Manning effect: every man pariticle have one charge
  integer :: man_s    !Manning effect: star chains
  integer :: Nq_salt_ions !Charged salt ions, which not include anions.
  real*8  :: Lx       !Length of cell in x direction
  real*8  :: Ly       !Length of cell in y direction
  real*8  :: Lz       !Distance of two plate
  real*8  :: Z_empty  !Empty space ratio of height and length in slab geometry
  real*8  :: ratio_xy !Rotio of length x and width y of the box
  real*8  :: sigmag   !Grafting density of brushes on the plate
  real*8  :: Beta     !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
  real*8  :: qq       !Charge of charged monomers
  real*8  :: qqi      !Charge of salt ions
  real*8  :: ion_ratio!Ratio of salt ions to the charge quantites of PE
  integer :: Nq_PE    !Charged monomers of PE
  real*8  :: R_bond   !length of chemical band

!##################end systems coefficient#################!

!##################running and Histogram###################!
  integer :: restart_or_continue  !restart or continue after breaking off
  integer :: uniform_or_random    !uniform grafted or random grafted
  integer :: StepNum0             !Steps of preheating
  integer :: StepNum              !Steps of running
  integer :: DeltaStep1           !step inteval, height
  integer :: DeltaStep2           !step inteval, histogram
  integer :: DeltaStep3           !step inteval, write data
  integer :: step                 !ordering number
  integer :: multistep            !each multistep recalculate the coulomb force
  real*8  :: dt                   !time of each move
  !
  !timing
  real*8  :: started              !time at starting
  real*8  :: finished             !time at finishing
  real*8  :: total_time=0         !total time of the simulation
  !
  !histogram
  integer :: SizeHist=500        !number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
  real*8, allocatable, dimension(:,:) :: pos    !array of position
  real*8, allocatable, dimension(:,:) :: vel    !array of velocity
  real*8, allocatable, dimension(:,:) :: acc    !array of accelaration
  integer, allocatable, dimension(:) :: charged_arm
!########################end arrays########################!


contains

  subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  post(pos or pos1), i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lz(used in period condition)
  !note:
  !  including period condition
  !-----------------------------------------!
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = pos(i,1:3) - pos(j,1:3)

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

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr


subroutine gauss_dist(mu, sigma, rnd)
  !---------------------------------------!
  !
  !---------------------------------------!
  implicit none
  real*8, intent(out) :: rnd
  real*8, intent(in) :: mu, sigma
  real*8 rnd1, rnd2
  
  call random_number(rnd1)
  call random_number(rnd2)
  rnd=sqrt(-2*log(rnd1))*cos(2*pi*rnd2)
  rnd=mu+rnd*sigma
end subroutine gauss_dist


end module global_variables
