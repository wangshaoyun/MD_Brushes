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
	integer :: uniform_or_random		!Uniform grafted or random grafted
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

end module global_variables
