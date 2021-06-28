module FatesInsectMemMod
        use FatesConstantsMod, only : r8 => fates_r8
	!use EDParamsMod             , only : insect_an

	implicit none
	integer,parameter insectType= 2  	! 1= MPB, 2=WPB 3=OtherType.

	! Here are input parameter for the MPBAttack subroutine
	! Need to move these parameters into the parameter file
	!real(r8),parameter :: an = insect_an                        

	integer,  parameter :: DomainSize = 2**8                    ! domain size 256 categories for physiological ages
	integer,  parameter :: numberInsectTypes = 1                ! number of insect types (currently only one-mountain pine beetle)
	integer,  parameter :: maxNumStages   = 20                  ! maximum number of stages for insect development
	integer,  parameter :: maxpft = 14					        ! maximum number of plant pfts (I've manually entered the current number, which is inelegant...)
						
						
    ! Here are input parameter for each species subroutine
    	select case (insectType)  ! 1= MPB, 2=WPB 3=OtherType. 
		case(1)
			real(r8), parameter :: an = -11.462723_r8			! controls tree loss rate as a function of beetle and tree density
			real(r8), parameter :: ab = -5.599606_r8			! controls proportion of beetles that attack
			real(r8), parameter :: alpha3 = -35.4185_r8		    ! The temperature in degrees centigrade at which only 50 % survival occurs
			real(r8), parameter :: Beta3 = 2.0_r8				! controls the rate of change of survival probability as a function of temperature 
			real(r8), parameter :: EndPopn= 40.0_r8             ! Generic Endemic population variable. Estimated at 40 for MPB
			real(r8), parameter :: FecMax = 54.67_r8			! Maximum fecundity per parent.
			real(r8), parameter :: FecMortR = .2526481_r8		! The rate of mortality (non-overwinter mortality) Egg to Adult.
			real(r8), parameter :: Mort_Fec = -18.0_r8		    ! The temperature at which ovipositing adults die
			real(r8), parameter :: Mort_EPT = -18.0_r8          ! The temperature at which Parents, Pupae and Eggs die
			real(r8), parameter :: Mort_Ads, = -18.0_r8         ! The temperature at which flying adults die
			real(r8), parameter :: FFTL = 17.53_r8				! The minimum temperature at which flight can occur
			real(r8), parameter :: FFTH = 42.00_r8				! The maximum temperature at which flight can occur
			real(r8), parameter :: FF1 = 25.00_r8 				! Flight is based on an equation with shape p(flight)=b0+(b1*x)+(b2*x**2)....+(b5*x**5) This is b0
			real(r8), parameter :: FF2 = -5.324e+00_r8          ! From above b1
			real(r8), parameter :: FF3 = 4.277e-01_r8	        ! From above b2
			real(r8), parameter :: FF4 = -1.633e-02_r8		    ! From above b3
			real(r8), parameter :: FF5 = 3.014e-04_r8           ! From above b4
			real(r8), parameter :: FF6 = -2.172e-06_r8          ! From above b5
			
		case(2)
			real(r8), parameter :: r1 = 0.2009 			        ! controls tree loss rate as a function of beetle and tree density
			real(r8), parameter :: x0 = 10.03 				    ! controls the intercept of the TDIA attack equation 
			real(r8), parameter :: x1 = 1.545			        ! The influence of drought on the number of beetles needed to kill a tree
			real(r8), parameter :: CWDvec = 2.0_r8				! This is currently a dummy variable for drought---should be replaced
			real(r8), parameter :: x2 = 0.506_r8                ! The influence of host class size (0-<31.6 cm and 1 for >31.6 cm)
			real(r8), parameter :: SizeFactor = 0.051_r8        ! The influence on attack success on choosing larger host
			real(r8), parameter :: EndPopn= 840.0_r8            ! Generic Endemic population variable. Estimated at 840 for WPB
			real(r8), parameter :: FecMax = 24_r8               ! Maximum fecundity per parent.
			real(r8), parameter :: FecMortR = .142857_r8        ! The rate of mortality (non-overwinter mortality) Egg to Adult.
			real(r8), parameter :: Mort_Fec = -12.2_r8          ! The temperature at which ovipositing adults die
			real(r8), parameter :: Mort_EPT = -15.0_r8          ! The temperature at which Parents, Pupae and Eggs die
			real(r8), parameter :: Mort_Ads= -12.2_r8           ! The temperature at which flying adults die
			real(r8), parameter :: FFTL = 18.6_r8               ! The minimum temperature at which flight can occur
			real(r8), parameter :: FFTH = 38.9_r8               ! The maximum temperature at which flight can occur
			real(r8), parameter :: FF1 = -1.29873473_r8         ! Flight is based on an equation with shape p(flight)=b0+(b1*x)+(b2*x**2)....+(b5*x**5) This is b0
			real(r8), parameter :: FF2 = 0.1090086_r8           ! From above b1
			real(r8), parameter :: FF3 = -0.0019475_r8	   	    ! From above b2
			real(r8), parameter :: FF4 = 0.0_r8                 ! From above b3
			real(r8), parameter :: FF5 =0.0_r8 					! From above b4
			real(r8), parameter :: FF6 =0.0_r8                  ! From above b5 
		case default
			fates_endrun("Missing InsectType Parameter") ! if there s not InsectType stop run 
		end select 
!-------------------------------------------------------------------------------------------------------------------------------
    ! Defining a site-level type from which to obtain variables.
    type ed_site_insect_type
        !variables for tracking insect dynamics
	
	! array to define the preference of insects (0-not preferred; 1-preferred)
	integer, allocatable :: InsectPFTPref(:,:)	

     ! Containers for the distributions of physiological age for each life stage. 
	! The PhysAge array holds physiological age distributions for each of the relevant stages of the beetles.
	real(r8), allocatable :: PhysAge(:,:)		               ! array to hold physiological age for eggs, L1, L2, L3, L4, P, T
	
	! Containers for densities of insects transitioning from one stage to another for each insect type.
	! The Transit array holds beetles trasitionaing into the egg life stage, into L1, into L2, into L3, 
	! into L4, into the pupal stage, into the teneral adult life stage (currently 7 stages).
	real(r8), allocatable :: Transit(:)
	
	!insect density for different insect types and stages
	! Note that for the mountain pine beetle densities (insect type 1--first row)
	! are expressed as (nper 225 m^2). This does not have to be the case for other insect species.
        real(r8), allocatable :: indensity(:,:) 	 
	
        ! Variables related to mountain pine beetle winter survival (these are specific to the mountain pine beetle)
        real(r8) :: ColdestT                         			! coldest winter temperature experienced to date (resets on a yearly basis)
	real(r8) :: FebInPopn                         			    ! a cencus of population that is used to decide whether insects attack trees
	
	! Maximum and minimum daily temperatures in degree C
	real(r8) :: MaxDailyT                         			
 	real(r8) :: MinDailyT                          			
	contains
           procedure  :: InitInsectSite
           procedure  :: ZeroInsectSite
    
    end type ed_site_insect_type

    contains

    !==========================================================================================================
    subroutine InitInsectSite(this)

        implicit none

        ! argument
        class(ed_site_insect_type), intent(inout) :: this
	
	! Set up the bins based on the insect type (as each has a different n stages. 
	select case (insectType)
		case(1) !MPB
			allocate(this%InsectPFTPref(1:numberInsectTypes, 1:maxpft))
			this%InsectPFTPref = 0 
				this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle

			allocate(this%PhysAge(1:DomainSize,1:7))
			this%PhysAge(1:DomainSize, 1:7) = 0.0_r8
			
			allocate(this%Transit(1:7))
			this%Transit(1:7) = 0.0_r8
			
			allocate(this%indensity(1:numberInsectTypes, 1:maxNumStages))	
			this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8

			this%ColdestT = 15.0_r8
			
			this%FebInPopn = 0.0_r8	
		case(2) ! WPB
			allocate(this%InsectPFTPref(1:numberInsectTypes, 1:maxpft))
			this%InsectPFTPref = 0 
				this%InsectPFTPref(1,2)= 1					! Resolve this for WPB!!!!!!!!!!!!!!!!!!!!!!!!
			
			allocate(this%PhysAge(1:DomainSize,1:6))
			this%PhysAge(1:DomainSize, 1:6) = 0.0_r8
			
			allocate(this%Transit(1:6))
			this%Transit(1:6) = 0.0_r8
			
			allocate(this%indensity(1:numberInsectTypes, 1:maxNumStages))	
			this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8

			this%ColdestT = 15.0_r8
		case default
			fates_endrun("Missing InsectType Parameter")
	
	end select		
	! As model runs typically start January 1, 
	! I have decided to initialize with non-reactive temperatures for insects.
	this%MaxDailyT = 0.0_r8
	this%MinDailyT = 0.0_r8
    end subroutine InitInsectSite
    
    !==========================================================================================================
    subroutine ZeroInsectSite(this)
        implicit none
        ! argument
        class(ed_site_insect_type), intent(inout) :: this
	select case (insectType)
		case(1) !MPB
			! initialize the preference, but will move to the parameter file later
			this%InsectPFTPref = 0 
			this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle
			this%PhysAge(1:DomainSize, 1:7) = 0.0_r8
			this%Transit(1:7) = 0.0_r8
			this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8
			this%ColdestT = 0.0_r8
			this%FebInPopn = 0.0_r8
		case(2) ! WPB
		! initialize the preference, but will move to the parameter file later
			this%InsectPFTPref = 0 
			this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle
			this%PhysAge(1:DomainSize, 1:6) = 0.0_r8
			this%Transit(1:6) = 0.0_r8
			this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8
			this%ColdestT = 0.0_r8
			this%FebInPopn = 0.0_r8
			! As model runs typically start January 1, 
			! I have decided to initialize with non-reactive temperatures for insects.
			this%MaxDailyT = 0.0_r8
			this%MinDailyT = 0.0_r8
		case default
			fates_endrun("Missing InsectType Parameter")
	end select 

    end subroutine ZeroInsectSite
       
end module FatesInsectMemMod
