module FatesInsectMemMod

    use FatesConstantsMod, only : r8 => fates_r8
    !use EDParamsMod             , only : insect_an

    implicit none

    ! Here are input parameter for the MPBAttack subroutine
    ! Need to move these parameters into the parameter file
    !real(r8),parameter :: an = insect_an                         	! controls tree loss rate as a function of beetle and tree density, default=-17.50628_r8_r8
    real(r8), parameter :: an = -17.50628_r8				! controls tree loss rate as a function of beetle and tree density
    real(r8), parameter :: ab = -12.06471_r8				! controls beetle loss rate as a function of beetle and tree density
    integer,  parameter :: DomainSize = 2**8                      	! domain size 256 categories for physiological ages
    integer,  parameter :: numberInsectTypes = 1                  	! number of insect types (currently only one-mountain pine beetle)
    integer,  parameter :: maxNumStages   = 20                    	! maximum number of stages for insect development
    integer,  parameter :: maxpft = 14					! maximum number of plant pfts (I've manually entered the current number, which is inelegant...)
    
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Defining a site-level type from which to obtain variables.
    type ed_site_insect_type
        !variables for tracking insect dynamics
	
	! array to define the preference of insects (0-not preferred; 1-preferred)
	integer, allocatable :: InsectPFTPref(:,:)	

        ! Containers for the distributions of physiological age for each life stage. 
	! The MPB_PhysAge array holds physiological age distributions for each of the relevant stages of the MPB.
	real(r8), allocatable :: MPB_PhysAge(:,:)		! array to hold physiological age for eggs, L1, L2, L3, L4, P, T
	
	! Containers for densities of insects transitioning from one stage to another for each insect type.
	! The MPB_Tranit array holds beetles trasitionaing into the egg life stage, into L1, into L2, into L3, 
	! into L4, into the pupal stage, into the teneral adult life stage (currently 7 stages).
	real(r8), allocatable :: MPB_Transit(:)
	
	!insect density for different insect types and stages
	! Note that for the mountain pine beetle densities (insect type 1--first row)
	! are expressed as (nper 225 m^2). This does not have to be the case for other insect species.
        real(r8), allocatable :: indensity(:,:) 	 
	
        ! Variables related to mountain pine beetle winter survival (these are specific to the mountain pine beetle)
        real(r8) :: PrS                         			! smallest probability of mountain pine beetle winter survival (unitless)
        real(r8) :: Ct                          			! The level of mountain pine beetle larval cold tolerance in the population (unitless).
	
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
	
	! initialize the preference, but will move to the parameter file later
	allocate(this%InsectPFTPref(1:numberInsectTypes, 1:maxpft))
	this%InsectPFTPref = 0 
        this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle

	allocate(this%MPB_PhysAge(1:DomainSize,1:7))
	this%MPB_PhysAge(1:DomainSize, 1:7) = 0.0_r8
	
	allocate(this%MPB_Transit(1:7))
	this%MPB_Transit(1:7) = 0.0_r8
	
	allocate(this%indensity(1:numberInsectTypes, 1:maxNumStages))	
	this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8

        this%PrS = 0.0_r8
        this%Ct = 0.0_r8
	
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
	
	! initialize the preference, but will move to the parameter file later
	this%InsectPFTPref = 0 
        this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle

	this%MPB_PhysAge(1:DomainSize, 1:7) = 0.0_r8
	
	this%MPB_Transit(1:7) = 0.0_r8
	
	this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8

        this%PrS = 0.0_r8
        this%Ct = 0.0_r8
	
	! As model runs typically start January 1, 
	! I have decided to initialize with non-reactive temperatures for insects.
	this%MaxDailyT = 0.0_r8
	this%MinDailyT = 0.0_r8

    end subroutine ZeroInsectSite
       
end module FatesInsectMemMod
