module FatesInsectMod
  use shr_kind_mod              , only : r8 => shr_kind_r8
  use FatesInterfaceMod         , only : bc_in_type
  use FatesInterfaceMod         , only : hlm_current_year, hlm_current_month, hlm_current_day, hlm_freq_day
  use EDtypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type
  use FatesInsectMemMod         , only : ed_site_insect_type, numberInsectTypes
  !use EDParamsMod               , only : insect_an

  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: insect_model

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: beetle_model	! calls the mountain pine beetle model
  private :: MPBSim2		! mountain pine beetle subroutine calls all of the affiliated MPB subroutines below
  private :: Ovipos		! mountain pine beetle oviposition subroutine
  private :: EPTDev		! development subroutine for mountain pine beetle eggs, pupae, and teneral adults
  private :: LarvDev		! mountain pine beetle larval development subroutine
  private :: AdDev		! mountain pine beetle adult development subroutine (unflown adults)
  private :: ConvolveSR		! Performs a convolution of two distributions (for insect dev.) using FFT
  private :: LnormPDF		! produces a lognormal distribution on a specified domain
  private :: RegnierFunc	! produces a hump shaped development rate curve for temperature-dependent insect development
  private :: FlightFunc		! temperature dependent flight probability for adult mountain pine beetles
  private :: MPBAttack		! simulates mountain pine beetle attack of various sized trees in a site.

contains

  !========================================================================
  subroutine insect_model(currentSite, bc_in)
    !
    ! !DESCRIPTION:
    !  Core of insect model, calling all insect population demography routines.
    !  The insect model is called at the site level, and all of the insect-
    !  related state variables are stored at the site level (see FatesInsectMemMod)
    !
    use FatesInterfaceMod    , only : bc_in_type
    use EDtypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target  :: currentSite
    type(bc_in_type)        , intent(in)	     :: bc_in

    !-----------------------------------------------------------------------
    ! calling the insect demography submodels (currently only the mountain
    ! pine beetle model is implemented). Later there will be calls to other
    ! insect models that may attack different plant functional types.
    !-----------------------------------------------------------------------
	
    call beetle_model(currentSite, bc_in)

  end subroutine insect_model

  !========================================================================
  subroutine beetle_model(currentSite, bc_in)
    !
    ! !DESCRIPTION:
    ! The mountain pine beetle model.
    !
    use FatesInsectMemMod    , only : an, ab, dd1			! these parameters will be passed using parameter file.
    use FatesInsectMemMod    , only : ed_site_insect_type
    use FatesInterfaceMod    , only : hlm_current_year, hlm_current_month, hlm_current_day, hlm_freq_day, bc_in_type
    use EDtypesMod           , only : ed_patch_type, ed_cohort_type

    ! !ARGUMENTS:
    type(ed_site_type)       , intent(inout), target  :: currentSite
    type(bc_in_type)         , intent(in)		:: bc_in
    !
    ! POINTERS TO GENERIC TYPES	
    type (ed_patch_type), pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    type(ed_site_insect_type), pointer :: si_insect
    integer :: iofp                         		! index fates patch age

    ! Temperature variables that drive the mountain pine beetle demography model. These are averages over all of the patches.
    ! at the specific site.
    real(r8) :: max_airTC                   	! maximum daily air temperature (degrees C) in the site at reference height
    real(r8) :: min_airTC                   	! minimum daily air temperature (degrees C) in the site at reference height

    !! Below are state variables that we track at the site level.

    ! Containers for the distributions of physiological age for each life stage. In the
    ! InitInsectSite subroutine these will be allocated with size equal to the domain size.
    real(r8) :: OE(2**8)           	! vector to hold physiological age distribution for eggs
    real(r8) :: OL1(2**8)          	! vector to hold physiological age distribution for first instar larvae
    real(r8) :: OL2(2**8)          	! vector to hold physiological age distribution for second instar larvae
    real(r8) :: OL3(2**8)          	! vector to hold physiological age distribution for third instar larvae
    real(r8) :: OL4(2**8)          	! vector to hold physiological age distribution for fourth instar larvae
    real(r8) :: OP(2**8)           	! vector to hold physiological age distribution for pupae
    real(r8) :: OT(2**8)           	! vector to hold physiological age distribution for teneral adults

    real(r8) :: NewEggstm1                  	! density of new eggs oviposited in the previous time step (t minus 1)
    real(r8) :: NewL1tm1                    	! density of new L1 in the previous time step (t minus 1)
    real(r8) :: NewL2tm1                    	! density of new L2 in the previous time step (t minus 1)
    real(r8) :: NewL3tm1                    	! density of new L3 in the previous time step (t minus 1)
    real(r8) :: NewL4tm1                    	! density of new L4 in the previous time step (t minus 1)
    real(r8) :: NewPtm1                     	! density of new pupae in the previous time step (t minus 1)
    real(r8) :: NewTtm1                     	! density of new teneral adults in the previous time step (t minus 1)

    real(r8) :: Fec                         	! the expected number of pre-eggs at each time per ha
    real(r8) :: E                           	! the expected number of eggs at each time per ha
    real(r8) :: L1                          	! the expected number of L1 at each time step per ha
    real(r8) :: L2                          	! the expected number of L2 at each time step per ha
    real(r8) :: L3                          	! the expected number of L3 at each time step per ha
    real(r8) :: L4                          	! the expected number of L4 at each time step per ha
    real(r8) :: P                          	! the expected number of pupae at each time step per ha
    real(r8) :: Te                          	! the expected number of tenerals at each time step per ha
    real(r8) :: A                           	! the expected number of flying adults at each time step ha
    real(r8) :: FA                          	! density of adults that initiated flight in the current time step per ha
    real(r8) :: Bt                		! beetles that remain in flight from the previous step per ha
    real(r8) :: Parents                     	! density of parent beetles in the current time step per ha

    ! related to the winter mortality model for mountain pine beetle:
    real(r8) :: ColdestT                       	! Coldest yearly temperature experienced to date.

    ! Current host tree densities for insects (in this case for mountain pine beetle) per ha
    ! Averaged over all patches within each site.
    real(r8) :: NtGEQ20				! initial susceptible host trees in the 20+ cm dbh size class

    ! I also make the equivalent container for the density of hosts prior to insect attack so that we can compute
    ! the proportion that died in the current step (daily time step).
    real(r8) :: Ntm1GEQ20            		! previous susceptible host trees in the 20+ cm dbh size class

    ! Here are variables that I use to decide whether to restart the mountain pine beetle population at endemic population levels
    real(r8) :: InPopn            		! current total population of insects within trees (if measured before they fly)
    real(r8) :: FebInPopn         		! current total population of insects estimated on Feb. first (before they would fly)
    real(r8), parameter :: EndMPBPopn = 40.0_r8 ! The minimum endemic parent mountain pine beetle population (male and female) per ha
    
    ! number of patches in the site
    integer :: NumPatches

    !----------------------------------------------------------------------------------------------------
    ! Grabbing the values of the state variables from currentSite

    ! The physiological age distributions
    OE = currentSite%si_insect%MPB_PhysAge(:,1)
    OL1 = currentSite%si_insect%MPB_PhysAge(:,2)
    OL2 = currentSite%si_insect%MPB_PhysAge(:,3)
    OL3 = currentSite%si_insect%MPB_PhysAge(:,4)
    OL4 = currentSite%si_insect%MPB_PhysAge(:,5)
    OP = currentSite%si_insect%MPB_PhysAge(:,6)
    OT = currentSite%si_insect%MPB_PhysAge(:,7)

    ! The transitioning individuals from one life stage to another.
    NewEggstm1 = currentSite%si_insect%MPB_Transit(1)
    NewL1tm1 = currentSite%si_insect%MPB_Transit(2)
    NewL2tm1 = currentSite%si_insect%MPB_Transit(3)
    NewL3tm1 = currentSite%si_insect%MPB_Transit(4)
    NewL4tm1 = currentSite%si_insect%MPB_Transit(5)
    NewPtm1 = currentSite%si_insect%MPB_Transit(6)
    NewTtm1 = currentSite%si_insect%MPB_Transit(7)

    ! The one in the row argumuent of the indensity array corresponds to mountain pine beetle
    ! (insect type 1). The number in the column argument of the array refers to the
    ! life stage: 1->Fec, 2->Eggs, 3->L1, 4->L2, 5->L4, 6->L4, 7->Pupae,
    ! 8->Teneral adults, 9->Adults, 10->Flown Adults, 11->flying beetles from prvious step (Bt),
    ! 12->Parents (flown adults that succesfully attacked trees that day--daily time step)
    ! Columns 13 to 20 of the indensity array are empty for the mountain pine beetle.
    Fec = currentSite%si_insect%indensity(1,1)
    E = currentSite%si_insect%indensity(1,2)
    L1 = currentSite%si_insect%indensity(1,3)
    L2 = currentSite%si_insect%indensity(1,4)
    L3 = currentSite%si_insect%indensity(1,5)
    L4 = currentSite%si_insect%indensity(1,6)
    P = currentSite%si_insect%indensity(1,7)
    Te = currentSite%si_insect%indensity(1,8)
    A = currentSite%si_insect%indensity(1,9)
    FA = currentSite%si_insect%indensity(1,10)
    Bt = currentSite%si_insect%indensity(1,11)
    Parents = CurrentSite%si_insect%indensity(1,12)

    ColdestT = currentSite%si_insect%ColdestT

    !----------------------------------------------------------------------------------------------------
    ! Calculate the site level average tree density in each of the size classes that we use in the 
    ! mountain pine beetle model across all patches. I then call the insect life cycle model at the site level.
    NtGEQ20 = 0.0_r8

    NumPatches = 0
    
    max_airTC = 0.0_r8
    min_airTC = 0.0_r8
    
    ! We cycle through the patches from oldest to youngest  
    currentPatch => currentSite%oldest_patch	! starting with the oldest 
    
    do while (associated(currentPatch))

	iofp = currentPatch%patchno             ! This is needed to get the relevant temperature variables from bc_in
    	currentCohort => currentPatch%tallest
	
	! Computing patch numbers
	NumPatches = NumPatches + 1
	
	! Computing mean temperature averaged across all patches (normalized later)
	max_airTC = max_airTC + (bc_in%tgcm_max_pa(iofp) - 273.15_r8 - 2.762601_r8)
    	min_airTC = min_airTC + (bc_in%tgcm_min_pa(iofp) - 273.15_r8 - 4.777561_r8)

   	do while(associated(currentCohort)) ! cycling through cohorts from tallest to shortest

        	! Below I compute the tree density per ha in each of the size classes
        	! used in the current version of the insect mortality model.

        	! Here is the 20+ cm dbh size class we use in the model. The fraction 
		! in parentheses ensures that when summed up over the whole site,
		! the density of trees will be the density per ha.
        	if(currentCohort%pft == 2 .and. currentCohort%dbh >= 20.0_r8)then
        		NtGEQ20 = NtGEQ20 + currentCohort%n*(currentPatch%area/10000.0_r8)
        	end if

        	currentCohort => currentCohort%shorter

    	end do ! This ends the cohort do loop
	
	currentPatch => currentPatch%younger
	
    end do	! Patch do loop
    
    ! Now completing the temperature averaging process.
    max_airTC = max_airTC/NumPatches
    min_airTC = min_airTC/NumPatches

    ! I record the number of trees in each of the size classes prior to attack.
    Ntm1GEQ20 = NtGEQ20
    
    ! In the case of beetle extinction, we re-initialize the parent beetle population with
    ! a small number (endemic beetle population level) of parent beetles. We count the
    ! population in February so that we know that none have flown yet, but if it is exceedingly
    ! small, we re-initialize with parents on July 14 of the same year.
    InPopn = Fec + E + L1 + L2 + L3 + L4 + P + Te + A

    if(hlm_current_month == 2 .and. hlm_current_day == 1) then
        FebInPopn = InPopn
    end if

    if(hlm_current_month == 7 .and. hlm_current_day == 21 .and. FebInPopn < EndMPBPopn) then
        ! The endemic mountain pine beetle population per hectare was estimated by Carroll et al
        ! to be 40 attacks (female beetles) per ha.
        Parents = EndMPBPopn
    end if
    
    ! Here's a hack to initialize the model with density of insects appropriate for the Rocky
    ! Mountain National Park.
    !if(hlm_current_year == 2000 .and. hlm_current_month == 7 .and. hlm_current_day == 21) then
        ! The intial number of parents in 2000 in Rocky Mountain National Park assuming that
	! each of 2.97418 killed trees per ha killed in 2000 was attacked by 173.216 females 
	! (our estimate of number of females per tree from Yellowstone dataset).
    !    Parents = 515.1756_r8
    !end if
    
    ! Here's a hack to initialize the model with density of insects appropriate for 
    ! Glacier National Park.
    if(hlm_current_year == 2005 .and. hlm_current_month == 7 .and. hlm_current_day == 21) then
        ! The intial number of parents in 2003 in Glacier National Park assuming that
	! each of 7.050161 killed trees per ha killed in 2005 was attacked by 173.216 females 
	! (our estimate of number of females per tree from Yellowstone dataset).
        Parents = 1221.201_r8
    end if
    
    ! Updating the coldest temperature
    if(min_airTC < ColdestT)then
    	ColdestT = Tmin
    end if
    
    ! Resetting the coldest temperature tracker on July 21 of each year:
    if(hlm_current_month == 7 .and. hlm_current_day == 21) then
       ColdestT = 15.0_r8
    end if

    !----------------------------------------------------------------------------------------------------
    ! Calling the full MPB simulation for the time step. 
    call MPBSim2(max_airTC, min_airTC, Parents, FA, OE, OL1, OL2, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, ColdestT, &
            NtGEQ20, Bt, an, ab, FebInPopn, EndMPBPopn, dd1)

    !----------------------------------------------------------------------------------------------------
    ! update the vegetation mortality.
    
    ! We cycle through the patches from oldest to youngest  
    currentPatch => currentSite%oldest_patch	! starting with the oldest 
    
    do while (associated(currentPatch))
    	currentCohort => currentPatch%tallest

    	! Note that insect mortality is greater than zero only if the beetle population is
    	! larger than the epidemic beetle population. Otherwise beetles only colonize trees
    	! that were already killed by other mortality causes so insect mortality is effectively zero.
    	do while(associated(currentCohort)) ! cycling through cohorts from tallest to shortest
        	! Below I compute the tree mortality rate (n/ha/year) in each of the size classes
        	! used in the current version of the insect mortality model.
	
        	! Here is the 20+ cm dbh size class we use in the model.
		! In each dbhclass we multiply the daily probability of mortality by 365.0_r8
		! to the mortality rate on a yearly basis.
        	if(FebInPopn > EndMPBPopn .and. currentCohort%pft == 2 .and. currentCohort%dbh >= &
			20.0_r8 .and. NtGEQ20 > 0.0_r8 .and. Ntm1GEQ20 > NtGEQ20)then
		
                		currentCohort%inmort = (1.0_r8 - NtGEQ20/Ntm1GEQ20)*365.0_r8	
			else
				currentCohort%inmort = 0.0_r8
        	end if

        	currentCohort => currentCohort%shorter

     	end do ! This ends the cohort do loop
	
	currentPatch => currentPatch%younger
	
    end do ! This ends the patch do loop

    !----------------------------------------------------------------------------------------------------
    !assign the updated values to the array for storage

     ! Mountain pine beetle densities in each life stage
     currentSite%si_insect%indensity(1,1) = Fec
     currentSite%si_insect%indensity(1,2) = E
     currentSite%si_insect%indensity(1,3) = L1
     currentSite%si_insect%indensity(1,4) = L2
     currentSite%si_insect%indensity(1,5) = L3
     currentSite%si_insect%indensity(1,6) = L4
     currentSite%si_insect%indensity(1,7) = P
     currentSite%si_insect%indensity(1,8) = Te
     currentSite%si_insect%indensity(1,9) = A
     currentSite%si_insect%indensity(1,10) = FA
     currentSite%si_insect%indensity(1,11) = Bt
     CurrentSite%si_insect%indensity(1,12) = Parents

     ! densities of individuals transitioning from one stage to another
     currentSite%si_insect%MPB_Transit(1) = NewEggstm1
     currentSite%si_insect%MPB_Transit(2) = NewL1tm1
     currentSite%si_insect%MPB_Transit(3) = NewL2tm1
     currentSite%si_insect%MPB_Transit(4) = NewL3tm1
     currentSite%si_insect%MPB_Transit(5) = NewL4tm1
     currentSite%si_insect%MPB_Transit(6) = NewPtm1
     currentSite%si_insect%MPB_Transit(7) = NewTtm1

     ! The physiological age distributions
     currentSite%si_insect%MPB_PhysAge(:,1) = OE
     currentSite%si_insect%MPB_PhysAge(:,2) = OL1
     currentSite%si_insect%MPB_PhysAge(:,3) = OL2
     currentSite%si_insect%MPB_PhysAge(:,4) = OL3
     currentSite%si_insect%MPB_PhysAge(:,5) = OL4
     currentSite%si_insect%MPB_PhysAge(:,6) = OP
     currentSite%si_insect%MPB_PhysAge(:,7) = OT

    ! The winter survival related variables.
    currentSite%si_insect%ColdestT = ColdestT
    
    ! Daily maximum and minimum temperatures for diagnostic purposes
    currentSite%si_insect%MaxDailyT = max_airTC
    currentSite%si_insect%MinDailyT = min_airTC

end subroutine beetle_model

!==================================================================================================
Subroutine MPBSim2(Tmax, Tmin, Parents, FA, OE, OL1, OL2, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, ColdestT, &
            NtGEQ20, Bt, an, ab, FebInPopn, EndMPBPopn, dd1)
    ! This subroutine simulates the demographic processes
    ! of the mountain pine beetle for a single time step including
    ! oviposition, the egg stage, the four larval instars,
    ! the pupal stage, the teneral adult stage,
    ! the adult stage, the flying adult stage, and then attack.

    implicit none

    ! Here are the input and output variables
    real(r8), intent(in) :: Tmax
    real(r8), intent(in) :: Tmin
    real(r8), intent(inout) :: Parents
    real(r8), intent(out) :: FA

    real(r8), intent(inout) :: OE(2**8)
    real(r8), intent(inout) :: OL1(2**8)
    real(r8), intent(inout) :: OL2(2**8)
    real(r8), intent(inout) :: OL3(2**8)
    real(r8), intent(inout) :: OL4(2**8)
    real(r8), intent(inout) :: OP(2**8)
    real(r8), intent(inout) :: OT(2**8)

    real(r8), intent(inout) :: NewEggstm1
    real(r8), intent(inout) :: NewL1tm1
    real(r8), intent(inout) :: NewL2tm1
    real(r8), intent(inout) :: NewL3tm1
    real(r8), intent(inout) :: NewL4tm1
    real(r8), intent(inout) :: NewPtm1
    real(r8), intent(inout) :: NewTtm1

    real(r8), intent(inout) :: Fec              ! the expected number of pre-eggs at each time
    real(r8), intent(inout) :: E                ! the expected number of eggs at each time
    real(r8), intent(inout) :: L1               ! the expected number of L1 at each time step
    real(r8), intent(inout) :: L2               ! the expected number of L2 at each time step
    real(r8), intent(inout) :: L3               ! the expected number of L3 at each time step
    real(r8), intent(inout) :: L4               ! the expected number of L4 at each time step
    real(r8), intent(inout) :: P                ! the expected number of pupae at each time step
    real(r8), intent(inout) :: Te               ! the expected number of tenerals at each time step
    real(r8), intent(inout) :: A                ! the expected number of flying adults at each time step

    ! The smallest probability of larval winter survival as a function of the lowest temperature to date.
    real(r8), intent(in) :: ColdestT

    ! input and output variables
    real(r8), intent(inout) :: NtGEQ20                ! initial susceptible host trees in the 20+ cm dbh size class
    real(r8), intent(inout) :: Bt                     ! beetles that remain in flight from the previous step

    ! input parameters
    real(r8), intent(in) :: an                        ! controls the tree loss rate
    real(r8), intent(in) :: ab                        ! controls the beetle loss rate
    real(r8), intent(in) :: FebInPopn                 ! February insect population
    real(r8), intent(in) :: EndMPBPopn 	      	      ! The endemic mountain pine beetle population (females per ha)
    real(r8), intent(in) :: dd1                       ! controls density dependent competition of juvenile mountain pine beetles

    !---------------------------------------------------------------------------------
    ! All of the parameters below are internal parameters (internal to the subroutine)

    ! iterator
    integer(kind = 4) :: j
    ! Defining the time step (1 day)

    real(r8), parameter :: deltat = 1.0 ! units are days

    !! Below I instantiate all of the rate constants. These
    !! are given in (from Regniere et al 2012) and were updated
    !! by personal communication with Jacques Regniere.

    ! parameters for oviposition rate
    real(r8), parameter :: sigma0 = 0.2458         ! controls rate variability
    real(r8), parameter :: TB0 = 4.6341            ! base temperature in degrees C
    real(r8), parameter :: DeltaB0 = 0.1
    real(r8), parameter :: TM0 = 27.7587           ! max temperature in degree C
    real(r8), parameter :: DeltaM0 = 3.0759
    real(r8), parameter :: omega0 = 0.3684
    real(r8), parameter :: psi0 = 0.005199

    ! parameters for development rate for eggs
    real(r8), parameter :: sigma1 = 0.1799         ! controls rate variability
    real(r8), parameter :: TB1 = 7.0               ! base temperature in degrees C
    real(r8), parameter :: DeltaB1 = 0.019297569
    real(r8), parameter :: TM1 = 30.0928           ! max temperature in degree C
    real(r8), parameter :: DeltaM1 = 4.4175
    real(r8), parameter :: omega1 = 0.2563
    real(r8), parameter :: psi1 = 0.02317

    ! parameters for development rate for L1 larvae
    real(r8), parameter :: sigma2 = 0.2911      ! controls rate variability
    real(r8), parameter :: TB2 = 3.5559
    real(r8), parameter :: DeltaB2 = 0.1
    real(r8), parameter :: TM2 = 29.2647
    real(r8), parameter :: DeltaM2 = 3.8227
    real(r8), parameter :: omega2 = 0.2398
    real(r8), parameter :: psi2 = 0.01082

    ! parameters for development rate for L2 larvae
    real(r8), parameter :: sigma3 = 0.3799       ! controls rate variability
    real(r8), parameter :: TB3 = 6.9598
    real(r8), parameter :: DeltaB3 = 0.097087379
    real(r8), parameter :: TM3 = 28.9047
    real(r8), parameter :: DeltaM3 = 3.0374
    real(r8), parameter :: omega3 = 0.3714
    real(r8), parameter :: psi3 = 0.01072

    ! parameters for development rate for L3 larvae
    real(r8), parameter :: sigma4 = 0.3868      ! controls rate variability
    real(r8), parameter :: TB4 = 6.8462
    real(r8), parameter :: DeltaB4 = 0.1
    real(r8), parameter :: TM4 = 28.7013
    real(r8), parameter :: DeltaM4 = 2.5359
    real(r8), parameter :: omega4 = 0.4399
    real(r8), parameter :: psi4 = 0.003892

    ! parameters for development rate for L4 larvae
    real(r8), parameter :: sigma5 = 0.3932       ! controls rate variability
    real(r8), parameter :: TB5 = 16.2464
    real(r8), parameter :: DeltaB5 = 0.039052889
    real(r8), parameter :: TM5 = 28.0
    real(r8), parameter :: DeltaM5 = 4.5504
    real(r8), parameter :: omega5 = 0.2593
    real(r8), parameter :: psi5 = 0.05034

    ! parameters for development rate for pupae
    real(r8), parameter :: sigma6 = 0.2998       ! controls rate variability
    real(r8), parameter :: TB6 = 5.63
    real(r8), parameter :: DeltaB6 = 0.10989011
    real(r8), parameter :: TM6 = 28.55
    real(r8), parameter :: DeltaM6 = 2.86
    real(r8), parameter :: omega6 = 0.1532
    real(r8), parameter :: psi6 = 0.02054

    ! parameters for development rate for teneral adults
    real(r8), parameter :: sigma7 = 0.5284       ! controls rate variability
    real(r8), parameter :: TB7 = 4.24
    real(r8), parameter :: DeltaB7 = 0.099967011
    real(r8), parameter :: TM7 = 35.0
    real(r8), parameter :: DeltaM7 = 7.1479
    real(r8), parameter :: omega7 = 0.1463
    real(r8), parameter :: psi7 = 0.01173
    
    ! Here are the parameters for the simplified model of larval mortality
    real(r8), parameter :: alpha3 = -35.4185	! The temperature in degrees centigrade at which only 50 % survival occurs
    real(r8), parameter :: Beta3 = 2.0		! controls the rate of change of survival probability as a function of temperature 

    ! Here are variables to hold the buffered under bark temperatures
    real(r8) :: Tmean     ! mean temperature at each time step in degrees Celcius
    real(r8) :: Tmin2     ! the buffered under-bark minimum temperature
    real(r8) :: Tmax2     ! the warmer under bark maximum temperature

    ! Variables that relate to the domain of the lognormal distribution
    integer(kind = 4), parameter :: n = 2**8    ! input variable. Must be specified
    real(r8), parameter :: Mx = 2.0
    real(r8), parameter :: Mn = 1.0e-20
    real(r8), parameter :: da = (Mx - Mn)/(2.0**8.0)

    real(r8) :: avec(n) = (/(Mn + j*da, j=0,n-1)/)          ! vector defining the domain

    ! parameters that change with each iteration
    real(r8) :: med0              ! median development rate in the pre-egg stage
    real(r8) :: med1              ! median development rate in egg stage
    real(r8) :: med2              ! median development rate in L1 stage
    real(r8) :: med3              ! median development rate in L2 stage
    real(r8) :: med4              ! median development rate in L3 stage
    real(r8) :: med5              ! median development rate in L4 stage
    real(r8) :: med6              ! median development rate in Pupa stage
    real(r8) :: med7              ! median development rate in teneral adult stage
    real(r8) :: mu1               ! mean of the log transformed random development rate in egg stage
    real(r8) :: mu2               ! mean of the log transformed random development rate in L1 stage
    real(r8) :: mu3               ! mean of the log transformed random development rate in L2 stage
    real(r8) :: mu4               ! mean of the log transformed random development rate in L3 stage
    real(r8) :: mu5               ! mean of the log transformed random development rate in L4 stage
    real(r8) :: mu6               ! mean of the log transformed random development rate in Pupa stage
    real(r8) :: mu7               ! mean of the log transformed random development rate in teneral adult stage

    ! New individuals that developed into the next life stage in the
    ! time step (these are each scalar values)
    real(r8) :: NewEggs
    real(r8) :: NewL1
    real(r8) :: NewL2
    real(r8) :: NewL3
    real(r8) :: NewL4
    real(r8) :: NewP
    real(r8) :: NewT
    real(r8) :: NewA

    !--------------------------------------------------------------------------------------------------

    !! We need to compute the mean phloem temperature according to the model of
    !! Bolstad, Bentz and Logan (1997).
    !! We compute mean phloem temperature by averaging maximum and minimum phloem temperature.
    !! Here I use the average temperature differential (6.6 degrees C)
    Tmean = 0.5_r8*(Tmax + Tmin) + 0.9_r8 + 6.6_r8*(Tmax - Tmin)/(2.0_r8*24.4_r8)
    Tmax2 = Tmax + 6.6_r8*(Tmax - Tmin)/24.4_r8
    Tmin2 = Tmin + 1.8_r8

    ! Computing the median development rate for each life stage in this time step
    call RegniereFunc(Tmean, TB0, DeltaB0, TM0, DeltaM0, omega0, psi0, med0)   ! for pre-eggs
    call RegniereFunc(Tmean, TB1, DeltaB1, TM1, DeltaM1, omega1, psi1, med1)   ! for eggs
    call RegniereFunc(Tmean, TB2, DeltaB2, TM2, DeltaM2, omega2, psi2, med2)   ! for L1
    call RegniereFunc(Tmean, TB3, DeltaB3, TM3, DeltaM3, omega3, psi3, med3)   ! for L2
    call RegniereFunc(Tmean, TB4, DeltaB4, TM4, DeltaM4, omega4, psi4, med4)   ! for L3
    call RegniereFunc(Tmean, TB5, DeltaB5, TM5, DeltaM5, omega5, psi5, med5)   ! for L4
    call RegniereFunc(Tmean, TB6, DeltaB6, TM6, DeltaM6, omega6, psi6, med6)   ! for pupae
    call RegniereFunc(Tmean, TB7, DeltaB7, TM7, DeltaM7, omega7, psi7, med7)   ! for teneral adults

    ! The mu parameter of the lognormal distribution is given by the natural logarithm of the median
    ! development rate.
    mu1 = 0.0_r8
    mu2 = 0.0_r8
    mu3 = 0.0_r8
    mu4 = 0.0_r8
    mu5 = 0.0_r8
    mu6 = 0.0_r8
    mu7 = 0.0_r8
    
    if(med1 > 0.0_r8) then
    	mu1 = log(med1*deltat)  ! for eggs
    end if 
    
    if(med2 > 0.0_r8) then
    	mu2 = log(med2*deltat)  ! for L1
    end if     
    
    if(med3 > 0.0_r8) then
    	mu3 = log(med3*deltat)  ! for L2
    end if 	

    if(med4 > 0.0_r8) then	
    	mu4 = log(med4*deltat)  ! for L3
    end if 	
	
    if(med5 > 0.0_r8) then	
    	mu5 = log(med5*deltat)  ! for L4
    end if 	
	
    if(med6 > 0.0_r8) then
    	mu6 = log(med6*deltat)  ! for pupae
    end if     
    
    if(med7 > 0.0_r8) then
    	mu7 = log(med7*deltat)  ! for teneral adults
    end if 

    !---------------------------------------------------------------------------------------------------------
    ! Now we can simulate each of the life stages by calling the appropriate subroutines

    ! Simulating oviposition:
    call Ovipos(Fec, Parents, med0, Tmin2, NewEggs)
    ! The output of this subroutine (NewEggs) is input for the next subroutine below.
    ! The Ovipos subroutine also updates the scalar value for fecundity.
    ! It takes as input the number of parents which comes from an initial value
    ! or from the MPBAttack subroutine called at the end of this sequence.

    ! Simulating egg development:
    call EPTDev(n, avec, med1, mu1, sigma1, Tmin2, NewEggs, NewEggstm1, OE, E, NewL1)
    ! The output of this subroutine (NewL1) is input for the next subroutine below.
    ! This updates the aging distribution (OE) and NewEggstm1 and
    ! outputs a scalar for the expected number of eggs.

    ! Simulating development of L1 larvae:
    call LarvDev(n, avec, med2, mu2, sigma2, NewL1, NewL1tm1, OL1, L1, NewL2)
    ! The output of this subroutine (NewL2) is input for the next subroutine below.
    ! This updates the aging distribution (OL1) and  NewL1tm1 and
    ! outputs a scalar for the expected number of first instar larvae (L1).

    ! Simulating development of L2 larvae:
    call LarvDev(n, avec, med3, mu3, sigma3, NewL2, NewL2tm1, OL2, L2, NewL3)
    ! The output of this subroutine (NewL3) is input for the next subroutine below.
    ! This updates the aging distribution (OL2) and  NewL2tm1 and
    ! outputs a scalar for the expected number of second instar larvae (L2).

    ! Simulating development of L3 larvae:
    call LarvDev(n, avec, med4, mu4, sigma4, NewL3, NewL3tm1, OL3, L3, NewL4)
    ! The output of this subroutine (NewL4) is input for the next subroutine below.
    ! This updates the aging distribution (OL3) and  NewL3tm1 and
    ! outputs a scalar for the expected number of third instar larvae (L3).

    ! Simulating development of L4 larvae:
    call LarvDev(n, avec, med5, mu5, sigma5, NewL4, NewL4tm1, OL4, L4, NewP)
    ! The output of this subroutine (NewP) is input for the next subroutine below.
    ! This updates the aging distribution (OL4) and  NewL4tm1 and
    ! outputs a scalar for the expected number of fourth instar larvae (L4).

    ! Applying larval mortality only as individuals exit the larval stage and
    ! develop into pupae because winter mortality depends on the coldest 
    ! temperature experienced over an individual's whole larval career. 
    ! Winter survival probability is modeled as a logistic curve function of 
    ! the coldest winter (air) temperature to date.
    NewP = NewP/(1.0_r8 + exp(-(ColdestT - alpha3)/Beta3))

    ! Simulating pupal development:
    call EPTDev(n, avec, med6, mu6, sigma6, Tmin2, NewP, NewPtm1, OP, P, NewT)
    ! The output of this subroutine (NewT) is input for the next subroutine below.
    ! This updates the aging distribution (OP) and  NewPtm1 and
    ! outputs a scalar for the expected number of pupae (P).

    ! Simulating teneral adult development:
    call EPTDev(n, avec, med7, mu7, sigma7, Tmin2, NewT, NewTtm1, OT, Te, NewA)
    ! The output of this subroutine (NewA) is input for the next subroutine below.
    ! This updates the aging distribution (OT) and  NewTtm1 and
    ! outputs a scalar for the expected number of teneral adults (Te).

    ! Simulating adult flight
    call AdSR(NewA, Tmin2, Tmax2, A, FA)
    ! This updates the expected number of adults (A) and flying adults (FA).
    
    ! Killing beetles that mistakenly remain in flight during cold temperatures.
    ! This prevents them from killing trees when they shouldn't be.
    if(Tmin <= -18.0)then
        FA = 0.0_r8
	Bt = 0.0_r8
    end if

    ! Simulating the attack of host trees
    call MPBAttack(NtGEQ20, Bt, FA, Parents, an, ab, FebInPopn, EndMPBPopn, dd1)
    ! This updates the density of trees in each of the size classes, and the density of beetles that remain in
    ! flight and outputs a number of parents that will start the oviposition process.
    
    contains
    !=================================================================================================================
subroutine MPBAttack(NtGEQ20, Bt, FA, Parents, an, ab, FebInPopn, EndMPBPopn, dd1)
    ! In this subroutine I solve the differential equations analytically.

    implicit none

    ! input and output variables.
    ! Tree density in size classes per ha
    real(r8), intent(inout) :: NtGEQ20              ! initial susceptible host trees in the 20+ cm dbh size class
    real(r8), intent(inout) :: Bt                   ! beetles that remain in flight from the previous step

    ! input variable
    real(r8), intent(in) :: FA                      ! Adults that just started to fly in this time step

    ! output variables
    real(r8), intent(out) :: Parents                ! the density of beetles that entered trees killed in this time step

    ! input parameters (dbh stands for tree diameter at breast height)
    real(r8), intent(in) :: an                      ! controls the tree loss rate
    real(r8), intent(in) :: ab                      ! controls the beetle loss rate 
    real(r8), intent(in) :: FebInPopn               ! February insect population
    real(r8), intent(in) :: EndMPBPopn              ! endemic mountain pine beetle population threshold
    real(r8), intent(in) :: dd1                     ! parameter controlling competition among juvenile beetles

    ! Here are internal variables and parameters
    real(r8) :: timestep = 1.0_r8         ! one day time step
    real(r8) :: Btp1                      ! an updated value for the beetles
    real(r8) :: Ntp1GEQ20                 ! updated susceptible host trees in the 20+ cm dbh size class
    real(r8) :: Atp1GEQ20                 ! updated parent beetles in the 20+ cm dbh size class
    real(r8) :: Itp1GEQ20                 ! updated infested trees in the 20+ cm dbh size class

    ! I add in the beetles that just started flying in the time step
    Bt = Bt + FA
    
    !---------------------------------------------------------------------------------------------
    ! Here I compute the analytic solutions

    ! Here's the solution for beetles
    Btp1 = Bt/(1.0_r8 + Bt*exp(ab)*timestep)

    ! Here's the analytic solution for attacking beetles
    Atp1GEQ20 = Bt*(Bt*exp(ab)*timestep/(1.0_r8 + Bt*exp(ab)*timestep))

    ! Here's the analytic solution for trees
    Ntp1GEQ20 = NtGEQ20*((1.0_r8 + Bt*exp(ab)*timestep)**(exp(an)/exp(ab)))*exp(-exp(an)*Bt*timestep)

    ! Here's the analytic solution for infested trees
    Itp1GEQ20 = NtGEQ20 - Ntp1GEQ20

    !------------------------------------------------------------------------------------------------
    
    ! Now I update all of the state variables. This depends on whether the population is endemic or not.
    ! when populations are in the endemic phase, they only attack weakened
    ! trees that are already functionally dead from other causes.
    ! 114363.64 is the surface area attacked by MPB in an average tree (from the Klein et al data)
    ! the formulation for Ricker type negative density dependence comes from Goodsman et al 2018.
    if(Itp1GEQ20 > 0.0_r8)then

        if(FebInPopn > EndMPBPopn)then
            Parents = Atp1GEQ20*exp(-dd1*sqrt(Atp1GEQ20/Itp1GEQ20/20.0*3.14159265359*((10.6)**2.0)/114363.64))
            Bt = Btp1
            NtGEQ20 = Ntp1GEQ20
            else
                ! Under the endemic scenario beetles do not kill trees.
                Parents = Atp1GEQ20*exp(-dd1*sqrt(Atp1GEQ20/Itp1GEQ20/20.0*3.14159265359*((10.6)**2.0)/114363.64))
                Bt = Btp1
                NtGEQ20 = NtGEQ20
        end if

    end if

end subroutine MPBAttack

!======================================================================================================
subroutine Ovipos(Fec, Parents, med, Tmn2, NewEggs)

    ! This subroutine simulate oviposition by parent beetles.
    ! The fecundity variable is the number of eggs remaining in the pool.
    ! med is median (temperature-dependent) oviposition rate.
    ! NewEggs are eggs that were laid in the time step.

    implicit none

    ! Here are the input and output variables
    real(r8), intent(inout) :: Fec            ! fecundity (eggs remaining to be laid)
    real(r8), intent(in) :: Parents           ! number of parents doing the ovipositing
    real(r8), intent(in) :: med               ! median oviposition rate
    real(r8), intent(in) :: Tmn2              ! minimum temperature under the bark
    real(r8), intent(out) :: NewEggs          ! Eggs laid in the time step

    ! internal parameters
    real(r8), parameter :: fmax = 54.66667_r8 ! Regniere et al 2012 estimate that 82 eggs are produced per female (2/3 of these are female)

    real(r8), parameter :: netp = 0.2526481_r8! This is one minus the net probability  of mortality from causes other
    					      ! than winter cold when the average winter mortality has been accounted for.

    ! Aplying winter mortality to egg laying adults
    if(Tmn2 <= -18.0)then
        Fec = 0.0_r8
    end if

    ! Computing new eggs. Note this has to be done before updating the
    ! Fec variable below.
    NewEggs = Fec*(1.0_r8 - exp(-med))*netp

    ! Simulating oviposition: (Fec represents the number of eggs remaining)
    ! each female lays an initial clutch of 82 eggs multiplied by two over three
    ! which is the proportion of the eggs that are female (the model only tracks
    ! female mountain pine beetles).
    Fec = Parents*fmax + Fec*exp(-med)

end subroutine Ovipos

!======================================================================================================
subroutine EPTDev(n, avec, med, mu, sigma, Tmn2, NewEPT, NewEPTtm1, OEPT, EPTcurrent, NewNext)

    ! This subroutine advances egg, pupal or teneral adult development and
    ! returns the number of individuals that move into the next stage
    ! (NewNext) in the time step as well as a count of the number
    ! of individuals currently in the stage (EPTcurrent).
    ! It can be used for the eggs, pupae or teneral adults, but it takes
    ! different life stage-dependent parameters (med, mu, sigma).

    ! This subroutine takes as input, n, which gives the size of the age domain,
    ! avec, which is the aging domain itself. The temperature-dependent
    ! median development rate (med) for the life stage, the corresponding
    ! mean of the log-normally distributed aging rate (mu),
    ! and the scale parameter of the log-normal distribution (sigma).
    ! Note that med, mu, and sigma are life stage-specific. The algorithm
    ! also takes as input the minimum daily temperature (Tmn2) which it uses
    ! to assign temperature-dependent mortality.
    ! The subroutine takes as input the number of new individuals in the
    ! current step (NewEPT), and from the previous step (NewEPTtm1).
    ! NewEPTtm1, OEPT is input and output and is the
    ! distribution of the accumulated physiological age in the life stage.

    implicit none

    ! Here are the input and output variables
    integer(kind = 4), intent(in) :: n              ! size of aging domain
    real(r8), intent(in) :: avec(n)           ! The aging domain
    real(r8), intent(in) :: med               ! median development rate
    real(r8), intent(in) :: mu                ! mean of the log transformed random development rate
    real(r8), intent(in) :: sigma             ! scale parameter of the log-normally distributed rate
    real(r8), intent(in) :: Tmn2              ! the buffered under-bark minimum temperature
    real(r8), intent(in) :: NewEPT            ! New individuals (just developed from previous stage)
    real(r8), intent(inout) :: NewEPTtm1      ! New individuals from the previous time step
    real(r8), intent(inout) :: OEPT(n)     ! Distribution of physiological age
    real(r8), intent(out) :: EPTcurrent       ! Number of individuals currently in the stage
    real(r8), intent(out) :: NewNext          ! New that just developed out of the current stage into the next one

    ! Here are variables related to the convolution
    real(r8) :: OldEPT(n)                  ! copy of the distribution of physiological age
    real(r8) :: AconvB(n)                  ! An array to hold the convolution result
    real(r8) PDF(n)                        ! An array to hold the aging kernel

    ! To compute the number of new individuals in the next stage,
    ! we need to know how many old individuals there were
    ! in the previous step.
    OldEPT = OEPT

    ! I use the -18 threshold to kill all eggs, pupae, or teneral
    ! adults as described in Regniere 2015.
    if(Tmn2 <= -18.0_r8)then
        OldEPT = 0.0
    end if

    if(med > 0.0_r8)then
        ! Computing the aging kernel
        call LnormPDF(avec, mu, sigma, PDF)
        PDF = PDF/sum(PDF)     ! ensuring that it sums to one

        ! Doing the convolution
        call ConvolveSR(PDF, OldEPT, AconvB)
        ! I assign all new individuals
        ! to the first age in the next stage.
        OEPT = AConvB + NewEPTtm1*PDF
        NewEPTtm1 = NewEPT

        ! Computing eggs
        EPTcurrent = sum(real(OEPT(1:128))) + NewEPT

        ! Computing NewNext: New individuals that developed into
        ! the next stage in this time step are individuals that
        ! exceeded the breakpoint of the current life stage (128).
        NewNext = sum(real(OEPT(129:n))) - sum(real(OldEPT(129:n)))

        else
            OEPT = OldEPT
            NewEPTtm1 = NewEPTtm1 + NewEPT
            EPTcurrent = sum(real(OEPT(1:128))) + NewEPTtm1

            ! Computing NewNext: if no development occurs there are no new individuals
            NewNext = 0.0_r8
        end if

        ! Just to make sure that silly things don't happen
        if(NewNext < 0.0 .or. isnan(NewNext))then
            NewNext = 0.0_r8
        end if

end subroutine EPTDev

!======================================================================================================
subroutine LarvDev(n, avec, med, mu, sigma, NewL, NewLtm1, OL, Lcurrent, NewNext)

    ! This subroutine advances larval development and returns
    ! the number of larvae in the current instar that advance
    ! to the next life stage (NewNext) in the time step as well as
    ! a count of the number larvae currently in the stage (Lcurrent).
    ! It can be used for any of the four mountain pine beetle larval
    ! instars but takes different instar-dependent parameters (med, mu, sigma).

    ! This subroutine takes as input, n, which gives the size of the age domain,
    ! avec, which is the aging domain itself. The temperature-dependent
    ! median development rate (med) for the larval stage
    ! (each instar has its own), the corresponding mean of the
    ! log-normally distributed aging rate (mu), and the scale parameter
    ! of the log-normal distribution (sigma). Note that mu and sigma
    ! are also larval instar-specific. The subroutine takes
    ! as input the number of new larvae in the current step (NewL)
    ! and from the previous step (NewLtm1). NewLtm1, and OL are input
    ! and output. OL is the distribution of the
    ! accumulated physiological age in the life stage.

    implicit none

    ! Here are the input and output variables
    integer(kind = 4), intent(in) :: n        ! size of aging domain
    real(r8), intent(in) :: avec(n)           ! The aging domain
    real(r8), intent(in) :: med               ! median development rate
    real(r8), intent(in) :: mu                ! mean of the log transformed random development rate
    real(r8), intent(in) :: sigma             ! scale parameter of the log-normally distributed rate
    real(r8), intent(in) :: NewL              ! New larvae (just developed from the previous stage)
    real(r8), intent(inout) :: NewLtm1        ! New larvae from the previous time step
    real(r8), intent(inout) :: OL(n)          ! Distribution of physiological age
    real(r8), intent(out) :: Lcurrent         ! Number of larvae
    real(r8), intent(out) :: NewNext          ! New individuals in the next stage (just developed out of current stage)

    ! Here are variables related to the convolution
    real(r8) :: OldL(n)                    ! copy of the distribution of physiological age
    real(r8) :: AconvB(n)                  ! An array to hold the convolution result
    real(r8) PDF(n)                        ! An array to hold the aging kernel

    ! To compute the number of new individuals in the next stage,
    ! we need to know how many old individuals there were
    ! in the previous step.
    OldL = OL

    if(med > 0.0_r8)then
        ! computing the aging kernel
        call LnormPDF(avec, mu, sigma, PDF)
        PDF = PDF/sum(PDF)     ! ensuring that it sums to one

        ! doing the convolution
        call ConvolveSR(PDF, OldL, AconvB)
        ! I assign all new individuals
        ! to the first age in the next stage.
        OL = AConvB + NewLtm1*PDF
        NewLtm1 = NewL

        ! computing Lcurrent
        Lcurrent = sum(real(OL(1:128))) + NewL

        ! Computing individuals that developed into the next stage
        ! in this time step: These are individuals that exceeded
        ! the breakpoint (128).
        NewNext = sum(real(OL(129:n))) - sum(real(OldL(129:n)))

    else
        OL = OldL
        NewLtm1 = NewLtm1 + NewL
        Lcurrent = sum(real(OL(1:128))) + NewLtm1

        ! Computing new individuals in the next stage
        NewNext = 0.0
    end if

    ! Just to make sure that silly things don't happen
    if(NewNext < 0.0 .or. isnan(NewNext))then
        NewNext = 0.0_r8
    end if

end subroutine LarvDev

!======================================================================================================
subroutine AdSR(NewA, Tmn2, Tmx2, Adtm1, FA)

    ! This subroutine keeps track of the number of mature adults
    ! that have not yet flown and outputs flown adults. This subroutine
    ! calls the PropFly subroutine to compute what proportion of adults
    ! fly as a function of maximum under-bark temperature.

    ! It takes as input a scalar (NewA) that comes from individuals that
    ! have developed out of the teneral adult stage (from the EPTDev subroutine),
    ! the scalar representing the minimum under bark temperature (Tmn2),
    ! the scalar representing the maximum under bark temperature(Tmx2),
    ! and a scalar number of adults in the previous time step Adtm1.
    ! It returns an updated scalar for the number of adults (Adtm1),
    ! and a scalar of the number of adults that flew in this timestep (FA).

    ! The number of adults that flew in the time step becomes input for a
    ! differential equation representing the attack of trees by the flying
    ! adults.

    implicit none

    ! Here are the input and output variables
    real(r8), intent(in) :: NewA              ! New adults
    real(r8), intent(in) :: Tmn2              ! minimum temperature under the bark
    real(r8), intent(in) :: Tmx2              ! maximum temperature under the bark
    real(r8), intent(inout) :: Adtm1          ! Number of adult individuals
    real(r8), intent(out) :: FA               ! Flying beetles

    ! An internal variable
    real(r8) :: PropFly                       ! The proportion of adult beetles that fly

    ! I use the -18 threshold to kill all adults as
    ! described in Regniere 2015
    if(Tmn2 <= -18.0_r8)then
        Adtm1 = 0.0_r8
    end if

    ! If it's warm enough, beetles fly. I'm using an empirical
    ! model based on the work of McCambridge (1971)
    call FlightFunc(Tmx2, PropFly)

    ! Computing adults that flew in this time step. Note that
    ! this needs to be done BEFORE updating the remaining adults.
    FA = Adtm1*PropFly

    ! Uptdating Adtm1
    Adtm1 = Adtm1*(1.0-PropFly) + NewA

end subroutine AdSR

!=============================================================================================================================
Subroutine ConvolveSR(ItemA, ItemB, AconvB)
    ! This subroutine convolves to arrays of the same dimensions (ItemA and ItemB)
    ! In this version I do the convolution using matrix multiplication rather than by
    ! using fast Fourier transforms as I had previously done even though the fft
    ! approach is more computationally efficient and elegant.

    implicit none

    real(r8), intent(in) :: ItemA(:)
    real(r8), intent(in) :: ItemB(:)
    real(r8), intent(out) :: AconvB(:)

    ! Parameters used in the convolution
    integer(kind = 4), parameter :: m = 2**8        ! input variable. Must be specified
    integer(kind = 4), parameter :: Twom = 2*m      ! Because I use a pad of size n, we double n

    ! Variables used in the convolution
    integer(kind = 4) :: i,j,k                      ! iterators
    real(r8) :: ItemBmatrix(Twom, Twom)       ! Here's a matrix for the item we are convolving
    real(r8) :: PaddedItemA(Twom)
    real(r8) :: Convolved(Twom)

    ! First I do the padding with zeros on the right hand side of the two input arrays
    PaddedItemA(1:m) = ItemA
    PaddedItemA(m+1:Twom) = 0.0_r8

    ItemBmatrix = 0.0_r8
    Convolved = 0.0_r8

    ! Now I fill in the elements of the matrix.
    do i = 1,m
       ItemBmatrix(i:(m + i - 1),i) = ItemB
    end do

    ! Now do the convolution using matrix multiplication.
    ! I do the matrix multiplication manually to minimize any
    ! dependencies on particular versions of Fortran.
    ! Because the ItemBmatrix is a diagonal matrix, we can
    ! do this more efficiently.
    do j = 1,Twom
        ! j loops over the rows
        ! we use the j index in the k do-loop below to
        ! capitalize on the diagonal nature of the ItemBmatrix
        do k = 1,j

            ! k loops over all of the columns
            Convolved(j) = Convolved(j) + ItemBmatrix(j,k)*PaddedItemA(k)

         end do ! ends the k loop over columns

    end do ! ends the j loop over rows

    ! To not loose individuals that developed past the threshold I add
    AconvB = Convolved(1:m)

    ! In this last funky step I take all of the individuals that developed into the
    ! padded region and add them to the last possible spot in the array so that they
    ! do not affect development in future steps.
    AconvB(m) = AconvB(m) + sum(Convolved(m+1:Twom))

End Subroutine ConvolveSR

!================================================================================================
Subroutine LnormPDF(x, mulog, sigmalog, PDF)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here is a subroutine that produces a lognormal probability density function
    ! over a range of values.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! I use the assumed-shape convention for the 1-D array
    implicit none

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: mulog
    real(r8), intent(in) :: sigmalog
    real(r8), intent(out) :: PDF(:)

    real(r8), parameter :: pi = 3.1415927

    PDF = 1.0/(x*sigmalog*sqrt(2.0*pi))*exp(-1.0/(2.0*sigmalog**2.0)*(log(x)-mulog)**2.0)

end Subroutine LnormPDF

!========================================================================================
subroutine RegniereFunc(TC, TB, DeltaB, TM, DeltaM, omega, psi, DevR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The RegniereFunc represents the Regniere function for
    ! temperature-dependent insect development (Regniere et al. 2012).
    ! TC is temperature in degrees Celcius.
    ! TB, DeltaB, TM,DeltaM, omega, psi are parameters
    ! The subroutine returns DevR, the median
    ! temperature-dependent development rate.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(r8), intent(in) :: TC
    real(r8), intent(in) :: TB
    real(r8), intent(in) :: DeltaB
    real(r8), intent(in) :: TM
    real(r8), intent(in) :: DeltaM
    real(r8), intent(in) :: omega
    real(r8), intent(in) :: psi
    real(r8), intent(out) :: DevR

    ! We start by defining it as zero in case the condition below does not hold
    DevR = 0.0_r8

    ! Computing the development rate
    if(TC >= TB .and. TC <= TM) then
        DevR = psi*(exp(omega*(TC - TB)) - (TM - TC)/(TM - TB)*exp(-omega*(TC - TB)/DeltaB) &
        - (TC - TB)/(TM - TB)*exp(omega*(TM - TB) - (TM - TC)/DeltaM))
    end if

end subroutine RegniereFunc

!===============================================================================
subroutine FlightFunc(TC, Flying)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Here's a subroutine that outputs the proportion of beetles that fly
    ! as a function of air temperature. The subroutine is based on an
    ! empirical fit of a polynomial function to the data of
    ! McCambridge (1971).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(r8), intent(in) :: TC
    real(r8), intent(out) :: Flying

    ! We start by defining it as zero in case the condition below does not hold
    Flying = 0.0_r8

    ! Computing the proportion that are flying rate
    if(TC >= 17.53_r8 .and. TC <= 42.00_r8) then
        Flying = 2.500e+01_r8 + (-5.324e+00_r8)*TC + (4.277e-01_r8)*(TC**2) + (-1.633e-02_r8)*(TC**3) + &
        (3.014e-04_r8)*(TC**4) + (-2.172e-06_r8)*(TC**5)
    end if

end subroutine FlightFunc

End Subroutine MPBSim2

end module FatesInsectMod
