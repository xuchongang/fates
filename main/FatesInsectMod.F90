module FatesInsectMod
  use shr_kind_mod              , only : r8 => shr_kind_r8
  use FatesInterfaceMod         , only : bc_in_type
  use FatesInterfaceMod         , only : hlm_current_month, hlm_current_day, hlm_freq_day
  use EDtypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type
  use FatesInsectMemMod         , only : ed_patch_insect_type, numberInsectTypes

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
  private :: RBMortSim		! Regniere and Bentz model of mountain pine beetle winter mortality
  private :: MPBAttack		! simulates mountain pine beetle attack of various sized trees in a patch.

contains

  !========================================================================
  subroutine insect_model(currentSite, bc_in)
    !
    ! !DESCRIPTION:
    !  Core of insect model, calling all insect population demography routines.
    !  Although the insect model is called at the site level, all of the insect-
    !  related state variables are stored at the patch level (see FatesInsectMemMod)
    !
    use EDTypesMod           , only : AREA
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target  :: currentSite
    type(bc_in_type)        , intent(in)             :: bc_in

    ! patch pointer	
    type (ed_patch_type), pointer :: currentPatch
    ! cohort pointer
    type (ed_cohort_type), pointer :: currentCohort  
    
    ! For each site we cycle through the patches from oldest to youngest  
    currentPatch => currentSite%oldest_patch	! starting with the oldest 
    do while (associated(currentPatch))

    	!-----------------------------------------------------------------------
    	! zero out the insect mortality for the current day

    	currentCohort => currentPatch%tallest

    	do while(associated(currentCohort))

     		currentCohort%inmort = 0.0_r8

     		currentCohort => currentCohort%shorter

    	end do

    	!-----------------------------------------------------------------------
    	! Calling the insect demography submodels (currently only the mountain
    	! pine beetle model is implemented). Later there will be calls to other
    	! insect models that may attack different plant functional types.

	call beetle_model(currentSite, currentPatch, bc_in)

    	!-----------------------------------------------------------------------
    	! Check the total insect mortality. If the proportion if insect caused mortality in the cohort for
    	! the day is larger than one, we reset it to one (one corresponds to 365.0_r8 for a 365 day long year)
    	! as mortality is expressed on a per year basis but implemented on a daily basis.

    	currentCohort => currentPatch%tallest

    	do while(associated(currentCohort))

     		if(currentCohort%inmort*hlm_freq_day > 1.0_r8) then
     			currentCohort%inmort = 1.0_r8/hlm_freq_day
     		end if

     		currentCohort => currentCohort%shorter

    	end do ! Cohort do loop
	
	currentPatch => currentPatch%younger
	
     end do	! Patch do loop

  end subroutine insect_model

  !========================================================================
  subroutine beetle_model(currentSite, currentPatch, bc_in)
    !
    ! !DESCRIPTION:
    ! The mountain pine beetle model.
    !
    use EDTypesMod           , only : AREA
    use FatesInsectMemMod    , only : delta1, an, bn, ab, bb		! these parameters will be passed using parameter file.
    use FatesInsectMemMod    , only : ed_patch_insect_type
    use FatesInterfaceMod    , only : hlm_current_month, hlm_current_day, hlm_freq_day

    ! !ARGUMENTS:
    type(ed_site_type)       , intent(inout), target  :: currentSite
    type(ed_patch_type)      , intent(inout), target  :: currentPatch
    type(bc_in_type)         , intent(in)             :: bc_in

    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_insect_type), pointer :: pa_insect
    integer :: iofp                         	! index fates patch age

    ! Temperature variables that drive the mountain pine beetle demography model.
    real(r8) :: max_airTC                   	! maximum daily air temperature (degrees C) in the patch at reference height
    real(r8) :: min_airTC                   	! minimum daily air temperature (degrees C) in the patch at reference height

    !! Below are state variables that we track at the patch level.

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

    real(r8) :: Fec                         	! the expected number of pre-eggs at each time per 225 m^2
    real(r8) :: E                           	! the expected number of eggs at each time per 225 m^2
    real(r8) :: L1                          	! the expected number of L1 at each time step per 225 m^2
    real(r8) :: L2                          	! the expected number of L2 at each time step per 225 m^2
    real(r8) :: L3                          	! the expected number of L3 at each time step per 225 m^2
    real(r8) :: L4                          	! the expected number of L4 at each time step per 225 m^2
    real(r8) :: P                          	! the expected number of pupae at each time step per 225 m^2
    real(r8) :: Te                          	! the expected number of tenerals at each time step per 225 m^2
    real(r8) :: A                           	! the expected number of flying adults at each time step per 225 m^2
    real(r8) :: FA                          	! density of adults that initiated flight in the current time step per 225 m^2
    real(r8) :: Bt                		! beetles that remain in flight from the previous step per 225 m^2
    real(r8) :: Parents                     	! density of parent beetles in the current time step per 225 m^2

    ! The smallest probability of larval winter survival as a function of the lowest temperature to date.
    real(r8) :: PrS                       	! probability of winter survival
    real(r8) :: Ct                        	! The level of larval cold tolerance in the population.
    integer :: counter				! duration of cold hardening in the RBMortsim subroutine

    ! Current host tree densities for insects (in this case for mountain pine beetle) per 225 m^2 (15 X 15 m gap).
    real(r8) :: Nt68				! initial susceptible host trees in the 5 to 8 inch dbh size class
    real(r8) :: Nt10                    	! initial susceptible host trees in the 8 to 10 inch dbh size class
    real(r8) :: Nt12              		! initial susceptible host trees in the 10 to 12 inch dbh size class
    real(r8) :: Nt14              		! initial susceptible host trees in the 12 to 14 inch dbh size class
    real(r8) :: Nt16s             		! initial susceptible host trees in the  14 inch or larger dbh size class

    ! I also make the equivalent container for the density of hosts prior to insect attack so that we can compute
    ! the proportion that died in the current step (daily time step).
    real(r8) :: Ntm168            		! previous susceptible host trees in the 5 to 8 inch dbh size class
    real(r8) :: Ntm110            		! previous susceptible host trees in the 8 to 10 inch dbh size class
    real(r8) :: Ntm112            		! previous susceptible host trees in the 10 to 12 inch dbh size class
    real(r8) :: Ntm114            		! previous susceptible host trees in the 12 to 14 inch dbh size class
    real(r8) :: Ntm116s           		! previous susceptible host trees in the  14 inch or larger dbh size class

    ! Here are variables that I use to decide whether to restart the mountain pine beetle population at endemic population levels
    real(r8) :: InPopn            		! current total population of insects within trees (if measured before they fly)
    real(r8) :: FebInPopn         		! current total population of insects estimated on Feb. first (before they would fly)
    real(r8), parameter :: EndMPBPopn = 0.684_r8! The minimum endemic parent mountain pine beetle population (male and female) per 225 m^2.

    !----------------------------------------------------------------------------------------------------
    ! Grabbing the values of the state variables from currentPatch

    ! The physiological age distributions
    OE = currentPatch%pa_insect%MPB_PhysAge(:,1)
    OL1 = currentPatch%pa_insect%MPB_PhysAge(:,2)
    OL2 = currentPatch%pa_insect%MPB_PhysAge(:,3)
    OL3 = currentPatch%pa_insect%MPB_PhysAge(:,4)
    OL4 = currentPatch%pa_insect%MPB_PhysAge(:,5)
    OP = currentPatch%pa_insect%MPB_PhysAge(:,6)
    OT = currentPatch%pa_insect%MPB_PhysAge(:,7)

    ! The transitioning individuals from one life stage to another.
    NewEggstm1 = currentPatch%pa_insect%MPB_Transit(1)
    NewL1tm1 = currentPatch%pa_insect%MPB_Transit(2)
    NewL2tm1 = currentPatch%pa_insect%MPB_Transit(3)
    NewL3tm1 = currentPatch%pa_insect%MPB_Transit(4)
    NewL4tm1 = currentPatch%pa_insect%MPB_Transit(5)
    NewPtm1 = currentPatch%pa_insect%MPB_Transit(6)
    NewTtm1 = currentPatch%pa_insect%MPB_Transit(7)

    ! The one in the row argumuent of the indensity array corresponds to mountain pine beetle
    ! (insect type 1). The number in the column argument of the array refers to the
    ! life stage: 1->Fec, 2->Eggs, 3->L1, 4->L2, 5->L4, 6->L4, 7->Pupae,
    ! 8->Teneral adults, 9->Adults, 10->Flown Adults, 11->flying beetles from prvious step (Bt),
    ! 12->Parents (flown adults that succesfully attacked trees that day--daily time step)
    ! Columns 13 to 20 of the indensity array are empty for the mountain pine beetle.
    Fec = currentPatch%pa_insect%indensity(1,1)
    E = currentPatch%pa_insect%indensity(1,2)
    L1 = currentPatch%pa_insect%indensity(1,3)
    L2 = currentPatch%pa_insect%indensity(1,4)
    L3 = currentPatch%pa_insect%indensity(1,5)
    L4 = currentPatch%pa_insect%indensity(1,6)
    P = currentPatch%pa_insect%indensity(1,7)
    Te = currentPatch%pa_insect%indensity(1,8)
    A = currentPatch%pa_insect%indensity(1,9)
    FA = currentPatch%pa_insect%indensity(1,10)
    Bt = currentPatch%pa_insect%indensity(1,11)
    Parents = CurrentPatch%pa_insect%indensity(1,12)

    PrS = currentPatch%pa_insect%PrS
    Ct = currentPatch%pa_insect%Ct
    counter = currentPatch%pa_insect%counter

    !----------------------------------------------------------------------------------------------------
    ! Calculate the density trees in each of the size classes that we use in the mountain pine beetle model
    ! for the current patch. I then call the insect life cycle model at the patch level.
    Nt68 = 0.0_r8
    Nt10 = 0.0_r8
    Nt12 = 0.0_r8
    Nt14 = 0.0_r8
    Nt16s = 0.0_r8

    iofp = currentPatch%patchno             ! This is needed to get the relevant temperature variables from bc_in
    currentCohort => currentPatch%tallest

    do while(associated(currentCohort)) ! cycling through cohorts from tallest to shortest

	! The first row of the InsectPFTPref array corresponds to insect type. The one specifies that
	! this submodel is for the mountain pine beetle.
        if(currentPatch%pa_insect%InsectPFTPref(1,currentCohort%pft) == 1) then  ! only attacks pfts preferred by mountain pine beetle
          ! Below I compute the tree density per 225 m^2 in each of the size classes
          ! used in the current version of the insect mortality model.

          ! Here is the 5-8 inch dbh size class we use in the model.
          if(currentCohort%dbh >= 12.7_r8 .and. currentCohort%dbh < 20.32_r8)then
              Nt68 = Nt68 + currentCohort%n/AREA*225.0_r8
          end if

          ! Here is the 8-10 inch dbh size class we use in the model.
          if(currentCohort%dbh >= 20.32_r8 .and. currentCohort%dbh < 25.4_r8)then
              Nt10 = Nt10 + currentCohort%n/AREA*225.0_r8
          end if

          ! Here is 10-12 inch dbh size class we use in the model.
          if(currentCohort%dbh >= 25.4_r8 .and. currentCohort%dbh < 30.48_r8)then
              Nt12 = Nt12 + currentCohort%n/AREA*225.0_r8
          end if

          ! Here is 12-14 inch dbh size class we use in the model.
          if(currentCohort%dbh >= 30.48_r8 .and. currentCohort%dbh < 35.56_r8)then
              Nt14 = Nt14 + currentCohort%n/AREA*225.0_r8
          end if

          ! Here is 14 inch dbh size class and larger we use in the model.
          if(currentCohort%dbh >= 35.56_r8)then
              Nt16s = Nt16s + currentCohort%n/AREA*225.0_r8
          end if

        endif

        currentCohort => currentCohort%shorter

    end do ! This ends the cohort do loop

    ! converting the minimum and maximum daily air temperatures
    ! from degrees K to degrees C.
    max_airTC = bc_in(currentSite)%tgcm_max_pa(iofp) - 273.15_r8
    min_airTC = bc_in(currentSite)%tgcm_min_pa(iofp) - 273.15_r8

    ! I record the number of trees in each of the size classes prior to attack.
    Ntm168 = Nt68
    Ntm110 = Nt10
    Ntm112 = Nt12
    Ntm114 = Nt14
    Ntm116s = Nt16s

    !----------------------------------------------------------------------------------------------------
    ! Calling the full MPB simulation for the time step. Note that the air
    ! temperatures are recorded at the patch level.
    call MPBSim2(max_airTC, min_airTC, Parents, FA, OE, OL1, OL2, &
        OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
        NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
        Fec, E, L1, L2, L3, L4, P, Te, A, PrS, Ct, counter, &
        Nt68, Nt10, Nt12, Nt14, Nt16s, Bt, &
        an, bn, ab, bb, delta1)

    ! In the case of beetle extinction, we re-initialize the parent beetle population with
    ! a small number (endemic beetle population level) of parent beetles. We count the
    ! population in February so that we know that none have flown yet, but if it is exceedingly
    ! small, we re-initialize with parents on July 14 of the same year.
    InPopn = Fec + E + L1 + L2 + L3 + L4 + P + Te + A

    if(hlm_current_month == 2 .and. hlm_current_day == 1) then
        FebInPopn = InPopn
    end if

    if(hlm_current_month == 7 .and. hlm_current_day == 14 .and. FebInPopn < EndMPBPopn) then
        ! The endemic mountain pine beetle population per hectare was estimated by Carroll et al
        ! to be 15.2 attacks (female beetles) beetles = 30.4 beetles including male and female.
        ! I convert this to a density of beetles per 225 m^2 (value stored in EndMPBPopn parameter).
        Parents = EndMPBPopn
    end if

    !----------------------------------------------------------------------------------------------------
    ! update the vegetation mortality
    currentCohort => currentPatch%tallest

    ! Note that insect mortality is greater than zero only if the beetle population is
    ! larger than the endemic beetle population. Otherwise beetles only colonize trees
    ! that were already killed by other mortality causes so insect mortality is effectively zero.
    do while(associated(currentCohort)) ! cycling through cohorts from tallest to shortest
        ! Below I compute the tree density per 225 m^2 in each of the size classes
        ! used in the current version of the insect mortality model.

        if(FebInPopn > EndInPopn .and. currentPatch%pa_insect%InsectPFTPref(1,currentCohort%pft) == 1) then
            ! Here is the 5-8 inch dbh size class we use in the model.
	    ! In each dbhclass we multiply the daily probability of mortality by 365.0_r8
	    ! to the mortality rate on a yearly basis.
            if(currentCohort%dbh >= 12.7_r8 .and. currentCohort%dbh < 20.32_r8)then
                currentCohort%inmort = (1.0_r8 - Nt68/Ntm168)*365.0_r8
            end if

            ! Here is the 8-10 inch dbh size class we use in the model.
            if(currentCohort%dbh >= 20.32_r8 .and. currentCohort%dbh < 25.4_r8)then
                currentCohort%inmort = (1.0_r8 - Nt10/Ntm110)*365.0_r8
            end if

            ! Here is 10-12 inch dbh size class we use in the model.
            if(currentCohort%dbh >= 25.4_r8 .and. currentCohort%dbh < 30.48_r8)then
                currentCohort%inmort = (1.0_r8 - Nt12/Ntm112)*365.0_r8
            end if

            ! Here is 12-14 inch dbh size class we use in the model.
            if(currentCohort%dbh >= 30.48_r8 .and. currentCohort%dbh < 35.56_r8)then
                currentCohort%inmort = (1.0_r8 - Nt14/Ntm114)*365.0_r8
            end if

            ! Here is 14 inch dbh size class and larger we use in the model.
            if(currentCohort%dbh >= 35.56_r8)then
                currentCohort%inmort = (1.0_r8 - Nt16s/Ntm116s)*365.0_r8
            end if

        end if

        currentCohort => currentCohort%shorter

     end do ! This ends the cohort do loop

    !----------------------------------------------------------------------------------------------------
    !assign the updated values to the array for storage

     ! Mountain pine beetle densities in each life stage
     currentPatch%pa_insect%indensity(1,1) = Fec
     currentPatch%pa_insect%indensity(1,2) = E
     currentPatch%pa_insect%indensity(1,3) = L1
     currentPatch%pa_insect%indensity(1,4) = L2
     currentPatch%pa_insect%indensity(1,5) = L3
     currentPatch%pa_insect%indensity(1,6) = L4
     currentPatch%pa_insect%indensity(1,7) = P
     currentPatch%pa_insect%indensity(1,8) = Te
     currentPatch%pa_insect%indensity(1,9) = A
     currentPatch%pa_insect%indensity(1,10) = FA
     currentPatch%pa_insect%indensity(1,11) = Bt
     CurrentPatch%pa_insect%indensity(1,12) = Parents

     ! densities of individuals transitioning from one stage to another
     currentPatch%pa_insect%MPB_Transit(1) = NewEggstm1
     currentPatch%pa_insect%MPB_Transit(2) = NewL1tm1
     currentPatch%pa_insect%MPB_Transit(3) = NewL2tm1
     currentPatch%pa_insect%MPB_Transit(4) = NewL3tm1
     currentPatch%pa_insect%MPB_Transit(4) = NewL4tm1
     currentPatch%pa_insect%MPB_Transit(4) = NewPtm1
     currentPatch%pa_insect%MPB_Transit(4) = NewTtm1

     ! The physiological age distributions
     currentPatch%pa_insect%MPB_PhysAge(:,1) = OE
     currentPatch%pa_insect%MPB_PhysAge(:,2) = OL1
     currentPatch%pa_insect%MPB_PhysAge(:,3) = OL2
     currentPatch%pa_insect%MPB_PhysAge(:,4) = OL3
     currentPatch%pa_insect%MPB_PhysAge(:,5) = OL4
     currentPatch%pa_insect%MPB_PhysAge(:,6) = OP
     currentPatch%pa_insect%MPB_PhysAge(:,7) = OT

    ! The winter survival parameters
    currentPatch%pa_insect%PrS = PrS
    currentPatch%pa_insect%Ct = Ct
    currentPatch%pa_insect%counter = counter

end subroutine beetle_model

!==================================================================================================
Subroutine MPBSim2(Tmax, Tmin, Parents, FA, OE, OL1, OL2, &
            OL3, OL4, OP, OT, NewEggstm1, NewL1tm1, &
            NewL2tm1, NewL3tm1, NewL4tm1, NewPtm1, NewTtm1, &
            Fec, E, L1, L2, L3, L4, P, Te, A, PrS, Ct, counter, &
            Nt68, Nt10, Nt12, Nt14, Nt16s, Bt, &
            an, bn, ab, bb, delta1)
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

    real(r8), intent(inout) :: Fec            ! the expected number of pre-eggs at each time
    real(r8), intent(inout) :: E                ! the expected number of eggs at each time
    real(r8), intent(inout) :: L1               ! the expected number of L1 at each time step
    real(r8), intent(inout) :: L2               ! the expected number of L2 at each time step
    real(r8), intent(inout) :: L3               ! the expected number of L3 at each time step
    real(r8), intent(inout) :: L4               ! the expected number of L4 at each time step
    real(r8), intent(inout) :: P                ! the expected number of pupae at each time step
    real(r8), intent(inout) :: Te               ! the expected number of tenerals at each time step
    real(r8), intent(inout) :: A              ! the expected number of flying adults at each time step

    ! The smallest probability of larval winter survival as a function of the lowest temperature to date.
    real(r8), intent(inout) :: PrS
    real(r8), intent(inout) :: Ct             ! The level of larval cold tolerance in the population.
    integer(kind = 4), intent(inout) :: counter     ! the duration of cold tolerance accumulation

    ! input and output variables
    real(r8), intent(inout) :: Nt68                   ! initial susceptible host trees in the 5 to 8 inch dbh size class
    real(r8), intent(inout) :: Nt10                   ! initial susceptible host trees in the 8 to 10 inch dbh size class
    real(r8), intent(inout) :: Nt12                   ! initial susceptible host trees in the 10 to 12 inch dbh size class
    real(r8), intent(inout) :: Nt14                   ! initial susceptible host trees in the 12 to 14 inch dbh size class
    real(r8), intent(inout) :: Nt16s                  ! initial susceptible host trees in the  14 inch or larger dbh size class
    real(r8), intent(inout) :: Bt                     ! beetles that remain in flight from the previous step

    ! input parameters (dbh stands for tree diameter at breast height)
    real(r8), intent(in) :: an                        ! controls the tree loss rate as a function of tree size class
    real(r8), intent(in) :: bn                        ! controls the tree loss rate as a function of tree size class
    real(r8), intent(in) :: ab                        ! controls the beetle loss rate as a function of tree size class
    real(r8), intent(in) :: bb                        ! controls the beetle loss rate as a function of tree size class
    real(r8), intent(in) :: delta1                  ! the beetle settling rate per hour estimated in Goodsman et al (2016)

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

    real(r8) :: PrSurvNew         ! a container for the new candidate for larval minimum survival probability

    !--------------------------------------------------------------------------------------------------

    ! I reset survival probability if there are no larvae.
    if(L1 + L2 + L3 + L4 < 0.001)then
        PrS = 1.0
        Ct = 0.0
        counter = 0
    end if

    !! We need to compute the mean phloem temperature according to the model of
    !! Bolstad, Bentz and Logan (1997).
    !! We compute mean phloem temperature by averaging maximum and minimum phloem temperature.
    !! Here I use the average temperature differential (6.6 degrees C)
    Tmean = 0.5*(Tmax + Tmin) + 0.9 + 6.6*(Tmax - Tmin)/(2.0*24.4)
    Tmax2 = Tmax + 6.6*(Tmax - Tmin)/24.4
    Tmin2 = Tmin + 1.8

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
    mu1 = log(med1*deltat)  ! for eggs
    mu2 = log(med2*deltat)  ! for L1
    mu3 = log(med3*deltat)  ! for L2
    mu4 = log(med4*deltat)  ! for L3
    mu5 = log(med5*deltat)  ! for L4
    mu6 = log(med6*deltat)  ! for pupae
    mu7 = log(med7*deltat)  ! for teneral adults

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

    ! Applying larval mortality. Calling this subroutine
    ! produces a new minimum survival probability estimate and a
    ! new estimate of the current level of cold-hardiness.
    call RBMortSim(Tmin2, Tmax2, PrSurvNew, Ct, counter)
    ! Updating the minimum survival probability estimate
    PrS = min(PrS, PrSurvNew)
    ! larval mortality is only applied as individuals exit the
    ! larval stage because it depends on the minimum survival
    ! probability they experienced over their larval career.
    NewP = NewP*PrS

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

    ! Simulating the attack of host trees
    call MPBAttack(Nt68, Nt10, Nt12, Nt14, Nt16s, Bt, FA, Parents, &
            an, bn, ab, bb, delta1)
    ! This updates the density of trees in each of the size classes, and the density of beetles that remain in
    ! flight and outputs a number of parents that will start the oviposition process.
    
    contains
    !=================================================================================================================
subroutine MPBAttack(Nt68, Nt10, Nt12, Nt14, Nt16s, Bt, FA, Parents, an, bn, ab, bb, delta1)
    ! In this subroutine I solve the differential equations using the Euler method with an exceedingly small time step.

    implicit none

    ! input and output variables.
    ! Tree density in size classes per 225 m^2 (15m X 15m gap).
    real(r8), intent(inout) :: Nt68                 ! initial susceptible host trees in the 5 to 8 inch dbh size class
    real(r8), intent(inout) :: Nt10                 ! initial susceptible host trees in the 8 to 10 inch dbh size class
    real(r8), intent(inout) :: Nt12                 ! initial susceptible host trees in the 10 to 12 inch dbh size class
    real(r8), intent(inout) :: Nt14                 ! initial susceptible host trees in the 12 to 14 inch dbh size class
    real(r8), intent(inout) :: Nt16s                ! initial susceptible host trees in the  14 inch or larger dbh size class
    real(r8), intent(inout) :: Bt                   ! beetles that remain in flight from the previous step

    ! input variable
    real(r8), intent(in) :: FA                      ! Adults that just started to fly in this time step

    ! output variables
    real(r8), intent(out) :: Parents                ! the density of beetles that entered trees killed in this time step

    ! input parameters (dbh stands for tree diameter at breast height)
    real(r8), intent(in) :: an                      ! controls the tree loss rate as a function of dbh class
    real(r8), intent(in) :: bn                      ! controls the tree loss rate as a function of dbh class
    real(r8), intent(in) :: ab                      ! controls the beetle loss rate as a function of dbh class
    real(r8), intent(in) :: bb                      ! controls the beetle loss rate as a function of dbh class
    real(r8), intent(in) :: delta1                  ! the beetle settling rate per hour estimated in Goodsman et al (2016)

    ! Here are internal variables

    ! First we need to keep track of how many trees there initially were in each size class so that we can
    ! compute the density of infested trees at the end of the one day time step.
    real(r8) :: Btp1                      ! an updated value for the beetles
    real(r8) :: Ntp168                    ! updated susceptible host trees in the 5 to 8 inch dbh size class
    real(r8) :: Ntp110                    ! updated susceptible host trees in the 8 to 10 inch dbh size class
    real(r8) :: Ntp112                    ! updated susceptible host trees in the 10 to 12 inch dbh size class
    real(r8) :: Ntp114                    ! updated susceptible host trees in the 12 to 14 inch dbh size class
    real(r8) :: Ntp116s                   ! updated susceptible host trees in the  14 inch or larger dbh size class

    real(r8) :: Pt68                      ! parent beetles the 5 to 8 inch dbh size class
    real(r8) :: Pt10                      ! parent beetles in the 8 to 10 inch dbh size class
    real(r8) :: Pt12                      ! parent beetles in the 10 to 12 inch dbh size class
    real(r8) :: Pt14                      ! parent beetles in the 12 to 14 inch dbh size class
    real(r8) :: Pt16s                     ! parent beetles in the 14 inch or larger dbh size class

    real(r8) :: Ptp168                    ! updated parent beetles the 5 to 8 inch dbh size class
    real(r8) :: Ptp110                    ! updated parent beetles in the 8 to 10 inch dbh size class
    real(r8) :: Ptp112                    ! updated parent beetles in the 10 to 12 inch dbh size class
    real(r8) :: Ptp114                    ! updated parent beetles in the 12 to 14 inch dbh size class
    real(r8) :: Ptp116s                   ! updated parent beetles in the 14 inch or larger dbh size class

    ! parameters and variables related to the integration routine
    integer(kind = 4), parameter :: timesteps = 8640  ! number of 10 second time steps in a 24 hour day
    integer(kind = 4) :: i                            ! the iterator
    real(r8) :: deltat = 1.0/8640.0

    ! I add in the beetles that just started flying in the time step
    Bt = Bt + FA

    ! initializing the parent beetles
    Pt68 = 0.0
    Pt10 = 0.0
    Pt12 = 0.0
    Pt14 = 0.0
    Pt16s = 0.0

    !---------------------------------------------------------------------------------------------
    ! Here I do the Euler integration
    do i = 1,timesteps
        ! This is the ODE for flying beetles
        Btp1 = Bt - exp(ab + bb*6.5)*Nt68*deltat - exp(ab + bb*9.0)*Bt*Nt10*deltat - exp(ab + bb*11.0)*Bt*Nt12*deltat &
        - exp(ab + bb*13.0)*Bt*Nt14*deltat - exp(ab + bb*15.0)*Bt*Nt16s*deltat - delta1*24.0*Bt*deltat
        Ntp168 = Nt68 - exp(an + bn*6.5)*Bt*Nt68*deltat           ! This is the discretized ODE for 5-8 inch DBH trees
        Ntp110 = Nt10 - exp(an + bn*9.0)*Bt*Nt10*deltat           ! This is the discretized ODE for 8-10 inch DBH trees
        Ntp112 = Nt12 - exp(an + bn*11.0)*Bt*Nt12*deltat          ! This is the discretized ODE for 10-12 inch DBH trees
        Ntp114 = Nt14 - exp(an + bn*13.0)*Bt*Nt14*deltat          ! This is the discretized ODE for 12-14 inch DBH trees
        Ntp116s = Nt16s - exp(an + bn*15.0)*Bt*Nt16s*deltat       ! This is the discretized ODE for 14 inch and larger DBH trees

        Ptp168 = Pt68 + exp(ab + bb*6.5)*Bt*Nt68*deltat         ! This is the discretized ODE for parent beetles in 5-8 inch DBH trees
        Ptp110 = Pt10 + exp(ab + bb*9.0)*Bt*Nt10*deltat         ! This is the discretized ODE for parent beetles in 8-10 inch DBH trees
        Ptp112 = Pt12 + exp(ab + bb*11.0)*Bt*Nt12*deltat        ! This is the discretized ODE for parent beetles in 10-12 inch DBH trees
        Ptp114 = Pt14 + exp(ab + bb*13.0)*Bt*Nt14*deltat        ! This is the discretized ODE for parent beetles in 12-14 inch DBH trees
        Ptp116s = Pt16s + exp(ab + bb*15.0)*Bt*Nt16s*deltat     ! This is the discretized ODE for parent beetles in 14 inch + trees

        ! Now I update all of the state variables
        Bt = Btp1
        Nt68 = Ntp168
        Nt10 = Ntp110
        Nt12 = Ntp112
        Nt14 = Ntp114
        Nt16s = Ntp116s

        Pt68 = Ptp168
        Pt10 = Ptp110
        Pt12 = Ptp112
        Pt14 = Ptp114
        Pt16s = Ptp116s

    end do
    !------------------------------------------------------------------------------------------------

    ! To prevent the algorithm from returning NaNs or negative values.
    if(isnan(Bt) .or. Bt < 0.0)then
        Bt = 0.0
    end if

    if(isnan(Nt68) .or. Nt68 < 0.0)then
        Nt68 = 0.0
    end if

    if(isnan(Nt10) .or. Nt10 < 0.0)then
        Nt10 = 0.0
    end if

    if(isnan(Nt12) .or. Nt12 < 0.0)then
        Nt12 = 0.0
    end if

    if(isnan(Nt14) .or. Nt14 < 0.0)then
        Nt14 = 0.0
    end if

    if(isnan(Nt16s) .or. Nt16s < 0.0)then
        Nt16s = 0.0
    end if

    if(isnan(Pt68) .or. Pt68 < 0.0)then
        Pt68 = 0.0
    end if

    if(isnan(Pt10) .or. Pt10 < 0.0)then
        Pt10 = 0.0
    end if

    if(isnan(Pt12) .or. Pt12 < 0.0)then
        Pt12 = 0.0
    end if

    if(isnan(Pt14) .or. Pt14 < 0.0)then
        Pt14 = 0.0
    end if

    if(isnan(Pt16s) .or. Pt16s < 0.0)then
        Pt16s = 0.0
    end if

    ! Now I compute the number of new parents in this one day time interval.
    Parents = Pt68 + Pt10 + Pt12 + Pt14 + Pt16s

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
    real(r8), parameter :: fmax = 82.0      ! Regniere et al 2012 estimate that 82 eggs are produced per female

    ! Aplying winter mortality to egg laying adults
     if(Tmn2 <= -18.0)then
        Fec = 0.0
    end if

    ! Computing new eggs. Note this has to be done before updating the
    ! Fec variable below.
    NewEggs = Fec*(1.0 - exp(-med))

    ! Simulating oviposition: (Fec represents the number of eggs remaining)
    ! Below I assume that half of the individuals that fly and
    ! lay eggs are females (after tree induced colonization mortality)
    ! and each female lays an initial clutch of 82 eggs.
    Fec = 0.5*Parents*fmax + Fec*exp(-med)

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
    if(Tmn2 <= -18.0)then
        OldEPT = 0.0
    end if

    if(med > 0.0)then
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
            NewNext = 0.0
        end if

        ! Just to make sure that silly things don't happen
        if(NewNext < 0.0 .or. isnan(NewNext))then
            NewNext = 0.0
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
    integer(kind = 4), intent(in) :: n              ! size of aging domain
    real(r8), intent(in) :: avec(n)           ! The aging domain
    real(r8), intent(in) :: med               ! median development rate
    real(r8), intent(in) :: mu                ! mean of the log transformed random development rate
    real(r8), intent(in) :: sigma             ! scale parameter of the log-normally distributed rate
    real(r8), intent(in) :: NewL              ! New larvae (just developed from the previous stage)
    real(r8), intent(inout) :: NewLtm1        ! New larvae from the previous time step
    real(r8), intent(inout) :: OL(n)       ! Distribution of physiological age
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

    if(med > 0.0)then
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
        NewNext = 0.0
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
    if(Tmn2 <= -18.0)then
        Adtm1 = 0.0
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
    PaddedItemA(m+1:Twom) = 0.0

    ItemBmatrix = 0.0
    Convolved = 0.0

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
    DevR = 0.0

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
    Flying = 0.0

    ! Computing the proportion that are flying rate
    if(TC >= 17.53 .and. TC <= 42.00) then
        Flying = 2.500e+01 + (-5.324e+00)*TC + (4.277e-01)*(TC**2) + (-1.633e-02)*(TC**3) + &
        (3.014e-04)*(TC**4) + (-2.172e-06)*(TC**5)
    end if

end subroutine FlightFunc

!===============================================================================
Subroutine RBMortSim(Tmx2, Tmn2, PrSurv, Ct, counter)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The RBMortSim subroutine implements the Regniere and Bentz (2007)
    ! larval mortality simulation model.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(r8), intent(in) :: Tmx2              ! Daily maximum under bark temperature
    real(r8), intent(in) :: Tmn2              ! Daily minimum under bark temperature
    real(r8), intent(out) :: PrSurv           ! Lowest temperature to date.
    real(r8), intent(inout) :: Ct             ! current accumulated level of cold hardiness
    integer(kind = 4), intent(inout) :: counter     ! steps up by one every time the function is called
                                                    ! as long as there are larvae present.

    ! below are internal variables and parameters
    real(r8) :: CtNew                         ! The new level of cold hardiness

    ! Defining the constants
    ! the "parameter" indicates that
    ! these are constants.

    ! Parameters of the winter larval mortality model
    ! from Table 2 in Regniere and Bentz 2007
    ! the parameters are broken into groups below by
    ! equation in the original manuscript

    ! parameters for distribution of supercooling points (equation 1)
    real(r8), parameter :: alpha1 = -9.8       ! mean supercooling point (SCP) in state 1 in degree C
    real(r8), parameter :: Beta1 = 2.26        ! spread of SCP in state 1
    real(r8), parameter :: alpha2 = -21.2      ! mean SCP in state 2
    real(r8), parameter :: Beta2 = 1.47        ! spread of SCP in state 2
    real(r8), parameter :: alpha3 = -32.3      ! mean SCP in state 3
    real(r8), parameter :: Beta3 = 2.42        ! spread of SCP in state 3

    ! parameters for the gain rate (equation 3)
    real(r8), parameter :: rhoG = 0.311        ! maximum gain rate
    real(r8), parameter :: sigmaG = 8.716      ! spread of the gain temperature response

    ! parameters for the changing optimal temperature for gains in cold hardiness (equation 5)
    real(r8), parameter :: muG = -5.0
    real(r8), parameter :: kappaG = -39.3

    ! parameters for the loss rate (equation 4)
    real(r8), parameter :: rhoL = 0.791
    real(r8), parameter :: sigmaL = 3.251

    ! parameters for the changing optimal temperature for loss of cold hardiness (equation 6)
    real(r8), parameter :: muL = 33.9
    real(r8), parameter :: kappaL = -32.7

    ! These are the breakpoints that dictate what proportion of individuals are in each SCP state
    real(r8), parameter :: LambdaZero = 0.254
    real(r8), parameter :: LambdaOne = 0.764

    ! Variables involved in iteration
    real(r8) :: T               ! mean temperature at each time step in degrees Celcius
    real(r8) :: Range1          ! maximum range of temperature under the bark

    ! Defining other variables that are recomputed in each time step
    real(r8) :: TG              ! optimum temperature for gaining cold hardiness
    real(r8) :: TL              ! optimum temperature for losing cold hardiness
    real(r8) :: Gt              ! gain of cold hardiness at time t
    real(r8) :: Lt              ! loss of cold hardiness at time t
    real(r8) :: p1              ! the proportion of individuals in stage 1
    real(r8) :: p2              ! the proportion of individuals in stage 2
    real(r8) :: p3              ! the proportion of individuals in stage 3

    !---------------------------------------------------------------------------
    ! Now for the calculations

    ! Computing the max and min under-bark temperature as well as the range
    ! using the Bolstad, Bentz, and Logan (1997) model
    T = (Tmn2 + Tmx2)/2.0 ! The mean temperature
    Range1 = Tmx2 - Tmn2

    ! Computing the optimum gain and loss temperatures (equations 5 and 6).
    TG = muG + kappaG*Ct
    TL = muL + kappaL*Ct

    ! Computing the gain and loss functions for this time step (equations 3 and 4)
    Gt = Range1*rhoG*exp(-(T - TG)/sigmaG)/(sigmaG*(1 + exp(-(T - TG)/sigmaG))**2.0)
    Lt = Range1*rhoL*exp(-(T - TL)/sigmaL)/(sigmaL*(1 + exp(-(T - TL)/sigmaL))**2.0)

    ! Computing the amount of cold hardiness achieved in the step (equation 7)
    if(counter <= 153 .and. Ct < 0.5)then
        CtNew = Ct + (1 - Ct)*Gt
    else
         CtNew = Ct + (1 - Ct)*Gt - Ct*Lt
    end if

    ! Computing the proportion of larvae that are in each SCP stage (equation 9)
    ! This is modified from equation 9 according to advice from Jacques Regniere
    p1 = max(0.0,min(1.0,(0.5-CtNew)/(0.5-LambdaZero)))
    p3 = max(0.0,min(1.0,(CtNew-0.5)/(LambdaOne-0.5)))
    p2 = 1 - p1 - p3

    ! Now computing the probability of surviving in each time step
    ! Due to supercooling points, the probability of survival is just the
    ! cumulative density function of the corresponding logistic distribution
    ! function
    PrSurv = p1/(1 + exp(-(Tmn2 - alpha1)/Beta1)) + &
    p2/(1 + exp(-(Tmn2 - alpha2)/Beta2)) + &
    p3/(1 + exp(-(Tmn2 - alpha3)/Beta3))

    ! Now I update the cold tolerance level
    Ct = CtNew

    ! now I update the counter
    counter = counter + 1

End Subroutine RBMortSim

End Subroutine MPBSim2

end module FatesInsectMod
