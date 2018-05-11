module EDMortalityFunctionsMod

  ! ============================================================================
  ! Functions that control mortality.
  ! ============================================================================

   use FatesConstantsMod     , only : r8 => fates_r8
   use FatesGlobals          , only : fates_log
   use EDPftvarcon           , only : EDPftvarcon_inst
   use EDTypesMod            , only : ed_cohort_type
   use EDTypesMod            , only : ed_site_type
   use EDTypesMod            , only : ed_patch_type
   use FatesConstantsMod     , only : itrue,ifalse
   use FatesAllometryMod     , only : bleaf
   use FatesAllometryMod     , only : storage_fraction_of_target
   use FatesInterfaceMod     , only : bc_in_type
   use FatesInterfaceMod     , only : hlm_use_ed_prescribed_phys
   use FatesInterfaceMod     , only : hlm_freq_day
   use EDLoggingMortalityMod , only : LoggingMortality_frac
   use EDParamsMod           , only : fates_mortality_disturbance_fraction
   use FatesInterfaceMod     , only : bc_in_type


   implicit none
   private
   
   
   public :: mortality_rates
   public :: Mortality_Derivative
   
   logical :: DEBUG_growth = .false.
   
   ! ============================================================================
   ! 10/30/09: Created by Rosie Fisher
   ! 02/20/18: Refactored Ryan Knox
   ! ============================================================================


contains



  subroutine mortality_rates( cohort_in,bc_in,cmort,hmort,bmort,frmort,d13cmort)

    ! ============================================================================
    !  Calculate mortality rates from carbon storage, hydraulic cavitation, 
    !  background and freezing
    ! ============================================================================
    
    use FatesConstantsMod,  only : tfrz => t_water_freeze_k_1atm
   

    type (ed_cohort_type), intent(in) :: cohort_in 
    type (bc_in_type), intent(in) :: bc_in
    real(r8),intent(out) :: bmort ! background mortality : Fraction per year
    real(r8),intent(out) :: cmort  ! carbon starvation mortality
    real(r8),intent(out) :: hmort  ! hydraulic failure mortality
    real(r8),intent(out) :: frmort ! freezing stress mortality
    real(r8),intent(out) :: d13cmort  ! d13c related drought induced mortality, Hang ZHOU

    real(r8) :: frac  ! relativised stored carbohydrate
    real(r8) :: b_leaf ! target leaf biomass kgC
    real(r8) :: hf_sm_threshold    ! hydraulic failure soil moisture threshold 
    real(r8) :: temp_dep           ! Temp. function (freezing mortality)
    real(r8) :: temp_in_C          ! Daily averaged temperature in Celcius
    real(r8),parameter :: frost_mort_buffer = 5.0_r8  ! 5deg buffer for freezing mortality
    
    
    ! ! Hang ZHOU
    ! real(r8), parameter :: d13c_critical = -20.0_r8 ! -20 Liang Wei, threshold
    ! real(r8), parameter :: d13c_mortrate = 0.6_r8  !Liang Wei define rate
    ! integer ::&
    !      yr,    &! year
    !      mon,   &! month
    !      day,   &! day of month
    !      tod     ! time of day (seconds past 0Z)
    ! real(r8) :: d13c_background = 0.0_r8

    logical, parameter :: test_zero_mortality = .false. ! Developer test which
                                                        ! may help to debug carbon imbalances
                                                        ! and the like

    if (hlm_use_ed_prescribed_phys .eq. ifalse) then
    
    
    ! ! Hang ZHOU, Liang wei, calculate background d13c, this may need further edits as people may start modeling from year 1
    ! call get_curr_date(yr, mon, day, tod)
    ! if (yr < 1740) then
    !   d13c_background = -6.429_r8
    ! else if (yr > 2019) then
    !   d13c_background = -9.000_r8
    ! else
    !   d13c_background = -6.429_r8 - 0.0060_r8 * exp(0.0217_r8 * (yr - 1740))
    ! endif

    ! 'Background' mortality (can vary as a function of 
    !  density as in ED1.0 and ED2.0, but doesn't here for tractability) 
    bmort = EDPftvarcon_inst%bmort(cohort_in%pft) 

    ! Proxy for hydraulic failure induced mortality. 
    hf_sm_threshold = EDPftvarcon_inst%hf_sm_threshold(cohort_in%pft)

    if(cohort_in%patchptr%btran_ft(cohort_in%pft) <= hf_sm_threshold)then 
       hmort = EDPftvarcon_inst%mort_scalar_hydrfailure(cohort_in%pft)
     else
       hmort = 0.0_r8
     endif 
    
    ! Carbon Starvation induced mortality.
    if ( cohort_in%dbh  >  0._r8 ) then
       call bleaf(cohort_in%dbh,cohort_in%pft,cohort_in%canopy_trim,b_leaf)
       call storage_fraction_of_target(b_leaf, cohort_in%bstore, frac)
       if( frac .lt. 1._r8) then
          cmort = max(0.0_r8,EDPftvarcon_inst%mort_scalar_cstarvation(cohort_in%pft) * &
               (1.0_r8 - frac))
       else
          cmort = 0.0_r8
       endif

    else
       write(fates_log(),*) 'dbh problem in mortality_rates', &
            cohort_in%dbh,cohort_in%pft,cohort_in%n,cohort_in%canopy_layer
    endif
    
    ! ! D13C related drought induced mortality
    ! ! some quick output of the daily weighted mean d13c flux for debugging

    ! !if((d13c_background - cohort_in%c13disc_acc) >= d13c_critical)then
    ! if((d13c_background - cohort_in%c13disc_acc) >= d13c_critical .and. cohort_in%c13disc_acc /= 0)then
    !    d13cmort = d13c_mortrate
    ! else
    !    d13cmort = 0.0_r8
    ! endif

    ! Hang ZHOU (joeyzhou1984@gmail.com) 2018-03-25
    ! turn off all the calcualtion of d13cmort here, and move to EDAccumulatedFluxesMod
    ! d13cmort would be caclulated there and passed into the cohort data-structure
    ! so when the static mode is used, EDGrowthFUnctionsMode is not called, d13cmort will still be calculated
    ! and when the static mode is not used, `mortality_rates` can still access the calculated d13cmort via the cohort data-structure
    
    d13cmort = cohort_in%d13cmort
    !Liang Wei, temp: output to log file to test
    !if (DEBUG_growth) write(fates_log(), *) 'MORTALITY I, c13disc_acc', cohort_in%c13disc_acc
    !if (DEBUG_growth) write(fates_log(), *) 'MORTALITY II, d13cmort', d13cmort
    !if (DEBUG_growth) write(fates_log(), *) 'MORTALITY III, d13c_background', d13c_background
    !if (DEBUG_growth) write(fates_log(), *) 'MORTALITY IV, year', yr

    !    Mortality due to cold and freezing stress (frmort), based on ED2 and:           
    !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions 
    !           of land-use change, CO2 fertilization, and climate variability to the    
    !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              
    !           doi: 10.1111/j.1365-2486.2006.01254.x                                    

    temp_in_C = bc_in%t_veg24_si - tfrz
    temp_dep  = max(0.0,min(1.0,1.0 - (temp_in_C - &
                EDPftvarcon_inst%freezetol(cohort_in%pft))/frost_mort_buffer) )
    frmort    = EDPftvarcon_inst%mort_scalar_coldstress(cohort_in%pft) * temp_dep


    !mortality_rates = bmort + hmort + cmort

    else ! i.e. hlm_use_ed_prescribed_phys is true
       if ( cohort_in%canopy_layer .eq. 1) then
          bmort = EDPftvarcon_inst%prescribed_mortality_canopy(cohort_in%pft)
       else
          bmort = EDPftvarcon_inst%prescribed_mortality_understory(cohort_in%pft)
       endif
       cmort  = 0._r8
       hmort  = 0._r8
       frmort = 0._r8
    endif

    if (test_zero_mortality) then
       cmort = 0.0_r8
       hmort = 0.0_r8
       frmort = 0.0_r8
       bmort = 0.0_r8
    end if
       
    return
 end subroutine mortality_rates

 ! ============================================================================

 subroutine Mortality_Derivative( currentSite, currentCohort, bc_in)

    !
    ! !DESCRIPTION:
    ! Calculate the change in number density per unit time from the contributing
    ! rates.  These rates are not disturbance-inducing rates (that is handled
    ! elsewhere).
    !
    ! !USES:

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_cohort_type),intent(inout), target :: currentCohort
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)
    real(r8) :: frmort   ! freezing mortality rate (fraction per year)
    real(r8) :: d13cmort  ! d13c related drought induced mortality, Hang ZHOU
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    integer  :: ipft              ! local copy of the pft index
    !----------------------------------------------------------------------

    ipft = currentCohort%pft
    
    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    call mortality_rates(currentCohort,bc_in,cmort,hmort,bmort,frmort,d13cmort)
    call LoggingMortality_frac(ipft, currentCohort%dbh, &
                               currentCohort%lmort_direct,                       &
                               currentCohort%lmort_collateral,                    &
                               currentCohort%lmort_infra )

    if (currentCohort%canopy_layer > 1)then 
       
       ! Include understory logging mortality rates not associated with disturbance
       dndt_logging = (currentCohort%lmort_direct     + &
                       currentCohort%lmort_collateral + &
                       currentCohort%lmort_infra)/hlm_freq_day

       currentCohort%dndt = -1.0_r8 * (cmort+hmort+bmort+frmort+dndt_logging) * currentCohort%n
    else
       currentCohort%dndt = -(1.0_r8 - fates_mortality_disturbance_fraction) &
            * (cmort+hmort+bmort+frmort) * currentCohort%n
    endif

    return

 end subroutine Mortality_Derivative

end module EDMortalityFunctionsMod
