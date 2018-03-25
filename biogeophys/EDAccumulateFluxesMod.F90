module EDAccumulateFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This routine accumulates NPP, GPP and respiration of each cohort over the course of each 24 hour period.
  ! The fluxes are stored per cohort, and the npp_tstep (etc) fluxes are calcualted in EDPhotosynthesis
  ! This routine cannot be in EDPhotosynthesis because EDPhotosynthesis is a loop and therefore would
  ! erroneously add these things up multiple times.
  ! Rosie Fisher. March 2014.
  !
  ! !USES:
  use FatesGlobals, only      : fates_endrun
  use FatesGlobals, only      : fates_log
  use shr_log_mod , only      : errMsg => shr_log_errMsg
  use FatesConstantsMod , only : r8 => fates_r8

  use clm_time_manager   ,  only : get_curr_date  !Liang Wei

  implicit none
  private
  !
  public :: AccumulateFluxes_ED

  logical :: DEBUG = .false.  ! for debugging this module

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

contains

  !------------------------------------------------------------------------------

  subroutine AccumulateFluxes_ED(nsites, sites, bc_in, bc_out, dt_time)

    !
    ! !DESCRIPTION:
    ! see above
    !
    ! !USES:

    use EDTypesMod        , only : ed_patch_type, ed_cohort_type, &
                                   ed_site_type, AREA
    use FatesInterfaceMod , only : bc_in_type,bc_out_type

    !
    ! !ARGUMENTS
    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    type(bc_in_type),   intent(in)            :: bc_in(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)
    real(r8),           intent(in)            :: dt_time  ! timestep interval
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer  :: ccohort ! current cohort
    type(ed_patch_type) , pointer  :: cpatch ! current patch
    integer :: iv !leaf layer
    integer :: c  ! clm/alm column
    integer :: s  ! ed site
    integer :: ifp ! index fates patch
    real(r8):: n_perm2

    !Liang Wei 12/21/2017 for d13C estimation
    ! type(ed_cohort_type), pointer    :: currentcohort
    ! real(r8) :: gpp_tstep   !formerly gpp_clm (kgC/indiv/timestep)
    ! real(r8) :: gpp_acc
    ! real(r8) :: c13disc_acc
    ! real(r8) :: c13disc_clm

    real(r8) :: d13cmort = 0.0_r8                    ! d13c related drought induced mortality, Hang ZHOU
    real(r8), parameter :: d13c_critical = -20.0_r8  ! -20 Liang Wei, threshold
    real(r8), parameter :: d13c_mortrate = 0.6_r8    !Liang Wei define rate
    integer ::&
         yr,    &! year
         mon,   &! month
         day,   &! day of month
         tod     ! time of day (seconds past 0Z)
    real(r8) :: d13c_background = 0.0_r8

    !----------------------------------------------------------------------

    do s = 1, nsites

       ifp = 0
       sites(s)%npp = 0.0_r8
       cpatch => sites(s)%oldest_patch
       do while (associated(cpatch))
          ifp = ifp+1

          if( bc_in(s)%filter_photo_pa(ifp) == 3 ) then
             ccohort => cpatch%shortest
             do while(associated(ccohort))

                ! Accumulate fluxes from hourly to daily values.
                ! _tstep fluxes are KgC/indiv/timestep _acc are KgC/indiv/day

                if ( DEBUG ) then

                   write(fates_log(),*) 'EDAccumFlux 64 ',ccohort%npp_tstep
                   write(fates_log(),*) 'EDAccumFlux 66 ',ccohort%gpp_tstep
                   write(fates_log(),*) 'EDAccumFlux 67 ',ccohort%resp_tstep

                endif


                ! Hang ZHOU,  Liang Wei (2018-03-25)
                d13cmort = 0.0_r8
                d13c_background = 0.0_r8

                !if((ccohort%gpp_acc + ccohort%gpp_clm) .eq. 0.0_r8) then
                if((ccohort%gpp_acc + ccohort%gpp_tstep) .eq. 0.0_r8) then
                  ccohort%c13disc_acc = 0.0_r8
                else
                   !ccohort%c13disc_acc  = ((ccohort%c13disc_acc * ccohort%gpp_acc) + (ccohort%c13disc_clm * ccohort%gpp_clm)) / &
                        !(ccohort%gpp_acc + ccohort%gpp_clm)
                    ccohort%c13disc_acc  = ((ccohort%c13disc_acc * ccohort%gpp_acc) + (ccohort%c13disc_clm * ccohort%gpp_tstep)) / &
                                      (ccohort%gpp_acc + ccohort%gpp_tstep)
                endif

                ccohort%npp_acc  = ccohort%npp_acc  + ccohort%npp_tstep
                ccohort%gpp_acc  = ccohort%gpp_acc  + ccohort%gpp_tstep
                ccohort%resp_acc = ccohort%resp_acc + ccohort%resp_tstep

                ! Hang ZHOU (2018-03-25)
                ! caclulate the d13cmort and saved to ccohor (current cohort)
                ! the related snippets stay in the `EDGrowthFunctionsMod/mortality_rates` before
                ! they are moved here to make sure d13cmore will still be calculated even the static option is used
                ! there maybe better way to organize the code structure and logic

                call get_curr_date(yr, mon, day, tod)
                if (yr < 1740) then
                  d13c_background = -6.429_r8
                else if (yr > 2019) then
                  d13c_background = -9.000_r8
                else
                  d13c_background = -6.429_r8 - 0.0060_r8 * exp(0.0217_r8 * (yr - 1740))
                endif
                if((d13c_background - ccohort%c13disc_acc) >= d13c_critical .and. ccohort%c13disc_acc /= 0)then
                   d13cmort = d13c_mortrate
                else
                   d13cmort = 0.0_r8
                endif
                ccohort%d13cmort = d13cmort

                !----- THE FOLLOWING IS ONLY IMPLEMENTED TEMPORARILY FOR B4B reproducibility
                !----- ALSO, THERE IS NO REASON TO USE THE ISNEW FLAG HERE
                if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                   n_perm2   = ccohort%n/AREA
                else
                   n_perm2   = 0.0_r8
                endif
                if ( .not. ccohort%isnew ) then
                   sites(s)%npp = sites(s)%npp + ccohort%npp_tstep * n_perm2 * 1.e3_r8 / dt_time
                endif

                do iv=1,ccohort%nv
                   if(ccohort%year_net_uptake(iv) == 999._r8)then ! note that there were leaves in this layer this year.
                      ccohort%year_net_uptake(iv) = 0._r8
                   end if
                   ccohort%year_net_uptake(iv) = ccohort%year_net_uptake(iv) + ccohort%ts_net_uptake(iv)
                enddo

                ccohort => ccohort%taller
             enddo ! while(associated(ccohort))
          end if
          cpatch => cpatch%younger
       end do  ! while(associated(cpatch))
    end do
    return

 end subroutine AccumulateFluxes_ED

end module EDAccumulateFluxesMod

