#include "fabm_driver.h"

!------------------------------------------------light model-----------------------------------------------------
!
! This routine is to simulate the photosynthetic active radiation (PAR) attenuation in the water column:
!  
!           Att = (attSW + attCHL*chl)*z
!
! Adapted for GOTM-FABM by BW on May 2020 in Dalhousie university    
!
! Please go to this link for more information about GoTM-Fabm coding:
! https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model
!-----------------------------------------------------------------------------------------------------------------

MODULE memg_bio_fennel_light

   USE fabm_types

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(type_base_model), PUBLIC :: type_memg_bio_fennel_light
      ! Variable identifiers [model input]
      TYPE (type_dependency_id) :: id_chl,id_thick 
      TYPE (type_surface_dependency_id) :: id_I_0
      
      ! Identifiers for diagnostic variable [model output]
      TYPE (type_diagnostic_variable_id) :: id_par

      ! Parameters
      real(rk) :: attSW            ! PAR attenuation coef due to sea water [m-1]
      real(rk) :: attCHL           ! PAR attenuation coef due to chlorophyll [m-1]
      real(rk) :: PARfrac          ! fraction of photosynthetically available radiation [dimensionless]

   CONTAINS
      PROCEDURE :: initialize
      ! Reference model procedures here.
      PROCEDURE :: do
      
   END TYPE type_memg_bio_fennel_light

CONTAINS

   SUBROUTINE initialize(self,configunit)
      CLASS (type_memg_bio_fennel_light), INTENT(INOUT), TARGET :: self
      INTEGER,                          INTENT(IN)              :: configunit
 
      ! Parameters
      real(rk) :: attSW            ! PAR attenuation coef due to sea water [m-1]
      real(rk) :: attCHL           ! PAR attenuation coef due to chlorophyll [m-1]
      real(rk) :: PARfrac          ! fraction of photosynthetically available radiation [dimensionless]
      
      namelist /memg_bio_fennel_light/ attSW, attCHL, PARfrac
!
! set default values for model parameters
! 
      attSW = 0.04_rk            
      attCHL = 0.02486_rk           
      PARfrac = 0.43_rk          

      ! Read the namelist
      if (configunit>0) read(configunit,nml=memg_bio_fennel_light,err=99,end=100)

      ! Store parameter values into parameters
      ! Note: all rates are provided in values per day and then converted to values per second
      call self%get_parameter(self%attSW,'attSW','m-1','PAR attenuation coef due to sea water',default=attSW)
      call self%get_parameter(self%attCHL,'attCHL','m-1','PAR attenuation coef due to chlorophyll',default=attCHL)
      call self%get_parameter(self%PARfrac,'PARfrac','-','fraction of photosynthetically available radiation',default=PARfrac)
 
      ! Register state variables
 
      ! Register environmental dependencies
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_thick, standard_variables%cell_thickness)
      call self%register_dependency(self%id_chl,standard_variable=type_bulk_standard_variable(name='chlorophyll',units='mg m-3'))
      
      ! Register the contribution of all state variables to total nitrogen
         
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_par,'par','W m-2', 'photosynthetic active radiation in water column',  &
         &      standard_variable=type_bulk_standard_variable(name='photosynthetic_active_radiation_in_water_column',units='W m-2'))

   return

99 call self%fatal_error('memg_bio_fennel_light','Error reading namelist memg_bio_fennel_light.')

100 call self%fatal_error('memg_bio_fennel_light','Namelist memg_bio_fennel_light was not found.')

   ! Register model parameters and variables here.
   END SUBROUTINE initialize

   SUBROUTINE do(self,_ARGUMENTS_DO_)
      class (type_memg_bio_fennel_light),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
   
      real(rk) :: chl,I_0,thick,par
      real(rk) :: Itop,AttFac,Att,ExpAtt  ! local variables
      
      ! Obtain surface solar shortwave radiation
      _GET_SURFACE_(self%id_I_0,I_0)
      Itop = I_0*self%PARfrac     ! surface PAR

      AttFac = 0.0_rk;

      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_chl,chl) ! chlorophyll concentration
         _GET_(self%id_thick,thick) ! layer thinkness

         IF (I_0 .gt. 0.0_rk) THEN
            Att = (AttFac + self%attSW + self%attCHL*chl)*thick
            ExpAtt = EXP(-Att)
            par = Itop*(1.0_rk-ExpAtt)/Att   ! average at cell center
            Itop=Itop*ExpAtt ! PAR at the cell bottom (also the next cell surface)
         ELSE 
            par = 0
         END IF
         _SET_DIAGNOSTIC_(self%id_par,par) ! PAR at layer centre
      _VERTICAL_LOOP_END_

   END SUBROUTINE do

END MODULE


