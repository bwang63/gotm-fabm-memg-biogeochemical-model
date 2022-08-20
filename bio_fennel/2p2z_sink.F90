#include "fabm_driver.h"

!------------------------------------------------- 2p2z_sink -----------------------------------------------------
!
! This routine is to calculate the sinking velocity and sinking flux of phytoplankton and detritus. Following PISCES
! (Aumont et al 2015), the sinking velocity of small detritus is assumed to be constant while the sinking velocity of
! large detritus increases linearly with the depth.
!
! wL = wL_min + (wL_max - wL_min)*max(0,depth-Zref)/Zmax
!
!
! Adapted for GOTM-FABM by Bin Wang on May 2020 in Dalhousie university
!
! references:
!  Aumont, O., Ethé, C., Tagliabue, A., Bopp, L., and Gehlen, M.: PISCES-v2: an ocean biogeochemical model for carbon
!       and ecosystem studies, Geosci. Model Dev., 8, 2465–2513, https://doi.org/10.5194/gmd-8-2465-2015, 2015.
!
! Please go to this link for more information about GoTM-Fabm coding:
! https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model
!
!-----------------------------------------------------------------------------------------------------------------

MODULE memg_bio_fennel_2p2z_sink

   USE fabm_types

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(type_base_model), PUBLIC :: type_memg_bio_fennel_2p2z_sink
      ! Variable identifiers
      TYPE (type_state_variable_id) :: id_phyS, id_phyL
      TYPE (type_state_variable_id) :: id_chlS, id_chlL
      TYPE (type_state_variable_id) :: id_SDeN, id_LDeN
      TYPE (type_state_variable_id) :: id_Opal, id_Calc
      
      
      TYPE (type_dependency_id)            :: id_dep, id_thick ! depth of the middle of cell, thickness of the cell  
      
      TYPE (type_diagnostic_variable_id)   :: id_wL ! sinking velocity of large detritus

      ! Parameters
      real(rk)   ::  wPS, wPL, wS ! sinking velocities (prescribed in fabm.ymal)
      real(rk)   ::  wL_min, wL_max ! minimum and maximum sinking velocity of large detritus
      
   CONTAINS
      PROCEDURE :: initialize
      ! Reference model procedures here.
      PROCEDURE :: get_vertical_movement
   END TYPE type_memg_bio_fennel_2p2z_sink

CONTAINS

   SUBROUTINE initialize(self,configunit)
      CLASS (type_memg_bio_fennel_2p2z_sink), INTENT(INOUT), TARGET :: self
      INTEGER,                          INTENT(IN)              :: configunit
      
      ! Parameters
      real(rk)   ::  wPS, wPL, wS ! sinking velocities (prescribed in fabm.ymal)
      real(rk)   ::  wL_min, wL_max ! minimum and maximum sinking velocity of large detritus
      
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      namelist /memg_bio_fennel_2p2z_sink/ wPS, wPL, wS, wL_min, wL_max
!
! set default values for model parameters
! 
      wPS  = 0.1_rk 
      wPL  = 0.1_rk
      wS  = 0.1_rk
      wL_min  = 10.0_rk
      wL_max  = 200.0_rk
      
      ! Read the namelist
      if (configunit>0) read(configunit,nml=memg_bio_fennel_2p2z_sink,err=99,end=100)
      
        ! Store parameter values into parameters
        ! Note: all sinking velocities are provided in values per day and then converted to values per second
        call self%get_parameter(self%wPS,'wPS','m d-1','Sinking velocity of small phytoplankton (<0 for sinking)',default=wPS,scale_factor=-1.0_rk*d_per_s)
        call self%get_parameter(self%wPL,'wPL','m d-1','Sinking velocity of large phytoplankton (<0 for sinking)',default=wPL,scale_factor=-1.0_rk*d_per_s)
        call self%get_parameter(self%wS,'wS','m d-1','Sinking velocity of samll detritus (<0 for sinking)',default=wS,scale_factor=-1.0_rk*d_per_s)
        call self%get_parameter(self%wL_min,'wL_min','m d-1','Minimum sinking velocity of large detritus (<0 for sinking)',default=wL_min,scale_factor=-1.0_rk*d_per_s)
        call self%get_parameter(self%wL_max,'wL_max','m d-1','Maximum sinking velocity of large detritus (<0 for sinking)',default=wL_max,scale_factor=-1.0_rk*d_per_s)

        ! Register state variables
        call self%register_state_variable(self%id_phyS,'phyS','mmol N m-3','small_phytoplankton', minimum=0.0_rk,   &
           &    standard_variable=type_bulk_standard_variable(name='small_phytoplankton',units='mmol N m-3'))
        call self%register_state_variable(self%id_phyL,'phyL','mmol N m-3','large_phytoplankton', minimum=0.0_rk,   &
           &    standard_variable=type_bulk_standard_variable(name='large_phytoplankton',units='mmol N m-3'))
     
        call self%register_state_variable(self%id_chlS,'chlS','mg m-3','small_chlorophyll', minimum=0.0_rk,    &
           &    standard_variable=type_bulk_standard_variable(name='small_chlorophyll',units='mg m-3'))
        call self%register_state_variable(self%id_chlL,'chlL','mg m-3','large_chlorophyll', minimum=0.0_rk,    &
           &    standard_variable=type_bulk_standard_variable(name='large_chlorophyll',units='mg m-3'))
      
        call self%register_state_variable(self%id_LDeN,'LDeN','mmol N m-3','LdetritusN', minimum=0.0_rk,   &
           &    standard_variable=type_bulk_standard_variable(name='LdetritusN',units='mmol N m-3'))
        call self%register_state_variable(self%id_SDeN,'SDeN','mmol N m-3','SdetritusN', minimum=0.0_rk,   &
           &    standard_variable=type_bulk_standard_variable(name='SdetritusN',units='mmol N m-3'))
         
        call self%register_state_variable(self%id_Opal,'Opal','mmol Si m-3','Opal', minimum=0.0_rk,   &
           &    standard_variable=type_bulk_standard_variable(name='Opal',units='mmol Si m-3'))
        call self%register_state_variable(self%id_Calc,'Calcite','mmol Ca m-3','Calcite',minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='Calcite',units='mmol Ca m-3'))
         
        ! Register environmental dependencies
        call self%register_dependency(self%id_dep, standard_variables%depth) ! Depth of middle of cell
        call self%register_dependency(self%id_thick, standard_variables%cell_thickness) ! Thickness of cell

        ! Register diagnostic variables
        call self%register_diagnostic_variable(self%id_wL,'wL','m d-1',  'Sinking velocity of large detritus (>0 for sinking)',   &
            &  standard_variable=type_bulk_standard_variable(name='wL',units='m d-1'))

   return

99 call self%fatal_error('memg_bio_fennel_2p2z_sink','Error reading namelist memg_bio_fennel_2p2z_sink.')

100 call self%fatal_error('memg_bio_fennel_2p2z_sink','Namelist memg_bio_fennel_2p2z_sink was not found.')

   END SUBROUTINE initialize

   SUBROUTINE get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_memg_bio_fennel_2p2z_sink),intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: phyS,phyL,chlS,chlL,LDeN,SDeN
      real(rk) :: wL
      real(rk) :: thick,depth
      real(rk) :: zref,zfact,zmax
      real(rk), parameter :: secs_per_day = 86400.0_rk
      
      _VERTICAL_LOOP_BEGIN_

         ! Obtain Dependency and State Variable Values for the specific depth of the loop
         _GET_(self%id_thick, thick) !Thickness of cell [m]
         _GET_(self%id_dep, depth) !Depth of middle of cell [m]

         ! Sinking velocity of large detritus is assumed to increase linearly with the depth
         ! based on Aumont et al 2015 GMD 
         ! in PISCES, this reference depth is defined as the maximum of mixed layer depth and
         ! the euphotic zone. 
         zref = 100.0_rk;
         zmax = 5000.0_rk;
         if (depth .ge. zref) then
            zfact = MAX( 0., depth - zref ) / zmax
            wL = self%wL_min + ( self%wL_max - self%wL_min ) * zfact 
         else
            wL = self%wL_min
         endif

         ! Send rates of change to FABM, using the GET_VERTICAL_MOV module from FABM
         _SET_VERTICAL_MOVEMENT_(self%id_phyS,self%wPS)
         _SET_VERTICAL_MOVEMENT_(self%id_chlS,self%wPS)
         
         _SET_VERTICAL_MOVEMENT_(self%id_phyL,self%wPL)
         _SET_VERTICAL_MOVEMENT_(self%id_chlL,self%wPL)
         
         _SET_VERTICAL_MOVEMENT_(self%id_SDeN,self%wS)
         _SET_VERTICAL_MOVEMENT_(self%id_LDeN,wL)
         
         _SET_VERTICAL_MOVEMENT_(self%id_Opal,wL)
         _SET_VERTICAL_MOVEMENT_(self%id_Calc,wL)
         
         _SET_DIAGNOSTIC_(self%id_wL,-wL*secs_per_day) ! PAR at layer centre

      _VERTICAL_LOOP_END_
   
   END SUBROUTINE get_vertical_movement

END MODULE