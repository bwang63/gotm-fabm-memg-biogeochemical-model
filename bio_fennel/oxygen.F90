#include "fabm_driver.h"

!------------------------------------------------oxygen model -----------------------------------------------------
!
! This routine is to simulate the source and sink fluxes of oxygen. 
!
! The source and sink terms of oxygen including:
!                                   (1) new primary production by NO3 (NewPP) 
!                                   (2) regenerate primary production by NH4 (RegPP)  
!                                   (3) nitrification (N_Nitrifi)
!                                   (4) zooplankton metabolism (Zmetabo)
!                                   (5) zooplankton excretion (Zexcret)
!                                   (6) remineralization of small detritus (RemineS)
!                                   (7) remineralization of large detritus (RemineL)
! were calculated by the nitrogen cycle model (e.g. 1p1z, or 2p2z) and were added into the oxygen pool by this oxygen
! model.
!
! The air-sea fluxes of oxygen (mmol O-2 day-2) were calcuated in this model by the subroutine do_surface
!  
!           air-sea flux = coef*(D0sat - DO)
!           coef = 0.31*wind_speed*wind_speed*sqrt(660/Schmidt_number)
!
! The bottm-water fluxes of oxygen (mmol O-2 day-2) were not included yet.
!
! Developed by KF
! Adapted for GOTM-FABM by BW on May 2020 in Dalhousie university    
!
! Please go to this link for more information about GoTM-Fabm coding:
! https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model
!-----------------------------------------------------------------------------------------------------------------

MODULE memg_bio_fennel_oxygen

   USE fabm_types

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(type_base_model), PUBLIC :: type_memg_bio_fennel_oxygen
      ! Variable identifiers
      TYPE (type_state_variable_id) :: id_oxy
      
      TYPE (type_dependency_id) :: id_temp,id_salt,id_cell
      TYPE (type_dependency_id) :: id_newPP,id_regPP,id_nitrifi,id_Zmetabo,id_Zexcret,id_remineS,id_remineL
      
      TYPE (type_horizontal_dependency_id) :: id_wind

      TYPE (type_surface_diagnostic_variable_id) :: id_oxyflx,id_oxysat

      ! Parameters
      
   CONTAINS
      PROCEDURE :: initialize
      ! Reference model procedures here.
      PROCEDURE :: do ! add source and sink fluxes into oxygen pool
      PROCEDURE :: do_surface ! calculate air-sea fluxes
      
   END TYPE

CONTAINS

   SUBROUTINE initialize(self,configunit)
      CLASS (type_memg_bio_fennel_oxygen), INTENT(INOUT), TARGET :: self
      INTEGER,                          INTENT(IN)              :: configunit
 
      ! Parameters
      real(rk) :: oxygen_initial      ! Initial oxygen Concentration

      namelist /memg_bio_fennel_oxygen/ oxygen_initial
!
! set default values for model parameters
!
      oxygen_initial = 300.0_rk
      
      ! Read the namelist
      if (configunit>0) read(configunit,nml=memg_bio_fennel_oxygen,err=99,end=100)

      ! Store parameter values into parameters
 
      ! Register state variables
      call self%register_state_variable(self%id_oxy,'oxy','mmol O m-3','dissolved_oxygen',oxygen_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='oxygen',units='mmol O m-3'))

      ! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)
      call self%register_dependency(self%id_cell, standard_variables%cell_thickness)
      
      call self%register_dependency(self%id_newPP, standard_variable=type_bulk_standard_variable(name='NewPP',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_regPP, standard_variable=type_bulk_standard_variable(name='RegPP',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_nitrifi, standard_variable=type_bulk_standard_variable(name='Nitrifi',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_Zmetabo, standard_variable=type_bulk_standard_variable(name='Zmetabo',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_Zexcret, standard_variable=type_bulk_standard_variable(name='Zexcret',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_remineS, standard_variable=type_bulk_standard_variable(name='RemineS',units='mmol N m-3 day-1'))
      call self%register_dependency(self%id_remineL, standard_variable=type_bulk_standard_variable(name='RemineL',units='mmol N m-3 day-1'))
         
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_oxyflx,'airseaflx','mmol O m-2 day-1','air sea flux of oxygen')
      call self%register_diagnostic_variable(self%id_oxysat,'oxysat','mmol O m-3 day-1','saturated oxygen concentration')

   return

99 call self%fatal_error('memg_bio_fennel_airseflx','Error reading namelist memg_bio_fennel_airseflx.')

100 call self%fatal_error('memg_bio_fennel_airseflx','Namelist memg_bio_fennel_airseflx was not found.')

   ! Register model parameters and variables here.
   END SUBROUTINE initialize

   SUBROUTINE do(self,_ARGUMENTS_DO_)
      class (type_memg_bio_fennel_oxygen),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
       
      real(rk) :: oxy 
      
      ! Local variables
      real(rk) :: N_NewProd,N_RegProd
      real(rk) :: N_Nitrifi
      real(rk) :: N_Zexcret,N_Zmetabo
      real(rk) :: N_RemineS,N_RemineL
      real(rk) :: d_oxy
      
      real(rk), parameter :: secs_pr_day = 86400.0_rk
      
      ! parameter associated with modeling oxygen
      real(rk), parameter :: rOxNO3= 8.625_rk       ! 138/16
      real(rk), parameter :: rOxNH4= 6.625_rk       ! 106/16
      
      _LOOP_BEGIN_
      
         ! Obtain concentration of oxygen and its source and sink fluxes
         _GET_(self%id_oxy,oxy)
         _GET_(self%id_newPP,N_NewProd)
         _GET_(self%id_regPP,N_RegProd)
         _GET_(self%id_nitrifi,N_Nitrifi)
         _GET_(self%id_Zmetabo,N_Zmetabo)
         _GET_(self%id_Zexcret,N_Zexcret)
         _GET_(self%id_remineS,N_RemineS)
         _GET_(self%id_remineL,N_RemineL)
         
         d_oxy=N_NewProd*rOxNO3+N_RegProd*rOxNH4   &
      &        -2.0_rk*N_Nitrifi &
      &        -rOxNH4*(N_Zmetabo+N_Zexcret) &
      &        -rOxNH4*(N_RemineS+N_RemineL)
         d_oxy = d_oxy/secs_pr_day
         
         _SET_ODE_(self%id_oxy,d_oxy)
      

      _LOOP_END_
   END SUBROUTINE do
   
   SUBROUTINE do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_memg_bio_fennel_oxygen),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
   
      real(rk) :: oxy,temp,salt,wind,zcel
      
      ! Local variables
      real(rk) :: cff2,cff3
      real(rk) :: u10squ,SchmidtN_Ox,TS,AA
      real(rk) :: O2_Flux,O2satu
      
      
      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk), parameter :: eps = 1.0e-20_rk
      
      ! parameter associated with modeling oxygen
      ! Schmidt number coefficients using the formulation of Wanninkhof (1992)
      real(rk), parameter :: A_O2 = 1953.4_rk
      real(rk), parameter :: B_O2 = 128.0_rk
      real(rk), parameter :: C_O2 = 3.9918_rk
      real(rk), parameter :: D_O2 = 0.050091_rk
      real(rk), parameter :: E_O2 = 0.0_rk
      ! Oxygen saturation coefficients
      real(rk), parameter :: OA0 = 2.00907_rk       
      real(rk), parameter :: OA1 = 3.22014_rk       
      real(rk), parameter :: OA2 = 4.05010_rk       
      real(rk), parameter :: OA3 = 4.94457_rk
      real(rk), parameter :: OA4 =-0.256847_rk
      real(rk), parameter :: OA5 = 3.88767_rk
      real(rk), parameter :: OB0 =-0.00624523_rk
      real(rk), parameter :: OB1 =-0.00737614_rk
      real(rk), parameter :: OB2 =-0.0103410_rk
      real(rk), parameter :: OB3 =-0.00817083_rk
      real(rk), parameter :: OC0 =-0.000000488682_rk
      real(rk), parameter :: l2mol = 1000.0_rk/22.3916_rk      ! liter to mol
      
      
      _HORIZONTAL_LOOP_BEGIN_
      
         ! Obtain concentration of biological variables.
         _GET_(self%id_oxy,oxy)

         ! Obtain environmental dependencies
         _GET_(self%id_temp,temp)
         _GET_(self%id_salt,salt)
         _GET_(self%id_cell,zcel)
         _GET_HORIZONTAL_(self%id_wind,wind)
         
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!

         
! ROMS also provided another coefficient values to calculate air-sea fluxes: cff2=dtdays*0.251_r8*24.0_r8/100.0_r8
         cff2=d_per_s*0.31_rk*24.0_rk/100.0_rk

!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
         u10squ=wind*wind
         SchmidtN_Ox=A_O2-temp*(B_O2-temp*(C_O2-temp*(D_O2-temp*E_O2))) ! Updated by AL
         cff3=cff2*u10squ*SQRT(660.0_rk/SchmidtN_Ox)
         
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
         TS=LOG((298.15_rk-temp)/(273.15_rk+temp))
         AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
     &          salt*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &          OC0*salt*salt
!
!  Convert from ml/l to mmol/m3.
!
         O2satu=l2mol*EXP(AA)

!
!  O2 gas exchange.
!
         O2_Flux=cff3*(O2satu-oxy)    
         
!
!-----------------------------------------------------------------------
!  Update the State variables and store the diagnostics
!-----------------------------------------------------------------------
!         
         _SET_SURFACE_EXCHANGE_(self%id_oxy,O2_Flux/zcel)
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_oxyflx,secs_pr_day*O2_Flux)
         _SET_SURFACE_DIAGNOSTIC_(self%id_oxysat,O2satu)

      _HORIZONTAL_LOOP_END_

   
   END SUBROUTINE do_surface

END MODULE
