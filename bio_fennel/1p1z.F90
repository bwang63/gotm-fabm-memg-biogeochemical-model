#include "fabm_driver.h"

!------------------------------------------------1p1z model-----------------------------------------------------
!
! This routine is a basic 1D version of the Fennel et at. (2006) ecosystem model. The model computes the sources and
! sinks for 2 nutrient pools (nitrate and ammonia), 1 phytoplankton pool, 1 zooplankton pool, and 2 detritus pools
! (large fast sinking detritus and small slow sinking detritus). Chlorophyll is also simulated as an independent
! variable. Photoacclimation is allowed based on Geider et al., (1997). 
!
! Developed by KF
! Adapted for GOTM-FABM by BW on May 2020 in Dalhousie university
!
! Reference:
!
!  Fennel, K., Wilkin, J., Levin, J., Moisan, J., O'Reilly, J., Haidvogel, D., 2006: 
!      Nitrogen cycling in the MidAtlantic Bight and implications for the North Atlantic nitrogen budget: Results
!      from a three-dimensional model. Global Biogeochemical Cycles 20, GB3007, doi:10.1029/2005GB002456.        
!
! Please go to this link for more information about GoTM-Fabm coding:
! https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model
!-----------------------------------------------------------------------------------------------------------------

MODULE memg_bio_fennel_1p1z

   USE fabm_types

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(type_base_model), PUBLIC :: type_memg_bio_fennel_1p1z
      ! Variable identifiers
      TYPE (type_state_variable_id) :: id_no3, id_nh4
      TYPE (type_state_variable_id) :: id_phy, id_chl
      TYPE (type_state_variable_id) :: id_zoo
      TYPE (type_state_variable_id) :: id_SDeN, id_LDeN
      
      TYPE (type_dependency_id) :: id_temp, id_par, id_oxy
      TYPE (type_horizontal_dependency_id) :: id_I_0

      TYPE (type_diagnostic_variable_id) :: id_newPP,id_regPP
      TYPE (type_diagnostic_variable_id) :: id_graz,id_Pmortal,id_coagP,id_Zmortal,id_Zmetabo,id_Zexcret
      TYPE (type_diagnostic_variable_id) :: id_assim,id_egest,id_coagD,id_remineS,id_remineL,id_nitrifi
      
      ! These diagnostic variables were set up here in order to check whether the model can correctly
      ! simulate source and sink terms. We don't need to setup all of these terms in a real application
      !TYPE (type_diagnostic_variable_id) :: id_bio1,id_bio2,id_bio3,id_bio4,id_bio5,id_bio6,id_bio7
      
      ! Parameters
      real(rk) :: attSW            ! PAR attenuation coef due to sea water [m-1]
      real(rk) :: attCHL           ! PAR attenuation coef due to chlorophyll [m-1]
      real(rk) :: PARfrac          ! fraction of photosynthetically available radiation [dimensionless]
      real(rk) :: I_thNH4          ! Radiation threshold for nitrification inhibition [Watts/m2]
      real(rk) :: D_p5NH4          ! Half-saturation radiation for nitrification inhibition [Watts/m2]
      real(rk) :: NitriR           ! Nitrification rate: oxidation of NH4 to NO3 [1/day]
      real(rk) :: K_NO3            ! Inverse half-saturation for phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: K_NH4            ! Inverse half-saturation for phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Vp0              ! Eppley temperature-limited growth parameter [nondimensional]
      real(rk) :: K_Phy            ! Zooplankton half-saturation constant for ingestion [1/day]
      real(rk) :: PhyCN            ! Phytoplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: PhyIS            ! Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]
      real(rk) :: PhyMR            ! Phytoplankton mortality rate [1/day]
      real(rk) :: Chl2C_m          ! Maximum chlorophyll to carbon ratio [mg Chl/mg C]
      real(rk) :: ZooAE_N          ! Zooplankton Nitrogen assimilation efficiency [nondimesnional]
      real(rk) :: ZooBM            ! Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooER            ! Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooGR            ! Zooplankton maximum growth rate [1/day]
      real(rk) :: ZooMR            ! Zooplankton mortality rate [1/day]
      real(rk) :: ZooCN            ! Zooplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: LDeRRN           ! Large detritus remineralization rate N-fraction [1/day]
      real(rk) :: SDeRRN           ! Small detritus remineralization rate N-fraction [1/day]
      real(rk) :: CoagR            ! Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN [1/day]
      real(rk) :: wP               ! Sinking velocity of phytoplankton [m/day]
      real(rk) :: wL               ! Sinking velocity of large detritus [m/day]
      real(rk) :: wS               ! Sinking velocity of small detritus [m/day]
      real(rk) :: PhyMin           ! Phytoplankton minimum threshold value [mmol N/m3]
      real(rk) :: ChlMin           ! Chlorophyll minimum threshold value [mmol N/m3]
      real(rk) :: ZooMin           ! Zooplankton minimum threshold value [mmol N/m3] 
         
      logical  :: oxy_model        ! whether to model oxygen
      
   CONTAINS
      PROCEDURE :: initialize
      ! Reference model procedures here.
      PROCEDURE :: do
      
   END TYPE

CONTAINS

   SUBROUTINE initialize(self,configunit)
      CLASS (type_memg_bio_fennel_1p1z), INTENT(INOUT), TARGET :: self
      INTEGER,                          INTENT(IN)              :: configunit
 
      ! Parameters
      real(rk) :: no3_initial      ! Initial Nitrate Concentration
      real(rk) :: nh4_initial      ! Initial Ammonia Concentration
      real(rk) :: phy_initial      ! Initial Phytoplankton Concentration
      real(rk) :: chl_initial      ! Initial Chlorophyll Concentration
      real(rk) :: zoo_initial      ! Initial Zooplankton Concentration
      real(rk) :: SDeN_initial     ! Initial small Detritus Concentration in nitrogen
      real(rk) :: LDeN_initial     ! Initial large Detritus Concentration in nitrogen
      
      real(rk) :: attSW            ! PAR attenuation coef due to sea water [m-1]
      real(rk) :: attCHL           ! PAR attenuation coef due to chlorophyll [m-1]
      real(rk) :: PARfrac          ! fraction of photosynthetically available radiation [dimensionless]
      real(rk) :: I_thNH4          ! Radiation threshold for nitrification inhibition [Watts/m2]
      real(rk) :: D_p5NH4          ! Half-saturation radiation for nitrification inhibition [Watts/m2]
      real(rk) :: NitriR           ! Nitrification rate: oxidation of NH4 to NO3 [1/day]
      real(rk) :: K_NO3            ! Inverse half-saturation for phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: K_NH4            ! Inverse half-saturation for phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Vp0              ! Eppley temperature-limited growth parameter [nondimensional]
      real(rk) :: K_Phy            ! Zooplankton half-saturation constant for ingestion [1/day]
      real(rk) :: PhyCN            ! Phytoplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: PhyIS            ! Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]
      real(rk) :: PhyMR            ! Phytoplankton mortality rate [1/day]
      real(rk) :: Chl2C_m          ! Maximum chlorophyll to carbon ratio [mg Chl/mg C]
      real(rk) :: ZooAE_N          ! Zooplankton Nitrogen assimilation efficiency [nondimesnional]
      real(rk) :: ZooBM            ! Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooER            ! Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooGR            ! Zooplankton maximum growth rate [1/day]
      real(rk) :: ZooMR            ! Zooplankton mortality rate [1/day]
      real(rk) :: ZooCN            ! Zooplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: LDeRRN           ! Large detritus remineralization rate N-fraction [1/day]
      real(rk) :: SDeRRN           ! Small detritus remineralization rate N-fraction [1/day]
      real(rk) :: CoagR            ! Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN [1/day]
      real(rk) :: wP               ! Sinking velocity of phytoplankton [m/day]
      real(rk) :: wL               ! Sinking velocity of large detritus [m/day]
      real(rk) :: wS               ! Sinking velocity of small detritus [m/day]
      real(rk) :: PhyMin           ! Phytoplankton minimum threshold value [mmol N/m3]
      real(rk) :: ChlMin           ! Chlorophyll minimum threshold value [mmol N/m3]
      real(rk) :: ZooMin           ! Zooplankton minimum threshold value [mmol N/m3] 
      
      logical  :: oxy_model        ! whether to model oxygen
      
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      namelist /memg_bio_fennel_1p1z/  no3_initial,nh4_initial, &
                                       chl_initial,phy_initial,zoo_initial,  &
                                       SDeN_initial,LDeN_initial,   &
                                       attSW, attCHL, PARfrac,   &
                                       I_thNH4, D_p5NH4, NitriR, K_NO3, K_NH4,   &
                                       Vp0, K_Phy, PhyCN, PhyIS, PhyMR,   & 
                                       Chl2C_m,   &
                                       ZooAE_N, ZooBM, ZooER, ZooGR, ZooMR, ZooCN,   &
                                       LDeRRN, SDeRRN, CoagR,   &
                                       wP, wL, wS,   &
                                       PhyMin, ChlMin, ZooMin
!
! set default values for initial condition
!
      no3_initial = 1.0_rk      
      nh4_initial = 0.1_rk     
      phy_initial = 0.1_rk     
      chl_initial = 0.1_rk     
      zoo_initial = 0.01_rk    
      SDeN_initial = 0.01_rk   
      LDeN_initial = 0.01_rk
!
! set default values for parameter
!
      attSW = 0.04_rk            
      attCHL = 0.02486_rk           
      PARfrac = 0.43_rk          
      I_thNH4 = 0.0095_rk         
      D_p5NH4 = 0.1_rk          
      NitriR = 0.2_rk          
      K_NO3 = 2.0_rk         
      K_NH4 = 2.0_rk           
      Vp0 = 0.69_rk           
      K_Phy = 0.5_rk   
      PhyCN = 6.625_rk         
      PhyIS = 0.025_rk          
      PhyMR = 0.072_rk            
      Chl2C_m = 0.0535_rk        
      ZooAE_N = 0.75_rk          
      ZooBM = 0.1_rk            
      ZooER = 0.1_rk           
      ZooGR = 0.75_rk         
      ZooMR = 0.025_rk
      ZooCN = 6.625_rk
      LDeRRN = 0.01_rk        
      SDeRRN = 0.03_rk         
      CoagR = 0.005_rk           
      wP = 0.1_rk               
      wL = 1.0_rk              
      wS = 0.1_rk              
      PhyMin = 0.000001_rk           
      ChlMin = 0.000001_rk            
      ZooMin = 0.000001_rk             

      oxy_model = .false.
      
      ! Read the namelist
      if (configunit>0) read(configunit,nml=memg_bio_fennel_1p1z,err=99,end=100)

      ! Store parameter values into parameters
      ! Note: all rates are provided in values per day
      call self%get_parameter(self%attSW,'attSW','m-1','PAR attenuation coef due to sea water',default=attSW)
      call self%get_parameter(self%attCHL,'attCHL','(mg Chl m2)-1','PAR attenuation coef due to chlorophyll',default=attCHL)
      call self%get_parameter(self%PARfrac,'PARfrac','-','fraction of photosynthetically available radiation',default=PARfrac)
      call self%get_parameter(self%I_thNH4,'I_thNH4','Watts.m-2','Radiation threshold for nitrification inhibition',default=I_thNH4)
      call self%get_parameter(self%D_p5NH4,'D_p5NH4','Watts.m-2','Half-saturation radiation for nitrification inhibition',default=D_p5NH4)
      call self%get_parameter(self%NitriR,'NitriR','d-1','Nitrification rate: oxidation of NH4 to NO3',default=NitriR)
      call self%get_parameter(self%K_NO3,'K_NO3','(mmol N m-3)-1','Inverse half-saturation for phytoplankton NO3 uptake',default=K_NO3)
      call self%get_parameter(self%K_NH4,'K_NH4','(mmol N m-3)-1','Inverse half-saturation for phytoplankton NH4 uptake',default=K_NH4)
      call self%get_parameter(self%Vp0,'Vp0','d-1','Eppley temperature-limited growth parameter',default=Vp0)
      call self%get_parameter(self%K_Phy,'K_Phy','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion',default=K_Phy)
      call self%get_parameter(self%PhyCN,'PhyCN','mmol C/mmol N','Phytoplankton Carbon:Nitrogen ratio',default=PhyCN)
      call self%get_parameter(self%PhyIS,'PhyIS','(Watts m-2 day)-1','Phytoplankton, initial slope of P-I curve',default=PhyIS)
      call self%get_parameter(self%PhyMR,'PhyMR','d-1','Phytoplankton mortality rate',default=PhyMR)
      call self%get_parameter(self%Chl2C_m,'Chl2C_m','mg Chl/mg C','Maximum chlorophyll to carbon ratio',default=Chl2C_m)
      call self%get_parameter(self%ZooAE_N,'ZooAE_N','-','Zooplankton Nitrogen assimilation efficiency',default=ZooAE_N)
      call self%get_parameter(self%ZooBM,'ZooBM','d-1','Zooplankton Basal metabolism',default=ZooBM)
      call self%get_parameter(self%ZooER,'ZooER','d-1','Zooplankton specific excretion rate',default=ZooER)
      call self%get_parameter(self%ZooGR,'ZooGR','d-1','Zooplankton maximum growth rate',default=ZooGR)
      call self%get_parameter(self%ZooMR,'ZooMR','d-1','Zooplankton mortality rate',default=ZooMR)
      call self%get_parameter(self%ZooCN,'ZooCN','mmol C/mmol N','Zooplankton Carbon:Nitrogen ratio',default=ZooCN)
      call self%get_parameter(self%LDeRRN,'LDeRRN','d-1','Large detritus remineralization rate N-fraction',default=LDeRRN)
      call self%get_parameter(self%SDeRRN,'SDeRRN','d-1','Small detritus remineralization rate N-fraction',default=SDeRRN)
      call self%get_parameter(self%CoagR,'CoagR','d-1','Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN',default=CoagR)
      call self%get_parameter(self%wP,'wP','m d-1','Sinking velocity of phytoplankton (<0 for sinking)',default=wP,scale_factor=-1.0_rk*d_per_s)
      call self%get_parameter(self%wL,'wL','m d-1','Sinking velocity of large detritus (<0 for sinking)',default=wL,scale_factor=-1.0_rk*d_per_s)
      call self%get_parameter(self%wS,'wS','m d-1','Sinking velocity of samll detritus (<0 for sinking)',default=wS,scale_factor=-1.0_rk*d_per_s)
      call self%get_parameter(self%PhyMin,'PhyMin','mmol N m-3','Phytoplankton minimum threshold value',default=PhyMin)
      call self%get_parameter(self%ChlMin,'ChlMin','mmol N m-3','Chlorophyll minimum threshold value',default=ChlMin)
      call self%get_parameter(self%ZooMin,'ZooMin','mmol N m-3','Zooplankton minimum threshold value',default=ZooMin)
      
      call self%get_parameter(self%oxy_model,'oxy_model','','Switch to modelling oxygen',default=oxy_model)
 
      ! Register state variables
      call self%register_state_variable(self%id_no3,'no3','mmol N m-3','nitrate',no3_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='nitrate',units='mmol N m-3'))
      call self%register_state_variable(self%id_nh4,'nh4','mmol N m-3','ammonia',nh4_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='ammonia',units='mmol N m-3'))
      call self%register_state_variable(self%id_phy,'phy','mmol N m-3','phytoplankton',phy_initial,minimum=0.0_rk,vertical_movement=self%wP,   &
         &    standard_variable=type_bulk_standard_variable(name='phytoplankton',units='mmol N m-3'))
      call self%register_state_variable(self%id_chl,'chl','mg m-3','chlorophyll',chl_initial,minimum=0.0_rk,vertical_movement=self%wP,     &
         &    standard_variable=type_bulk_standard_variable(name='chlorophyll',units='mg m-3'))
      call self%register_state_variable(self%id_zoo,'zoo','mmol N m-3','zooplankton',zoo_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='zooplankton',units='mmol N m-3'))
      call self%register_state_variable(self%id_LDeN,'LDeN','mmol N m-3','LdetritusN',LDeN_initial,minimum=0.0_rk,vertical_movement=self%wL,   &
         &    standard_variable=type_bulk_standard_variable(name='LdetritusN',units='mmol N m-3'))
      call self%register_state_variable(self%id_SDeN,'SDeN','mmol N m-3','SdetritusN',SDeN_initial,minimum=0.0_rk,vertical_movement=self%wS,   &
         &    standard_variable=type_bulk_standard_variable(name='SdetritusN',units='mmol N m-3'))
      

      ! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_shortwave_flux)
      
      if (self%oxy_model) then
         call self%register_dependency(self%id_oxy, standard_variable=type_bulk_standard_variable(name='oxygen',units='mmol O m-3')) ! calculated by the oxygen model
      end if
      
      call self%register_dependency(self%id_par, standard_variable=type_bulk_standard_variable(name='photosynthetic_active_radiation_in_water_column',units='W m-2')) ! calculated by the light model
      
      ! Register the contribution of all state variables to total nitrogen
         
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_newPP,'NewPP','mmol N m-3 day-1','new primary production by NO3',   &
         &    standard_variable=type_bulk_standard_variable(name='NewPP',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_regPP,'RegPP','mmol N m-3 day-1','regenerate primary production by NH4',   &
         &    standard_variable=type_bulk_standard_variable(name='RegPP',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_graz,'Graz','mmol N m-3 day-1','grazing by zooplankton',   &
         &    standard_variable=type_bulk_standard_variable(name='Graz',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Pmortal,'Pmortal','mmol N m-3 day-1','phytoplankton mortality',   &
         &    standard_variable=type_bulk_standard_variable(name='Pmortal',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_coagP,'CoagP','mmol N m-3 day-1','coaggregation of phytoplankton',   &
         &    standard_variable=type_bulk_standard_variable(name='CoagP',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Zmortal,'Zmortal','mmol N m-3 day-1','zooplankton mortality',   &
         &    standard_variable=type_bulk_standard_variable(name='Zmortal',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Zmetabo,'Zmetabo','mmol N m-3 day-1','zooplankton metabolism',   &
         &    standard_variable=type_bulk_standard_variable(name='Zmetabo',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Zexcret,'Zexcret','mmol N m-3 day-1','zooplankton excretion',   &
         &    standard_variable=type_bulk_standard_variable(name='Zexcret',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_assim,'Assim','mmol N m-3 day-1','zooplankton assimilation',   &
         &    standard_variable=type_bulk_standard_variable(name='Assim',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_egest,'Egest','mmol N m-3 day-1','zooplankton egestion',   &
         &    standard_variable=type_bulk_standard_variable(name='Egest',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_coagD,'CoagD','mmol N m-3 day-1','coaggregation of small detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='CoagD',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_remineS,'RemineS','mmol N m-3 day-1','remineralization of small detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='RemineS',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_remineL,'RemineL','mmol N m-3 day-1','remineralization of large detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='RemineL',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_nitrifi,'Nitrifi','mmol N m-3 day-1','nitrification from NH4 into NO3',   &
         &    standard_variable=type_bulk_standard_variable(name='Nitrifi',units='mmol N m-3 day-1'))
      
      
      ! These diagnostic variables were set up here in order to check whether the model can correctly
      ! simulate source and sink terms. We don't need to setup all of these terms in a real application
      !call self%register_diagnostic_variable(self%id_bio1,'bio1',' ','')
      !call self%register_diagnostic_variable(self%id_bio2,'bio2',' ','')
      !call self%register_diagnostic_variable(self%id_bio3,'bio3',' ','')
      !call self%register_diagnostic_variable(self%id_bio4,'bio4',' ','')
      !call self%register_diagnostic_variable(self%id_bio5,'bio5',' ','')
      !call self%register_diagnostic_variable(self%id_bio6,'bio6',' ','')
      !call self%register_diagnostic_variable(self%id_bio7,'bio7',' ','')
     
      
   return

99 call self%fatal_error('memg_bio_fennel_1p1z','Error reading namelist memg_bio_fennel_1p1z.')

100 call self%fatal_error('memg_bio_fennel_1p1z','Namelist memg_bio_fennel_1p1z was not found.')

   ! Register model parameters and variables here.
   END SUBROUTINE initialize

   SUBROUTINE do(self,_ARGUMENTS_DO_)
      class (type_memg_bio_fennel_1p1z),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
   
      real(rk) :: no3,nh4,phy,zoo,LDeN,SDeN,chl
      real(rk) :: temp,par,oxy,I_0 
      
      ! Local variables
      real(rk) :: Chl2C,Vp,Epp,t_PPmax
      real(rk) :: inhNH4,L_NH4,L_NO3,LTOT
      real(rk) :: cff,cff1,cff2,cff3,cff4,cff5
      real(rk) :: fac1,fac2,fac3
      
      real(rk) :: N_NewProd,N_RegProd,Chl_Prod
      real(rk) :: N_Nitrifi
      real(rk) :: N_Graz,Chl_Graz
      real(rk) :: N_Assim,N_Egest,N_Pmortal,Chl_Pmortal
      real(rk) :: N_Zmortal,N_Zexcret,N_Zmetabo
      real(rk) :: N_CoagP,N_CoagD,Chl_Coag
      real(rk) :: N_RemineS,N_RemineL
      real(rk) :: d_no3,d_nh4,d_phy,d_zoo,d_LDeN,d_SDeN,d_chl
      
      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk), parameter :: eps = 1.0e-20_rk
      
      ! parameter associated with modeling oxygen
      real(rk), parameter :: rOxNO3= 8.625_rk       ! 138/16
      real(rk), parameter :: rOxNH4= 6.625_rk       ! 106/16
      
      _LOOP_BEGIN_
      
         ! Obtain concentration of biological variables.
         _GET_(self%id_no3,no3)
         _GET_(self%id_nh4,nh4)
         _GET_(self%id_phy,phy)
         _GET_(self%id_chl,chl)
         _GET_(self%id_zoo,zoo)
         _GET_(self%id_LDeN,LDeN)
         _GET_(self%id_SDeN,SDeN)
           
         ! Obtain environmental dependencies
         _GET_(self%id_temp,temp)
         
         if (self%oxy_model) then
            _GET_(self%id_oxy,oxy)
         end if
         
         _GET_(self%id_par,par)
         _GET_HORIZONTAL_(self%id_I_0,I_0)
         
         ! These diagnostic variables were set up here in order to check whether the model can correctly
         ! simulate source and sink terms. We don't need to setup these terms in a real application
         ! In my case, I firstly used model outputs of these 7 biological variables to do an offline
         ! calculation of all biological fluxes by a well-calibrated matlab version model. However, the
         ! difference between the GoTM fluxes and the matlab fluxes, e.g. newPP can be as large as 2%. I
         ! realized this is because GoTM used concentrations before adding source/sink fluxes while matlab
         ! used concentrations after that. Therefore, I also set those concentrations before adding source/
         ! sink fluxes as output. By using these concentrations, the difference in fluxes between GoTM and
         ! matlab is within 2e-7. That suggested my GoTM version of bio_Fennel model could calculate all source
         ! and sink terms correctly.
         !_SET_DIAGNOSTIC_(self%id_bio1,no3)
         !_SET_DIAGNOSTIC_(self%id_bio2,nh4)
         !_SET_DIAGNOSTIC_(self%id_bio3,phy)
         !_SET_DIAGNOSTIC_(self%id_bio4,chl)
         !_SET_DIAGNOSTIC_(self%id_bio5,zoo)
         !_SET_DIAGNOSTIC_(self%id_bio6,LDeN)
         !_SET_DIAGNOSTIC_(self%id_bio7,SDeN)
         
         
         IF (I_0 .gt. 0.0_rk) THEN
!-----------------------------------------------------------------------
! Compute ration of Chlorophyll-a to phytoplankton, [mg Chla / (mg C)]
!-----------------------------------------------------------------------
            cff = self%PhyCN*12.0_rk
            Chl2C = MIN(chl/(phy*cff+eps), self%Chl2C_m)

!-----------------------------------------------------------------------
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!-----------------------------------------------------------------------
            Vp = self%Vp0*0.59_rk*(1.066_rk**temp)
            fac1 = par*self%PhyIS
            Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
            t_PPmax=Epp*fac1
            
!-----------------------------------------------------------------------
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!-----------------------------------------------------------------------
            cff1=nh4*self%K_NH4
            cff2=no3*self%K_NO3   
            inhNH4=1.0_rk/(1.0_rk+cff1)
            L_NH4=cff1/(1.0_rk+cff1)
            L_NO3=cff2*inhNH4/(1.0_rk+cff2)
            LTOT=L_NO3+L_NH4
         
            !  Nitrate and ammonium uptake by Phytoplankton.
            fac1=d_per_s*t_PPmax
            cff4=fac1*self%K_NO3*inhNH4/(1.0_rk+cff2)*phy
            cff5=fac1*self%K_NH4/(1.0_rk+cff1)*phy
            
            ! In ROMS, bio_fennel deals with some sink terms in implicit way
            ! The gotm version of bio_fennel model trys to deal with the
            ! source and sink terms the same as in ROMS. no3 and nh4 values
            ! will be changed in order to calculate some implicit sink terms
            ! but not be saved into the output files.
            no3=no3/(1.0_rk+cff4)
            nh4=nh4/(1.0_rk+cff5)

            N_NewProd=no3*cff4
            N_RegProd=nh4*cff5
            
            phy = phy+N_NewProd+N_RegProd
            
            ! increse in chlorophyll because of phytoplankton growth 
            Chl_Prod = (d_per_s*t_PPmax*t_PPmax*LTOT*LTOT*self%Chl2C_m*chl)/  &
     &                 (self%PhyIS*MAX(Chl2C,eps)*par+eps)
            
            chl=chl+Chl_Prod
            
            ! increase in oxygen because of phytoplankton growth
            if (self%oxy_model) then
               oxy=oxy+N_NewProd*rOxNO3+N_RegProd*rOxNH4
            end if

!-----------------------------------------------------------------------
!  The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
!  only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
! 
!          NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!          NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
! 
!  Note that the entire process has a total loss of two moles of O2 per
!  mole of NH4. If we were to resolve NO2 profiles, this is where we
!  would change the code to split out the differential effects of the
!  two different bacteria types. If OXYGEN is defined, nitrification is
!  inhibited at low oxygen concentrations using a Michaelis-Menten term.
!-----------------------------------------------------------------------
            if (self%oxy_model) then
               fac2=MAX(oxy,0.0_rk)     ! O2 max
               fac3=MAX(fac2/(3.0_rk+fac2),0.0_rk) ! MM for O2 dependence
               fac1=d_per_s*self%NitriR*fac3
            else
               fac1=d_per_s*self%NitriR
            end if
            cff1=(par-self%I_thNH4)/(self%D_p5NH4+par-2.0_rk*self%I_thNH4)      
            cff2=1.0_rk-MAX(0.0_rk,cff1)
            cff3=fac1*cff2
            nh4=nh4/(1.0_rk+cff3)
            N_Nitrifi=nh4*cff3
            no3=no3+N_Nitrifi
            
            if (self%oxy_model) then
               oxy=oxy-2.0_rk*N_Nitrifi
            end if
            
         ELSE ! If I0=0, no phytoplankton growth and nitrification ocurs at maximum rate (NitriR)
            N_NewProd = 0.0_rk
            N_RegProd = 0.0_rk
            Chl_Prod = 0.0_rk
            
            cff3=d_per_s*self%NitriR
            nh4=nh4/(1.0_rk+cff3)
            N_Nitrifi=nh4*cff3
            no3=no3+N_Nitrifi
            
            if (self%oxy_model) then
               oxy=oxy-2.0_rk*N_Nitrifi
            end if
         END IF

!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
         fac1=d_per_s*self%ZooGR
         cff2=d_per_s*self%PhyMR
!
! Phytoplankton grazing by zooplankton.
!
         cff1=fac1*zoo*phy/(self%K_Phy+phy*phy)
         cff3=1.0_rk/(1.0_rk+cff1)
         phy=phy*cff3
         chl=chl*cff3
         
         N_Graz=cff1*phy
         Chl_Graz=cff1*chl
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
         N_Assim=cff1*phy*self%ZooAE_N
         N_Egest=cff1*phy*(1.0_rk-self%ZooAE_N)
         
         zoo = zoo+N_Assim
         SDeN = SDeN+N_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
         N_Pmortal=cff2*MAX(phy-self%PhyMin,0.0_rk)
         phy=phy-N_Pmortal
         
         Chl_Pmortal=cff2*MAX(chl-self%ChlMin,0.0_rk)
         chl=chl-Chl_Pmortal
         
         SDeN=SDeN+N_Pmortal
         
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
         cff1=d_per_s*self%ZooBM
         fac2=d_per_s*self%ZooMR
         fac3=d_per_s*self%ZooER

         fac1=fac3*phy*phy/(self%K_Phy+phy*phy)                               
         cff2=fac2*zoo
         cff3=fac1*self%ZooAE_N
         zoo=zoo/(1.0_rk+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
         N_Zmortal=cff2*zoo
         N_Zexcret=cff3*zoo
         
         nh4=nh4+N_Zexcret
         SDeN=SDeN+N_Zmortal

!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
         N_Zmetabo=cff1*MAX(zoo-self%ZooMin,0.0_rk)
         zoo=zoo-N_Zmetabo
         nh4=nh4+N_Zmetabo
         if (self%oxy_model) then
            oxy=oxy-rOxNH4*(N_Zmetabo+N_Zexcret)
         end if
!         
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
         fac1=d_per_s*self%CoagR
         cff1=fac1*(SDeN+phy)
         cff2=1.0_rk/(1.0_rk+cff1)
         phy=phy*cff2
         chl=chl*cff2
         SDeN=SDeN*cff2
         
         N_CoagP=phy*cff1
         Chl_Coag=Chl*cff1
         N_CoagD=SDeN*cff1
         
         LDeN=LDeN+N_CoagP+N_CoagD
         
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
         if (self%oxy_model) then
            fac1=MAX(oxy-6.0_rk,0.0_rk) ! O2 off max
            fac2=MAX(fac1/(3.0_rk+fac1),0.0_rk) ! MM for O2 dependence
            
            cff1=d_per_s*self%SDeRRN*fac2
            cff2=1.0_rk/(1.0_rk+cff1)
            cff3=d_per_s*self%LDeRRN*fac2
            cff4=1.0_rk/(1.0_rk+cff3)
            
            SDeN=SDeN*cff2
            LDeN=LDeN*cff4
            
            N_RemineS=SDeN*cff1
            N_RemineL=LDeN*cff3
            
            nh4=nh4+N_RemineS+N_RemineL
            oxy=oxy-(N_RemineS+N_RemineL)*rOxNH4
         else
            cff1=d_per_s*self%SDeRRN
            cff2=1.0_rk/(1.0_rk+cff1)
            cff3=d_per_s*self%LDeRRN
            cff4=1.0_rk/(1.0_rk+cff3)
         
            SDeN=SDeN*cff2
            LDeN=LDeN*cff4
         
            N_RemineS=SDeN*cff1
            N_RemineL=LDeN*cff3
         
            nh4=nh4+N_RemineS+N_RemineL
         end if

!
!-----------------------------------------------------------------------
!  Sum all contributor to change rate:
!-----------------------------------------------------------------------
!
         d_no3 = -N_NewProd+N_Nitrifi
         d_nh4 = -N_RegProd-N_Nitrifi+N_Zmetabo+N_Zexcret+N_RemineS+N_RemineL
         d_phy = N_NewProd+N_RegProd-N_Graz-N_Pmortal-N_CoagP
         d_chl = Chl_Prod-Chl_Graz-Chl_Pmortal-Chl_Coag
         d_zoo = N_Assim-N_Zmortal-N_Zexcret-N_Zmetabo
         d_SDeN = N_Egest+N_Pmortal+N_Zmortal-N_CoagD-N_RemineS
         d_LDeN = N_CoagP+N_CoagD-N_RemineL   
         
         
         
         ! Send rates of change to FABM.
         _SET_ODE_(self%id_no3,d_no3)
         _SET_ODE_(self%id_nh4,d_nh4)
         _SET_ODE_(self%id_phy,d_phy)
         _SET_ODE_(self%id_chl,d_chl)
         _SET_ODE_(self%id_zoo,d_zoo)
         _SET_ODE_(self%id_LDeN,d_LDeN)
         _SET_ODE_(self%id_SDeN,d_SDeN)
         
         ! Provide diagnostic variables to FABM.
         _SET_DIAGNOSTIC_(self%id_newPP,secs_pr_day*N_NewProd)
         _SET_DIAGNOSTIC_(self%id_regPP,secs_pr_day*N_RegProd)
         
         ! These diagnostic variables were set up here in order to check whether the model can correctly
         ! simulate source and sink terms. We don't need to setup all of these terms in a real application
         ! But please be careful that some fluxes are useful in the oxygen model
         _SET_DIAGNOSTIC_(self%id_graz,secs_pr_day*N_Graz)
         _SET_DIAGNOSTIC_(self%id_Pmortal,secs_pr_day*N_Pmortal)
         _SET_DIAGNOSTIC_(self%id_coagP,secs_pr_day*N_CoagP)
         _SET_DIAGNOSTIC_(self%id_Zmortal,secs_pr_day*N_Zmortal)
         _SET_DIAGNOSTIC_(self%id_Zmetabo,secs_pr_day*N_Zmetabo)
         _SET_DIAGNOSTIC_(self%id_Zexcret,secs_pr_day*N_Zexcret)
         _SET_DIAGNOSTIC_(self%id_assim,secs_pr_day*N_Assim)
         _SET_DIAGNOSTIC_(self%id_egest,secs_pr_day*N_Egest)
         _SET_DIAGNOSTIC_(self%id_coagD,secs_pr_day*N_CoagD)
         _SET_DIAGNOSTIC_(self%id_remineS,secs_pr_day*N_RemineS)
         _SET_DIAGNOSTIC_(self%id_remineL,secs_pr_day*N_RemineL)
         _SET_DIAGNOSTIC_(self%id_nitrifi,secs_pr_day*N_Nitrifi)
         
      _LOOP_END_

   
   END SUBROUTINE do

END MODULE
