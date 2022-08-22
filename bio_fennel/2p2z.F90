#include "fabm_driver.h"

!--------------------------------------------------- 2p2z -------------------------------------------------------
!
! This model is an update version of the biogeochemical model described in Laurent et al., (2021) to include 2 
! different parameterization schemes of large POC. 
!
! (1) Ballast scheme:
! A part of large POC is protected from remineralization by minerals, including opal and CaCO3.
! 
! Following nemuro model (Kishi et al 2007) and BEC model (Lima et al 2014), the opal is produced by mortality,
! grazing, and aggregation of large phytoplankton
! The Si:N ratio of large phytoplankton is assumed constant.
!
! The calcite carbonate (CaCO3) production is modelled as a function of large detritus production following MEDUSA
! model (Yool et al 2011). The ratio between calcite and carbon production is dependent on the latitude.
!
! (2) WLin scheme
! The sinking velocity of large POC is assumed to increase linearly with the depth
!
! Adapted for GOTM-FABM by Bin Wang on May 2020 in Dalhousie university
!
! references:
!  Laurent, A., Fennel, K., Kuhn, A.: An observation-based evaluation and ranking of historical Earth system model
!      simulations in the northwest North Atlantic Ocean, Biogeosciences, 18, 1803–1822,
!      https://doi.org/10.5194/bg-18-1803-2021, 2021
!
!  Kishi, M. J., et  all, 2007: Nemuro - a lower trophic level model for the North Pacific marine ecosystem,
!      Ecological Modelling, 202, 12-25.
!
!  Lima, I. D., Lam, P. J., and Doney, S. C.: Dynamics of particulate organic carbon flux in a global ocean model,
!      Biogeosciences, 11, 1177–1198, https://doi.org/10.5194/bg-11-1177-2014, 2014.
!
!  Yool, A., Popova, E. E., and Anderson, T. R.: Medusa-1.0: a new intermediate complexity plankton ecosystem model
!      for the global domain, Geosci. Model Dev., 4, 381–417, https://doi.org/10.5194/gmd-4-381-2011, 2011
!
! Please go to this website for more information about GoTM-Fabm coding:
! https://github.com/fabm-model/fabm/wiki/Developing-a-new-biogeochemical-model
!-----------------------------------------------------------------------------------------------------------------

MODULE memg_bio_fennel_2p2z

   USE fabm_types

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(type_base_model), PUBLIC :: type_memg_bio_fennel_2p2z
      ! Variable identifiers
      TYPE (type_state_variable_id) :: id_no3, id_nh4
      TYPE (type_state_variable_id) :: id_phyS, id_phyL
      TYPE (type_state_variable_id) :: id_chlS, id_chlL
      TYPE (type_state_variable_id) :: id_zooS, id_zooL
      TYPE (type_state_variable_id) :: id_SDeN, id_LDeN
      TYPE (type_state_variable_id) :: id_Opal, id_Calc
      TYPE (type_state_variable_id) :: id_reg_no3, id_pre_no3
      
      TYPE (type_dependency_id) :: id_temp, id_par, id_oxy
      TYPE (type_dependency_id) :: id_dep ! depth of the middle of cell 
      TYPE (type_horizontal_dependency_id) :: id_I_0, id_lat

      TYPE (type_diagnostic_variable_id) :: id_newPP,id_regPP
      TYPE (type_diagnostic_variable_id) :: id_Zmetabo,id_Zexcret,id_remineS,id_remineL,id_nitrifi,id_fragment
      
      TYPE (type_diagnostic_variable_id) :: id_newPP_Ps,id_newPP_Pl
      
      TYPE (type_diagnostic_variable_id) :: id_fcaco3
      
      ! Parameters
      real(rk) :: I_thNH4          ! Radiation threshold for nitrification inhibition [Watts/m2]
      real(rk) :: D_p5NH4          ! Half-saturation radiation for nitrification inhibition [Watts/m2]
      real(rk) :: NitriR           ! Nitrification rate: oxidation of NH4 to NO3 [1/day]
      real(rk) :: Ps_K_NO3         ! Inverse half-saturation for small phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: Pl_K_NO3         ! Inverse half-saturation for large phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: Ps_K_NH4         ! Inverse half-saturation for small phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Pl_K_NH4         ! Inverse half-saturation for large phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Ps_Vp0           ! Eppley temperature-limited growth parameter for small phytoplankton [d-1]
      real(rk) :: Pl_Vp0           ! Eppley temperature-limited growth parameter for large phytoplankton [d-1]
      real(rk) :: K_ZsPs           ! Zooplankton half-saturation constant for ingestion [1/day] (phyS --> zooS)
      real(rk) :: K_ZlPs           ! Zooplankton half-saturation constant for ingestion [1/day] (phyS --> zooL)
      real(rk) :: K_ZsPl           ! Zooplankton half-saturation constant for ingestion [1/day] (phyL --> zooS)
      real(rk) :: K_ZlPl           ! Zooplankton half-saturation constant for ingestion [1/day] (phyL --> zooL)
      real(rk) :: K_ZlZs           ! Zooplankton half-saturation constant for ingestion [1/day] (zooS --> zooL)
      real(rk) :: K_ZlDs           ! Zooplankton half-saturation constant for ingestion [1/day] (SDeN --> zooL)
      real(rk) :: PhyCN            ! Phytoplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: DiatSiN          ! Diatom (large Phytoplankton) Silica:Nitrogen ratio [mmol Si/mmol N]
      real(rk) :: Ps_PhyIS         ! Small Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]
      real(rk) :: Pl_PhyIS         ! Large Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]      
      real(rk) :: PhySMR           ! Small Phytoplankton mortality rate [1/day]
      real(rk) :: PhyLMR           ! Large Phytoplankton mortality rate [1/day]
      real(rk) :: Ps_Chl2C_m       ! Maximum chlorophyll to carbon ratio for small phytoplankton [mg Chl/mg C]
      real(rk) :: Pl_Chl2C_m       ! Maximum chlorophyll to carbon ratio for large phytoplankton [mg Chl/mg C]     
      real(rk) :: ZooSAE_N         ! Small Zooplankton Nitrogen assimilation efficiency [nondimesnional]
      real(rk) :: ZooLAE_N         ! Large Zooplankton Nitrogen assimilation efficiency [nondimesnional]  
      real(rk) :: ZooSBM           ! Small Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooLBM           ! Large Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooSER           ! Small Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooLER           ! Large Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooSPsGR         ! Zooplankton maximum growth rate [1/day] (phyS --> zooS)
      real(rk) :: ZooSPlGR         ! Zooplankton maximum growth rate [1/day] (phyL --> zooS)
      real(rk) :: ZooLPsGR         ! Zooplankton maximum growth rate [1/day] (phyS --> zooL)
      real(rk) :: ZooLPlGR         ! Zooplankton maximum growth rate [1/day] (phyL --> zooL)
      real(rk) :: ZooLZsGR         ! Zooplankton maximum growth rate [1/day] (zooS --> zooL)
      real(rk) :: ZooLDsGR         ! Zooplankton maximum growth rate [1/day] (SDeN --> zooL)
      real(rk) :: ZooL_GrInPs      ! Large Zooplankton inhibition coefficient for grazing on Ps [(mmol N m-3)-1]
      real(rk) :: ZooS_GrInPl      ! Small Zooplankton inhibition coefficient for grazing on Pl [(mmol N m-3)-1]
      real(rk) :: ZooL_GrInDs      ! Large Zooplankton inhibition coefficient for grazing on Ds [(mmol N m-3)-1]
      real(rk) :: ZooSMR           ! Small Zooplankton mortality rate [1/day]
      real(rk) :: ZooLMR           ! Large Zooplankton mortality rate [1/day]
      real(rk) :: ZooCN            ! Zooplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: LDeRRN           ! Large detritus remineralization rate N-fraction [1/day]
      real(rk) :: FragRN           ! Large detritus fragmentation rate N-fraction [1/day]
      real(rk) :: SDeRRN           ! Small detritus remineralization rate N-fraction [1/day]
      real(rk) :: CoagR            ! Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN [1/day]
      ! real(rk) :: wP               ! Sinking velocity of phytoplankton [m/day]
      ! real(rk) :: wPS              ! Sinking velocity of small phytoplankton [m/day]
      ! real(rk) :: wPL              ! Sinking velocity of large phytoplankton [m/day]
      ! real(rk) :: wL               ! Sinking velocity of large detritus [m/day]
      ! real(rk) :: wS               ! Sinking velocity of small detritus [m/day]
      real(rk) :: OpalPR           ! Opal protection ratio [mmol N/mmol Si]
      real(rk) :: OpalDR           ! Opal dissolution rate Si-fraction [1/day]
      real(rk) :: fcaco3_0         ! Equatorial Calcite:organic C ratio [mmol Ca/mmol C]
      real(rk) :: fcaco3_90        ! Polar Calcite:organic C ratio [mmol Ca/mmol C]
      real(rk) :: CalcPR           ! Calcite (CaCO3) protection ratio [mmol N/mmol Ca]
      real(rk) :: CalcDR           ! Calcite (CaCO3) dissolution rate Ca-fraction [1/day]
      real(rk) :: PhyMin           ! Phytoplankton minimum threshold value [mmol N/m3]
      real(rk) :: ChlMin           ! Chlorophyll minimum threshold value [mmol N/m3]
      real(rk) :: ZooMin           ! Zooplankton minimum threshold value [mmol N/m3]
      
      real(rk) :: Kz               ! The half saturation constant for ramping up remineralization rate of LDeN
         
      logical  :: oxy_model        ! whether to model oxygen
      logical  :: TEMP_RATES       ! whether the phytoplankton and zooplankton mortality rates,
                                   ! zooplankton grazing rate,
                                   ! zooplankton metabolism and excrete rates,
                                   ! are temperature dependent
      
   CONTAINS
      PROCEDURE :: initialize
      ! Reference model procedures here.
      PROCEDURE :: do
      PROCEDURE :: do_bottom ! instance remineralization for the bottom layer
      
   END TYPE

CONTAINS

   SUBROUTINE initialize(self,configunit)
      CLASS (type_memg_bio_fennel_2p2z), INTENT(INOUT), TARGET :: self
      INTEGER,                          INTENT(IN)              :: configunit
 
      ! Default initial values
      real(rk) :: no3_initial      ! Initial Nitrate Concentration
      real(rk) :: nh4_initial      ! Initial Ammonia Concentration
      
      real(rk) :: phyS_initial     ! Initial small Phytoplankton Concentration
      real(rk) :: phyL_initial     ! Initial large Phytoplankton Concentration
     
      real(rk) :: chlS_initial     ! Initial small Chlorophyll Concentration
      real(rk) :: chlL_initial     ! Initial large Chlorophyll Concentration
      
      real(rk) :: zooS_initial     ! Initial small Zooplankton Concentration
      real(rk) :: zooL_initial     ! Initial large Zooplankton Concentration
      
      real(rk) :: SDeN_initial     ! Initial small Detritus Concentration in nitrogen
      real(rk) :: LDeN_initial     ! Initial large Detritus Concentration in nitrogen
      
      real(rk) :: Opal_initial     ! Initial Opal Concentration in silica
      real(rk) :: Calc_initial     ! Initial CaCO3 Concentration in Ca
      
      real(rk) :: reg_no3_initial      ! Initial regenerated Nitrate Concentration
      real(rk) :: pre_no3_initial      ! Initial preformed Nitrate Concentration
      
      ! Parameters
      real(rk) :: I_thNH4          ! Radiation threshold for nitrification inhibition [Watts/m2]
      real(rk) :: D_p5NH4          ! Half-saturation radiation for nitrification inhibition [Watts/m2]
      real(rk) :: NitriR           ! Nitrification rate: oxidation of NH4 to NO3 [1/day]
      real(rk) :: Ps_K_NO3         ! Inverse half-saturation for small phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: Pl_K_NO3         ! Inverse half-saturation for large phytoplankton NO3 uptake [1/(mmol N m-3)]
      real(rk) :: Ps_K_NH4         ! Inverse half-saturation for small phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Pl_K_NH4         ! Inverse half-saturation for large phytoplankton NH4 uptake [1/(mmol N m-3)]
      real(rk) :: Ps_Vp0           ! Eppley temperature-limited growth parameter for small phytoplankton [d-1]
      real(rk) :: Pl_Vp0           ! Eppley temperature-limited growth parameter for large phytoplankton [d-1]
      real(rk) :: K_ZsPs           ! Zooplankton half-saturation constant for ingestion [1/day] (phyS --> zooS)
      real(rk) :: K_ZlPs           ! Zooplankton half-saturation constant for ingestion [1/day] (phyS --> zooL)
      real(rk) :: K_ZsPl           ! Zooplankton half-saturation constant for ingestion [1/day] (phyL --> zooS)
      real(rk) :: K_ZlPl           ! Zooplankton half-saturation constant for ingestion [1/day] (phyL --> zooL)
      real(rk) :: K_ZlZs           ! Zooplankton half-saturation constant for ingestion [1/day] (zooS --> zooL)
      real(rk) :: K_ZlDs           ! Zooplankton half-saturation constant for ingestion [1/day] (SDeN --> zooL)
      real(rk) :: PhyCN            ! Phytoplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: DiatSiN          ! Diatom (large Phytoplankton) Silica:Nitrogen ratio [mmol Si/mmol N]
      real(rk) :: Ps_PhyIS         ! Small Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]
      real(rk) :: Pl_PhyIS         ! Large Phytoplankton, initial slope of P-I curve [mg_C/(mg_Chl Watts m-2 day)]      
      real(rk) :: PhySMR           ! Small Phytoplankton mortality rate [1/day]
      real(rk) :: PhyLMR           ! Large Phytoplankton mortality rate [1/day]
      real(rk) :: Ps_Chl2C_m       ! Maximum chlorophyll to carbon ratio for small phytoplankton [mg Chl/mg C]
      real(rk) :: Pl_Chl2C_m       ! Maximum chlorophyll to carbon ratio for large phytoplankton [mg Chl/mg C]     
      real(rk) :: ZooSAE_N         ! Small Zooplankton Nitrogen assimilation efficiency [nondimesnional]
      real(rk) :: ZooLAE_N         ! Large Zooplankton Nitrogen assimilation efficiency [nondimesnional]  
      real(rk) :: ZooSBM           ! Small Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooLBM           ! Large Zooplankton Basal metabolism [1/day]
      real(rk) :: ZooSER           ! Small Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooLER           ! Large Zooplankton specific excretion rate [1/day]
      real(rk) :: ZooSPsGR         ! Zooplankton maximum growth rate [1/day] (phyS --> zooS)
      real(rk) :: ZooSPlGR         ! Zooplankton maximum growth rate [1/day] (phyL --> zooS)
      real(rk) :: ZooLPsGR         ! Zooplankton maximum growth rate [1/day] (phyS --> zooL)
      real(rk) :: ZooLPlGR         ! Zooplankton maximum growth rate [1/day] (phyL --> zooL)
      real(rk) :: ZooLZsGR         ! Zooplankton maximum growth rate [1/day] (zooS --> zooL)
      real(rk) :: ZooLDsGR         ! Zooplankton maximum growth rate [1/day] (SDeN --> zooL)
      real(rk) :: ZooL_GrInPs      ! Large Zooplankton inhibition coefficient for grazing on Ps [(mmol N m-3)-1]
      real(rk) :: ZooS_GrInPl      ! Small Zooplankton inhibition coefficient for grazing on Pl [(mmol N m-3)-1]
      real(rk) :: ZooL_GrInDs      ! Large Zooplankton inhibition coefficient for grazing on Ds [(mmol N m-3)-1]
      real(rk) :: ZooSMR           ! Small Zooplankton mortality rate [1/day]
      real(rk) :: ZooLMR           ! Large Zooplankton mortality rate [1/day]
      real(rk) :: ZooCN            ! Zooplankton Carbon:Nitrogen ratio [mmol C/mmol N]
      real(rk) :: LDeRRN           ! Large detritus remineralization rate N-fraction [1/day]
      real(rk) :: FragRN           ! Large detritus fragmentation rate N-fraction [1/day]
      real(rk) :: SDeRRN           ! Small detritus remineralization rate N-fraction [1/day]
      real(rk) :: CoagR            ! Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN [1/day]
      ! real(rk) :: wP               ! Sinking velocity of phytoplankton [m/day]
      ! real(rk) :: wPS              ! Sinking velocity of small phytoplankton [m/day]
      ! real(rk) :: wPL              ! Sinking velocity of large phytoplankton [m/day]
      ! real(rk) :: wL               ! Sinking velocity of large detritus [m/day]
      ! real(rk) :: wS               ! Sinking velocity of small detritus [m/day]
      real(rk) :: OpalPR           ! Opal protection ratio [mmol N/mmol Si]
      real(rk) :: OpalDR           ! Opal dissolution rate Si-fraction [1/day]
      real(rk) :: fcaco3_0         ! Equatorial Calcite:organic C ratio [mmol Ca/mmol C]
      real(rk) :: fcaco3_90        ! Polar Calcite:organic C ratio [mmol Ca/mmol C]
      real(rk) :: CalcPR           ! Calcite (CaCO3) protection ratio [mmol N/mmol Ca]
      real(rk) :: CalcDR           ! Calcite (CaCO3) dissolution rate Ca-fraction [1/day]
      real(rk) :: PhyMin           ! Phytoplankton minimum threshold value [mmol N/m3]
      real(rk) :: ChlMin           ! Chlorophyll minimum threshold value [mmol N/m3]
      real(rk) :: ZooMin           ! Zooplankton minimum threshold value [mmol N/m3]
      
      real(rk) :: Kz               ! The half saturation constant for ramping up remineralization rate of LDeN
         
      logical  :: oxy_model        ! whether to model oxygen
      logical  :: TEMP_RATES       ! whether the phytoplankton and zooplankton mortality rates,
                                   ! zooplankton grazing rate,
                                   ! zooplankton metabolism and excrete rates,
                                   ! are temperature dependent
      
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      namelist /memg_bio_fennel_2p2z/  no3_initial,nh4_initial, &
                                       phyS_initial,phyL_initial,    &
                                       chlS_initial,chlL_initial,    &
                                       zooS_initial,zooL_initial,    &
                                       SDeN_initial,LDeN_initial,   &
                                       Opal_initial,Calc_initial,   &
                                       reg_no3_initial,pre_no3_initial,   &
                                       I_thNH4,D_p5NH4,NitriR,    &
                                       Ps_K_NO3,Pl_K_NO3,Ps_K_NH4,Pl_K_NH4,    &
                                       Ps_Vp0,Pl_Vp0,    &
                                       K_ZsPs,K_ZlPs,K_ZsPl,K_ZlPl,K_ZlZs,K_ZlDs,    &
                                       PhyCN,Ps_PhyIS,Pl_PhyIS,PhySMR,PhyLMR,Ps_Chl2C_m,Pl_Chl2C_m,    &
                                       ZooSAE_N,ZooLAE_N,ZooSBM,ZooLBM,ZooSER,ZooLER,    &
                                       ZooSPsGR,ZooSPlGR,ZooLPsGR,ZooLPlGR,ZooLZsGR,ZooLDsGR,    &
                                       ZooL_GrInPs,ZooS_GrInPl,ZooL_GrInDs,ZooSMR,ZooLMR,ZooCN,    &
                                       LDeRRN,FragRN,SDeRRN,CoagR,   &
                                       DiatSiN, OpalPR, OpalDR, &
                                       fcaco3_0,fcaco3_90,CalcPR,CalcDR, &
                                       PhyMin,ChlMin,ZooMin,Kz
                                       
!
! set default values for initial condition
!
      no3_initial = 1.0_rk      
      nh4_initial = 0.001_rk     
      
      phyS_initial = 0.001_rk
      phyL_initial = 0.001_rk
      
      chlS_initial = 0.001_rk
      chlL_initial = 0.001_rk
      
      zooS_initial = 0.001_rk
      zooL_initial = 0.001_rk
      
      SDeN_initial = 0.001_rk   
      LDeN_initial = 0.001_rk
      
      Opal_initial = 0.001_rk
      Calc_initial = 0.001_rk
      
      reg_no3_initial = 0.001_rk
      pre_no3_initial = 0.001_rk
      
!
! set default values for parameter
!
      I_thNH4 = 0.0095_rk 
      D_p5NH4 = 0.1_rk      
      NitriR = 0.2_rk       
      Ps_K_NO3 = 2.0_rk   
      Pl_K_NO3 = 2.0_rk         
      Ps_K_NH4 = 2.0_rk         
      Pl_K_NH4 = 2.0_rk        
      Ps_Vp0 = 1.9711_rk         
      Pl_Vp0 = 1.9054_rk          
      K_ZsPs = 2.0_rk          
      K_ZlPs = 2.0_rk          
      K_ZsPl = 2.0_rk         
      K_ZlPl = 2.0_rk           
      K_ZlZs = 2.0_rk           
      K_ZlDs = 2.0_rk            
      PhyCN = 6.625_rk
      DiatSiN = 2.0_rk
      Ps_PhyIS = 0.0405_rk       
      Pl_PhyIS = 0.0393_rk             
      PhySMR = 0.2377_rk         
      PhyLMR = 0.1169_rk          
      Ps_Chl2C_m = 0.0328_rk       
      Pl_Chl2C_m = 0.0386_rk          
      ZooSAE_N = 0.75_rk        
      ZooLAE_N = 0.75_rk          
      ZooSBM = 0.0886_rk
      ZooLBM = 0.0886_rk           
      ZooSER = 0.0886_rk           
      ZooLER = 0.0886_rk        
      ZooSPsGR = 6.6761_rk        
      ZooSPlGR = 6.6761_rk        
      ZooLPsGR = 3.33805_rk  
      ZooLPlGR = 1.1126_rk    
      ZooLZsGR = 6.6761_rk              
      ZooLDsGR = 0.0_rk         
      ZooL_GrInPs = 3.010_rk     
      ZooS_GrInPl = 3.010_rk         
      ZooL_GrInDs = 4.605_rk      
      ZooSMR = 0.0224_rk        
      ZooLMR = 0.0224_rk         
      ZooCN = 6.625_rk           
      LDeRRN = 0.01_rk
      FragRN = 0.01_rk
      SDeRRN = 0.4_rk          
      CoagR = 0.0023_rk
      ! wP  = 0.1_rk
      ! wPS = 0.1_rk              
      ! wPL = 0.5_rk             
      ! wL = 5.0_rk              
      ! wS = 0.1_rk
      OpalPR = 0.02_rk
      OpalDR = 0.0225_rk
      fcaco3_0 = 0.10_rk
      fcaco3_90 = 0.02_rk
      CalcPR = 0.088_rk
      CalcDR = 0.0_rk
      PhyMin = 0.00000001_rk           
      ChlMin = 0.00000001_rk           
      ZooMin = 0.00000001_rk
      
      Kz = 0.0_rk
         
      oxy_model = .false.        
      TEMP_RATES  = .true.     
      
      ! Read the namelist
      if (configunit>0) read(configunit,nml=memg_bio_fennel_2p2z,err=99,end=100)

      ! Store parameter values into parameters
      call self%get_parameter(self%I_thNH4,'I_thNH4','Watts.m-2','Radiation threshold for nitrification inhibition',default=I_thNH4)
      call self%get_parameter(self%D_p5NH4,'D_p5NH4','Watts.m-2','Half-saturation radiation for nitrification inhibition',default=D_p5NH4)
      call self%get_parameter(self%NitriR,'NitriR','d-1','Nitrification rate: oxidation of NH4 to NO3',default=NitriR)
      call self%get_parameter(self%Ps_K_NO3,'Ps_K_NO3','(mmol N m-3)-1','Inverse half-saturation for small phytoplankton NO3 uptake',default=Ps_K_NO3)
      call self%get_parameter(self%Pl_K_NO3,'Pl_K_NO3','(mmol N m-3)-1','Inverse half-saturation for large phytoplankton NO3 uptake',default=Pl_K_NO3)
      call self%get_parameter(self%Ps_K_NH4,'Ps_K_NH4','(mmol N m-3)-1','Inverse half-saturation for small phytoplankton NH4 uptake',default=Ps_K_NH4)
      call self%get_parameter(self%Pl_K_NH4,'Pl_K_NH4','(mmol N m-3)-1','Inverse half-saturation for large phytoplankton NH4 uptake',default=Pl_K_NH4)
      call self%get_parameter(self%Ps_Vp0,'Ps_Vp0','d-1','Eppley temperature-limited growth parameter for small phytoplankton',default=Ps_Vp0)
      call self%get_parameter(self%Pl_Vp0,'Pl_Vp0','d-1','Eppley temperature-limited growth parameter for large phytoplankton',default=Pl_Vp0)
      call self%get_parameter(self%K_ZsPs,'K_ZsPs','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (phyS --> zooS)',default=K_ZsPs)
      call self%get_parameter(self%K_ZlPs,'K_ZlPs','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (phyS --> zooL)',default=K_ZlPs)
      call self%get_parameter(self%K_ZsPl,'K_ZsPl','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (phyL --> zooS)',default=K_ZsPl)
      call self%get_parameter(self%K_ZlPl,'K_ZlPl','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (phyL --> zooL)',default=K_ZlPl)
      call self%get_parameter(self%K_ZlZs,'K_ZlZs','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (zooS --> zooL)',default=K_ZlZs)
      call self%get_parameter(self%K_ZlDs,'K_ZlDs','(millimole_N m-3)2','Zooplankton half-saturation constant for ingestion (SDeN --> zooL)',default=K_ZlDs)
      call self%get_parameter(self%PhyCN,'PhyCN','mmol C/mmol N','Phytoplankton Carbon:Nitrogen ratio',default=PhyCN)
      call self%get_parameter(self%DiatSiN,'DiatSiN','mmol Si/mmol N','Diatom (large Phytoplankton) Silica:Nitrogen ratio',default=DiatSiN)
      call self%get_parameter(self%Ps_PhyIS,'Ps_PhyIS','(Watts m-2 day)-1','Small Phytoplankton, initial slope of P-I curve',default=Ps_PhyIS)
      call self%get_parameter(self%Pl_PhyIS,'Pl_PhyIS','(Watts m-2 day)-1','Large Phytoplankton, initial slope of P-I curve',default=Pl_PhyIS)
      call self%get_parameter(self%PhySMR,'PhySMR','d-1','Small Phytoplankton mortality rate',default=PhySMR)
      call self%get_parameter(self%PhyLMR,'PhyLMR','d-1','Large Phytoplankton mortality rate',default=PhyLMR)
      call self%get_parameter(self%Ps_Chl2C_m,'Ps_Chl2C_m','mg Chl/mg C','Maximum chlorophyll to carbon ratio for small phytoplankton',default=Ps_Chl2C_m)
      call self%get_parameter(self%Pl_Chl2C_m,'Pl_Chl2C_m','mg Chl/mg C','Maximum chlorophyll to carbon ratio for large phytoplankton',default=Pl_Chl2C_m)
      call self%get_parameter(self%ZooSAE_N,'ZooSAE_N','-','Small Zooplankton Nitrogen assimilation efficiency',default=ZooSAE_N)
      call self%get_parameter(self%ZooLAE_N,'ZooLAE_N','-','Large Zooplankton Nitrogen assimilation efficiency',default=ZooLAE_N)
      call self%get_parameter(self%ZooSBM,'ZooSBM','d-1','Small Zooplankton Basal metabolism',default=ZooSBM)
      call self%get_parameter(self%ZooLBM,'ZooLBM','d-1','Large Zooplankton Basal metabolism',default=ZooLBM)
      call self%get_parameter(self%ZooSER,'ZooSER','d-1','Small Zooplankton specific excretion rate',default=ZooSER)
      call self%get_parameter(self%ZooLER,'ZooLER','d-1','Large Zooplankton specific excretion rate',default=ZooLER)
      call self%get_parameter(self%ZooSPsGR,'ZooSPsGR','d-1','Zooplankton maximum growth rate (phyS --> zooS)',default=ZooSPsGR)
      call self%get_parameter(self%ZooSPlGR,'ZooSPlGR','d-1','Zooplankton maximum growth rate (phyL --> zooS)',default=ZooSPlGR)
      call self%get_parameter(self%ZooLPsGR,'ZooLPsGR','d-1','Zooplankton maximum growth rate (phyS --> zooL)',default=ZooLPsGR)
      call self%get_parameter(self%ZooLPlGR,'ZooLPlGR','d-1','Zooplankton maximum growth rate (phyL --> zooL)',default=ZooLPlGR)
      call self%get_parameter(self%ZooLZsGR,'ZooLZsGR','d-1','Zooplankton maximum growth rate (zooS --> zooL)',default=ZooLZsGR)
      call self%get_parameter(self%ZooLDsGR,'ZooLDsGR','d-1','Zooplankton maximum growth rate (SDeN --> zooL)',default=ZooLDsGR)
      call self%get_parameter(self%ZooL_GrInPs,'ZooL_GrInPs','(mmol N m-3)-1','Large Zooplankton inhibition coefficient for grazing on Ps',default=ZooL_GrInPs)
      call self%get_parameter(self%ZooS_GrInPl,'ZooS_GrInPl','(mmol N m-3)-1','Small Zooplankton inhibition coefficient for grazing on Pl',default=ZooS_GrInPl)
      call self%get_parameter(self%ZooL_GrInDs,'ZooL_GrInDs','(mmol N m-3)-1','Large Zooplankton inhibition coefficient for grazing on Ds',default=ZooL_GrInDs)
      call self%get_parameter(self%ZooSMR,'ZooSMR','d-1','Small Zooplankton mortality rate',default=ZooSMR)
      call self%get_parameter(self%ZooLMR,'ZooLMR','d-1','Large Zooplankton mortality rate',default=ZooLMR)
      call self%get_parameter(self%ZooCN,'ZooCN','mmol C/mmol N','Zooplankton Carbon:Nitrogen ratio',default=ZooCN)
      call self%get_parameter(self%LDeRRN,'LDeRRN','d-1','Large detritus remineralization rate N-fraction',default=LDeRRN)
      call self%get_parameter(self%FragRN,'FragRN','d-1','Large detritus fragmentation rate N-fraction',default=FragRN)
      call self%get_parameter(self%SDeRRN,'SDeRRN','d-1','Small detritus remineralization rate N-fraction',default=SDeRRN)
      call self%get_parameter(self%CoagR,'CoagR','d-1','Coagulation rate: aggregation rate of SDeN + Phy ==> LDeN',default=CoagR)
      ! call self%get_parameter(self%wP,'wP','m d-1','Sinking velocity of small phytoplankton (<0 for sinking)',default=wP,scale_factor=-1.0_rk*d_per_s)
      ! call self%get_parameter(self%wPS,'wPS','m d-1','Sinking velocity of small phytoplankton (<0 for sinking)',default=wPS,scale_factor=-1.0_rk*d_per_s)
      ! call self%get_parameter(self%wPL,'wPL','m d-1','Sinking velocity of large phytoplankton (<0 for sinking)',default=wPL,scale_factor=-1.0_rk*d_per_s)
      ! call self%get_parameter(self%wL,'wL','m d-1','Sinking velocity of large detritus (<0 for sinking)',default=wL,scale_factor=-1.0_rk*d_per_s)
      ! call self%get_parameter(self%wS,'wS','m d-1','Sinking velocity of samll detritus (<0 for sinking)',default=wS,scale_factor=-1.0_rk*d_per_s)
      call self%get_parameter(self%OpalPR,'OpalPR','mmol N/mmol Si','Opal protection ratio',default=OpalPR)
      call self%get_parameter(self%OpalDR,'OpalDR','d-1','Opal dissolution rate Si-fraction',default=OpalDR)
      call self%get_parameter(self%fcaco3_0,'fcaco3_0','mmol Ca/mmol C','Equatorial Calcite:organic C ratio',default=fcaco3_0)
      call self%get_parameter(self%fcaco3_90,'fcaco3_90','mmol Ca/mmol C','Polar Calcite:organic C ratio',default=fcaco3_90)
      call self%get_parameter(self%CalcPR,'CalcPR','mmol N/mmol Ca','Calcite (CaCO3) protection ratio',default=CalcPR)
      call self%get_parameter(self%CalcDR,'CalcDR','d-1','Calcite (CaCO3) dissolution rate Ca-fraction',default=CalcDR)
      call self%get_parameter(self%PhyMin,'PhyMin','mmol N m-3','Phytoplankton minimum threshold value',default=PhyMin)
      call self%get_parameter(self%ChlMin,'ChlMin','mmol N m-3','Chlorophyll minimum threshold value',default=ChlMin)
      call self%get_parameter(self%ZooMin,'ZooMin','mmol N m-3','Zooplankton minimum threshold value',default=ZooMin)
      
      call self%get_parameter(self%Kz,'Kz','m','The half saturation constant for ramping up remineralization rate of LDeN',default=Kz)
      
      call self%get_parameter(self%oxy_model,'oxy_model','','Switch to modelling oxygen',default=oxy_model)
      call self%get_parameter(self%TEMP_RATES,'TEMP_RATES','','Switch to temperature dependent rate',default=TEMP_RATES)


      ! Register state variables
      call self%register_state_variable(self%id_no3,'no3','mmol N m-3','nitrate',no3_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='nitrate',units='mmol N m-3'))
      call self%register_state_variable(self%id_nh4,'nh4','mmol N m-3','ammonia',nh4_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='ammonia',units='mmol N m-3'))
      
      call self%register_state_variable(self%id_phyS,'phyS','mmol N m-3','small_phytoplankton',phyS_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='small_phytoplankton',units='mmol N m-3'))
      call self%register_state_variable(self%id_phyL,'phyL','mmol N m-3','large_phytoplankton',phyL_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='large_phytoplankton',units='mmol N m-3'))
     
      call self%register_state_variable(self%id_chlS,'chlS','mg m-3','small_chlorophyll',chlS_initial,minimum=0.0_rk,    &
         &    standard_variable=type_bulk_standard_variable(name='small_chlorophyll',units='mg m-3'))
      call self%register_state_variable(self%id_chlL,'chlL','mg m-3','large_chlorophyll',chlL_initial,minimum=0.0_rk,    &
         &    standard_variable=type_bulk_standard_variable(name='large_chlorophyll',units='mg m-3'))
         
      call self%register_state_variable(self%id_zooS,'zooS','mmol N m-3','small_zooplankton',zooS_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='small_zooplankton',units='mmol N m-3'))
      call self%register_state_variable(self%id_zooL,'zooL','mmol N m-3','large_zooplankton',zooL_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='large_zooplankton',units='mmol N m-3'))

      call self%register_state_variable(self%id_LDeN,'LDeN','mmol N m-3','LdetritusN',LDeN_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='LdetritusN',units='mmol N m-3'))
      call self%register_state_variable(self%id_SDeN,'SDeN','mmol N m-3','SdetritusN',SDeN_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='SdetritusN',units='mmol N m-3'))
         
      call self%register_state_variable(self%id_Opal,'Opal','mmol Si m-3','Opal',Opal_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='Opal',units='mmol Si m-3'))
      call self%register_state_variable(self%id_Calc,'Calcite','mmol Ca m-3','Calcite',Calc_initial,minimum=0.0_rk,   &
         &    standard_variable=type_bulk_standard_variable(name='Calcite',units='mmol Ca m-3'))
         
      call self%register_state_variable(self%id_reg_no3,'reg_no3','mmol N m-3','regenerate_nitrate',reg_no3_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='regenerate_nitrate',units='mmol N m-3'))
      call self%register_state_variable(self%id_pre_no3,'pre_no3','mmol N m-3','preformed_nitrate', pre_no3_initial,minimum=0.0_rk,no_precipitation_dilution=.true.,   &
         &    standard_variable=type_bulk_standard_variable(name='preformed_nitrate',units='mmol N m-3'))
         
      ! Register environmental dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_par, standard_variable=type_bulk_standard_variable(name='photosynthetic_active_radiation_in_water_column',units='W m-2')) ! calculated by the light model
      call self%register_dependency(self%id_lat, standard_variables%latitude)
      call self%register_dependency(self%id_dep, standard_variables%depth) ! Depth of middle of cell
        
      if (self%oxy_model) then
         call self%register_dependency(self%id_oxy, standard_variable=type_bulk_standard_variable(name='oxygen',units='mmol O m-3')) ! calculated by the oxygen model
      end if
      
      ! Register the contribution of all state variables to total chlorophyll
      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='chlorophyll',units="mg m-3",aggregate_variable=.true.),self%id_chlS)
      call self%add_to_aggregate_variable(type_bulk_standard_variable(name='chlorophyll',units="mg m-3",aggregate_variable=.true.),self%id_chlL)
      
      ! Register the contribution of all state variables to total phytoplankton
      ! call self%add_to_aggregate_variable(type_bulk_standard_variable(name='phytoplankton',units="mmol N m-3",aggregate_variable=.true.),self%id_phyS)
      ! call self%add_to_aggregate_variable(type_bulk_standard_variable(name='phytoplankton',units="mmol N m-3",aggregate_variable=.true.),self%id_phyL)
      
      ! Register the contribution of all state variables to total zooplankotn
      ! call self%add_to_aggregate_variable(type_bulk_standard_variable(name='zooplankton',units="mmol N m-3",aggregate_variable=.true.),self%id_zooS)
      ! call self%add_to_aggregate_variable(type_bulk_standard_variable(name='zooplankton',units="mmol N m-3",aggregate_variable=.true.),self%id_zooL)
      
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_newPP,'NewPP','mmol N m-3 day-1','new primary production by NO3',   &
         &    standard_variable=type_bulk_standard_variable(name='NewPP',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_regPP,'RegPP','mmol N m-3 day-1','regenerate primary production by NH4',   &
         &    standard_variable=type_bulk_standard_variable(name='RegPP',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Zmetabo,'Zmetabo','mmol N m-3 day-1','zooplankton metabolism',   &
         &    standard_variable=type_bulk_standard_variable(name='Zmetabo',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_Zexcret,'Zexcret','mmol N m-3 day-1','zooplankton excretion',   &
         &    standard_variable=type_bulk_standard_variable(name='Zexcret',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_remineS,'RemineS','mmol N m-3 day-1','remineralization of small detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='RemineS',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_fragment,'Fragment','mmol N m-3 day-1','fragmentation of large detritus into small detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='Fragment',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_remineL,'RemineL','mmol N m-3 day-1','remineralization of large detritus',   &
         &    standard_variable=type_bulk_standard_variable(name='RemineL',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_nitrifi,'Nitrifi','mmol N m-3 day-1','nitrification from NH4 into NO3',   &
         &    standard_variable=type_bulk_standard_variable(name='Nitrifi',units='mmol N m-3 day-1'))
      
      call self%register_diagnostic_variable(self%id_newPP_Ps,'NewPP_Ps','mmol N m-3 day-1','new primary production by NO3 of small phytoplankton',   &
         &    standard_variable=type_bulk_standard_variable(name='NewPP_Ps',units='mmol N m-3 day-1'))
      call self%register_diagnostic_variable(self%id_newPP_Pl,'NewPP_Pl','mmol N m-3 day-1','new primary production by NO3 of large phytoplankton',   &
         &    standard_variable=type_bulk_standard_variable(name='NewPP_Pl',units='mmol N m-3 day-1'))
      
      call self%register_diagnostic_variable(self%id_fcaco3,'fcaco3','mmol Ca/mmol C','Calcite:organic C ratio',   &
         &    standard_variable=type_bulk_standard_variable(name='fcaco3',units='mmol Ca/mmol C'))
         

   return

99 call self%fatal_error('memg_bio_fennel_2p2z','Error reading namelist memg_bio_fennel_2p2z.')

100 call self%fatal_error('memg_bio_fennel_2p2z','Namelist memg_bio_fennel_2p2z was not found.')

   ! Register model parameters and variables here.
   END SUBROUTINE initialize

   SUBROUTINE do(self,_ARGUMENTS_DO_)
      class (type_memg_bio_fennel_2p2z),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
   
      real(rk) :: no3,nh4,phyS,phyL,chlS,chlL,zooS,zooL,LDeN,SDeN,Opal,Calcite, reg_no3, pre_no3
      real(rk) :: temp,par,oxy,I_0,lat,depth 
      
      ! Local variables
      real(rk) :: Ps_Chl2C,Pl_Chl2C,Vp,Epp,t_PPmax
      real(rk) :: inhNH4,L_NH4,L_NO3,LTOT
      real(rk) :: cff,cff1,cff2,cff3,cff4,cff5
      real(rk) :: fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8
      real(rk) :: ppref
      real(rk) :: switch_ZlDs
      
      real(rk) :: N_NewProd,N_NewProd_Ps,N_NewProd_Pl
      real(rk) :: N_RegProd,N_RegProd_Ps,N_RegProd_Pl
      real(rk) :: Chl_Prod,Chl_Prod_Ps,Chl_Prod_Pl
      real(rk) :: N_Nitrifi
      real(rk) :: N_Graz,N_Graz_ZsPs,N_Graz_ZlPs,N_Graz_ZsPl,N_Graz_ZlPl
      real(rk) :: Chl_Graz,Chl_Graz_ZsPs,Chl_Graz_ZlPs,Chl_Graz_ZsPl,Chl_Graz_ZlPl
      real(rk) :: N_Assim,N_Assim_Zs,N_Assim_Zl
      real(rk) :: N_Assim_ZsPs,N_Assim_ZlPs,N_Assim_ZsPl,N_Assim_ZlPl,N_Assim_ZlZs,N_Assim_ZlDs
      real(rk) :: N_Egest,N_Egest_Zs,N_Egest_Zl
      real(rk) :: N_Egest_ZsPs,N_Egest_ZlPs,N_Egest_ZsPl,N_Egest_ZlPl,N_Egest_ZlZs,N_Egest_ZlDs
      real(rk) :: N_Pmortal,N_Pmortal_Ps,N_Pmortal_Pl
      real(rk) :: Chl_Pmortal,Chl_Pmortal_Ps,Chl_Pmortal_Pl
      real(rk) :: N_Zmortal,N_Zsmortal,N_Zlmortal
      real(rk) :: N_Zexcret,N_Zsexcret,N_Zlexcret
      real(rk) :: N_Zmetabo,N_Zsmetabo,N_Zlmetabo
      real(rk) :: N_CoagP,N_CoagD,Chl_Coag
      real(rk) :: N_RemineS,N_RemineL,N_Fragment
      real(rk) :: Opal_Dissolve,Calc_Dissolve
      real(rk) :: fcaco3
      real(rk) :: d_no3,d_nh4
      real(rk) :: d_phyS,d_phyL
      real(rk) :: d_zooS,d_zooL
      real(rk) :: d_chlS,d_chlL
      real(rk) :: d_LDeN,d_SDeN
      real(rk) :: d_Opal,d_Calc
      real(rk) :: d_reg_no3,d_pre_no3
      ! free part and protected part of LDeN
      real(rk) :: fLDeN,pLDeN
      
      real(rk) :: zfact ! coefficient to ramp remineraization rate rL = rL*zfact
      
      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk), parameter :: eps = 1.0e-20_rk
      
      ! parameters associated with temperature dependent rates
      real(rk), parameter :: Tcoef0 = 0.59_rk 
      real(rk), parameter :: Tcoef1 = 1.066_rk
      
      ! parameters associated with modeling oxygen
      real(rk), parameter :: rOxNO3= 8.625_rk       ! 138/16
      real(rk), parameter :: rOxNH4= 6.625_rk       ! 106/16
      
      ! the remineralization rate of LDeN will be ramped to its defined value at zref
      real(rk), parameter :: zref= 100.0_rk  ! 
      
      _LOOP_BEGIN_
      
         ! Obtain concentration of biological variables.
         _GET_(self%id_no3,no3)
         _GET_(self%id_nh4,nh4)
         
         _GET_(self%id_phyS,phyS)
         _GET_(self%id_phyL,phyL)
         
         _GET_(self%id_chlS,chlS)
         _GET_(self%id_chlL,chlL)
         
         _GET_(self%id_zooS,zooS)
         _GET_(self%id_zooL,zooL)
         
         _GET_(self%id_LDeN,LDeN)
         _GET_(self%id_SDeN,SDeN)
         
         _GET_(self%id_Opal,Opal)
         _GET_(self%id_Calc,Calcite)
           
         _GET_(self%id_reg_no3,reg_no3)
         _GET_(self%id_pre_no3,pre_no3)
         
         ! Obtain environmental dependencies
         _GET_(self%id_temp,temp)
         _GET_(self%id_dep, depth) !Depth of middle of cell [m]
         
         if (self%oxy_model) then
            _GET_(self%id_oxy,oxy)
         end if
         
         _GET_(self%id_par,par)
         _GET_HORIZONTAL_(self%id_I_0,I_0)
         _GET_HORIZONTAL_(self%id_lat,lat) ! latitude [degree_north]
         
         if (self%ZooL_GrInDs .le. 0.0_rk) then
            switch_ZlDs = 0.0_rk
         else
            switch_ZlDs = 1.0_rk
         endif
         
!-----------------------------------------------------------------------
! Calculate the latitude dependent rain ratio, fcaco3 [mmol Ca/mmol C]
! following Yool et al (2011).
!-----------------------------------------------------------------------
         fcaco3=self%fcaco3_90+(self%fcaco3_0-self%fcaco3_90)*(90.0_rk-abs(lat))/90.0_rk
         fcaco3 = fcaco3*self%PhyCN
         
         IF (I_0 .gt. 0.0_rk) THEN
!-----------------------------------------------------------------------
! Compute ration of Chlorophyll-a to phytoplankton, [mg Chla / (mg C)]
!-----------------------------------------------------------------------
            cff = self%PhyCN*12.0_rk
            Ps_Chl2C = MIN(chlS/(phyS*cff+eps), self%Ps_Chl2C_m)
            Pl_Chl2C = MIN(chlL/(phyL*cff+eps), self%Pl_Chl2C_m)
 
! small phytoplankton
!-----------------------------------------------------------------------
!  Temperature-limited and light-limited growth rate 
!-----------------------------------------------------------------------
            Vp=self%Ps_Vp0*Tcoef0*(Tcoef1**temp)
            fac1=par*self%Ps_PhyIS
            Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
            t_PPmax=Epp*fac1
!-----------------------------------------------------------------------
!  Nutrient-limitation terms
!-----------------------------------------------------------------------
            cff1=nh4*self%Ps_K_NH4
            cff2=no3*self%Ps_K_NO3
            inhNH4=1.0_rk/(1.0_rk+cff1)
            L_NH4=cff1/(1.0_rk+cff1)
            L_NO3=cff2*inhNH4/(1.0_rk+cff2)
            LTOT=L_NO3+L_NH4
            !  Nitrate and ammonium uptake by small Phytoplankton.
            fac1=d_per_s*t_PPmax
            cff4=fac1*self%Ps_K_NO3*inhNH4/(1.0_rk+cff2)*phyS
            cff5=fac1*self%Ps_K_NH4/(1.0_rk+cff1)*phyS
               
            no3=no3/(1.0_rk+cff4)
            nh4=nh4/(1.0_rk+cff5)
            N_NewProd_Ps=no3*cff4
            N_RegProd_Ps=nh4*cff5
               
            phyS=phyS+N_NewProd_Ps+N_RegProd_Ps

            ! increse in chlorophyll because of small phytoplankton growth 
            Chl_Prod_Ps = (d_per_s*t_PPmax*t_PPmax*LTOT*LTOT*self%Ps_Chl2C_m*chlS)/  &
     &                    (self%Ps_PhyIS*MAX(Ps_Chl2C,eps)*par+eps)
            
            chlS=chlS+Chl_Prod_Ps
               
! large phytoplankton
!-----------------------------------------------------------------------
!  Temperature-limited and light-limited growth rate 
!-----------------------------------------------------------------------
            Vp=self%Pl_Vp0*Tcoef0*(Tcoef1**temp)
            fac1=par*self%Pl_PhyIS
            Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
            t_PPmax=Epp*fac1
!-----------------------------------------------------------------------
!  Nutrient-limitation terms
!-----------------------------------------------------------------------
            cff1=nh4*self%Pl_K_NH4
            cff2=no3*self%Pl_K_NO3 
            inhNH4=1.0_rk/(1.0_rk+cff1)
            L_NH4=cff1/(1.0_rk+cff1)
            L_NO3=cff2*inhNH4/(1.0_rk+cff2)
            LTOT=L_NO3+L_NH4
            !  Nitrate and ammonium uptake by large Phytoplankton.
            fac1=d_per_s*t_PPmax
            cff4=fac1*self%Pl_K_NO3*inhNH4/(1.0_rk+cff2)*phyL
            cff5=fac1*self%Pl_K_NH4/(1.0_rk+cff1)*phyL
               
            no3=no3/(1.0_rk+cff4)
            nh4=nh4/(1.0_rk+cff5)
            N_NewProd_Pl=no3*cff4
            N_RegProd_Pl=nh4*cff5
               
            phyL=phyL+N_NewProd_Pl+N_RegProd_Pl

            ! increse in chlorophyll because of large phytoplankton growth 
            Chl_Prod_Pl = (d_per_s*t_PPmax*t_PPmax*LTOT*LTOT*self%Pl_Chl2C_m*chlL)/  &
     &                    (self%Pl_PhyIS*MAX(Pl_Chl2C,eps)*par+eps)
            
            chlL=chlL+Chl_Prod_Pl
               
            ! sum up fluxes by small and large phytoplankton
            N_NewProd=N_NewProd_Ps+N_NewProd_Pl
            N_RegProd=N_RegProd_Ps+N_RegProd_Pl
            Chl_Prod = Chl_Prod_Ps+Chl_Prod_Pl
               
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
            N_NewProd_Ps = 0.0_rk
            N_NewProd_Pl = 0.0_rk
            N_RegProd_Ps = 0.0_rk
            N_RegProd_Pl = 0.0_rk
            Chl_Prod_Ps = 0.0_rk
            Chl_Prod_Pl = 0.0_rk
            
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
         if (self%TEMP_RATES) then
            fac1=d_per_s*self%ZooSPsGR*Tcoef0*(Tcoef1**temp) ! phyS --> zooS
            fac2=d_per_s*self%ZooSPlGR*Tcoef0*(Tcoef1**temp) ! phyL --> zooS
            fac3=d_per_s*self%ZooLPsGR*Tcoef0*(Tcoef1**temp) ! phyS --> zooL
            fac4=d_per_s*self%ZooLPlGR*Tcoef0*(Tcoef1**temp) ! phyL --> zooL
            fac5=d_per_s*self%ZooLZsGR*Tcoef0*(Tcoef1**temp) ! zooS --> zooL
            fac6=d_per_s*self%ZooLDsGR*Tcoef0*(Tcoef1**temp) ! SDeN --> zooL
         else
            fac1=d_per_s*self%ZooSPsGR ! phyS --> zooS
            fac2=d_per_s*self%ZooSPlGR ! phyL --> zooS
            fac3=d_per_s*self%ZooLPsGR ! phyS --> zooL
            fac4=d_per_s*self%ZooLPlGR ! phyL --> zooL
            fac5=d_per_s*self%ZooLZsGR ! zooS --> zooL
            fac6=d_per_s*self%ZooLDsGR ! SDeN --> zooL
         end if   ! endif TEMP_RATES
         
!
! small phytoplankton grazing by small zooplankton (1).
!
         cff1=fac1*zooS*phyS/(self%K_ZsPs+phyS*phyS)
         cff3=1.0_rk/(1.0_rk+cff1)
         phyS=cff3*phyS
         chlS=cff3*chlS
              
         N_Graz_ZsPs=cff1*phyS
         Chl_Graz_ZsPs=cff1*chlS
              
         N_Assim_ZsPs=cff1*phyS*self%ZooSAE_N
         N_Egest_ZsPs=cff1*phyS*(1.0_rk-self%ZooSAE_N)
            
!
! small phytoplankton grazing by large zooplankton (2).
!
         ppref = EXP(-self%ZooL_GrInPs*(phyL+zooS+SDeN*switch_ZlDs))
         cff1=ppref*fac3*zooL*phyS/(self%K_ZlPs+phyS*phyS)
         cff4=1.0_rk/(1.0_rk+cff1)
         phyS=cff4*phyS
         chlS=cff4*chlS
              
         N_Graz_ZlPs=cff1*phyS
         Chl_Graz_ZlPs=cff1*chlS
              
         N_Assim_ZlPs=cff1*phyS*self%ZooLAE_N
         N_Egest_ZlPs=cff1*phyS*(1.0_rk-self%ZooLAE_N)
              
!
! large phytoplankton grazing by small zooplankton (3).
!
         ppref = EXP(-self%ZooS_GrInPl*phyS)
         cff1=ppref*fac2*zooS*phyL/(self%K_ZsPl+phyL*phyL)
         cff4=1.0_rk/(1.0_rk+cff1)
         phyL=cff4*phyL
         chlL=cff4*chlL
              
         N_Graz_ZsPl=cff1*phyL
         Chl_Graz_ZsPl=cff1*chlL
              
         N_Assim_ZsPl=cff1*phyL*self%ZooSAE_N
         N_Egest_ZsPl=cff1*phyL*(1.0_rk-self%ZooSAE_N)         
!
! large phytoplankton grazing by large zooplankton (4).
!             
         cff1=fac4*zooL*phyL/(self%K_ZlPl+phyL*phyL)
         cff4=1.0_rk/(1.0_rk+cff1)
         phyL=cff4*phyL
         chlL=cff4*chlL
              
         N_Graz_ZlPl=cff1*phyL
         Chl_Graz_ZlPl=cff1*chlL
              
         N_Assim_ZlPl=cff1*phyL*self%ZooLAE_N
         N_Egest_ZlPl=cff1*phyL*(1.0_rk-self%ZooLAE_N)
              
! 
! small zooplankton grazing by large zooplankton (5).
!               
         cff1=fac5*zooL*zooS/(self%K_ZlZs+zooS*zooS)
         cff4=1.0_rk/(1.0_rk+cff1)
         zooS=cff4*zooS
              
         N_Assim_ZlZs=cff1*zooS*self%ZooLAE_N
         N_Egest_ZlZs=cff1*zooS*(1.0_rk-self%ZooLAE_N)
              
!
! small detritus grazing by large zooplankton (6).
!
         ppref = EXP(-self%ZooL_GrInDs*(phyS+phyL+zooS))
         cff1=ppref*fac6*zooL*SDeN/(self%K_ZlDs+SDeN*SDeN)
         cff4=1.0_rk/(1.0_rk+cff1)
         SDeN=cff4*SDeN
         N_Assim_ZlDs=cff1*SDeN*self%ZooLAE_N
         N_Egest_ZlDs=cff1*SDeN*(1.0_rk-self%ZooLAE_N)
            
!
! sum up all assimilation and egestion fluxes.
!
         N_Assim_Zs=N_Assim_ZsPs+N_Assim_ZsPl 
         N_Egest_Zs=N_Egest_ZsPs+N_Egest_ZsPl 
         zooS=zooS+N_Assim_Zs
     
         N_Assim_Zl=N_Assim_ZlPs+N_Assim_ZlPl+N_Assim_ZlZs+N_Assim_ZlDs
         N_Egest_Zl=N_Egest_ZlPs+N_Egest_ZlPl+N_Egest_ZlZs+N_Egest_ZlDs
              
         zooL=zooL+N_Assim_Zl
         SDeN=SDeN+N_Egest_Zs+N_Egest_Zl
     
         N_Assim = N_Assim_Zs+N_Assim_Zl     
         N_Egest = N_Egest_Zs+N_Egest_Zl
         
         Opal = Opal + (N_Graz_ZsPl+N_Graz_ZlPl)*self%DiatSiN

!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
         if (self%TEMP_RATES) then
            fac1=d_per_s*self%PhySMR*Tcoef0*(Tcoef1**temp)
            fac2=d_per_s*self%PhyLMR*Tcoef0*(Tcoef1**temp)
         else
            fac1=d_per_s*self%PhySMR
            fac2=d_per_s*self%PhyLMR
         end if

!
! Small phytoplankton mortality (1).
!
         N_Pmortal_Ps=fac1*MAX(PhyS-self%PhyMin,0.0_rk)
         phyS=phyS-N_Pmortal_Ps
              
         Chl_Pmortal_Ps=fac1*MAX(chlS-self%ChlMin,0.0_rk)
         chlS=chlS-Chl_Pmortal_Ps
              
!
! Large phytoplankton mortality (2).
!
         N_Pmortal_Pl=fac2*MAX(PhyL-self%PhyMin,0.0_rk)
         phyL=phyL-N_Pmortal_Pl
              
         Chl_Pmortal_Pl=fac2*MAX(chlL-self%ChlMin,0.0_rk)
         chlL=chlL-Chl_Pmortal_Pl

!
! sum up all phytoplankton mortality fluxes.
!
         N_Pmortal=N_Pmortal_Ps+N_Pmortal_Pl
         Chl_Pmortal=Chl_Pmortal_Ps+Chl_Pmortal_Pl
         
         SDeN=SDeN+N_Pmortal
         
         Opal=Opal+N_Pmortal_Pl*self%DiatSiN
         
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
         if (self%TEMP_RATES) then
            fac2=d_per_s*self%ZooSBM*Tcoef0*(Tcoef1**temp) !metabolism
            fac3=d_per_s*self%ZooSMR*Tcoef0*(Tcoef1**temp) !mortality
            fac4=d_per_s*self%ZooSER*Tcoef0*(Tcoef1**temp) !excretion
          
            fac5=d_per_s*self%ZooLBM*Tcoef0*(Tcoef1**temp) !metabolism
            fac6=d_per_s*self%ZooLMR*Tcoef0*(Tcoef1**temp) !mortality
            fac7=d_per_s*self%ZooLER*Tcoef0*(Tcoef1**temp) !excretion
         else
            fac2=d_per_s*self%ZooSBM !metabolism
            fac3=d_per_s*self%ZooSMR !mortality
            fac4=d_per_s*self%ZooSER !excretion
          
            fac5=d_per_s*self%ZooLBM !metabolism
            fac6=d_per_s*self%ZooLMR !mortality
            fac7=d_per_s*self%ZooLER !excretion
         end if     

!
! Small zooplankton mortality 
!
         cff2=fac3*zooS
         zooS=zooS/(1.0_rk+cff2) 
         N_Zsmortal=cff2*zooS
         SDeN=SDeN+N_Zsmortal
         
!
! Small zooplankton excretion of small and large phytoplankton
!
         fac1=fac4*phyS*phyS/(self%K_ZsPs+phyS*phyS)
         ppref = EXP(-self%ZooS_GrInPl*phyS)
         fac3=ppref*fac4*phyL*phyL/(self%K_ZsPl+phyL*phyL)
         cff3=(fac1+fac3)*self%ZooSAE_N
         zooS=zooS/(1.0_rk+cff3)
         N_Zsexcret=cff3*zooS 

!  Small zooplankton basal metabolism (limited by a zooplankton minimum).
!
         N_Zsmetabo=fac2*MAX(zooS-self%zooMin,0.0_rk)
         zooS=zooS-N_Zsmetabo
         nh4=nh4+N_Zsmetabo 

!
! Large zooplankton mortality
!
         cff2=fac6*zooL
         zooL=zooL/(1.0_rk+cff2)
         N_Zlmortal=cff2*zooL
         LDeN=LDeN+N_Zlmortal
         
         Calcite=Calcite+N_Zlmortal*fcaco3
         
!
! Large zooplankton excretion of small and large phytoplankton and small 
! zooplankton and small detritus
!
         ppref = EXP(-self%ZooL_GrInPs*(phyL+zooS+SDeN*switch_ZlDs))
         fac1=ppref*fac7*phyS*phyS/(self%K_ZlPs+phyS*phyS)
         fac2=fac7*phyL*phyL/(self%K_ZlPl+phyL*phyL)
         fac3=fac7*zooS*zooS/(self%K_ZlZs+zooS*zooS)
         ppref = EXP(-self%ZooL_GrInDs*(phyS+phyL+zooS))
         fac4=fac7*SDeN*SDeN/(self%K_ZlDs+SDeN*SDeN)
         cff3=(fac1+fac2+fac3+fac4)*self%ZooLAE_N
         zooL=zooL/(1.0_rk+cff3)
         N_Zlexcret=cff3*zooL

! Large zooplankton basal metabolism (limited by a zooplankton minimum).
!
         N_Zlmetabo=fac5*MAX(zooL-self%ZooMin,0.0_rk)
         zooL=zooL-N_Zlmetabo
         nh4=nh4+N_Zlmetabo

!
! sum up all fluxes.
!
         N_Zmortal= N_Zsmortal+N_Zlmortal
         N_Zmetabo= N_Zsmetabo+N_Zlmetabo
         N_Zexcret= N_Zsexcret+N_Zlexcret
         nh4=nh4+N_Zexcret

         if (self%oxy_model) then
            oxy=oxy-rOxNH4*(N_Zmetabo+N_Zexcret)
         end if
!         
!-----------------------------------------------------------------------
!  Coagulation of large phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
         fac1=d_per_s*self%CoagR
         cff1=fac1*(SDeN+phyL)
         cff2=1.0_rk/(1.0_rk+cff1)
         phyL=phyL*cff2
         chlL=chlL*cff2
         SDeN=SDeN*cff2
         
         N_CoagP=phyL*cff1
         Chl_Coag=ChlL*cff1
         N_CoagD=SDeN*cff1
         
         LDeN=LDeN+N_CoagP+N_CoagD
         Opal=Opal+N_CoagP*self%DiatSiN
         Calcite=Calcite+(N_CoagP+N_CoagD)*fcaco3
         
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!

! only free part of LDeN is subject to remineralization and fragmentation
         pLDeN = MIN(Opal*self%OpalPR+Calcite*self%CalcPR,LDeN)
         fLDeN = MAX(LDeN-pLDeN,0.0_rk)
         
         zfact = depth/(depth+self%Kz)*(zref+self%Kz)/zref
         zfact = MIN(zfact,1.0_rk)
         
         if (self%TEMP_RATES) then
            fac3=d_per_s*self%SDeRRN*Tcoef0*(Tcoef1**temp) 
            fac4=d_per_s*self%LDeRRN*Tcoef0*(Tcoef1**temp)*zfact 
         else
            fac3=d_per_s*self%SDeRRN 
            fac4=d_per_s*self%LDeRRN*zfact  
         end if
         
         
         if (self%oxy_model) then
            fac1=MAX(oxy-6.0_rk,0.0_rk) ! O2 off max
            fac2=MAX(fac1/(3.0_rk+fac1),0.0_rk) ! MM for O2 dependence
            
            cff1=fac3*fac2
            cff2=1.0_rk/(1.0_rk+cff1)
            cff3=fac4*fac2
            cff4=1.0_rk/(1.0_rk+cff3)
            
            SDeN=SDeN*cff2
            fLDeN=fLDeN*cff4
            
            N_RemineS=SDeN*cff1
            N_RemineL=fLDeN*cff3
            
            nh4=nh4+N_RemineS+N_RemineL
            oxy=oxy-(N_RemineS+N_RemineL)*rOxNH4
         else
            cff1=fac3
            cff2=1.0_rk/(1.0_rk+cff1)
            cff3=fac4
            cff4=1.0_rk/(1.0_rk+cff3)
         
            SDeN=SDeN*cff2
            fLDeN=fLDeN*cff4
         
            N_RemineS=SDeN*cff1
            N_RemineL=fLDeN*cff3
         
            nh4=nh4+N_RemineS+N_RemineL
         end if

!-----------------------------------------------------------------------
!  Dissolution of opal into Si(OH)4 (note: we didn't simulate Si(OH)4)
!-----------------------------------------------------------------------
!
         fac3=d_per_s*self%OpalDR 
         cff1=fac3
         cff2=1.0_rk/(1.0_rk+cff1)
         Opal=Opal*cff2
         Opal_Dissolve=Opal*cff1

!-----------------------------------------------------------------------
!  Dissolution of CaCO3
!-----------------------------------------------------------------------
!
         fac3=d_per_s*self%CalcDR 
         cff1=fac3
         cff2=1.0_rk/(1.0_rk+cff1)
         Calcite=Calcite*cff2
         Calc_Dissolve=Calcite*cff1
         
!-----------------------------------------------------------------------
!  Large Detritus disaggregate into small Detritus, fragmentation
!-----------------------------------------------------------------------
!

         fac3=d_per_s*self%FragRN 
         cff1=fac3
         cff2=1.0_rk/(1.0_rk+cff1)
         fLDeN=fLDeN*cff2
         N_Fragment=fLDeN*cff1
         SDeN=SDeN+N_Fragment

!-----------------------------------------------------------------------
!  The true large Detritus (sum of free and protected part)
!-----------------------------------------------------------------------
!
         LDeN = fLDeN+pLDeN
!
!-----------------------------------------------------------------------
!  Sum all contributor to change rate:
!-----------------------------------------------------------------------
!
         d_no3 = -N_NewProd+N_Nitrifi
         d_nh4 = -N_RegProd-N_Nitrifi+N_Zexcret+N_Zmetabo+N_RemineS+N_RemineL
         
         d_phyS = N_NewProd_Ps+N_RegProd_Ps-N_Graz_ZsPs-N_Graz_ZlPs-N_Pmortal_Ps
         d_phyL = N_NewProd_Pl+N_RegProd_Pl-N_Graz_ZsPl-N_Graz_ZlPl-N_Pmortal_Pl-N_CoagP
         
         d_chlS = Chl_Prod_Ps-Chl_Graz_ZsPs-Chl_Graz_ZlPs-Chl_Pmortal_Ps
         d_chlL = Chl_Prod_Pl-Chl_Graz_ZsPl-Chl_Graz_ZlPl-Chl_Pmortal_Pl-Chl_Coag
         
         d_zooS = -N_Assim_ZlZs-N_Egest_ZlZs+N_Assim_Zs-N_Zsmortal-N_Zsexcret-N_Zsmetabo
         d_zooL = N_Assim_Zl-N_Zlmortal-N_Zlexcret-N_Zlmetabo
         
         d_SDeN = -N_Assim_ZlDs-N_Egest_ZlDs+N_Egest_Zs+N_Egest_Zl &
      &           +N_Pmortal+N_Zsmortal-N_RemineS-N_CoagD+N_Fragment
         d_LDeN = N_Zlmortal+N_CoagP+N_CoagD-N_RemineL-N_Fragment
         
         d_Opal = (N_Graz_ZsPl+N_Graz_ZlPl+N_Pmortal_Pl+N_CoagP)*self%DiatSiN-Opal_Dissolve
         d_Calc = (N_Zlmortal+N_CoagP+N_CoagD)*fcaco3-Calc_Dissolve
         
         d_pre_no3 = -N_NewProd
         d_reg_no3 = N_Nitrifi
         
         ! Send rates of change to FABM.
         _SET_ODE_(self%id_no3,d_no3)
         _SET_ODE_(self%id_nh4,d_nh4)
         
         _SET_ODE_(self%id_phyS,d_phyS)
         _SET_ODE_(self%id_phyL,d_phyL)
         
         _SET_ODE_(self%id_chlS,d_chlS)
         _SET_ODE_(self%id_chlL,d_chlL)
         
         _SET_ODE_(self%id_zooS,d_zooS)
         _SET_ODE_(self%id_zooL,d_zooL)
         
         _SET_ODE_(self%id_LDeN,d_LDeN)
         _SET_ODE_(self%id_SDeN,d_SDeN)
         
         _SET_ODE_(self%id_Opal,d_Opal)
         _SET_ODE_(self%id_Calc,d_Calc)
         
         _SET_ODE_(self%id_pre_no3,d_pre_no3)
         _SET_ODE_(self%id_reg_no3,d_reg_no3)
         
         ! Provide diagnostic variables to FABM.
         ! those fluxes will be useful in the oxygen model
         _SET_DIAGNOSTIC_(self%id_newPP,secs_pr_day*N_NewProd)
         _SET_DIAGNOSTIC_(self%id_regPP,secs_pr_day*N_RegProd)
         _SET_DIAGNOSTIC_(self%id_Zmetabo,secs_pr_day*N_Zmetabo)
         _SET_DIAGNOSTIC_(self%id_Zexcret,secs_pr_day*N_Zexcret)
         _SET_DIAGNOSTIC_(self%id_remineS,secs_pr_day*N_RemineS)
         _SET_DIAGNOSTIC_(self%id_remineL,secs_pr_day*N_RemineL)
         _SET_DIAGNOSTIC_(self%id_fragment,secs_pr_day*N_Fragment)
         _SET_DIAGNOSTIC_(self%id_nitrifi,secs_pr_day*N_Nitrifi)
         
         _SET_DIAGNOSTIC_(self%id_newPP_Ps,secs_pr_day*N_NewProd_Ps)
         _SET_DIAGNOSTIC_(self%id_newPP_Pl,secs_pr_day*N_NewProd_Pl)
         
         _SET_DIAGNOSTIC_(self%id_fcaco3,fcaco3/self%PhyCN)
                  
      _LOOP_END_

   
   END SUBROUTINE do
   
   SUBROUTINE do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_memg_bio_fennel_2p2z),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
      real(rk) :: phyS,phyL,chlS,chlL,zooS,zooL,LDeN,SDeN
      real(rk) :: Opal,Calcite
      real(rk) :: decay_no3
      real(rk) :: decay_phyS,decay_phyL
      real(rk) :: decay_zooS,decay_zooL
      real(rk) :: decay_chlS,decay_chlL
      real(rk) :: decay_LDeN,decay_SDeN
      real(rk) :: decay_Opal,decay_Calc
      real(rk) :: decay_reg_no3
      
      real(rk) :: fac,cff1,cff2
      
      real(rk), parameter :: secs_pr_day = 86400.0_rk
      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      ! we assume that all organic N remineralized quickly (1 hour) at bottom
      real(rk), parameter :: decay_rate = 24.0_rk

      _HORIZONTAL_LOOP_BEGIN_
      
      ! Obtain concentration of biological variables.
         _GET_(self%id_phyS,phyS)
         _GET_(self%id_phyL,phyL)
         
         _GET_(self%id_chlS,chlS)
         _GET_(self%id_chlL,chlL)
         
         _GET_(self%id_zooS,zooS)
         _GET_(self%id_zooL,zooL)
         
         _GET_(self%id_LDeN,LDeN)
         _GET_(self%id_SDeN,SDeN)
         
         _GET_(self%id_Opal,Opal)
         _GET_(self%id_Calc,Calcite)
         
         fac=d_per_s*decay_rate
!
! all organic N decay into NH4 quickly.
!
         decay_phyS=-fac*MAX(PhyS-self%PhyMin,0.0_rk)
         decay_chlS=-fac*MAX(chlS-self%ChlMin,0.0_rk)
         decay_phyL=-fac*MAX(PhyL-self%PhyMin,0.0_rk)
         decay_chlL=-fac*MAX(chlL-self%ChlMin,0.0_rk)
         
         decay_zooS=-fac*MAX(zooS-self%PhyMin,0.0_rk)
         decay_zooL=-fac*MAX(zooL-self%PhyMin,0.0_rk)
         decay_SDeN=-fac*MAX(SDeN-self%PhyMin,0.0_rk)
         decay_LDeN=-fac*MAX(LDeN-self%PhyMin,0.0_rk)
         
         decay_no3 = -(decay_phyS+decay_phyL+decay_zooS+decay_zooL+decay_SDeN+decay_LDeN)
         
         decay_Opal=-fac*MAX(Opal-self%PhyMin,0.0_rk)
         decay_Calc=-fac*MAX(Calcite-self%PhyMin,0.0_rk)
         
         ! Send rates of change to FABM.
         _SET_BOTTOM_EXCHANGE_(self%id_no3,decay_no3)
         
         _SET_BOTTOM_EXCHANGE_(self%id_phyS,decay_phyS)
         _SET_BOTTOM_EXCHANGE_(self%id_phyL,decay_phyL)
         
         _SET_BOTTOM_EXCHANGE_(self%id_chlS,decay_chlS)
         _SET_BOTTOM_EXCHANGE_(self%id_chlL,decay_chlL)
         
         _SET_BOTTOM_EXCHANGE_(self%id_zooS,decay_zooS)
         _SET_BOTTOM_EXCHANGE_(self%id_zooL,decay_zooL)
         
         _SET_BOTTOM_EXCHANGE_(self%id_LDeN,decay_LDeN)
         _SET_BOTTOM_EXCHANGE_(self%id_SDeN,decay_SDeN)
         
         _SET_BOTTOM_EXCHANGE_(self%id_Opal,decay_Opal)
         _SET_BOTTOM_EXCHANGE_(self%id_Calc,decay_Calc)
         
         _SET_BOTTOM_EXCHANGE_(self%id_reg_no3,decay_reg_no3)
         
      _HORIZONTAL_LOOP_END_

   END SUBROUTINE do_bottom

END MODULE
