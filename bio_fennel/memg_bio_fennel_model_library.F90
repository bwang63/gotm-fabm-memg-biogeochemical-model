module memg_bio_fennel_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   ! Add new bio_fennel modules here
   use memg_bio_fennel_1p1z
   use memg_bio_fennel_2p2z
   use memg_bio_fennel_2p2z_sink
   use memg_bio_fennel_light
   use memg_bio_fennel_oxygen

   
   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: memg_bio_fennel_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         ! Add new bio_fennel models here
         case ('1p1z'); allocate(type_memg_bio_fennel_1p1z::model)
         case ('2p2z'); allocate(type_memg_bio_fennel_2p2z::model)
         case ('2p2z_sink'); allocate(type_memg_bio_fennel_2p2z_sink::model)
         case ('light'); allocate(type_memg_bio_fennel_light::model)
         case ('oxygen'); allocate(type_memg_bio_fennel_oxygen::model)
      end select
   end subroutine create

end module memg_bio_fennel_model_library
