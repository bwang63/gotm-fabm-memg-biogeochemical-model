module memg_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   ! Add use statements for new models here
   use memg_bio_fennel_model_library

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: initialize      
      procedure :: create
   end type

   type (type_factory),save,target,public :: memg_model_factory

contains

   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
      ! Add additional child model factories here
      call self%add(memg_bio_fennel_model_factory,'bio_fennel')

      ! Go through default initializaton steps.
      ! This also allows newly added child model factories to initialize.
      call self%type_base_model_factory%initialize()
   end subroutine initialize
   
   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         ! Add case statements for new models here
         case default
            call self%type_base_model_factory%create(name,model)

      end select

   end subroutine create

end module