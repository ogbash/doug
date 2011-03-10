!> Parameters for the modules of the components.
module mparameters

  ! Distribution component

  !> Use elemental input module as the distribution component
  integer, parameter :: DISTRIBUTION_TYPE_ELEMENTAL = 1
  !> Use assembled input module as the distribution component
  integer, parameter :: DISTRIBUTION_TYPE_ASSEMBLED = 2
  !> Use structured mesh generator module as the distribution component
  integer, parameter :: DISTRIBUTION_TYPE_STRUCTURED = 3

end module mparameters
