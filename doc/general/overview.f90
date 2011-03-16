!> \page p_overview DOUG overview
!!
!! DOUG has currently has 2 executables: \c doug_geom and \c doug_aggr, which can read different input formats and use
!!  different preconditioners:
!! - the \c doug_geom executable can read elemental input and use geometric coarse preconditioner,
!! - the \c doug_aggr executable can read assembled input, generate mesh locally, use aggregation and aggregation based coarse preconditioners.
!!
!! They occasionally need to be combined, but it requires more refactoring.
!!
