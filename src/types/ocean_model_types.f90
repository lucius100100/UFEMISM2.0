MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  type type_ocean_matrix

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S                           ! [PSU]             Salinity

    ! Time fields for interpolation
    REAL(dp) :: t0      ! Start time for interpolation
    REAL(dp) :: t1      ! End time for interpolation

    ! Timeframes for interpolation
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe0  ! Ocean state at t0
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe1  ! Ocean state at t1

  end type type_ocean_matrix

  TYPE type_ocean_model
    ! The ocean model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S                           ! [PSU]             Salinity

    ! Secondary data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_draft                     ! [degrees Celsius] Temperature at ice base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_freezing_point            ! [degrees Celsius] Pressure freezing point of water

    type(type_ocean_matrix) :: matrix

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_ocean_model

CONTAINS

END MODULE ocean_model_types