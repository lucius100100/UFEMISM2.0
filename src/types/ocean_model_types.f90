MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE :: type_timeframe
  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: T                            ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: S                            ! [PSU]             Salinity
  
  END TYPE type_timeframe

  TYPE type_ocean_matrix

    ! Time fields for interpolation
    REAL(dp) :: t0                                                         ! Start time for interpolation
    REAL(dp) :: t1                                                         ! End time for interpolation     

    ! Timeframes for interpolation
    TYPE(type_timeframe), ALLOCATABLE :: timeframe0                        ! Ocean state at t0
    TYPE(type_timeframe), ALLOCATABLE :: timeframe1                        ! Ocean state at t1
    
  END TYPE type_ocean_matrix

  TYPE type_ocean_model
    ! The ocean model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S                           ! [PSU]             Salinity

    ! Secondary data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_draft                     ! [degrees Celsius] Temperature at ice base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_freezing_point            ! [degrees Celsius] Pressure freezing point of water

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

    ! Matrix ocean model
    TYPE(type_ocean_matrix)                 :: matrix                      ! Data structure to handle multiple timeframes for interpolation

  END TYPE type_ocean_model

CONTAINS

END MODULE ocean_model_types