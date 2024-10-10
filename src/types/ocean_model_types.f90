MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_timeframe
  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: T                            ! Temperature array for this timeframe
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: S                            ! Salinity array for this timeframe
  
  END TYPE type_timeframe

  TYPE type_ocean_matrix_interpolation

    ! Time fields for interpolation
    REAL(dp) :: t0                                                         ! Start time for interpolation
    REAL(dp) :: t1                                                         ! End time for interpolation     

    ! Timeframes containing temperature and salinity data
    TYPE(type_timeframe) :: timeframe0                                     ! LGM
    TYPE(type_timeframe) :: timeframe1                                     ! PI                      

  END TYPE type_ocean_matrix_interpolation

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

  END TYPE type_ocean_model

CONTAINS

END MODULE ocean_model_types