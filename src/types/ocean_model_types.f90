MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_ocean_matrix_interpolation

    ! Time fields for interpolation
    REAL(dp) :: t0                                                         ! Start time for interpolation
    REAL(dp) :: t1                                                         ! End time for interpolation     

    ! Temperature and salinity arrays at each timeframe
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: T0  ! Temperature at t0
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: S0  ! Salinity at t0
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: T1  ! Temperature at t1
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: S1  ! Salinity at t1    

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