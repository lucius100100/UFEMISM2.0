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
    REAL(dp) :: t2      ! In case of more than two snapshots
    REAL(dp) :: t3      
    REAL(dp) :: t4      
    REAL(dp) :: t5      
    REAL(dp) :: t6      
    REAL(dp) :: t7      
    REAL(dp) :: t8      
    REAL(dp) :: t9      
    REAL(dp) :: t10      

    ! Timeframes for interpolation
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe0  ! Ocean state at t0
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe1  ! Ocean state at t1
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe2  ! Ocean state at t2
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe3  ! Ocean state at t3
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe4  ! Ocean state at t4
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe5  ! Ocean state at t5
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe6  ! Ocean state at t6
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe7  ! Ocean state at t7
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe8  ! Ocean state at t8
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe9  ! Ocean state at t9
    TYPE(type_ocean_matrix), ALLOCATABLE :: timeframe10  ! Ocean state at t10

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