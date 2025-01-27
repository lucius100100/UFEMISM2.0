MODULE climate_model_types

  ! The different data types used in the climate modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_timeframe
    
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T2m                         ! [K]      Monthly 2-m air temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Precip                      ! [m.w.e.] Monthly precipitation

  END TYPE type_timeframe

  TYPE type_climate_matrix_interpolation

    ! Time fields for interpolation
    REAL(dp)                                :: t0                          ! Start time for interpolation
    REAL(dp)                                :: t1                          ! End time for interpolation     

    ! Timeframes containing temperature and salinity data
    TYPE(type_timeframe)                    :: timeframe0                  ! PI
    TYPE(type_timeframe)                    :: timeframe1                  ! LGM     

  END TYPE type_climate_matrix_interpolation

  TYPE type_climate_model
    ! The climate model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T2m                         ! [K]      Monthly 2-m air temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: Precip                      ! [m.w.e.] Monthly precipitation

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

    ! Matrix component
    TYPE(type_climate_matrix_interpolation)   :: matrix

  END TYPE type_climate_model

CONTAINS

END MODULE climate_model_types