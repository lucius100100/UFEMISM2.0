MODULE climate_matrix

    ! Timeframe0 = PI
    ! Timeframe1 = LGM

    ! Matrix climate models
  
  ! ====================
  ! ===== Preamble =====
  ! ====================
  
    USE precisions                                             , ONLY: dp
    USE mpi_basic                                              , ONLY: par, sync
    USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
    USE model_configuration                                    , ONLY: C
    USE parameters
    USE mesh_types                                             , ONLY: type_mesh
    USE ice_model_types                                        , ONLY: type_ice_model
    USE climate_model_types                                    , ONLY: type_climate_model
    USE region_types                                           , ONLY: type_model_region
    USE netcdf_input                                           , ONLY: read_field_from_file_2D_monthly
    USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
    USE reference_geometry_types                               , ONLY: type_reference_geometry
    USE grid_types                                             , ONLY: type_grid
    USE mesh_data_smoothing                                    , ONLY: smooth_Gaussian_2D

    IMPLICIT NONE
  
  CONTAINS
  
  ! =========================
  ! ===== Main routines =====
  ! =========================

  ! TO DO, needs more elaborate implementation of matrix method as shown in for example Scherrenberg et al., 2024

    SUBROUTINE linear_time_interpolation(mesh, climate, time)
      ! Linear interpolation between two climate snapshots

      IMPLICIT NONE

      ! In/output variables:
      TYPE(type_mesh),                        INTENT(IN)    :: mesh
      TYPE(type_climate_model),               INTENT(INOUT) :: climate
      REAL(dp),                               INTENT(IN)    :: time

      ! Local variables:
      CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'linear_time_interpolation'
      REAL(dp)                                              :: wt0
      INTEGER                                               :: vi, m

      ! Add routine to path
      CALL init_routine( routine_name)

      ! Ensure matrix timeframes are allocated
          IF (.NOT. ALLOCATED(climate%matrix%timeframe0%T2m) .OR. &
          .NOT. ALLOCATED(climate%matrix%timeframe1%T2m) .OR. &
          .NOT. ALLOCATED(climate%matrix%timeframe0%Precip) .OR. &
          .NOT. ALLOCATED(climate%matrix%timeframe1%Precip)) THEN
        CALL crash('Climate matrix timeframes are not allocated.')
      END IF

      ! Check for division by 0 error
      IF (ABS(climate%matrix%t1 - climate%matrix%t0) < 1e-5_dp) THEN
        CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
      END IF

      ! Calculate weights for linear interpolation
      wt0 = (time - climate%matrix%t1) / (climate%matrix%t0 - climate%matrix%t1)

      IF (par%master) THEN
        print *, "Interpolation weight climate =", wt0
      END IF

      ! Apply linear interpolation
      DO vi = mesh%vi1, mesh%vi2
          DO m = 1, 12
            climate%T2m(vi,m) = wt0 * climate%matrix%timeframe0%T2m(vi,m) + (1.0_dp - wt0) * climate%matrix%timeframe1%T2m(vi,m)
            climate%Precip(vi,m) = wt0 * climate%matrix%timeframe0%Precip(vi,m) + (1.0_dp - wt0) * climate%matrix%timeframe1%Precip(vi,m)
          END DO
      END DO

      ! Finalise routine path
      CALL finalise_routine( routine_name)

    END SUBROUTINE linear_time_interpolation

    SUBROUTINE run_climate_model_matrix( mesh, ice, climate, time, region_name)
      ! Calculate the climate using an interpolating matrix climate scheme
    
      IMPLICIT NONE
    
      ! In/output variables:
      TYPE(type_mesh),                        INTENT(IN)    :: mesh
      TYPE(type_ice_model),                   INTENT(IN)    :: ice
      TYPE(type_climate_model),               INTENT(INOUT) :: climate
      CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
      REAL(dp),                               INTENT(IN)    :: time
    
      ! Local variables:
      CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_matrix'

      ! Add routine to path
      CALL init_routine( routine_name) 

      ! Perform time interpolation
      SELECT CASE (TRIM(C%choice_climate_model_matrix))
      CASE('linear_time')
        CALL linear_time_interpolation(mesh, climate, time)
      CASE DEFAULT
        CALL crash('Unknown choice_climate_model_matrix' // TRIM(C%choice_climate_model_matrix))
      END SELECT

      ! Finalise routine path
      CALL finalise_routine( routine_name)
    
    END SUBROUTINE run_climate_model_matrix
    
    SUBROUTINE initialise_climate_model_matrix( mesh, climate, region_name)
      ! Initialise the climate matrix model
    
      IMPLICIT NONE
    
      ! In/output variables:
      TYPE(type_mesh),                        INTENT(IN)    :: mesh
      TYPE(type_climate_model),               INTENT(INOUT) :: climate
      CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    
      ! Local variables:
      CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_matrix'
      CHARACTER(LEN=256)                                    :: filename1, filename2, filename_insolation, filename_tas, filename_sos
      INTEGER                                               :: i, j, vi, m

      ! Add routine to path
      CALL init_routine( routine_name)
    
      ! Print to terminal
      IF (par%master)  WRITE(*,"(A)") '     Initialising matrix climate model "' // &
        colour_string( TRIM( C%choice_climate_model_matrix),'light blue') // '"...'

      ! Start and ending of simulation
      climate%matrix%t0 = REAL(C%start_time_of_run, dp)             ! PI
      climate%matrix%t1 = REAL(C%end_time_of_run, dp)               ! LGM

      ! Allocate memory for climate timeframes T2m and Precip array if not already allocated
      ! timeframe0%T2m
      IF (.NOT. ALLOCATED(climate%matrix%timeframe0%T2m)) THEN
        ALLOCATE(climate%matrix%timeframe0%T2m(mesh%vi1:mesh%vi2, 12))
        climate%matrix%timeframe0%T2m = 0._dp
      END IF
    
      ! timeframe0%Precip
      IF (.NOT. ALLOCATED(climate%matrix%timeframe0%Precip)) THEN
        ALLOCATE(climate%matrix%timeframe0%Precip(mesh%vi1:mesh%vi2, 12))
        climate%matrix%timeframe0%Precip = 0._dp
      END IF
      
      ! timeframe1%T2m
      IF (.NOT. ALLOCATED(climate%matrix%timeframe1%T2m)) THEN
        ALLOCATE(climate%matrix%timeframe1%T2m(mesh%vi1:mesh%vi2, 12))
        climate%matrix%timeframe0%T2m = 0._dp
      END IF
    
      ! timeframe1%Precip
      IF (.NOT. ALLOCATED(climate%matrix%timeframe1%Precip)) THEN
        ALLOCATE(climate%matrix%timeframe1%Precip(mesh%vi1:mesh%vi2, 12))
        climate%matrix%timeframe0%Precip = 0._dp
      END IF

      ! Construct filenames
      filename1            = TRIM(C%filename_climate_matrix_base1)  ! Snapshot PI
      filename2            = TRIM(C%filename_climate_matrix_base2)  ! Snapshot LGM

      ! Read PI (timeframe0) from netCDF
      CALL read_field_from_file_2D_monthly(filename1, 'T2m',   mesh, climate%matrix%timeframe0%T2m)
      CALL read_field_from_file_2D_monthly(filename1, 'Precip',mesh, climate%matrix%timeframe0%Precip)
    
      ! Read LGM (timeframe1) from netCDF
      CALL read_field_from_file_2D_monthly(filename2, 'T2m',   mesh, climate%matrix%timeframe1%T2m)
      CALL read_field_from_file_2D_monthly(filename2, 'Precip',mesh, climate%matrix%timeframe1%Precip)
    
      ! Ensure correct model choice and load in files if necessary
      SELECT CASE (TRIM(C%choice_climate_model_matrix))
      CASE('linear_time')
        ! Global weight variables or time, no specific intialisation necessary

      CASE DEFAULT
        CALL crash('Unknown choice_climate_model_matrix' // TRIM(C%choice_climate_model_matrix))
      END SELECT

      ! Prescribe initial climate state (start from PI snapshot)
      DO vi = mesh%vi1, mesh%vi2
        DO m = 1, 12
          climate%T2m(vi,m)    = climate%matrix%timeframe0%T2m(vi,m)
          climate%Precip(vi,m) = climate%matrix%timeframe0%Precip(vi,m)
        END DO
      END DO

      ! Print snapshot info
      IF (par%master) THEN
        WRITE(*,*) 'Snapshot timeframe0 T2m min/max:', MINVAL(climate%matrix%timeframe0%T2m), &
                                                  MAXVAL(climate%matrix%timeframe0%T2m)
        WRITE(*,*) 'Snapshot timeframe1 T2m min/max:', MINVAL(climate%matrix%timeframe1%T2m), &
                                                  MAXVAL(climate%matrix%timeframe1%T2m)
        WRITE(*,*) 'Snapshot timeframe0 Precip min/max:', MINVAL(climate%matrix%timeframe0%Precip), &
                                                     MAXVAL(climate%matrix%timeframe0%Precip)
        WRITE(*,*) 'Snapshot timeframe1 Precip min/max:', MINVAL(climate%matrix%timeframe1%Precip), &
                                                     MAXVAL(climate%matrix%timeframe1%Precip)
      END IF

      ! Finalise routine path
      CALL finalise_routine( routine_name)
    
    END SUBROUTINE initialise_climate_model_matrix

END MODULE climate_matrix