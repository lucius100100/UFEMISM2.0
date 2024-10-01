MODULE ocean_matrix

    ! Matrix ocean models
  
  ! ===== Preamble =====
  ! ====================
  
    USE precisions                                             , ONLY: dp
    USE mpi_basic                                              , ONLY: par, sync
    USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
    USE model_configuration                                    , ONLY: C
    USE parameters
    USE mesh_types                                             , ONLY: type_mesh
    USE ice_model_types                                        , ONLY: type_ice_model
    USE ocean_model_types                                      , ONLY: type_ocean_model
    USE netcdf_input                                           , ONLY: read_field_from_file_3D_ocean
    USE netcdf_basic                                           , ONLY: field_name_options_T_ocean, field_name_options_S_ocean
  
    IMPLICIT NONE
  
  CONTAINS
  
  ! ===== Main routines =====
  ! =========================

  SUBROUTINE count_allocated_timeframes(ocean, num_timeframes)
    ! Count number of timeframes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ocean_model),               INTENT(IN)    :: ocean
    INTEGER,                              INTENT(OUT)   :: num_timeframes

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'count_allocated_timeframes'
    INTEGER, PARAMETER                                  :: MAX_TIMEFRAMES = 11

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Count 
    num_timeframes = 0
    IF (ALLOCATED(ocean%matrix%timeframe0%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe1%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe2%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe3%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe4%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe5%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe6%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe7%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe8%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe9%T)) num_timeframes = num_timeframes + 1
    IF (ALLOCATED(ocean%matrix%timeframe10%T)) num_timeframes = num_timeframes + 1
    ! Can always add more if needed, for now capped at 10 timeframes

    ! Don't exceed amount of supported number of timeframes
    IF (num_timeframes > MAX_TIMEFRAMES) THEN
      CALL crash('Number of timeframes exceeds the coded limit.')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE count_allocated_timeframes

  SUBROUTINE fill_times_array(ocean, times)
    ! Fill times array for amount of timeframes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ocean_model),               INTENT(IN)    :: ocean
    REAL(dp), DIMENSION(:),               INTENT(OUT)   :: times
      
    ! Local variables:
    INTEGER                                             :: n
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'fill_times_array'

    ! Add routine to path
    CALL init_routine( routine_name)

    n = 0
    IF (ALLOCATED(ocean%matrix%timeframe0%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t0
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe1%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t1
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe2%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t2
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe3%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t3
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe4%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t4
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe5%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t5
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe6%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t6
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe7%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t7
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe8%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t8
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe9%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t9
    END IF
    IF (ALLOCATED(ocean%matrix%timeframe10%T)) THEN
        n = n + 1
        times(n) = ocean%matrix%t10
    END IF
    ! Can always add more if needed, for now capped at 10 timeframes

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE fill_times_array

  SUBROUTINE update_ocean_matrix_timeframes(mesh, ocean, region_name, time)
    ! Update the ocean matrix timeframes
    !
    ! Use an interpolating matrix ocean scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_ocean_matrix_timeframes'
    INTEGER                                               :: year0, year1
    CHARACTER(LEN=256)                                    :: filename1, filename2
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine two closest years to the current time
    year0 = FLOOR(time)
    year1 = CEILING(time)

    ! Set time for the two snapshots
    ocean%matrix%t0 = REAL(year0, dp)
    ocean%matrix%t1 = REAL(year1, dp)

    ! Allocate memory for all possible timeframes
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe0)) THEN
      ALLOCATE(ocean%matrix%timeframe0)
    END IF
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe1)) THEN
        ALLOCATE(ocean%matrix%timeframe1)
    END IF
    ! Can always add more if needed, for now capped at 10 timeframes

    ! Construct filenames for the two ocean snapshots
    filename1 = TRIM(C%filename_ocean_matrix_base1) // TRIM(region_name) // '_' // TRIM(ADJUSTL(CHAR(year0))) // '.nc'
    filename2 = TRIM(C%filename_ocean_matrix_base2) // TRIM(region_name) // '_' // TRIM(ADJUSTL(CHAR(year1))) // '.nc'

    ! Read the ocean snapshots
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_T_ocean, mesh, ocean%matrix%timeframe0%T)
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_S_ocean, mesh, ocean%matrix%timeframe0%S)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_T_ocean, mesh, ocean%matrix%timeframe1%T)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_S_ocean, mesh, ocean%matrix%timeframe1%S)
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_ocean_matrix_timeframes

  SUBROUTINE linear_interpolation(mesh, ocean, time)
    ! Linear interpolation 
    !
    ! FIX linear interpolation for more than 2 timeframes

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'linear_interpolation'
    REAL(dp)                                              :: wt0, wt1
    INTEGER                                               :: i, j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate weights for linear interpolation
    wt0 = (ocean%matrix%t1 - time) / (ocean%matrix%t1 - ocean%matrix%t0)
    wt1 = 1.0_dp - wt0

    ! Apply linear interpolation
    DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
            ocean%T(i,j) = wt0 * ocean%matrix%timeframe0%T(i,j) + wt1 * ocean%matrix%timeframe1%T(i,j)
            ocean%S(i,j) = wt0 * ocean%matrix%timeframe0%S(i,j) + wt1 * ocean%matrix%timeframe1%S(i,j)
        END DO
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE linear_interpolation

  SUBROUTINE polynomial_interpolation(mesh, ocean, time, num_timeframes, times, weights)
    ! Polynomial interpolation (Lagrange)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                          INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                   INTENT(INOUT) :: ocean
    REAL(dp),                                 INTENT(IN)    :: time
    INTEGER,                                  INTENT(IN)    :: num_timeframes
    REAL(dp), DIMENSION(:),                   INTENT(IN)    :: times
    REAL(dp), DIMENSION(:),                   INTENT(OUT)   :: weights

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'polynomial_interpolation'
    INTEGER                                                 :: n, m, i, j, polynomial_order

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Polynomial order based on the number of timeframes
    polynomial_order = num_timeframes - 1

    ! Calculate weights for polynomial interpolation
    DO n = 1, num_timeframes
        weights(n) = 1.0_dp
        DO m = 1, num_timeframes
            IF (m /= n) THEN
                weights(n) = weights(n) * (time - times(m)) / (times(n) - times(m))
            END IF
        END DO
    END DO

    ! Apply polynomial interpolation
    DO i = mesh%vi1, mesh%vi2
      DO j = 1, C%nz_ocean
          ocean%T(i, j) = 0.0_dp
          ocean%S(i, j) = 0.0_dp

          ! Sum the contributions from each available timeframe
          n = 0
          IF (ALLOCATED(ocean%matrix%timeframe0%T)) THEN
              n = n + 1
              IF (n <= num_timeframes) THEN
                  ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe0%T(i, j)
                  ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe0%S(i, j)
              END IF
          END IF
          IF (ALLOCATED(ocean%matrix%timeframe1%T)) THEN
              n = n + 1
              IF (n <= num_timeframes) THEN
                  ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe1%T(i, j)
                  ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe1%S(i, j)
              END IF
          END IF
          IF (ALLOCATED(ocean%matrix%timeframe2%T)) THEN
              n = n + 1
              IF (n <= num_timeframes) THEN
                  ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe2%T(i, j)
                  ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe2%S(i, j)
              END IF
          END IF
          IF (ALLOCATED(ocean%matrix%timeframe3%T)) THEN
            n = n + 1
            IF (n <= num_timeframes) THEN
                ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe3%T(i, j)
                ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe3%S(i, j)
            END IF
          END IF
          ! Limited polynomial order due to computational time
      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE polynomial_interpolation

  SUBROUTINE run_ocean_model_matrix( mesh, ice, ocean, time, region_name)
    ! Calculate the ocean
    !
    ! Use an interpolating matrix ocean scheme, choice between linear and polynomial
  
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_matrix'
    INTEGER                                               :: num_timeframes, required_timeframes
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: times, weights

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Update timeframes if necessary
    IF (time < ocean%matrix%t0 .OR. time > ocean%matrix%t1) THEN
      CALL update_ocean_matrix_timeframes(mesh, ocean, region_name, time)
    END IF

    ! Get the number of available timeframes
    CALL count_allocated_timeframes(ocean, num_timeframes)

    ! Arrays for timeframe and weight recognition
    ALLOCATE(times(num_timeframes))
    ALLOCATE(weights(num_timeframes))

    ! Fill the times array with the available timeframes
    CALL fill_times_array(ocean, times)
    
    ! Minimum required amount of timeframes for each interpolation method
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear') THEN
      required_timeframes = 2
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial') THEN
        required_timeframes = 2  ! Minimum required for polynomial interpolation (capped to 4 for now)
    ELSE
        CALL crash('Unknown interpolation method: ' // TRIM(C%choice_ocean_model_matrix))
        RETURN
    END IF

    ! Check if enough timeframes are available for the selected method
    IF (num_timeframes < required_timeframes) THEN
        CALL crash('Insufficient timeframes for interpolation method.')
        RETURN
    END IF

    ! Check for division by 0 error
    IF (ABS(ocean%matrix%t1 - ocean%matrix%t0) < 1e-8_dp) THEN
      CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
    END IF

    ! Perform time interpolation
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear') THEN
      CALL linear_interpolation(mesh, ocean, time)
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial') THEN
        CALL polynomial_interpolation(mesh, ocean, time, num_timeframes, times, weights)
    ELSE
        CALL crash('Unknown choice_ocean_model_matrix "' // TRIM(C%choice_ocean_model_matrix) // '"')
    END IF

    ! Clean up arrays
    DEALLOCATE(times)
    DEALLOCATE(weights)

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE run_ocean_model_matrix
  
  SUBROUTINE initialise_ocean_model_matrix( mesh, ocean, region_name)
    ! Initialise the ocean matrix model
    !
    ! Use an matrix ocean scheme
  
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_matrix'
    CHARACTER(LEN=256)                                    :: filename_ocean_snapshot
    CHARACTER(LEN=256)                                    :: filename_ocean_matrix_base1, filename_ocean_matrix_base2

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising matrix ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_matrix),'light blue') // '"...'
  
    ! Run the chosen matrix ocean model
    IF (C%choice_ocean_model_matrix == 'snapshot') THEN
      ! Read single-time data from external file
  
      ! Determine which ocean model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
      ELSEIF (region_name == 'EAS') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
      ELSEIF (region_name == 'GRL') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
      ELSEIF (region_name == 'ANT') THEN
        filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF
  
      ! Fill in  main variables
      CALL read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean, mesh, ocean%T)
      CALL read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean, mesh, ocean%S)
  
    ELSE
      CALL crash('unknown choice_ocean_model_matrix "' // TRIM( C%choice_ocean_model_matrix) // '"')
    END IF
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_ocean_model_matrix

END MODULE ocean_matrix