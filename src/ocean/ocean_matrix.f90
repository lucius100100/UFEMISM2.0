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
    USE ocean_model_types                                      , ONLY: type_ocean_model, type_ocean_matrix_interpolation
    USE netcdf_input                                           , ONLY: read_field_from_file_3D_ocean
    USE netcdf_basic                                           , ONLY: field_name_options_T_ocean, field_name_options_S_ocean
  
    IMPLICIT NONE
  
  CONTAINS
  
  ! ===== Main routines =====
  ! =========================

  SUBROUTINE update_ocean_matrix_timeframes(mesh, ocean, matrix, region_name, time)
    ! Update the ocean matrix timeframes

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(INOUT) :: matrix
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'update_ocean_matrix_timeframes'
    CHARACTER(LEN=256)                                    :: filename1, filename2
    INTEGER                                               :: ndepth
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: depth

    ! Add routine to path
    CALL init_routine( routine_name)

    ! For now, hardcode the years
    matrix%t0 = -21000.0_dp   ! LGM
    matrix%t1 = 0.0_dp        ! PI

    ! For now, hardcode the depth 
    !ndepth = 11
    !ALLOCATE(depth(ndepth))
    !depth = (/0.0_dp, 150.0_dp, 300.0_dp, 450.0_dp, 600.0_dp, 750.0_dp, 900.0_dp, &
              !1050.0_dp, 1200.0_dp, 1350.0_dp, 1500.0_dp/)

    ! Allocate memory for timeframes' T and S array if not already allocated
    IF (.NOT. ALLOCATED(matrix%timeframe0%T)) THEN
      ALLOCATE(matrix%timeframe0%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ALLOCATE(matrix%timeframe0%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
    END IF
    IF (.NOT. ALLOCATED(matrix%timeframe1%T)) THEN
      ALLOCATE(matrix%timeframe1%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ALLOCATE(matrix%timeframe1%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))      
    END IF

    ! Construct filenames for the two ocean snapshots
    filename1 = TRIM(C%filename_ocean_matrix_base1)
    filename2 = TRIM(C%filename_ocean_matrix_base2)

    ! Read the ocean snapshots
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_T_ocean, mesh, matrix%timeframe0%T)
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_S_ocean, mesh, matrix%timeframe0%S)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_T_ocean, mesh, matrix%timeframe1%T)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_S_ocean, mesh, matrix%timeframe1%S)
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_ocean_matrix_timeframes

  SUBROUTINE linear_interpolation(mesh, ocean, matrix, time)
    ! Linear interpolation between two ocean snapshots

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'linear_interpolation'
    REAL(dp)                                              :: wt0, wt1
    INTEGER                                               :: i, j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate weights for linear interpolation
    wt0 = (matrix%t1 - time) / (matrix%t1 - matrix%t0)
    wt1 = 1.0_dp - wt0

    ! Apply linear interpolation
    DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
            ocean%T(i,j) = wt0 * matrix%timeframe0%T(i,j) + wt1 * matrix%timeframe1%T(i,j)
            ocean%S(i,j) = wt0 * matrix%timeframe0%S(i,j) + wt1 * matrix%timeframe1%S(i,j)
        END DO
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE linear_interpolation

  !SUBROUTINE polynomial_interpolation(mesh, ocean, time, num_timeframes, times, weights)
    ! Polynomial interpolation (Lagrange)

    !IMPLICIT NONE

    ! In/output variables:
    !TYPE(type_mesh),                          INTENT(IN)    :: mesh
    !TYPE(type_ocean_model),                   INTENT(INOUT) :: ocean
    !REAL(dp),                                 INTENT(IN)    :: time
    !INTEGER,                                  INTENT(IN)    :: num_timeframes
    !REAL(dp), DIMENSION(:),                   INTENT(IN)    :: times
    !REAL(dp), DIMENSION(:),                   INTENT(OUT)   :: weights

    ! Local variables:
    !CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'polynomial_interpolation'
    !INTEGER                                               :: n, m, i, j, polynomial_order

    ! Add routine to path
    !CALL init_routine( routine_name)

    ! Polynomial order based on the number of timeframes
    !polynomial_order = num_timeframes - 1

    ! Calculate weights for polynomial interpolation
    !DO n = 1, num_timeframes
        !weights(n) = 1.0_dp
        !DO m = 1, num_timeframes
            !IF (m /= n) THEN
                !weights(n) = weights(n) * (time - times(m)) / (times(n) - times(m))
            !END IF
        !END DO
    !END DO

    ! Apply polynomial interpolation
    !DO i = mesh%vi1, mesh%vi2
      !DO j = 1, C%nz_ocean
          !ocean%T(i, j) = 0.0_dp
          !ocean%S(i, j) = 0.0_dp

          ! Sum the contributions from each available timeframe
          !n = 0
          !IF (ALLOCATED(matrix%timeframe0%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * matrix%timeframe0%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * matrix%timeframe0%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(matrix%timeframe1%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * matrix%timeframe1%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * matrix%timeframe1%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(matrix%timeframe2%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * matrix%timeframe2%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * matrix%timeframe2%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(matrix%timeframe3%T)) THEN
            !n = n + 1
            !IF (n <= num_timeframes) THEN
                !ocean%T(i, j) = ocean%T(i, j) + weights(n) * matrix%timeframe3%T(i, j)
                !ocean%S(i, j) = ocean%S(i, j) + weights(n) * matrix%timeframe3%S(i, j)
            !END IF
          !END IF
          ! Limited polynomial order due to computational time
      !END DO
    !END DO

    ! Finalise routine path
    !CALL finalise_routine( routine_name)
  
  !END SUBROUTINE polynomial_interpolation

  SUBROUTINE run_ocean_model_matrix( mesh, ice, ocean, time, region_name)
    ! Calculate the ocean
    !
    ! Use an interpolating matrix ocean scheme
  
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_matrix'
    TYPE(type_ocean_matrix_interpolation)                 :: matrix
    INTEGER                                               :: required_timeframes

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Update timeframes if necessary
    IF (time < matrix%t0 .OR. time > matrix%t1) THEN
      CALL update_ocean_matrix_timeframes(mesh, ocean, matrix, region_name, time)
    END IF
    
    PRINT *, 'After timeframe update:'
    PRINT *, 'Time:', time
    PRINT *, 'matrix%t0:', matrix%t0
    PRINT *, 'matrix%t1:', matrix%t1    

    ! Minimum required amount of timeframes for each interpolation method
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear') THEN
      required_timeframes = 2
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial') THEN
        required_timeframes = 3 
    ELSE
        CALL crash('Unknown interpolation method: ' // TRIM(C%choice_ocean_model_matrix))
        RETURN
    END IF

    ! Check for division by 0 error
    !IF (ABS(matrix%t1 - matrix%t0) < 1e-5_dp) THEN
      !CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
    !END IF

    ! Perform time interpolation
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear') THEN
      CALL linear_interpolation(mesh, ocean, matrix, time)
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial') THEN
        CALL crash('Polynomial interpolation not implemented yet')
    ELSE
        CALL crash('Unknown choice_ocean_model_matrix' // TRIM(C%choice_ocean_model_matrix))
    END IF

    PRINT *, 'After interpolation:'
    PRINT *, 'Time:', time
    PRINT *, 'matrix%t0:', matrix%t0
    PRINT *, 'matrix%t1:', matrix%t1

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE run_ocean_model_matrix
  
  SUBROUTINE initialise_ocean_model_matrix( mesh, ocean, region_name)
    ! Initialise the ocean matrix model
  
    IMPLICIT NONE
  
    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_matrix'
    TYPE(type_ocean_matrix_interpolation)                 :: matrix

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising matrix ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_matrix),'light blue') // '"...'

    ! Run the chosen matrix ocean model
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear') THEN
      ! For linear interpolation, no initialisation is required here
      ! Timeframes will be updated when `run_ocean_model_matrix` is called
    ELSE
      CALL crash('Unknown choice_ocean_model_matrix: "' // TRIM(C%choice_ocean_model_matrix) // '"')
    END IF
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_ocean_model_matrix

END MODULE ocean_matrix