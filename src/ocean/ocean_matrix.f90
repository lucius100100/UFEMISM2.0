MODULE ocean_matrix

    ! matrix ocean models
  
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
    real(dp)                                              :: wt0, wt1
    INTEGER                                               :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Update timeframes if necessary
    IF (time < ocean%matrix%t0 .OR. time > ocean%matrix%t1) THEN
      CALL update_ocean_matrix_timeframes(mesh, ocean, region_name, time)
    END IF

    ! Check for division by 0 error
    IF (ABS(ocean%matrix%t1 - ocean%matrix%t0) < 1e-6_dp) THEN
      CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
    END IF

    ! Perform time interpolation
    ! Linear interpolation
    IF     (C%choice_ocean_model_matrix == 'linear') THEN
      wt0 = (ocean%matrix%t1 - time) / (ocean%matrix%t1 - ocean%matrix%t0)
      wt1 = 1._dp - wt0

      ! Apply interpolation
      DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
          ocean%T(i,j) = wt0 * ocean%matrix%timeframe0%T(i,j) + wt1 * ocean%matrix%timeframe1%T(i,j)
          ocean%S(i,j) = wt0 * ocean%matrix%timeframe0%S(i,j) + wt1 * ocean%matrix%timeframe1%S(i,j)
        END DO
      END DO

    ELSE IF (C%choice_ocean_model_matrix == 'polynomial') THEN
      ! FIX for later
    ELSE
      CALL crash('unknown choice_ocean_model_matrix "' // TRIM( C%choice_ocean_model_matrix) // '"')
    END IF

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
    INTEGER                                               :: year0, year1
    CHARACTER(LEN=256)                                    :: filename0, filename1
    
    ! Determine two closest years to the current time
    year0 = FLOOR(time)
    year1 = CEILING(time)

    ! Set time for the two snapshots
    ocean%matrix%t0 = REAL(year0, dp)
    ocean%matrix%t1 = REAL(year1, dp)

    ! Allocate memory for timeframe0 and timeframe1
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe0)) THEN
      ALLOCATE(ocean%matrix%timeframe0)
    END IF
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe1)) THEN
        ALLOCATE(ocean%matrix%timeframe1)
    END IF

    ! Construct filenames for the two ocean snapshots
    filename0 = TRIM(C%filename_ocean_matrix_base) // TRIM(region_name) // '_' // TRIM(ADJUSTL(CHAR(year0))) // '.nc'
    filename1 = TRIM(C%filename_ocean_matrix_base) // TRIM(region_name) // '_' // TRIM(ADJUSTL(CHAR(year1))) // '.nc'

    ! Read the ocean snapshots
    CALL read_field_from_file_3D_ocean(filename0, field_name_options_T_ocean, mesh, ocean%matrix%timeframe0%T)
    CALL read_field_from_file_3D_ocean(filename0, field_name_options_S_ocean, mesh, ocean%matrix%timeframe0%S)
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_T_ocean, mesh, ocean%matrix%timeframe1%T)
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_S_ocean, mesh, ocean%matrix%timeframe1%S)
  
  END SUBROUTINE update_ocean_matrix_timeframes

END MODULE ocean_matrix