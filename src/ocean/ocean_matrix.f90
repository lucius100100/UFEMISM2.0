MODULE ocean_matrix

  ! Concept:
  ! w_GHG = w_CO2 + w_CH4 + w_N2O
  ! w_e = w_GHG + w_insolation
  ! w_f = w_TS_relation + w_circulation ! Probably going to ignore this, too complicated to implement / outside scope of research / too little impact on outcome

  ! List of forcing:
  ! GHG concentrations (CO2, CH4, N2O) https://www.ipcc.ch/site/assets/uploads/2018/03/TAR-06.pdf
  ! Insolation variation
  ! Ocean circulation (?) https://cp.copernicus.org/articles/19/1081/2023/cp-19-1081-2023.pdf  CH4 and Antarctic Circumpolar Current (ACC) (?)
  ! Freshwater flux
  ! T and S relationship (empirical relation?)

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

  SUBROUTINE interpolation_with_insolation(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on insolation values

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_insolation'
    REAL(dp)                                              :: w_ins
    INTEGER                                               :: i, j
    REAL(dp)                                              :: ins_current, ins_PI, ins_LGM
    !REAL(dp), PARAMETER                                   :: ins_PI = 440.0_dp
    !REAL(dp), PARAMETER                                   :: ins_LGM = 70.0_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
    CALL get_insolation(time, ins_current)
    CALL get_insolation(0.0_dp, ins_PI) ! Should actually be 100.0_dp, but insolation solution is in ka
    CALL get_insolation(21000.0_dp, ins_LGM)

    ! Compute w_ins
    w_ins = (ins_current - ins_LGM) / (ins_PI - ins_LGM)

    ! Clamp between cutoff values if enabled
    IF (C%clamp_weights) THEN
        w_ins = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_ins))
    END IF

    ! Apply interpolation using w_ins
    DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
            ocean%T(i,j) = w_ins * ocean%matrix%timeframe1%T(i,j) + (1.0_dp - w_ins) * ocean%matrix%timeframe0%T(i,j)
            ocean%S(i,j) = w_ins * ocean%matrix%timeframe1%S(i,j) + (1.0_dp - w_ins) * ocean%matrix%timeframe0%S(i,j)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_insolation

  SUBROUTINE get_insolation(time, ins_current)
    ! Get insolation values at the given time

    IMPLICIT NONE

    ! In- and output variables
    REAL(dp), INTENT(IN)                                  :: time
    REAL(dp), INTENT(OUT)                                 :: ins_current

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'get_insolation'
    INTEGER, PARAMETER                                    :: num_points = 22
    REAL(dp), DIMENSION(num_points)                       :: age_data
    REAL(dp), DIMENSION(num_points)                       :: ins_data
    INTEGER                                               :: i
    LOGICAL                                               :: found

    ! Add routine to path
    CALL init_routine( routine_name)

    ! FIX
    ! Hardcoded values for now, should be reading in from csv
    
    ! Berger, A; Loutre, Marie-France (1999)
    ! Initialize data arrays
    age_data = (/ &
    0.0_dp, 1000.0_dp, 2000.0_dp, 3000.0_dp, 4000.0_dp, 5000.0_dp, 6000.0_dp, 7000.0_dp, 8000.0_dp, 9000.0_dp, 10000.0_dp, 11000.0_dp, &
    12000.0_dp, 13000.0_dp, 14000.0_dp, 15000.0_dp, 16000.0_dp, 17000.0_dp, 18000.0_dp, 19000.0_dp, 20000.0_dp, 21000.0_dp /)
    
    ins_data = (/ &
	  426.76_dp, 430.12_dp, 434.69_dp, 440.20_dp, 446.28_dp, 452.48_dp, 458.31_dp, 463.29_dp, 467.00_dp, 469.12_dp, &
    469.44_dp, 467.92_dp, 464.67_dp, 459.95_dp, 454.12_dp, 447.62_dp, 440.92_dp, 434.50_dp, 428.77_dp, 424.07_dp, &
    420.64_dp, 418.62_dp /)

    ! Initialize found flag
    found = .FALSE.
    
    ! Search for the time in age_data
    DO i = 1, num_points
        IF (ABS(age_data(i) - time) < 1e-3_dp) THEN
            ins_current = ins_data(i)
            found = .TRUE.
            EXIT
        END IF
    END DO
    
    IF (.NOT. found) THEN
        CALL crash('Time value not found in age_data.')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation

  SUBROUTINE linear_time_interpolation(mesh, ocean, matrix, time)
    ! Linear interpolation between two ocean snapshots

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'linear_time_interpolation'
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
    
    ! Check for division by 0 error
    !IF (ABS(ocean%matrix%t1 - ocean%matrix%t0) < 1e-5_dp) THEN
      !CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
    !END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE linear_time_interpolation

  SUBROUTINE interpolation_with_GHG_basic(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on GHG

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_GHG_basic'
    REAL(dp)                                              :: w_GHG, w_GHG_CO2, w_GHG_CH4, w_GHG_N2O
    INTEGER                                               :: i, j
    REAL(dp)                                              :: CO2_current, CH4_current, N2O_current
    REAL(dp)                                              :: CO2_PI, CH4_PI, N2O_PI
    REAL(dp)                                              :: CO2_LGM, CH4_LGM, N2O_LGM
    CHARACTER(LEN=256)                                    :: GHG_inclusion

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
    CALL get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)
    CALL get_GHG_concentrations(100.0_dp, CO2_PI, CH4_PI, N2O_PI)
    CALL get_GHG_concentrations(21000.0_dp, CO2_LGM, CH4_LGM, N2O_LGM)

    ! Get config settings
    GHG_inclusion    = TRIM(C%choice_ghg_inclusion)

    ! Compute weights based on GHG concentration ratios
    ! w_CO2
    IF (GHG_inclusion == 'CO2' .OR. GHG_inclusion == 'CO2_CH4' .OR. GHG_inclusion == 'CO2_CH4_N2O') THEN
      w_GHG_CO2 = (CO2_current - CO2_LGM) / (CO2_PI - CO2_LGM)
    ELSE
        w_GHG_CO2 = 0.0_dp
    END IF

    ! w_CH4
    IF (GHG_inclusion == 'CO2_CH4' .OR. GHG_inclusion == 'CO2_CH4_N2O') THEN
        w_GHG_CH4 = (CH4_current - CH4_LGM) / (CH4_PI - CH4_LGM)
    ELSE
        w_GHG_CH4 = 0.0_dp
    END IF

    ! w_N2O
    IF (GHG_inclusion == 'CO2_CH4_N2O') THEN
        w_GHG_N2O = (N2O_current - N2O_LGM) / (N2O_PI - N2O_LGM)
    ELSE
        w_GHG_N2O = 0.0_dp
    END IF

    ! Combine weights based on selected GHGs
    IF (GHG_inclusion == 'CO2') THEN
        w_GHG = w_GHG_CO2
    ELSE IF (GHG_inclusion == 'CO2_CH4') THEN
        w_GHG = (w_GHG_CO2 + w_GHG_CH4) / 2.0_dp
    ELSE IF (GHG_inclusion == 'CO2_CH4_N2O') THEN
        w_GHG = (w_GHG_CO2 + w_GHG_CH4 + w_GHG_N2O) / 3.0_dp
    ELSE
        CALL crash('Unknown choice_ghg_inclusion: "' // TRIM(GHG_inclusion) // '"')
    END IF

    ! Clamp between cutoff values if enabled
    IF (C%clamp_weights) THEN
        w_GHG = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_GHG))
    END IF

    ! Apply interpolation using w_GHG
    DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
            ocean%T(i,j) = w_GHG * ocean%matrix%timeframe1%T(i,j) + (1.0_dp - w_GHG) * ocean%matrix%timeframe0%T(i,j)
            ocean%S(i,j) = w_GHG * ocean%matrix%timeframe1%S(i,j) + (1.0_dp - w_GHG) * ocean%matrix%timeframe0%S(i,j)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_GHG_basic

  SUBROUTINE interpolation_with_GHG_radiative(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on GHG

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_GHG_radiative'
    REAL(dp)                                              :: w_GHG
    INTEGER                                               :: i, j
    REAL(dp)                                              :: CO2_current, CH4_current, N2O_current
    REAL(dp)                                              :: CO2_PI, CH4_PI, N2O_PI
    REAL(dp)                                              :: CO2_LGM, CH4_LGM, N2O_LGM
    REAL(dp)                                              :: DeltaF_CO2, DeltaF_CH4, DeltaF_N2O
    REAL(dp)                                              :: DeltaF_CO2_PI, DeltaF_CH4_PI, DeltaF_N2O_PI
    REAL(dp)                                              :: DeltaF_CO2_LGM, DeltaF_CH4_LGM, DeltaF_N2O_LGM
    REAL(dp)                                              :: DeltaF_total, DeltaF_PI, DeltaF_LGM
    CHARACTER(LEN=256)                                    :: CO2_relationship
    CHARACTER(LEN=256)                                    :: GHG_inclusion

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
    CALL get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)
    CALL get_GHG_concentrations(100.0_dp, CO2_PI, CH4_PI, N2O_PI)
    CALL get_GHG_concentrations(21000.0_dp, CO2_LGM, CH4_LGM, N2O_LGM)

    ! Get config settings
    GHG_inclusion    = TRIM(C%choice_ghg_inclusion)
    CO2_relationship = TRIM(C%choice_CO2_relationship)

    ! Initialize DeltaFs
    DeltaF_total = 0.0_dp
    DeltaF_PI    = 0.0_dp
    DeltaF_LGM   = 0.0_dp

    ! Compute CO2 forcing
    IF (ghg_inclusion == 'CO2' .OR. ghg_inclusion == 'CO2_CH4' .OR. ghg_inclusion == 'CO2_CH4_N2O') THEN
      CALL compute_CO2_forcing(CO2_current, CO2_LGM, CO2_relationship, DeltaF_CO2)
      CALL compute_CO2_forcing(CO2_PI,     CO2_LGM, CO2_relationship, DeltaF_CO2_PI)
      CALL compute_CO2_forcing(CO2_LGM,    CO2_LGM, CO2_relationship, DeltaF_CO2_LGM) 

      DeltaF_total = DeltaF_total + DeltaF_CO2
      DeltaF_PI    = DeltaF_PI    + DeltaF_CO2_PI
      DeltaF_LGM   = DeltaF_LGM   + DeltaF_CO2_LGM
    END IF

    ! Include CH4 forcing if selected
    IF (ghg_inclusion == 'CO2_CH4' .OR. ghg_inclusion == 'CO2_CH4_N2O') THEN
      CALL compute_CH4_forcing(CH4_current, CH4_LGM, N2O_LGM, DeltaF_CH4)
      CALL compute_CH4_forcing(CH4_PI,      CH4_LGM, N2O_LGM, DeltaF_CH4_PI)
      CALL compute_CH4_forcing(CH4_LGM,     CH4_LGM, N2O_LGM, DeltaF_CH4_LGM)

      DeltaF_total = DeltaF_total + DeltaF_CH4
      DeltaF_PI    = DeltaF_PI    + DeltaF_CH4_PI
      DeltaF_LGM   = DeltaF_LGM   + DeltaF_CH4_LGM
    END IF

    ! Include N2O forcing if selected
    IF (ghg_inclusion == 'CO2_CH4_N2O') THEN
      CALL compute_N2O_forcing(N2O_current, N2O_LGM, CH4_LGM, DeltaF_N2O)
      CALL compute_N2O_forcing(N2O_PI,      N2O_LGM, CH4_LGM, DeltaF_N2O_PI)
      CALL compute_N2O_forcing(N2O_LGM,     N2O_LGM, CH4_LGM, DeltaF_N2O_LGM)

      DeltaF_total = DeltaF_total + DeltaF_N2O
      DeltaF_PI    = DeltaF_PI    + DeltaF_N2O_PI
      DeltaF_LGM   = DeltaF_LGM   + DeltaF_N2O_LGM
    END IF

    ! Compute w_GHG
    w_GHG = (DeltaF_total - DeltaF_LGM) / (DeltaF_PI - DeltaF_LGM)

    ! Clamp between cutoff values if enabled
    IF (C%clamp_weights) THEN
        w_GHG = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_GHG))
    END IF

    ! Apply interpolation using w_GHG
    DO i = mesh%vi1, mesh%vi2
        DO j = 1, C%nz_ocean
            ocean%T(i,j) = w_GHG * ocean%matrix%timeframe1%T(i,j) + (1.0_dp - w_GHG) * ocean%matrix%timeframe0%T(i,j)
            ocean%S(i,j) = w_GHG * ocean%matrix%timeframe1%S(i,j) + (1.0_dp - w_GHG) * ocean%matrix%timeframe0%S(i,j)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_GHG_radiative

  SUBROUTINE get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)
    ! Get GHG concentrations at the given time

    IMPLICIT NONE

    ! In- and output variables
    REAL(dp), INTENT(IN)                                  :: time
    REAL(dp), INTENT(OUT)                                 :: CO2_current
    REAL(dp), INTENT(OUT)                                 :: CH4_current
    REAL(dp), INTENT(OUT)                                 :: N2O_current

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'get_GHG_concentrations'
    INTEGER, PARAMETER                                    :: num_points = 373
    REAL(dp), DIMENSION(num_points)                       :: age_data
    REAL(dp), DIMENSION(num_points)                       :: CO2_data
    REAL(dp), DIMENSION(num_points)                       :: CH4_data
    REAL(dp), DIMENSION(num_points)                       :: N2O_data
    INTEGER                                               :: i
    !LOGICAL                                               :: found
    INTEGER                                               :: idx_low, idx_high
    REAL(dp)                                              :: t_low, t_high, fraction

    ! Add routine to path
    CALL init_routine( routine_name)

    ! FIX
    ! Hardcoded values for now, should be reading in from csv
    ! https://stackoverflow.com/questions/8828377/reading-data-from-txt-file-in-fortran
    
    ! Initialize data arrays
    age_data = (/ &
    0.0_dp, 100.0_dp, 137.0_dp, 148.0_dp, 212.0_dp, 268.0_dp, 279.0_dp, 280.0_dp, 395.0_dp, 400.0_dp, &
    404.0_dp, 485.0_dp, 513.0_dp, 559.0_dp, 572.0_dp, 672.0_dp, 677.0_dp, 754.0_dp, 768.0_dp, 769.0_dp, &
    775.0_dp, 874.0_dp, 877.0_dp, 950.0_dp, 951.0_dp, 1029.0_dp, 1060.0_dp, 1116.0_dp, 1153.0_dp, 1201.0_dp, &
    1233.0_dp, 1305.0_dp, 1350.0_dp, 1402.0_dp, 1453.0_dp, 1489.0_dp, 1552.0_dp, 1571.0_dp, 1638.0_dp, 1673.0_dp, &
    1733.0_dp, 1734.0_dp, 1812.0_dp, 1855.0_dp, 1931.0_dp, 1982.0_dp, 2057.0_dp, 2069.0_dp, 2128.0_dp, 2134.0_dp, &
    2212.0_dp, 2262.0_dp, 2334.0_dp, 2366.0_dp, 2433.0_dp, 2453.0_dp, 2456.0_dp, 2536.0_dp, 2552.0_dp, 2604.0_dp, &
    2657.0_dp, 2728.0_dp, 2741.0_dp, 2804.0_dp, 2806.0_dp, 2902.0_dp, 2977.0_dp, 3040.0_dp, 3053.0_dp, 3116.0_dp, &
    3133.0_dp, 3215.0_dp, 3230.0_dp, 3336.0_dp, 3363.0_dp, 3450.0_dp, 3453.0_dp, 3523.0_dp, 3531.0_dp, 3621.0_dp, &
    3622.0_dp, 3626.0_dp, 3714.0_dp, 3721.0_dp, 3790.0_dp, 3842.0_dp, 3910.0_dp, 3946.0_dp, 4004.0_dp, 4036.0_dp, &
    4040.0_dp, 4096.0_dp, 4097.0_dp, 4161.0_dp, 4249.0_dp, 4257.0_dp, 4324.0_dp, 4331.0_dp, 4374.0_dp, 4447.0_dp, &
    4480.0_dp, 4526.0_dp, 4573.0_dp, 4666.0_dp, 4667.0_dp, 4703.0_dp, 4737.0_dp, 4766.0_dp, 4843.0_dp, 4874.0_dp, &
    4984.0_dp, 5004.0_dp, 5051.0_dp, 5094.0_dp, 5160.0_dp, 5170.0_dp, 5274.0_dp, 5282.0_dp, 5370.0_dp, 5387.0_dp, &
    5469.0_dp, 5470.0_dp, 5476.0_dp, 5562.0_dp, 5580.0_dp, 5657.0_dp, 5682.0_dp, 5711.0_dp, 5716.0_dp, 5855.0_dp, &
    5868.0_dp, 5998.0_dp, 6026.0_dp, 6039.0_dp, 6058.0_dp, 6131.0_dp, 6163.0_dp, 6263.0_dp, 6318.0_dp, 6354.0_dp, &
    6397.0_dp, 6434.0_dp, 6470.0_dp, 6545.0_dp, 6617.0_dp, 6622.0_dp, 6694.0_dp, 6713.0_dp, 6783.0_dp, 6838.0_dp, &
    6922.0_dp, 6941.0_dp, 7028.0_dp, 7045.0_dp, 7112.0_dp, 7158.0_dp, 7234.0_dp, 7258.0_dp, 7320.0_dp, 7375.0_dp, &
    7413.0_dp, 7489.0_dp, 7507.0_dp, 7590.0_dp, 7600.0_dp, 7691.0_dp, 7716.0_dp, 7731.0_dp, 7751.0_dp, 7781.0_dp, &
    7804.0_dp, 7876.0_dp, 7928.0_dp, 7961.0_dp, 7986.0_dp, 7990.0_dp, 8041.0_dp, 8050.0_dp, 8096.0_dp, 8121.0_dp, &
    8148.0_dp, 8173.0_dp, 8181.0_dp, 8191.0_dp, 8211.0_dp, 8250.0_dp, 8281.0_dp, 8298.0_dp, 8327.0_dp, 8374.0_dp, &
    8387.0_dp, 8393.0_dp, 8412.0_dp, 8475.0_dp, 8476.0_dp, 8477.0_dp, 8488.0_dp, 8579.0_dp, 8610.0_dp, 8653.0_dp, &
    8702.0_dp, 8784.0_dp, 8802.0_dp, 8869.0_dp, 8904.0_dp, 8973.0_dp, 9016.0_dp, 9092.0_dp, 9140.0_dp, 9232.0_dp, &
    9241.0_dp, 9317.0_dp, 9347.0_dp, 9480.0_dp, 9536.0_dp, 9545.0_dp, 9597.0_dp, 9615.0_dp, 9721.0_dp, 9739.0_dp, &
    9807.0_dp, 9909.0_dp, 9946.0_dp, 9983.0_dp, 10024.0_dp, 10088.0_dp, 10133.0_dp, 10209.0_dp, 10240.0_dp, 10294.0_dp, &
    10298.0_dp, 10388.0_dp, 10417.0_dp, 10514.0_dp, 10527.0_dp, 10613.0_dp, 10621.0_dp, 10682.0_dp, 10744.0_dp, 10805.0_dp, &
    10824.0_dp, 10827.0_dp, 10895.0_dp, 10933.0_dp, 10935.0_dp, 11014.0_dp, 11038.0_dp, 11087.0_dp, 11119.0_dp, 11136.0_dp, &
    11201.0_dp, 11236.0_dp, 11241.0_dp, 11278.0_dp, 11338.0_dp, 11342.0_dp, 11392.0_dp, 11436.0_dp, 11469.0_dp, 11484.0_dp, &
    11558.0_dp, 11580.0_dp, 11631.0_dp, 11635.0_dp, 11676.0_dp, 11727.0_dp, 11765.0_dp, 11819.0_dp, 11865.0_dp, 11896.0_dp, &
    11958.0_dp, 12050.0_dp, 12085.0_dp, 12122.0_dp, 12167.0_dp, 12188.0_dp, 12294.0_dp, 12371.0_dp, 12411.0_dp, 12496.0_dp, &
    12545.0_dp, 12642.0_dp, 12760.0_dp, 12849.0_dp, 12942.0_dp, 12994.0_dp, 13090.0_dp, 13116.0_dp, 13241.0_dp, 13253.0_dp, &
    13400.0_dp, 13440.0_dp, 13542.0_dp, 13571.0_dp, 13653.0_dp, 13727.0_dp, 13804.0_dp, 13890.0_dp, 13937.0_dp, 13948.0_dp, &
    14089.0_dp, 14196.0_dp, 14303.0_dp, 14353.0_dp, 14550.0_dp, 14671.0_dp, 14725.0_dp, 14868.0_dp, 14890.0_dp, 15012.0_dp, &
    15035.0_dp, 15175.0_dp, 15233.0_dp, 15274.0_dp, 15438.0_dp, 15491.0_dp, 15570.0_dp, 15635.0_dp, 15742.0_dp, 15786.0_dp, &
    15886.0_dp, 15957.0_dp, 16071.0_dp, 16073.0_dp, 16230.0_dp, 16260.0_dp, 16391.0_dp, 16452.0_dp, 16548.0_dp, 16659.0_dp, &
    16738.0_dp, 16870.0_dp, 16898.0_dp, 16910.0_dp, 17093.0_dp, 17111.0_dp, 17323.0_dp, 17375.0_dp, 17494.0_dp, 17565.0_dp, &
    17730.0_dp, 17809.0_dp, 17855.0_dp, 17943.0_dp, 18163.0_dp, 18285.0_dp, 18425.0_dp, 18541.0_dp, 18797.0_dp, 18798.0_dp, &
    18828.0_dp, 18868.0_dp, 18921.0_dp, 19221.0_dp, 19347.0_dp, 19366.0_dp, 19509.0_dp, 19585.0_dp, 19591.0_dp, 19597.0_dp, &
    19748.0_dp, 19871.0_dp, 19988.0_dp, 20009.0_dp, 20168.0_dp, 20197.0_dp, 20355.0_dp, 20357.0_dp, 20502.0_dp, 20578.0_dp, &
    20748.0_dp, 20848.0_dp, 21000.0_dp /)

    CO2_data = (/ &
    280.4_dp, 280.4_dp, 280.4_dp, 279.9381679_dp, 277.251145_dp, 274.9_dp, 277.9_dp, 277.9103448_dp, 279.1_dp, 280.6555556_dp, &
    281.9_dp, 277.7_dp, 278.9864865_dp, 281.1_dp, 281.2265487_dp, 282.2_dp, 282.0719512_dp, 280.1_dp, 279.9065041_dp, 279.8926829_dp, &
    279.8097561_dp, 278.4414634_dp, 278.4_dp, 276.6_dp, 276.6227273_dp, 278.3954545_dp, 279.1_dp, 278.2569892_dp, 277.7_dp, 278.3_dp, &
    278.7_dp, 277.9_dp, 277.4_dp, 278.3087379_dp, 279.2_dp, 279.4909091_dp, 280.0_dp, 279.7569767_dp, 278.9_dp, 278.8263158_dp, &
    278.7_dp, 278.6911392_dp, 278.0_dp, 277.602521_dp, 276.9_dp, 276.8190476_dp, 276.7_dp, 276.7_dp, 276.7_dp, 276.7642857_dp, &
    277.6_dp, 277.7229508_dp, 277.9_dp, 276.6070707_dp, 273.9_dp, 274.8708738_dp, 275.0165049_dp, 278.9_dp, 278.0529412_dp, 275.3_dp, &
    275.0435484_dp, 274.7_dp, 274.9666667_dp, 276.2589744_dp, 276.3_dp, 274.6_dp, 275.4443709_dp, 276.1536424_dp, 276.3_dp, 273.1_dp, &
    273.2545455_dp, 274.0_dp, 274.1239669_dp, 275.0_dp, 274.6307692_dp, 273.4410256_dp, 273.4_dp, 273.0_dp, 272.8787879_dp, 271.5151515_dp, &
    271.5_dp, 271.6575758_dp, 275.1242424_dp, 275.4_dp, 274.9_dp, 273.5133333_dp, 271.7_dp, 271.6617021_dp, 271.6_dp, 272.0173913_dp, &
    272.0695652_dp, 272.8_dp, 272.78_dp, 271.5_dp, 271.2840491_dp, 271.2644172_dp, 271.1_dp, 270.82_dp, 269.1_dp, 269.5820755_dp, &
    269.8_dp, 270.6408602_dp, 271.5_dp, 270.9276923_dp, 270.9215385_dp, 270.7_dp, 269.9444444_dp, 269.3_dp, 268.8009259_dp, 268.6_dp, &
    269.6153846_dp, 269.8_dp, 268.6511111_dp, 267.6_dp, 265.3_dp, 265.2912281_dp, 265.2_dp, 265.4_dp, 267.6_dp, 267.3273585_dp, &
    266.0122642_dp, 265.9962264_dp, 265.9_dp, 265.5_dp, 264.5905263_dp, 260.7_dp, 263.2423729_dp, 266.1915254_dp, 266.7_dp, 265.5_dp, &
    265.2909091_dp, 263.2_dp, 262.8585366_dp, 262.7_dp, 262.3902174_dp, 261.2_dp, 261.1757576_dp, 261.1_dp, 260.0725275_dp, 259.4_dp, &
    260.85125_dp, 262.1_dp, 262.3594595_dp, 262.9_dp, 258.1_dp, 258.0739583_dp, 257.6989583_dp, 257.6_dp, 260.232_dp, 262.3_dp, &
    262.8708738_dp, 263.0_dp, 260.7_dp, 260.2345238_dp, 258.4_dp, 259.0409836_dp, 260.1_dp, 260.1837209_dp, 260.4_dp, 259.9860215_dp, &
    259.7_dp, 259.2957447_dp, 259.2_dp, 260.8_dp, 260.6811881_dp, 259.6_dp, 259.5166667_dp, 259.4666667_dp, 259.4_dp, 259.3_dp, &
    259.0578947_dp, 258.3_dp, 259.6684211_dp, 260.5368421_dp, 261.1947368_dp, 261.3_dp, 260.79_dp, 260.7_dp, 261.0862595_dp, 261.2961832_dp, &
    261.5229008_dp, 261.7328244_dp, 261.8_dp, 261.52_dp, 260.96_dp, 259.868_dp, 259.0_dp, 259.304717_dp, 259.8245283_dp, 260.6669811_dp, &
    260.9_dp, 260.8666667_dp, 260.7611111_dp, 260.4111111_dp, 260.4055556_dp, 260.4_dp, 260.2813725_dp, 259.3_dp, 260.4310811_dp, 262.0_dp, &
    262.6358779_dp, 263.7_dp, 263.7211765_dp, 263.8_dp, 264.2711538_dp, 265.2_dp, 263.5378151_dp, 260.6_dp, 260.9_dp, 263.0_dp, &
    263.0847059_dp, 263.8_dp, 263.8821918_dp, 264.2465753_dp, 264.4_dp, 264.3704918_dp, 264.2_dp, 264.1709677_dp, 264.0_dp, 263.8744186_dp, &
    263.4_dp, 265.7_dp, 265.3_dp, 264.9_dp, 265.9152381_dp, 267.5_dp, 267.2768595_dp, 266.9_dp, 266.5717647_dp, 266.0_dp, &
    265.9707317_dp, 265.3121951_dp, 265.1_dp, 267.3045455_dp, 267.6_dp, 265.0382979_dp, 264.8_dp, 264.8_dp, 264.8_dp, 265.0_dp, &
    265.2590909_dp, 265.3_dp, 264.4_dp, 264.1_dp, 264.1024691_dp, 264.2_dp, 264.2986301_dp, 264.5_dp, 264.1734694_dp, 264.0_dp, &
    263.0_dp, 265.2_dp, 264.4380952_dp, 258.8_dp, 260.8_dp, 260.4_dp, 255.4_dp, 253.9_dp, 253.8_dp, 253.3810811_dp, &
    251.3144144_dp, 250.7_dp, 249.7727273_dp, 249.7_dp, 251.1_dp, 250.7_dp, 248.4695652_dp, 245.3_dp, 245.3_dp, 245.3_dp, &
    246.6_dp, 243.2_dp, 241.7902778_dp, 240.3_dp, 239.7939759_dp, 239.5578313_dp, 238.3658635_dp, 237.5_dp, 237.532_dp, 237.6_dp, &
    236.4589041_dp, 234.2_dp, 238.3_dp, 237.810989_dp, 237.3_dp, 237.5108108_dp, 237.9_dp, 237.8483444_dp, 237.6_dp, 237.5276382_dp, &
    236.641206_dp, 236.4_dp, 239.2_dp, 239.0432432_dp, 238.6_dp, 238.6_dp, 238.6_dp, 238.8986111_dp, 239.0618056_dp, 239.1_dp, &
    234.8898592_dp, 231.6949296_dp, 228.5_dp, 228.4797571_dp, 228.4_dp, 226.8097143_dp, 226.1_dp, 225.32_dp, 225.2_dp, 224.5_dp, &
    224.239819_dp, 222.6561086_dp, 222.0_dp, 221.8_dp, 221.0_dp, 220.9598485_dp, 220.9_dp, 220.3331395_dp, 219.4_dp, 217.75_dp, &
    214.0_dp, 211.5320856_dp, 207.5695187_dp, 207.5_dp, 207.6679144_dp, 207.7_dp, 204.425_dp, 202.9_dp, 201.926087_dp, 200.8_dp, &
    198.7033175_dp, 195.2_dp, 195.0489627_dp, 194.9842324_dp, 193.9970954_dp, 193.9_dp, 191.5712121_dp, 191.0_dp, 189.4342105_dp, 188.5_dp, &
    188.5_dp, 188.5_dp, 188.7402985_dp, 189.2_dp, 187.7847953_dp, 187.0_dp, 187.875_dp, 188.6_dp, 189.3135889_dp, 189.3163763_dp, &
    189.4_dp, 192.3_dp, 188.3_dp, 188.5816901_dp, 188.7_dp, 188.7117284_dp, 188.8_dp, 189.18159_dp, 189.2117155_dp, 189.241841_dp, &
    190.0_dp, 188.975_dp, 188.0_dp, 188.0233333_dp, 188.2_dp, 195.0_dp, 191.2701639_dp, 191.2229508_dp, 187.8_dp, 187.5219512_dp, &
    186.9_dp, 186.7479087_dp, 186.51673_dp /)

    CH4_data = (/ &
    907.0_dp, 808.8243243_dp, 772.4993243_dp, 761.7_dp, 682.7_dp, 676.3588235_dp, 675.1132353_dp, 675.0_dp, 683.9125_dp, 684.3_dp, &
    682.8840708_dp, 654.2115044_dp, 644.3_dp, 676.0322034_dp, 685.0_dp, 667.5714286_dp, 666.7_dp, 677.5307692_dp, 679.5_dp, 674.3285714_dp, &
    643.3_dp, 647.0_dp, 646.0922078_dp, 624.0025974_dp, 623.7_dp, 634.7_dp, 637.9068966_dp, 643.7_dp, 637.1705882_dp, 628.7_dp, &
    626.0230769_dp, 620.0_dp, 625.242268_dp, 631.3_dp, 631.5344828_dp, 631.7_dp, 625.2463415_dp, 623.3_dp, 613.25_dp, 608.0_dp, &
    612.9180328_dp, 613.0_dp, 613.1933884_dp, 613.3_dp, 607.5551181_dp, 603.7_dp, 608.0103448_dp, 608.7_dp, 602.3461538_dp, 601.7_dp, &
    598.8359375_dp, 597.0_dp, 593.3307692_dp, 591.7_dp, 595.2733333_dp, 596.34_dp, 596.5_dp, 592.5_dp, 591.7_dp, 585.9057143_dp, &
    580.0_dp, 578.5630952_dp, 578.3_dp, 577.0_dp, 576.950289_dp, 574.5641618_dp, 572.7_dp, 576.0_dp, 574.8397849_dp, 569.2172043_dp, &
    567.7_dp, 579.1969072_dp, 581.3_dp, 574.4458647_dp, 572.7_dp, 575.2892857_dp, 575.3785714_dp, 577.4619048_dp, 577.7_dp, 578.0_dp, &
    578.54_dp, 580.7_dp, 569.0_dp, 569.0164063_dp, 569.178125_dp, 569.3_dp, 568.2538462_dp, 567.7_dp, 572.8555556_dp, 575.7_dp, &
    563.0_dp, 570.1719298_dp, 570.3_dp, 573.2473684_dp, 577.3_dp, 567.7_dp, 569.5108108_dp, 569.7_dp, 566.7344828_dp, 561.7_dp, &
    561.7_dp, 561.7_dp, 562.1333333_dp, 562.9907801_dp, 563.0_dp, 563.1542857_dp, 563.3_dp, 564.5037736_dp, 567.7_dp, 566.2269504_dp, &
    561.0_dp, 561.5970149_dp, 563.0_dp, 561.3016807_dp, 558.694958_dp, 558.3_dp, 554.0285714_dp, 553.7_dp, 563.1704762_dp, 565.0_dp, &
    565.9879518_dp, 566.0_dp, 566.1254545_dp, 567.9236364_dp, 568.3_dp, 569.3568627_dp, 569.7_dp, 584.3_dp, 584.0038217_dp, 575.7700637_dp, &
    575.0_dp, 565.9493671_dp, 564.0_dp, 570.90625_dp, 581.0_dp, 577.5238095_dp, 576.0_dp, 576.6451613_dp, 577.0_dp, 578.3670886_dp, &
    580.0_dp, 583.9027397_dp, 587.7_dp, 597.2230263_dp, 606.3651316_dp, 607.0_dp, 595.3_dp, 597.0078652_dp, 603.3_dp, 604.0913669_dp, &
    605.3_dp, 605.5626016_dp, 606.7650407_dp, 607.0_dp, 608.600885_dp, 609.7_dp, 607.876_dp, 607.3_dp, 613.1290598_dp, 618.3_dp, &
    612.7666667_dp, 601.7_dp, 603.3216216_dp, 610.7990991_dp, 611.7_dp, 610.1310345_dp, 609.7_dp, 617.0_dp, 621.0_dp, 621.9622642_dp, &
    622.7_dp, 622.8741935_dp, 623.0_dp, 623.3_dp, 617.7_dp, 617.1909091_dp, 610.7_dp, 611.1909091_dp, 613.7_dp, 620.0_dp, &
    597.7_dp, 612.3_dp, 601.9444444_dp, 589.0_dp, 614.0_dp, 624.3_dp, 629.4666667_dp, 632.3_dp, 605.7_dp, 632.7_dp, &
    645.2210526_dp, 651.0_dp, 640.3_dp, 641.875_dp, 641.9_dp, 641.8_dp, 640.7_dp, 647.1147541_dp, 649.3_dp, 654.1608696_dp, &
    659.7_dp, 654.78_dp, 653.7_dp, 656.0647059_dp, 657.3_dp, 658.3473214_dp, 659.0_dp, 658.0806452_dp, 657.5_dp, 657.3178218_dp, &
    657.3_dp, 666.4056604_dp, 670.0_dp, 664.7_dp, 664.9584615_dp, 665.0_dp, 663.7371429_dp, 663.3_dp, 672.1903226_dp, 673.7_dp, &
    676.5251208_dp, 680.7628019_dp, 682.3_dp, 682.774359_dp, 683.3_dp, 679.1899083_dp, 676.3_dp, 680.5616822_dp, 682.3_dp, 684.8137931_dp, &
    685.0_dp, 678.0_dp, 682.6031746_dp, 698.0_dp, 695.3737374_dp, 678.0_dp, 679.1594203_dp, 688.0_dp, 676.3422535_dp, 664.8725352_dp, &
    661.3_dp, 661.0486486_dp, 655.3513514_dp, 652.1675676_dp, 652.0_dp, 662.7378641_dp, 666.0_dp, 674.8925926_dp, 680.7_dp, 679.6688525_dp, &
    675.7262295_dp, 673.6032787_dp, 673.3_dp, 684.180198_dp, 701.8237624_dp, 703.0_dp, 678.4577465_dp, 656.8605634_dp, 640.6626761_dp, 633.3_dp, &
    654.3_dp, 657.4342466_dp, 664.7_dp, 661.5447761_dp, 629.2037313_dp, 588.9746269_dp, 559.0_dp, 510.4_dp, 469.0_dp, 468.8168182_dp, &
    468.4504545_dp, 467.9068182_dp, 467.7_dp, 471.1485437_dp, 475.3427184_dp, 477.3_dp, 482.7_dp, 469.274359_dp, 462.3_dp, 465.0910448_dp, &
    466.7_dp, 490.7266447_dp, 519.9549342_dp, 542.0_dp, 596.5172414_dp, 627.0_dp, 639.8262295_dp, 643.3_dp, 645.7635036_dp, 646.0_dp, &
    661.7_dp, 659.7350877_dp, 654.7245614_dp, 653.3_dp, 649.3051282_dp, 645.7_dp, 632.9453988_dp, 618.7_dp, 613.0_dp, 614.2302632_dp, &
    630.0_dp, 611.7_dp, 588.9369427_dp, 578.3_dp, 507.9251572_dp, 464.7_dp, 463.6857868_dp, 461.0_dp, 460.7760479_dp, 459.5341317_dp, &
    459.3_dp, 473.3_dp, 474.7060606_dp, 475.7_dp, 457.0327189_dp, 451.0_dp, 451.3840278_dp, 451.7_dp, 454.7470199_dp, 456.0_dp, &
    462.0701754_dp, 466.38_dp, 473.3_dp, 472.6132075_dp, 418.7_dp, 420.3024845_dp, 427.3_dp, 418.5191083_dp, 404.7_dp, 405.8417143_dp, &
    406.6542857_dp, 408.012_dp, 408.3_dp, 381.3_dp, 376.8690073_dp, 376.4331719_dp, 371.3_dp, 375.0707602_dp, 383.7_dp, 378.4652542_dp, &
    366.3_dp, 369.46_dp, 371.3_dp, 372.4428571_dp, 375.3_dp, 373.2977099_dp, 371.0_dp, 373.3946381_dp, 378.6793566_dp, 378.7_dp, &
    377.7780142_dp, 376.5486998_dp, 374.9198582_dp, 365.7_dp, 375.5193103_dp, 377.0_dp, 362.2551111_dp, 354.4186667_dp, 353.8_dp, 354.2692857_dp, &
    366.0796429_dp, 375.7_dp, 375.3608696_dp, 375.3_dp, 377.8274566_dp, 378.2884393_dp, 380.8_dp, 380.7775785_dp, 379.1520179_dp, 378.3_dp, &
    378.3_dp, 378.3_dp, 378.3_dp /)

    N2O_data = (/ &
    284.5_dp, 272.3378378_dp, 267.8378378_dp, 266.5_dp, 268.7_dp, 264.0058824_dp, 263.0838235_dp, 263.0_dp, 267.7916667_dp, 268.0_dp, &
    267.6566372_dp, 260.7035398_dp, 258.3_dp, 265.0830508_dp, 267.0_dp, 263.1904762_dp, 263.0_dp, 265.2597826_dp, 265.6706522_dp, 265.7_dp, &
    266.7_dp, 261.7_dp, 261.7233766_dp, 262.2922078_dp, 262.3_dp, 261.0_dp, 263.7436782_dp, 268.7_dp, 267.1338624_dp, 265.1021164_dp, &
    263.747619_dp, 260.7_dp, 261.7670103_dp, 263.0_dp, 264.7586207_dp, 266.0_dp, 264.6939024_dp, 264.3_dp, 267.1901961_dp, 268.7_dp, &
    268.0114754_dp, 268.0_dp, 267.5487603_dp, 267.3_dp, 265.1456693_dp, 263.7_dp, 267.6655172_dp, 268.3_dp, 268.3_dp, 268.3_dp, &
    267.934375_dp, 267.7_dp, 267.9076923_dp, 268.0_dp, 270.3103448_dp, 271.0_dp, 270.8787879_dp, 267.6464646_dp, 267.0_dp, 265.167619_dp, &
    263.3_dp, 267.272619_dp, 268.0_dp, 268.3_dp, 268.2965318_dp, 268.1300578_dp, 268.0_dp, 264.0_dp, 263.6784946_dp, 262.1204301_dp, &
    261.7_dp, 262.2072165_dp, 262.3_dp, 261.2639098_dp, 261.0_dp, 261.3_dp, 261.262963_dp, 260.3987654_dp, 260.3_dp, 260.9631579_dp, &
    260.9705263_dp, 261.0_dp, 256.3_dp, 256.7046875_dp, 260.69375_dp, 263.7_dp, 262.3923077_dp, 261.7_dp, 261.8933333_dp, 262.0_dp, &
    264.7_dp, 270.8894737_dp, 271.0_dp, 270.4526316_dp, 269.7_dp, 268.7_dp, 263.8108108_dp, 263.3_dp, 265.3017241_dp, 268.7_dp, &
    264.6481013_dp, 259.0_dp, 261.9542857_dp, 267.8_dp, 267.6943662_dp, 263.8915493_dp, 260.3_dp, 258.9320755_dp, 255.3_dp, 256.3992908_dp, &
    260.3_dp, 260.3430108_dp, 260.444086_dp, 260.5365591_dp, 260.6784946_dp, 260.7_dp, 264.9714286_dp, 265.3_dp, 261.4447619_dp, 260.7_dp, &
    266.0_dp, 265.9072072_dp, 265.3504505_dp, 257.3702703_dp, 255.7_dp, 261.7392157_dp, 263.7_dp, 263.7_dp, 263.5407643_dp, 259.1140127_dp, &
    258.7_dp, 256.7253165_dp, 256.3_dp, 258.0875_dp, 260.7_dp, 256.2504762_dp, 254.3_dp, 253.6548387_dp, 253.3_dp, 255.8974684_dp, &
    259.0_dp, 259.3547945_dp, 259.7_dp, 259.3546053_dp, 259.0230263_dp, 259.0_dp, 265.3_dp, 263.2505618_dp, 255.7_dp, 257.1244604_dp, &
    259.3_dp, 258.9910569_dp, 257.5764228_dp, 257.3_dp, 260.5017699_dp, 262.7_dp, 263.46_dp, 263.7_dp, 260.6794872_dp, 258.0_dp, &
    259.3333333_dp, 262.0_dp, 261.3513514_dp, 258.3603604_dp, 258.0_dp, 256.4310345_dp, 256.0_dp, 256.1704545_dp, 256.3977273_dp, 256.7386364_dp, &
    257.0_dp, 258.1612903_dp, 259.0_dp, 257.539823_dp, 256.4336283_dp, 256.2566372_dp, 254.0_dp, 253.9411215_dp, 253.6401869_dp, 253.4766355_dp, &
    253.3_dp, 254.7705882_dp, 255.2411765_dp, 255.8294118_dp, 257.0058824_dp, 259.3_dp, 258.0922078_dp, 257.4298701_dp, 256.3_dp, 256.3317568_dp, &
    256.3405405_dp, 256.3445946_dp, 256.3574324_dp, 256.4_dp, 256.4461538_dp, 256.4923077_dp, 257.0_dp, 260.7295082_dp, 262.0_dp, 260.5978261_dp, &
    259.0_dp, 264.166_dp, 265.3_dp, 265.5627451_dp, 265.7_dp, 260.7714286_dp, 257.7_dp, 257.883871_dp, 258.0_dp, 259.8217822_dp, &
    260.0_dp, 260.2150943_dp, 260.3_dp, 256.0_dp, 257.7230769_dp, 258.0_dp, 256.2914286_dp, 255.7_dp, 259.9741935_dp, 260.7_dp, &
    260.568599_dp, 260.3714976_dp, 260.3_dp, 257.7858974_dp, 255.0_dp, 257.9357798_dp, 260.0_dp, 259.2897196_dp, 259.0_dp, 268.5896552_dp, &
    269.3_dp, 270.0_dp, 270.3912698_dp, 271.7_dp, 271.2535354_dp, 268.3_dp, 268.8449275_dp, 273.0_dp, 269.0704225_dp, 265.2042254_dp, &
    264.0_dp, 264.2567568_dp, 270.0765766_dp, 273.3288288_dp, 273.5_dp, 269.8184466_dp, 268.7_dp, 267.0666667_dp, 266.0_dp, 265.7213115_dp, &
    264.6557377_dp, 264.0819672_dp, 264.0_dp, 265.8316832_dp, 268.8019802_dp, 269.0_dp, 267.6971831_dp, 266.5507042_dp, 265.6908451_dp, 265.3_dp, &
    253.3_dp, 254.0232877_dp, 255.7_dp, 255.380597_dp, 252.1067164_dp, 248.0343284_dp, 245.0_dp, 245.0_dp, 245.0_dp, 244.1968182_dp, &
    242.5904545_dp, 240.2068182_dp, 239.3_dp, 237.9463415_dp, 236.3_dp, 236.6307087_dp, 238.3_dp, 241.1957265_dp, 242.7_dp, 244.3492537_dp, &
    245.3_dp, 247.6611842_dp, 250.5335526_dp, 252.7_dp, 248.8517241_dp, 246.7_dp, 250.6344262_dp, 251.7_dp, 263.8350365_dp, 265.0_dp, &
    266.7_dp, 264.5947368_dp, 259.2263158_dp, 257.7_dp, 259.8025641_dp, 261.7_dp, 261.7_dp, 261.7_dp, 257.4_dp, 258.1671053_dp, &
    268.0_dp, 271.0_dp, 268.7509554_dp, 267.7_dp, 250.5399371_dp, 240.0_dp, 239.3695431_dp, 237.7_dp, 237.0413174_dp, 233.3886228_dp, &
    232.7_dp, 240.0_dp, 241.3474747_dp, 242.3_dp, 226.9580645_dp, 222.0_dp, 234.2340278_dp, 244.3_dp, 229.2066225_dp, 223.0_dp, &
    221.4210526_dp, 220.3_dp, 226.0_dp, 226.0666667_dp, 231.3_dp, 229.8465839_dp, 223.5_dp, 224.8598726_dp, 227.0_dp, 227.5842105_dp, &
    228.0_dp, 232.95_dp, 234.0_dp, 256.0_dp, 223.3_dp, 223.926087_dp, 231.3_dp, 233.1245614_dp, 237.3_dp, 232.9978814_dp, &
    223.0_dp, 234.5656_dp, 241.3_dp, 244.7285714_dp, 253.3_dp, 230.2038168_dp, 203.7_dp, 208.0032258_dp, 217.5_dp, 217.5176887_dp, &
    218.0483491_dp, 218.7558962_dp, 219.6933962_dp, 225.0_dp, 241.5103448_dp, 244.0_dp, 243.347032_dp, 243.0_dp, 242.85_dp, 242.7_dp, &
    231.1270073_dp, 221.7_dp, 238.9108696_dp, 242.0_dp, 242.3198276_dp, 242.3781609_dp, 242.695977_dp, 242.7_dp, 266.5167421_dp, 279.0_dp, &
    257.1518519_dp, 244.3_dp, 233.8079646_dp /)

    ! Initialize found flag
    !found = .FALSE.
    
    ! Find the two indices such that age_data(idx_low) <= time <= age_data(idx_high)
    idx_low = -1
    idx_high = -1

    IF (time <= age_data(1)) THEN
        ! Time is before the first data point
        idx_low = 1
        idx_high = 1
    ELSE IF (time >= age_data(num_points)) THEN
        ! Time is after the last data point
        idx_low = num_points
        idx_high = num_points
    ELSE
        ! Time is within the data range
        DO i = 1, num_points - 1
            IF (age_data(i) <= time .AND. time <= age_data(i+1)) THEN
                idx_low = i
                idx_high = i + 1
                EXIT
            END IF
        END DO
    END IF

    IF (idx_low == -1 .OR. idx_high == -1) THEN
        CALL crash('Time value not within the range of age_data.')
    END IF

    t_low = age_data(idx_low)
    t_high = age_data(idx_high)

    IF (t_high == t_low) THEN
        fraction = 0.0_dp
    ELSE
        fraction = (time - t_low) / (t_high - t_low)
    END IF

    ! Interpolate GHG concentrations
    CO2_current = CO2_data(idx_low) + fraction * (CO2_data(idx_high) - CO2_data(idx_low))
    CH4_current = CH4_data(idx_low) + fraction * (CH4_data(idx_high) - CH4_data(idx_low))
    N2O_current = N2O_data(idx_low) + fraction * (N2O_data(idx_high) - N2O_data(idx_low))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_GHG_concentrations

  SUBROUTINE compute_CO2_forcing(CO2_current, CO2_ref, CO2_relationship, DeltaF_CO2)
    ! Compute radiative forcing due to CO2 using specified relationship

    IMPLICIT NONE

    ! Input variables
    REAL(dp), INTENT(IN)                                  :: CO2_current       ! Current CO2 concentration (ppm)
    REAL(dp), INTENT(IN)                                  :: CO2_ref           ! LGM CO2 concentration (ppm)
    CHARACTER(LEN=256), INTENT(IN)                        :: CO2_relationship  ! Relationship to use ('relationship1', etc.)
    REAL(dp), INTENT(OUT)                                 :: DeltaF_CO2        ! Radiative forcing due to CO2

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_CO2_forcing'
    REAL(dp)                                              :: alpha, beta
    REAL(dp)                                              :: g_current, g_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute DeltaF_CO2 based on selected relationship
    SELECT CASE (CO2_relationship)

    CASE ('relationship1')
      ! ∆F= α ln(C/C0)
      alpha = 5.35_dp
      DeltaF_CO2 = alpha * LOG(CO2_current / CO2_ref)
      
    CASE ('relationship2')
      ! ∆F= α ln(C/C0) + β(√C − √C0)
      alpha = 4.841_dp
      beta = 0.0906_dp
      DeltaF_CO2 = alpha * LOG(CO2_current / CO2_ref) + beta * (SQRT(CO2_current) - SQRT(CO2_ref))
    
    CASE ('relationship3')
      ! ∆F= α(g(C)–g(C0)), where g(C) = ln(1 + 1.2C + 0.005C² + 1.4 × 10⁻⁶C³)
      alpha = 3.35_dp
      g_current = LOG(1.0_dp + 1.2_dp*CO2_current + 0.005_dp*CO2_current**2 + 1.4e-6_dp*CO2_current**3)
      g_ref = LOG(1.0_dp + 1.2_dp*CO2_ref + 0.005_dp*CO2_ref**2 + 1.4e-6_dp*CO2_ref**3)
      DeltaF_CO2 = alpha * g_current - g_ref 
    
    CASE DEFAULT
      CALL crash('Unknown CO2 relationship: ' // CO2_relationship)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_CO2_forcing

  FUNCTION f_overlap(M, N) RESULT(f)
    ! Function radiative forcing CH4 and N2O
    ! f(M,N) = 0.47 * ln[1 + 2.01e-5 * (M*N)^0.75 + 5.31e-15 * M * (M*N)^1.52]

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: M, N

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'f_overlap'
    REAL(dp)                                              :: f
    REAL(dp)                                              :: term1, term2, argument

    ! Add routine to path
    CALL init_routine( routine_name)

    term1 = 2.01e-5_dp * (M * N)**0.75_dp
    term2 = 5.31e-15_dp * M * (M * N)**1.52_dp
    argument = 1.0_dp + term1 + term2
    f = 0.47_dp * LOG(argument)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END FUNCTION f_overlap

  SUBROUTINE compute_CH4_forcing(M, M0, N0, DeltaF_CH4)
    ! Compute radiative forcing contributions from CH4 and N2O
    
    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: M            ! CH4 concentration (ppb)
    REAL(dp), INTENT(IN)                                  :: M0           ! LGM CH4 concentration (ppb)
    REAL(dp), INTENT(IN)                                  :: N0           ! LGM N2O concentration (ppb)
    REAL(dp), INTENT(OUT)                                 :: DeltaF_CH4   ! Radiative forcing due to CH4

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_CH4_N2O_forcing'
    REAL(dp), PARAMETER                                   :: alpha_CH4 = 0.036_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute DeltaF_CH4
    ! ∆F= α(√M–√M0)–(f(M,N0)–f(M0,N0))
    DeltaF_CH4 = alpha_CH4 * (SQRT(M) - SQRT(M0)) - (f_overlap(M, N0) - f_overlap(M0, N0))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_CH4_forcing

  SUBROUTINE compute_N2O_forcing(N, N0, M0, DeltaF_N2O)
    ! Compute radiative forcing contributions from N2O
    
    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: M0           ! LGM CH4 concentration (ppb)
    REAL(dp), INTENT(IN)                                  :: N            ! N2O concentration (ppb)
    REAL(dp), INTENT(IN)                                  :: N0           ! LGM N2O concentration (ppb)
    REAL(dp), INTENT(OUT)                                 :: DeltaF_N2O   ! Radiative forcing due to N2O

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_CH4_N2O_forcing'
    REAL(dp), PARAMETER                                   :: alpha_N2O = 0.12_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Compute DeltaF_N2O
    ! ∆F= α(√N–√N0)–(f(M0,N)–f(M0,N 0))
    DeltaF_N2O = alpha_N2O * (SQRT(N) - SQRT(N0)) - (f_overlap(M0, N) - f_overlap(M0, N0))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_N2O_forcing

  !SUBROUTINE polynomial_time_interpolation(mesh, ocean, time, num_timeframes, times, weights)
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
    !CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'polynomial_time_interpolation'
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
          !IF (ALLOCATED(ocean%matrix%timeframe0%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe0%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe0%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(ocean%matrix%timeframe1%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe1%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe1%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(ocean%matrix%timeframe2%T)) THEN
              !n = n + 1
              !IF (n <= num_timeframes) THEN
                  !ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe2%T(i, j)
                  !ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe2%S(i, j)
              !END IF
          !END IF
          !IF (ALLOCATED(ocean%matrix%timeframe3%T)) THEN
            !n = n + 1
            !IF (n <= num_timeframes) THEN
                !ocean%T(i, j) = ocean%T(i, j) + weights(n) * ocean%matrix%timeframe3%T(i, j)
                !ocean%S(i, j) = ocean%S(i, j) + weights(n) * ocean%matrix%timeframe3%S(i, j)
            !END IF
          !END IF
          ! Limited polynomial order due to computational time
      !END DO
    !END DO

    ! Finalise routine path
    !CALL finalise_routine( routine_name)
  
  !END SUBROUTINE polynomial_time_interpolation

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

    ! Minimum required amount of timeframes for each linear time interpolation method
    !IF (TRIM(C%choice_ocean_model_matrix) == 'linear_time') THEN
      !required_timeframes = 2
    !ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial_time') THEN
        !required_timeframes = 3 
    !ELSE
    !END IF

    ! Perform time interpolation
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear_time') THEN
      CALL linear_time_interpolation(mesh, ocean, matrix, time)
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial_time') THEN
        CALL crash('Polynomial interpolation not implemented yet')
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'GHG_radiative_based') THEN
      CALL interpolation_with_GHG_radiative(mesh, ocean, matrix, time)
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'GHG_based') THEN
      CALL interpolation_with_GHG_basic(mesh, ocean, matrix, time)
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'insolation') THEN
      CALL interpolation_with_insolation(mesh, ocean, matrix, time)
    ELSE
        CALL crash('Unknown choice_ocean_model_matrix' // TRIM(C%choice_ocean_model_matrix))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE run_ocean_model_matrix
  
  SUBROUTINE initialise_ocean_model_matrix( mesh, ocean, region_name)
    ! Initialise the ocean matrix model
  
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_matrix'
    CHARACTER(LEN=256)                                    :: filename1, filename2
    INTEGER                                               :: i, j
    !INTEGER                                               :: ndepth
    !REAL(dp), DIMENSION(:), ALLOCATABLE                   :: depth

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising matrix ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_matrix),'light blue') // '"...'

    ! Start and ending of simulation
    ocean%matrix%t0 = REAL(C%start_time_of_run, dp)     ! LGM
    ocean%matrix%t1 = REAL(C%end_time_of_run, dp)       ! PI

    ! Possibility to hardcode the depth, should then be called in read_field_from_file_3D_ocean command
    !ndepth = 11
    !ALLOCATE(depth(ndepth))
    !depth = (/0.0_dp, 150.0_dp, 300.0_dp, 450.0_dp, 600.0_dp, 750.0_dp, 900.0_dp, &
              !1050.0_dp, 1200.0_dp, 1350.0_dp, 1500.0_dp/)

    ! Allocate memory for timeframes' T and S array if not already allocated
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe0%T)) THEN
      ALLOCATE(ocean%matrix%timeframe0%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ALLOCATE(ocean%matrix%timeframe0%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
    END IF
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe1%T)) THEN
      ALLOCATE(ocean%matrix%timeframe1%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ALLOCATE(ocean%matrix%timeframe1%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))      
    END IF

    ! Construct filenames for the two ocean snapshots
    filename1 = TRIM(C%filename_ocean_matrix_base1) ! LGM
    filename2 = TRIM(C%filename_ocean_matrix_base2) ! PI

    ! Read the ocean snapshots
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_T_ocean, mesh, ocean%matrix%timeframe0%T)
    CALL read_field_from_file_3D_ocean(filename1, field_name_options_S_ocean, mesh, ocean%matrix%timeframe0%S)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_T_ocean, mesh, ocean%matrix%timeframe1%T)
    CALL read_field_from_file_3D_ocean(filename2, field_name_options_S_ocean, mesh, ocean%matrix%timeframe1%S)

    ! Ensure correct model choice
    IF (TRIM(C%choice_ocean_model_matrix) == 'linear_time') THEN
      ! Check
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'polynomial_time') THEN
      CALL crash('Polynomial interpolation not implemented yet')
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'GHG_radiative_based') THEN
      ! Check
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'GHG_based') THEN
      ! Check
    ELSE IF (TRIM(C%choice_ocean_model_matrix) == 'insolation') THEN
      ! Check
    ELSE
      CALL crash('Unknown choice_ocean_model_matrix' // TRIM(C%choice_ocean_model_matrix))
    END IF

    ! Prescribe initial ocean state
    DO i = mesh%vi1, mesh%vi2
      DO j = 1, C%nz_ocean
        ocean%T(i,j) = ocean%matrix%timeframe0%T(i,j)
        ocean%S(i,j) = ocean%matrix%timeframe0%S(i,j)
      END DO
    END DO
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_ocean_model_matrix

END MODULE ocean_matrix