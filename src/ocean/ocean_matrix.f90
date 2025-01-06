MODULE ocean_matrix

    ! Timeframe0 = PI
    ! Timeframe1 = LGM

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
    USE region_types                                           , ONLY: type_model_region
    USE netcdf_input                                           , ONLY: read_field_from_file_3D_ocean, read_field_from_file_2D
    USE netcdf_basic                                           , ONLY: field_name_options_T_ocean, field_name_options_S_ocean
    USE ocean_utilities                                        , ONLY: debug_ocean_matrix_state, initialise_ocean_vertical_grid
    USE mesh_utilities                                         , ONLY: extrapolate_Gaussian
    USE reference_geometry_types                               , ONLY: type_reference_geometry
    USE grid_types                                             , ONLY: type_grid
    USE mesh_data_smoothing                                    , ONLY: smooth_Gaussian_2D

    IMPLICIT NONE
  
  CONTAINS
  
  ! ===== Main routines =====
  ! =========================
  
  SUBROUTINE get_sea_level_values(time, sea_level_current)
    ! Get sea level values at the given time, limited to 30 ka for now, can be extended to 798 ka.
    ! Spratt, R. M., & Lisiecki, L. E. (2016). A Late Pleistocene sea level stack. Climate of the Past, 12(4), 1079-1092.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: time
    REAL(dp), INTENT(OUT)                                 :: sea_level_current

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'get_sea_level_values'
    INTEGER, PARAMETER                                    :: num_points = 31
    REAL(dp), DIMENSION(num_points)                       :: age_data, sea_level_data
    INTEGER                                               :: i
    !LOGICAL                                               :: found
    INTEGER                                               :: idx_low, idx_high
    REAL(dp)                                              :: t_low, t_high, fraction

    ! Add routine to path
    CALL init_routine( routine_name)

    ! FIX
    ! Hardcoded values for now, should be reading in from csv / txt / dat
    ! https://stackoverflow.com/questions/8828377/reading-data-from-txt-file-in-fortran
    
    ! Initialize data arrays
    age_data = (/ &
    0._dp, 1000._dp	, 2000._dp, 3000._dp, 4000._dp, 5000._dp, 6000._dp, 7000._dp, 8000._dp, 9000._dp, 10000._dp, &
    11000._dp, 12000._dp, 13000._dp, 14000._dp, 15000._dp, 16000._dp, 17000._dp, 18000._dp, 19000._dp, 20000._dp, &
    21000._dp, 22000._dp, 23000._dp, 24000._dp, 25000._dp, 26000._dp, 27000._dp, 28000._dp, 29000._dp, 30000._dp /)
  
    sea_level_data = (/ &
    8.49_dp, 7.63_dp, 4.01_dp, 4.35_dp, 3.13_dp, 0._dp, -4.01_dp, -6.11_dp, -9.09_dp, -15.83_dp, -24.59_dp, -35.85_dp, -51.05_dp, &
    -66.3_dp, -76.64_dp, -86.57_dp, -98.06_dp, -107.3_dp, -113.01_dp, -116.68_dp, -117.56_dp, -120.01_dp, -125.82_dp, -128.72_dp, &
    -130.0_dp, -126.89_dp, -122.4_dp, -118.28_dp, -115.19_dp, -110.87_dp, -105.78_dp /)

    ! Initialize found flag
    !found = .FALSE.
    
    ! Find the two indices such that age_data(idx_low) <= time <= age_data(idx_high)
    idx_low = -1
    idx_high = -1

    ! Edge cases
    IF (time <= age_data(1)) THEN
        ! Time is before the first data point
        idx_low = 1
        idx_high = 1
    ELSE IF (time >= age_data(num_points)) THEN
        ! Time is after the last data point
        idx_low = num_points
        idx_high = num_points
    ! Main data search
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

    ! Interpolate sea level values to match to model runtime
    sea_level_current = sea_level_data(idx_low) + fraction * (sea_level_data(idx_high) - sea_level_data(idx_low))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_sea_level_values

  SUBROUTINE interpolation_with_d18O(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on GHG

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_d18O'
    REAL(dp)                                              :: w_d18O
    INTEGER                                               :: i, j, vi, k
    REAL(dp)                                              :: d18O_current, d18O_PI, d18O_LGM
    REAL(dp)                                              :: scale_T, scale_S

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Retrieve scaling factors from configuration
    scale_T = C%scale_d18O_temperature
    scale_S = C%scale_d18O_salinity

    ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
    CALL get_d18O_values(time, d18O_current)    ! d18O during runtime
    CALL get_d18O_values(0.0_dp, d18O_PI)     ! PI
    CALL get_d18O_values(21000.0_dp, d18O_LGM)  ! LGM

    print *, "d18O current = ", d18O_current

    ! Compute weight based on d18O value ratios
    w_d18O = (d18O_current - d18O_LGM) / (d18O_PI - d18O_LGM)

    ! Clamp between cutoff values if enabled
    IF (C%clamp_weights) THEN
      w_d18O = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_d18O))
    END IF

    print *, "Interpolation weight calculated = ", w_d18O

    ! Apply interpolation using w_d18O and scaling
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, C%nz_ocean
        ocean%T(vi,k) = w_d18O * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_d18O) * ocean%matrix%timeframe1%T(vi,k) + scale_T * d18O_current
        ocean%S(vi,k) = w_d18O * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_d18O) * ocean%matrix%timeframe1%S(vi,k) + scale_S * d18O_current
      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_d18O

  SUBROUTINE get_d18O_values(time, d18O_current)
    ! Get d18O values at the given time
    ! Stenni, B., et al. (2006), EPICA Dome C Stable Isotope Data to 44.8 KYrBP.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: time
    REAL(dp), INTENT(OUT)                                 :: d18O_current

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'get_d18O_values'
    INTEGER, PARAMETER                                    :: num_points = 914
    REAL(dp), DIMENSION(num_points)                       :: age_data, delta18O_data
    INTEGER                                               :: i
    !LOGICAL                                               :: found
    INTEGER                                               :: idx_low, idx_high
    REAL(dp)                                              :: t_low, t_high, fraction

    ! Add routine to path
    CALL init_routine( routine_name)

    ! FIX
    ! Hardcoded values for now, should be reading in from csv / txt
    ! https://stackoverflow.com/questions/8828377/reading-data-from-txt-file-in-fortran
    
    ! Initialize data arrays
    age_data = (/ &
    42.5_dp,     51.15_dp,     60.08_dp,     69.36_dp,     78.58_dp,     87.77_dp,     97.06_dp,     106.43_dp,     116.11_dp,     125.83_dp,     &
    135.47_dp,     145.35_dp,     155.61_dp,     165.93_dp,     176.29_dp,     186.69_dp,     197.1_dp,     207.53_dp,     218.05_dp,     228.94_dp,     &
    239.97_dp,     251.22_dp,     262.57_dp,     274.02_dp,     285.58_dp,     297.21_dp,     309.07_dp,     321.01_dp,     333.03_dp,     345.1_dp,     &
    357.25_dp,     369.47_dp,     381.79_dp,     394.22_dp,     406.75_dp,     419.28_dp,     431.81_dp,     444.41_dp,     457.05_dp,     469.8_dp,     &
    482.6_dp,     495.46_dp,     508.31_dp,     521.13_dp,     534.03_dp,     546.96_dp,     559.81_dp,     572.64_dp,     585.5_dp,     598.38_dp,     &
    611.3_dp,     624.29_dp,     637.4_dp,     650.6_dp,     663.89_dp,     677.33_dp,     690.84_dp,     704.44_dp,     718.03_dp,     731.51_dp,     &
    745.07_dp,     758.82_dp,     772.25_dp,     785.27_dp,     798.34_dp,     811.44_dp,     824.63_dp,     837.83_dp,     850.99_dp,     864.17_dp,     &
    877.48_dp,     890.81_dp,     904.19_dp,     917.7_dp,     931.33_dp,     945.09_dp,     958.89_dp,     972.67_dp,     986.48_dp,     1000.49_dp,     &
    1014.52_dp,     1028.55_dp,     1042.56_dp,     1056.55_dp,     1070.59_dp,     1084.65_dp,     1098.58_dp,     1112.49_dp,     1126.5_dp,     1140.56_dp,     &
    1154.68_dp,     1168.75_dp,     1182.72_dp,     1196.69_dp,     1210.64_dp,     1224.57_dp,     1238.51_dp,     1252.64_dp,     1266.81_dp,     1281.03_dp,     &
    1295.24_dp,     1309.44_dp,     1323.77_dp,     1338.25_dp,     1352.78_dp,     1367.34_dp,     1381.92_dp,     1396.54_dp,     1411.37_dp,     1426.34_dp,     &
    1441.6_dp,     1456.89_dp,     1472.17_dp,     1487.48_dp,     1502.84_dp,     1518.4_dp,     1534.05_dp,     1549.76_dp,     1565.49_dp,     1581.25_dp,     &
    1597.07_dp,     1613.02_dp,     1629.07_dp,     1645.24_dp,     1661.36_dp,     1677.42_dp,     1693.33_dp,     1709.29_dp,     1725.75_dp,     1742.27_dp,     &
    1758.7_dp,     1775.06_dp,     1791.33_dp,     1807.75_dp,     1824.29_dp,     1840.75_dp,     1857.18_dp,     1873.65_dp,     1890.15_dp,     1906.71_dp,     &
    1923.21_dp,     1939.55_dp,     1956.13_dp,     1973.01_dp,     1989.91_dp,     2006.82_dp,     2023.8_dp,     2040.81_dp,     2057.91_dp,     2075.06_dp,     &
    2092.34_dp,     2109.61_dp,     2126.85_dp,     2144.21_dp,     2161.68_dp,     2179.29_dp,     2196.95_dp,     2232.37_dp,     2250.26_dp,     2268.33_dp,     &
    2286.75_dp,     2305.2_dp,     2323.67_dp,     2342.16_dp,     2360.64_dp,     2379.11_dp,     2397.55_dp,     2415.9_dp,     2434.21_dp,     2452.44_dp,     &
    2470.63_dp,     2488.72_dp,     2506.69_dp,     2524.56_dp,     2542.31_dp,     2560.01_dp,     2577.64_dp,     2595.25_dp,     2612.79_dp,     2630.34_dp,     &
    2647.89_dp,     2665.39_dp,     2682.81_dp,     2700.1_dp,     2717.31_dp,     2734.42_dp,     2751.49_dp,     2768.53_dp,     2785.61_dp,     2802.81_dp,     &
    2820.15_dp,     2837.69_dp,     2855.29_dp,     2872.92_dp,     2890.5_dp,     2908.1_dp,     2925.99_dp,     2943.98_dp,     2962.05_dp,     2980.13_dp,     &
    2998.19_dp,     3016.33_dp,     3034.58_dp,     3052.76_dp,     3070.92_dp,     3107.93_dp,     3126.72_dp,     3145.59_dp,     3164.59_dp,     3183.64_dp,     &
    3202.78_dp,     3221.95_dp,     3241.15_dp,     3260.36_dp,     3279.56_dp,     3298.74_dp,     3317.89_dp,     3336.98_dp,     3356.02_dp,     3374.98_dp,     &
    3393.86_dp,     3412.65_dp,     3431.35_dp,     3449.98_dp,     3468.51_dp,     3486.99_dp,     3505.38_dp,     3523.74_dp,     3542.02_dp,     3560.27_dp,     &
    3578.48_dp,     3596.67_dp,     3614.85_dp,     3633.06_dp,     3651.28_dp,     3669.58_dp,     3687.92_dp,     3706.39_dp,     3724.93_dp,     3743.6_dp,     &
    3762.37_dp,     3781.26_dp,     3800.26_dp,     3819.33_dp,     3838.48_dp,     3857.68_dp,     3876.92_dp,     3896.17_dp,     3915.43_dp,     3934.67_dp,     &
    3953.9_dp,     3973.1_dp,     3992.27_dp,     4011.42_dp,     4030.55_dp,     4049.67_dp,     4068.8_dp,     4087.93_dp,     4107.06_dp,     4126.2_dp,     &
    4145.34_dp,     4164.48_dp,     4183.62_dp,     4202.77_dp,     4221.94_dp,     4241.13_dp,     4260.38_dp,     4279.66_dp,     4298.99_dp,     4318.33_dp,     &
    4337.71_dp,     4357.1_dp,     4376.5_dp,     4395.91_dp,     4415.33_dp,     4434.78_dp,     4454.26_dp,     4473.8_dp,     4493.37_dp,     4512.98_dp,     &
    4532.62_dp,     4552.3_dp,     4572.01_dp,     4591.74_dp,     4611.52_dp,     4631.33_dp,     4651.2_dp,     4671.08_dp,     4691.02_dp,     4710.96_dp,     &
    4730.9_dp,     4750.8_dp,     4770.67_dp,     4790.45_dp,     4810.18_dp,     4829.82_dp,     4849.43_dp,     4869.0_dp,     4888.56_dp,     4908.12_dp,     &
    4927.69_dp,     4947.29_dp,     4966.9_dp,     4986.51_dp,     5006.11_dp,     5025.7_dp,     5045.25_dp,     5064.78_dp,     5084.27_dp,     5103.75_dp,     &
    5123.23_dp,     5142.72_dp,     5162.23_dp,     5181.79_dp,     5201.38_dp,     5221.02_dp,     5240.7_dp,     5260.47_dp,     5280.32_dp,     5300.29_dp,     &
    5320.39_dp,     5340.63_dp,     5361.0_dp,     5381.47_dp,     5402.0_dp,     5422.56_dp,     5443.07_dp,     5463.52_dp,     5483.83_dp,     5504.04_dp,     &
    5524.09_dp,     5544.01_dp,     5563.83_dp,     5583.52_dp,     5603.13_dp,     5622.63_dp,     5642.06_dp,     5661.39_dp,     5680.66_dp,     5699.86_dp,     &
    5719.05_dp,     5738.22_dp,     5757.42_dp,     5776.66_dp,     5796.01_dp,     5815.42_dp,     5834.98_dp,     5854.61_dp,     5874.36_dp,     5894.17_dp,     &
    5914.04_dp,     5933.93_dp,     5953.83_dp,     5973.72_dp,     5993.6_dp,     6013.48_dp,     6033.36_dp,     6053.26_dp,     6073.18_dp,     6093.14_dp,     &
    6113.14_dp,     6133.19_dp,     6153.31_dp,     6173.47_dp,     6193.72_dp,     6214.0_dp,     6234.37_dp,     6254.77_dp,     6275.23_dp,     6295.7_dp,     &
    6316.2_dp,     6336.69_dp,     6357.18_dp,     6377.66_dp,     6398.15_dp,     6418.67_dp,     6439.22_dp,     6459.83_dp,     6480.5_dp,     6501.24_dp,     &
    6522.05_dp,     6542.92_dp,     6563.83_dp,     6584.76_dp,     6605.66_dp,     6626.56_dp,     6647.39_dp,     6668.2_dp,     6688.94_dp,     6709.63_dp,     &
    6730.25_dp,     6750.8_dp,     6771.29_dp,     6791.69_dp,     6812.05_dp,     6832.39_dp,     6852.75_dp,     6873.18_dp,     6893.69_dp,     6914.33_dp,     &
    6935.11_dp,     6956.03_dp,     6977.14_dp,     6998.38_dp,     7019.83_dp,     7041.38_dp,     7063.14_dp,     7084.98_dp,     7106.96_dp,     7128.98_dp,     &
    7151.04_dp,     7172.44_dp,     7193.25_dp,     7213.1_dp,     7232.52_dp,     7251.85_dp,     7271.16_dp,     7290.43_dp,     7309.69_dp,     7328.92_dp,     &
    7348.12_dp,     7367.28_dp,     7386.35_dp,     7405.36_dp,     7424.21_dp,     7443.0_dp,     7461.62_dp,     7480.17_dp,     7498.61_dp,     7517.02_dp,     &
    7535.39_dp,     7553.76_dp,     7572.13_dp,     7590.52_dp,     7608.92_dp,     7627.34_dp,     7645.79_dp,     7664.27_dp,     7682.78_dp,     7701.34_dp,     &
    7719.95_dp,     7738.62_dp,     7757.34_dp,     7776.11_dp,     7794.96_dp,     7813.86_dp,     7832.89_dp,     7852.0_dp,     7871.28_dp,     7890.67_dp,     &
    7910.21_dp,     7929.87_dp,     7949.63_dp,     7969.5_dp,     7989.43_dp,     8009.43_dp,     8029.45_dp,     8049.51_dp,     8069.55_dp,     8089.58_dp,     &
    8109.57_dp,     8129.52_dp,     8149.4_dp,     8169.25_dp,     8189.04_dp,     8208.8_dp,     8228.51_dp,     8248.17_dp,     8267.75_dp,     8287.26_dp,     &
    8306.67_dp,     8325.97_dp,     8345.19_dp,     8364.28_dp,     8383.31_dp,     8402.23_dp,     8421.12_dp,     8439.95_dp,     8458.78_dp,     8477.61_dp,     &
    8496.47_dp,     8515.36_dp,     8534.3_dp,     8553.28_dp,     8572.28_dp,     8591.3_dp,     8610.32_dp,     8629.32_dp,     8648.27_dp,     8667.18_dp,     &
    8686.02_dp,     8704.79_dp,     8723.49_dp,     8742.12_dp,     8760.71_dp,     8779.25_dp,     8797.78_dp,     8816.3_dp,     8834.84_dp,     8853.41_dp,     &
    8872.0_dp,     8890.63_dp,     8909.29_dp,     8927.97_dp,     8946.65_dp,     8965.31_dp,     8983.91_dp,     9002.46_dp,     9020.9_dp,     9039.27_dp,     &
    9057.53_dp,     9075.71_dp,     9093.82_dp,     9111.85_dp,     9129.84_dp,     9147.75_dp,     9165.62_dp,     9183.42_dp,     9201.21_dp,     9218.98_dp,     &
    9236.74_dp,     9272.3_dp,     9290.1_dp,     9307.9_dp,     9325.69_dp,     9343.44_dp,     9361.14_dp,     9378.75_dp,     9396.29_dp,     9413.73_dp,     &
    9431.08_dp,     9448.36_dp,     9465.55_dp,     9482.69_dp,     9499.78_dp,     9516.87_dp,     9533.97_dp,     9551.1_dp,     9568.32_dp,     9585.61_dp,     &
    9603.0_dp,     9620.48_dp,     9638.02_dp,     9655.62_dp,     9673.24_dp,     9690.88_dp,     9708.51_dp,     9726.12_dp,     9743.7_dp,     9761.24_dp,     &
    9778.72_dp,     9796.15_dp,     9813.51_dp,     9830.82_dp,     9848.08_dp,     9865.31_dp,     9882.5_dp,     9899.69_dp,     9916.89_dp,     9934.13_dp,     &
    9951.42_dp,     9968.8_dp,     9986.25_dp,     10003.82_dp,     10021.46_dp,     10039.22_dp,     10057.02_dp,     10074.91_dp,     10092.84_dp,     10110.83_dp,     &
    10128.83_dp,     10146.85_dp,     10164.82_dp,     10182.76_dp,     10200.59_dp,     10218.35_dp,     10235.97_dp,     10253.52_dp,     10270.98_dp,     10288.43_dp,     &
    10305.88_dp,     10323.41_dp,     10341.0_dp,     10358.75_dp,     10376.57_dp,     10394.57_dp,     10412.64_dp,     10430.85_dp,     10449.13_dp,     10467.52_dp,     &
    10485.99_dp,     10504.54_dp,     10523.16_dp,     10541.81_dp,     10560.5_dp,     10579.2_dp,     10597.92_dp,     10616.65_dp,     10635.39_dp,     10654.16_dp,     &
    10672.94_dp,     10691.75_dp,     10710.57_dp,     10729.42_dp,     10748.29_dp,     10767.17_dp,     10786.07_dp,     10804.99_dp,     10823.95_dp,     10842.95_dp,     &
    10862.01_dp,     10881.1_dp,     10900.21_dp,     10919.32_dp,     10938.39_dp,     10957.44_dp,     10976.42_dp,     10995.37_dp,     11014.28_dp,     11033.17_dp,     &
    11052.05_dp,     11070.91_dp,     11089.76_dp,     11108.59_dp,     11127.42_dp,     11146.27_dp,     11165.12_dp,     11184.0_dp,     11202.88_dp,     11221.76_dp,     &
    11240.62_dp,     11259.48_dp,     11278.36_dp,     11297.28_dp,     11316.3_dp,     11335.37_dp,     11354.59_dp,     11373.87_dp,     11393.26_dp,     11412.7_dp,     &
    11432.2_dp,     11451.77_dp,     11471.37_dp,     11491.08_dp,     11510.86_dp,     11530.82_dp,     11550.9_dp,     11571.23_dp,     11591.76_dp,     11612.58_dp,     &
    11633.69_dp,     11655.05_dp,     11676.74_dp,     11698.61_dp,     11720.81_dp,     11743.14_dp,     11765.77_dp,     11788.55_dp,     11811.57_dp,     11834.78_dp,     &
    11858.22_dp,     11881.91_dp,     11905.78_dp,     11929.93_dp,     11954.21_dp,     11978.72_dp,     12003.33_dp,     12028.11_dp,     12053.0_dp,     12078.05_dp,     &
    12103.29_dp,     12128.7_dp,     12154.38_dp,     12180.22_dp,     12206.38_dp,     12232.68_dp,     12259.31_dp,     12286.11_dp,     12313.24_dp,     12340.62_dp,     &
    12368.27_dp,     12396.19_dp,     12424.29_dp,     12452.57_dp,     12480.93_dp,     12509.36_dp,     12537.81_dp,     12566.25_dp,     12594.67_dp,     12623.06_dp,     &
    12651.38_dp,     12679.64_dp,     12707.81_dp,     12735.92_dp,     12763.94_dp,     12791.91_dp,     12819.78_dp,     12847.61_dp,     12875.34_dp,     12903.0_dp,     &
    12930.58_dp,     12958.07_dp,     12985.48_dp,     13012.78_dp,     13040.04_dp,     13067.22_dp,     13094.4_dp,     13121.6_dp,     13148.85_dp,     13176.18_dp,     &
    13203.6_dp,     13231.11_dp,     13258.71_dp,     13286.37_dp,     13314.12_dp,     13341.92_dp,     13369.8_dp,     13397.74_dp,     13425.79_dp,     13453.93_dp,     &
    13482.17_dp,     13510.44_dp,     13538.74_dp,     13566.98_dp,     13595.17_dp,     13623.22_dp,     13651.2_dp,     13679.0_dp,     13706.68_dp,     13734.2_dp,     &
    13761.55_dp,     13788.75_dp,     13815.75_dp,     13842.64_dp,     13869.36_dp,     13896.0_dp,     13922.49_dp,     13948.89_dp,     13975.12_dp,     14001.2_dp,     &
    14027.1_dp,     14052.8_dp,     14078.36_dp,     14103.75_dp,     14129.08_dp,     14154.39_dp,     14179.77_dp,     14205.33_dp,     14231.11_dp,     14257.22_dp,     &
    14283.72_dp,     14310.57_dp,     14337.92_dp,     14365.55_dp,     14393.72_dp,     14422.13_dp,     14451.02_dp,     14480.1_dp,     14509.52_dp,     14539.07_dp,     &
    14568.76_dp,     14598.46_dp,     14628.16_dp,     14657.77_dp,     14687.32_dp,     14716.77_dp,     14746.17_dp,     14775.53_dp,     14804.9_dp,     14834.3_dp,     &
    14863.78_dp,     14893.35_dp,     14923.11_dp,     14952.99_dp,     14983.13_dp,     15013.42_dp,     15044.04_dp,     15074.87_dp,     15106.08_dp,     15137.57_dp,     &
    15169.4_dp,     15201.58_dp,     15233.99_dp,     15266.7_dp,     15299.53_dp,     15332.56_dp,     15365.66_dp,     15398.87_dp,     15432.13_dp,     15465.46_dp,     &
    15498.86_dp,     15532.31_dp,     15565.85_dp,     15599.46_dp,     15633.21_dp,     15667.05_dp,     15701.14_dp,     15735.38_dp,     15769.92_dp,     15804.71_dp,     &
    15839.78_dp,     15875.2_dp,     15910.84_dp,     15946.86_dp,     15983.05_dp,     16019.58_dp,     16056.24_dp,     16093.15_dp,     16130.2_dp,     16167.43_dp,     &
    16204.86_dp,     16242.45_dp,     16280.33_dp,     16318.38_dp,     16356.81_dp,     16395.45_dp,     16434.55_dp,     16473.92_dp,     16513.78_dp,     16554.02_dp,     &
    16594.67_dp,     16635.82_dp,     16677.29_dp,     16719.32_dp,     16761.63_dp,     16804.49_dp,     16847.57_dp,     16891.08_dp,     16934.79_dp,     16978.78_dp,     &
    17022.98_dp,     17067.34_dp,     17111.99_dp,     17156.8_dp,     17202.02_dp,     17247.46_dp,     17293.42_dp,     17339.66_dp,     17386.41_dp,     17433.53_dp,     &
    17481.07_dp,     17529.01_dp,     17577.22_dp,     17625.77_dp,     17674.45_dp,     17723.34_dp,     17772.3_dp,     17821.39_dp,     17870.56_dp,     17919.82_dp,     &
    17969.18_dp,     18018.62_dp,     18068.15_dp,     18117.72_dp,     18167.35_dp,     18217.01_dp,     18266.72_dp,     18316.47_dp,     18366.28_dp,     18416.14_dp,     &
    18466.03_dp,     18515.94_dp,     18565.85_dp,     18615.72_dp,     18665.58_dp,     18715.41_dp,     18765.23_dp,     18815.05_dp,     18864.88_dp,     18914.74_dp,     &
    18964.61_dp,     19014.49_dp,     19064.38_dp,     19114.27_dp,     19164.14_dp,     19213.98_dp,     19263.75_dp,     19313.46_dp,     19363.08_dp,     19412.63_dp,     &
    19462.1_dp,     19511.57_dp,     19561.05_dp,     19610.68_dp,     19660.38_dp,     19710.33_dp,     19760.39_dp,     19810.64_dp,     19860.96_dp,     19911.34_dp,     &
    19961.65_dp,     20011.89_dp,     20061.92_dp,     20111.83_dp,     20161.5_dp,     20211.08_dp,     20260.48_dp,     20309.8_dp,     20358.99_dp,     20408.1_dp,     &
    20457.11_dp,     20506.06_dp,     20554.98_dp,     20603.89_dp,     20652.81_dp,     20701.77_dp,     20750.75_dp,     20799.76_dp,     20848.78_dp,     20897.81_dp,     &
    20946.88_dp,     20995.98_dp,     21045.16_dp,     21094.39_dp /)
  
    delta18O_data = (/ &
    -50.12_dp,     -49.16_dp,     -48.14_dp,     -50.35_dp,     -50.94_dp,     -50.6_dp,     -50.52_dp,     -51.07_dp,     -49.97_dp,     -50.89_dp,     &
    -51.6_dp,     -50.45_dp,     -50.81_dp,     -49.98_dp,     -50.62_dp,     -51.07_dp,     -51.72_dp,     -50.11_dp,     -52.51_dp,     -50.22_dp,     &
    -49.92_dp,     -49.86_dp,     -50.44_dp,     -50.7_dp,     -51.7_dp,     -50.98_dp,     -50.72_dp,     -51.68_dp,     -51.17_dp,     -51.26_dp,     &
    -50.64_dp,     -51.5_dp,     -51.56_dp,     -51.89_dp,     -51.0_dp,     -50.55_dp,     -50.58_dp,     -51.93_dp,     -50.5_dp,     -50.49_dp,     &
    -51.06_dp,     -51.24_dp,     -50.23_dp,     -52.52_dp,     -51.03_dp,     -50.94_dp,     -49.54_dp,     -50.14_dp,     -51.18_dp,     -49.53_dp,     &
    -51.2_dp,     -51.18_dp,     -49.33_dp,     -49.82_dp,     -50.43_dp,     -51.31_dp,     -51.06_dp,     -51.1_dp,     -51.02_dp,     -51.49_dp,     &
    -51.68_dp,     -51.58_dp,     -50.7_dp,     -50.82_dp,     -50.44_dp,     -50.65_dp,     -50.24_dp,     -49.93_dp,     -49.7_dp,     -49.57_dp,     &
    -51.58_dp,     -50.27_dp,     -51.79_dp,     -49.99_dp,     -51.45_dp,     -51.04_dp,     -51.16_dp,     -51.28_dp,     -50.41_dp,     -51.06_dp,     &
    -50.98_dp,     -51.55_dp,     -50.61_dp,     -49.76_dp,     -50.36_dp,     -51.06_dp,     -51.2_dp,     -50.17_dp,     -51.09_dp,     -49.6_dp,     &
    -51.35_dp,     -50.9_dp,     -49.91_dp,     -49.89_dp,     -50.1_dp,     -51.15_dp,     -50.29_dp,     -50.31_dp,     -50.3_dp,     -50.8_dp,     &
    -50.77_dp,     -49.77_dp,     -50.47_dp,     -50.0_dp,     -50.9_dp,     -49.43_dp,     -50.25_dp,     -50.89_dp,     -50.76_dp,     -51.0_dp,     &
    -51.03_dp,     -50.98_dp,     -49.87_dp,     -51.44_dp,     -50.93_dp,     -50.99_dp,     -50.65_dp,     -50.03_dp,     -51.38_dp,     -50.48_dp,     &
    -50.38_dp,     -51.65_dp,     -51.3_dp,     -50.27_dp,     -51.61_dp,     -50.43_dp,     -50.47_dp,     -51.28_dp,     -51.11_dp,     -50.64_dp,     &
    -50.55_dp,     -50.16_dp,     -51.68_dp,     -50.68_dp,     -51.41_dp,     -51.45_dp,     -50.2_dp,     -51.44_dp,     -49.36_dp,     -50.68_dp,     &
    -50.65_dp,     -50.89_dp,     -50.56_dp,     -50.77_dp,     -51.29_dp,     -51.34_dp,     -50.33_dp,     -51.36_dp,     -49.56_dp,     -51.17_dp,     &
    -51.13_dp,     -50.3_dp,     -51.1_dp,     -50.86_dp,     -51.46_dp,     -51.04_dp,     -50.89_dp,     -49.96_dp,     -50.31_dp,     -50.82_dp,     &
    -50.38_dp,     -52.07_dp,     -51.68_dp,     -50.99_dp,     -51.17_dp,     -51.28_dp,     -52.22_dp,     -51.22_dp,     -50.53_dp,     -51.03_dp,     &
    -51.04_dp,     -50.28_dp,     -50.85_dp,     -49.87_dp,     -50.64_dp,     -50.23_dp,     -50.94_dp,     -50.26_dp,     -49.78_dp,     -50.29_dp,     &
    -50.47_dp,     -51.01_dp,     -50.5_dp,     -51.27_dp,     -49.59_dp,     -49.89_dp,     -49.91_dp,     -49.9_dp,     -50.48_dp,     -49.63_dp,     &
    -49.79_dp,     -49.59_dp,     -50.5_dp,     -51.17_dp,     -51.6_dp,     -50.09_dp,     -50.42_dp,     -50.77_dp,     -51.24_dp,     -50.58_dp,     &
    -51.07_dp,     -50.71_dp,     -50.57_dp,     -50.64_dp,     -50.18_dp,     -49.46_dp,     -51.71_dp,     -51.17_dp,     -51.64_dp,     -50.18_dp,     &
    -51.2_dp,     -51.98_dp,     -50.28_dp,     -51.09_dp,     -51.41_dp,     -50.03_dp,     -50.99_dp,     -51.75_dp,     -50.66_dp,     -50.85_dp,     &
    -50.88_dp,     -50.72_dp,     -50.53_dp,     -50.31_dp,     -51.5_dp,     -49.7_dp,     -50.46_dp,     -50.91_dp,     -50.56_dp,     -51.01_dp,     &
    -49.88_dp,     -49.96_dp,     -50.42_dp,     -50.02_dp,     -49.99_dp,     -51.44_dp,     -49.54_dp,     -49.63_dp,     -51.79_dp,     -49.7_dp,     &
    -50.65_dp,     -50.22_dp,     -51.32_dp,     -51.31_dp,     -49.93_dp,     -50.61_dp,     -50.85_dp,     -50.62_dp,     -50.85_dp,     -50.18_dp,     &
    -51.14_dp,     -51.07_dp,     -49.91_dp,     -49.97_dp,     -51.36_dp,     -49.67_dp,     -49.82_dp,     -50.82_dp,     -50.75_dp,     -50.31_dp,     &
    -50.03_dp,     -51.83_dp,     -49.96_dp,     -49.49_dp,     -50.47_dp,     -49.4_dp,     -50.86_dp,     -51.31_dp,     -50.46_dp,     -50.22_dp,     &
    -50.45_dp,     -50.83_dp,     -51.17_dp,     -50.23_dp,     -49.51_dp,     -50.86_dp,     -50.0_dp,     -51.07_dp,     -50.86_dp,     -51.23_dp,     &
    -50.44_dp,     -50.41_dp,     -49.99_dp,     -50.88_dp,     -50.64_dp,     -50.62_dp,     -50.72_dp,     -49.84_dp,     -50.91_dp,     -51.38_dp,     &
    -50.85_dp,     -50.82_dp,     -51.24_dp,     -50.73_dp,     -49.94_dp,     -50.38_dp,     -50.53_dp,     -50.49_dp,     -50.14_dp,     -49.86_dp,     &
    -50.93_dp,     -51.0_dp,     -50.54_dp,     -50.67_dp,     -50.24_dp,     -51.33_dp,     -50.47_dp,     -50.56_dp,     -50.12_dp,     -50.04_dp,     &
    -50.05_dp,     -50.35_dp,     -50.77_dp,     -50.56_dp,     -51.75_dp,     -50.74_dp,     -49.85_dp,     -50.46_dp,     -49.86_dp,     -50.61_dp,     &
    -50.42_dp,     -50.16_dp,     -52.12_dp,     -51.23_dp,     -51.41_dp,     -50.56_dp,     -51.56_dp,     -52.2_dp,     -51.06_dp,     -50.48_dp,     &
    -49.54_dp,     -50.65_dp,     -50.49_dp,     -51.11_dp,     -50.57_dp,     -50.54_dp,     -51.11_dp,     -50.07_dp,     -51.01_dp,     -49.35_dp,     &
    -51.17_dp,     -49.8_dp,     -49.76_dp,     -51.2_dp,     -50.14_dp,     -49.96_dp,     -50.15_dp,     -50.8_dp,     -51.11_dp,     -51.02_dp,     &
    -51.7_dp,     -50.88_dp,     -50.02_dp,     -50.55_dp,     -50.83_dp,     -51.12_dp,     -50.36_dp,     -50.05_dp,     -51.6_dp,     -49.8_dp,     &
    -51.25_dp,     -50.95_dp,     -50.89_dp,     -50.42_dp,     -51.01_dp,     -51.37_dp,     -50.35_dp,     -50.7_dp,     -51.54_dp,     -51.72_dp,     &
    -50.62_dp,     -51.53_dp,     -51.11_dp,     -50.79_dp,     -49.92_dp,     -50.48_dp,     -51.65_dp,     -51.44_dp,     -49.77_dp,     -50.84_dp,     &
    -51.54_dp,     -51.18_dp,     -52.15_dp,     -50.63_dp,     -51.57_dp,     -51.97_dp,     -50.01_dp,     -50.06_dp,     -50.51_dp,     -51.47_dp,     &
    -51.34_dp,     -51.16_dp,     -51.24_dp,     -50.88_dp,     -49.92_dp,     -50.44_dp,     -49.65_dp,     -50.97_dp,     -51.06_dp,     -51.25_dp,     &
    -50.85_dp,     -50.79_dp,     -51.62_dp,     -51.01_dp,     -50.87_dp,     -51.23_dp,     -51.34_dp,     -50.84_dp,     -51.63_dp,     -52.63_dp,     &
    -51.42_dp,     -51.16_dp,     -51.66_dp,     -51.42_dp,     -51.61_dp,     -50.76_dp,     -49.67_dp,     -50.94_dp,     -51.93_dp,     -51.13_dp,     &
    -51.73_dp,     -51.57_dp,     -50.46_dp,     -51.43_dp,     -51.61_dp,     -51.01_dp,     -50.31_dp,     -51.09_dp,     -49.76_dp,     -50.83_dp,     &
    -49.77_dp,     -51.34_dp,     -51.72_dp,     -50.87_dp,     -50.51_dp,     -49.92_dp,     -51.65_dp,     -51.16_dp,     -49.75_dp,     -49.92_dp,     &
    -51.95_dp,     -51.21_dp,     -51.4_dp,     -51.44_dp,     -50.57_dp,     -50.01_dp,     -50.76_dp,     -51.06_dp,     -51.12_dp,     -51.34_dp,     &
    -51.17_dp,     -51.96_dp,     -51.09_dp,     -52.79_dp,     -51.13_dp,     -50.76_dp,     -51.46_dp,     -51.77_dp,     -52.31_dp,     -52.24_dp,     &
    -51.22_dp,     -51.58_dp,     -50.83_dp,     -51.2_dp,     -51.44_dp,     -51.26_dp,     -51.26_dp,     -52.03_dp,     -51.66_dp,     -51.82_dp,     &
    -51.47_dp,     -51.0_dp,     -50.79_dp,     -50.83_dp,     -52.39_dp,     -50.82_dp,     -49.95_dp,     -50.88_dp,     -51.25_dp,     -50.23_dp,     &
    -51.2_dp,     -51.11_dp,     -50.81_dp,     -51.43_dp,     -51.41_dp,     -50.91_dp,     -51.3_dp,     -50.53_dp,     -51.79_dp,     -51.51_dp,     &
    -50.88_dp,     -50.5_dp,     -50.96_dp,     -50.69_dp,     -50.77_dp,     -51.0_dp,     -50.04_dp,     -50.28_dp,     -51.29_dp,     -50.62_dp,     &
    -50.76_dp,     -50.57_dp,     -50.97_dp,     -50.51_dp,     -50.62_dp,     -52.09_dp,     -51.56_dp,     -50.13_dp,     -51.1_dp,     -49.96_dp,     &
    -49.66_dp,     -50.36_dp,     -50.77_dp,     -51.36_dp,     -50.21_dp,     -50.27_dp,     -50.34_dp,     -50.51_dp,     -50.14_dp,     -50.2_dp,     &
    -49.85_dp,     -49.74_dp,     -50.41_dp,     -51.69_dp,     -50.29_dp,     -50.95_dp,     -50.05_dp,     -50.49_dp,     -49.3_dp,     -50.08_dp,     &
    -50.91_dp,     -50.13_dp,     -50.41_dp,     -49.08_dp,     -50.07_dp,     -49.95_dp,     -49.63_dp,     -49.99_dp,     -48.91_dp,     -49.53_dp,     &
    -51.94_dp,     -50.06_dp,     -50.17_dp,     -51.19_dp,     -49.62_dp,     -50.55_dp,     -49.93_dp,     -50.48_dp,     -50.81_dp,     -49.63_dp,     &
    -50.36_dp,     -50.4_dp,     -49.94_dp,     -49.72_dp,     -50.0_dp,     -49.65_dp,     -50.81_dp,     -49.84_dp,     -48.89_dp,     -49.43_dp,     &
    -49.71_dp,     -50.16_dp,     -49.77_dp,     -50.14_dp,     -50.68_dp,     -51.06_dp,     -49.5_dp,     -49.99_dp,     -50.01_dp,     -49.82_dp,     &
    -51.15_dp,     -49.68_dp,     -50.89_dp,     -50.25_dp,     -50.32_dp,     -50.98_dp,     -49.65_dp,     -49.42_dp,     -49.01_dp,     -49.58_dp,     &
    -49.8_dp,     -48.06_dp,     -50.59_dp,     -49.88_dp,     -50.16_dp,     -50.42_dp,     -50.34_dp,     -50.19_dp,     -49.9_dp,     -48.94_dp,     &
    -50.32_dp,     -50.45_dp,     -50.33_dp,     -50.65_dp,     -50.15_dp,     -49.89_dp,     -49.71_dp,     -49.36_dp,     -49.9_dp,     -50.77_dp,     &
    -50.29_dp,     -49.79_dp,     -49.53_dp,     -50.36_dp,     -50.96_dp,     -50.11_dp,     -49.61_dp,     -49.91_dp,     -49.38_dp,     -49.69_dp,     &
    -50.8_dp,     -50.46_dp,     -49.86_dp,     -50.53_dp,     -50.75_dp,     -50.34_dp,     -49.28_dp,     -49.25_dp,     -49.43_dp,     -49.85_dp,     &
    -51.21_dp,     -50.07_dp,     -49.85_dp,     -50.23_dp,     -48.57_dp,     -48.69_dp,     -50.9_dp,     -50.26_dp,     -50.62_dp,     -50.11_dp,     &
    -49.62_dp,     -49.18_dp,     -49.62_dp,     -48.82_dp,     -49.27_dp,     -49.92_dp,     -50.34_dp,     -50.13_dp,     -50.27_dp,     -50.32_dp,     &
    -49.47_dp,     -50.14_dp,     -49.46_dp,     -50.23_dp,     -49.49_dp,     -50.46_dp,     -49.86_dp,     -50.17_dp,     -49.29_dp,     -50.36_dp,     &
    -51.37_dp,     -50.83_dp,     -51.01_dp,     -50.74_dp,     -50.99_dp,     -50.77_dp,     -51.49_dp,     -51.51_dp,     -50.66_dp,     -50.82_dp,     &
    -50.66_dp,     -51.67_dp,     -51.37_dp,     -51.59_dp,     -52.94_dp,     -51.84_dp,     -51.67_dp,     -52.4_dp,     -51.36_dp,     -51.08_dp,     &
    -52.03_dp,     -50.97_dp,     -52.54_dp,     -52.82_dp,     -52.11_dp,     -53.19_dp,     -51.78_dp,     -51.35_dp,     -52.36_dp,     -52.43_dp,     &
    -52.29_dp,     -53.51_dp,     -53.11_dp,     -53.21_dp,     -52.83_dp,     -53.3_dp,     -52.24_dp,     -52.63_dp,     -52.42_dp,     -53.3_dp,     &
    -53.63_dp,     -52.52_dp,     -52.98_dp,     -52.04_dp,     -52.57_dp,     -52.95_dp,     -53.03_dp,     -52.72_dp,     -52.17_dp,     -52.54_dp,     &
    -52.68_dp,     -53.11_dp,     -52.87_dp,     -52.52_dp,     -51.9_dp,     -52.28_dp,     -52.78_dp,     -52.14_dp,     -51.91_dp,     -52.43_dp,     &
    -52.62_dp,     -53.43_dp,     -53.17_dp,     -52.63_dp,     -53.14_dp,     -52.36_dp,     -51.38_dp,     -52.6_dp,     -52.25_dp,     -52.97_dp,     &
    -53.06_dp,     -53.22_dp,     -53.15_dp,     -52.16_dp,     -52.91_dp,     -52.54_dp,     -51.48_dp,     -53.0_dp,     -51.54_dp,     -52.68_dp,     &
    -52.83_dp,     -52.42_dp,     -51.79_dp,     -51.24_dp,     -51.01_dp,     -51.57_dp,     -52.3_dp,     -52.19_dp,     -52.44_dp,     -51.54_dp,     &
    -51.15_dp,     -51.45_dp,     -51.69_dp,     -52.33_dp,     -49.94_dp,     -50.48_dp,     -50.66_dp,     -51.08_dp,     -50.87_dp,     -51.27_dp,     &
    -50.93_dp,     -53.05_dp,     -51.82_dp,     -51.02_dp,     -52.08_dp,     -52.15_dp,     -51.68_dp,     -52.49_dp,     -51.93_dp,     -52.97_dp,     &
    -52.74_dp,     -52.14_dp,     -52.66_dp,     -52.43_dp,     -51.76_dp,     -51.95_dp,     -52.19_dp,     -51.75_dp,     -51.98_dp,     -52.19_dp,     &
    -53.21_dp,     -50.88_dp,     -52.41_dp,     -52.82_dp,     -52.81_dp,     -51.84_dp,     -52.16_dp,     -52.35_dp,     -52.12_dp,     -52.59_dp,     &
    -52.8_dp,     -53.34_dp,     -53.01_dp,     -53.14_dp,     -53.06_dp,     -53.31_dp,     -53.27_dp,     -52.65_dp,     -52.87_dp,     -53.02_dp,     &
    -53.44_dp,     -52.93_dp,     -53.75_dp,     -52.57_dp,     -52.87_dp,     -53.0_dp,     -53.15_dp,     -52.93_dp,     -53.16_dp,     -53.38_dp,     &
    -53.59_dp,     -52.94_dp,     -53.96_dp,     -53.38_dp,     -53.33_dp,     -53.83_dp,     -54.25_dp,     -53.93_dp,     -53.31_dp,     -53.44_dp,     &
    -53.47_dp,     -53.38_dp,     -54.49_dp,     -53.61_dp,     -53.37_dp,     -53.79_dp,     -53.93_dp,     -53.73_dp,     -54.21_dp,     -55.16_dp,     &
    -54.19_dp,     -54.53_dp,     -54.19_dp,     -53.92_dp,     -55.5_dp,     -54.47_dp,     -55.11_dp,     -55.04_dp,     -55.19_dp,     -55.66_dp,     &
    -54.66_dp,     -54.81_dp,     -55.27_dp,     -54.27_dp,     -55.3_dp,     -54.4_dp,     -56.04_dp,     -55.21_dp,     -55.19_dp,     -55.76_dp,     &
    -55.29_dp,     -55.73_dp,     -56.03_dp,     -56.08_dp,     -56.07_dp,     -56.22_dp,     -54.84_dp,     -55.54_dp,     -56.24_dp,     -55.49_dp,     &
    -55.5_dp,     -55.91_dp,     -56.39_dp,     -56.23_dp,     -55.67_dp,     -55.86_dp,     -55.55_dp,     -55.83_dp,     -55.96_dp,     -55.91_dp,     &
    -55.97_dp,     -56.15_dp,     -56.03_dp,     -56.69_dp,     -55.4_dp,     -55.75_dp,     -55.04_dp,     -56.25_dp,     -55.63_dp,     -55.98_dp,     &
    -55.87_dp,     -55.83_dp,     -56.15_dp,     -55.51_dp,     -56.57_dp,     -55.9_dp,     -55.16_dp,     -56.71_dp,     -56.27_dp,     -55.09_dp,     &
    -56.07_dp,     -55.87_dp,     -54.67_dp,     -56.07_dp,     -55.52_dp,     -56.27_dp,     -55.69_dp,     -56.23_dp,     -56.18_dp,     -56.55_dp,     &
    -55.96_dp,     -56.46_dp,     -56.38_dp,     -55.6_dp,     -55.6_dp,     -55.24_dp,     -55.95_dp,     -56.74_dp,     -56.37_dp,     -55.53_dp,     &
    -55.51_dp,     -55.99_dp,     -55.59_dp,     -56.36_dp,     -55.88_dp,     -55.88_dp,     -56.32_dp,     -56.31_dp,     -56.67_dp,     -55.73_dp,     &
    -55.6_dp,     -55.84_dp,     -56.0_dp,     -56.24_dp /)

    ! Initialize found flag
    !found = .FALSE.
    
    ! Find the two indices such that age_data(idx_low) <= time <= age_data(idx_high)
    idx_low = -1
    idx_high = -1

    ! Edge cases
    IF (time <= age_data(1)) THEN
        ! Time is before the first data point
        idx_low = 1
        idx_high = 1
    ELSE IF (time >= age_data(num_points)) THEN
        ! Time is after the last data point
        idx_low = num_points
        idx_high = num_points
    ! Main data search
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

    ! Interpolate d18O values to match to model runtime
    d18O_current = delta18O_data(idx_low) + fraction * (delta18O_data(idx_high) - delta18O_data(idx_low))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_d18O_values

  SUBROUTINE apply_anomaly_fields(mesh, ocean, matrix, time)
    ! Apply climatological annual mean anomaly fields progressively from PI to LGM over 21,000 years
  
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'apply_anomaly_fields'
    REAL(dp)                                              :: w_global
    INTEGER                                               :: i, j, vi, k
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: tas, sos
    REAL(dp)                                              :: wt0, wt1
  
    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Calculate weights for linear interpolation
    wt0 = (time - ocean%matrix%t1) / (ocean%matrix%t0 - ocean%matrix%t1)

    print *, "Interpolation weight calculated = ", wt0

    ! Apply linear interpolation
    DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
            ocean%T(vi,k) = wt0 * ocean%matrix%timeframe0%T(vi,k) + wt1 * ocean%matrix%timeframe1%T(vi,k)
            ocean%S(vi,k) = wt0 * ocean%matrix%timeframe0%S(vi,k) + wt1 * ocean%matrix%timeframe1%S(vi,k)
        END DO
    END DO
  
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE apply_anomaly_fields  

  SUBROUTINE interpolation_with_insolation_and_GHG_radiative(mesh, ocean, matrix, time, grid_smooth)
    ! Linear interpolation forcing based on insolation values

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_insolation_and_GHG_radiative'
    INTEGER                                               :: i, j, vi, k
    ! Insolation
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: ins_current, w_ins
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: w_ins_smooth, w_ice
    CHARACTER(LEN=256)                                    :: choice_insolation_forcing, filename_insolation
    REAL(dp)                                              :: time_to_read, r, w_ins_av, sum_local
    LOGICAL                                               :: apply_forcing
    INTEGER                                               :: count_pts
    ! GHG
    REAL(dp)                                              :: w_GHG
    REAL(dp)                                              :: CO2_current, CH4_current, N2O_current
    REAL(dp)                                              :: CO2_PI, CH4_PI, N2O_PI
    REAL(dp)                                              :: CO2_LGM, CH4_LGM, N2O_LGM
    REAL(dp)                                              :: DeltaF_CO2, DeltaF_CH4, DeltaF_N2O
    REAL(dp)                                              :: DeltaF_CO2_PI, DeltaF_CH4_PI, DeltaF_N2O_PI
    REAL(dp)                                              :: DeltaF_CO2_LGM, DeltaF_CH4_LGM, DeltaF_N2O_LGM
    REAL(dp)                                              :: DeltaF_total, DeltaF_PI, DeltaF_LGM
    CHARACTER(LEN=256)                                    :: CO2_relationship, GHG_inclusion
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: w_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ======================================= w_ins =======================================
    ! Retrieve configuration settings
    choice_insolation_forcing = TRIM(C%choice_insolation_forcing)
    filename_insolation       = TRIM(C%filename_insolation)

    ! Initialize apply_forcing flag based on forcing choice
    IF (choice_insolation_forcing == 'none') THEN
      apply_forcing = .FALSE.
    ELSE
      apply_forcing = .TRUE.
    END IF

    IF (apply_forcing) THEN

      ! Allocate arrays
      IF (.NOT. ALLOCATED(ins_current)) THEN
        ALLOCATE(ins_current(mesh%vi1:mesh%vi2))
      END IF
      IF (.NOT. ALLOCATED(w_ins)) THEN
        ALLOCATE(w_ins(mesh%vi1:mesh%vi2))
      END IF
  
      ! Determine the time to read based on the forcing choice
      SELECT CASE (choice_insolation_forcing)
        CASE ('static')
          time_to_read = C%static_insolation_time
        CASE ('realistic')
          time_to_read = time
        CASE DEFAULT
          CALL crash('Unknown choice_insolation_forcing: "' // choice_insolation_forcing // '"')
      END SELECT
  
      ! Read insolation field at the determined time
      CALL read_field_from_file_2D(TRIM(filename_insolation), 'Q_TOA', mesh, ins_current, time_to_read = time_to_read)

      ! Compute w_ins for each mesh point
      DO vi = mesh%vi1, mesh%vi2
        ! w_ins = (I_current - I_LGM) / (I_PI - I_LGM)
        w_ins(vi) = (ins_current(vi) - ocean%matrix%Q_TOA_LGM(vi)) / (ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi))

        ! Handle division by zero or very small denominators
        IF (ABS(ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi)) < 1E-5_dp) THEN
          w_ins(vi) = 0.0_dp
        END IF
      END DO

      ! Smooth insolation field
      IF (.NOT. ALLOCATED(w_ins_smooth)) THEN
        ALLOCATE(w_ins_smooth(mesh%vi1:mesh%vi2))
      END IF

      w_ins_smooth = w_ins

      r = 200000._dp
      CALL smooth_Gaussian_2D(mesh, grid_smooth, w_ins_smooth, r)

      ! Calculate average
      sum_local = 0.0_dp
      count_pts = mesh%vi2 - mesh%vi1 + 1

      DO vi = mesh%vi1, mesh%vi2
        sum_local = sum_local + w_ins(vi)
      END DO

      w_ins_av = sum_local / REAL(count_pts, dp)

      ! Combine factors
      IF (.NOT. ALLOCATED(w_ice)) THEN
        ALLOCATE(w_ice(mesh%vi1:mesh%vi2))
      END IF

      DO vi = mesh%vi1, mesh%vi2
        w_ice(vi) = (1.0_dp * w_ins(vi) + &
                     3.0_dp * w_ins_smooth(vi) + &
                     3.0_dp * w_ins_av) / 7.0_dp
      END DO
  
      ! ======================================= w_GHG =======================================

      ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
      CALL get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)
      CALL get_GHG_concentrations(0.0_dp, CO2_PI, CH4_PI, N2O_PI)
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

      ! ======================================= w_tot =======================================

      IF (.NOT. ALLOCATED(w_tot)) THEN
        ALLOCATE(w_tot(mesh%vi1:mesh%vi2))
      END IF
      
      DO vi = mesh%vi1, mesh%vi2
        w_tot(vi) = (w_ice(vi) + w_GHG) / 2
      END DO

      ! Clamp weights between cutoff values if enabled
      IF (C%clamp_weights) THEN
        DO vi = mesh%vi1, mesh%vi2
          w_tot(vi) = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_tot(vi)))
        END DO
      END IF

      ! Apply interpolation using w_tot
      DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
          ocean%T(vi,k) = w_tot(vi) * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_tot(vi)) * ocean%matrix%timeframe1%T(vi,k)
          ocean%S(vi,k) = w_tot(vi) * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_tot(vi)) * ocean%matrix%timeframe1%S(vi,k)
        END DO
      END DO
  
      ! Deallocate arrays
      DEALLOCATE(ins_current, w_ins, w_tot, w_ins_smooth, w_ice)
  
    ELSE
      ! If forcing is 'none', do not modify the ocean state
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_insolation_and_GHG_radiative

  SUBROUTINE interpolation_with_insolation_and_GHG(mesh, ocean, matrix, time, grid_smooth)
    ! Linear interpolation forcing based on insolation values

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_insolation_and_GHG'
    INTEGER                                               :: i, j, vi, k
    ! Insolation
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: ins_current, w_ins
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: w_ins_smooth, w_ice
    CHARACTER(LEN=256)                                    :: choice_insolation_forcing, filename_insolation
    REAL(dp)                                              :: time_to_read, r, w_ins_av, sum_local
    LOGICAL                                               :: apply_forcing
    INTEGER                                               :: count_pts
    ! GHG
    REAL(dp)                                              :: w_GHG, w_GHG_CO2, w_GHG_CH4, w_GHG_N2O
    REAL(dp)                                              :: CO2_current, CH4_current, N2O_current
    REAL(dp)                                              :: CO2_PI, CH4_PI, N2O_PI
    REAL(dp)                                              :: CO2_LGM, CH4_LGM, N2O_LGM
    CHARACTER(LEN=256)                                    :: GHG_inclusion
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: w_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ======================================= w_ins =======================================
    ! Retrieve configuration settings
    choice_insolation_forcing = TRIM(C%choice_insolation_forcing)
    filename_insolation       = TRIM(C%filename_insolation)

    ! Initialize apply_forcing flag based on forcing choice
    IF (choice_insolation_forcing == 'none') THEN
      apply_forcing = .FALSE.
    ELSE
      apply_forcing = .TRUE.
    END IF

    IF (apply_forcing) THEN

      ! Allocate arrays
      IF (.NOT. ALLOCATED(ins_current)) THEN
        ALLOCATE(ins_current(mesh%vi1:mesh%vi2))
      END IF
      IF (.NOT. ALLOCATED(w_ins)) THEN
        ALLOCATE(w_ins(mesh%vi1:mesh%vi2))
      END IF
  
      ! Determine the time to read based on the forcing choice
      SELECT CASE (choice_insolation_forcing)
        CASE ('static')
          time_to_read = C%static_insolation_time
        CASE ('realistic')
          time_to_read = time
        CASE DEFAULT
          CALL crash('Unknown choice_insolation_forcing: "' // choice_insolation_forcing // '"')
      END SELECT
  
      ! Read insolation field at the determined time
      CALL read_field_from_file_2D(TRIM(filename_insolation), 'Q_TOA', mesh, ins_current, time_to_read = time_to_read)

      ! Compute w_ins for each mesh point
      DO vi = mesh%vi1, mesh%vi2
        ! w_ins = (I_current - I_LGM) / (I_PI - I_LGM)
        w_ins(vi) = (ins_current(vi) - ocean%matrix%Q_TOA_LGM(vi)) / (ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi))

        ! Handle division by zero or very small denominators
        IF (ABS(ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi)) < 1E-5_dp) THEN
          w_ins(vi) = 0.0_dp
        END IF
      END DO

      ! Smooth insolation field
      IF (.NOT. ALLOCATED(w_ins_smooth)) THEN
        ALLOCATE(w_ins_smooth(mesh%vi1:mesh%vi2))
      END IF

      w_ins_smooth = w_ins

      r = 200000._dp
      CALL smooth_Gaussian_2D(mesh, grid_smooth, w_ins_smooth, r)

      ! Calculate average
      sum_local = 0.0_dp
      count_pts = mesh%vi2 - mesh%vi1 + 1

      DO vi = mesh%vi1, mesh%vi2
        sum_local = sum_local + w_ins(vi)
      END DO

      w_ins_av = sum_local / REAL(count_pts, dp)

      ! Combine factors
      IF (.NOT. ALLOCATED(w_ice)) THEN
        ALLOCATE(w_ice(mesh%vi1:mesh%vi2))
      END IF

      DO vi = mesh%vi1, mesh%vi2
        w_ice(vi) = (1.0_dp * w_ins(vi) + &
                     3.0_dp * w_ins_smooth(vi) + &
                     3.0_dp * w_ins_av) / 7.0_dp
      END DO
  
      ! ======================================= w_GHG =======================================

      ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
      CALL get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)  ! Concentration during runtime
      CALL get_GHG_concentrations(0.0_dp, CO2_PI, CH4_PI, N2O_PI)             ! PI
      CALL get_GHG_concentrations(21000.0_dp, CO2_LGM, CH4_LGM, N2O_LGM)        ! LGM
  
      ! Get config settings
      GHG_inclusion = TRIM(C%choice_ghg_inclusion)
  
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

      ! ======================================= w_tot =======================================
  
      IF (.NOT. ALLOCATED(w_tot)) THEN
        ALLOCATE(w_tot(mesh%vi1:mesh%vi2))
      END IF

      DO vi = mesh%vi1, mesh%vi2
        w_tot(vi) = (w_ice(vi) + w_GHG) / 2
      END DO

      ! Clamp weights between cutoff values if enabled
      IF (C%clamp_weights) THEN
        DO vi = mesh%vi1, mesh%vi2
          w_tot(vi) = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_tot(vi)))
        END DO
      END IF

      ! Apply interpolation using w_tot
      DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
          ocean%T(vi,k) = w_tot(vi) * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_tot(vi)) * ocean%matrix%timeframe1%T(vi,k)
          ocean%S(vi,k) = w_tot(vi) * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_tot(vi)) * ocean%matrix%timeframe1%S(vi,k)
        END DO
      END DO
  
      ! Deallocate arrays
      DEALLOCATE(ins_current, w_ins, w_tot, w_ins_smooth, w_ice)
  
    ELSE
      ! If forcing is 'none', do not modify the ocean state
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_insolation_and_GHG

  SUBROUTINE interpolation_with_insolation(mesh, ocean, matrix, time, grid_smooth)
    ! Linear interpolation forcing based on insolation values

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_insolation'
    INTEGER                                               :: i, j, vi, k, count_pts
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: ins_current, w_ins
    REAL(dp), DIMENSION(:), ALLOCATABLE                   :: w_ins_smooth, w_ice
    CHARACTER(LEN=256)                                    :: choice_insolation_forcing, filename_insolation
    REAL(dp)                                              :: time_to_read, r, w_ins_av, sum_local
    LOGICAL                                               :: apply_forcing    

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Retrieve configuration settings
    choice_insolation_forcing = TRIM(C%choice_insolation_forcing)
    filename_insolation       = TRIM(C%filename_insolation)

    ! Initialize apply_forcing flag based on forcing choice
    IF (choice_insolation_forcing == 'none') THEN
      apply_forcing = .FALSE.
    ELSE
      apply_forcing = .TRUE.
    END IF

    IF (apply_forcing) THEN

      ! Allocate arrays
      IF (.NOT. ALLOCATED(ins_current)) THEN
        ALLOCATE(ins_current(mesh%vi1:mesh%vi2))
      END IF
      IF (.NOT. ALLOCATED(w_ins)) THEN
        ALLOCATE(w_ins(mesh%vi1:mesh%vi2))
      END IF
  
      ! Determine the time to read based on the forcing choice
      SELECT CASE (choice_insolation_forcing)
        CASE ('static')
          time_to_read = C%static_insolation_time
        CASE ('realistic')
          time_to_read = time
        CASE DEFAULT
          CALL crash('Unknown choice_insolation_forcing: "' // choice_insolation_forcing // '"')
      END SELECT
  
      ! Read insolation field at the determined time
      CALL read_field_from_file_2D(TRIM(filename_insolation), 'Q_TOA', mesh, ins_current, time_to_read = time_to_read)

      ! Compute w_ins for each mesh point
      DO vi = mesh%vi1, mesh%vi2
        ! w_ins = (I_current - I_LGM) / (I_PI - I_LGM)
        w_ins(vi) = (ins_current(vi) - ocean%matrix%Q_TOA_LGM(vi)) / (ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi))

        ! Handle division by zero or very small denominators
        IF (ABS(ocean%matrix%Q_TOA_PI(vi) - ocean%matrix%Q_TOA_LGM(vi)) < 1E-5_dp) THEN
          w_ins(vi) = 0.0_dp
        END IF
      END DO

      ! Smooth insolation field
      IF (.NOT. ALLOCATED(w_ins_smooth)) THEN
        ALLOCATE(w_ins_smooth(mesh%vi1:mesh%vi2))
      END IF

      w_ins_smooth = w_ins

      r = 200000._dp
      CALL smooth_Gaussian_2D(mesh, grid_smooth, w_ins_smooth, r)

      ! Calculate average
      sum_local = 0.0_dp
      count_pts = mesh%vi2 - mesh%vi1 + 1

      DO vi = mesh%vi1, mesh%vi2
        sum_local = sum_local + w_ins(vi)
      END DO

      w_ins_av = sum_local / REAL(count_pts, dp)

      ! Combine factors
      IF (.NOT. ALLOCATED(w_ice)) THEN
        ALLOCATE(w_ice(mesh%vi1:mesh%vi2))
      END IF

      DO vi = mesh%vi1, mesh%vi2
        w_ice(vi) = (1.0_dp * w_ins(vi) + &
                     3.0_dp * w_ins_smooth(vi) + &
                     3.0_dp * w_ins_av) / 7.0_dp
      END DO

      ! Clamp weights between cutoff values if enabled
      IF (C%clamp_weights) THEN
        DO vi = mesh%vi1, mesh%vi2
          w_ice(vi) = MAX(C%clamp_cutoff_low, MIN(C%clamp_cutoff_high, w_ice(vi)))
        END DO
      END IF

      ! Apply interpolation using w_ins
      DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
          ocean%T(vi,k) = w_ice(vi) * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_ice(vi)) * ocean%matrix%timeframe1%T(vi,k)
          ocean%S(vi,k) = w_ice(vi) * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_ice(vi)) * ocean%matrix%timeframe1%S(vi,k)
        END DO
      END DO
  
      ! Deallocate arrays
      DEALLOCATE(ins_current, w_ins, w_ins_smooth, w_ice)
  
    ELSE
      ! If forcing is 'none', do not modify the ocean state
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_insolation

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
    INTEGER                                               :: i, j, vi, k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check for division by 0 error
    IF (ABS(ocean%matrix%t1 - ocean%matrix%t0) < 1e-5_dp) THEN
      CALL crash('t0 and t1 are too close or identical, interpolation cannot be performed.')
    END IF

    ! Calculate weights for linear interpolation
    wt0 = (time - ocean%matrix%t1) / (ocean%matrix%t0 - ocean%matrix%t1)

    print *, "Interpolation weight calculated =", wt0

    ! Apply linear interpolation
    DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
            ocean%T(vi,k) = wt0 * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - wt0) * ocean%matrix%timeframe1%T(vi,k)
            ocean%S(vi,k) = wt0 * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - wt0) * ocean%matrix%timeframe1%S(vi,k)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE linear_time_interpolation

  SUBROUTINE interpolation_with_GHG_basic(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on GHG

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_GHG_basic'
    REAL(dp)                                              :: w_GHG, w_GHG_CO2, w_GHG_CH4, w_GHG_N2O
    INTEGER                                               :: i, j, vi, k
    REAL(dp)                                              :: CO2_current, CH4_current, N2O_current
    REAL(dp)                                              :: CO2_PI, CH4_PI, N2O_PI
    REAL(dp)                                              :: CO2_LGM, CH4_LGM, N2O_LGM
    CHARACTER(LEN=256)                                    :: GHG_inclusion

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get current, LGM, and PI GHG concentrations (time must match entry in age_data)
    CALL get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)  ! Concentration during runtime
    CALL get_GHG_concentrations(0.0_dp, CO2_PI, CH4_PI, N2O_PI)             ! PI
    CALL get_GHG_concentrations(21000.0_dp, CO2_LGM, CH4_LGM, N2O_LGM)        ! LGM

    ! Get config settings
    GHG_inclusion = TRIM(C%choice_ghg_inclusion)

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

    print *, "Interpolation weight calculated = ", w_GHG

    ! Apply interpolation using w_GHG
    DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
            ocean%T(vi,k) = w_GHG * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_GHG) * ocean%matrix%timeframe1%T(vi,k)
            ocean%S(vi,k) = w_GHG * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_GHG) * ocean%matrix%timeframe1%S(vi,k)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_GHG_basic

  SUBROUTINE interpolation_with_GHG_radiative(mesh, ocean, matrix, time)
    ! Linear interpolation forcing based on GHG

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'interpolation_with_GHG_radiative'
    REAL(dp)                                              :: w_GHG
    INTEGER                                               :: i, j, vi, k
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
    CALL get_GHG_concentrations(0.0_dp, CO2_PI, CH4_PI, N2O_PI)
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
    DO vi = mesh%vi1, mesh%vi2
        DO k = 1, C%nz_ocean
            ocean%T(vi,k) = w_GHG * ocean%matrix%timeframe0%T(vi,k) + (1.0_dp - w_GHG) * ocean%matrix%timeframe1%T(vi,k)
            ocean%S(vi,k) = w_GHG * ocean%matrix%timeframe0%S(vi,k) + (1.0_dp - w_GHG) * ocean%matrix%timeframe1%S(vi,k)
        END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE interpolation_with_GHG_radiative

  SUBROUTINE get_GHG_concentrations(time, CO2_current, CH4_current, N2O_current)
    ! Get GHG concentrations at the given time, linearly interpolated concentrations to match age data
    ! CH4 & N2O: EPICA Dome C – Nitrous Oxide and Methane Data (Spahni, R.)
    ! CO2:       EPICA Dome C - 800KYr CO2 Data (Luthi, D.)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), INTENT(IN)                                  :: time
    REAL(dp), INTENT(OUT)                                 :: CO2_current
    REAL(dp), INTENT(OUT)                                 :: CH4_current
    REAL(dp), INTENT(OUT)                                 :: N2O_current

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'get_GHG_concentrations'
    INTEGER, PARAMETER                                    :: num_points = 373
    REAL(dp), DIMENSION(num_points)                       :: age_data, CO2_data, CH4_data, N2O_data
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

    ! Edge cases
    IF (time <= age_data(1)) THEN
        ! Time is before the first data point
        idx_low = 1
        idx_high = 1
    ELSE IF (time >= age_data(num_points)) THEN
        ! Time is after the last data point
        idx_low = num_points
        idx_high = num_points
    ! Main data search
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

    ! Interpolate GHG concentrations to match to model runtime
    CO2_current = CO2_data(idx_low) + fraction * (CO2_data(idx_high) - CO2_data(idx_low))
    CH4_current = CH4_data(idx_low) + fraction * (CH4_data(idx_high) - CH4_data(idx_low))
    N2O_current = N2O_data(idx_low) + fraction * (N2O_data(idx_high) - N2O_data(idx_low))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_GHG_concentrations

  SUBROUTINE compute_CO2_forcing(CO2_current, CO2_ref, CO2_relationship, DeltaF_CO2)
    ! Compute radiative forcing due to CO2 using specified relationship
    ! IPCC TAR-06 (2018)

    IMPLICIT NONE

    ! In/output variables:
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

    CASE ('relationship1') ! IPCC (1990)
      ! ∆F= α ln(C/C0)
      alpha = 5.35_dp
      DeltaF_CO2 = alpha * LOG(CO2_current / CO2_ref)
      
    CASE ('relationship2') ! Shi (1992)
      ! ∆F= α ln(C/C0) + β(√C − √C0)
      alpha = 4.841_dp
      beta = 0.0906_dp
      DeltaF_CO2 = alpha * LOG(CO2_current / CO2_ref) + beta * (SQRT(CO2_current) - SQRT(CO2_ref))
    
    CASE ('relationship3') ! WMO (1999)
      ! ∆F= α(g(C)–g(C0)), where g(C) = ln(1 + 1.2C + 0.005C² + 1.4 × 10⁻⁶C³)
      alpha = 3.35_dp
      g_current = LOG(1.0_dp + 1.2_dp*CO2_current + 0.005_dp*CO2_current**2 + 1.4e-6_dp*CO2_current**3)
      g_ref = LOG(1.0_dp + 1.2_dp*CO2_ref + 0.005_dp*CO2_ref**2 + 1.4e-6_dp*CO2_ref**3)
      DeltaF_CO2 = alpha * (g_current - g_ref) 
    
    CASE DEFAULT
      CALL crash('Unknown CO2 relationship: ' // CO2_relationship)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_CO2_forcing

  FUNCTION f_overlap(M, N) RESULT(f)
    ! Function radiative forcing CH4 and N2O
    ! f(M,N) = 0.47 * ln[1 + 2.01e-5 * (M*N)^0.75 + 5.31e-15 * M * (M*N)^1.52]
    ! IPCC TAR-06 (2018)

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
    ! IPCC TAR-06 (2018)
    
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
    ! IPCC TAR-06 (2018)
    
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

  SUBROUTINE run_ocean_model_matrix( mesh, ice, ocean, time, region_name, grid_smooth)
    ! Calculate the ocean
    ! Use an interpolating matrix ocean scheme
  
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth
  
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_matrix'
    TYPE(type_ocean_matrix_interpolation)                 :: matrix
    INTEGER                                               :: i, j, vi, k
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2)                :: mask_ocean
    REAL(dp)                                              :: max_ocean_size, sigma
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: T_field, S_field
    REAL(dp)                                              :: sea_level_current
    TYPE(type_reference_geometry)                         :: refgeo_PD

    ! Add routine to path
    CALL init_routine( routine_name) 

    ! Perform time interpolation
    SELECT CASE (TRIM(C%choice_ocean_model_matrix))
    CASE('linear_time')
      CALL linear_time_interpolation(mesh, ocean, matrix, time)
    CASE('GHG_radiative')
      CALL interpolation_with_GHG_radiative(mesh, ocean, matrix, time)
    CASE('GHG')
      CALL interpolation_with_GHG_basic(mesh, ocean, matrix, time)
    CASE('insolation')
      CALL interpolation_with_insolation(mesh, ocean, matrix, time, grid_smooth)
    CASE('d18O')
      CALL interpolation_with_d18O(mesh, ocean, matrix, time)
    CASE('anomaly_field')
      CALL apply_anomaly_fields(mesh, ocean, matrix, time)
    CASE('insolation+GHG')
      CALL interpolation_with_insolation_and_GHG(mesh, ocean, matrix, time, grid_smooth)
    CASE('insolation+GHG_radiative')
      CALL interpolation_with_insolation_and_GHG_radiative(mesh, ocean, matrix, time, grid_smooth)
    CASE DEFAULT
      CALL crash('Unknown choice_ocean_model_matrix' // TRIM(C%choice_ocean_model_matrix))
    END SELECT

    ! == Ocean Extrapolation after interpolation ==
    IF(C%enable_jourdain) THEN

      ! Initialize mask
      mask_ocean = 0
      
      ! Set basic mask (3 = horizontally extrapolated points, 2 = seed points, 1 = extrapolation points, 0 = ignore)
      DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_icefree_ocean(vi)) THEN
              mask_ocean(vi) = 2      ! Open ocean: use as seed
          ELSEIF (ice%mask_floating_ice(vi)) THEN
              mask_ocean(vi) = 1      ! Ice cavity: needs extrapolation
          END IF
      END DO
      
      ! Calculate extrapolation radius
      max_ocean_size = MINVAL(mesh%R)
      DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_floating_ice(vi)) THEN
              max_ocean_size = MAX(max_ocean_size, mesh%R(vi))
          END IF
      END DO
      sigma = max_ocean_size / 3._dp
      
      ! Print output
      IF (par%master) THEN
          WRITE(*,*) 'Ocean extrapolation runtime:'
          WRITE(*,*) '  - Number of seed points:', COUNT(mask_ocean == 2)
          WRITE(*,*) '  - Number of cavity points:', COUNT(mask_ocean == 1)
          WRITE(*,*) '  - Sigma:', sigma
      END IF
      
      ! Step 1: horizontal extrapolation into shelf cavities
      DO k = 1, C%nz_ocean
          ! Store current layer values
          T_field = ocean%T(:,k)
          S_field = ocean%S(:,k)
          
          ! Perform extrapolation
          CALL extrapolate_Gaussian(mesh, mask_ocean, T_field, sigma)
          CALL extrapolate_Gaussian(mesh, mask_ocean, S_field, sigma)
          
          ! Update only cavity points
          DO vi = mesh%vi1, mesh%vi2
              IF (mask_ocean(vi) == 1) THEN
                  ocean%T(vi,k) = T_field(vi)
                  ocean%S(vi,k) = S_field(vi)
                  ! After extrapolation mark as seed for vertical extrapolation
                  mask_ocean(vi) = 3  
              END IF
          END DO
      END DO

      ! Step 2: vertical extrapolation into sill-blocked shelf cavities
      DO vi = mesh%vi1, mesh%vi2
          IF (mask_ocean(vi) == 3) THEN
              DO k = 2, C%nz_ocean
                  ocean%T(vi,k) = ocean%T(vi,k-1)
                  ocean%S(vi,k) = ocean%S(vi,k-1)
              END DO
          END IF
      END DO

    END IF

    ! === Sea level ===
    SELECT CASE (TRIM(C%choice_sealevel_model))
      CASE ('fixed')
        ! Fixed sea level
        sea_level_current = C%fixed_sealevel

      CASE ('prescribed')
        ! Sea-level prescribed from external record file
        CALL get_sea_level_values(time, sea_level_current)

      CASE ('eustatic')
        ! Eustatic sea level
        CALL crash('Sea level initialisation: eustatic method not implement yet!')

      CASE ('SELEN')
        ! Sea level from SELEN
        CALL crash('Sea level initialisation: SELEN method not implement yet!')

      CASE DEFAULT
        ! Unknown case
        CALL crash('unknown choice_sealevel_model "' // &
                    TRIM( C%choice_sealevel_model) // '"!')

    END SELECT

    IF (par%master) THEN
      WRITE(*,*) 'Current sea level:', sea_level_current
    END IF

    ! Increase or decrease water column
    ! UNSURE this is correct, probably not
    DO k = 1, C%nz_ocean
      C%z_ocean(k) = sea_level_current + (k-1)*C%ocean_vertical_grid_dz
    END DO

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
    CHARACTER(LEN=256)                                    :: filename1, filename2, filename_insolation, filename_tas, filename_sos
    INTEGER                                               :: i, j, vi, k
    REAL(dp)                                              :: time_to_read_ins_PI, time_to_read_ins_LGM, time_to_read_tas, time_to_read_sos
    TYPE(type_ocean_matrix_interpolation)                 :: matrix
    INTEGER,  DIMENSION(mesh%vi1:mesh%vi2)                :: mask_ocean
    REAL(dp)                                              :: max_ocean_size, sigma
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: T_field, S_field
    TYPE(type_ice_model)                                  :: ice
    REAL(dp)                                              :: scale_snapshot_LGM_T, scale_snapshot_LGM_S, scale_snapshot_PI_T, scale_snapshot_PI_S

    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising matrix ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_matrix),'light blue') // '"...'

    ! Start and ending of simulation
    ocean%matrix%t0 = REAL(C%start_time_of_run, dp)             ! PI
    ocean%matrix%t1 = REAL(C%end_time_of_run, dp)               ! LGM

    ! Allocate memory for ocean and oceanic timeframes T and S array if not already allocated
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe0%T) .AND. .NOT. ALLOCATED(ocean%matrix%timeframe0%S)) THEN        ! PI
      ALLOCATE(ocean%matrix%timeframe0%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ALLOCATE(ocean%matrix%timeframe0%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ocean%matrix%timeframe0%T = 0._dp
      ocean%matrix%timeframe0%S = 0._dp
    END IF 
    IF (.NOT. ALLOCATED(ocean%matrix%timeframe1%T) .AND. .NOT. ALLOCATED(ocean%matrix%timeframe1%S)) THEN        ! LGM
      ALLOCATE(ocean%matrix%timeframe1%T(mesh%vi1:mesh%vi2, 1:C%nz_ocean))                                                 
      ALLOCATE(ocean%matrix%timeframe1%S(mesh%vi1:mesh%vi2, 1:C%nz_ocean))
      ocean%matrix%timeframe1%T = 0._dp
      ocean%matrix%timeframe1%S = 0._dp      
    END IF    

    ! Construct filenames
    filename1            = TRIM(C%filename_ocean_matrix_base1)  ! Snapshot PI
    filename2            = TRIM(C%filename_ocean_matrix_base2)  ! Snapshot LGM
    filename_tas         = TRIM(C%filename_tas)                 ! Anomaly field T
    filename_sos         = TRIM(C%filename_sos)                 ! Anomaly field S

    ! Set up scaling option
    scale_snapshot_LGM_T = C%scale_snapshot_LGM_T
    scale_snapshot_LGM_S = C%scale_snapshot_LGM_S
    scale_snapshot_PI_T = C%scale_snapshot_PI_T
    scale_snapshot_PI_S = C%scale_snapshot_PI_S

    ! Read the ocean snapshots
    CALL read_field_from_file_3D_ocean(filename1, 't_an', mesh, ocean%matrix%timeframe0%T) ! Ocean T PI
    CALL read_field_from_file_3D_ocean(filename1, 's_an', mesh, ocean%matrix%timeframe0%S) ! Ocean S PI
    CALL read_field_from_file_3D_ocean(filename2, 't_an', mesh, ocean%matrix%timeframe1%T) ! Ocean T LGM
    CALL read_field_from_file_3D_ocean(filename2, 's_an', mesh, ocean%matrix%timeframe1%S) ! Ocean S LGM

    ! Ensure correct model choice and load in files if necessary
    SELECT CASE (TRIM(C%choice_ocean_model_matrix))
    CASE('linear_time', 'GHG_radiative', 'GHG', 'd18O')
      ! Global weight variables or time, no specific intialisation necessary

    CASE('anomaly_field')
      ! Allocate arrays for anomaly fields

      IF (.NOT. ALLOCATED(ocean%matrix%tas) .AND. .NOT. ALLOCATED(ocean%matrix%sos)) THEN
        ALLOCATE(ocean%matrix%tas(mesh%vi1:mesh%vi2))
        ALLOCATE(ocean%matrix%sos(mesh%vi1:mesh%vi2))
      END IF 

      ! Read anomaly fields
      CALL read_field_from_file_2D(filename_tas, 'tas', mesh, ocean%matrix%tas) ! Ocean T
      CALL read_field_from_file_2D(filename_sos, 'sos', mesh, ocean%matrix%sos) ! Ocean S 
      WRITE(*, *) 'Lowest T anomaly:', MINVAL(ocean%matrix%tas)
      WRITE(*, *) 'Highest T anomaly:', MAXVAL(ocean%matrix%tas)
      WRITE(*, *) 'Lowest S anomaly:', MINVAL(ocean%matrix%sos)
      WRITE(*, *) 'Highest S anomaly:', MAXVAL(ocean%matrix%sos)

      ! Determine timeframe0 to be PI (LGM + anomaly)
      DO vi = mesh%vi1, mesh%vi2
        ocean%matrix%timeframe0%T(vi,1) = ocean%matrix%timeframe1%T(vi,1) - ocean%matrix%tas(vi)
        ocean%matrix%timeframe0%S(vi,1) = ocean%matrix%timeframe1%S(vi,1) - ocean%matrix%sos(vi)
      END DO

      WRITE(*, *) 'Anomaly field min/max:'
      PRINT *, "T min/max:", MINVAL(ocean%matrix%tas), MAXVAL(ocean%matrix%tas)
      DO k = 1, C%nz_ocean
        PRINT *, "Layer timeframe0", k, "T min/max:", MINVAL(ocean%matrix%timeframe0%T(:,k)), MAXVAL(ocean%matrix%timeframe0%T(:,k))
      END DO
      
    CASE('insolation', 'insolation+GHG', 'insolation+GHG_radiative')

      CALL initialise_insolation_fields(mesh, ocean, matrix)

    CASE DEFAULT
      CALL crash('Unknown choice_ocean_model_matrix' // TRIM(C%choice_ocean_model_matrix))
    END SELECT

    ! Prescribe initial ocean state (start from PI snapshot)
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, C%nz_ocean
        ocean%T(vi,k) = MAX(-3.0_dp, ocean%matrix%timeframe0%T(vi,k) + scale_snapshot_PI_T)
        ocean%S(vi,k) = ocean%matrix%timeframe0%S(vi,k) + scale_snapshot_PI_S
      END DO
    END DO

    ! Options to scale snapshot LGM
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, C%nz_ocean
        ocean%matrix%timeframe1%T(vi,k) = MAX(-3.0_dp, ocean%matrix%timeframe1%T(vi,k) + scale_snapshot_LGM_T)
        ocean%matrix%timeframe1%S(vi,k) = ocean%matrix%timeframe1%S(vi,k) + scale_snapshot_LGM_S
      END DO
    END DO

    ! Print snapshot info
    IF (par%master) THEN
      WRITE(*, *) 'Snapshot min/max over depth:'
      DO k = 1, C%nz_ocean
        PRINT *, "Layer timeframe0", k, "T min/max:", MINVAL(ocean%matrix%timeframe0%T(:,k)), MAXVAL(ocean%matrix%timeframe0%T(:,k))
      END DO  
      DO k = 1, C%nz_ocean
        PRINT *, "Layer timeframe1", k, "T min/max:", MINVAL(ocean%matrix%timeframe1%T(:,k)), MAXVAL(ocean%matrix%timeframe1%T(:,k))
      END DO   
      DO k = 1, C%nz_ocean
        PRINT *, "Layer ocean%T", k, "T min/max:", MINVAL(ocean%T(:,k)), MAXVAL(ocean%T(:,k))
      END DO   
      WRITE(*, *) 'Snapshot difference:'
      WRITE(*, *) 'Max T difference:', MAXVAL(ABS(ocean%matrix%timeframe0%T - ocean%matrix%timeframe1%T))
      WRITE(*, *) 'Max S difference:', MAXVAL(ABS(ocean%matrix%timeframe0%S - ocean%matrix%timeframe1%S))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_ocean_model_matrix

  SUBROUTINE initialise_insolation_fields(mesh, ocean, matrix)
    ! Initialize the insolation fields Q_TOA_PI and Q_TOA_LGM for the ocean model.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ocean_matrix_interpolation),  INTENT(IN)    :: matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_insolation_fields'
    CHARACTER(LEN=256)                                    :: filename_insolation
    REAL(dp)                                              :: time_to_read_ins_PI, time_to_read_ins_LGM

    ! Begin initialization routine
    CALL init_routine(routine_name)

    filename_insolation  = TRIM(C%filename_insolation) 
    time_to_read_ins_PI  = 0_dp                                ! Insolation PI
    time_to_read_ins_LGM = -21000_dp                           ! Insolation LGM

    ! Allocate arrays for insolation at PI and LGM
    IF (.NOT. ALLOCATED(ocean%matrix%Q_TOA_PI)) THEN
      ALLOCATE(ocean%matrix%Q_TOA_PI(mesh%vi1:mesh%vi2))
    END IF
    IF (.NOT. ALLOCATED(ocean%matrix%Q_TOA_LGM)) THEN
      ALLOCATE(ocean%matrix%Q_TOA_LGM(mesh%vi1:mesh%vi2))
    END IF

    ! Read fields at PI and LGM, starts at 2000 AD (Laskar solution starts at 2000 AD and is in 1ka timesteps)
    CALL read_field_from_file_2D(filename_insolation, 'Q_TOA', mesh, ocean%matrix%Q_TOA_PI, time_to_read_ins_PI)   ! PI  
    CALL read_field_from_file_2D(filename_insolation, 'Q_TOA', mesh, ocean%matrix%Q_TOA_LGM, time_to_read_ins_LGM) ! LGM 

    WRITE(*, *) 'Lowest Q_TOA_PI:', MINVAL(ocean%matrix%Q_TOA_PI)
    WRITE(*, *) 'Highest Q_TOA_PI:', MAXVAL(ocean%matrix%Q_TOA_PI)
    WRITE(*, *) 'Lowest Q_TOA_LGM:', MINVAL(ocean%matrix%Q_TOA_LGM)
    WRITE(*, *) 'Highest Q_TOA_LGM:', MAXVAL(ocean%matrix%Q_TOA_LGM)

    ! Finalize the initialization routine
    CALL finalise_routine(routine_name)

  END SUBROUTINE initialise_insolation_fields

END MODULE ocean_matrix