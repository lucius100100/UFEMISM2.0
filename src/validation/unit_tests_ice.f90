MODULE unit_tests_ice

  ! Unit tests for different ice velocity solvers.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE scalar_types                                           , ONLY: type_regional_scalars
  USE reference_geometries                                   , ONLY: type_reference_geometry, initialise_reference_geometries_raw, initialise_reference_geometries_on_model_mesh
  USE mesh_creation                                          , ONLY: create_mesh_from_gridded_geometry, write_mesh_success
  USE ice_model_main                                         , ONLY: initialise_ice_model, calc_zeta_gradients
  USE ice_velocity_main                                      , ONLY: initialise_velocity_solver, solve_stress_balance
  USE mesh_operators                                         , ONLY: calc_3D_matrix_operators_mesh
  USE netcdf_basic                                           , ONLY: create_new_netcdf_file_for_writing, open_existing_netcdf_file_for_writing, close_netcdf_file
  USE netcdf_output                                          , ONLY: add_zeta_dimension_to_file, setup_mesh_in_netcdf_file, add_field_mesh_dp_3D_b_notime, &
                                                                     write_to_field_multopt_mesh_dp_3D_b_notime, add_field_mesh_dp_2D_b_notime, &
                                                                     write_to_field_multopt_mesh_dp_2D_b_notime, add_time_dimension_to_file, add_field_mesh_dp_2D, &
                                                                     add_field_mesh_dp_2D_b, write_time_to_file, write_to_field_multopt_mesh_dp_2D, &
                                                                     write_to_field_multopt_mesh_dp_2D_b
  USE mesh_utilities                                         , ONLY: find_containing_triangle, integrate_over_domain, average_over_domain
  USE ice_thickness                                          , ONLY: calc_Hi_tplusdt
  USE math_utilities                                         , ONLY: ice_surface_elevation
  USE analytical_solutions                                   , ONLY: Halfar_dome

  IMPLICIT NONE

! ===== Global variables =====
! ============================

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE run_all_ice_unit_tests
    ! Run all unit tests for the different ice velocity solvers

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'run_all_ice_unit_tests'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run all ice-dynamics-related unit tests
    CALL test_ice_velocities_Halfar_dome
    CALL test_ISMIP_HOM_all
    CALL test_thickness_evolution_Halfar_dome

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_all_ice_unit_tests

  ! ===== Ice velocities in the Halfar dome geometry =====
  ! ======================================================

  SUBROUTINE test_ice_velocities_Halfar_dome
    ! Generate a simple Halfar dome geometry and calculate ice velocities
    ! with the SIA, DIVA, and BPA to check if the results make sense.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ice_velocities_Halfar_dome'
    TYPE(type_mesh)                                                    :: mesh
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars
    TYPE(type_reference_geometry)                                      :: refgeo_init, refgeo_PD, refgeo_GIAeq
    CHARACTER(LEN=256)                                                 :: region_name, mesh_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_SIA , v_3D_b_SIA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_DIVA, v_3D_b_DIVA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_BPA , v_3D_b_BPA
    LOGICAL                                                            :: found_errors_SIA, found_errors_DIVA, found_errors_BPA
    REAL(dp), PARAMETER                                                :: R_check         = 440E3_dp  ! At this distance from the ice divide, uabs_surf should be equal to...
    REAL(dp), PARAMETER                                                :: uabs_surf_check = 105._dp   ! ...this value...
    REAL(dp), PARAMETER                                                :: uabs_surf_tol   = 20._dp    ! ...within this tolerance.
    INTEGER                                                            :: ti
    REAL(dp)                                                           :: R, uabs_surf_target, uabs_surf_SIA, uabs_surf_DIVA, uabs_surf_BPA
    CHARACTER(LEN=256)                                                 :: filename
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Model domain
    region_name = 'ANT'

    ! Set model configuration for this experiment
    CALL set_config_for_Halfar_dome

  ! Initialise the model
  ! ====================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Create mesh from gridded initial geometry data
    mesh_name = 'Halfar_mesh'
    CALL create_mesh_from_gridded_geometry( region_name, mesh_name, &
      refgeo_init%grid_raw, &
      refgeo_init%Hi_grid_raw, &
      refgeo_init%Hb_grid_raw, &
      refgeo_init%Hs_grid_raw, &
      refgeo_init%SL_grid_raw, &
      C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
      mesh)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Initialise the ice model
    C%choice_stress_balance_approximation = 'SIA'
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Also initialise DIVA and BPA solvers
    C%choice_stress_balance_approximation = 'DIVA'
    CALL initialise_velocity_solver( mesh, ice, region_name)
    C%choice_stress_balance_approximation = 'BPA'
    CALL initialise_velocity_solver( mesh, ice, region_name)

    ! Calculate necessary 3-D matrix operators
    CALL calc_zeta_gradients( mesh, ice)
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

  ! Calculate ice velocities
  ! ========================

    ! Allocate memory
    ALLOCATE( u_3D_b_SIA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_SIA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_DIVA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_DIVA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_BPA(  mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_BPA(  mesh%ti1:mesh%ti2, mesh%nz))

    ! Calculate velocities

    ! SIA
    IF (par%master) WRITE(0,*) '   Calculating ice velocities for the Halfar dome with the ' // colour_string( 'SIA','light blue') // '...'
    C%choice_stress_balance_approximation = 'SIA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_SIA = ice%u_3D_b
    v_3D_b_SIA = ice%v_3D_b

    ! DIVA
    IF (par%master) WRITE(0,*) '   Calculating ice velocities for the Halfar dome with the ' // colour_string( 'DIVA','light blue') // '...'
    C%choice_stress_balance_approximation = 'DIVA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_DIVA = ice%u_3D_b
    v_3D_b_DIVA = ice%v_3D_b

    ! BPA
    IF (par%master) WRITE(0,*) '   Calculating ice velocities for the Halfar dome with the ' // colour_string( 'BPA','light blue') // '...'
    C%choice_stress_balance_approximation = 'BPA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_BPA = ice%u_3D_b
    v_3D_b_BPA = ice%v_3D_b

    ! Validate results - check if the surface velocity follows the expected
    ! linear increase away from the ice divide

    found_errors_SIA  = .FALSE.
    found_errors_DIVA = .FALSE.
    found_errors_BPA  = .FALSE.

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate distance from the ice divide
      R = NORM2( mesh%TriGC( ti,:))

      ! Only check up to a certain distance (near the margin, the SIA and BPA become a bit "noisy")
      IF (R > R_check) CYCLE

      ! Calculate target surface velocity
      uabs_surf_target = uabs_surf_check * R / R_check

      ! Calculate actual surface velocities
      uabs_surf_SIA  = SQRT( u_3D_b_SIA(  ti,1)**2 + v_3D_b_SIA(  ti,1)**2)
      uabs_surf_DIVA = SQRT( u_3D_b_DIVA( ti,1)**2 + v_3D_b_DIVA( ti,1)**2)
      uabs_surf_BPA  = SQRT( u_3D_b_BPA(  ti,1)**2 + v_3D_b_BPA(  ti,1)**2)

      ! Compare
      IF (ABS( uabs_surf_SIA  - uabs_surf_target) > uabs_surf_tol) found_errors_SIA  = .TRUE.
      IF (ABS( uabs_surf_DIVA - uabs_surf_target) > uabs_surf_tol) found_errors_DIVA = .TRUE.
      IF (ABS( uabs_surf_BPA  - uabs_surf_target) > uabs_surf_tol) found_errors_BPA  = .TRUE.

    END DO ! DO ti = mesh%ti1, mesh%ti2

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_SIA , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_DIVA, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_BPA , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    IF (.NOT. (found_errors_SIA .OR. found_errors_DIVA .OR. found_errors_BPA)) THEN
      IF (par%master) CALL happy('validated SIA, DIVA, and BPA solvers for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_SIA) THEN
      IF (par%master) CALL warning('found errors in SIA solver for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_DIVA) THEN
      IF (par%master) CALL warning('found errors in DIVA solver for diagnostic velocities in the Halfar dome geometry')
    END IF
    IF (found_errors_BPA) THEN
      IF (par%master) CALL warning('found errors in BPA solver for diagnostic velocities in the Halfar dome geometry')
    END IF

  ! Write results to NetCDF
  ! =======================

    ! Create a NetCDF output file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add a zeta dimension for the 3-D ice velocities
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the ice velocities as fields
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_SIA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_SIA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_DIVA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_DIVA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_BPA' )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_BPA' )

    ! Write the velocities to the file
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_SIA' , u_3D_b_SIA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_SIA' , v_3D_b_SIA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_DIVA', u_3D_b_DIVA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_DIVA', v_3D_b_DIVA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_BPA' , u_3D_b_BPA )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_BPA' , v_3D_b_BPA )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( u_3D_b_SIA)
    DEALLOCATE( v_3D_b_SIA)
    DEALLOCATE( u_3D_b_DIVA)
    DEALLOCATE( v_3D_b_DIVA)
    DEALLOCATE( u_3D_b_BPA)
    DEALLOCATE( v_3D_b_BPA)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ice_velocities_Halfar_dome

  SUBROUTINE set_config_for_Halfar_dome
    ! Set the config for the Halfar dome experiment

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'set_config_for_Halfar_dome'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Model domain
    C%lambda_M_ANT                          = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
    C%phi_M_ANT                             = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
    C%beta_stereo_ANT                       = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
    C%xmin_ANT                              = -800E3_dp                        ! Western  boundary     of the Antarctica domain [m]
    C%xmax_ANT                              =  800E3_dp                        ! Eastern  boundary     of the Antarctica domain [m]
    C%ymin_ANT                              = -800E3_dp                        ! Southern boundary     of the Antarctica domain [m]
    C%ymax_ANT                              =  800E3_dp                        ! Northern boundary     of the Antarctica domain [m]

    ! == Reference geometries (initial, present-day, and GIA equilibrium)

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                         = 2._dp                            ! Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    C%remove_Lake_Vostok                    = .FALSE.

    ! == Initial geometry
    C%choice_refgeo_init_ANT                = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    C%choice_refgeo_init_idealised          = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_init_idealised              = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == Present-day geometry
    C%choice_refgeo_PD_ANT                  = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    C%choice_refgeo_PD_idealised            = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_PD_idealised                = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == GIA equilibrium geometry
    C%choice_refgeo_GIAeq_ANT               = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    C%choice_refgeo_GIAeq_idealised         = 'Halfar'                         ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_GIAeq_idealised             = 5000._dp                         ! Resolution of square grid used for idealised present-day geometry

    ! == Parameters for idealised geometries
    C%refgeo_idealised_Halfar_H0            = 3000._dp                         ! Suggested value: 3000 m
    C%refgeo_idealised_Halfar_R0            = 500E3_dp                         ! Suggested value: 500E3 m

    ! == Mesh generation

    ! How to set up the initial mesh
    C%choice_initial_mesh_ANT               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'

    ! Resolutions for different parts of the ice sheet
    C%maximum_resolution_uniform            = 100e3_dp                         ! [m]          Maximum resolution for the entire domain
    C%maximum_resolution_grounded_ice       = 20e3_dp                          ! [m]          Maximum resolution for grounded ice
    C%maximum_resolution_floating_ice       = 100e3_dp                         ! [m]          Maximum resolution for floating ice
    C%maximum_resolution_grounding_line     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    C%grounding_line_width                  = 100e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    C%maximum_resolution_calving_front      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    C%calving_front_width                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    C%maximum_resolution_ice_front          = 10e3_dp                          ! [m]          Maximum resolution for the ice front
    C%ice_front_width                       = 5e3_dp                           ! [m]          Width of the band around the ice front that should get this resolution
    C%maximum_resolution_coastline          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    C%coastline_width                       = 100e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

    ! Regions of interest
    C%choice_regions_of_interest            = ''                               ! Regions of interest where other (higher) resolutions apply. Separated by double vertical bars "||", e.g. "PineIsland||Thwaites"

    ! Advanced geometry parameters
    C%do_singlecore_mesh_creation           = .TRUE.                           !              Whether or not to use only a single core for mesh generation (for better reproducibility)
    C%alpha_min                             = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle
    C%nit_Lloyds_algorithm                  = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement
    C%mesh_resolution_tolerance             = 1.25_dp                          ! [-]          Factors the target resolution for trangle-size requirement. 1=strict, use >1 to avoid unnecesarily high resolution

    ! Memory
    C%nC_mem                                = 32                               ! [-]          How many columns of memory should be allocated for connectivity lists

    ! == The scaled vertical coordinate zeta

    C%choice_zeta_grid                      = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    C%nz                                    = 12                               ! The number of vertical layers to use

    ! == Ice dynamics - velocity

    ! General
    C%choice_stress_balance_approximation   = 'SIA'                            ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA"
    C%n_flow                                = 3._dp                            ! Exponent in Glen's flow law
    C%m_enh_sheet                           = 1._dp                            ! Ice flow enhancement factor for grounded ice
    C%m_enh_shelf                           = 1._dp                            ! Ice flow enhancement factor for floating ice
    C%choice_hybrid_SIASSA_scheme           = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    C%do_GL_subgrid_friction                = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    C%subgrid_friction_exponent             = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta
    C%do_include_SSADIVA_crossterms         = .TRUE.                           ! Whether or not to include the gradients of the effective viscosity (the "cross-terms") in the solution of the SSA/DIVA

    ! Initialisation
    C%choice_initial_velocity_ANT           = 'zero'

    ! Some parameters for numerically solving the stress balance
    C%SIA_maximum_diffusivity               = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    C%visc_it_norm_dUV_tol                  = 5E-6_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    C%visc_it_nit                           = 500._dp                          ! Maximum number of effective viscosity iterations
    C%visc_it_relax                         = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    C%epsilon_sq_0                          = 1E-10_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    C%visc_eff_min                          = 1E0_dp                           ! Minimum value for effective viscosity
    C%beta_max                              = 1E20_dp                          ! Maximum value for basal friction coefficient
    C%vel_max                               = 5000._dp                         ! Velocities are limited to this value
    C%stress_balance_PETSc_rtol             = 1E-6_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    C%stress_balance_PETSc_abstol           = 1E-4_dp                          ! PETSc solver - stop criterion, absolute difference

    ! == Ice dynamics - sliding law

    ! Sliding laws
    C%choice_sliding_law                    = 'no_sliding'                     ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"

    ! == Ice dynamics - boundary conditions

    C%BC_u_west                             = 'zero'                           ! Boundary conditions for the ice velocity field at the domain border
    C%BC_u_east                             = 'zero'                           ! Allowed choices: "infinite", "zero", "zero"
    C%BC_u_south                            = 'zero'
    C%BC_u_north                            = 'zero'
    C%BC_v_west                             = 'zero'
    C%BC_v_east                             = 'zero'
    C%BC_v_south                            = 'zero'
    C%BC_v_north                            = 'zero'
    C%BC_H_west                             = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    C%BC_H_east                             = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    C%BC_H_south                            = 'zero'
    C%BC_H_north                            = 'zero'

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    C%choice_ice_rheology                   = 'uniform'                        ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
    C%uniform_flow_factor                   = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_config_for_Halfar_dome

  ! ===== ISMIP-HOM =====
  ! =====================

  SUBROUTINE test_ISMIP_HOM_all
    ! Run and validate all ISMIP-HOM experiments
    !
    ! The Ice-Sheet Model Intercomparison Project for Higher-Order Models (ISMIP-HOM)
    ! consists of 6 experiments, called A to E. Experiments A, B, C, and D concern
    ! diagnostic velocities in an idealised geometry. Experiment E consists of diagnostic
    ! velocities in a realistic 1-D geometry (which is skipped here, as UFEMISM doesn't
    ! lend itself well to 1-D experiments). Experiment F consists of finding a steady-state
    ! ice-sheet in an idealised geometry.
    !
    ! See also: Pattyn et al., 2008: Benchmark experiments for higher-order and full-Stokes
    !           ice sheet models (ISMIP-HOM), The Cryosphere 2, 95-108

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ISMIP_HOM_all'
    REAL(dp), DIMENSION(6), PARAMETER                                  :: Ls = [160E3_dp, 80E3_dp, 40E3_dp, 20E3_dp, 10E3_dp, 5E3_dp]
    INTEGER                                                            :: li
    REAL(dp)                                                           :: L

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Generate a new mesh for each scale, then run all 6 experiments on that same mesh
    DO li = 1, 6

      ! Experiment scale
      L = Ls( li)

      ! Run the ISMIP-HOM experiments
      CALL test_ISMIP_HOM_A( L)
      CALL test_ISMIP_HOM_C( L)

    END DO ! DO li = 1, 6

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ISMIP_HOM_all

  SUBROUTINE test_ISMIP_HOM_A( L)
    ! Run and validate ISMIP-HOM experiment A at this length scale

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                                            INTENT(IN)    :: L

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ISMIP_HOM_A'
    TYPE(type_mesh)                                                    :: mesh
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars
    TYPE(type_reference_geometry)                                      :: refgeo_init, refgeo_PD, refgeo_GIAeq
    CHARACTER(LEN=256)                                                 :: region_name, mesh_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_SIASSA, v_3D_b_SIASSA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_DIVA  , v_3D_b_DIVA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_BPA   , v_3D_b_BPA
    LOGICAL                                                            :: found_errors_SIASSA, found_errors_DIVA, found_errors_BPA
    REAL(dp), PARAMETER                                                :: u_surf_min_check = 1.8_dp
    REAL(dp), PARAMETER                                                :: u_surf_min_tol   = 0.5_dp
    REAL(dp), PARAMETER                                                :: u_surf_max_check = 100.0_dp
    REAL(dp), PARAMETER                                                :: u_surf_max_tol   = 25.0_dp
    INTEGER                                                            :: ti
    REAL(dp), DIMENSION(2)                                             :: p
    REAL(dp)                                                           :: u_surf_min_SIASSA, u_surf_min_DIVA, u_surf_min_BPA
    REAL(dp)                                                           :: u_surf_max_SIASSA, u_surf_max_DIVA, u_surf_max_BPA
    CHARACTER(LEN=256)                                                 :: filename
    CHARACTER(LEN=256)                                                 :: filename_ext
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Filename extension for experiment scale
    IF     (L == 160E3_dp) THEN
      filename_ext = '_160km'
    ELSEIF (L == 80E3_dp) THEN
      filename_ext = '_80km'
    ELSEIF (L == 40E3_dp) THEN
      filename_ext = '_40km'
    ELSEIF (L == 20E3_dp) THEN
      filename_ext = '_20km'
    ELSEIF (L == 10E3_dp) THEN
      filename_ext = '_10km'
    ELSEIF (L == 5E3_dp) THEN
      filename_ext = '_5km'
    ELSE
      CALL crash('invalid ISMIP-HOM length scale L = {dp_01} km', dp_01 = L)
    END IF

    region_name = 'ANT'

    ! Set up configuration for this experiment
    CALL set_config_for_ISMIP_HOM( L)

    C%choice_refgeo_init_idealised          = 'ISMIP-HOM_A'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_refgeo_PD_idealised            = 'ISMIP-HOM_A'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_refgeo_GIAeq_idealised         = 'ISMIP-HOM_A'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_sliding_law                    = 'no_sliding'                     ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"

  ! Initialise the model
  ! ====================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Create mesh from gridded initial geometry data
    mesh_name = 'ISMIP_HOM_A' // TRIM( filename_ext) // '_mesh'
    CALL create_mesh_from_gridded_geometry( region_name, mesh_name, &
      refgeo_init%grid_raw, &
      refgeo_init%Hi_grid_raw, &
      refgeo_init%Hb_grid_raw, &
      refgeo_init%Hs_grid_raw, &
      refgeo_init%SL_grid_raw, &
      C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
      mesh)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Initialise the ice model
    C%choice_stress_balance_approximation = 'SIA/SSA'
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Also initialise DIVA and BPA solvers
    C%choice_stress_balance_approximation = 'DIVA'
    CALL initialise_velocity_solver( mesh, ice, region_name)
    C%choice_stress_balance_approximation = 'BPA'
    CALL initialise_velocity_solver( mesh, ice, region_name)

    ! Calculate necessary 3-D matrix operators
    CALL calc_zeta_gradients( mesh, ice)
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

  ! Calculate ice velocities
  ! ========================

    ! Allocate memory
    ALLOCATE( u_3D_b_SIASSA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_SIASSA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_DIVA(   mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_DIVA(   mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_BPA(    mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_BPA(    mesh%ti1:mesh%ti2, mesh%nz))

    ! Calculate velocities

    ! SIA/SSA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM A with the ' // colour_string( 'SIA/SSA','light blue') // '...'
    C%choice_stress_balance_approximation = 'SIA/SSA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_SIASSA = ice%u_3D_b
    v_3D_b_SIASSA = ice%v_3D_b

    ! DIVA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM A with the ' // colour_string( 'DIVA','light blue') // '...'
    C%choice_stress_balance_approximation = 'DIVA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_DIVA = ice%u_3D_b
    v_3D_b_DIVA = ice%v_3D_b

    ! BPA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM A with the ' // colour_string( 'BPA','light blue') // '...'
    C%choice_stress_balance_approximation = 'BPA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_BPA = ice%u_3D_b
    v_3D_b_BPA = ice%v_3D_b

  ! ===== Validate results =====
  ! ============================

    found_errors_SIASSA = .FALSE.
    found_errors_DIVA   = .FALSE.
    found_errors_BPA    = .FALSE.

    ! Look, I'm not going to manually define desired min/max/tol velocities for 4 experiments,
    ! each at 6 scales, for 3 different stress balance approximations.
    ! Instead, we'll look only at the 16 km version, where all three approximations give
    ! roughly the same answer.

    IF (L == 160E3_dp) THEN

    ! == Minimum surface velocity
    ! ===========================

      ! Find the triangle that should contain the minimum velocity
      p = [0.25_dp, 0.25_dp] * L
      ti = 1
      CALL find_containing_triangle( mesh, p, ti)

      ! IF this triangle lies in this process' domain, obtain velocities and compare to target
      IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

        ! Obtain surface velocities
        u_surf_min_SIASSA = u_3D_b_SIASSA( ti,1)
        u_surf_min_DIVA   = u_3D_b_DIVA(   ti,1)
        u_surf_min_BPA    = u_3D_b_BPA(    ti,1)

        ! Compare to target
        IF (ABS( u_surf_min_SIASSA - u_surf_min_check) > u_surf_min_tol) found_errors_SIASSA = .TRUE.
        IF (ABS( u_surf_min_DIVA   - u_surf_min_check) > u_surf_min_tol) found_errors_DIVA   = .TRUE.
        IF (ABS( u_surf_min_BPA    - u_surf_min_check) > u_surf_min_tol) found_errors_BPA    = .TRUE.

      END IF ! IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

    ! == Maximum surface velocity
    ! ===========================

      ! Find the triangle that should contain the maximum velocity
      p = [-0.25_dp, 0.25_dp] * L
      ti = 1
      CALL find_containing_triangle( mesh, p, ti)

      ! IF this triangle lies in this process' domain, obtain velocities and compare to target
      IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

        ! Obtain surface velocities
        u_surf_max_SIASSA = u_3D_b_SIASSA( ti,1)
        u_surf_max_DIVA   = u_3D_b_DIVA(   ti,1)
        u_surf_max_BPA    = u_3D_b_BPA(    ti,1)

        ! Compare to target
        IF (ABS( u_surf_max_SIASSA - u_surf_max_check) > u_surf_max_tol) found_errors_SIASSA = .TRUE.
        IF (ABS( u_surf_max_DIVA   - u_surf_max_check) > u_surf_max_tol) found_errors_DIVA   = .TRUE.
        IF (ABS( u_surf_max_BPA    - u_surf_max_check) > u_surf_max_tol) found_errors_BPA    = .TRUE.

      END IF ! IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

    END IF ! IF (L == 160E3_dp) THEN

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_SIASSA, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_DIVA  , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_BPA   , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    IF (.NOT. (found_errors_SIASSA .OR. found_errors_DIVA .OR. found_errors_BPA)) THEN
      IF (par%master) CALL happy('validated SIA/SSA, DIVA, and BPA solvers for diagnostic velocities in the ISMIP-HOM A geometry')
    END IF
    IF (found_errors_SIASSA) THEN
      IF (par%master) CALL warning('found errors in SIA/SSA solver for diagnostic velocities in the ISMIP-HOM A geometry')
    END IF
    IF (found_errors_DIVA) THEN
      IF (par%master) CALL warning('found errors in DIVA solver for diagnostic velocities in the ISMIP-HOM A geometry')
    END IF
    IF (found_errors_BPA) THEN
      IF (par%master) CALL warning('found errors in BPA solver for diagnostic velocities in the ISMIP-HOM A geometry')
    END IF

  ! Write results to NetCDF
  ! =======================

    ! Create a NetCDF output file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // TRIM( filename_ext) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add a zeta dimension for the 3-D ice velocities
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the ice velocities as fields
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_SIASSA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_SIASSA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_DIVA'  )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_DIVA'  )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_BPA'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_BPA'   )

    ! Write the velocities to the file
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_SIASSA', u_3D_b_SIASSA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_SIASSA', v_3D_b_SIASSA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_DIVA'  , u_3D_b_DIVA  )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_DIVA'  , v_3D_b_DIVA  )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_BPA'   , u_3D_b_BPA   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_BPA'   , v_3D_b_BPA   )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( u_3D_b_SIASSA)
    DEALLOCATE( v_3D_b_SIASSA)
    DEALLOCATE( u_3D_b_DIVA)
    DEALLOCATE( v_3D_b_DIVA)
    DEALLOCATE( u_3D_b_BPA)
    DEALLOCATE( v_3D_b_BPA)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ISMIP_HOM_A

  SUBROUTINE test_ISMIP_HOM_C( L)
    ! Run and validate ISMIP-HOM experiment C at this length scale

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                                            INTENT(IN)    :: L

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_ISMIP_HOM_C'
    TYPE(type_mesh)                                                    :: mesh
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars
    TYPE(type_reference_geometry)                                      :: refgeo_init, refgeo_PD, refgeo_GIAeq
    CHARACTER(LEN=256)                                                 :: region_name, mesh_name
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_SIASSA, v_3D_b_SIASSA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_DIVA  , v_3D_b_DIVA
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                            :: u_3D_b_BPA   , v_3D_b_BPA
    LOGICAL                                                            :: found_errors_SIASSA, found_errors_DIVA, found_errors_BPA
    REAL(dp), PARAMETER                                                :: u_surf_min_check = 10.5_dp
    REAL(dp), PARAMETER                                                :: u_surf_min_tol   = 2.0_dp
    REAL(dp), PARAMETER                                                :: u_surf_max_check = 125.0_dp
    REAL(dp), PARAMETER                                                :: u_surf_max_tol   = 10.0_dp
    INTEGER                                                            :: ti
    REAL(dp), DIMENSION(2)                                             :: p
    REAL(dp)                                                           :: u_surf_min_SIASSA, u_surf_min_DIVA, u_surf_min_BPA
    REAL(dp)                                                           :: u_surf_max_SIASSA, u_surf_max_DIVA, u_surf_max_BPA
    CHARACTER(LEN=256)                                                 :: filename
    CHARACTER(LEN=256)                                                 :: filename_ext
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Filename extension for experiment scale
    IF     (L == 160E3_dp) THEN
      filename_ext = '_160km'
    ELSEIF (L == 80E3_dp) THEN
      filename_ext = '_80km'
    ELSEIF (L == 40E3_dp) THEN
      filename_ext = '_40km'
    ELSEIF (L == 20E3_dp) THEN
      filename_ext = '_20km'
    ELSEIF (L == 10E3_dp) THEN
      filename_ext = '_10km'
    ELSEIF (L == 5E3_dp) THEN
      filename_ext = '_5km'
    ELSE
      CALL crash('invalid ISMIP-HOM length scale L = {dp_01} km', dp_01 = L)
    END IF

    region_name = 'ANT'

    ! Set up configuration for this experiment
    CALL set_config_for_ISMIP_HOM( L)

    C%choice_refgeo_init_idealised          = 'ISMIP-HOM_C'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_refgeo_PD_idealised            = 'ISMIP-HOM_C'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_refgeo_GIAeq_idealised         = 'ISMIP-HOM_C'                    ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%choice_sliding_law                    = 'idealised'                      ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
    C%choice_idealised_sliding_law          = 'ISMIP-HOM_C'                    ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"

  ! Initialise the model
  ! ====================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Create mesh from gridded initial geometry data
    mesh_name = 'ISMIP_HOM_C' // TRIM( filename_ext) // '_mesh'
    CALL create_mesh_from_gridded_geometry( region_name, mesh_name, &
      refgeo_init%grid_raw, &
      refgeo_init%Hi_grid_raw, &
      refgeo_init%Hb_grid_raw, &
      refgeo_init%Hs_grid_raw, &
      refgeo_init%SL_grid_raw, &
      C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
      mesh)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Initialise the ice model
    C%choice_stress_balance_approximation = 'SIA/SSA'
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Also initialise DIVA and BPA solvers
    C%choice_stress_balance_approximation = 'DIVA'
    CALL initialise_velocity_solver( mesh, ice, region_name)
    C%choice_stress_balance_approximation = 'BPA'
    CALL initialise_velocity_solver( mesh, ice, region_name)

    ! Calculate necessary 3-D matrix operators
    CALL calc_zeta_gradients( mesh, ice)
    CALL calc_3D_matrix_operators_mesh( mesh, ice)

  ! Calculate ice velocities
  ! ========================

    ! Allocate memory
    ALLOCATE( u_3D_b_SIASSA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_SIASSA( mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_DIVA(   mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_DIVA(   mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( u_3D_b_BPA(    mesh%ti1:mesh%ti2, mesh%nz))
    ALLOCATE( v_3D_b_BPA(    mesh%ti1:mesh%ti2, mesh%nz))

    ! Calculate velocities

    ! SIA/SSA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM C with the ' // colour_string( 'SIA/SSA','light blue') // '...'
    C%choice_stress_balance_approximation = 'SIA/SSA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_SIASSA = ice%u_3D_b
    v_3D_b_SIASSA = ice%v_3D_b

    ! DIVA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM C with the ' // colour_string( 'DIVA','light blue') // '...'
    C%choice_stress_balance_approximation = 'DIVA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_DIVA = ice%u_3D_b
    v_3D_b_DIVA = ice%v_3D_b

    ! BPA
    IF (par%master) WRITE(0,*) '  Calculating ice velocities for ISMIP-HOM C with the ' // colour_string( 'BPA','light blue') // '...'
    C%choice_stress_balance_approximation = 'BPA'
    CALL solve_stress_balance( mesh, ice)
    u_3D_b_BPA = ice%u_3D_b
    v_3D_b_BPA = ice%v_3D_b

  ! ===== Validate results =====
  ! ============================

    found_errors_SIASSA = .FALSE.
    found_errors_DIVA   = .FALSE.
    found_errors_BPA    = .FALSE.

    ! Look, I'm not going to manually define desired min/max/tol velocities for 4 experiments,
    ! each at 6 scales, for 3 different stress balance approximations.
    ! Instead, we'll look only at the 16 km version, where all three approximations give
    ! roughly the same answer.

    IF (L == 160E3_dp) THEN

    ! == Minimum surface velocity
    ! ===========================

      ! Find the triangle that should contain the minimum velocity
      p = [0.25_dp, 0.25_dp] * L
      ti = 1
      CALL find_containing_triangle( mesh, p, ti)

      ! IF this triangle lies in this process' domain, obtain velocities and compare to target
      IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

        ! Obtain surface velocities
        u_surf_min_SIASSA = u_3D_b_SIASSA( ti,1)
        u_surf_min_DIVA   = u_3D_b_DIVA(   ti,1)
        u_surf_min_BPA    = u_3D_b_BPA(    ti,1)

        ! Compare to target
        IF (ABS( u_surf_min_SIASSA - u_surf_min_check) > u_surf_min_tol) found_errors_SIASSA = .TRUE.
        IF (ABS( u_surf_min_DIVA   - u_surf_min_check) > u_surf_min_tol) found_errors_DIVA   = .TRUE.
        IF (ABS( u_surf_min_BPA    - u_surf_min_check) > u_surf_min_tol) found_errors_BPA    = .TRUE.

      END IF ! IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

    ! == Maximum surface velocity
    ! ===========================

      ! Find the triangle that should contain the maximum velocity
      p = [-0.25_dp, 0.25_dp] * L
      ti = 1
      CALL find_containing_triangle( mesh, p, ti)

      ! IF this triangle lies in this process' domain, obtain velocities and compare to target
      IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

        ! Obtain surface velocities
        u_surf_max_SIASSA = u_3D_b_SIASSA( ti,1)
        u_surf_max_DIVA   = u_3D_b_DIVA(   ti,1)
        u_surf_max_BPA    = u_3D_b_BPA(    ti,1)

        ! Compare to target
        IF (ABS( u_surf_max_SIASSA - u_surf_max_check) > u_surf_max_tol) found_errors_SIASSA = .TRUE.
        IF (ABS( u_surf_max_DIVA   - u_surf_max_check) > u_surf_max_tol) found_errors_DIVA   = .TRUE.
        IF (ABS( u_surf_max_BPA    - u_surf_max_check) > u_surf_max_tol) found_errors_BPA    = .TRUE.

      END IF ! IF (ti >= mesh%ti1 .AND. ti < mesh%ti2) THEN

    END IF ! IF (L == 160E3_dp) THEN

    ! If no errors occurred, we are happy
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_SIASSA, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_DIVA  , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, found_errors_BPA   , 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

    IF (.NOT. (found_errors_SIASSA .OR. found_errors_DIVA .OR. found_errors_BPA)) THEN
      IF (par%master) CALL happy('validated SIA/SSA, DIVA, and BPA solvers for diagnostic velocities in the ISMIP-HOM C geometry')
    END IF
    IF (found_errors_SIASSA) THEN
      IF (par%master) CALL warning('found errors in SIA/SSA solver for diagnostic velocities in the ISMIP-HOM C geometry')
    END IF
    IF (found_errors_DIVA) THEN
      IF (par%master) CALL warning('found errors in DIVA solver for diagnostic velocities in the ISMIP-HOM C geometry')
    END IF
    IF (found_errors_BPA) THEN
      IF (par%master) CALL warning('found errors in BPA solver for diagnostic velocities in the ISMIP-HOM C geometry')
    END IF

  ! Write results to NetCDF
  ! =======================

    ! Create a NetCDF output file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // TRIM( filename_ext) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add a zeta dimension for the 3-D ice velocities
    CALL add_zeta_dimension_to_file( filename, ncid, mesh%zeta)

    ! Add all the ice velocities as fields
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_SIASSA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_SIASSA')
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_DIVA'  )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_DIVA'  )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'u_3D_b_BPA'   )
    CALL add_field_mesh_dp_3D_b_notime( filename, ncid, 'v_3D_b_BPA'   )

    ! Write the velocities to the file
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_SIASSA', u_3D_b_SIASSA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_SIASSA', v_3D_b_SIASSA)
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_DIVA'  , u_3D_b_DIVA  )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_DIVA'  , v_3D_b_DIVA  )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'u_3D_b_BPA'   , u_3D_b_BPA   )
    CALL write_to_field_multopt_mesh_dp_3D_b_notime( mesh, filename, ncid, 'v_3D_b_BPA'   , v_3D_b_BPA   )

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Clean up after yourself
    DEALLOCATE( u_3D_b_SIASSA)
    DEALLOCATE( v_3D_b_SIASSA)
    DEALLOCATE( u_3D_b_DIVA)
    DEALLOCATE( v_3D_b_DIVA)
    DEALLOCATE( u_3D_b_BPA)
    DEALLOCATE( v_3D_b_BPA)

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_ISMIP_HOM_C

  SUBROUTINE set_config_for_ISMIP_HOM( L)
    ! Set the config for the ISMIP-HOM experiments
    !
    ! NOTE: general settings only, the reference geometries and sliding law
    !       still need to be set for the specific A-F experiments.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                                            INTENT(IN)    :: L

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'set_config_for_ISMIP_HOM'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Antarctica
    C%lambda_M_ANT                          = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
    C%phi_M_ANT                             = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
    C%beta_stereo_ANT                       = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
    C%xmin_ANT                              = -L                               ! Western  boundary     of the Antarctica domain [m]
    C%xmax_ANT                              =  L                               ! Eastern  boundary     of the Antarctica domain [m]
    C%ymin_ANT                              = -L                               ! Southern boundary     of the Antarctica domain [m]
    C%ymax_ANT                              =  L                               ! Northern boundary     of the Antarctica domain [m]

  ! == Reference geometries (initial, present-day, and GIA equilibrium)
  ! ===================================================================

    ! Some pre-processing stuff for reference ice geometry
    C%refgeo_Hi_min                         = 2.0_dp                           ! Remove ice thinner than this value in the reference ice geometry. Particularly useful for BedMachine Greenland, which somehow covers the entire tundra with half a meter of ice...
    C%remove_Lake_Vostok                    = .TRUE.

    ! == Initial geometry
    ! ===================

    C%choice_refgeo_init_ANT                = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_init == 'idealised'
    C%choice_refgeo_init_idealised          = ''                               ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_init_idealised              = L / 10._dp                       ! Resolution of square grid used for idealised present-day geometry

    ! == Present-day geometry
    ! =======================

    C%choice_refgeo_PD_ANT                  = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_PD == 'idealised'
    C%choice_refgeo_PD_idealised            = ''                               ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_PD_idealised                = L / 10._dp                       ! Resolution of square grid used for idealised present-day geometry

    ! == GIA equilibrium geometry
    ! ===========================

    C%choice_refgeo_GIAeq_ANT               = 'idealised'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", or "read_from_file"
    ! Idealised geometry when choice_refgeo_GIAeq == 'idealised'
    C%choice_refgeo_GIAeq_idealised         = ''                               ! Choice of idealised present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
    C%dx_refgeo_GIAeq_idealised             = L / 10._dp                       ! Resolution of square grid used for idealised present-day geometry

    ! == Parameters for idealised geometries
    ! ======================================

    C%refgeo_idealised_ISMIP_HOM_L          = L                                ! Suggested value: 5E3 - 160E3 m

  ! == Mesh generation
  ! ==================

    ! How to set up the initial mesh
    C%choice_initial_mesh_ANT               = 'calc_from_initial_geometry'     ! Options: 'calc_from_initial_geometry', 'read_from_file'

    ! Resolutions for different parts of the ice sheet
    C%maximum_resolution_uniform            = L / 20._dp                       ! [m]          Maximum resolution for the entire domain
    C%maximum_resolution_grounded_ice       = 100e3_dp                         ! [m]          Maximum resolution for grounded ice
    C%maximum_resolution_floating_ice       = 100e3_dp                         ! [m]          Maximum resolution for floating ice
    C%maximum_resolution_grounding_line     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    C%grounding_line_width                  = 100e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    C%maximum_resolution_calving_front      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    C%calving_front_width                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    C%maximum_resolution_ice_front          = 100e3_dp                         ! [m]          Maximum resolution for the ice front
    C%ice_front_width                       = 100e3_dp                         ! [m]          Width of the band around the ice front that should get this resolution
    C%maximum_resolution_coastline          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    C%coastline_width                       = 100e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

    ! Regions of interest
    C%choice_regions_of_interest            = ''                               ! Regions of interest where other (higher) resolutions apply. Separated by double vertical bars "||", e.g. "PineIsland||Thwaites"
    C%ROI_maximum_resolution_uniform        = 100e3_dp                         ! [m]          Maximum resolution for the entire domain
    C%ROI_maximum_resolution_grounded_ice   = 50e3_dp                          ! [m]          Maximum resolution for grounded ice
    C%ROI_maximum_resolution_floating_ice   = 20e3_dp                          ! [m]          Maximum resolution for floating ice
    C%ROI_maximum_resolution_grounding_line = 5e3_dp                           ! [m]          Maximum resolution for the grounding line
    C%ROI_grounding_line_width              = 5e3_dp                           ! [m]          Width of the band around the grounding line that should get this resolution
    C%ROI_maximum_resolution_calving_front  = 10e3_dp                          ! [m]          Maximum resolution for the calving front
    C%ROI_calving_front_width               = 10e3_dp                          ! [m]          Width of the band around the calving front that should get this resolution
    C%ROI_maximum_resolution_ice_front      = 20e3_dp                          ! [m]          Maximum resolution for the ice front
    C%ROI_ice_front_width                   = 20e3_dp                          ! [m]          Width of the band around the ice front that should get this resolution
    C%ROI_maximum_resolution_coastline      = 50e3_dp                          ! [m]          Maximum resolution for the coastline
    C%ROI_coastline_width                   = 50e3_dp                          ! [m]          Width of the band around the coastline that should get this resolution

    ! Advanced geometry parameters
    C%do_singlecore_mesh_creation           = .TRUE.                           !              Whether or not to use only a single core for mesh generation (for better reproducibility)
    C%alpha_min                             = 0.4363_dp                        ! [radians]    Smallest allowed internal triangle angle
    C%nit_Lloyds_algorithm                  = 3                                ! [-]          Number of iterations of Lloyds algorithm to be applied after refinement
    C%mesh_resolution_tolerance             = 1.25_dp                          ! [-]          Factors the target resolution for trangle-size requirement. 1=strict, use >1 to avoid unnecesarily high resolution

    ! Memory
    C%nC_mem                                = 32                               ! [-]          How many columns of memory should be allocated for connectivity lists

  ! == The scaled vertical coordinate zeta
  ! ======================================

    C%choice_zeta_grid                      = 'regular'                        ! The type of vertical grid to use; can be "regular", "irregular_log", "old_15_layer_zeta"
    C%nz                                    = 12                               ! The number of vertical layers to use
    C%zeta_irregular_log_R                  = 10._dp                           ! Ratio between surface and base layer spacings

  ! == Ice dynamics - velocity
  ! ==========================

    ! General
    C%choice_stress_balance_approximation   = 'SSA'                            ! Choice of stress balance approximation: "none" (= no flow, though geometry can still change due to mass balance), "SIA", "SSA", "SIA/SSA", "DIVA", "BPA"
    C%n_flow                                = 3.0_dp                           ! Exponent in Glen's flow law
    C%m_enh_sheet                           = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
    C%m_enh_shelf                           = 1.0_dp                           ! Ice flow enhancement factor for floating ice
    C%choice_hybrid_SIASSA_scheme           = 'add'                            ! Choice of scheme for combining SIA and SSA velocities in the hybrid approach
    C%do_GL_subgrid_friction                = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
    C%subgrid_friction_exponent             = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta
    C%do_include_SSADIVA_crossterms         = .TRUE.                           ! Whether or not to include the gradients of the effective viscosity (the "cross-terms") in the solution of the SSA/DIVA

    ! Initialisation
    C%choice_initial_velocity_ANT           = 'zero'

    ! Some parameters for numerically solving the stress balance
    C%SIA_maximum_diffusivity               = 1E5_dp                           ! Limit the diffusivity in the SIA to this value
    C%visc_it_norm_dUV_tol                  = 5E-6_dp                          ! Stop criterion for the viscosity iteration: the L2-norm of successive velocity solutions should be smaller than this number
    C%visc_it_nit                           = 500                              ! Maximum number of effective viscosity iterations
    C%visc_it_relax                         = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
    C%epsilon_sq_0                          = 1E-15_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
    C%visc_eff_min                          = 1E4_dp                           ! Minimum value for effective viscosity
    C%beta_max                              = 1E20_dp                          ! Maximum value for basal friction coefficient
    C%vel_max                               = 5000._dp                         ! Velocities are limited to this value
    C%stress_balance_PETSc_rtol             = 1E-5_dp                          ! PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
    C%stress_balance_PETSc_abstol           = 1E-3_dp                          ! PETSc solver - stop criterion, absolute difference

  ! == Ice dynamics - sliding law
  ! =============================

    ! Sliding laws
    C%choice_sliding_law                    = 'idealised'                      ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Budd", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
    C%choice_idealised_sliding_law          = ''                               ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"

    ! Exponents
    C%slid_Weertman_m                       = 3._dp                            ! Exponent in Weertman sliding law
    C%slid_Budd_q_plastic                   = 0.3_dp                           ! Scaling exponent in Budd sliding law
    C%slid_ZI_p                             = 5._dp                            ! Velocity exponent used in the Zoet-Iverson sliding law

    ! Stability
    C%slid_delta_v                          = 1.0E-3_dp                        ! Normalisation parameter to prevent errors when velocity is zero
    C%slid_Budd_u_threshold                 = 100._dp                          ! Threshold velocity in Budd sliding law
    C%slid_ZI_ut                            = 200._dp                          ! (uniform) transition velocity used in the Zoet-Iverson sliding law [m/yr]

  ! == Ice dynamics - boundary conditions
  ! =====================================

    C%BC_u_west                             = 'periodic_ISMIP-HOM'             ! Boundary conditions for the ice velocity field at the domain border
    C%BC_u_east                             = 'periodic_ISMIP-HOM'             ! Allowed choices: "infinite", "zero", "periodic_ISMIP-HOM"
    C%BC_u_south                            = 'periodic_ISMIP-HOM'
    C%BC_u_north                            = 'periodic_ISMIP-HOM'
    C%BC_v_west                             = 'periodic_ISMIP-HOM'
    C%BC_v_east                             = 'periodic_ISMIP-HOM'
    C%BC_v_south                            = 'periodic_ISMIP-HOM'
    C%BC_v_north                            = 'periodic_ISMIP-HOM'
    C%BC_H_west                             = 'zero'                           ! Boundary conditions for ice thickness at the domain boundary
    C%BC_H_east                             = 'zero'                           ! Allowed choices:  "infinite", "zero", "ISMIP_HOM_F"
    C%BC_H_south                            = 'zero'
    C%BC_H_north                            = 'zero'

    ! Rheological model (relating Glen's flow parameter to ice temperature)
    C%choice_ice_rheology                   = 'uniform'                        ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
    C%uniform_flow_factor                   = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")

    ! Add routine to path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_config_for_ISMIP_HOM

  ! ===== Ice thickness evolution in the Halfar dome geometry =====
  ! ===============================================================

  SUBROUTINE test_thickness_evolution_Halfar_dome
    ! Generate a simple Halfar dome geometry and calculate ice thickness evolution
    ! for a few thousand years to see if it matches the analytical solution.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                      :: routine_name = 'test_thickness_evolution_Halfar_dome'
    TYPE(type_mesh)                                                    :: mesh
    TYPE(type_ice_model)                                               :: ice
    TYPE(type_regional_scalars)                                        :: scalars
    TYPE(type_reference_geometry)                                      :: refgeo_init, refgeo_PD, refgeo_GIAeq
    CHARACTER(LEN=256)                                                 :: region_name, mesh_name
    REAL(dp), PARAMETER                                                :: tstart    =     0._dp   ! [yr] Start time of simulation
    REAL(dp), PARAMETER                                                :: tstop     = 5000._dp    ! [yr] End   time of simulation
    REAL(dp)                                                           :: dt                      ! [yr] Time step
    REAL(dp), PARAMETER                                                :: dt_output =  100._dp    ! [yr] Time step
    REAL(dp)                                                           :: time, time_next_output
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                            :: SMB, BMB
    INTEGER                                                            :: vi
    REAL(dp)                                                           :: ice_volume_init, ice_volume
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                            :: Hi_an, dHi_sq
    REAL(dp)                                                           :: MSE_Hi, RMSE_Hi
    CHARACTER(LEN=256)                                                 :: filename
    INTEGER                                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to screen that we're doing this experiment
    IF (par%master) WRITE(0,*) '   Calculating ice thickness evolution for the Halfar dome...'

    ! Model domain
    region_name = 'ANT'

    ! Set model configuration for this experiment
    CALL set_config_for_Halfar_dome

    ! Specific for this experiment
    C%maximum_resolution_uniform            = 40e3_dp                          ! [m]          Maximum resolution for the entire domain
    C%maximum_resolution_grounded_ice       = 100e3_dp                         ! [m]          Maximum resolution for grounded ice
    C%maximum_resolution_floating_ice       = 100e3_dp                         ! [m]          Maximum resolution for floating ice
    C%maximum_resolution_grounding_line     = 100e3_dp                         ! [m]          Maximum resolution for the grounding line
    C%grounding_line_width                  = 100e3_dp                         ! [m]          Width of the band around the grounding line that should get this resolution
    C%maximum_resolution_calving_front      = 100e3_dp                         ! [m]          Maximum resolution for the calving front
    C%calving_front_width                   = 100e3_dp                         ! [m]          Width of the band around the calving front that should get this resolution
    C%maximum_resolution_ice_front          = 100e3_dp                         ! [m]          Maximum resolution for the ice front
    C%ice_front_width                       = 100e3_dp                         ! [m]          Width of the band around the ice front that should get this resolution
    C%maximum_resolution_coastline          = 100e3_dp                         ! [m]          Maximum resolution for the coastline
    C%coastline_width                       = 100e3_dp                         ! [m]          Width of the band around the coastline that should get this resolution

  ! Initialise the model
  ! ====================

    ! Initialise all the reference geometries on their raw input grids
    CALL initialise_reference_geometries_raw( region_name, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Create mesh from gridded initial geometry data
    mesh_name = 'Halfar_mesh'
    CALL create_mesh_from_gridded_geometry( region_name, mesh_name, &
      refgeo_init%grid_raw, &
      refgeo_init%Hi_grid_raw, &
      refgeo_init%Hb_grid_raw, &
      refgeo_init%Hs_grid_raw, &
      refgeo_init%SL_grid_raw, &
      C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
      mesh)

    ! Write the mesh creation success message to the terminal
    CALL write_mesh_success( mesh)

    ! Remap reference geometries from their raw input grids to the model mesh
    CALL initialise_reference_geometries_on_model_mesh( region_name, mesh, refgeo_init, refgeo_PD, refgeo_GIAeq)

    ! Initialise the ice model
    C%choice_stress_balance_approximation = 'SIA'
    CALL initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, scalars, region_name)
    ice%A_flow_3D = C%uniform_flow_factor
    ALLOCATE( ice%beta_b( mesh%vi1:mesh%vi2))

    ! Calculate initial ice volume
    CALL integrate_over_domain( mesh, ice%Hi, ice_volume_init)

    ! Allocate memory for (unused) SMB and BMB
    ALLOCATE( SMB(   mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( BMB(   mesh%vi1:mesh%vi2), source = 0._dp)
    ALLOCATE( Hi_an( mesh%vi1:mesh%vi2), source = 0._dp)

  ! Create a nice output file
  ! =========================

    ! Create a NetCDF output file
    filename = TRIM( C%output_dir) // TRIM( routine_name) // '_output.nc'
    CALL create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( filename, ncid)

    ! Add data fields
    CALL add_field_mesh_dp_2D(   filename, ncid, 'Hi_an'  , long_name = 'Ice thickness (analytical solution)', units = 'm')
    CALL add_field_mesh_dp_2D(   filename, ncid, 'Hi'     , long_name = 'Ice thickness', units = 'm')
    CALL add_field_mesh_dp_2D_b( filename, ncid, 'u_vav_b', long_name = 'Vertically averaged velocity in x-direction', units = 'm/yr')
    CALL add_field_mesh_dp_2D_b( filename, ncid, 'v_vav_b', long_name = 'Vertically averaged velocity in y-direction', units = 'm/yr')

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! Integrate the model through time
  ! ================================

    time             = tstart
    time_next_output = tstart
    dt               = 0.5_dp
    DO WHILE (time < tstop)

      ! Solve the stress balance to calculate ice velocities
      CALL solve_stress_balance( mesh, ice)

      ! After 100 years we can go for a larger time step
      IF (time >= 100._dp) dt = 1._dp
      ! After 1000 years we can go for an even larger time step
      IF (time >= 1000._dp) dt = 2.5_dp
      ! After 2000 years we can go for an even even larger time step
      IF (time >= 2000._dp) dt = 5._dp

      ! Calculate ice thickness in next time step
      CALL calc_Hi_tplusdt( mesh, ice%Hi, ice%u_vav_b, ice%v_vav_b, SMB, BMB, dt, ice%dHi_dt, ice%Hi_tplusdt)

      ! Write to output
      IF (time >= time_next_output - dt / 2._dp) THEN

        ! Update timer
        time_next_output = time_next_output + dt_output

        ! Open NetCDF file
        CALL open_existing_netcdf_file_for_writing(     filename, ncid)

        ! Write current time to file
        CALL write_time_to_file(                        filename, ncid, time)

        ! Calculate analytical solution
        DO vi = mesh%vi1, mesh%vi2
          ! Calculate analytical solution
          CALL Halfar_dome( C%uniform_flow_factor, C%n_flow, C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
            mesh%V( vi,1), mesh%V( vi,2), time, Hi_an( vi))
        END DO ! DO vi = mesh%vi1, mesh%vi2

        ! Write model data to file
        CALL write_to_field_multopt_mesh_dp_2D(   mesh, filename, ncid, 'Hi_an'  , Hi_an      )
        CALL write_to_field_multopt_mesh_dp_2D(   mesh, filename, ncid, 'Hi'     , ice%Hi     )
        CALL write_to_field_multopt_mesh_dp_2D_b( mesh, filename, ncid, 'u_vav_b', ice%u_vav_b)
        CALL write_to_field_multopt_mesh_dp_2D_b( mesh, filename, ncid, 'v_vav_b', ice%v_vav_b)

        ! Close the file
        CALL close_netcdf_file( ncid)

      END IF ! IF (time >= time_next_output) THEN

      ! Advance to next time step
      ice%Hi = ice%Hi_tplusdt
      DO vi = mesh%vi1, mesh%vi2
        ice%Hs( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))
      END DO
      time   = time + dt

    END DO ! DO WHILE (time < tstop)

    ! Open NetCDF file
    CALL open_existing_netcdf_file_for_writing(     filename, ncid)

    ! Write current time to file
    CALL write_time_to_file(                        filename, ncid, time)

    ! Calculate analytical solution
    DO vi = mesh%vi1, mesh%vi2
      ! Calculate analytical solution
      CALL Halfar_dome( C%uniform_flow_factor, C%n_flow, C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
        mesh%V( vi,1), mesh%V( vi,2), time, Hi_an( vi))
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Write model data to file
    CALL write_to_field_multopt_mesh_dp_2D(   mesh, filename, ncid, 'Hi_an'  , Hi_an      )
    CALL write_to_field_multopt_mesh_dp_2D(   mesh, filename, ncid, 'Hi'     , ice%Hi     )
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, filename, ncid, 'u_vav_b', ice%u_vav_b)
    CALL write_to_field_multopt_mesh_dp_2D_b( mesh, filename, ncid, 'v_vav_b', ice%v_vav_b)

    ! Close the file
    CALL close_netcdf_file( ncid)

  ! == Validation
  ! =============

    ! Calculate integrated ice volume and compare to initial volume (should be identical)
    CALL integrate_over_domain( mesh, ice%Hi, ice_volume)
    IF (ABS( 1._dp - ice_volume / ice_volume_init) > 1E-5_dp) THEN
      IF (par%master) CALL warning('Halfar dome experiment violates conservation of mass')
    END IF

    ! Calculate (squared) error in modelled ice thickness
    ALLOCATE( dHi_sq( mesh%vi1:mesh%vi2))
    dHi_sq = (ice%Hi - Hi_an)**2

    ! Calculate root-mean-square-error in modelled ice thickness
    CALL average_over_domain( mesh, dHi_sq, MSE_Hi)
    RMSE_Hi = SQRT( MSE_Hi)
    IF (RMSE_Hi > 0.2E2_dp) THEN
      IF (par%master) CALL warning('Halfar dome experiment deviates further than expected from analytical solution')
    END IF

    ! Clean up after yourself
    DEALLOCATE( SMB   )
    DEALLOCATE( BMB   )
    DEALLOCATE( Hi_an )
    DEALLOCATE( dHi_sq)

    ! Finalise routine
    CALL finalise_routine( routine_name)

  END SUBROUTINE test_thickness_evolution_Halfar_dome

END MODULE unit_tests_ice
