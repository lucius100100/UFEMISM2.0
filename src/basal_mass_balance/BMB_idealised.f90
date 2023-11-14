MODULE BMB_idealised

  ! Idealised BMB models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  USE mesh_refinement                                        , ONLY: calc_polygon_Pine_Island_Glacier, calc_polygon_Thwaites_Glacier, &
                                                                     calc_polygon_Tijn_test_ISMIP_HOM_A
  USE math_utilities                                         , ONLY: is_in_polygon, is_in_polygons
  USE reference_geometry_types                               , ONLY: type_reference_geometry

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_idealised( mesh, ice, BMB, time, region_name, refgeo)
    ! Calculate the basal mass balance
    !
    ! Use an idealised BMB scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen idealised BMB model
    SELECT CASE (C%choice_BMB_model_idealised)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_idealised "' // TRIM( C%choice_BMB_model_idealised) // '"')
      CASE ('MISMIP+')
        CALL run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
      CASE ('MISMIPplus')
        CALL run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
      CASE ('Tijn_ROI_Thwaites')
        CALL run_BMB_model_idealised_Tijn_ROI_Thwaites( mesh, ice, BMB, time, region_name, refgeo)
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised

  SUBROUTINE run_BMB_model_idealised_MISMIPplus( mesh, ice, BMB, time)
    ! The schematic basal melt used in the MISMIPplus experiments

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_idealised_MISMIPplus'
    INTEGER                                            :: vi
    REAL(dp)                                           :: zd, cavity_thickness

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    BMB%BMB = 0._dp

    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_floating_ice( vi)) THEN

        zd = ice%Hs( vi) - ice%Hi( vi)
        cavity_thickness = MAX( 0._dp, zd - ice%Hb( vi))

        ! Cornford et al. (2020), Eq. 7
        BMB%BMB( vi) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)

      END IF ! IF (ice%mask_floating_ice( vi)) THEN
    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised_MISMIPplus

  SUBROUTINE run_BMB_model_idealised_Tijn_ROI_Thwaites( mesh, ice, BMB, time, region_name, refgeo)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    TYPE(type_reference_geometry),          INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_idealised_Tijn_ROI_Thwaites'
    INTEGER                                            :: vi
    CHARACTER(LEN=256)                                             :: all_names_ROI, name_ROI
    INTEGER                                                        :: i
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                        :: poly_ROI
    INTEGER                                                        :: ti
    REAL(dp), DIMENSION(2)                                         :: p

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Go over all listed regions of interest
    all_names_ROI = C%choice_regions_of_interest

    DO WHILE (.TRUE.)

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      IF (i == 0) THEN
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      ELSE
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      END IF


      SELECT CASE (region_name)
        CASE DEFAULT
          CALL crash('unknown region name "' // region_name // '"!')
        CASE ('NAM')
          ! North america

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('EAS')
          ! Eurasia

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('GRL')
          ! Greenland

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              ! Don't need to do anything
              EXIT
            CASE ('Thwaites')
              ! Don't need to do anything
              EXIT
          END SELECT

        CASE ('ANT')

          SELECT CASE (name_ROI)
            CASE DEFAULT
              CALL crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
            CASE ('')
              ! Don't need to do anything
              EXIT
            CASE ('PineIsland')
              CALL calc_polygon_Pine_Island_Glacier( poly_ROI)
            CASE ('Thwaites')
              CALL calc_polygon_Thwaites_Glacier( poly_ROI)
            CASE ('Tijn_test_ISMIP_HOM_A')
              CALL calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
          END SELECT

      END SELECT

      ! Find all vertices that lie within this region of interest
      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_floating_ice( vi) .OR. ice%mask_icefree_ocean( vi) .OR. ice%mask_gl_gr( vi)) THEN
          p = mesh%V( vi,:)
          IF (is_in_polygon( poly_ROI, p)) THEN
            BMB%BMB( vi) = -50._dp
          ELSE
            BMB%BMB( vi) = (ice%Hi( vi) - refgeo%Hi( vi)) / -5._dp
          END IF
        ELSE
          BMB%BMB( vi) = 0._dp
        END IF
      END DO

      ! Clean up after yourself
      DEALLOCATE( poly_ROI)

      ! If no names are left, we are finished
      IF (all_names_ROI == '') EXIT

    END DO ! DO WHILE (.TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_idealised_Tijn_ROI_Thwaites

  SUBROUTINE initialise_BMB_model_idealised( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Use an idealised BMB scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '   Initialising idealised BMB model "' // &
      colour_string( TRIM( C%choice_BMB_model_idealised),'light blue') // '"...'

    ! Initialise the chosen idealised BMB model
    SELECT CASE (C%choice_BMB_model_idealised)
      CASE DEFAULT
        CALL crash('unknown choice_BMB_model_idealised "' // TRIM( C%choice_BMB_model_idealised) // '"')
      CASE ('MISMIP+')
        ! No need to do anything
      CASE ('Tijn_ROI_Thwaites')
        ! No need to do anything
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_idealised

END MODULE BMB_idealised
