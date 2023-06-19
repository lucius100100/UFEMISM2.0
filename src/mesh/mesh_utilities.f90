MODULE mesh_utilities

  ! Generally useful functions used in mesh creation and updating.

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE reallocate_mod
  USE mesh_types                                             , ONLY: type_mesh
  USE math_utilities                                         , ONLY: geometric_center, is_in_triangle, lies_on_line_segment, circumcenter, &
                                                                     line_from_points, line_line_intersection, encroaches_upon, crop_line_to_domain

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

! == Finding the vertices of a vertex' Voronoi cell
  SUBROUTINE calc_Voronoi_cell( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)
    ! Find the points spanning the Voronoi cell of vertex vi of the mesh.
    ! Sorted counted-clockwise, with no double entries.
    !
    ! Point Vor( i) corresponds to the circumcentre of triangle Vor_ti( i).
    !
    ! The line connecting Vor( i) and Vor( i+1) is shared with the Voronoi cell
    ! of vertex Vor_vi( i).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(IN)    :: dx
    REAL(dp), DIMENSION(mesh%nC_mem,2),  INTENT(OUT)   :: Vor
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_vi
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_ti
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:

    ! Vertices lying on the domain border or corners are best treated separately
    IF (mesh%VBI( vi) == 0) THEN
      ! Free vertex
      CALL calc_Voronoi_cell_free(   mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)
    ELSE
      ! Border vertex
      CALL calc_Voronoi_cell_border( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)
    END IF

  END SUBROUTINE calc_Voronoi_cell

  SUBROUTINE calc_Voronoi_cell_free( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)
    ! Find the points spanning the Voronoi cell of vertex vi of the mesh%
    ! Sorted counted-clockwise, with no double entries.
    !
    ! Point Vor( i) corresponds to the circumcentre of triangle Vor_ti( i).
    !
    ! The line connecting Vor( i) and Vor( i+1) is shared with the Voronoi cell
    ! of vertex Vor_vi( i).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(IN)    :: dx
    REAL(dp), DIMENSION(mesh%nC_mem,2),  INTENT(OUT)   :: Vor
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_vi
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_ti
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: iti,ti,vj,n1,n2,vori,vori_prev,vori_next,vk
    REAL(dp), DIMENSION(2)                             :: cc,cc_prev,cc_next,cc_prev0,cc_next0,cc1,cc2
    LOGICAL                                            :: is_valid_line

    ! Initialise
    Vor    = 0._dp
    Vor_vi = 0
    Vor_ti = 0
    nVor   = 0

    ! List the circumcentres of all triangles surrounding vi
    ! as points spanning the Voronoi cell

    DO iti = 1, mesh%niTri( vi)

      ti  = mesh%iTri( vi,iti)

      ! Find vertex vj such that triangle ti lies to the right of the line vi-vj
      vj = 0
      DO n1 = 1, 3
        n2 = n1 + 1
        IF (n2 == 4) n2 = 1
        IF (mesh%Tri( ti,n2) == vi) THEN
          vj = mesh%Tri( ti,n1)
          EXIT
        END IF
      END DO
      ! Safety
      IF (vj == 0) CALL crash('calc_Voronoi_cell_free - couldnt find vertex vj in triangle ti!')

      ! List the new Voronoi cell-spanning point
      nVor = nVor + 1
      Vor(     nVor,:) = mesh%Tricc( ti,:)
      Vor_vi(  nVor  ) = vj
      Vor_ti(  nVor  ) = ti

    END DO

    ! Correct any circumcentres that lie outside of the mesh domain
    vori = 1
    DO WHILE (vori < nVor)

      vori_prev = vori - 1
      IF (vori_prev == 0   ) vori_prev = nVor
      vori_next = vori + 1
      IF (vori_next >  nVor) vori_next = 1

      cc      = Vor( vori     ,:)
      cc_prev = Vor( vori_prev,:)
      cc_next = Vor( vori_next,:)

      IF (cc( 1) < mesh%xmin .OR. cc( 1) > mesh%xmax .OR. cc( 2) < mesh%ymin .OR. cc( 2) > mesh%ymax) THEN
        ! This triangle circumcentre is outside of the mesh domain

        ! Find the two points on the domain border that should actually be part
        ! of the Voronoi cell boundary
        CALL crop_line_to_domain( cc_prev, cc, mesh%xmin - dx, mesh%xmax + dx, mesh%ymin - dx, mesh%ymax + dx, mesh%tol_dist, cc_prev0, cc1, is_valid_line)
        CALL crop_line_to_domain( cc_next, cc, mesh%xmin - dx, mesh%xmax + dx, mesh%ymin - dx, mesh%ymax + dx, mesh%tol_dist, cc_next0, cc2, is_valid_line)

        ! The previous and next vertex, and the edge connecting them
        vj  = Vor_vi(  vori_prev)
        vk  = Vor_vi(  vori)

        ! Split the invalid Voronoi cell boundary point
        nVor = nVor + 1
        Vor(     1:nVor,1) = [Vor(     1:vori-1,1), cc1( 1), cc2( 1), Vor(     vori+1:nVor-1,1)]
        Vor(     1:nVor,2) = [Vor(     1:vori-1,2), cc1( 2), cc2( 2), Vor(     vori+1:nVor-1,2)]
        Vor_vi(  1:nVor  ) = [Vor_vi(  1:vori-1  ), 0      , vk     , Vor_vi(  vori+1:nVor-1  )]
        Vor_ti(  1:nVor  ) = [Vor_ti(  1:vori-1  ), 0      , 0      , Vor_ti(  vori+1:nVor-1  )]

        ! Move to the next to-be-checked point
        vori = vori + 2

      ELSE
        ! This triangle circumcentre is inside the mesh domain; move the next one

        vori = vori + 1

      END IF

    END DO

  END SUBROUTINE calc_Voronoi_cell_free

  SUBROUTINE calc_Voronoi_cell_border( mesh, vi, dx, Vor, Vor_vi, Vor_ti, nVor)
    ! Find the points spanning the Voronoi cell of vertex vi of the mesh.
    ! Sorted counted-clockwise, with no double entries.
    !
    ! Point Vor( i) corresponds to the circumcentre of triangle Vor_ti( i).
    !
    ! The line connecting Vor( i) and Vor( i+1) is shared with the Voronoi cell
    ! of vertex Vor_vi( i).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,                             INTENT(IN)    :: vi
    REAL(dp),                            INTENT(IN)    :: dx
    REAL(dp), DIMENSION(mesh%nC_mem,2),  INTENT(OUT)   :: Vor
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_vi
    INTEGER,  DIMENSION(mesh%nC_mem  ),  INTENT(OUT)   :: Vor_ti
    INTEGER,                             INTENT(OUT)   :: nVor

    ! Local variables:
    INTEGER                                            :: ci,vj,ti,iti,tj,n,n2
    REAL(dp), DIMENSION(2)                             :: p,q,pp,qq
    LOGICAL                                            :: is_valid_line

    ! Initialise
    Vor    = 0._dp
    Vor_vi = 0
    Vor_ti = 0
    nVor   = 0

    ! For each neighbouring vertex vj, list the circumcentre of the triangle
    ! lying right of the line vi-vj. Except, we skip the first neighbouring
    ! vertex, as DO that one the area to the right of the line vi-vj is
    ! outside of the mesh domain

    DO ci = 2, mesh%nC( vi)

      vj  = mesh%C( vi,ci)

      ! Find the triangle ti lying right of the line vi-vj
      ti = 0
      DO iti = 1, mesh%niTri( vi)
        tj = mesh%iTri( vi,iti)
        DO n = 1, 3
          n2 = n + 1
          IF (n2 == 4) n2 = 1
          IF (mesh%Tri( tj,n) == vj .AND. mesh%Tri( tj,n2) == vi) THEN
            ti = tj
            EXIT
          END IF
        END DO
        IF (ti > 0) EXIT
      END DO

      ! Safety
      IF (ti == 0) CALL crash('couldnt find triangle right of the line vi-vj!')

      ! List the circumcentre of this triangle
      nVor = nVor + 1
      Vor(     nVor,:) = mesh%Tricc( ti,:)
      Vor_vi(  nVor  ) = vj
      Vor_ti(  nVor  ) = ti

    END DO ! DO ci = 2, mesh%nC( vi)

    ! Fix the first point

    IF (Vor( 1,1) < mesh%xmin - dx .OR. Vor( 1,1) > mesh%xmax + dx .OR. &
        Vor( 1,2) < mesh%ymin - dx .OR. Vor( 1,2) > mesh%ymax + dx) THEN
      ! The first circumcentre lies outside the mesh domain; crop it

      p = Vor( 1,:)
      q = Vor( 2,:)

      CALL crop_line_to_domain( p, q, mesh%xmin - dx, mesh%xmax + dx, mesh%ymin - dx, mesh%ymax + dx, mesh%tol_dist, pp, qq, is_valid_line)

      Vor(    1,:) = pp
      Vor_ti( 1  ) = 0

    ELSE
      ! The first circumcentre lies inside the mesh domain; add its projection
      ! on the domain boundary as an additional point

      nVor = nVor + 1
      Vor(     2:nVor,:) = Vor(     1:nVor-1,:)
      Vor_vi(  2:nVor  ) = Vor_vi(  1:nVor-1  )
      Vor_ti(  2:nVor  ) = Vor_ti(  1:nVor-1  )

      Vor_vi( 1) = mesh%C(    vi,1)
      Vor_ti( 1) = mesh%iTri( vi,1)

      IF     (mesh%VBI( vi) == 1 .OR. mesh%VBI( vi) == 2) THEN
        ! vi is on the northern border, or on the northeast corner; its first neighbour lies on the northern border
        Vor(    1,:) = [Vor( 2,1), mesh%ymax + dx]
      ELSEIF (mesh%VBI( vi) == 3 .OR. mesh%VBI( vi) == 4) THEN
        ! vi is on the eastern border, or on the southeast corner; its first neighbour lies on the eastern border
        Vor(    1,:) = [mesh%xmax + dx, Vor( 2,2)]
      ELSEIF (mesh%VBI( vi) == 5 .OR. mesh%VBI( vi) == 6) THEN
        ! vi is on the southern border, or on the southwest corner; its first neighbour lies on the southern border
        Vor(    1,:) = [Vor( 2,1), mesh%ymin - dx]
      ELSEIF (mesh%VBI( vi) == 7 .OR. mesh%VBI( vi) == 8) THEN
        ! vi is on the western border, or on the northwest corner; its first neighbour lies on the western border
        Vor(    1,:) = [mesh%xmin - dx, Vor( 2,2)]
      END IF

    END IF

    ! Fix the last point

    IF (Vor( nVor,1) < mesh%xmin - dx .OR. Vor( nVor,1) > mesh%xmax + dx .OR. &
        Vor( nVor,2) < mesh%ymin - dx .OR. Vor( nVor,2) > mesh%ymax + dx) THEN
      ! The last circumcentre lies outside the mesh domain; crop it

      p = Vor( nVor-1,:)
      q = Vor( nVor  ,:)

      CALL crop_line_to_domain( p, q, mesh%xmin - dx, mesh%xmax + dx, mesh%ymin - dx, mesh%ymax + dx, mesh%tol_dist, pp, qq, is_valid_line)

      Vor( nVor,:) = qq

    ELSE
      ! The last circumcentre lies inside the mes domain; add its projection
      ! on the domain boundary as an additional point

      nVor = nVor + 1

      Vor_vi( nVor) = mesh%C(    vi, mesh%nC( vi))
      Vor_ti( nVor) = mesh%iTri( vi, mesh%nC( vi))

      IF     (mesh%VBI( vi) == 2 .OR. mesh%VBI( vi) == 3) THEN
        ! vi is on the eastern border, or on the northeast corner; its last neighbour lies on the eastern border
        Vor(    nVor,:) = [mesh%xmax + dx, Vor( nVor-1,2)]
      ELSEIF (mesh%VBI( vi) == 4 .OR. mesh%VBI( vi) == 5) THEN
        ! vi is on the southern border, or on the southeast corner; its last neighbour lies on the southern border
        Vor(    nVor,:) = [Vor( nVor-1,1), mesh%ymin - dx]
      ELSEIF (mesh%VBI( vi) == 6 .OR. mesh%VBI( vi) == 7) THEN
        ! vi is on the western border, or on the southwest corner; its last neighbour lies on the western border
        Vor(    nVor,:) = [mesh%xmin - dx, Vor( nVor-1,2)]
      ELSEIF (mesh%VBI( vi) == 8 .OR. mesh%VBI( vi) == 1) THEN
        ! vi is on the northern border, or on the northwest corner; its last neighbour lies on the northern border
        Vor(    nVor,:) = [Vor( nVor-1,1), mesh%ymax + dx]
      END IF

    END IF

    ! In the case of the four corners, add the corners themselves as well
    IF     (mesh%VBI( vi) == 2) THEN
      ! Northeast corner
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax + dx, mesh%ymax + dx]
    ELSEIF (mesh%VBI( vi) == 4) THEN
      ! Southeast corner
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmax + dx, mesh%ymin - dx]
    ELSEIF (mesh%VBI( vi) == 6) THEN
      ! Southwest corner
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin - dx, mesh%ymin - dx]
    ELSEIF (mesh%VBI( vi) == 8) THEN
      ! Northwest corner
      nVor = nVor + 1
      Vor( nVor,:) = [mesh%xmin - dx, mesh%ymax + dx]
    END IF

  END SUBROUTINE calc_Voronoi_cell_border

  SUBROUTINE find_shared_Voronoi_boundary( mesh, ei, cc1, cc2)
    ! Return the endpoints of the shared Voronoi cell boundary represented by edge ei

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: ei
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: cc1, cc2

    ! Local variables
    REAL(dp)                                      :: dx
    INTEGER                                       :: til,tir,ti
    REAL(dp), DIMENSION(2)                        :: cc1_raw, cc2_raw
    LOGICAL                                       :: is_valid_line

    dx = ((mesh%xmax - mesh%xmin) + (mesh%ymax - mesh%ymin)) / 100._dp

    til = mesh%ETri( ei,1)
    tir = mesh%ETri( ei,2)

    IF (mesh%EBI( ei) == 0) THEN
      ! This is an internal edge, with two adjacent triangles

      cc1_raw = mesh%Tricc( til,:)
      cc2_raw = mesh%Tricc( tir,:)

      ! Crop the line to the mesh domain
      CALL crop_line_to_domain( cc1_raw, cc2_raw, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, mesh%tol_dist, cc1, cc2, is_valid_line)

      ! Safety
      IF (.NOT. is_valid_line) CALL crash('find_shared_Voronoi_boundary - couldnt crop shared Voronoi boundary!')

    ELSE ! IF (mesh%EBI( ei) == 0) THEN
      ! This edge lies on the domain border, so the two vertices share only a single triangle

      IF     (til > 0 .AND. tir == 0) THEN
        ti = til
      ELSEIF (tir > 0 .AND. til == 0) THEN
        ti = tir
      ELSE
        CALL crash('find_shared_Voronoi_boundary - something is seriously wrong with EBI!')
      END IF

      cc1 = mesh%Tricc( ti,:)

      ! Get the projection of this circumcentre on the domain border
      IF     (mesh%EBI( ei) == 1) THEN
        ! Northern border
        cc2 = [cc1( 1), mesh%ymax + dx]
      ELSEIF (mesh%EBI( ei) == 3) THEN
        ! Eastern border
        cc2 = [mesh%xmax + dx, cc1( 2)]
      ELSEIF (mesh%EBI( ei) == 5) THEN
        ! Southern border
        cc2 = [cc1( 1), mesh%ymin - dx]
      ELSEIF (mesh%EBI( ei) == 7) THEN
        ! Western border
        cc2 = [mesh%xmin - dx, cc1( 2)]
      END IF

      ! The circumcentre itself may not lie outside of the domain
      cc1( 1) = MIN( mesh%xmax, MAX( mesh%xmin, cc1( 1) ))
      cc1( 2) = MIN( mesh%ymax, MAX( mesh%ymin, cc1( 2) ))

    END IF ! IF (mesh%EBI( ei) == 0) THEN

  END SUBROUTINE find_shared_Voronoi_boundary

! == Some basic geometrical operations

  SUBROUTINE list_border_vertices_all( mesh, nvi_border, lvi_border)
    ! List all border vertices of the mesh, ordered counter-clockwise,
    ! starting from vertex 1 in the southwest corner.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it

    ! Initialise
    lvi_border = 0
    nvi_border = 0

    ! Trace the border
    vi = 1
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! If we've reached the southwest corner again, stop
      IF (vi == 1) EXIT

    END DO

  END SUBROUTINE list_border_vertices_all

  SUBROUTINE list_border_vertices_south( mesh, nvi_border, lvi_border)
    ! List all southern border vertices of the mesh, ordered east to west.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! South
    vi_start = 2
    vi_stop  = 1

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi, mesh%nC( vi))

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_south

  SUBROUTINE list_border_vertices_east( mesh, nvi_border, lvi_border)
    ! List all eastern border vertices of the mesh, ordered south to north.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! East
    vi_start = 2
    vi_stop  = 3

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_east

  SUBROUTINE list_border_vertices_north( mesh, nvi_border, lvi_border)
    ! List all northern border vertices of the mesh, ordered east to west.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! North
    vi_start = 3
    vi_stop  = 4

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi,1)

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_north

  SUBROUTINE list_border_vertices_west( mesh, nvi_border, lvi_border)
    ! List all western border vertices of the mesh, ordered south to north.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    INTEGER,                      INTENT(OUT)     :: nvi_border
    INTEGER,  DIMENSION(mesh%nV), INTENT(OUT)     :: lvi_border

    ! Local variables:
    INTEGER                                       :: vi, it, vi_start, vi_stop

    ! West
    vi_start = 1
    vi_stop  = 4

    ! Initialise
    lvi_border = 0
    lvi_border( 1) = vi_start
    nvi_border = 1

    ! Start in the northwest corner
    vi = vi_start

    ! Trace the border
    it = 0
    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nV) CALL crash('border trace got stuck!')

      ! Move to next vertex along the border
      vi = mesh%C( vi, mesh%nC( vi))

      ! List current vertex
      nvi_border = nvi_border + 1
      lvi_border( nvi_border) = vi

      ! If we've reached the southwest corner, stop
      IF (vi == vi_stop) EXIT

    END DO

  END SUBROUTINE list_border_vertices_west

  SUBROUTINE encroaches_upon_any( mesh, p, isso, vi_encroached, vj_encroached)
    ! Check if p encroaches upon a border edge.
    ! If so, return the indices [vi,vj] of the vertices spanning that edge.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p
    LOGICAL,                    INTENT(OUT)       :: isso
    INTEGER,                    INTENT(OUT)       :: vi_encroached, vj_encroached

    ! Local variables:
    INTEGER                                       :: nvi_border
    INTEGER,  DIMENSION(mesh%nV)                  :: lvi_border
    INTEGER                                       :: i, j, vi, vj
    REAL(dp), DIMENSION(2)                        :: pa, pb

    isso          = .FALSE.
    vi_encroached = 0
    vj_encroached = 0

    ! List border vertices
    CALL list_border_vertices_all( mesh, nvi_border, lvi_border)

    ! Check if p encroaches on any of the border edges
    DO i = 1, nvi_border

      j = i + 1
      IF (i == nvi_border) j = 1

      vi = lvi_border( i)
      vj = lvi_border( j)

      pa = mesh%V( vi,:)
      pb = mesh%V( vj,:)

      IF (encroaches_upon( pa, pb, p, mesh%tol_dist)) THEN
        isso          = .TRUE.
        vi_encroached = vi
        vj_encroached = vj
        EXIT
      END IF

    END DO ! DO i = 1, nvi_border

  END SUBROUTINE encroaches_upon_any

  PURE FUNCTION is_border_edge( mesh, vi, vj) RESULT(isso)
    ! Determine whether or not the edge between two vertices is a border edge

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: vi, vj

    ! Local variables:
    LOGICAL                                       :: isso

    IF (mesh%VBI( vi) == 0 .OR. mesh%VBI( vj) == 0) THEN
      isso = .FALSE.
      RETURN
    END IF

    isso = .FALSE.

    IF (mesh%VBI( vi) == 1) THEN
      IF (mesh%VBI( vj) == 8 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 2) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 2) THEN
      IF (mesh%VBI( vj) == 8 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 4) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 3) THEN
      IF (mesh%VBI( vj) == 2 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 4) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 4) THEN
      IF (mesh%VBI( vj) == 2 .OR. &
          mesh%VBI( vj) == 3 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 6) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 5) THEN
      IF (mesh%VBI( vj) == 4 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 6) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 6) THEN
      IF (mesh%VBI( vj) == 4 .OR. &
          mesh%VBI( vj) == 5 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 8) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 7) THEN
      IF (mesh%VBI( vj) == 6 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 8) THEN
        isso = .TRUE.
      END IF
    ELSEIF (mesh%VBI( vi) == 8) THEN
      IF (mesh%VBI( vj) == 6 .OR. &
          mesh%VBI( vj) == 7 .OR. &
          mesh%VBI( vj) == 1 .OR. &
          mesh%VBI( vj) == 2) THEN
        isso = .TRUE.
      END IF
    END IF

  END FUNCTION is_border_edge

  PURE FUNCTION is_walltowall( mesh, ti) RESULT(isso)
   ! Determine whether or not a triangle is "wall to wall"
   ! (i.e. contains vertices lying on opposite domain borders)

    IMPLICIT NONE

    ! In/output variables:
   TYPE(type_mesh),          INTENT(IN)          :: mesh
   INTEGER,                  INTENT(IN)          :: ti
   LOGICAL                                       :: isso

    ! Local variables:
   LOGICAL                                       :: has_north, has_south, has_east, has_west
   INTEGER                                       :: n, vi

   has_north = .FALSE.
   has_south = .FALSE.
   has_east  = .FALSE.
   has_west  = .FALSE.

   DO n = 1, 3
     vi = mesh%Tri( ti,n)
     IF (mesh%VBI( vi) == 1) has_north = .TRUE.
     IF (mesh%VBI( vi) == 3) has_east  = .TRUE.
     IF (mesh%VBI( vi) == 5) has_south = .TRUE.
     IF (mesh%VBI( vi) == 7) has_west  = .TRUE.
   END DO

   isso = .FALSE.
   IF (has_north .AND. has_south) isso = .TRUE.
   IF (has_west  .AND. has_east ) isso = .TRUE.

  END FUNCTION is_walltowall

  SUBROUTINE update_triangle_circumcenter( mesh, ti)
    ! Calculate the circumcenter of mesh triangle ti

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: v1, v2, v3, cc

    v1 = mesh%V( mesh%Tri( ti,1),:)
    v2 = mesh%V( mesh%Tri( ti,2),:)
    v3 = mesh%V( mesh%Tri( ti,3),:)
    cc = circumcenter( v1, v2, v3)

    ! If find_circumcenter yields infinity, it's because p and q have the
    ! same y-coordinate. Rearrange vertices in triangle matrix (maintaining
    ! counter-clockwise orientation) and try again.

    IF (cc( 1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      mesh%Tri( ti,:) = [mesh%Tri( ti,2), mesh%Tri( ti,3), mesh%Tri( ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:)
      cc = circumcenter( v1, v2, v3)
    END IF

    IF (cc( 1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      mesh%Tri( ti,:) = [mesh%Tri( ti,2), mesh%Tri( ti,3), mesh%Tri( ti,1)]
      v1 = mesh%V( mesh%Tri( ti,1),:)
      v2 = mesh%V( mesh%Tri( ti,2),:)
      v3 = mesh%V( mesh%Tri( ti,3),:)
      cc = circumcenter( v1, v2, v3)
    END IF

    IF (cc(1) > (mesh%xmax - mesh%xmin) * 1E5_dp .OR. cc( 2) > (mesh%ymax - mesh%ymin) * 1E5_dp) THEN
      CALL warning('update_triangle_circumcenter: triangle  doesn yield a valid circumcenter!')
    END IF

    mesh%Tricc( ti,:) = cc

  END SUBROUTINE update_triangle_circumcenter

! == Tools for handling the triangle refinement stack

  SUBROUTINE add_triangle_to_refinement_stack_first( mesh, ti)
    ! Add triangle ti to the top of the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    IF (mesh%refinement_map( ti) == 0) THEN
      mesh%refinement_map( ti) = 1
      mesh%refinement_stackN = mesh%refinement_stackN + 1
      mesh%refinement_stack( 2: mesh%nTri_mem) = mesh%refinement_stack( 1: mesh%nTri_mem-1)
      mesh%refinement_stack( 1) = ti
    END IF

  END SUBROUTINE add_triangle_to_refinement_stack_first

  SUBROUTINE add_triangle_to_refinement_stack_last( mesh, ti)
    ! Add triangle ti to the end of the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    IF (mesh%refinement_map( ti) == 0) THEN
      mesh%refinement_map( ti) = 1
      mesh%refinement_stackN = mesh%refinement_stackN + 1
      mesh%refinement_stack( mesh%refinement_stackN) = ti
    END IF

  END SUBROUTINE add_triangle_to_refinement_stack_last

  SUBROUTINE remove_triangle_from_refinement_stack( mesh, ti)
    ! Remove triangle ti from the refinement stack

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    INTEGER,                    INTENT(IN)        :: ti

    ! Local variables:
    INTEGER                                       :: i

    ! Safety
    IF (ti == 0) CALL crash('received ti=0!')

    mesh%refinement_map( ti) = 0

    i = 1
    DO WHILE (i <= mesh%refinement_stackN)

      IF (mesh%refinement_stack( i) == ti) THEN
        mesh%refinement_stack( 1:mesh%refinement_stackN) = [mesh%refinement_stack( 1:i-1), mesh%refinement_stack( i+1:mesh%refinement_stackN), 0]
        mesh%refinement_stackN = mesh%refinement_stackN - 1
      ELSE
        i = i+1
      END IF

    END DO

  END SUBROUTINE remove_triangle_from_refinement_stack

! == Some basic search operations on a mesh

  SUBROUTINE find_containing_triangle( mesh, p, ti_in)
    ! Find the triangle containing the point p. First do a "linear search":
    ! Start at initial guess ti_in. Check all neighbours of ti_in, find the one
    ! closest to p, select that one as the new ti_in. Repeat until all neighbours
    ! of ti_in are further away from p than ti_in itself.
    ! Then (if needed) search outward from there using a flood-fill algorithm.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: ti_in

    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: qq, rr, ss
    INTEGER                                       :: ncycle, t_prev
    REAL(dp), DIMENSION(2)                        :: gcti, gctc
    REAL(dp)                                      :: d, dc, dcmin
    INTEGER                                       :: tc, tcmin
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: map, stack1, stack2
    INTEGER                                       :: stackN1, stackN2
    LOGICAL                                       :: FoundIt
    INTEGER                                       :: n, ti, n2, tin

    ! If p lies outside the mesh domain, throw an error
    IF (p( 1) < mesh%xmin .OR. p( 1) > mesh%xmax .OR. p( 2) < mesh%ymin .OR. p( 2) > mesh%ymax) THEN
      CALL crash('p lies outside mesh domain!')
    END IF

    ! See if the initial guess is correct.
    qq = mesh%V( mesh%Tri( ti_in,1),:)
    rr = mesh%V( mesh%Tri( ti_in,2),:)
    ss = mesh%V( mesh%Tri( ti_in,3),:)
    IF (is_in_triangle( qq,rr,ss,p)) RETURN

    ! If not, start with a linear search.

  ! == Linear search ==
  ! ===================

    ncycle = 0
    t_prev = ti_in
    DO WHILE (ncycle < mesh%nTri)

      qq = mesh%V( mesh%Tri( ti_in,1),:)
      rr = mesh%V( mesh%Tri( ti_in,2),:)
      ss = mesh%V( mesh%Tri( ti_in,3),:)
      gcti = geometric_center( qq,rr,ss)
      d = NORM2( gcti - p)

      dcmin = d + 10._dp
      tcmin = 0
      DO n = 1, 3
        tc   = mesh%TriC( ti_in,n)
        IF (tc == 0)      CYCLE ! This triangle neighbour doesn't exist
        IF (tc == t_prev) CYCLE ! This triangle neighbour is the one we just came from
        qq = mesh%V( mesh%Tri( tc,1),:)
        rr = mesh%V( mesh%Tri( tc,2),:)
        ss = mesh%V( mesh%Tri( tc,3),:)
        gctc = geometric_center( qq,rr,ss)
        dc = NORM2( gctc - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          tcmin = tc
        END IF
      END DO

      IF (dcmin < d) THEN
        t_prev = ti_in
        ti_in = tcmin
      ELSE
        EXIT
      END IF

    END DO ! DO WHILE (ncycle < mesh%nTri)

    ! Check if the result from the linear search is correct.
    qq = mesh%V(mesh%Tri( ti_in,1),:)
    rr = mesh%V(mesh%Tri( ti_in,2),:)
    ss = mesh%V(mesh%Tri( ti_in,3),:)
    IF (is_in_triangle( qq, rr, ss, p)) RETURN
    IF (lies_on_line_segment( qq, rr, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( rr, ss, p, mesh%tol_dist)) RETURN
    IF (lies_on_line_segment( ss, qq, p, mesh%tol_dist)) RETURN

    ! It's not. Perform a flood-fill style outward search.

  ! == Flood-fill search ==
  ! =======================

    ALLOCATE( map(    mesh%nTri), source = 0)
    ALLOCATE( stack1( mesh%nTri), source = 0)
    ALLOCATE( stack2( mesh%nTri), source = 0)

    map( ti_in)  = 1 ! We checked that one.

    ! Add ti_in's neighbours to the stack.
    stackN1 = 0
    DO n = 1, 3
      IF (mesh%TriC( ti_in,n) > 0) THEN
        stackN1 = stackN1 + 1
        stack1( stackN1) = mesh%TriC( ti_in,n)
      END IF
    END DO

    FoundIt = .FALSE.
    DO WHILE (.NOT. FoundIt)
      ! Check all triangles in the stack. If they're not it, add their
      ! non-checked neighbours to the new stack.

      stack2  = 0
      stackN2 = 0

      DO n = 1, stackN1
        ti = stack1( n)
        qq = mesh%V( mesh%Tri( ti,1),:)
        rr = mesh%V( mesh%Tri( ti,2),:)
        ss = mesh%V( mesh%Tri( ti,3),:)
        IF (is_in_triangle( qq, rr, ss, p)) THEN

          ! Found it!
          FoundIt = .TRUE.
          ti_in = ti
          EXIT

        ELSE ! if (is_in_triangle(ti,p))
          ! Did not find it. And add this triangle's non-checked neighbours to the new stack.

          DO n2 = 1, 3
            tin = mesh%TriC( ti,n2)
            IF (tin == 0)       CYCLE ! This neighbour doesn't exist.
            IF (map( tin) == 1) CYCLE ! This neighbour has already been checked or is already in the stack.
            stackN2 = stackN2 + 1
            stack2( stackN2) = tin
            map( tin) = 1
          END DO

        END IF ! IF (is_in_triangle(q, r, s, p, tol)) THEN
      END DO ! DO n = 1, mesh%triStackN1

      ! Cycle stacks.
      stack1  = stack2
      stackN1 = stackN2

      ! If no more non-checked neighbours could be found, terminate and throw an error.
      IF (stackN2 == 0 .AND. .NOT. FoundIt) THEN
        CALL crash('couldnt find triangle containing p = [{dp_01}, {dp_02}]', dp_01 = p(1), dp_02 = p(2))
      END IF

    END DO ! DO WHILE (.NOT. FoundIt)

    ! Clean up after yourself
    DEALLOCATE( map   )
    DEALLOCATE( stack1)
    DEALLOCATE( stack2)

  END SUBROUTINE find_containing_triangle

  SUBROUTINE find_containing_vertex( mesh, p, vi)
    ! Find the vertex whose Voronoi cell contains the point p, using a "linear search"
    ! Start at initial guess vi. Check all neighbours of vi, find the one
    ! closest to p, select that one as the new vi. Repeat until all neighbours
    ! of vi are further away from p than vi itself.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(INOUT)       :: vi

    ! Local variables:
    INTEGER                                       :: ncycle, vi_prev, ci, vc, vcmin
    REAL(dp)                                      :: d, dc, dcmin


    ncycle = 0
    vi_prev = vi
    DO WHILE (ncycle < mesh%nV)

      d = NORM2( mesh%V( vi,:) - p)

      dcmin = d + 10._dp
      vcmin = 0
      DO ci = 1, mesh%nC( vi)
        vc = mesh%C( vi,ci)
        IF (vc == vi_prev) CYCLE ! This is the neighbour we just came from
        dc = NORM2( mesh%V( vc,:) - p)
        IF (dc < dcmin) THEN
          dcmin = dc
          vcmin = vc
        END IF
      END DO

      IF (dcmin < d) THEN
        vi_prev = vi
        vi = vcmin
      ELSE
        RETURN
      END IF

    END DO ! DO WHILE (ncycle < mesh%nV)

    ! If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
    CALL crash('couldnt find closest vertex!')

  END SUBROUTINE find_containing_vertex

  PURE FUNCTION is_in_Voronoi_cell( mesh, p, vi) RESULT( isso)
    ! If the point p lies closer to vertex vi than to any other vertex, then
    ! by definition it lies inside the Voronoi cell of vertex vi.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    REAL(dp), DIMENSION(  2), INTENT(IN)          :: p
    INTEGER,                  INTENT(IN)          :: vi

    ! Local variables:
    LOGICAL                                       :: isso
    REAL(dp)                                      :: dist_vi
    INTEGER                                       :: vvi, vj

    isso = .TRUE.

    dist_vi = NORM2( mesh%V( vi,:) - p)

    DO vvi = 1, mesh%nC( vi)
      vj = mesh%C( vi,vvi)
      IF (NORM2( mesh%V( vj,:) - p) < dist_vi) THEN
        isso = .FALSE.
        RETURN
      END IF
    END DO

  END FUNCTION is_in_Voronoi_cell

! == Spatial integration / averaging

  SUBROUTINE integrate_over_domain( mesh, d, int_d)
    ! Calculate the integral int_d over the model domain of a 2-D data field d

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)          :: mesh
    REAL(dp), DIMENSION( mesh%vi1:mesh%vi2), INTENT(IN)          :: d
    REAL(dp),                                INTENT(OUT)         :: int_d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'integrate_over_domain'
    INTEGER                                                      :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Integrate over process domain
    int_d = 0._dp
    DO vi = mesh%vi1, mesh%vi2
      int_d = int_d + d( vi) * mesh%A( vi)
    END DO

    ! Reduce over processes to find total domain integral
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, int_d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE integrate_over_domain

  SUBROUTINE average_over_domain( mesh, d, av_d)
    ! Calculate the average int_d over the model domain of a 2-D data field d

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)          :: mesh
    REAL(dp), DIMENSION( mesh%vi1:mesh%vi2), INTENT(IN)          :: d
    REAL(dp),                                INTENT(OUT)         :: av_d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                                :: routine_name = 'average_over_domain'
    REAL(dp)                                                     :: int_d

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Integrate over the domain
    CALL integrate_over_domain( mesh, d, int_d)

    ! Divide by domain area to find the average
    av_d = int_d / ((mesh%xmax - mesh%xmin) * (mesh%ymax - mesh%ymin))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE average_over_domain

! == Flood-fill operations

  SUBROUTINE extend_group_single_iteration_a( mesh, map, stack, stackN)
    ! Extend the group of vertices described by the provided map and stack
    ! outward by a single flood-fill iteration

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(mesh%nV),        INTENT(INOUT) :: map, stack
    INTEGER,                             INTENT(INOUT) :: stackN

    ! Local variables:
    INTEGER                                            :: n,i,vi,ci,vj

    ! Safety
    IF (stackN == 0) CALL crash('extend_group_single_iteration_a - needs at least one vertex to start!')

    n = stackN

    DO i = 1, n

      vi = stack( i)

      ! Safety
      IF (map( vi) == 0) CALL crash('extend_group_single_iteration_a - map and stack do not match!')

      ! Add all non-listed neighbours of vi
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        IF (map( vj) == 0) THEN
          map( vj) = 1
          stackN = stackN + 1
          stack( stackN) = vj
        END IF
      END DO

    END DO ! DO i = 1, n

  END SUBROUTINE extend_group_single_iteration_a

  SUBROUTINE extend_group_single_iteration_b( mesh, map, stack, stackN)
    ! Extend the group of triangles described by the provided map and stack
    ! outward by a single flood-fill iteration

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    INTEGER,  DIMENSION(mesh%nTri),      INTENT(INOUT) :: map, stack
    INTEGER,                             INTENT(INOUT) :: stackN

    ! Local variables:
    INTEGER                                            :: n,i,ti,n2,tj

    ! Safety
    IF (stackN == 0) CALL crash('extend_group_single_iteration_b - needs at least one vertex to start!')

    n = stackN

    DO i = 1, n

      ti = stack( i)

      ! Safety
      IF (map( ti) == 0) CALL crash('extend_group_single_iteration_b - map and stack do not match!')

      ! Add all non-listed neighbours of vi
      DO n2 = 1, 3
        tj = mesh%TriC( ti,n2)
        IF (tj == 0) CYCLE
        IF (map( tj) == 0) THEN
          map( tj) = 1
          stackN = stackN + 1
          stack( stackN) = tj
        END IF
      END DO

    END DO ! DO i = 1, n

  END SUBROUTINE extend_group_single_iteration_b

  SUBROUTINE extend_group_single_iteration_c( mesh, map, stack, stackN)
    ! Extend the group of edges described by the provided map and stack
    ! outward by a single flood-fill iteration

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    INTEGER,  DIMENSION(mesh%nE),        INTENT(INOUT) :: map, stack
    INTEGER,                             INTENT(INOUT) :: stackN

    ! Local variables:
    INTEGER                                            :: n,i,ei,vi,vj,vl,vr,ci,eil,ejl,eir,ejr

    ! Safety
    IF (stackN == 0) CALL crash('extend_group_single_iteration_c - needs at least one triangles to start!')

    n = stackN

    DO i = 1, n

      ei = stack( i)

      ! Safety
      IF (map( ei) == 0) CALL crash('extend_group_single_iteration_c - map and stack do not match!')

      ! Find the (up to) four nearby edges

      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)
      vl = mesh%EV( ei,3)
      vr = mesh%EV( ei,4)

      eil = 0
      ejl = 0
      eir = 0
      ejr = 0

      IF (vl > 0) THEN

        DO ci = 1, mesh%nC( vi)
          IF (mesh%C( vi,ci) == vl) eil = mesh%VE( vi,ci)
        END DO
        DO ci = 1, mesh%nC( vj)
          IF (mesh%C( vj,ci) == vl) ejl = mesh%VE( vj,ci)
        END DO

        ! Safety
        IF (eil == 0 .OR. ejl == 0) CALL crash('inconsistency in mesh%VE!')

        IF (map( eil) == 0) THEN
          map( eil) = 1
          stackN = stackN + 1
          stack( stackN) = eil
        END IF

        IF (map( ejl) == 0) THEN
          map( ejl) = 1
          stackN = stackN + 1
          stack( stackN) = ejl
        END IF

      END IF ! IF (vl > 0) THEN

      IF (vr > 0) THEN

        DO ci = 1, mesh%nC( vi)
          IF (mesh%C( vi,ci) == vr) eir = mesh%VE( vi,ci)
        END DO
        DO ci = 1, mesh%nC( vj)
          IF (mesh%C( vj,ci) == vr) ejr = mesh%VE( vj,ci)
        END DO

        ! Safety
        IF (eir == 0 .OR. ejr == 0) CALL crash('inconsistency in mesh%VE!')

        IF (map( eir) == 0) THEN
          map( eir) = 1
          stackN = stackN + 1
          stack( stackN) = eir
        END IF

        IF (map( ejr) == 0) THEN
          map( ejr) = 1
          stackN = stackN + 1
          stack( stackN) = ejr
        END IF

      END IF ! IF (vl > 0) THEN

    END DO ! DO i = 1, n

  END SUBROUTINE extend_group_single_iteration_c

! == Contours and polygons for mesh generation

  SUBROUTINE calc_mesh_contour_as_line( mesh, d, f, line, mask)
    ! Calculate a contour line at level f for data d on a mesh.
    ! Generate the contour line in UFEMISM line-segment format (i.e. unordered
    ! individual line segments).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d
    REAL(dp),                                INTENT(IN)    :: f
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: line
    LOGICAL,  DIMENSION(:    ), OPTIONAL,    INTENT(IN)    :: mask

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_mesh_contour_as_line'
    REAL(dp), PARAMETER                                    :: tol = 1E-5_dp
    REAL(dp), DIMENSION(:    ), ALLOCATABLE                :: d_scaled
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE                :: mask_loc
    INTEGER                                                :: n
    INTEGER                                                :: ei, vi, vj
    REAL(dp)                                               :: di, dj
    REAL(dp), DIMENSION(2)                                 :: p,q

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( d,1) /= mesh%nV) CALL crash('wrong vector size!')

    ! Trivial case: if all values of d are greater or smaller than f,
    ! the contour line is empty
    IF (MINVAL( d) >= f .OR. MAXVAL( d) <= f) THEN
      ALLOCATE( line( 0,0))
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Shift d so the contour lies at d_scaled = 0
    ALLOCATE( d_scaled( mesh%nV))
    d_scaled = d - f

    ! Set the mask to optionally skip certain grid cells
    ALLOCATE( mask_loc( mesh%nV_loc))
    IF (PRESENT( mask)) THEN
      mask_loc = mask
    ELSE
      mask_loc =  .TRUE.
    END IF

    ! Allocate memory for the line segments
    ALLOCATE( line( mesh%nE, 4))
    n = 0

    ! Go over all edges; if the contour lies on it, add it to the line segment list
    DO ei = 1, mesh%nE

      ! The four vertices adjacent to this edge
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)

      ! Skip this edge if told so
      IF (.NOT. mask_loc( vi) .AND. .NOT. mask_loc( vj)) CYCLE

      ! The values of d on these vertices
      di = d( vi)
      dj = d( vj)

      ! If the product is negative, the contour lies in between vi and vj
      IF (di * dj <= 0._dp) THEN
        ! Find the shared Voronoi cell boundary between vi and vj
        CALL find_shared_Voronoi_boundary( mesh, ei, p, q)
        ! Add this shared boundary as a line segment
        n = n + 1
        line( n,:) = [p(1), p(2), q(1), q(2)]
      END IF

    END DO ! DO ei = 1, mesh%nE

    ! Crop memory
    CALL reallocate( line, n, 4)

    ! Clean up after yourself
    DEALLOCATE( d_scaled)
    DEALLOCATE( mask_loc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_contour_as_line

  SUBROUTINE calc_mesh_mask_as_polygons( mesh, mask, poly_mult)
    ! Calculate a set of polygon enveloping all TRUE-valued mask cells

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    LOGICAL,  DIMENSION(:    ),              INTENT(IN)    :: mask
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: poly_mult

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_mesh_mask_as_polygons'
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE                :: mask_loc
    INTEGER                                                :: vi
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: poly
    INTEGER                                                :: n_poly, n_tot

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( mask,1) /= mesh%nV) CALL crash('incorrect data dimensions!')

    ! Trivial case for no TRUE values at all
    IF (.NOT. ANY( mask)) THEN
      ALLOCATE( poly_mult( 0,0))
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Make a local copy of the logical mask
    ALLOCATE( mask_loc( mesh%nV))
    mask_loc = mask

    ! Initialise poly_mult and poly
    ALLOCATE( poly_mult( mesh%nE,2))
    ALLOCATE( poly(      mesh%nE,2))
    n_tot = 0

    ! Calculate polygons for all TRUE regions of the mask
    DO vi = 1, mesh%nV

      IF (mask_loc( vi)) THEN
        ! Found a seed for a TRUE region

        ! Calculate a polygon enveloping this TRUE region, and
        ! remove the region from the mask
        CALL calc_mesh_mask_as_polygon( mesh, mask_loc, vi, poly, n_poly)

        ! Add this polygon to poly_mult
        poly_mult( n_tot+1,1) = REAL( n_poly,dp)
        poly_mult( n_tot+1,2) = 0._dp
        poly_mult( n_tot+2:n_tot+1+n_poly,:) = poly( 1:n_poly,:)
        n_tot = n_tot + 1 + n_poly

      END IF ! IF (mask_loc( vi)) THEN

    END DO

    ! Crop memory
    CALL reallocate( poly_mult, n_tot, 2)

    ! Clean up after yourself
    DEALLOCATE( mask_loc)
    DEALLOCATE( poly)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_mask_as_polygons

  SUBROUTINE calc_mesh_mask_as_polygon( mesh, mask, vi0, poly, n_poly)
    ! Calculate a polygon enveloping the set of TRUE-valued mask vertices around vi0,
    ! and remove that set of vertices from the mask

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    LOGICAL,  DIMENSION(:    ),              INTENT(INOUT) :: mask
    INTEGER,                                 INTENT(IN)    :: vi0
    REAL(dp), DIMENSION(:,:  ),              INTENT(OUT)   :: poly
    INTEGER,                                 INTENT(OUT)   :: n_poly

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'calc_mesh_mask_as_polygon'
    INTEGER,  DIMENSION(:    ), ALLOCATABLE                :: map, stack
    INTEGER                                                :: stackN
    INTEGER                                                :: vi, ci, vj, ei0, ei, it, cii, ci2, vk, ck, ck2, it2, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( mask,1) /= mesh%nV) CALL crash('incorrect data dimensions!')
    IF (.NOT. mask( vi0)) CALL crash('seed at vi0 is not TRUE!')

    ! Use a flood-fill algorithm to find the map of same-valued grid cells around mask cell vi0
    ALLOCATE( map(   mesh%nV), source = 0)
    ALLOCATE( stack( mesh%nV), source = 0)

    map( vi0) = 1
    stackN = 1
    stack( 1) = vi0

    DO WHILE (stackN > 0)

      ! Take the last element from the stack
      vi = stack( stackN)
      stackN = stackN - 1

      ! Mark it as mapped
      map( vi) = 2

      ! Remove it from the input mask
      mask( vi) = .FALSE.

      ! Add all non-mapped, non-stacked, TRUE-valued neighbours to the stack
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        ! IF this neighbour lies outside the TRUE region, store the connection
        ! as a starting point for the outline tracer
        IF (map( vj) == 0 .AND. .NOT. mask( vj) .AND. .NOT. (mesh%VBI( vi) > 0 .AND. mesh%VBI( vj) > 0)) ei0 = mesh%VE( vi,ci)
        ! IF this neighbour lies inside the TRUE region and isn't
        ! marked yet, add it to the stack and mark it
        IF (map( vj) == 0 .AND. mask( vj)) THEN
          ! Add this neighbour to the stack
          stackN = stackN + 1
          stack( stackN) = vj
          ! Mark this neighbour on the map as stacked
          map( vj) = 1
        END IF
      END DO

    END DO ! DO DO WHILE (stackN > 0)
    ! Safety
    IF (ei0 == 0) CALL crash('couldnt find starting edge!')

    ! Starting at the edge we found earlier, trace the outline of the TRUE region
    ei     = ei0
    n_poly = 0
    it     = 0

    DO WHILE (.TRUE.)

      ! Safety
      it = it + 1
      IF (it > mesh%nE) CALL crash('outline tracer got stuck!')

      ! Find the vertices vi and vj spanning the current edge ei,
      ! sorted such that vi lies inside the TRUE region and vj lies outside of it.
      IF     (map( mesh%EV( ei,1)) == 2 .AND. map( mesh%EV( ei,2)) == 0) THEN
        vi  = mesh%EV(   ei,1)
        vj  = mesh%EV(   ei,2)
        til = mesh%ETri( ei,1)
        tir = mesh%ETri( ei,2)
      ELSEIF (map( mesh%EV( ei,1)) == 0 .AND. map( mesh%EV( ei,2)) == 2) THEN
        vi  = mesh%EV(   ei,2)
        vj  = mesh%EV(   ei,1)
        til = mesh%ETri( ei,2)
        tir = mesh%ETri( ei,1)
      ELSE
        ! Apparently this edge doesn't cross the border of the TRUE region
        CALL crash('found non-border edge!')
      END IF

      ! Add this edge to the polygon
      n_poly = n_poly + 1
      IF (tir == 0) THEN
        ! vi-vj is a border edge?
        IF (mesh%VBI( vi) > 0 .AND. mesh%VBI( vj) > 0) THEN
          poly( n_poly,:) = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp
        ELSE
          CALL crash('expected vi-vj to be a border edge!')
        END IF
      ELSE
        poly( n_poly,:) = mesh%Tricc( tir,:)
      END IF

      ! Find the index ci such that VE( vi,ci) = ei
      ci = 0
      DO cii = 1, mesh%nC( vi)
        IF (mesh%VE( vi,cii) == ei) THEN
          ci = cii
          EXIT
        END IF
      END DO
      ! Safety
      IF (ci == 0) CALL crash('couldnt find connection index ci such that VE( vi,ci) = ei!')

      IF (mesh%VBI( vi) > 0 .AND. ci == mesh%nC( vi)) THEN
        ! Special case DO when the tracer reaches the domain border

        ! Move along the border until we find the other END of the TRUE region

        ! Add the last Voronoi cell boundary section
        n_poly = n_poly + 1
        poly( n_poly,:) = (mesh%V( vi,:) + mesh%V( vj,:)) / 2._dp

        ! Add the section from the Voronoi cell boundary to vi
        n_poly = n_poly + 1
        poly( n_poly,:) = mesh%V( vi,:)

        ! Add all border sections until we find the other END of the TRUE region
        it2 = 0
        DO WHILE (map( mesh%C( vi,1)) == 2)
          ! Safety
          it2 = it2 + 1
          IF (it2 > mesh%nV) CALL crash('outline tracer DO mesh border got stuck!')
          ! Move to next border vertex
          vi = mesh%C( vi,1)
          ! Add section to polygon
          n_poly = n_poly + 1
          poly( n_poly,:) = mesh%V( vi,:)
        END DO

        ! The next edge
        ei = mesh%VE( vi,1)

      ELSE ! IF (mesh%VBI( vi) > 0 .AND. ci == mesh%nC( vi))
        ! Regular case DO the domain interior

        ! The index ci2 of the next edge originating from vi, counterclockwise from ci
        ci2 = ci + 1
        IF (ci2 > mesh%nC( vi)) ci2 = 1

        ! The vertex that this edge points to
        vk = mesh%C( vi,ci2)

        IF (map( vk) == 2) THEN
          ! IF vk also lies inside the TRUE region, move to its Voronoi cell boundary

          ! Find the connection index ck such that C( vk,ci) = vi
          ck = 0
          DO cii = 1, mesh%nC( vk)
            IF (mesh%C( vk,cii) == vi) THEN
              ck = cii
              EXIT
            END IF
          END DO
          ! Safety
          IF (ck == 0) CALL crash('couldnt find connection index ck such that C( vk,ci) = vi')

          ! The index ck2 of the next edge originating from vk, counterclockwise from ck
          ck2 = ck + 1
          IF (ck2 > mesh%nC( vk)) ck2 = 1

          ! The next edge
          ei = mesh%VE( vk,ck2)

        ELSE ! IF (map( vk) == 2)
          ! IF vk lies outside the TRUE region, move to the next edge along vi

          ! The next edge
          ei = mesh%VE( vi,ci2)

        END IF ! IF (map( vk) == 2)

      END IF ! IF  IF (mesh%VBI( vi) > 0 .AND. ci == mesh%nC( vi))

      ! IF we've reached the starting point again, stop
      IF (it > 1) THEN
        IF (NORM2( poly( n_poly,:) - poly( 1,:)) < mesh%tol_dist) THEN
          EXIT
        END IF
      END IF

    END DO ! DO DO WHILE (.TRUE.)

    ! Clean up after yourself
    DEALLOCATE( map)
    DEALLOCATE( stack)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_mesh_mask_as_polygon

! == ISMIP-HOM periodic boundary conditions

  SUBROUTINE find_ti_copy_ISMIP_HOM_periodic( mesh, ti, ti_copy, wti_copy)
    ! Periodic boundary conditions in the ISMIP-HOM experiments are implemented by
    ! taking advantage of the fact that u(x,y) = u(x+L/2,y+L/2)
    !
    ! Velocities at the boundary can therefore be set equal to the interior value
    ! diagonally across from the boundary point (displaced by [L/2,L/2])
    !
    ! This routine finds the interior triangle to copy velocities from

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)              :: mesh
    INTEGER,                             INTENT(IN)              :: ti
    INTEGER,  DIMENSION(mesh%nC_mem),    INTENT(OUT)             :: ti_copy
    REAL(dp), DIMENSION(mesh%nC_mem),    INTENT(OUT)             :: wti_copy

    ! Local variables:
    REAL(dp), DIMENSION(2)                                       :: gc, p
    INTEGER                                                      :: vi, iti, tj
    REAL(dp)                                                     :: dist

    ! This triangle's geometric centre
    gc = mesh%TriGC( ti,:)

    ! The point where we want to copy the previous velocity solution
    IF (gc( 1) > 0._dp) THEN
      p( 1) = gc( 1) - C%refgeo_idealised_ISMIP_HOM_L / 2._dp
    ELSE
      p( 1) = gc( 1) + C%refgeo_idealised_ISMIP_HOM_L / 2._dp
    END IF
    IF (gc( 2) > 0._dp) THEN
      p( 2) = gc( 2) - C%refgeo_idealised_ISMIP_HOM_L / 2._dp
    ELSE
      p( 2) = gc( 2) + C%refgeo_idealised_ISMIP_HOM_L / 2._dp
    END IF

    ! The vertex whose Voronoi cell contains this point
    vi = 5
    CALL find_containing_vertex( mesh, p, vi)

    ! Weighted average over the triangles surrounding this vertex
    ti_copy  = 0
    wti_copy = 0._dp

    DO iti = 1, mesh%niTri( vi)
      tj = mesh%iTri( vi,iti)
      dist = NORM2( p - mesh%TriGC( tj,:))
      ti_copy(  iti) = tj
      wti_copy( iti) = 1._dp / dist**2
    END DO

    ! Normalise weights
    wti_copy( 1:mesh%niTri( vi)) = wti_copy( 1:mesh%niTri( vi)) / SUM( wti_copy( 1:mesh%niTri( vi)))

  END SUBROUTINE find_ti_copy_ISMIP_HOM_periodic

! == Diagnostic tools

  SUBROUTINE check_if_meshes_are_identical( mesh1, mesh2, isso)
    ! Check if two meshes are identical

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh1, mesh2
    LOGICAL,                             INTENT(OUT)   :: isso

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_if_meshes_are_identical'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: vi,ci,iti,ti,n

    ! Add routine to path
    CALL init_routine( routine_name)

    isso = .TRUE.

    ! Size
    IF (mesh1%nV /= mesh2%nV .OR. mesh1%nTri /= mesh2%nTri) THEN
      isso = .FALSE.
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Vertex coordinates
    DO vi = 1, mesh1%nV
      IF (NORM2( mesh1%V( vi,:) - mesh2%V( vi,:)) > tol) THEN
        isso = .FALSE.
        CALL finalise_routine( routine_name)
        RETURN
      END IF
    END DO

    ! Vertex-to-vertex connectivity
    DO vi = 1, mesh1%nV
      IF (mesh1%nC( vi) /= mesh2%nC( vi)) THEN
        isso = .FALSE.
        CALL finalise_routine( routine_name)
        RETURN
      END IF
      DO ci = 1, mesh1%nC( vi)
        IF (mesh1%C( vi,ci) /= mesh2%C( vi,ci)) THEN
          isso = .FALSE.
          CALL finalise_routine( routine_name)
          RETURN
        END IF
      END DO
    END DO

    ! Vertex-to-triangle connectivity
    DO vi = 1, mesh1%nV
      IF (mesh1%niTri( vi) /= mesh2%niTri( vi)) THEN
        isso = .FALSE.
        CALL finalise_routine( routine_name)
        RETURN
      END IF
      DO iti = 1, mesh1%niTri( vi)
        IF (mesh1%iTri( vi,iti) /= mesh2%iTri( vi,iti)) THEN
          isso = .FALSE.
          CALL finalise_routine( routine_name)
          RETURN
        END IF
      END DO
    END DO

    ! Triangle-to-vertex connectivity
    DO ti = 1, mesh1%nTri
      DO n = 1, 3
        IF (mesh1%Tri( ti,n) /= mesh2%Tri( ti,n)) THEN
          isso = .FALSE.
          CALL finalise_routine( routine_name)
          RETURN
        END IF
      END DO
    END DO

    ! Triangle circumcenter coordinates
    DO ti = 1, mesh1%nTri
      IF (NORM2( mesh1%Tricc( ti,:) - mesh2%Tricc( ti,:)) > tol) THEN
        isso = .FALSE.
        CALL finalise_routine( routine_name)
        RETURN
      END IF
    END DO

    ! Triangle-to-triangle connectivity
    DO ti = 1, mesh1%nTri
      DO n = 1, 3
        IF (mesh1%TriC( ti,n) /= mesh2%TriC( ti,n)) THEN
          isso = .FALSE.
          CALL finalise_routine( routine_name)
          RETURN
        END IF
      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_if_meshes_are_identical

  SUBROUTINE write_mesh_to_screen( mesh)
    ! Write a mesh to the screen. Best not to do this with large meshes.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh

    ! Local variables:
    INTEGER                                       :: vi, ti

    WRITE(0,*) '============================================================================'
    WRITE(0,*) ''
    WRITE(0,*) ' name = ', mesh%name
    WRITE(0,*) ' xmin = ', mesh%xmin
    WRITE(0,*) ' xmax = ', mesh%xmax
    WRITE(0,*) ' ymin = ', mesh%ymin
    WRITE(0,*) ' ymax = ', mesh%ymax
    WRITE(0,*) ''
    WRITE(0,*) ' vi    nC             C            niTri         iTri         VBI                x              y'
    DO vi = 1, mesh%nV
      WRITE(0,'(A,I3,A,I3,A,6I3,A,I3,A,6I3,A,I3,A,F12.1,A,F12.1)') &
      ' ', vi, '   ', mesh%nC(vi), '    ', mesh%C(vi,1:6), '    ', mesh%niTri(vi), '    ', mesh%iTri(vi,1:6), '    ', mesh%VBI(vi), &
      '    ', mesh%V(vi,1), '    ', mesh%V(vi,2)
    END DO

    WRITE(0,*) ' ti       Tri         TriC'
    DO ti = 1, mesh%nTri
      WRITE(0,'(A,I3,A,3I3,A,3I3)') &
      ' ', ti, '   ', mesh%Tri(ti,:), '    ', mesh%TriC(ti,:)
    END DO
    WRITE(0,*) '============================================================================'

  END SUBROUTINE write_mesh_to_screen

  SUBROUTINE write_mesh_to_text_file( mesh, filename)
    ! Write a mesh to a text file.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    CHARACTER(LEN=*),         INTENT(IN)          :: filename

    ! Local variables:
    INTEGER                                       :: vi, ci, ti, fp

    ! Create a new text file
    OPEN(newUNIT  = fp, FILE = TRIM( filename), STATUS = 'REPLACE')

    ! Header
    WRITE(UNIT = fp, FMT = '(A)')       ' Mesh data'
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  xmin    = ', mesh%xmin
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  xmax    = ', mesh%xmax
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  ymin    = ', mesh%ymin
    WRITE(UNIT = fp, FMT = '(A,F14.4)') '  ymax    = ', mesh%ymax
    WRITE(UNIT = fp, FMT = '(A,I2)')    '  nC_mem = ', mesh%nC_mem
    WRITE(UNIT = fp, FMT = '(A,I6)')    '  nV      = ', mesh%nV
    WRITE(UNIT = fp, FMT = '(A,I6)')    '  nTri    = ', mesh%nTri
    WRITE(UNIT = fp, FMT = '(A)')       ''
    WRITE(UNIT = fp, FMT = '(A,I6,A,I3,A)')       'Vertex data: ', mesh%nV, ' rows, ', 2 + 1 + mesh%nC_mem + 1 + mesh%nC_mem + 1, ' columns'
    WRITE(UNIT = fp, FMT = '(A)')       ''

    ! Vertex data
    WRITE(UNIT = fp, FMT = '(A)')       'V  nC  C  niTri  iTri  VBI'
    DO vi = 1, mesh%nV
      WRITE(UNIT = fp, FMT = '(2F24.14,I3)', ADVANCE = 'NO') mesh%V(vi,1), mesh%V(vi,2), mesh%nC(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%C(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%niTri(vi)
      DO ci = 1, mesh%nC_mem
        WRITE(UNIT = fp, FMT = '(I6)', ADVANCE = 'NO') mesh%iTri(vi,ci)
      END DO
      WRITE(UNIT = fp, FMT = '(I3)', ADVANCE = 'NO') mesh%VBI(vi)
      WRITE(UNIT = fp, FMT = '(A)') ''
    END DO
    WRITE(UNIT = fp, FMT = '(A)')       ''

    ! Triangle data
    WRITE(UNIT = fp, FMT = '(A)')       'Tri  TriC'
    DO ti = 1, mesh%nTri
      WRITE(UNIT = fp, FMT = '(6I6)') mesh%Tri(ti,1), mesh%Tri(ti,2), mesh%Tri(ti,3), mesh%TriC(ti,1), mesh%TriC(ti,2), mesh%TriC(ti,3)
    END DO

    ! Close the text file
    CLOSE(UNIT = fp)

  END SUBROUTINE write_mesh_to_text_file

  SUBROUTINE check_mesh( mesh)
    ! Check if the mesh data is self-consistent

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh

    ! Local variables:
    INTEGER                                       :: vi, ci, vc, ci2, vc2, iti, iti2, ti, n, v1, v2, v3, ti2, n2
    LOGICAL                                       :: FoundIt
    INTEGER                                       :: tj, tivia, tivib, tivic, tjvia, tjvib, tjvic

    ! == V
    ! =============================================================
    DO vi = 1, mesh%nV
      IF (mesh%V(vi,1) < mesh%xmin - mesh%tol_dist .OR. mesh%V(vi,1) > mesh%xmax + mesh%tol_dist .OR. &
          mesh%V(vi,2) < mesh%ymin - mesh%tol_dist .OR. mesh%V(vi,2) > mesh%ymax + mesh%tol_dist) THEN
        CALL warning('vertex {int_01} outside mesh domain! (x = [{dp_01}, {dp_02}, {dp_03}], y = [{dp_04}, {dp_05}, {dp_06}]', &
          int_01 = vi, dp_01 = mesh%xmin, dp_02 = mesh%V( vi,1), dp_03 = mesh%xmax, dp_04 = mesh%ymin, dp_05 = mesh%V( vi,2), dp_06 = mesh%ymax)
      END IF
    END DO

    ! == nC
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        IF (mesh%C(vi,ci) == 0) CALL warning('vertex {int_01} has fewer connections than nC says!', int_01 = vi)
      END DO
      DO ci = mesh%nC(vi)+1, mesh%nC_mem
        IF (mesh%C(vi,ci) > 0)  CALL warning('vertex {int_01} has more connections than nC says!', int_01 = vi)
      END DO
    END DO

    ! == C
    ! =============================================================
    DO vi = 1, mesh%nV
      DO ci = 1, mesh%nC(vi)
        vc = mesh%C(vi,ci)
        FoundIt = .FALSE.
        DO ci2 = 1, mesh%nC(vc)
          vc2 = mesh%C(vc,ci2)
          IF (vc2==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) CALL warning('vertex {int_01} is connected to {int_02}, but not the other way round!', int_01 = vi, int_02 = vc)
      END DO
    END DO

    ! == niTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        IF (mesh%iTri(vi,iti) == 0) CALL warning('vertex {int_01} has fewer iTriangles than niTri says!', int_01 = vi)
      END DO
      DO iti = mesh%niTri(vi)+1, mesh%nC_mem
        IF (mesh%iTri(vi,iti) > 0)  CALL warning('vertex {int_01} has more iTriangles than niTri says!', int_01 = vi)
      END DO
    END DO

    ! == iTri
    ! =============================================================
    DO vi = 1, mesh%nV
      DO iti = 1, mesh%niTri(vi)
        ti = mesh%iTri(vi,iti)
        FoundIt = .FALSE.
        DO n = 1, 3
          IF (mesh%Tri(ti,n)==vi) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) CALL warning('vertex {int_01} lists triangle {int_02} in iTri, but that triangle {int_03} doesnt contain vertex {int_04}!', &
          int_01 = vi, int_02 = ti, int_03 = ti, int_04 = vi)
      END DO

      IF (mesh%VBI(vi) == 0) THEN

        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          n = 0
          DO iti = 1, mesh%niTri(vi)
            ti = mesh%iTri(vi,iti)
            DO iti2 = 1, mesh%niTri(vc)
              ti2 = mesh%iTri(vc,iti2)
              IF (ti==ti2) THEN
                n = n+1
                EXIT
              END IF
            END DO
          END DO
          IF (.NOT. (n==2)) CALL warning('non-edge vertices {int_01} and {int_02} share {int_03} triangles', int_01 = vi, int_02 = vc, int_03 = n)
        END DO

      ELSE ! IF (mesh%VBI(vi) == 0) THEN

        DO ci = 1, mesh%nC(vi)
          vc = mesh%C(vi,ci)
          IF (mesh%VBI(vc)==0) THEN

            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF (.NOT. (n==2)) CALL warning('edge vertex {int_01} and non-edge vertex {int_02} share {int_03} triangles', int_01 = vi, int_02 = vc, int_03 = n)

          ELSE ! IF (mesh%VBI(vc)==0) THEN

            n = 0
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO iti2 = 1, mesh%niTri(vc)
                ti2 = mesh%iTri(vc,iti2)
                IF (ti==ti2) THEN
                  n = n+1
                  EXIT
                END IF
              END DO
            END DO
            IF (.NOT. is_border_edge( mesh, vi, vc)) CYCLE
            IF (.NOT. (n==1)) CALL warning('edge vertices {int_01} and {int_02} share {int_03} triangles', int_01 = vi, int_02 = vc,int_03 = n)

          END IF
        END DO

      END IF
    END DO

    ! == border index
    ! =============================================================
    DO vi = 1, mesh%nV

      IF (mesh%VBI(vi) == 0) THEN

        IF (mesh%V(vi,1) <= mesh%xmin .OR. mesh%V(vi,1) >= mesh%xmax .OR. mesh%V(vi,2) <= mesh%ymin .OR. mesh%V(vi,2) >= mesh%ymax) THEN
          CALL warning('vertex {int_01} has border index 0 but lies on or beyond the mesh domain border!', int_01 = vi)
        END IF

        ! First and last neighbours must be connected
        vc  = mesh%C(vi,1)
        vc2 = mesh%C(vi,mesh%nC(vi))

        FoundIt = .FALSE.
        DO ci = 1, mesh%nC(vc)
          IF (mesh%C(vc,ci)==vc2) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) CALL warning('vertex {int_01} has border index 0, but its first and last neighbours are not connected!', int_01 = vi)

      ELSEIF (mesh%VBI(vi) == 1) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 1 but does not lie on the N border!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==8 .OR. mesh%VBI(vc)==1)) THEN
          CALL warning('vertex {int_01} has border index 1 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==1 .OR. mesh%VBI(vc)==2)) THEN
          CALL warning('vertex {int_01} has border index 1 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 2) THEN

        IF (.NOT. vi==3) CALL warning('vertex {int_01} is listed as NE corner!', int_01 = vi)

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 2 but does not lie on the NE corner!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==8 .OR. mesh%VBI(vc)==1)) THEN
          CALL warning('vertex {int_01} has border index 2 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==3 .OR. mesh%VBI(vc)==4)) THEN
          CALL warning('vertex {int_01} has border index 2 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 3) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 3 but does not lie on the E border!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==2 .OR. mesh%VBI(vc)==3)) THEN
          CALL warning('vertex {int_01} has border index 3 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==3 .OR. mesh%VBI(vc)==4)) THEN
          CALL warning('vertex {int_01} has border index 3 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 4) THEN

        IF (.NOT. vi==2) CALL warning('vertex {int_01} is listed as SE corner!', int_01 = vi)

        IF (ABS(mesh%V(vi,1) - mesh%xmax) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 4 but does not lie on the SE corner!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==2 .OR. mesh%VBI(vc)==3)) THEN
          CALL warning('vertex {int_01} has border index 4 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==5 .OR. mesh%VBI(vc)==6)) THEN
          CALL warning('vertex {int_01} has border index 4 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 5) THEN

        IF (ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 5 but does not lie on the S border!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==4 .OR. mesh%VBI(vc)==5)) THEN
          CALL warning('vertex {int_01} has border index 5 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==5 .OR. mesh%VBI(vc)==6)) THEN
          CALL warning('vertex {int_01} has border index 5 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 6) THEN

        IF (.NOT. vi==1) CALL warning('vertex {int_01} is listed as SW corner!', int_01 = vi)

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymin) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 6 but does not lie on the SW corner!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==4 .OR. mesh%VBI(vc)==5)) THEN
          CALL warning('vertex {int_01} has border index 6 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==7 .OR. mesh%VBI(vc)==8)) THEN
          CALL warning('vertex {int_01} has border index 6 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 7) THEN

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 7 but does not lie on the W border!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==6 .OR. mesh%VBI(vc)==7)) THEN
          CALL warning('vertex {int_01} has border index 7 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==7 .OR. mesh%VBI(vc)==8)) THEN
          CALL warning('vertex {int_01} has border index 7 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      ELSEIF (mesh%VBI(vi) == 8) THEN

        IF (.NOT. vi==4) CALL warning('vertex {int_01} is listed as NW corner!', int_01 = vi)

        IF (ABS(mesh%V(vi,1) - mesh%xmin) > mesh%tol_dist .OR. ABS(mesh%V(vi,2) - mesh%ymax) > mesh%tol_dist) THEN
          CALL warning('vertex {int_01} has border index 8 but does not lie on the NW corner!', int_01 = vi)
        END IF
        vc = mesh%C(vi,1)
        IF (.NOT. (mesh%VBI(vc)==6 .OR. mesh%VBI(vc)==7)) THEN
          CALL warning('vertex {int_01} has border index 8 but its first connection doesnt have a matching border index!', int_01 = vi)
        END IF
        vc = mesh%C(vi,mesh%nC(vi))
        IF (.NOT. (mesh%VBI(vc)==1 .OR. mesh%VBI(vc)==2)) THEN
          CALL warning('vertex {int_01} has border index 8 but its last connection doesnt have a matching border index!', int_01 = vi)
        END IF

      END IF

    END DO

    ! == Tri
    ! =============================================================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        vi = mesh%Tri(ti,n)
        FoundIt = .FALSE.
        DO iti = 1, mesh%niTri(vi)
          IF (mesh%iTri(vi,iti) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) CALL warning('triangle {int_01} contains vertex {int_02}, but that vertex doesnt list ti as an iTri!', int_01 = ti, int_02 = vi)
      END DO

      v1 = mesh%Tri(ti,1)
      v2 = mesh%Tri(ti,2)
      v3 = mesh%Tri(ti,3)

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v2) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) CALL warning('triangle {int_01} contains unconnected vertices {int_02} and {int_03}!', int_01 = ti, int_02 = v1, int_03 = v2)

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v1)
        vc = mesh%C(v1,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) CALL warning('triangle {int_01} contains unconnected vertices {int_02} and {int_03}!', int_01 = ti, int_02 = v1, int_03 = v3)

      FoundIt = .FALSE.
      DO ci = 1, mesh%nC(v2)
        vc = mesh%C(v2,ci)
        IF (vc==v3) THEN
          FoundIt = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. FoundIt) CALL warning('triangle {int_01} contains unconnected vertices {int_02} and {int_03}!', int_01 = ti, int_02 = v2, int_03 = v3)

      tivia = mesh%Tri( ti,1)
      tivib = mesh%Tri( ti,2)
      tivic = mesh%Tri( ti,3)
      DO tj = ti+1, mesh%nTri
        tjvia = mesh%Tri( tj,1)
        tjvib = mesh%Tri( tj,2)
        tjvic = mesh%Tri( tj,3)
        IF (ANY( [tjvia,tjvib,tjvic] == tivia) .AND. &
            ANY( [tjvia,tjvib,tjvic] == tivib) .AND. &
            ANY( [tjvia,tjvib,tjvic] == tivic)) CALL warning('triangles {int_01} and {int_02} are made up of the same vertices!', int_01 = ti, int_02 = tj)
      END DO
    END DO

    ! == TriC
    ! =============================================================

    DO ti = 1, mesh%nTri
      DO n = 1, 3
        ti2 = mesh%TriC(ti,n)
        IF (ti2 == 0) THEN
          CYCLE
        END IF
        FoundIt = .FALSE.
        DO n2 = 1, 3
          IF (mesh%TriC(ti2,n2) == ti) THEN
            FoundIt = .TRUE.
            EXIT
          END IF
        END DO
        IF (.NOT. FoundIt) CALL warning('triangle {int_01} is connected to {int_02}, but not the other way round!', int_01 = ti, int_02 = ti2)
      END DO
    END DO

  END SUBROUTINE check_mesh

  SUBROUTINE check_if_triangle_already_exists( mesh, ti)
    ! Check if any duplicate triangles exist

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN)          :: mesh
    INTEGER,                  INTENT(IN)          :: ti

    ! Local variables:
    INTEGER                                       :: via, vib, vic, tj

    via = mesh%Tri( ti,1)
    vib = mesh%Tri( ti,2)
    vic = mesh%Tri( ti,3)

    DO tj = 1, mesh%nTri
      IF (tj == ti) CYCLE
      IF (ANY( via == mesh%Tri( tj,:)) .AND. ANY( vib == mesh%Tri( tj,:)) .AND. ANY( vic == mesh%Tri( tj,:))) THEN
        CALL crash('duplicate triangles detected at ti = {int_01}, tj = {int_02}', int_01 = ti, int_02 = tj)
      END IF
    END DO

  END SUBROUTINE check_if_triangle_already_exists

END MODULE mesh_utilities
