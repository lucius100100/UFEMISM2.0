module assertions_unit_tests_mesh

  ! The assertions/unit tests for meshes.

  use assertions_unit_tests_basic, only: ASSERTION, UNIT_TEST, process_test_result
  use precisions, only: dp
  use mesh_types, only: type_mesh
  use control_resources_and_error_messaging, only: warning

  implicit none

  private

  public :: test_tol_mesh, test_mesh_is_self_consistent

contains

  !> Test if two meshes are identical to within tolerance
  subroutine test_tol_mesh( mesh1, mesh2, tol_dist, test_mode, message)
    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh1, mesh2
    real(dp),         intent(in   ) :: tol_dist
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result
    real(dp), dimension(:), allocatable :: dist

    test_result = .true.

    test_result = test_result .and. mesh1%nV          == mesh2%nV
    test_result = test_result .and. mesh1%nTri        == mesh2%nTri
    test_result = test_result .and. (mesh1%lambda_M    - mesh2%lambda_M   ) <= tol_dist
    test_result = test_result .and. (mesh1%phi_M       - mesh2%phi_M      ) <= tol_dist
    test_result = test_result .and. (mesh1%beta_stereo - mesh2%beta_stereo) <= tol_dist
    test_result = test_result .and. (mesh1%xmin        - mesh2%xmin       ) <= tol_dist
    test_result = test_result .and. (mesh1%xmax        - mesh2%xmax       ) <= tol_dist
    test_result = test_result .and. (mesh1%ymin        - mesh2%ymin       ) <= tol_dist
    test_result = test_result .and. (mesh1%ymax        - mesh2%ymax       ) <= tol_dist

    ! Vertex data
    if (mesh1%nV == mesh2%nV) then

      dist = hypot( mesh1%V(1:mesh1%nV,1) - mesh2%V(1:mesh2%nV,1), &
                    mesh1%V(1:mesh1%nV,2) - mesh2%V(1:mesh2%nV,2))
      test_result = test_result .and. all(dist <= tol_dist)
      test_result = test_result .and. all( mesh1%nC   ( 1:mesh1%nV  ) == mesh2%nC   ( 1:mesh2%nV  ))
      test_result = test_result .and. all( mesh1%C    ( 1:mesh1%nV,:) == mesh2%C    ( 1:mesh2%nV,:))
      test_result = test_result .and. all( mesh1%niTri( 1:mesh1%nV  ) == mesh2%niTri( 1:mesh2%nV  ))
      test_result = test_result .and. all( mesh1%iTri ( 1:mesh1%nV,:) == mesh2%iTri ( 1:mesh2%nV,:))
      test_result = test_result .and. all( mesh1%VBI  ( 1:mesh1%nV  ) == mesh2%VBI  ( 1:mesh2%nV  ))

    else
      test_result = .false.
    end if

    ! Triangle data
    if (mesh1%nTri == mesh2%nTri) then

      dist = hypot( mesh1%Tricc(1:mesh1%nTri,1) - mesh2%Tricc(1:mesh2%nTri,1), &
                    mesh1%Tricc(1:mesh1%nTri,2) - mesh2%Tricc(1:mesh2%nTri,2))
      test_result = test_result .and. all(dist <= tol_dist)
      test_result = test_result .and. all( mesh1%Tri ( 1:mesh1%nTri,:) == mesh2%Tri ( 1:mesh2%nTri,:))
      test_result = test_result .and. all( mesh1%TriC( 1:mesh1%nTri,:) == mesh2%TriC( 1:mesh2%nTri,:))

    else
      test_result = .false.
    end if

    call process_test_result( test_mode, test_result, message)

  end subroutine test_tol_mesh

  !> Test if a mesh is self-consistent
  subroutine test_mesh_is_self_consistent( mesh, test_mode, message)
    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh
    integer,          intent(in   ) :: test_mode
    character(len=*), intent(in   ) :: message
    ! Local variables:
    logical :: test_result

    test_result = is_self_consistent( mesh)

    call process_test_result( test_mode, test_result, message)

  end subroutine test_mesh_is_self_consistent

  !> Check if a mesh is self-consistent.
  function is_self_consistent( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res

    res = .true.

    res = res .and. is_self_consistent_vertex_locations( mesh)
    res = res .and. is_self_consistent_vertex_duplicates( mesh)
    res = res .and. is_self_consistent_nC( mesh)
    res = res .and. is_self_consistent_C( mesh)
    res = res .and. is_self_consistent_niTri( mesh)
    res = res .and. is_self_consistent_iTri( mesh)
    res = res .and. is_self_consistent_shared_triangles( mesh)
    res = res .and. is_self_consistent_border_VBI( mesh)
    res = res .and. is_self_consistent_border_V( mesh)
    res = res .and. is_self_consistent_Tri_iTri( mesh)
    res = res .and. is_self_consistent_Tri_C( mesh)
    res = res .and. is_self_consistent_triangle_duplicates( mesh)
    res = res .and. is_self_consistent_TriC( mesh)

  end function is_self_consistent

  !> Check if all vertices lie inside the mesh domain [xmin,xmax],[ymin,ymax]
  function is_self_consistent_vertex_locations( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res

    res = &
      all( mesh%V(:,1) >= mesh%xmin) .and. &
      all( mesh%V(:,1) <= mesh%xmax) .and. &
      all( mesh%V(:,2) >= mesh%ymin) .and. &
      all( mesh%V(:,2) <= mesh%xmax)

    if (.not. res) then
      call warning('Mesh inconsistency: found vertices outside the mesh domain')
    end if

  end function is_self_consistent_vertex_locations

  !> Check if any vertices are duplicate (i.e. coincide within the tolerance distance)
  function is_self_consistent_vertex_duplicates( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,vj

    res = .true.
    do vi = 1, mesh%nV-1
      do vj = vi+1, mesh%nV
        res = res .and. .not. hypot( mesh%V(vi,1) - mesh%V(vj,1), mesh%V(vi,2) - mesh%V(vj,2)) <= mesh%tol_dist
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: found duplicate vertices')
    end if

  end function is_self_consistent_vertex_duplicates

  !> Check if nC matches the number of entries in C
  function is_self_consistent_nC( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,ci

    res = .true.
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        res = res .and. mesh%C( vi,ci) > 0
      end do
      do ci = mesh%nC( vi)+1, mesh%nC_mem
        res = res .and. mesh%C( vi,ci) == 0
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: nC does not match number of entries in C')
    end if

  end function is_self_consistent_nC

  !> Check if vertices listed in C are connected both ways
  function is_self_consistent_C( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,ci,vj,cj,vk
    logical :: found_return_connection

    res = .true.
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        found_return_connection = .false.
        do cj = 1, mesh%nC( vj)
          vk = mesh%C( vj,cj)
          if (vk == vi) then
            found_return_connection = .true.
            exit
          end if
        end do
        res = res .and. found_return_connection
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: some vertices list neighbours in C that do not list a return connection')
    end if

  end function is_self_consistent_C

  !> Check if niTri matches the number of entries in iTri
  function is_self_consistent_niTri( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,iti

    res = .true.
    do vi = 1, mesh%nV
      do iti = 1, mesh%niTri( vi)
        res = res .and. mesh%iTri( vi,iti) > 0
      end do
      do iti = mesh%niTri( vi)+1, mesh%nC_mem
        res = res .and. mesh%iTri( vi,iti) == 0
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: niTri does not match number of entries in iTri')
    end if

  end function is_self_consistent_niTri

  !> Check if all iTriangles have a corresponding entry in the triangle-vertex connectivity list Tri
  function is_self_consistent_iTri( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,iti,ti,n,vj
    logical :: found_it

    res = .true.
    do vi = 1, mesh%nV
      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        found_it = .false.
        do n = 1, 3
          vj = mesh%Tri( ti,n)
          if (vj == vi) then
            found_it = .true.
            exit
          end if
        end do
        res = res .and. found_it
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: some vertices list triangles in iTri that do not list a return connection')
    end if

  end function is_self_consistent_iTri

  !> Check if connected vertices share the expected number of triangles
  function is_self_consistent_shared_triangles( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi,ci,vj,iti,itj,ti,tj,n_shared

    res = .true.
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi, ci)
        ! Loop over the iTriangles of both vi and vj and count the number of shared triangles
        n_shared = 0
        do iti = 1, mesh%niTri( vi)
          ti = mesh%iTri( vi,iti)
          do itj = 1, mesh%niTri( vj)
            tj = mesh%iTri( vj,itj)
            if (tj == ti) then
              n_shared = n_shared + 1
            end if
          end do
        end do
        ! If vi and vj both lie on the same border, they should share 1 triangle; otherwise, 2.
        if (is_border_edge( mesh,vi,vj)) then
          ! Both vi and vj lie on the border
          res = res .and. (n_shared == 1)
          if (n_shared /= 1) write(0,*) 'whaa'
        else
          ! Either vi or vj lies in the interior
          res = res .and. (n_shared == 2)
          if (n_shared /= 2) write(0,*) 'whoo'
        end if
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: connected vertices have an incorrect number of shared triangles')
    end if

  end function is_self_consistent_shared_triangles

  !> Check border indices of connected border vertices
  function is_self_consistent_border_VBI( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi, VBI, VBI_prev, VBI_prevprev, VBI_next, VBI_nextnext, vj_first, VBI_first, vj_last, VBI_last

    res = .true.

    ! Check if the first and last neighbours of border vertices
    ! are border vertices themselves too (and the correct ones at that)
    do vi = 1, mesh%nV

      VBI = mesh%VBI( vi)

      if (VBI > 0) then

        VBI_prev = VBI - 1
        if (VBI_prev == 0) VBI_prev = 8
        VBI_prevprev = VBI_prev - 1
        if (VBI_prevprev == 0) VBI_prevprev = 8

        VBI_next = VBI + 1
        if (VBI_next == 9) VBI_next = 1
        VBI_nextnext = VBI_next + 1
        if (VBI_nextnext == 9) VBI_nextnext = 1

        vj_first = mesh%C( vi,1)
        vj_last  = mesh%C( vi, mesh%nC( vi))
        VBI_first = mesh%VBI( vj_first)
        VBI_last  = mesh%VBI( vj_last)

        select case (mesh%VBI(vi))
        case (1,3,5,7)
          ! Border vertex (excluding corners)
          res = res .and. (VBI_first == VBI .or. VBI_first == VBI_prev)
          res = res .and. (VBI_last  == VBI .or. VBI_last  == VBI_next)
        case (2,4,6,8)
          ! Corner vertex
          res = res .and. (VBI_first == VBI_prev .or. VBI_first == VBI_prevprev)
          res = res .and. (VBI_last  == VBI_next .or. VBI_last  == VBI_nextnext)
        end select

      end if

    end do

    if (.not. res) then
      call warning('Mesh inconsistency: connected border vertices have incorrect border indices')
    end if

  end function is_self_consistent_border_VBI

  !> Check locations of border vertices
  function is_self_consistent_border_V( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: vi

    res = .true.

    ! Check if border vertices really lie on the border
    do vi = 1, mesh%nV
      select case (mesh%VBI( vi))
      case default
        ! Unrecognised border index
        res = .false.
      case (0)
        ! Free (interior) vertex
        res = res .and. &
          mesh%V( vi,1) > mesh%xmin .and. &
          mesh%V( vi,1) < mesh%xmax .and. &
          mesh%V( vi,2) > mesh%ymin .and. &
          mesh%V( vi,2) < mesh%ymax
      case (1)
        ! North border
        res = res .and. &
          mesh%V( vi,1) > mesh%xmin .and. &
          mesh%V( vi,1) < mesh%xmax .and. &
          abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
      case (5)
        ! South border
        res = res .and. &
          mesh%V( vi,1) > mesh%xmin .and. &
          mesh%V( vi,1) < mesh%xmax .and. &
          abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
      case (3)
        ! East border
        res = res .and. &
          mesh%V( vi,2) > mesh%ymin .and. &
          mesh%V( vi,2) < mesh%ymax .and. &
          abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist
      case (7)
        ! West border
        res = res .and. &
          mesh%V( vi,2) > mesh%ymin .and. &
          mesh%V( vi,2) < mesh%ymax .and. &
          abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist
      case (2)
        ! Northeast corner
        res = res .and. &
          abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist .and. &
          abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
      case (4)
        ! Southeast corner
        res = res .and. &
          abs( mesh%V( vi,1) - mesh%xmax) < mesh%tol_dist .and. &
          abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
      case (6)
        ! Southwest corner
        res = res .and. &
          abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist .and. &
          abs( mesh%V( vi,2) - mesh%ymin) < mesh%tol_dist
      case (8)
        ! Northwest corner
        res = res .and. &
          abs( mesh%V( vi,1) - mesh%xmin) < mesh%tol_dist .and. &
          abs( mesh%V( vi,2) - mesh%ymax) < mesh%tol_dist
      end select
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: border vertices do not lie on the border')
    end if

  end function is_self_consistent_border_V

  !> Check if vertices listed in Tri list those triangles in iTri
  function is_self_consistent_Tri_iTri( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: ti,n,vi,iti,tj
    logical :: found_it

    res = .true.
    do ti = 1, mesh%nTri
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        found_it = .false.
        do iti = 1, mesh%niTri( vi)
          tj = mesh%iTri( vi,iti)
          if (tj == ti) then
            found_it = .true.
            exit
          end if
        end do
        res = res .and. found_it
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: vertices listed in Tri do not list those triangles in iTri')
    end if

  end function is_self_consistent_Tri_iTri

  !> Check if vertices listed in Tri are connected
  function is_self_consistent_Tri_C( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: ti,n1,n2,vi,vj,ci,vk
    logical :: are_connected

    res = .true.
    do ti = 1, mesh%nTri
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1
        vi = mesh%Tri( ti,n1)
        vj = mesh%Tri( ti,n2)
        are_connected = .false.
        do ci = 1, mesh%nC( vi)
          vk = mesh%C( vi,ci)
          if (vk == vj) then
            are_connected = .true.
            exit
          end if
        end do
        res = res .and. are_connected
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: vertices listed in Tri are not connected')
    end if

  end function is_self_consistent_Tri_C

  !> Check if any triangles are duplicate (i.e. consist of the same three vertices)
  function is_self_consistent_triangle_duplicates( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: ti,tj,n,vi,n_match

    res = .true.
    do ti = 1, mesh%nTri-1
      do tj = ti+1, mesh%nTri
        n_match = 0
        do n = 1, 3
          vi = mesh%Tri( ti,n)
          if (any(mesh%Tri( tj,:) == vi)) n_match = n_match + 1
        end do
        res = res .and. n_match >= 0 .and. n_match <= 2
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: found duplicate triangles')
    end if

  end function is_self_consistent_triangle_duplicates

  !> Check if triangles listed in TriC are connected both ways
  function is_self_consistent_TriC( mesh) result( res)

    ! In/output variables:
    type(type_mesh), intent( in) :: mesh
    logical                      :: res
    ! Local variables:
    integer :: ti,n,tj,n2,tk
    logical :: found_it

    res = .true.
    do ti = 1, mesh%nTri
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (tj == 0) cycle
        found_it = .false.
        do n2 = 1, 3
          tk = mesh%TriC( tj,n2)
          if (tk == ti) then
            found_it = .true.
            exit
          end if
        end do
        res = res .and. found_it
      end do
    end do

    if (.not. res) then
      call warning('Mesh inconsistency: triangles list neighbours in TriC that do not list a return connection')
    end if

  end function is_self_consistent_TriC

  !> Check if the edge between two vertices is a border edge
  function is_border_edge( mesh, vi, vj) result(res)

    ! In/output variables:
    type(type_mesh), intent(in) :: mesh
    integer,         intent(in) :: vi, vj
    logical                     :: res

    res = .false.

    if (mesh%VBI( vi) == 0 .or. mesh%VBI( vj) == 0) then
      ! At least one of the vertices is an interior vertex; trivial answer
      return
    end if

    if (mesh%VBI( vi) == 1) then
      if (mesh%VBI( vj) == 8 .or. &
          mesh%VBI( vj) == 1 .or. &
          mesh%VBI( vj) == 2) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 2) then
      if (mesh%VBI( vj) == 8 .or. &
          mesh%VBI( vj) == 1 .or. &
          mesh%VBI( vj) == 3 .or. &
          mesh%VBI( vj) == 4) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 3) then
      if (mesh%VBI( vj) == 2 .or. &
          mesh%VBI( vj) == 3 .or. &
          mesh%VBI( vj) == 4) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 4) then
      if (mesh%VBI( vj) == 2 .or. &
          mesh%VBI( vj) == 3 .or. &
          mesh%VBI( vj) == 5 .or. &
          mesh%VBI( vj) == 6) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 5) then
      if (mesh%VBI( vj) == 4 .or. &
          mesh%VBI( vj) == 5 .or. &
          mesh%VBI( vj) == 6) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 6) then
      if (mesh%VBI( vj) == 4 .or. &
          mesh%VBI( vj) == 5 .or. &
          mesh%VBI( vj) == 7 .or. &
          mesh%VBI( vj) == 8) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 7) then
      if (mesh%VBI( vj) == 6 .or. &
          mesh%VBI( vj) == 7 .or. &
          mesh%VBI( vj) == 8) then
        res = .true.
      end if
    elseif (mesh%VBI( vi) == 8) then
      if (mesh%VBI( vj) == 6 .or. &
          mesh%VBI( vj) == 7 .or. &
          mesh%VBI( vj) == 1 .or. &
          mesh%VBI( vj) == 2) then
        res = .true.
      end if
    end if

  end function is_border_edge

end module assertions_unit_tests_mesh
