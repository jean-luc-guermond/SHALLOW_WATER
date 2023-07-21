MODULE lin_alg_tools
CONTAINS

  SUBROUTINE Axb(a,b,x)
    USE matrix_type
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: a
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: b
    REAL(KIND=8), DIMENSION(:), INTENT(OUT):: x
    INTEGER :: i, ps, pe
    DO i = 1, SIZE(a%ia)-1
       ps = a%ia(i)
       pe = a%ia(i+1)-1
       x(i) = SUM(a%aa(ps:pe)*b(a%ja(ps:pe)))
    END DO
  END SUBROUTINE Axb

  SUBROUTINE compute_mass(mesh,mass)
    USE matrix_type
    USE def_type_mesh
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)
  END SUBROUTINE compute_mass

  SUBROUTINE compute_stiffness(mesh,stiff)
    USE matrix_type
    USE def_type_mesh
    USE st_matrix
    USE fem_s_M
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: stiff
    CALL st_csr(mesh%jj, stiff%ia, stiff%ja)
    ALLOCATE(stiff%aa(SIZE(stiff%ja)))
    stiff%aa = 0.d0
    CALL qs_11_M (mesh, 1.d0, stiff%ia, stiff%ja, stiff%aa)
  END SUBROUTINE compute_stiffness

  SUBROUTINE lumped_mass(mesh,mass,lumped)
    USE matrix_type
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                     :: mesh
    TYPE(matrice_bloc)                  :: mass
    REAL(KIND=8), DIMENSION(:), POINTER :: lumped
    INTEGER :: i, np
    ALLOCATE(lumped(mesh%np))
    np = SIZE(mass%ia)-1
    ALLOCATE(lumped(np))
    DO i = 1, np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO
  END SUBROUTINE lumped_mass

  SUBROUTINE diag_mat(ia,ja,diag)
    IMPLICIT NONE
    INTEGER, DIMENSION(:)          :: ia,ja
    INTEGER, DIMENSION(:), POINTER :: diag
    INTEGER :: np, i, p
    np = SIZE(ia)-1
    ALLOCATE(diag(np))
    DO i = 1, np
       DO p = ia(i), ia(i+1) - 1
          IF (i==ja(p)) THEN
             diag(i) = p
             EXIT
          END IF
       END DO
    END DO
  END SUBROUTINE diag_mat

END MODULE lin_alg_tools
