MODULE fct
  USE matrix_type
  USE input_data
  USE problem_setup
  USE space_dim
CONTAINS
  SUBROUTINE FCT_positivity(ulow,hdry,maxn,minn,lumped,mat,lij)
    IMPLICIT NONE
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    TYPE(matrice_bloc),         INTENT(OUT) :: lij
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ulow, lumped
    LOGICAL,      DIMENSION(:), INTENT(IN)  :: hdry
    REAL(KIND=8), DIMENSION(:)              :: maxn, minn
    REAL(KIND=8), DIMENSION(SIZE(ulow))     :: Qminus, Pminus, Rminus
    REAL(KIND=8) :: fij, Pminus_small
    INTEGER      :: i, j, p
    qminus=MIN(lumped*(minn-ulow),0.d0)
    Pminus = 0.d0
    DO i = 1, SIZE(ulow)
       !===Go back to low-order if water height too small
       IF (hdry(i) .OR. minn(i) .LE. inputs%htiny ) THEN
          Rminus(i) = 0.d0
          CYCLE
       END IF
       
       Pminus_small = -inputs%htiny*lumped(i)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          IF (i==j) THEN
             CYCLE
          END IF
          fij = mat%aa(p)
          IF (fij.LT.0.d0) THEN
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       IF (pminus(i).GE.Pminus_small) THEN
          rminus(i) = 1.d0
       ELSE
          Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
       END IF
    END DO

    DO i = 1, SIZE(ulow)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          IF (i==j) THEN
             CYCLE
          END IF
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij%aa(p) = Rminus(j)
          ELSE
             lij%aa(p) = Rminus(i)
          END IF
       END DO
    END DO
    !WRITE(*,*) ' VERIF fct', SUM(lij%aa*mat%aa)
  END SUBROUTINE FCT_positivity

  SUBROUTINE FCT_generic(ulow,hdry,maxn,minn,mat,lumped,lij)
    IMPLICIT NONE
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    TYPE(matrice_bloc),         INTENT(OUT) :: lij
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ulow, lumped
    LOGICAL,      DIMENSION(:), INTENT(IN)  :: hdry
    REAL(KIND=8), DIMENSION(:)              :: maxn, minn
    REAL(KIND=8), DIMENSION(SIZE(ulow))     :: Qminus, Pminus, Rminus, Qplus, Pplus, Rplus
    REAL(KIND=8) :: fij, Pminus_small, Pplus_small
    INTEGER      :: i, j, p
    qminus = MIN(lumped*(minn-ulow),0.d0)
    qplus  = MAX(lumped*(maxn-ulow),0.d0)
    Pminus = 0.d0
    Pplus  = 0.d0
    DO i = 1, SIZE(ulow)
       !===Go back to low-order if water height too small
       IF (hdry(i)) THEN
          Rminus(i) = 0.d0
          Rplus(i)  = 0.d0
          CYCLE
       END IF
       Pminus_small = -inputs%htiny*lumped(i)
       Pplus_small  =  inputs%htiny*lumped(i)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          IF (i==j) THEN
             CYCLE
          END IF
          fij = mat%aa(p)
          IF (fij.LT.0.d0) THEN
             Pminus(i) = Pminus(i) + fij
          ELSE
             Pplus(i)  = Pplus(i) + fij
          END IF
       END DO
       IF (pminus(i).GE.Pminus_small) THEN
          rminus(i) = 1.d0
       ELSE
          Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
       END IF
       IF (Pplus(i).LE.Pplus_small) THEN
          Rplus(i) = 1.d0
       ELSE
          Rplus(i) =  MIN(Qplus(i)/Pplus(i),1.d0)
       END IF
    END DO

    DO i = 1, SIZE(ulow)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          IF (i==j) THEN
             CYCLE
          END IF
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij%aa(p) = MIN(Rminus(j),Rplus(i))
          ELSE
             lij%aa(p) = MIN(Rminus(i),Rplus(j))
          END IF
       END DO
    END DO
    !WRITE(*,*) ' VERIF fct', SUM(lij%aa*mat%aa)
  END SUBROUTINE FCT_generic
  
  SUBROUTINE maxmin(un,mat,maxn,minn)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: maxn, minn
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
    INTEGER      :: i
    DO i = 1, SIZE(un)
       maxn(i) = MAXVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
       minn(i) = MINVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
    END DO
  END SUBROUTINE maxmin

END MODULE fct
