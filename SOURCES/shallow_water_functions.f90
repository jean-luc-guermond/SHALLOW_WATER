MODULE shallow_water_functions
  USE mesh_handling
  USE space_dim
  USE input_data
  PUBLIC :: compute_velocity, flux, compute_one_over_h, compute_lambda_vacc

  PRIVATE
CONTAINS
  FUNCTION compute_velocity(un) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,1),k_dim) :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,1)) :: one_over_h, h
    REAL(KIND=8) :: small_h
    INTEGER :: k
    h = un(:,1)
    one_over_h = compute_one_over_h(h)
    DO k = 1, k_dim
       vv(:,k) = un(:,k+1)*one_over_h
    END DO
  END FUNCTION  compute_velocity

  FUNCTION flux(un) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim) :: vv
    REAL(KIND=8), DIMENSION(mesh%np,k_dim) :: velocity
    INTEGER :: k, ks
    velocity = compute_velocity(un)
    DO k = 1, k_dim
       DO ks = 1, inputs%syst_size
          vv(:,ks,k) = velocity(:,k)*un(:,ks)
       END DO
    END DO
  END FUNCTION flux

  FUNCTION compute_one_over_h(un) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: vv
    vv = 2*max(un,0.d0)/(un**2+max(un,inputs%htiny)**2)
  END FUNCTION compute_one_over_h

  SUBROUTINE compute_lambda_vacc(ul,ur,vell,velr,nij,lambda)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: lambda
    REAL(KIND=8), DIMENSION(k_dim), INTENT(IN)  :: nij, vell, velr
    REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN)   :: ur, ul
    REAL(KIND=8) :: ht, vl, vr, hl, hr, lbdl, lbdr, sql, sqr
    REAL(KIND=8) :: fh, a, c, Delta, x0, hmin, hmax, vmin, vmax, sqrmin, sqrmax
    REAL(KIND=8) :: ovhl, ovhr
    hl = max(ul(1),inputs%htiny)
    hr = max(ur(1),inputs%htiny)
    ovhl=1.d0/hl
    ovhr=1.d0/hr
    vl =  SUM(vell*nij)
    vr =  SUM(velr*nij)
    sql = SQRT(inputs%gravity*ABS(hl))
    sqr = SQRT(inputs%gravity*ABS(hr))

    x0=(2.d0*SQRT(2.d0)-1.d0)**2
    IF(hl.LE.hr) THEN
       hmin = hl
       vmin = vl
       hmax = hr
       vmax = vr
    ELSE
       hmin = hr
       vmin = vr
       hmax = hl
       vmax = vl
    END IF

    !===Case 1
    fh = phi(x0*hmin,vl,hl,vr,hr)
    IF (0.d0.LE.fh) THEN
       ht = min(x0*hmin,(MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity))
       lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
       lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
       lambda = MAX(ABS(lbdl),ABS(lbdr))
       RETURN
    END IF
    !===Case 2
    fh = phi(x0*hmax,vl,hl,vr,hr)
    IF (0.d0.LE.fh) THEN
       sqrmin = SQRT(hmin)
       sqrmax = SQRT(hmax)
       a = 1/(2.d0*SQRT(2.d0))
       c = -hmin*a -sqrmin*sqrmax + sqrmin*(vr-vl)/(2*sqrt(inputs%gravity))
       Delta = hmin-4*a*c
       IF (Delta<0.d0) THEN
          WRITE(*,*) ' BUG in compute_lambda_vacc'
          STOP
       END IF
       ht = min(x0*hmax,((-sqrmin+sqrt(Delta))/(2*a))**2)
       lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
       lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
       lambda = MAX(ABS(lbdl),ABS(lbdr))
       RETURN
    END IF
    !===Case 3
    sqrmin = SQRT(hmin)
    sqrmax = SQRT(hmax)
    ht = sqrmin*sqrmax*(1+sqrt(2/inputs%gravity)*(vl-vr)/(sqrmin+sqrmax))
    lbdl = vl - sql*SQRT((1+max((ht-hl)/(2*hl),0.d0))*(1+max((ht-hl)/hl,0.d0)))
    lbdr = vr + sqr*SQRT((1+max((ht-hr)/(2*hr),0.d0))*(1+max((ht-hr)/hr,0.d0)))
    lambda = MAX(ABS(lbdl),ABS(lbdr))
    RETURN
  END SUBROUTINE compute_lambda_vacc

  FUNCTION phi(h,ul,hl,ur,hr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: h, ul, hl, ur, hr
    REAL(KIND=8)             :: vv, fl, fr
    IF (h>hl) THEN
       fl = (h-hl)*SQRT((inputs%gravity/2)*(h+hl)/(h*hl))
    ELSE
       fl = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hl))
    END IF
    IF (h>hr) THEN
       fr = (h-hr)*SQRT((inputs%gravity/2)*(h+hr)/(h*hr))
    ELSE
       fr = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hr))
    END IF
    vv = fl + fr + ur - ul
  END FUNCTION phi

END MODULE shallow_water_functions
