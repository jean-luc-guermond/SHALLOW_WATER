MODULE IDP_limiting_shallow_water
  USE input_data
  USE mesh_handling
  USE matrix_type
  USE shallow_water_functions
  USE problem_setup
  USE lin_alg_tools
  PUBLIC :: compute_bounds_and_low_flux, compute_bounds_and_low_flux_and_source, convex_limiting_proc
  PRIVATE
  LOGICAL :: once_bounds=.TRUE., once_limiting=.TRUE.
  REAL(KIND=8), DIMENSION(:), POINTER   :: vel2max, hmin, hmax, Kinmax, Q2max
  REAL(KIND=8), POINTER, DIMENSION(:)   :: hsmall 
  TYPE(matrice_bloc)                    :: mc_minus_ml, lij, stiff
  TYPE(matrice_bloc), DIMENSION(k_dim+1):: fctmat
  REAL(KIND=8), DIMENSION(:),   POINTER :: relaxi
  INTEGER, PARAMETER                    :: it_max_limiting = 2
  REAL(KIND=8), PARAMETER :: epsilon_hsmall = 0.01d0
CONTAINS
  SUBROUTINE compute_bounds_and_low_flux(dt,velocity,dijL,cij,FluxijL,lumped,diag,un,opt_utest)
    IMPLICIT NONE
    REAL(KIND=8)                                      :: dt
    TYPE(matrice_bloc)                                :: dijL
    TYPE(matrice_bloc), DIMENSION(k_dim)              :: cij
    TYPE(matrice_bloc), DIMENSION(:)                  :: FluxijL
    REAL(KIND=8), DIMENSION(mesh%np,k_dim)            :: velocity
    REAL(KIND=8), DIMENSION(mesh%np)                  :: lumped
    INTEGER,      DIMENSION(mesh%np)                  :: diag
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: un, utest
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), OPTIONAL :: opt_utest
    REAL(KIND=8), DIMENSION(inputs%syst_size)               :: ubarij, usource, usourceL
    REAL(KIND=8), DIMENSION(mesh%np)                        :: overh
    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji, overhij, lambdai, xx
    INTEGER :: comp, d, i, j, p

    IF (once_bounds) THEN
       ALLOCATE(hmin(mesh%np),hmax(mesh%np),vel2max(mesh%np),kinmax(mesh%np),Q2max(mesh%np))
       once_bounds=.FALSE.
    END IF
    overh = compute_one_over_h(un(:,1))

    !===initialize bounds
    DO i = 1, mesh%np
       vel2max(i) =  SUM(velocity(i,:)**2)
       Q2max(i) = SUM(un(i,2:k_dim+1)**2)
       kinmax(i) = Q2max(i)*overh(i)
    END DO
    hmax = un(:,1)
    hmin = hmax

    DO i = 1, mesh%np
       lambdai = 1.d0/(dijL%ia(i+1) - 1.d0 - dijL%ia(i)) !===Diagonal term is removed
       usource = 0.d0
       usourceL= 0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          DO comp = 1, inputs%syst_size
             DO d = 1, k_dim !===Taking care of the gas dynamics flux.  
                usource(comp) = usource(comp) - 2*cij(d)%aa(p)*velocity(i,d)*un(i,comp)*ratij
             END DO
             usource(comp) = usource(comp) - 2*dijL%aa(p)*un(i,comp)*ratij
          END DO
       END DO

       usourceL = 0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          ubarij = 0.d0
          DO comp = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                ubarij(comp) = ubarij(comp) &
                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji - velocity(i,d)*un(i,comp)*ratij)
                xx = xx &
                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji + velocity(i,d)*un(i,comp)*ratij)    
             END DO
             fluxijL(comp)%aa(p) = xx+ dijL%aa(p)*(un(j,comp)*ratji-un(i,comp)*ratij)
             IF (comp.NE.1) THEN !===Hydrostatic pressure
                ubarij(comp) = ubarij(comp) &
                     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2)
                fluxijL(comp)%aa(p) = fluxijL(comp)%aa(p) &
                     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2 &
                     + un(i,1)**2)
             END IF
             ubarij(comp) = ubarij(comp)/(2*dijL%aa(p)) + (un(j,comp)*ratji + un(i,comp)*ratij)/2
          END DO

          IF (i.NE.j) THEN
             ubarij = ubarij + dt*usource/lumped(i)
             !TESTTTTTTTTTTTTTTTT
             usourceL = usourceL + (2*dt*dijL%aa(p)/lumped(i))*ubarij
             !TESTTTTTTTTTTTTTTTT
          ELSE
             ubarij = un(i,:)+dt*usource/lumped(i)
             !TESTTTTTTTTTTTTTTTT
             usourceL = usourceL + ubarij*(1.d0+2*dt*dijL%aa(diag(i))/lumped(i))
             !TESTTTTTTTTTTTTTTTT
          END IF

          hmin(i) = MIN(hmin(i), ubarij(1))
          hmax(i) = MAX(hmax(i), ubarij(1))
          overhij = 2*MAX(ubarij(1),0.d0)/(ubarij(1)**2+MAX(ubarij(1),inputs%htiny)**2)
          vel2max(i) = MAX(vel2max(i),SUM(ubarij(2:k_dim+1)**2)*overhij**2)
          Kinmax(i) = MAX(kinmax(i),SUM(ubarij(2:k_dim+1)**2)*overhij)
          Q2max(i) = MAX(Q2max(i),SUM(ubarij(2:k_dim+1)**2))
       END DO

       IF (PRESENT(opt_utest)) THEN
          opt_utest(i,:) = usourceL
       END IF

    END DO
    !TESTTTTTTTTTT
    !WRITE(*,*) ' MAXVAL(vel2max)', MAXVAL(vel2max)
    if (PRESENT(opt_utest)) THEN
       DO i = 1, mesh%np
          IF (opt_utest(i,1) - hmin(i).lt. -1.d-5) THEN
             WRITE(*,*) ' BUG 1', opt_utest(i,1),  hmin(i), opt_utest(i,1)-hmin(i)
             stop
          ELSE IF (hmax(i) - opt_utest(i,1) .lt. -1.d-5) THEN
             WRITE(*,*) ' BUG 2', hmax(i), opt_utest(i,1), hmax(i)- opt_utest(i,1)
             !stop
          END IF
          !IF (psi_func(opt_utest(i,:),vel2max(i))< -1.d-10*vel2max(i)) THEN
          !   WRITE(*,*) 'BUG compute_bounds_and_low_flux', psi_func(opt_utest(i,:),vel2max(i)),vel2max(i), i
          !   !STOP
          !END IF
       END DO
    END if
    !TESTTTTTTTTTTTTTTT
  END SUBROUTINE compute_bounds_and_low_flux

  SUBROUTINE compute_bounds_and_low_flux_and_source(dt,velocity,dijL,cij,FluxijL, SiL, lumped,diag,un,opt_utest)
    IMPLICIT NONE
    REAL(KIND=8)                                      :: dt
    TYPE(matrice_bloc)                                :: dijL
    TYPE(matrice_bloc), DIMENSION(k_dim)              :: cij
    TYPE(matrice_bloc), DIMENSION(:)                  :: FluxijL
    REAL(KIND=8), DIMENSION(mesh%np,k_dim)            :: velocity
    REAL(KIND=8), DIMENSION(mesh%np)                  :: lumped
    INTEGER,      DIMENSION(mesh%np)                  :: diag
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: un, utest, SiL
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), OPTIONAL :: opt_utest
    REAL(KIND=8), DIMENSION(inputs%syst_size)               :: ubarij, usource, usourceL
    REAL(KIND=8), DIMENSION(mesh%np)                        :: overh
    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji, overhij, lambdai, xx
    INTEGER :: comp, d, i, j, p

    IF (once_bounds) THEN
       ALLOCATE(hmin(mesh%np),hmax(mesh%np),vel2max(mesh%np),kinmax(mesh%np),Q2max(mesh%np))
       once_bounds=.FALSE.
    END IF
    overh = compute_one_over_h(un(:,1))

    !===initialize bounds
    DO i = 1, mesh%np
       vel2max(i) =  SUM(velocity(i,:)**2)
       Q2max(i) = SUM(un(i,2:k_dim+1)**2)
       kinmax(i) = Q2max(i)*overh(i)
    END DO
    hmax = un(:,1)
    hmin = hmax

    SiL =0.d0
    DO i = 1, mesh%np
       lambdai = 1.d0/(dijL%ia(i+1) - 1.d0 - dijL%ia(i)) !===Diagonal term is removed
       usource = 0.d0
       usourceL= 0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          DO comp = 1, inputs%syst_size
             DO d = 1, k_dim !===Taking care of the gas dynamics flux.  
                usource(comp) = usource(comp) - 2*cij(d)%aa(p)*velocity(i,d)*un(i,comp)*ratij
             END DO
             usource(comp) = usource(comp) - 2*dijL%aa(p)*un(i,comp)*ratij
          END DO
       END DO

       usourceL = 0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          ubarij = 0.d0
          DO comp = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                ubarij(comp) = ubarij(comp) &
                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji - velocity(i,d)*un(i,comp)*ratij)
                xx = xx &
                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji + velocity(i,d)*un(i,comp)*ratij)    
             END DO
             fluxijL(comp)%aa(p) = xx + dijL%aa(p)*(un(j,comp)*ratji-un(i,comp)*ratij)
             IF (comp.NE.1) THEN !===Hydrostatic pressure
                ubarij(comp) = ubarij(comp) &
                     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2)
                !fluxijL(comp)%aa(p) = fluxijL(comp)%aa(p) &
                !     - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)**2/2 + un(i,1)**2/2)
                !SiL(i,comp) = SiL(i,comp) - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2 &
                !     -un(j,1)**2/2 + un(i,1)**2/2)
                fluxijL(comp)%aa(p) = fluxijL(comp)%aa(p) &
                     - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)*un(i,1))
                SiL(i,comp) = SiL(i,comp) - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2 &
                     -un(j,1)*un(i,1))
             END IF
             ubarij(comp) = ubarij(comp)/(2*dijL%aa(p)) + (un(j,comp)*ratji + un(i,comp)*ratij)/2
          END DO

          IF (i.NE.j) THEN
             ubarij = ubarij + dt*usource/lumped(i)
             !TESTTTTTTTTTTTTTTTT
             usourceL = usourceL + (2*dt*dijL%aa(p)/lumped(i))*ubarij
             !TESTTTTTTTTTTTTTTTT
          ELSE
             ubarij = un(i,:)+dt*usource/lumped(i)
             !TESTTTTTTTTTTTTTTTT
             usourceL = usourceL + ubarij*(1.d0+2*dt*dijL%aa(diag(i))/lumped(i))
             !TESTTTTTTTTTTTTTTTT
          END IF
          hmin(i) = MIN(hmin(i), ubarij(1))
          hmax(i) = MAX(hmax(i), ubarij(1))
          overhij = 2*MAX(ubarij(1),0.d0)/(ubarij(1)**2+MAX(ubarij(1),inputs%htiny)**2)
          vel2max(i) = MAX(vel2max(i),SUM(ubarij(2:k_dim+1)**2)*overhij**2)
       END DO

       IF (PRESENT(opt_utest)) THEN
          opt_utest(i,:) = usourceL
       END IF

    END DO
    !TESTTTTTTTTTT
    !WRITE(*,*) ' MAXVAL(vel2max)', MAXVAL(vel2max)
    if (PRESENT(opt_utest)) THEN
       DO i = 1, mesh%np
          IF (opt_utest(i,1) - hmin(i).lt. -1.d-5) THEN
             WRITE(*,*) ' BUG 1', opt_utest(i,1),  hmin(i), opt_utest(i,1)-hmin(i)
             stop
          ELSE IF (hmax(i) - opt_utest(i,1) .lt. -1.d-5) THEN
             WRITE(*,*) ' BUG 2', hmax(i), opt_utest(i,1), hmax(i)- opt_utest(i,1)
             !stop
          END IF
          !IF (psi_func(opt_utest(i,:),vel2max(i))< -1.d-10*vel2max(i)) THEN
          !   WRITE(*,*) 'BUG compute_bounds_and_low_flux', psi_func(opt_utest(i,:),vel2max(i)),vel2max(i), i
          !   !STOP
          !END IF
       END DO
    END if
    !TESTTTTTTTTTTTTTTT
  END SUBROUTINE compute_bounds_and_low_flux_and_source


  SUBROUTINE convex_limiting_proc(velocity,un,ulow,unext,FluxijH,FluxijL,&
       mass,lumped,diag,opt_src_shift)
    USE space_dim
    USE st_matrix
    USE CSR_transpose
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim)   :: velocity
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: un, ulow, unext
    TYPE(matrice_bloc)                       :: mass
    TYPE(matrice_bloc), DIMENSION(:)         :: FluxijL, FluxijH
    REAL(KIND=8), DIMENSION(mesh%np)         :: lumped
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size),OPTIONAL:: opt_src_shift
    REAL(KIND=8), DIMENSION(mesh%np)         :: vel2, q2, kin, lim
    INTEGER,      DIMENSION(mesh%np)         :: diag
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: du
    INTEGER :: k, it, i, ps, pe
    REAL(KIND=8) :: x
    INTEGER :: mm(1)
    REAL(KIND=8) :: DD

    IF (once_limiting) THEN
       once_limiting=.FALSE.
       !===Mass - lumped
       CALL duplicate(mass,mc_minus_ml)
       mc_minus_ml%aa = mass%aa
       mc_minus_ml%aa(diag) = mc_minus_ml%aa(diag) - lumped
       !===lij
       CALL duplicate(mass,lij)
       !===fctmat
       DO k = 1, inputs%syst_size
          CALL duplicate(mass,fctmat(k))
       END DO
       !===stiff (for relaxation)
       CALL compute_stiffness(mesh,stiff)

       !===Define relaxi for use in relaxation
       ALLOCATE(relaxi(mesh%np))
       DD = SUM(lumped) !==DD is measure of domain
       IF (k_dim==1) THEN
          relaxi = (lumped/DD)*SQRT(lumped/DD)
       ELSE IF (k_dim==2) THEN
          relaxi = (lumped/DD)**(0.75)
       END IF
    END IF

    IF (inputs%if_relax_bounds) THEN
       x = 1.d0
       CALL relax(0.d0,-x,un(:,1),hmin)
       CALL relax(1.d30,x,un(:,1),hmax)
       DO i = 1, mesh%np
          vel2(i) = SUM(velocity(i,:)**2) !=== v**2
          q2(i) = SUM(un(i,2:k_dim)**2)   !=== q**2=v**2*h**2
          kin(i) = SQRT(vel2(i)*q2(i))    !=== v*q=v**2*h
       END DO
       x=1.d0
       CALL relax(1.d30,x,vel2,vel2max)
       CALL relax(1.d30,x,kin,kinmax)
       CALL relax(1.d30,x,Q2,Q2max)
       !hmin = 0.d0
       !hmax = 1.d10
    END IF

    !===Time increment if consistent matrix is used
    IF (inputs%if_lumped) THEN
       du = 0.d0
    ELSE
       du = unext-un
    END IF

    !===Limiting matrix, viscosity + mass matrix correction
    CALL compute_fct_matrix_full(du,fctmat,FluxijH,FluxijL,mc_minus_ml)

    !===Convex limiting
    DO it = 1, it_max_limiting !and more is good for large CFL
       lij%aa = 1.d0
       !===Limit density
       CALL LOCAL_limit(ulow,hmax,hmin,fctmat(1),lumped,1,diag,lij)!===Works best

       !===Limit
       !CALL limit_v2(ulow,vel2max,lumped)              !===(V^2)_max - Q^2/h^2
       !CALL quadratic_limiting(ulow,lumped,vel2max)    !===(V_^2)_max h^2 - Q^2
       CALL quadratic_kin_energy_limiting(ulow,lumped) !=== (V.Q)_max h - Q^2
       !CALL quadratic_Q2_limiting(ulow,lumped)         !=== (Q^2)_max - Q^2

       !===Tranpose lij
       CALL transpose_op(lij,'min')
       !===TEST
       !lij%aa = 1.d0
       !lij%aa(diag)=0.d0
       !===TEST
       CALL update_fct_full(ulow,ulow,fctmat,lij,lumped)
       DO k = 1, inputs%syst_size
          fctmat(k)%aa = (1-lij%aa)*fctmat(k)%aa
       END DO
       IF (present(opt_Src_shift)) THEN
          lij%aa(diag)=1.d0
          DO i = 1, mesh%np
             ps = lij%ia(i)
             pe = lij%ia(i+1)-1
             opt_Src_shift(i,:) = opt_Src_shift(i,:)*MINVAL(lij%aa(ps:pe))
          END DO
       END IF
    END DO
    unext = ulow
    !===End of computation

  END SUBROUTINE convex_limiting_proc

  SUBROUTINE compute_fct_matrix_full(du,fctmatrix,FluxijH,FluxijL,mc_minus_ml)
    IMPLICIT NONE
    TYPE(matrice_bloc)                           :: mc_minus_ml
    TYPE(matrice_bloc), DIMENSION(:)             :: FluxijL, FluxijH
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)       :: fctmatrix
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: du
    INTEGER :: i, j, k, p
    DO k = 1, inputs%syst_size
       fctmatrix(k)%aa = inputs%dt*(FluxijH(k)%aa - FluxijL(k)%aa) 
    END DO
    DO i = 1, mesh%np
       DO p = mc_minus_ml%ia(i), mc_minus_ml%ia(i+1) - 1
          j = mc_minus_ml%ja(p)
          DO k = 1, inputs%syst_size
             fctmatrix(k)%aa(p) = fctmatrix(k)%aa(p) - mc_minus_ml%aa(p)*(du(j,k)-du(i,k))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix_full

  SUBROUTINE update_fct_full(ulow,unext,fctmat,lij,lumped)
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: unext
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)   :: fctmat
    TYPE(matrice_bloc)                                :: lij
    REAL(KIND=8), DIMENSION(mesh%np)                  :: lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: x
    INTEGER :: i, k, p
    DO i = 1, mesh%np
       x = 0.d0
       DO p = lij%ia(i), lij%ia(i+1) - 1
          DO k = 1, inputs%syst_size
             x(k) = x(k) + lij%aa(p)*fctmat(k)%aa(p)
          END DO
       END DO
       DO k = 1, inputs%syst_size
          unext(i,k) = ulow(i,k) + x(k)/lumped(i)
       END DO
    END DO
  END SUBROUTINE update_fct_full

  SUBROUTINE LOCAL_limit(ulow,maxn,minn,mat,mass,comp,diag,lij)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: mass, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)   :: mat
    INTEGER,                    INTENT(IN)   :: comp
    INTEGER, DIMENSION(:),     INTENT(IN)    :: diag
    TYPE(matrice_bloc),         INTENT(OUT)  :: lij
    REAL(KIND=8) :: maxni, minni, ui, uij, xij, lambdai, usmall
    INTEGER      :: i, j, p
    LOGICAL, SAVE :: once=.true.
  
    IF (once) THEN
       ALLOCATE(hsmall(mesh%np))
       hsmall = epsilon_hsmall*nondim_meshsize_loc*inputs%max_water_h
       once = .false.
    END IF

    !===Compute lij
    usmall = inputs%htiny
    lij%aa = 1.d0
    DO i = 1, mesh%np
       lambdai = 1.d0/(mat%ia(i+1) - 1.d0 - mat%ia(i))  !===Skip diagonal term
       !lambdai = 1.d0/(mat%ia(i+1) - mat%ia(i)) !===Do not skip diagonal term
       maxni = maxn(i)
       minni = minn(i)
       !write(*,*) minni, ulow(i,comp), maxni
       ui = ulow(i,comp)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          IF(i==j) CYCLE !===Skip diagonal term
          xij = mat%aa(p)/(mass(i)*lambdai)
          uij = ui + xij
          !IF (uij.LE. inputs%htiny) THEN
          IF (uij.LE.hsmall(i)) THEN
             lij%aa(p) = 0.d0 !threshold(uij/hsmall(i))
          ELSE IF (uij<minni) THEN
             lij%aa(p) = MIN(ABS(minni - ui)/(ABS(xij)+usmall),1.d0)
          ELSE IF (uij>maxni) THEN
             lij%aa(p) = MIN(ABS(maxni - ui)/(ABS(xij)+usmall),1.d0)
          END IF
       END DO
    END DO

    !lij%aa(diag) = 0.d0 !===Skip diagonal term
  END SUBROUTINE LOCAL_limit

  FUNCTION threshold(x) RESULT(g)
    USE mesh_handling
    USE sub_plot
    IMPLICIT NONE

    REAL(KIND=8)  :: x, z, t, zp, relu, f, g
    REAL(KIND=8), PARAMETER :: x0 = 0.25d0, x1=SQRT(3.d0)*x0 !x0=0.05 good for P1 (quadratic threshold)
!!$    REAL(KIND=8), PARAMETER :: x0 = 0.015d0 !x0=0.1 (cubic threshold)

!!$    REAL(KIND=8), DIMENSION(mesh%np)  :: xx, gg   
!!$    integer :: i
!!$    do i = 1, mesh%np
!!$       xx(i) = (i-1.d0)/mesh%np -1+2*x0
!!$       z = xx(i)-x0
!!$       zp = xx(i)-2*x0
!!$       relu = (zp+ABS(zp))/2
!!$       f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
!!$       gg(i) = max((f + 2*x0**3)/(4*x0**3),0.d0)
!!$       if (xx(i).le.0.d0) gg(i) =0.d0
!!$    end do
!!$    CALL plot_1d(xx,gg,'threshold1.plt') 
!!$    stop
    !===Quadratic threshold
    x = x - 1+2*x0
    if (x.le.0.d0) THEN
       g =0.d0
    ELSE
       z = x-x0   
       zp = x-2*x0
       relu = (zp+ABS(zp))/2
       f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
       g = (f + 2*x0**3)/(4*x0**3)
    END iF

    !===Cubic threshold
!!$    relu = ((x-2*x0)+abs(x-2*x0))/2
!!$    t = x/(2*x0)
!!$    g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)

    !CALL plot_1d(x,g,'threshold2.plt') 
    !stop
    RETURN
  END FUNCTION threshold

  SUBROUTINE quadratic_limiting(ulow,lumped,v2max)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: v2max, lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ur, Pij
    REAL(KIND=8) :: lambdai, psir, lr, a, b, c, delta, coeff, max_v2max, min_psi
    INTEGER      :: i, j, p, k
    !=== BUG max_v2max = inputs%gravity*(inputs%max_water_h)**2 !===g (h_0)^2
    max_v2max = inputs%gravity*(inputs%max_water_h) !===g (h_0)
    min_psi = max_v2max*inputs%htiny**2 !===(g h_0^2)*(htiny)^2
    DO i = 1, mesh%np
       lambdai = 1.d0/(lij%ia(i+1) - 1.d0 - lij%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       IF (v2max(i)*ulow(i,1)**2 .LE. min_psi) THEN !===Zero kinetic energy
          lij%aa(lij%ia(i):lij%ia(i+1) - 1) = 0.d0  !===Important for well-balancing
          CYCLE
       END IF
       DO p = lij%ia(i), lij%ia(i+1) - 1
          j =  lij%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(i,k) + lr*Pij(k)
          END DO
          psir = v2max(i)*ur(1)**2 - SUM(ur(2:k_dim+1)**2) !=== (V_max)^2 H^2 - ||Q||^2
          IF (psir.GE.0.d0) THEN
             CYCLE
          END IF

          a = -SUM(Pij(2:k_dim+1)**2) !===a is negative
          b = 2*(v2max(i)*Pij(1)*ulow(i,1) - SUM(Pij(2:k_dim+1)*ulow(i,2:k_dim+1)))
          c = v2max(i)*ulow(i,1)**2 - SUM(ulow(i,2:k_dim+1)**2)
          delta = b**2 - 4*a*c
          IF (delta<0.d0 .OR. a .GE. 0.d0) THEN
             CYCLE
          ELSE
             lr = MAX(0.d0,(-b-SQRT(delta))/(2*a))
             lij%aa(p) = MIN(lr,lij%aa(p))
          END IF
       END DO
    END DO
  END SUBROUTINE quadratic_limiting

  SUBROUTINE quadratic_kin_energy_limiting(ulow,lumped)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ur, Pij
    REAL(KIND=8) :: lambdai, psir, lr, a, b, c, delta, coeff, v2max, min_psi
    INTEGER      :: i, j, p, k
    v2max = inputs%gravity*inputs%max_water_h !===g * h_0
    min_psi = v2max*inputs%htiny**2 !===(g h_0)*(htiny)^2
    DO i = 1, mesh%np
       lambdai = 1.d0/(lij%ia(i+1) - 1.d0 - lij%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       IF (kinmax(i)*ulow(i,1) .LE. min_psi) THEN !===Zero kinetic energy
          lij%aa(lij%ia(i):lij%ia(i+1) - 1) = 0.d0  !===Important for well-balancing
          CYCLE
       END IF
       DO p = lij%ia(i), lij%ia(i+1) - 1
          j =  lij%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(i,k) + lr*Pij(k)
          END DO
          psir = kinmax(i)*ur(1) - SUM(ur(2:k_dim+1)**2) !=== Kmax H - ||Q||^2
          IF (psir.GE.0.d0) THEN
             CYCLE
          END IF

          a = -SUM(Pij(2:k_dim+1)**2) !===a is negative
          b = kinmax(i)*Pij(1) - 2*SUM(Pij(2:k_dim+1)*ulow(i,2:k_dim+1))
          c = kinmax(i)*ulow(i,1) - SUM(ulow(i,2:k_dim+1)**2)
          delta = b**2 - 4*a*c
          IF (delta<0.d0 .OR. a .GE. 0.d0) THEN
             CYCLE
          ELSE
             lr = MAX(0.d0,(-b-SQRT(delta))/(2*a))
             lij%aa(p) = MIN(lr,lij%aa(p))
          END IF
       END DO
    END DO
  END SUBROUTINE quadratic_kin_energy_limiting

  SUBROUTINE quadratic_Q2_limiting(ulow,lumped)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ur, Pij
    REAL(KIND=8) :: lambdai, psir, lr, a, b, c, delta, coeff, v2max, min_psi
    INTEGER      :: i, j, p, k
    v2max = inputs%gravity*inputs%max_water_h !===g * h_0
    min_psi = v2max*inputs%htiny**2 !===(g h_0)*(htiny)^2
    DO i = 1, mesh%np
       lambdai = 1.d0/(lij%ia(i+1) - 1.d0 - lij%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       IF (Q2max(i) .LE. min_psi) THEN !===Zero kinetic energy
          lij%aa(lij%ia(i):lij%ia(i+1) - 1) = 0.d0  !===Important for well-balancing
          CYCLE
       END IF
       DO p = lij%ia(i), lij%ia(i+1) - 1
          j =  lij%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(i,k) + lr*Pij(k)
          END DO
          psir = Q2max(i) - SUM(ur(2:k_dim+1)**2) !=== Q2max - ||Q||^2
          IF (psir.GE.0.d0) THEN
             CYCLE
          END IF

          a = -SUM(Pij(2:k_dim+1)**2) !===a is negative
          b = -2*SUM(Pij(2:k_dim+1)*ulow(i,2:k_dim+1))
          c = Q2max(i) - SUM(ulow(i,2:k_dim+1)**2)
          delta = b**2 - 4*a*c
          IF (delta<0.d0 .OR. a.GE.0.d0) THEN
             CYCLE
          ELSE
             lr = MAX(0.d0,(-b-SQRT(delta))/(2*a))
             lij%aa(p) = MIN(lr,lij%aa(p))
          END IF
       END DO
    END DO
  END SUBROUTINE quadratic_Q2_limiting


  SUBROUTINE relax(global_bound,pm,func,bound)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np)       :: func
    REAL(KIND=8),               INTENT(IN) :: global_bound
    REAL(KIND=8),               INTENT(IN) :: pm
    REAL(KIND=8), DIMENSION(:)             :: bound
    REAL(KIND=8), DIMENSION(mesh%np)       :: alpha, curvature
    INTEGER      :: i, j, p
    REAL(KIND=8) :: norm
    alpha = 0.d0
    DO i = 1, mesh%np
       norm = 0.d0
       DO p = stiff%ia(i), stiff%ia(i+1) - 1
          j = stiff%ja(p)
          !IF (i==j) CYCLE
          alpha(i) = alpha(i) + stiff%aa(p)*(func(j) - func(i))
          norm = norm + ABS(stiff%aa(p))
       END DO
       alpha(i) = 2*alpha(i)/norm
    END DO
    SELECT CASE(inputs%limiter_type)
    CASE('avg') !==Average
       curvature = 0.d0
       DO i = 1, mesh%np
          DO p = stiff%ia(i), stiff%ia(i+1) - 1
             j = stiff%ja(p)
             IF (i==j) CYCLE
             curvature(i) = curvature(i) + alpha(j) + alpha(i)
          END DO
          curvature(i) = curvature(i)/(2*(stiff%ia(i+1)-stiff%ia(i)-1))
       END DO
    CASE ('minmod') !===Minmod
       curvature = alpha    
       DO i = 1, mesh%np
          DO p = stiff%ia(i), stiff%ia(i+1) - 1
             j = stiff%ja(p)
             IF (i==j) CYCLE
             IF (curvature(i)*alpha(j).LE.0.d0) THEN
                curvature(i) = 0.d0
             ELSE IF (ABS(curvature(i)) > ABS(alpha(j))) THEN
                curvature(i) = alpha(j)
             END IF
          END DO
       END DO
    CASE ('none')
       curvature=alpha
    CASE DEFAULT
       WRITE(*,*) ' BUG in relax'
       STOP
    END SELECT

    IF (pm<0) THEN
       !bound = MAX(global_bound, (1-relaxi)*bound, bound + pm*ABS(curvature))
       bound = MAX((1-relaxi)*bound, bound + pm*ABS(curvature))
       !bound = MAX(global_bound,bound + pm*ABS(curvature))
       !bound = (1-relaxi)*bound
    ELSE
       !bound = MIN(global_bound, (1+relaxi)*bound, bound + pm*ABS(curvature))
       bound = MIN( (1+relaxi)*bound, bound + pm*ABS(curvature))
       !bound = MIN(global_bound, bound + pm*ABS(curvature))
       !bound = (1+relaxi)*bound
    END IF

  END SUBROUTINE RELAX

  SUBROUTINE limit_v2(ulow,v2m,lumped)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: v2m
    REAL(KIND=8), DIMENSION(mesh%np)                     :: lumped
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ul, ur, Pij
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, psi_small
    INTEGER      :: i, j, p, k
    psi_small = inputs%gravity*inputs%htiny !===g (h_0)
    !===
    DO i = 1, mesh%np
       lambdai = 1.d0/(lij%ia(i+1) - 1 - lij%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))

       DO p = lij%ia(i), lij%ia(i+1) - 1
          j =  lij%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p)
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(i,k) + lr*Pij(k) !===Density must be positive
          END DO
          if (ur(1)<0.d0) then
             lij%aa(p) = 0.d0  !===CFL is too large
             CYCLE
          end if
          psir = psi_func(ur,v2m(i))
          IF (psir.GE.-psi_small) THEN
             lij%aa(p) = lij%aa(p)
             CYCLE
          END IF
          ll = 0.d0
          ul = ulow(i,:)
          psil = psi_func(ul,v2m(i))
          DO WHILE (ABS(psil-psir) .GT. psi_small)
             llold = ll
             lrold = lr
             ll = ll - psil*(lr-ll)/(psir-psil)
             lr = lr - psir/psi_prime_func(Pij,lr,psir,ur,v2m(i))
             IF (ll.GE.lr) THEN
                ll = lr !lold
                EXIT
             END IF
             IF (ll< llold) THEN
                ll = llold
                EXIT
             END IF
             IF (lr > lrold) THEN
                lr = lrold
                EXIT
             END IF
             ul = ulow(i,:) + ll*Pij
             ur = ulow(i,:) + lr*Pij
             psil = psi_func(ul,v2m(i))
             psir = psi_func(ur,v2m(i))
          END DO
          IF (psir.GE.-psi_small) THEN
             lij%aa(p) = lr
          ELSE
             lij%aa(p) = ll
          END IF
       END DO
    END DO
  CONTAINS
    FUNCTION psi_func(u,v2m) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: u
      REAL(KIND=8),                     INTENT(IN) :: v2m
      REAL(KIND=8)                                 :: psi
      psi = v2m - SUM(u(2:inputs%syst_size)**2)*(2/(u(1)**2 + max(u(1),inputs%htiny)**2))
    END FUNCTION psi_func
    FUNCTION psi_prime_func(Pij,lr,psir,ur,V2m) RESULT(psi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: ur, Pij
    REAL(KIND=8),                              INTENT(IN) :: lr, psir, V2m
    REAL(KIND=8)                                 :: psi
    REAL(KIND=8), DIMENSION(inputs%syst_size)    :: up
    REAL(KIND=8), PARAMETER :: eps=1.d-7
    up = ur - eps*Pij
    psi = (psir-psi_func(up,V2m))/eps
  END FUNCTION psi_prime_func
  END SUBROUTINE limit_v2

!!$
!!$     SUBROUTINE compute_bounds_and_low_modified_flux(dt,velocity,dijL,cij,FluxijL, SiL, lumped,diag,un,opt_utest)
!!$    IMPLICIT NONE
!!$    REAL(KIND=8)                                      :: dt
!!$    TYPE(matrice_bloc)                                :: dijL
!!$    TYPE(matrice_bloc), DIMENSION(k_dim)              :: cij
!!$    TYPE(matrice_bloc), DIMENSION(:)                  :: FluxijL
!!$    REAL(KIND=8), DIMENSION(mesh%np,k_dim)            :: velocity
!!$    REAL(KIND=8), DIMENSION(mesh%np)                  :: lumped
!!$    INTEGER,      DIMENSION(mesh%np)                  :: diag
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: un, utest, SiL
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), OPTIONAL :: opt_utest
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size)               :: ubarij, usource, usourceL
!!$    REAL(KIND=8), DIMENSION(mesh%np)                        :: overh
!!$    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji, overhij, lambdai, xx
!!$    INTEGER :: comp, d, i, j, p
!!$
!!$    IF (once_bounds) THEN
!!$       ALLOCATE(hmin(mesh%np),hmax(mesh%np),vel2max(mesh%np))
!!$       once_bounds=.FALSE.
!!$    END IF
!!$    overh = compute_one_over_h(un(:,1))
!!$
!!$    !===initialize bounds
!!$    DO i = 1, mesh%np
!!$       vel2max(i) =  SUM(velocity(i,:)**2)
!!$    END DO
!!$    hmax = un(:,1)
!!$    hmin = hmax
!!$
!!$    SiL =0.d0
!!$    DO i = 1, mesh%np
!!$       lambdai = 1.d0/(dijL%ia(i+1) - 1.d0 - dijL%ia(i)) !===Diagonal term is removed
!!$       usource = 0.d0
!!$       usourceL= 0.d0
!!$       DO p = dijL%ia(i), dijL%ia(i+1) - 1
!!$          j = dijL%ja(p)
!!$          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
!!$          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
!!$          ratij = Hstarij*overh(i)
!!$          ratji = Hstarji*overh(j)
!!$          DO comp = 1, inputs%syst_size
!!$             DO d = 1, k_dim !===Taking care of the gas dynamics flux.  
!!$                usource(comp) = usource(comp) - 2*cij(d)%aa(p)*velocity(i,d)*un(i,comp)*ratij
!!$             END DO
!!$             usource(comp) = usource(comp) - 2*dijL%aa(p)*un(i,comp)*ratij
!!$          END DO
!!$       END DO
!!$
!!$       usourceL = 0.d0
!!$       DO p = dijL%ia(i), dijL%ia(i+1) - 1
!!$          j = dijL%ja(p)
!!$          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
!!$          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
!!$          ratij = Hstarij*overh(i)
!!$          ratji = Hstarji*overh(j)
!!$          ubarij = 0.d0
!!$          DO comp = 1, inputs%syst_size
!!$             xx = 0.d0
!!$             DO d = 1, k_dim
!!$                ubarij(comp) = ubarij(comp) &
!!$                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji - velocity(i,d)*un(i,comp)*ratij)
!!$                xx = xx &
!!$                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp)*ratji + velocity(i,d)*un(i,comp)*ratij)    
!!$             END DO
!!$             fluxijL(comp)%aa(p) = xx + dijL%aa(p)*(un(j,comp)*ratji-un(i,comp)*ratij)
!!$             IF (comp.NE.1) THEN !===Hydrostatic pressure
!!$                ubarij(comp) = ubarij(comp) &
!!$                     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2)
!!$                fluxijL(comp)%aa(p) = fluxijL(comp)%aa(p) &
!!$                     - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)**2/2 + un(i,1)**2/2)
!!$                SiL(i,comp) = SiL(i,comp) - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)*ratji)**2/2 - (un(i,1)*ratij)**2/2 &
!!$                     -un(j,1)**2/2 + un(i,1)**2/2)
!!$             END IF
!!$             ubarij(comp) = ubarij(comp)/(2*dijL%aa(p)) + (un(j,comp)*ratji + un(i,comp)*ratij)/2
!!$          END DO
!!$
!!$          IF (i.NE.j) THEN
!!$             ubarij = ubarij + dt*usource/lumped(i)
!!$             !TESTTTTTTTTTTTTTTTT
!!$             usourceL = usourceL + (2*dt*dijL%aa(p)/lumped(i))*ubarij
!!$             !TESTTTTTTTTTTTTTTTT
!!$          ELSE
!!$             ubarij = un(i,:)+dt*usource/lumped(i)
!!$             !TESTTTTTTTTTTTTTTTT
!!$             usourceL = usourceL + ubarij*(1.d0+2*dt*dijL%aa(diag(i))/lumped(i))
!!$             !TESTTTTTTTTTTTTTTTT
!!$          END IF
!!$          hmin(i) = MIN(hmin(i), ubarij(1))
!!$          hmax(i) = MAX(hmax(i), ubarij(1))
!!$          overhij = 2*MAX(ubarij(1),0.d0)/(ubarij(1)**2+MAX(ubarij(1),inputs%htiny)**2)
!!$          vel2max(i) = MAX(vel2max(i),SUM(ubarij(2:k_dim+1)**2)*overhij**2)
!!$       END DO
!!$
!!$       IF (PRESENT(opt_utest)) THEN
!!$          opt_utest(i,:) = usourceL
!!$       END IF
!!$
!!$    END DO
!!$    !TESTTTTTTTTTT
!!$    !WRITE(*,*) ' MAXVAL(vel2max)', MAXVAL(vel2max)
!!$    if (PRESENT(opt_utest)) THEN
!!$       DO i = 1, mesh%np
!!$          IF (opt_utest(i,1) - hmin(i).lt. -1.d-5) THEN
!!$             WRITE(*,*) ' BUG 1', opt_utest(i,1),  hmin(i), opt_utest(i,1)-hmin(i)
!!$             stop
!!$          ELSE IF (hmax(i) - opt_utest(i,1) .lt. -1.d-5) THEN
!!$             WRITE(*,*) ' BUG 2', hmax(i), opt_utest(i,1), hmax(i)- opt_utest(i,1)
!!$             !stop
!!$          END IF
!!$          !IF (psi_func(opt_utest(i,:),vel2max(i))< -1.d-10*vel2max(i)) THEN
!!$          !   WRITE(*,*) 'BUG compute_bounds_and_low_flux', psi_func(opt_utest(i,:),vel2max(i)),vel2max(i), i
!!$          !   !STOP
!!$          !END IF
!!$       END DO
!!$    END if
!!$    !TESTTTTTTTTTTTTTTT
!!$  END SUBROUTINE compute_bounds_and_low_modified_flux
!!$
!!$  SUBROUTINE newton_limiting(ulow,lumped)
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: ulow
!!$    REAL(KIND=8), DIMENSION(mesh%np)                     :: lumped
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ul, ur, Pij
!!$    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold
!!$    REAL(KIND=8) :: max_v2max, min_psi
!!$    INTEGER      :: i, j, p, k, it=0
!!$
!!$    max_v2max = inputs%gravity*(inputs%max_water_h)**2 !===g (h_0)^2
!!$    min_psi = max_v2max*inputs%htiny**2 !===(g h_0^2)*(htiny)^2
!!$
!!$    DO i = 1, mesh%np
!!$
!!$       !TESTTTTTTTTTT
!!$       !IF (psi_func(ulow(i,:),vel2max(i))<0.d0) THEN
!!$       !   WRITE(*,*) 'newton_limiting', psi_func(ulow(i,:),vel2max(i)), i
!!$       !END IF
!!$       !TESTTTTTTTTTT
!!$
!!$       lambdai = 1.d0/(lij%ia(i+1) - 1 - lij%ia(i))
!!$       coeff = 1.d0/(lambdai*lumped(i))
!!$
!!$       DO p = lij%ia(i), lij%ia(i+1) - 1
!!$          j =  lij%ja(p)
!!$          IF (i==j) THEN
!!$             lij%aa(p) = 0.d0
!!$             CYCLE
!!$          END IF
!!$          lr = lij%aa(p)
!!$          DO k = 1, inputs%syst_size
!!$             Pij(k) = fctmat(k)%aa(p)*coeff
!!$             ur(k) = ulow(i,k) + lr*Pij(k) 
!!$          END DO
!!$          IF (ur(1)<0.d0) THEN !===Water must be positive
!!$             lij%aa(p) = 0.d0
!!$             !write(*,*) 'Negative water', ur(1)
!!$             CYCLE
!!$          END IF
!!$          psir = psi_func(ur,vel2max(i))
!!$          IF (psir.GE.-min_psi) THEN
!!$             lij%aa(p) = lij%aa(p)
!!$             CYCLE
!!$          END IF
!!$          ll = 0.d0
!!$          ul = ulow(i,:)
!!$          psil = psi_func(ul,vel2max(i))
!!$          DO WHILE (psil-psir .GT. min_psi)
!!$             it = it +1
!!$             !TEST!!!!!!!!!!!!!!!!!!
!!$             !IF (i==154 .AND. it>1000) THEN
!!$             !   write(*,*) psi_func(ulow(i,:),vel2max(i)), psil, ulow(i,1)
!!$             !   ll =0
!!$             !   DO k = 1, 100
!!$             !      lrold = ll+(lr-ll)*(k-1)/99.d0
!!$             !      ur = ulow(i,:) + lrold*Pij
!!$             !      write(12,*) lrold, psi_func(ur,vel2max(i))
!!$             !   END DO
!!$             !   STOP
!!$             !   write(*,*) i, psil, psir, psi_prime_func(Pij,lr,psir,ur,Vel2max(i))
!!$             !END IF
!!$  
!!$             !TEST!!!!!!!!!!!!!!!!!!!1
!!$             llold = ll
!!$             lrold = lr
!!$             ll = ll - psil*(lr-ll)/(psir-psil)
!!$             lr = lr - psir/psi_prime_func(Pij,lr,psir,ur,Vel2max(i))
!!$
!!$
!!$             IF (ll.GE.lr) THEN
!!$                ll = lr !lold
!!$                EXIT
!!$             END IF
!!$             IF (ll< llold) THEN
!!$                ll = llold
!!$                EXIT
!!$             END IF
!!$             IF (lr > lrold) THEN
!!$                lr = lrold
!!$                EXIT
!!$             END IF
!!$             ul = ulow(i,:) + ll*Pij
!!$             ur = ulow(i,:) + lr*Pij
!!$             psil = psi_func(ul,vel2max(i))
!!$             psir = psi_func(ur,vel2max(i))
!!$          END DO
!!$          IF (psir.GE.-min_psi) THEN
!!$             lij%aa(p) = lr
!!$          ELSE
!!$             lij%aa(p) = ll
!!$          END IF
!!$       END DO
!!$    END DO
!!$
!!$    !TEST
!!$    !DO i = 1, mesh%np
!!$    !   ur = ulow(i,:)
!!$    !   DO p = lij%ia(i), lij%ia(i+1) - 1
!!$    !      DO k =1, inputs%syst_size
!!$    !        ur(k) =  ur(k) + lij%aa(p)*fctmat(k)%aa(p)/(lumped(i)*lambdai)
!!$    !      END DO
!!$    !   END DO
!!$    !   psir = psi_func(ur,vel2max(i))
!!$    !   if (psir<-min_psi) THEN
!!$    !      write(*,*) ' Negative Psi', i, psir
!!$    !      STOP
!!$    !   END if
!!$    END DO
!!$    !TEST
!!$    !CONTAINS
!!$  END SUBROUTINE newton_limiting
!!$
!!$  FUNCTION psi_func(u,v2max) RESULT(psi)
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: u
!!$    REAL(KIND=8),                     INTENT(IN) :: v2max
!!$    REAL(KIND=8)                                 :: psi, overh
!!$    overh = 2*MAX(u(1),0.d0)/(u(1)**2+MAX(u(1),inputs%htiny)**2)
!!$    psi = v2max - SUM(u(2:inputs%syst_size)**2)*overh**2
!!$  END FUNCTION psi_func
!!$
!!$  FUNCTION psi_prime_func(Pij,lr,psir,ur,V2max) RESULT(psi)
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN) :: ur, Pij
!!$    REAL(KIND=8),                              INTENT(IN) :: lr, psir, V2max
!!$    REAL(KIND=8)                                 :: psi
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size)    :: up
!!$    REAL(KIND=8), PARAMETER :: eps=1.d-7
!!$    up = ur - eps*Pij
!!$    psi = (psir-psi_func(up,V2max))/eps
!!$  END FUNCTION psi_prime_func

END MODULE IDP_limiting_shallow_water
