MODULE IDP_update_shallow_water
  USE matrix_type
  USE space_dim
  USE mesh_handling
  USE pardiso_solve
  USE Butcher_tableau
  USE input_data

  PUBLIC:: IDP_construct_shallow_water_matrices, full_step_ERK, bdy, euler, compute_dt
  INTEGER,  PUBLIC :: isolve_shallow_water_pardiso
  TYPE(BT), PUBLIC :: ERK
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: lumped
  PRIVATE

  INTEGER, DIMENSION(:), POINTER      :: diag
  TYPE(matrice_bloc)                  :: mass
  TYPE(matrice_bloc), DIMENSION(k_dim)  :: cij
  TYPE(matrice_bloc)                    :: dijL, dijH, muijL, muijH, stiff
  TYPE(matrice_bloc), DIMENSION(:), POINTER  :: Fluxij
  TYPE(matrice_bloc), DIMENSION(:), POINTER  :: FluxijL, FluxijH
  TYPE(matrice_bloc), DIMENSION(:,:), POINTER:: Fluxij_High
  REAL(KIND=8), DIMENSION(:,:,:), POINTER :: SiH
  REAL(KIND=8), DIMENSION(:,:),   POINTER :: SiL
  REAL(KIND=8), DIMENSION(:,:), POINTER :: velocity

CONTAINS 

  SUBROUTINE IDP_construct_shallow_water_matrices
    USE st_matrix
    USE mesh_handling
    USE fem_s_M
    USE lin_alg_tools
    USE CSR_transpose
    IMPLICIT NONE
    INTEGER :: m, p, ni, nj, i, j, d, k, l
    REAL(KIND=8), DIMENSION(k_dim) :: x
    LOGICAL, DIMENSION(mesh%np) :: if_udotn_zero


    !===Mass
    CALL compute_mass(mesh,mass)

    !===lumped
    CALL lumped_mass(mesh,mass,lumped)

    !===diag
    CALL diag_mat(mass%ia,mass%ja,diag)

    !===stiff (for smoothness indicator)
    CALL compute_stiffness(mesh,stiff)

    !===fluxij
    ALLOCATE(fluxij(inputs%syst_size))
    DO k = 1, inputs%syst_size
       CALL duplicate(mass,fluxij(k))
    END DO

    !===FluxijL
    ALLOCATE(fluxijL(inputs%syst_size))
    DO k = 1, inputs%syst_size
       CALL duplicate(mass,fluxijL(k))
    END DO

    !===FluxijH
    ALLOCATE(fluxijH(inputs%syst_size))
    DO k = 1, inputs%syst_size
       CALL duplicate(mass,fluxijH(k))
    END DO

    !===fluxij_High
    ALLOCATE(fluxij_High(inputs%syst_size,ERK%s))
    DO l = 1, ERK%s
       DO k = 1, inputs%syst_size
          CALL duplicate(mass,fluxij_High(k,l))
       END DO
    END DO

    !===dijL
    CALL duplicate(mass,dijL)

    !===dijH
    CALL duplicate(mass,dijH)

    !===muijL
    CALL duplicate(mass,muijL)

    !===muijH
    CALL duplicate(mass,muijH)

    !===cij = \int_K \GRAD(\phi_j) \phi_i \dif x
    if_udotn_zero = .FALSE.
    if_udotn_zero(udotn_js_D) = .TRUE.
    if_udotn_zero(h_js_D) = .FALSE.
    DO d = 1, k_dim
       CALL duplicate(mass,cij(d))
    END DO
    DO m = 1, mesh%me
       DO ni = 1, mesh%gauss%n_w  
          i = mesh%jj(ni, m)
          IF (if_udotn_zero(i)) THEN
             DO nj = 1, mesh%gauss%n_w  
                j = mesh%jj(nj, m)
                DO d = 1, k_dim
                   x(d) = -SUM(mesh%gauss%dw(d,ni,:,m) * mesh%gauss%ww(nj,:)*mesh%gauss%rj(:,m))
                END DO
                DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                   IF (cij(1)%ja(p) == j) THEN  
                      DO d = 1, k_dim
                         cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                      END DO
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             DO nj = 1, mesh%gauss%n_w  
                j = mesh%jj(nj, m)
                DO d = 1, k_dim
                   x(d) = SUM(mesh%gauss%dw(d,nj,:,m) * mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
                END DO
                DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                   IF (cij(1)%ja(p) == j) THEN  
                      DO d = 1, k_dim
                         cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                      END DO
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          END IF
       ENDDO
    ENDDO

    !===Velocity
    ALLOCATE(velocity(mesh%np,k_dim))
    ALLOCATE(SiH(mesh%np,inputs%syst_size,ERK%s))
    ALLOCATE(SiL(mesh%np,inputs%syst_size))

  END SUBROUTINE IDP_construct_shallow_water_matrices

  SUBROUTINE divide_by_lumped(if_lumped,rk)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: if_lumped
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: rk
    REAL(KIND=8), DIMENSION(mesh%np) :: ff
    INTEGER :: k
    IF (if_lumped) THEN
       DO k = 1, inputs%syst_size
          rk(:,k) = rk(:,k)/lumped
       END DO
    ELSE
       DO k = 1, inputs%syst_size
          CALL solve_pardiso(mass%aa,mass%ia,mass%ja,rk(:,k),ff,isolve_shallow_water_pardiso)
          isolve_shallow_water_pardiso=ABS(isolve_shallow_water_pardiso)
          rk(:,k) = ff
       END DO
    END IF
  END SUBROUTINE divide_by_lumped

  SUBROUTINE full_step_ERK(un_in)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: un_in
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,ERK%s+1) :: un
    INTEGER  :: stage
    un(:,:,1) = un_in
    DO stage = 2, ERK%s+1
       !CALL one_stage_ERK(stage,un)
       CALL one_stage_modified_ERK(stage,un)
    END DO
    un_in = un(:,:,ERK%s+1)
  END SUBROUTINE full_step_ERK

  SUBROUTINE one_stage_ERK(stage,un)
    USE shallow_water_functions
    USE IDP_limiting_shallow_water
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: stage
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,ERK%s+1) :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: rk, ulow, uhigh
    INTEGER :: comp, stage_prime
    REAL(KIND=8) :: time_stage
    INTEGER :: l
    stage_prime = ERK%lp_of_l(stage)
    time_stage = inputs%time+ERK%C(stage)*inputs%dt

    SELECT CASE(inputs%method_type)
    CASE('galerkin')
       RETURN
    CASE('viscous')
       velocity = compute_velocity(un(:,:,stage_prime))
       CALL compute_dijL(un(:,:,stage_prime))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijL)
       !===TEST
       CALL compute_bounds_and_low_flux(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL,lumped,diag,un(:,:,stage_prime))
       CALL sum_flux(FluxijL,rk)
       CALL divide_by_lumped(.true.,rk)
       un(:,:,stage) = un(:,:,stage_prime)+ERK%inc_C(stage)*inputs%dt*rk
       CALL bdy(un(:,:,stage),time_stage)!===Enforce boundary condition
       RETURN

    CASE('high')
       !===Best result obtained without smoothness indicator (both low order and high order)
       !===Low-order update
       velocity = compute_velocity(un(:,:,stage_prime))
       CALL compute_dijL(un(:,:,stage_prime))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijL)
       !CALL compute_bounds_and_low_flux(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL,&
       !     lumped,diag,un(:,:,stage_prime),opt_utest=uhigh)
       !===TEST
       CALL compute_bounds_and_low_flux(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL,&
            lumped,diag,un(:,:,stage_prime))
       DO comp = 1, inputs%syst_size
          fluxijL(comp)%aa = ERK%inc_C(stage)*fluxijL(comp)%aa
       END DO
       CALL sum_flux(FluxijL,rk)
       CALL divide_by_lumped(.true.,rk)
       ulow = un(:,:,stage_prime)+inputs%dt*rk
       !===TEST
       !WRITE(*,*) ' test', MAXVAL(ABS(ulow-uhigh))
       !STOP
       !===TEST

       !===High-order update
       CALL compute_inviscid_high_order_flux(un(:,:,stage-1),fluxij)
       DO comp = 1, inputs%syst_size
          fluxij_High(comp,stage-1)%aa = fluxij(comp)%aa
          fluxijH(comp)%aa=0.d0
          DO l = 1, stage-1
             fluxijH(comp)%aa = fluxijH(comp)%aa+ERK%MatRK(stage,l)*fluxij_High(comp,l)%aa
          END DO
       END DO
       CALL compute_muijL
       !CALL entropy_residual(un(:,:,stage-1)) !===Entropy residual at stage-1
       !CALL smoothness(un(:,1,stage_prime),dijL,dijH)
       CALL entropy_commutator(un(:,:,stage-1))
       
       CALL add_visc_to_high_flux(dijH,fluxijH,ERK%inc_C(stage),un(:,:,stage_prime))
       CALL sum_flux(fluxijH,rk)
       CALL divide_by_lumped(inputs%if_lumped,rk)
       un(:,:,stage) = un(:,:,stage_prime)+inputs%dt*rk !===Unlimited high-order solution
       IF (inputs%if_convex_limiting) THEN
          CALL convex_limiting_proc(velocity,un(:,:,stage_prime),ulow,un(:,:,stage),&
               FluxijH,FluxijL,mass,lumped,diag)
       END IF
       CALL bdy(un(:,:,stage),time_stage)!===Enforce boundary condition
       RETURN
    END SELECT
  END SUBROUTINE one_stage_ERK

    SUBROUTINE one_stage_modified_ERK(stage,un)
    USE shallow_water_functions
    USE IDP_limiting_shallow_water
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: stage
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,ERK%s+1) :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: rk, ulow, uhigh, SL, Source_shift
    INTEGER :: comp, stage_prime
    REAL(KIND=8) :: time_stage
    INTEGER :: l
    stage_prime = ERK%lp_of_l(stage)
    time_stage = inputs%time+ERK%C(stage)*inputs%dt
    SELECT CASE(inputs%method_type)
    CASE('galerkin')
       RETURN
    CASE('viscous')
       velocity = compute_velocity(un(:,:,stage_prime))
       CALL compute_dijL(un(:,:,stage_prime))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijL)
       !===TEST
       CALL compute_bounds_and_low_flux(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL,lumped,diag,un(:,:,stage_prime))
       CALL sum_flux(FluxijL,rk)
       CALL divide_by_lumped(.true.,rk)
       un(:,:,stage) = un(:,:,stage_prime)+ERK%inc_C(stage)*inputs%dt*rk

       CALL bdy(un(:,:,stage),time_stage)!===Enforce boundary condition
       RETURN

    CASE('test')
       velocity = compute_velocity(un(:,:,stage_prime))
       CALL compute_dijL(un(:,:,stage_prime))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijL)
       !===TEST
       CALL compute_bounds_and_low_flux_and_source(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL,SL, &
            lumped,diag,un(:,:,stage_prime))
       CALL sum_flux(FluxijL,rk)
       rk = rk + SL
       CALL divide_by_lumped(.true.,rk)
       un(:,:,stage) = un(:,:,stage_prime)+ERK%inc_C(stage)*inputs%dt*rk 

       CALL bdy(un(:,:,stage),time_stage)!===Enforce boundary condition
       RETURN

    CASE('high')
       !===Low-order update
       velocity = compute_velocity(un(:,:,stage_prime))
       CALL compute_dijL(un(:,:,stage_prime))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijL)
       !===TEST
       CALL compute_bounds_and_low_flux_and_source(ERK%inc_C(stage)*inputs%dt,velocity,dijL,cij,FluxijL, SL, &
            lumped,diag,un(:,:,stage_prime))
       DO comp = 1, inputs%syst_size
          fluxijL(comp)%aa = ERK%inc_C(stage)*fluxijL(comp)%aa
          SL(:,comp)=ERK%inc_C(stage)*SL(:,comp)
       END DO
       CALL sum_flux(FluxijL,rk)
       rk = rk + SL
       CALL divide_by_lumped(.true.,rk)
       ulow = un(:,:,stage_prime)+inputs%dt*rk 

       !===High-order update
       CALL compute_inviscid_high_order_flux_and_source(un(:,:,stage-1),fluxij, SiH(:,:,stage-1))
       Source_shift = 0.d0
       DO comp = 1, inputs%syst_size
          fluxij_High(comp,stage-1)%aa = fluxij(comp)%aa
          fluxijH(comp)%aa=0.d0
          DO l = 1, stage-1
             fluxijH(comp)%aa     = fluxijH(comp)%aa     + ERK%MatRK(stage,l)*fluxij_High(comp,l)%aa
             Source_shift(:,comp) = Source_shift(:,comp) + ERK%MatRK(stage,l)*SiH(:,comp,l)
          END DO
       END DO
       CALL entropy_residual(un(:,:,stage-1)) !===Entropy residual at stage-1
       !CALL entropy_commutator(un(:,:,stage-1))
       !===TEST
       !CALL smoothness(un(:,1,stage_prime),dijL,dijH)
       !===TEST
       CALL add_visc_to_high_flux(dijH,fluxijH,ERK%inc_C(stage),un(:,:,stage_prime))
       CALL sum_flux(fluxijH,rk)
       rk = rk + Source_shift
       CALL divide_by_lumped(inputs%if_lumped,rk)
       un(:,:,stage) = un(:,:,stage_prime)+inputs%dt*rk !===Unlimited high-order solution

       IF (inputs%if_convex_limiting) THEN
          Source_shift = Source_shift - SL
          CALL convex_limiting_proc(velocity,un(:,:,stage_prime),ulow,un(:,:,stage),&
               FluxijH,FluxijL,mass,lumped,diag)
          CALL divide_by_lumped(.true.,source_shift)
          un(:,:,stage) = un(:,:,stage) + inputs%dt*source_shift
       END IF
       CALL bdy(un(:,:,stage),time_stage)!===Enforce boundary condition
       RETURN
    END SELECT
  END SUBROUTINE one_stage_modified_ERK

  SUBROUTINE compute_inviscid_high_order_flux(un,fluxij)
    USE problem_setup
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
    TYPE(matrice_bloc), DIMENSION(:)                   :: fluxij
    REAL(KIND=8) :: xx
    INTEGER :: i,j, p, comp, d
    DO i = 1, mesh%np
       DO p = dijH%ia(i), dijH%ia(i+1) - 1
          j = dijH%ja(p)
          DO comp = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp))
             END DO
             fluxij(comp)%aa(p) = xx 
             IF (comp.NE.1) THEN !===Hydrostatic pressure
                !fluxij(comp)%aa(p) = fluxij(comp)%aa(p)   &
                !     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)**2 + un(i,1)**2)/2+un(i,1)*(bath(j)-bath(i)))
                !===Reference
                fluxij(comp)%aa(p) = fluxij(comp)%aa(p)   &
                     - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)*un(i,1)+un(i,1)*(bath(j)-bath(i)))
                !===Reference
             END IF
          END DO
       END DO
    END DO
  END SUBROUTINE compute_inviscid_high_order_flux

  SUBROUTINE compute_inviscid_high_order_flux_and_source(un,fluxij,Si)
    USE problem_setup
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un, Si
    TYPE(matrice_bloc), DIMENSION(:)                   :: fluxij
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: ss
    REAL(KIND=8) :: xx
    INTEGER :: i,j, p, comp, d
    DO i = 1, mesh%np
       ss =0.d0
       DO p = dijH%ia(i), dijH%ia(i+1) - 1
          j = dijH%ja(p)
          DO comp = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp))
             END DO
             fluxij(comp)%aa(p) = xx 
             IF (comp.NE.1) THEN !===Hydrostatic pressure
                !fluxij(comp)%aa(p) = fluxij(comp)%aa(p)   &
                !     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)**2 + un(i,1)**2)/2)
                !===Reference
                fluxij(comp)%aa(p) = fluxij(comp)%aa(p)   &
                     - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)*un(i,1))
                ss(comp) = ss(comp) - cij(comp-1)%aa(p)*inputs%gravity*(un(i,1)*(bath(j)-bath(i)))
                !===Reference
             END IF
          END DO
       END DO
       Si(i,:) = ss 
    END DO
  END SUBROUTINE compute_inviscid_high_order_flux_and_source

  SUBROUTINE sum_flux(Fluxij,rk)
    IMPLICIT NONE
    TYPE(matrice_bloc), DIMENSION(:) :: Fluxij
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: rk
    INTEGER :: i, ps, pe, l, comp
    DO i = 1, SIZE(dijL%ia)-1
       ps = dijL%ia(i)
       pe = dijL%ia(i+1)-1
       DO comp = 1, inputs%syst_size
          rk(i,comp) = SUM(fluxij(comp)%aa(ps:pe))
       END DO
    END DO
  END SUBROUTINE sum_flux

  SUBROUTINE add_visc_to_low_flux(dij,fluxij,scalar,uu)
    USE shallow_water_functions
    USE problem_setup
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: dij
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: uu
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)   :: fluxij
    REAL(KIND=8), DIMENSION(mesh%np)                  :: overh
    REAL(KIND=8) :: scalar, Hstarij, Hstarji, ratij, ratji
    INTEGER :: i, j, p, comp
    overh = compute_one_over_h(uu(:,1))
    DO i = 1, SIZE(dij%ia)-1
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          Hstarij = MAX(0.d0,uu(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,uu(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          DO comp = 1, inputs%syst_size
             fluxij(comp)%aa(p) = fluxij(comp)%aa(p)+scalar*dij%aa(p)*(uu(j,comp)*ratji-uu(i,comp)*ratij)
          END DO
       END DO
    END DO
  END SUBROUTINE add_visc_to_low_flux

  SUBROUTINE add_visc_to_high_flux(dij,fluxij,scalar,uu)
    USE shallow_water_functions
    USE problem_setup
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(IN) :: dij
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: uu
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size)   :: fluxij
    REAL(KIND=8), DIMENSION(mesh%np)                  :: overh
    REAL(KIND=8) :: scalar, Hstarij, Hstarji, ratij, ratji
    INTEGER :: i, j, p, comp
    overh = compute_one_over_h(uu(:,1))
    DO i = 1, SIZE(dij%ia)-1
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          Hstarij = MAX(0.d0,uu(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,uu(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          DO comp = 1, inputs%syst_size
             !===Reference
             fluxij(comp)%aa(p) = fluxij(comp)%aa(p)+scalar*dij%aa(p)*(uu(j,comp)*ratji-uu(i,comp)*ratij)
          END DO
       END DO
    END DO
  END SUBROUTINE add_visc_to_high_flux

  SUBROUTINE entropy_residual(un)
    USE mesh_handling
    USE pardiso_solve
    USE problem_setup
    USE shallow_water_functions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)       :: un
    REAL(KIND=8), DIMENSION(mesh%np)                        :: scal, res, ent, rescale, maxn, minn
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim) :: ff
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)       :: Entprime
    REAL(KIND=8), DIMENSION(mesh%np,k_dim)                  :: ent_flux
    REAL(KIND=8) :: xx, yy, aa, bb, small_ent, small_rescale
    INTEGER :: k, i, j, p

    !small_ent = inputs%gravity*inputs%htiny*inputs%htiny
    !small_ent = inputs%gravity*inputs%max_water_h*inputs%htiny
    small_ent = 1.d-2*inputs%gravity*inputs%max_water_h*inputs%max_water_h
    !===TEST
    !small_ent = 1.d-10*inputs%gravity*inputs%max_water_h*inputs%max_water_h
    !===TEST
    small_rescale = MAXVAL(ABS(cij(1)%aa))*small_ent
    
    scal = inputs%gravity*un(:,1)**2
    !scal = scal + inputs%gravity*un(:,1)*bath !===eta+ghz ***(entropy with bathymetry)
    DO k = 1, k_dim
       scal = scal + 0.5d0*velocity(:,k)*un(:,k+1)
    END DO
    ent =  scal - 0.5d0*inputs%gravity*un(:,1)**2 !===|v|^2 h/2 + g h^2/2

    !===TEST
    DO i = 1, mesh%np
       maxn(i) = MAXVAL(ent(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       minn(i) = MINVAL(ent(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO
    !===TEST

    DO k = 1, k_dim
       ent_flux(:,k) = velocity(:,k)*scal !===v (|v|^2 h/2 + g h^2); includes vhgz if entropy with bathymetry
    END DO

    Entprime(:,1) = inputs%gravity*un(:,1) !===-|v|^2/2 + g h
    !Entprime(:,1) = Entprime(:,1) + inputs%gravity*bath !=== + g z ***(entropy with bathymetry)
    DO k = 1, k_dim
       Entprime(:,1) = Entprime(:,1) -0.5d0*velocity(:,k)**2
       Entprime(:,k+1) = velocity(:,k)  !===v
    END DO

    ff = flux(un)
    DO k = 1, k_dim
       ff(:,k+1,k) = ff(:,k+1,k) + 0.5d0*inputs%gravity*un(:,1)**2
    END DO

    DO i = 1, mesh%np
       aa=0.d0
       bb=0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          !yy = inputs%gravity*un(1,i)*bath(j) !===g h z ***(entropy with bathymetry)
          DO k = 1, k_dim
             !bb = bb - cij(k)%aa(p)*Entprime(k+1,i)*yy !=== h_i v.grad(z) ***(entropy with bathymetry)
             aa = aa + cij(k)%aa(p)*ent_flux(j,k)
             bb = bb - cij(k)%aa(p)*SUM(Entprime(i,1:k_dim+1)*ff(j,1:k_dim+1,k))
          END DO
       END DO
       res(i) = aa + bb
       rescale(i) = ABS(aa)+ABS(bb)
    END DO

    res = ABS(res)/(rescale+small_rescale)

    !===TEST
    !rescale =MAX(ABS(maxn-minn)/2, small_rescale)
    !res = ABS(res)/rescale
    !===TEST

 
    !DO i = 1, mesh%np
    !   IF (un(i,1).LE.inputs%htiny) res(i) = 1.d0 
    !END DO

    !===Thresholding is bad without limiting, subcritical test
    !res = threshold(res)
    !===Thresholding is bad without limiting, subcritical test

    DO i = 1, mesh%np
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          IF (i==j) CYCLE
          !===Average is not robust wrt dry states in paraboloid test, but good for accuracy in subcritical test
          !dijH%aa(p) = dijL%aa(p)*(res(i)+res(j))/2
          !dijH%aa(p)  = dijL%aa(p)*max(res(i),res(j))
          !dijH%aa(p) = muijL%aa(p)*max(res(i),res(j))
          !dijH%aa(p) = muijL%aa(p)*(res(i)+res(j))/2
          !===!TEST
          dijH%aa(p) = dijL%aa(p)*min(1.d0,max(res(i),res(j)))
          !===TEST
       END DO
    END DO

    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),ABS(res),'res.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, ABS(res), 'res.plt')
       END SELECT
    END IF

  END SUBROUTINE entropy_residual

  FUNCTION threshold(x) RESULT(g)
    USE mesh_handling
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np)  :: x, z, t, zp, relu, f, g
    !REAL(KIND=8), PARAMETER :: x0 = 0.1d0, x1=SQRT(3.d0)*x0 !x0=0.05 good for P1 (qudratic threshold)
    REAL(KIND=8), PARAMETER :: x0 = 0.2d0 !x0=0.1 (cubic threshold)
    !integer :: i
    !do i = 1, mesh%np
    !   x(i) = (i-1.d0)/mesh%np
    !end do
!!$    !===Quadratic threshold
!!$    z = x-x0
!!$    zp = x-2*x0
!!$    relu = (zp+abs(zp))/2
!!$    f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
!!$    g = (f + 2*x0**3)/(4*x0**3)
!!$    CALL plot_1d(x,g,'threshold1.plt') 
!!$
!!$    !===Cubic threshold
    relu = ((x-2*x0)+abs(x-2*x0))/2
    t = x/(2*x0)
    g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)

    !CALL plot_1d(x,g,'threshold2.plt') 
    !stop
    RETURN
  END FUNCTION threshold

  SUBROUTINE compute_dijL(un)
    USE mesh_handling
    USE shallow_water_functions
    USE CSR_transpose
    USE problem_setup
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: norm_cij, lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij, vell, velr
    REAL(KIND=8), DIMENSION(inputs%syst_size)     :: ul, ur
    REAL(KIND=8), DIMENSION(mesh%np)             :: overh
    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji

    !===Viscosity using compute_lambda
    overh = compute_one_over_h(un(:,1))   !===Fix July 6, 2023
    DO i = 1, mesh%np
       vell = velocity(i,:)
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          ul = un(i,:)*ratij
          ur = un(j,:)*ratji
          velr = velocity(j,:)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p) !=== definition of cij is same as in the paper
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             CALL compute_lambda_vacc(ul,ur,vell,velr,nij,lambda)
             dijL%aa(p) = norm_cij*lambda
          ELSE
             dijL%aa(p) = 0.d0
          END IF
       END DO
    END DO
    CALL transpose_op(dijL,'max')

    DO i = 1, mesh%np
       dijL%aa(diag(i)) = -SUM(dijL%aa(dijL%ia(i):dijL%ia(i+1)-1))
    END DO
    RETURN
  END SUBROUTINE compute_dijL

  SUBROUTINE compute_muijL
    USE mesh_handling
    USE CSR_transpose
    IMPLICIT NONE
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8), DIMENSION(k_dim)                :: ur, ul
    !===Viscosity using speed only
    DO i = 1, mesh%np
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p)
             END DO
             ul = velocity(i,:)
             ur = velocity(j,:)
             !!lambda=MAX(MAX(-SUM(nij*ul),0.d0),MAX(SUM(nij*ur),0.d0))
             lambda=MAX(ABS(SUM(nij*ul)),ABS(SUM(nij*ur)))
             muijL%aa(p) = lambda
          ELSE
             muijL%aa(p) = 0.d0
          END IF
       END DO
    END DO
    CALL transpose_op(muijL,'max')
    DO i = 1, mesh%np
       muijL%aa(diag(i)) = -SUM(muijL%aa(dijL%ia(i):dijL%ia(i+1)-1))
    END DO
    RETURN
  END SUBROUTINE compute_muijL

  SUBROUTINE COMPUTE_DT(u0)
    USE input_data
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: u0
    CALL compute_dijL(u0)
    inputs%dt = 0.5d0*inputs%CFL/MAXVAL(ABS(dijL%aa(diag))/lumped)
    !===Rescale time step
    inputs%dt = inputs%dt*ERK%s
    !===Rescale time step
  END SUBROUTINE COMPUTE_DT

  SUBROUTINE smoothness(waterh,dij_in,dij_out)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np)       :: waterh
    TYPE(matrice_bloc)                     :: dij_in, dij_out
    REAL(KIND=8), DIMENSION(mesh%np)       :: alpha, num, denom
    INTEGER      :: i, j, p
    num = 0.d0
    denom = 0.d0
    DO i = 1, mesh%np
       DO p = stiff%ia(i), stiff%ia(i+1) - 1
          j = stiff%ja(p)
          IF (i==j) CYCLE
          num(i) = num(i)     + stiff%aa(p)*(waterh(i) - waterh(j))
          denom(i) = denom(i) + abs(stiff%aa(p))*ABS(abs(waterh(i) - waterh(j))+inputs%htiny)
       END DO
       alpha(i) = ABS(num(i)/denom(i))**2
    END DO
    DO i = 1, mesh%np
       DO p = stiff%ia(i), stiff%ia(i+1) - 1
          j = stiff%ja(p)
          IF (i==j) THEN
             dij_out%aa(p) =0.d0
          ELSE
             !dij_out%aa(p) = ((alpha(i)+alpha(j))/2)*dij_in%aa(p)
             dij_out%aa(p) = max(alpha(i),alpha(j))*dij_in%aa(p)
          END IF
       END DO
    END DO
    DO i = 1, mesh%np
       dij_out%aa(diag(i)) = -SUM(dij_out%aa(dij_out%ia(i):dij_out%ia(i+1)-1))
    END DO
  END SUBROUTINE smoothness

  SUBROUTINE entropy_commutator(un)
    USE mesh_handling
    USE pardiso_solve
    USE shallow_water_functions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)       :: un
    REAL(KIND=8), DIMENSION(mesh%np)                        :: scal, res, ent, rescale, minn, maxn
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim) :: ff
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)       :: Entprime
    REAL(KIND=8), DIMENSION(mesh%np,k_dim)                  :: ent_flux
    REAL(KIND=8) :: xx, yy, small_ent, small_rescale
    INTEGER :: k, i, j, p

    small_ent = 1.d-10*inputs%gravity*inputs%max_water_h*inputs%max_water_h
    small_rescale = MAXVAL(ABS(cij(1)%aa))*small_ent

    scal = inputs%gravity*un(:,1)**2
    !scal = scal + inputs%gravity*un(:,1)*bath !===eta+ghz ***(entropy with bathymetry)
    DO k = 1, k_dim
       scal = scal + 0.5d0*velocity(:,k)*un(:,k+1)
    END DO
    ent =  scal - 0.5d0*inputs%gravity*un(:,1)**2 !===|v|^2 h/2 + g h^2/2

    DO i = 1, mesh%np
       maxn(i) = MAXVAL(ent(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
       minn(i) = MINVAL(ent(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO

    DO k = 1, k_dim
       ent_flux(:,k) = velocity(:,k)*scal !===v (v|^2 h/2 + g h^2); includes vhgz if entropy with bathymetry
    END DO

    Entprime(:,1) = inputs%gravity*un(:,1) !===-|v|^2/2 + g h 
    !Entprime(:,1) = Entprime(:,1) + inputs%gravity*bath !=== + g z ***(entropy with bathymetry)
    DO k = 1, k_dim 
       Entprime(:,1) = Entprime(:,1) -0.5d0*velocity(:,k)**2 
       Entprime(:,k+1) = velocity(:,k)  !===v
    END DO

    ff = flux(un)
    DO k = 1, k_dim
       ff(:,k+1,k) = ff(:,k+1,k) + 0.5d0*inputs%gravity*un(:,1)**2
    END DO

    res  = 0.d0
    DO i = 1, mesh%np
       xx=0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          !yy = inputs%gravity*un(i,1)*bath(j) !===g h z ***(entropy with bathymetry)
          DO k = 1, k_dim
             !xx = xx - cij(k)%aa(p)*Entprime(i,k+1)*yy !=== h_i v.grad(z) ***(entropy with bathymetry)
             xx = xx + cij(k)%aa(p)*(ent_flux(j,k)-SUM(Entprime(i,:)*ff(j,:,k)))
          END DO
       END DO
       res(i) = res(i) + xx
    END DO

    !rescale =MAX(ABS(maxn-minn)/2, inputs%gravity*inputs%htiny**2)
    !rescale =MAX(ABS(maxn-minn)/2,inputs%epsilon_htiny*maxn)
    rescale =MAX(ABS(maxn-minn)/2, small_rescale)
    res = ABS(res)/rescale

    !res = threshold(res)

    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          dijH%aa(p) = dijL%aa(p)*min(1.d0,max(res(i),res(j)))
       END DO
    END DO
    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),ABS(res),'res.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, ABS(res), 'res.plt')
       END SELECT
    END IF

  END SUBROUTINE entropy_commutator








  !==========LEFTOVER

  SUBROUTINE euler(un,unext)
    USE mesh_handling
    USE shallow_water_functions
    USE fct
    USE pardiso_solve
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: rk

    velocity = compute_velocity(un)

    !===Compute first-order viscosity
    CALL compute_dijL(un)
    CALL smb_1(un,rk)
    CALL divide_by_lumped(.true.,rk)

    !===Compute First-Order solution
    unext = un+inputs%dt*rk

    IF (inputs%method_type=='viscous') THEN
       CALL check_Hmin(unext)
       RETURN
    END IF
  END SUBROUTINE euler

  SUBROUTINE smb_1(un,rk)
    USE mesh_handling
    USE shallow_water_functions
    USE problem_setup
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim)        :: vv
    REAL(KIND=8), DIMENSION(mesh%np)                               :: overh
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, Hstarij, Hstarji, ratij, ratji

    overh = compute_one_over_h(un(:,1))   !===Fix July 6, 2023

    vv=flux(un)
    rk=0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))

          !===Fix July 6, 2023
          ratij = Hstarij*overh(i)
          ratji = Hstarji*overh(j)
          !===End

!!$          IF (un(i,1).LE.inputs%htiny) THEN !Important for long-time WB
!!$             ratij=0.d0
!!$          ELSE 
!!$             ratij = Hstarij/un(i,1)
!!$          END IF
!!$          IF (un(j,1).LE.inputs%htiny) THEN !Important for long-time WB
!!$             ratji=0.d0
!!$          ELSE
!!$             ratji = Hstarji/un(j,1)
!!$          END IF

          DO k = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(j,k,d)*ratji + vv(i,k,d)*ratij)
             END DO
             rk(i,k) = rk(i,k) + xx + dijL%aa(p)*(un(j,k)*ratji-un(i,k)*ratij)
          END DO
          DO k = 1, k_dim
             rk(i,k+1) = rk(i,k+1) &
                  - 0.5d0*inputs%gravity*(Hstarji**2 - Hstarij**2)*cij(k)%aa(p)
          END DO
       END DO
    END DO

    IF (inputs%if_friction)  THEN
       DO k = 1, k_dim
          rk(:,k+1)  = rk(:,k+1) - lumped*friction(k,un)
       END DO
    END IF
  END SUBROUTINE smb_1

  SUBROUTINE bdy(uu,t)
    USE problem_setup
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: uu
    REAL(KIND=8), INTENT(IN) :: t
    IF (SIZE(h_js_D).NE.0)  uu(h_js_D,1)  = sol_anal(1,mesh%rr(:,h_js_D),t)
    IF (SIZE(ux_js_D).NE.0) uu(ux_js_D,2) = sol_anal(2,mesh%rr(:,ux_js_D),t)
    IF (k_dim==2) THEN
       IF (SIZE(ux_js_D).NE.0) uu(uy_js_D,3) = sol_anal(3,mesh%rr(:,uy_js_D),t)
    END IF
  END SUBROUTINE bdy

  SUBROUTINE check_Hmin(h)
    USE mesh_handling
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: h
    INTEGER :: i
    LOGICAL, DIMENSION(mesh%np) :: check
    check=.TRUE.
    check(js_D)=.FALSE.
    DO i = 1, mesh%np
       IF (.NOT.check(i)) CYCLE
       IF (h(i,1)<0.d0) THEN
          WRITE(*,*) 'Min h<0, STOP', h(1,i)
          velocity(ux_js_D,1) = 0.d0
          IF (k_dim==2) velocity(uy_js_D,2) = 0.d0
          WRITE(*,*) 'MAXVAL(vel)', MAXVAL(ABS(velocity(:,1))), MAXVAL(ABS(velocity(:,k_dim)))
          CALL plot_scalar_field(mesh%jj, mesh%rr, h(:,1), 'h.plt')
          CALL plot_scalar_field(mesh%jj, mesh%rr, velocity(:,1), 'vx.plt')
          if (k_dim==2) CALL plot_scalar_field(mesh%jj, mesh%rr, velocity(:,2), 'vy.plt')
          STOP
       END IF
    END DO
  END SUBROUTINE check_Hmin

!!$    SUBROUTINE compute_inviscid_high_order_modified_flux(un,fluxij,Si)
!!$    USE problem_setup
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un, Si
!!$    TYPE(matrice_bloc), DIMENSION(:)                   :: fluxij
!!$    REAL(KIND=8) :: xx
!!$    INTEGER :: i,j, p, comp, d
!!$    Si = 0.d0
!!$    DO i = 1, mesh%np
!!$       DO p = dijH%ia(i), dijH%ia(i+1) - 1
!!$          j = dijH%ja(p)
!!$          DO comp = 1, inputs%syst_size
!!$             xx = 0.d0
!!$             DO d = 1, k_dim
!!$                xx = xx - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp))
!!$             END DO
!!$             fluxij(comp)%aa(p) = xx 
!!$             IF (comp.NE.1) THEN !===Hydrostatic pressure
!!$               fluxij(comp)%aa(p) = fluxij(comp)%aa(p)   &
!!$                    - cij(comp-1)%aa(p)*inputs%gravity*(un(j,1)**2 + un(i,1)**2)/2
!!$               Si(i,comp) = Si(i,comp) - cij(comp-1)%aa(p)*inputs%gravity*un(i,1)*(bath(j)-bath(i)) 
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$  END SUBROUTINE compute_inviscid_high_order_modified_flux


!!$  SUBROUTINE compute_high_order_flux(un,fluxij)
!!$    USE problem_setup
!!$    USE shallow_water_functions
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
!!$    TYPE(matrice_bloc), DIMENSION(:)                   :: fluxij
!!$    INTEGER :: i,j, p, comp, d
!!$    REAL(KIND=8) ::  Hstarij, Hstarji, ratij, ratji
!!$    REAL(KIND=8), DIMENSION(mesh%np) :: overh
!!$    !TEST
!!$    !overh = compute_one_over_h(un(:,1))   !===Fix July 6, 2023
!!$    !TEST
!!$    DO i = 1, mesh%np
!!$       DO p = dijH%ia(i), dijH%ia(i+1) - 1
!!$          j = dijH%ja(p)
!!$          !TEST
!!$          !Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
!!$          !Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
!!$          !ratij = Hstarij*overh(i)
!!$          !ratji = Hstarji*overh(j)
!!$          !TEST
!!$          DO comp = 1, inputs%syst_size
!!$             DO d = 1, k_dim
!!$                fluxij(comp)%aa(p) = &
!!$                     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp)) &
!!$                     + dijH%aa(p)*(un(j,comp)-un(i,comp))
!!$                !fluxij(comp)%aa(p) = &
!!$                !     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp)) &
!!$                !    + muijH%aa(p)*(un(j,comp)-un(i,comp))
!!$                !fluxij(comp)%aa(p) = &
!!$                !     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp)) &
!!$                !     + dijH%aa(p)*(un(j,comp)*ratji-un(i,comp)*ratij)
!!$                !fluxij(comp)%aa(p) = &
!!$                !     - cij(d)%aa(p)*(velocity(j,d)*un(j,comp) + velocity(i,d)*un(i,comp)) &
!!$                !     + muijH%aa(p)*(un(j,comp)*ratji-un(i,comp)*ratij)
!!$             END DO
!!$             IF (comp.NE.1) THEN !===Hydrostatic pressure
!!$                fluxij(comp)%aa(p) = fluxij(comp)%aa(p) &
!!$                     - cij(comp-1)%aa(p)*inputs%gravity*((un(j,1)**2 - un(i,1)**2)/2+un(i,1)*(bath(j)-bath(i)))       
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$  END SUBROUTINE compute_high_order_flux

!!$  SUBROUTINE low_flux(un,fluxij)
!!$    USE mesh_handling
!!$    USE shallow_water_functions
!!$    USE problem_setup
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
!!$    TYPE(matrice_bloc), DIMENSION(:)                               :: Fluxij
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim)        :: vv
!!$    REAL(KIND=8), DIMENSION(mesh%np)                               :: overh
!!$    INTEGER :: d, i, j, comp, p
!!$    REAL(KIND=8) :: xx, Hstarij, Hstarji, ratij, ratji
!!$
!!$    overh = compute_one_over_h(un(:,1))   !===Fix July 6, 2023
!!$
!!$    vv=flux(un)
!!$    DO i = 1, mesh%np
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
!!$          Hstarji = MAX(0.d0,un(j,1)+bath(j)-MAX(bath(i),bath(j)))
!!$          !===Fix July 6, 2023
!!$          ratij = Hstarij*overh(i)
!!$          ratji = Hstarji*overh(j)
!!$          !===End
!!$          DO comp = 1, inputs%syst_size
!!$             xx = 0.d0
!!$             DO d = 1, k_dim
!!$                xx = xx - cij(d)%aa(p)*(vv(j,comp,d)*ratji + vv(i,comp,d)*ratij)
!!$             END DO
!!$             fluxij(comp)%aa(p) =  xx
!!$          END DO
!!$          !DO comp = 1, k_dim
!!$          !   fluxij(comp+1)%aa(p) = fluxij(comp+1)%aa(p) &
!!$          !        - 0.5d0*inputs%gravity*(Hstarji**2 - Hstarij**2)*cij(comp)%aa(p)
!!$          !END DO
!!$          DO comp = 1, k_dim
!!$             fluxij(comp+1)%aa(p) = fluxij(comp+1)%aa(p) &
!!$                  - 0.5d0*inputs%gravity*(Hstarji**2 + Hstarij**2)*cij(comp)%aa(p)
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    !===FIX ME. MUST CREATE SEPARATE SOURCE 
!!$    !IF (inputs%if_friction)  THEN
!!$    !   DO comp = 1, k_dim
!!$    !      fluxij(comp+1)%aa(diag) = fluxij(comp+1)%aa(diag) - lumped*friction(comp,un)
!!$    !   END DO
!!$    !END IF
!!$  END SUBROUTINE low_flux

!!$  SUBROUTINE low_sources(un,Si)
!!$    USE mesh_handling
!!$    USE shallow_water_functions
!!$    USE problem_setup
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: Si
!!$    REAL(KIND=8), DIMENSION(mesh%np)                   :: overh
!!$    INTEGER :: d, i, j, comp, p
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size) :: xx
!!$    REAL(KIND=8) :: Hstarij, ratij
!!$
!!$    overh = compute_one_over_h(un(:,1))
!!$    
!!$    Si = 0.d0
!!$    DO i = 1, mesh%np
!!$       xx = 0.d0
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          Hstarij = MAX(0.d0,un(i,1)+bath(i)-MAX(bath(i),bath(j)))
!!$          ratij = Hstarij*overh(i)
!!$          DO comp = 1, k_dim
!!$             xx(comp+1) = xx(comp+1) + inputs%gravity*(un(i,1)*ratij)**2*cij(comp)%aa(p)
!!$          END DO
!!$       END DO
!!$       Si(i,:) = Si(i,:) + xx 
!!$    END DO
!!$    
!!$    !===FIX ME. MUST CREATE SEPARATE SOURCE 
!!$    IF (inputs%if_friction)  THEN
!!$       DO comp = 1, k_dim
!!$          Si(:,comp+1) = Si(:,comp+1) - lumped*friction(comp,un)
!!$       END DO
!!$    END IF
!!$  END SUBROUTINE low_sources

!!$  SUBROUTINE high_flux(un,fluxij)
!!$    USE mesh_handling
!!$    USE shallow_water_functions
!!$    USE problem_setup
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
!!$    TYPE(matrice_bloc), DIMENSION(:)                               :: Fluxij
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,k_dim)        :: vv
!!$    INTEGER :: d, i, j, comp, p
!!$    REAL(KIND=8) :: xx
!!$
!!$    vv=flux(un)
!!$    DO i = 1, mesh%np
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          !===End
!!$          DO comp = 1, inputs%syst_size
!!$             xx = 0.d0
!!$             DO d = 1, k_dim
!!$                xx = xx - cij(d)%aa(p)*(vv(j,comp,d) + vv(i,comp,d))
!!$             END DO
!!$             fluxij(comp)%aa(p) =  xx
!!$          END DO
!!$          !DO comp = 1, k_dim 
!!$          !   fluxij(comp+1)%aa(p) = fluxij(comp+1)%aa(p) &
!!$          !        - inputs%gravity*un(j,1)*un(i,1)*cij(comp)%aa(p) !&
!!$                     !- inputs%gravity* bath(j)*un(i,1)*cij(comp)%aa(p)
!!$          !END DO
!!$          DO comp = 1, k_dim
!!$             fluxij(comp+1)%aa(p) = fluxij(comp+1)%aa(p) &
!!$                  - inputs%gravity*(un(j,1)**2/2 + un(i,1)**2/2)*cij(comp)%aa(p) !&
!!$                     !- inputs%gravity* bath(j)*un(i,1)*cij(comp)%aa(p)
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    !===FIX ME. MUST CREATE SEPARATE SOURCE 
!!$    !IF (inputs%if_friction)  THEN
!!$    !   DO comp = 1, k_dim
!!$    !      fluxij(comp+1)%aa(diag) = fluxij(comp+1)%aa(diag) - lumped*friction(comp,un)
!!$    !   END DO
!!$    !END IF
!!$  END SUBROUTINE high_flux
!!$
!!$  SUBROUTINE high_sources(un,Si)
!!$    USE mesh_handling
!!$    USE problem_setup
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
!!$    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: Si
!!$    INTEGER :: d, i, j, comp, p
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size) :: xx
!!$
!!$    Si = 0.d0
!!$    DO i = 1, mesh%np
!!$       xx = 0.d0
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          DO comp = 1, k_dim 
!!$             xx(comp+1) = xx(comp+1) - inputs%gravity*bath(j)*un(i,1)*cij(comp)%aa(p) !&
!!$                  !+ inputs%gravity*(un(j,1)**2-un(i,1)**2)*cij(comp)%aa(p)
!!$          END DO
!!$       END DO
!!$       Si(i,:) = Si(i,:) + xx 
!!$    END DO
!!$    
!!$    IF (inputs%if_friction)  THEN
!!$       DO comp = 1, k_dim
!!$          Si(:,comp+1) = Si(:,comp+1) - lumped*friction(comp,un)
!!$       END DO
!!$    END IF
!!$  END SUBROUTINE high_sources



END MODULE IDP_update_shallow_water
