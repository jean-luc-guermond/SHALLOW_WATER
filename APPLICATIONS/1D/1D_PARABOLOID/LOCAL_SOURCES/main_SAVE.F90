PROGRAM shallow_water
  USE space_dim
  USE input_data
  USE mesh_handling
  USE IDP_update_shallow_water
  USE IDP_start_shallow_water
  USE problem_setup
  USE sub_plot
  USE fem_tn
  !USE mesh_interpolation
  USE timing_tools
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rk, un, ui, uo, uinit
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: hmovie
  INTEGER                                   :: it, it_max
  REAL(KIND=8)                              :: tps, to, q2, q3
  LOGICAL :: once=.TRUE.

  CALL IDP_start_shallow_water_sub

  ALLOCATE(rk(inputs%syst_size,mesh%np),un(inputs%syst_size,mesh%np),&
       ui(inputs%syst_size,mesh%np),uo(inputs%syst_size,mesh%np),hmovie(mesh%np),&
       uinit(inputs%syst_size,mesh%np))
  inputs%time = 0.d0
  
  CALL initialize(un) 
  uinit = un

  CALL plot_1d(mesh%rr(1,:), bath, 'bath.plt')
  !CALL plot_1d(mesh%rr(1,:), uinit(1,:), 'hinit.plt')
  !CALL plot_1d(mesh%rr(1,:), uinit(1,:)+bath, 'hpluszinit.plt')
  !CALL plot_1d(mesh%rr(1,:), uinit(2,:), 'qxinit.plt')
  CALL COMPUTE_DT(un)

  WRITE(*,*) ' Begin Mass ', SUM(un(1,:)*lumped)

  it_max = INT(inputs%Tfinal/inputs%dt)
  IF (it_max==0) THEN
     it_max = 1
     inputs%dt = inputs%Tfinal
  ELSE
     inputs%dt = inputs%Tfinal/it_max
  END IF

  !DO it = 1, it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     CALL COMPUTE_DT(un)
     !write(*,*) 'time ', inputs%time, inputs%dt,' Mass1', SUM(un(1,:)*lumped)

     to = inputs%time
     uo = un  !t
     !===Step 1
     CALL euler(uo,un) !t
     !write(*,*) 'time ', inputs%time, inputs%dt, ' Mass2', SUM(uo(1,:)*lumped), SUM(un(1,:)*lumped)
     if (ABS(SUM(uo(1,:)*lumped)-SUM(un(1,:)*lumped))/SUM(uo(1,:)*lumped).ge.1d-7) STOP
     CALL bdy(un,inputs%time+inputs%dt) !t+dt

     !===Step 2
     inputs%time=to+inputs%dt
     CALL euler(un,ui) !t+dt
     !write(*,*) 'time ', inputs%time, inputs%dt, ' Mass3', SUM(un(1,:)*lumped), SUM(ui(1,:)*lumped)
     un = (3*uo+ ui)/4
     CALL bdy(un,inputs%time+inputs%dt/2) !t+dt/2
     !write(*,*) 'time ', inputs%time, inputs%dt, ' Mass4', SUM(ui(1,:)*lumped), SUM(un(1,:)*lumped)
     
     !===Step 3
     inputs%time =  to + inputs%dt/2
     CALL euler(un,ui) !t+dt/2
    ! write(*,*) 'time ', inputs%time, inputs%dt, ' Mass5', SUM(un(1,:)*lumped), SUM(ui(1,:)*lumped)
     un = (uo+ 2*ui)/3
     CALL bdy(un,inputs%time+inputs%dt) !t+dt
     inputs%time = to + inputs%dt
     !write(*,*) 'time ', inputs%time, inputs%dt, ' Mass6', SUM(ui(1,:)*lumped), SUM(un(1,:)*lumped)
     
     !===Monitor convergence
     hmovie = sol_anal(1,mesh%rr,inputs%time)
     CALL ns_l1 (mesh, hmovie-un(1,:), q3)
     CALL ns_l1 (mesh, hmovie, q2)
     WRITE(11,*) 1/inputs%time, q3/q2
        
     IF (once) THEN
        tps = user_time()
        once=.FALSE.
     END IF
     
  END DO
  tps = user_time() - tps
  WRITE(*,*) ' total time ', tps, 'Time per time step and dof', tps/(it_max*mesh%np), it_max
  WRITE(*,*) ' End Mass ', SUM(un(1,:)*lumped)
  CALL compute_errors

  CALL plot_1d(mesh%rr(1,:), un(1,:)+bath, 'HplusZ.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:)-sol_anal(1,mesh%rr,inputs%time), 'errh.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:), 'h.plt')
  CALL plot_1d(mesh%rr(1,:), un(2,:), 'qx.plt')
  CALL plot_1d(mesh%rr(1,:), un(2,:)-sol_anal(2,mesh%rr,inputs%time), 'errqx.plt')

CONTAINS

  SUBROUTINE bdy(uu,t)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: uu
    REAL(KIND=8), INTENT(IN) :: t
    IF (SIZE(h_js_D).NE.0)  uu(1,h_js_D)  = sol_anal(1,mesh%rr(:,h_js_D),t)
    IF (SIZE(ux_js_D).NE.0) uu(2,ux_js_D) = sol_anal(2,mesh%rr(:,ux_js_D),t)
  END SUBROUTINE bdy


  SUBROUTINE compute_errors
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: hh
    REAL(KIND=8), DIMENSION(k_dim,mesh%gauss%l_G*mesh%me):: rr_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me)  :: uexact
    REAL(KIND=8), DIMENSION(mesh%np)                 :: zero
    REAL(KIND=8) :: err, norm

    CALL r_gauss(rr_gauss) 
    uexact = sol_anal(1,rr_gauss,inputs%time)
    zero = 0.d0
    CALL ns_anal_l1(mesh, un(1,:), uexact, err)
    CALL ns_anal_l1(mesh, zero, uexact, norm)
    WRITE(*,*) ' Relative L1 error on h (Gaussian)', err/norm, err
    !CALL ns_l1 (mesh, hh-un(1,:), err)
    !CALL ns_l1 (mesh, hh, norm)
    !WRITE(*,*) ' Relative L1 error on h', err/norm

    CALL ns_anal_0(mesh, un(1,:), uexact, err)
    CALL ns_anal_0(mesh, zero, uexact, norm)
    WRITE(*,*) ' Relative L2 error on h (Gaussian)', err/norm, err

    hh = un(1,:)-sol_anal(1,mesh%rr,inputs%time)
    err = MAXVAL(ABS(hh))
    norm = MAXVAL(ABS(sol_anal(1,mesh%rr,inputs%time)))
    WRITE(*,*) ' Relative Linfty error on h       ', err/norm, err

    uexact = sol_anal(2,rr_gauss,inputs%time)
    CALL ns_anal_l1(mesh, un(2,:), uexact, err)
    CALL ns_anal_l1(mesh, zero, uexact, norm)

    WRITE(*,*) ' Relative L1 error on Q (Gaussian)', (err)/(norm), (err)

    uexact = sol_anal(2,rr_gauss,inputs%time)
    CALL ns_anal_0(mesh, un(2,:), uexact, err)
    CALL ns_anal_0(mesh, zero, uexact, norm)
    WRITE(*,*) ' Relative L2 error on Q (Gaussian)', (err)/(norm), (err)

    hh = ABS(un(2,:)-sol_anal(2,mesh%rr,inputs%time))
    err = MAXVAL(hh)
    hh = ABS(sol_anal(2,mesh%rr,inputs%time))
    norm = MAXVAL(hh)
    WRITE(*,*) ' Relative Linfty error on Q       ', err/norm, err
  END SUBROUTINE compute_errors

  SUBROUTINE r_gauss(rr_gauss)
    USE Gauss_points
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: rr_gauss
    INTEGER :: index, k, l, m, n
    REAL(KIND=8) :: s 
    index = 0
    DO m = 1, mesh%me
       DO l = 1,  mesh%gauss%l_G; index = index + 1
          DO k = 1, mesh%gauss%k_d
             s = 0
             DO n=1, mesh%gauss%n_w
                s = s + mesh%gauss%ww(n,l)*mesh%rr(k,mesh%jj(n,m))
             END DO
             rr_gauss(k,index) = s
          END DO
       END DO
    END DO

  END SUBROUTINE r_gauss

END PROGRAM
