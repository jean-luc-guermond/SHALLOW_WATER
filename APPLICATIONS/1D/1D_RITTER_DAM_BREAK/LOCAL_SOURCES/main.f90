PROGRAM shallow_water
  USE space_dim
  USE input_data
  USE mesh_handling
  USE IDP_update_shallow_water
  USE IDP_start_shallow_water
  USE shallow_water_functions
  USE problem_setup
  USE sub_plot
  USE fem_tn
  !USE mesh_interpolation
  USE timing_tools
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rk, un, ui, uo, uinit
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: hmovie
  INTEGER                                   :: it, it_max, it_time=0
  REAL(KIND=8)                              :: tps, to, q2, q3, err=0.d0, norm=0.d0
  LOGICAL :: once=.TRUE.

  CALL IDP_start_shallow_water_sub

  ALLOCATE(rk(mesh%np,inputs%syst_size),un(mesh%np,inputs%syst_size),&
       ui(mesh%np,inputs%syst_size),uo(mesh%np,inputs%syst_size),hmovie(mesh%np),&
       uinit(mesh%np,inputs%syst_size))
  inputs%time = 1.d0

  CALL initialize(un)
 
  uinit = un

  CALL plot_1d(mesh%rr(1,:), bath, 'bath.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(:,1), 'hinit.plt')
  !CALL plot_1d(mesh%rr(1,:), uinit(1,:)+bath, 'hpluszinit.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(:,2), 'qxinit.plt')
  CALL COMPUTE_DT(un)

  WRITE(*,*) ' Begin Mass ', SUM(un(:,1)*lumped)

  it_max = INT(inputs%Tfinal/inputs%dt)
  IF (it_max==0) THEN
     it_max = 1
     inputs%dt = inputs%Tfinal
  ELSE
     inputs%dt = inputs%Tfinal/it_max
  END IF

  !DO it = 1, 100!it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     CALL COMPUTE_DT(un)
     !write(*,*) 'time ', inputs%time, inputs%dt,' Mass1', SUM(un(:,1)*lumped)
     !IF (inputs%time + inputs%dt>inputs%Tfinal) THEN
     !   inputs%dt=inputs%Tfinal-inputs%time
     !END IF
     CALL full_step_ERK(un)
     inputs%time = inputs%time + inputs%dt
     it_time = it_time + 1

     !===Monitor convergence
     hmovie = sol_anal(2,mesh%rr,inputs%time)
     CALL ns_l1 (mesh, hmovie-un(:,2), q3)
     IF (once) THEN
        tps = user_time()
        once=.FALSE.
     END IF
  END DO
  tps = user_time() - tps

  CALL write_restart(un)

  WRITE(*,*) ' total time ', tps, 'Time per time step and dof', tps/(it_time*mesh%np), it_time
  WRITE(*,*) ' End Mass ', SUM(un(:,1)*lumped)
  CALL compute_errors

  CALL plot_1d(mesh%rr(1,:), un(:,1)+bath, 'HplusZ.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,1)-sol_anal(1,mesh%rr,inputs%time), 'errh.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,1), 'h.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,2), 'qx.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,2)*compute_one_over_h(un(:,1)), 'vx.plt')
  hmovie = sol_anal(2,mesh%rr,inputs%time)*compute_one_over_h(sol_anal(1,mesh%rr,inputs%time))
  CALL plot_1d(mesh%rr(1,:), sol_anal(1,mesh%rr,inputs%time), 'hexact.plt')
  CALL plot_1d(mesh%rr(1,:), hmovie, 'vx_exact.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,2)-sol_anal(2,mesh%rr,inputs%time), 'errqx.plt')

CONTAINS

  SUBROUTINE write_restart(un)                                                  
    IMPLICIT NONE                                              
    REAL(KIND=8), DIMENSION(:,:) :: un                                                         
    OPEN(unit = 10, file = 'restart.final.'//inputs%file_name, form = 'unformatted', status = 'unknown')
    WRITE(10) inputs%time, un                                                    
    WRITE(*,*) ' inputs%time at checkpoint', inputs%time                            
    CLOSE(10)                                                                       
  END SUBROUTINE write_restart

  SUBROUTINE compute_errors
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: hh
    REAL(KIND=8), DIMENSION(mesh%np)  :: uexact
    REAL(KIND=8) :: err, norm

    uexact = sol_anal(1,mesh%rr,inputs%time)

    CALL ns_l1(mesh, un(:,1)-uexact, err)
    CALL ns_l1(mesh, uexact, norm)
    WRITE(*,*) ' Relative L1 error on h    ', err/norm, err

    CALL ns_0(mesh, un(:,1)-uexact, err)
    CALL ns_0(mesh, uexact, norm)
    WRITE(*,*) ' Relative L2 error on h    ', err/norm, err

    hh = un(:,1)-sol_anal(1,mesh%rr,inputs%time)
    err = MAXVAL(ABS(hh))
    norm = MAXVAL(ABS(sol_anal(1,mesh%rr,inputs%time)))
    WRITE(*,*) ' Relative Linfty error on h', err/norm, err

    uexact = sol_anal(2,mesh%rr,inputs%time)
    CALL ns_l1(mesh, un(:,2)-uexact, err)
    CALL ns_l1(mesh, uexact, norm)

    WRITE(*,*) ' Relative L1 error on Q    ', (err)/(norm), (err)

    uexact = sol_anal(2,mesh%rr,inputs%time)
    CALL ns_0(mesh, un(:,2)- uexact, err)
    CALL ns_0(mesh, uexact, norm)
    WRITE(*,*) ' Relative L2 error on Q    ', (err)/(norm), (err)

    hh = ABS(un(:,2)-sol_anal(2,mesh%rr,inputs%time))
    err = MAXVAL(hh)
    norm = MAXVAL(ABS(uexact))
    WRITE(*,*) ' Relative Linfty error on Q', err/norm, err
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

END PROGRAM shallow_water
