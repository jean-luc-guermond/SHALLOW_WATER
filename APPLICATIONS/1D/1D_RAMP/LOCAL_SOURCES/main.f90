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
  REAL(KIND=8)                              :: tps, to, q2, q3
  LOGICAL :: once=.TRUE.

  CALL IDP_start_shallow_water_sub

  ALLOCATE(rk(mesh%np,inputs%syst_size),un(mesh%np,inputs%syst_size),&
       ui(mesh%np,inputs%syst_size),uo(mesh%np,inputs%syst_size),hmovie(mesh%np),&
       uinit(mesh%np,inputs%syst_size))
  inputs%time = 0.d0

  CALL initialize(un) 
  uinit = un

  CALL plot_1d(mesh%rr(1,:), bath, 'bath.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,1), 'hinit.plt')
  CALL COMPUTE_DT(un)
  write(*,*) inputs%dt

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
     uo = un
     write(*,*) 'time ', inputs%time, inputs%dt
     CALL full_step_ERK(un)
     inputs%time = inputs%time + inputs%dt
     it_time = it_time + 1

     WRITE(11,*) inputs%time, MAXVAL(ABS(uo(:,1)-un(:,1)))
     IF (once) THEN
        tps = user_time()
        once=.FALSE.
     END IF
  END DO
  tps = user_time() - tps

  CALL write_restart(un)

  WRITE(*,*) ' total time ', tps, 'Time per time step and dof', tps/(it_time*mesh%np), it_max
  WRITE(*,*) ' End Mass ', SUM(un(:,1)*lumped)

  CALL plot_1d(mesh%rr(1,:), un(:,1)+bath, 'HplusZ.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,1), 'h.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,2), 'qx.plt')
  CALL plot_1d(mesh%rr(1,:), un(:,2)*compute_one_over_h(un(:,1)), 'vx.plt')

CONTAINS

  SUBROUTINE write_restart(un)                                                  
    IMPLICIT NONE                                              
    REAL(KIND=8), DIMENSION(:,:) :: un                                                         
    OPEN(unit = 10, file = 'restart.final.'//inputs%file_name, form = 'unformatted', status = 'unknown')
    WRITE(10) inputs%time, un                                                    
    WRITE(*,*) ' inputs%time at checkpoint', inputs%time                            
    CLOSE(10)                                                                       
  END SUBROUTINE write_restart


END PROGRAM shallow_water
