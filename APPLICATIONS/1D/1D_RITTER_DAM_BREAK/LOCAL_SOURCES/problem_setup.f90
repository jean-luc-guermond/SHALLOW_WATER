!===============================================================================
! This is the problem_setup module. To set up a new test case,
! you should:
!    1. Define the bathymetry in bath_function()
!    2. Define the parameters in read_in_parameters()
!    3. Set up the initial profile in initial_profile()
!===============================================================================
MODULE problem_setup
  USE space_dim
  USE input_data
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: bath
  !===For problem specific constants
  TYPE my_parameters
     REAL(KIND=8) ::  hL, x0, L
  END TYPE my_parameters
  TYPE(my_parameters), PUBLIC :: parameters
CONTAINS

   !===This is subroutine is used for initializing the problem. No need to modify
  SUBROUTINE initialize(un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: un
    INTEGER :: k

    !===Define constants
    inputs%gravity=9.81d0
    parameters%hL = 0.005d0
    parameters%x0 = 5.d0
    parameters%L = 10.d0
    

    ALLOCATE (bath(mesh%np))

    !===Initialize bathymetry using bath_function
    bath = bath_function(mesh%rr)

    !===Initialization of un
    DO k = 1, inputs%syst_size
       un(:, k) = sol_anal(k, mesh%rr, inputs%time)
    END DO

    !===If we restart, we initialize time and solution from file
    IF (inputs%if_restart) THEN
       OPEN(unit=10, file='restart.'//inputs%file_name, form = 'unformatted', status = 'unknown')
       READ (10) inputs%time, un
    END IF

    !===Define max water height as pointer for everyone to read
    inputs%max_water_h = MAXVAL(un(:, 1))

    !===Then define htiny. Here htiny = max_x(h(x,0))*epsilon
    inputs%htiny = inputs%epsilon_tiny * inputs%max_water_h

    !===Sanity check
    IF (inputs%max_water_h .LT. 1.d-14) THEN
       WRITE (*, *) ' --- BUG in initialize(): Forgot to initialize max_water_h !!! ', &
            inputs%max_water_h
       WRITE (*, *) ' --- Aborting! '
       STOP
    ELSE IF (MINVAL(un(:, 1)) .LT. 0.d0) THEN
       WRITE (*, *) ' --- BUG in initialize(): Initial negative water depth !!! ', &
            MINVAL(un(:, 1))
       WRITE (*, *) ' --- Aborting! '
       !STOP
    END IF
  END SUBROUTINE initialize


  !===This function is used to define the bathymetry.
  FUNCTION bath_function(rr) RESULT(bath_func)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: x
    REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: bath_func
    INTEGER :: i
    bath_func = 0.d0
  END FUNCTION bath_function

  FUNCTION friction(k,un) RESULT(rk)
     USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
    INTEGER                                            :: k
    REAL(KIND=8), DIMENSION(mesh%np)  :: rk
    rk = 0.d0
  END FUNCTION  friction

  !===This function is used to define the initial conditions.
  FUNCTION sol_anal(k, rr, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr, 2)) :: vv
    REAL(KIND=8) :: xa, xb
    INTEGER :: i
    SELECT CASE(k)
    CASE(1)
       IF (t<1.d-14) THEN
          DO i = 1, SIZE(rr,2)
             IF (rr(1,i)<parameters%x0) THEN
                vv(i) = parameters%hL
             ELSE
                vv(i) = 0.d0
             END IF
          END DO
          RETURN
       END IF
       xa = parameters%x0 - t*SQRT(inputs%gravity*parameters%hL)
       xb = parameters%x0 + 2*t*SQRT(inputs%gravity*parameters%hL)
       DO i = 1, SIZE(rr,2)
          IF (rr(1,i)<xa) THEN
             vv(i) = parameters%hL
          ELSE IF (rr(1,i) < xb) THEN
             vv(i) = (4/(9*inputs%gravity))*(SQRT(inputs%gravity*parameters%hL) - (rr(1,i)-parameters%x0)/(2*t))**2
          ELSE
             vv(i) = 0.d0
          END IF
       END DO
    CASE(2)
       IF (t.LE.1.d-14) THEN
          vv =0.d0
       ELSE
          xa = parameters%x0 - t*SQRT(inputs%gravity*parameters%hL)
          xb = parameters%x0 + 2*t*SQRT(inputs%gravity*parameters%hL)
          DO i = 1, SIZE(rr,2)
             IF (rr(1,i)<xa) THEN
                vv(i) = 0.d0
             ELSE IF (rr(1,i) < xb) THEN
                vv(i) = (2.d0/3)*((rr(1,i)-parameters%x0)/t + SQRT(inputs%gravity*parameters%hL))&
                     * (4/(9*inputs%gravity))*(SQRT(inputs%gravity*parameters%hL) - (rr(1,i)-parameters%x0)/(2*t))**2
             ELSE
                vv(i) = 0.d0
             END IF
          END DO
       END IF
    CASE DEFAULT
       vv = 0.d0
    END SELECT
    RETURN
  END FUNCTION sol_anal
   END MODULE problem_setup
