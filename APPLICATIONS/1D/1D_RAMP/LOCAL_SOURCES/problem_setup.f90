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
   !===For problem specific constants that are read in 'parameters.in' file
   TYPE my_parameters
      REAL(KIND=8) :: hL, hR, slope
   END TYPE my_parameters
   TYPE(my_parameters), PUBLIC :: parameters
CONTAINS

   !===This function is used to define the bathymetry.
   FUNCTION bath_function(rr) RESULT(bath_func)
      USE mesh_handling
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: x
      REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: bath_func

      !===Define x coordinate
      x = rr(1, :)

      !===Define bathymetry
      bath_func = MAX(parameters%slope *x , 0.d0)
    END FUNCTION bath_function
    
    FUNCTION friction(k,un) RESULT(rk)
      USE mesh_handling
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: un
      INTEGER                                            :: k
      REAL(KIND=8), DIMENSION(mesh%np)  :: rk
      rk = 0.d0
    END FUNCTION  friction

   FUNCTION sol_anal(k, rr, t) RESULT(vv)
      USE mesh_handling
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: k
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND=8), INTENT(IN) :: t
      REAL(KIND=8), DIMENSION(SIZE(rr, 2)) :: vv, x
      
      !===Define x coordinate
      x = rr(1, :)

      !===Define initial condition for each component
      ! IF (t .LE. 1.d-14) THEN
      !===
      SELECT CASE (k)
      CASE (1) ! water depth h
         WHERE (x .LE. 0.d0)
            vv = MAX(0.d0, parameters%hL - bath_function(rr))
         ELSE WHERE
            vv = MAX(0.d0, parameters%hR - bath_function(rr))
         END WHERE
      CASE (2) ! x momentum h * u
         vv = 0.d0
      END SELECT
    END FUNCTION sol_anal

   !===This is subroutine is used for initializing the problem. No need to modify
   SUBROUTINE initialize(un)
      USE mesh_handling
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: un
      INTEGER :: k
      ALLOCATE (bath(mesh%np))
      
      !===Define constants
      inputs%gravity=9.81d0
      parameters%hL = 2.d0
      parameters%hR = 1.95d0
      parameters%slope = 1.d0

      !===Initialize bathymetry using bath_function
      bath = bath_function(mesh%rr)

      !===Initialization of un
      DO k = 1, inputs%syst_size
         un(:,k) = sol_anal(k, mesh%rr, inputs%time)
      END DO

      !===If we restart, we initialize time and solution from file
      IF (inputs%if_restart) THEN
         OPEN(unit=10, file='restart.'//inputs%file_name, form = 'unformatted', status = 'unknown')
         READ (10) inputs%time, un
      END IF

      !===Define max water height as pointer for everyone to read
      inputs%max_water_h = MAXVAL(un(:,1))

      !===Then define htiny. Here htiny = max_x(h(x,0))*epsilon
      inputs%htiny = inputs%epsilon_tiny * inputs%max_water_h
      !un(:,1) = MAX(un(:,1) , inputs%htiny) !===Bad fix

      !===Sanity check
      IF (inputs%max_water_h .LT. 1.d-14) THEN
         WRITE (*, *) ' --- BUG in initialize(): Forgot to initialize max_water_h !!! ', &
            inputs%max_water_h
         WRITE (*, *) ' --- Aborting! '
         STOP
      ELSE IF (MINVAL(un(1, :)) .LT. 0.d0) THEN
         WRITE (*, *) ' --- BUG in initialize(): Initial negative water depth !!! ', &
            MINVAL(un(:,1))
         WRITE (*, *) ' --- Aborting! '
         STOP
      END IF
   END SUBROUTINE initialize
END MODULE problem_setup
