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
     REAL(KIND=8) :: a, h0, hL, Bx
  END TYPE my_parameters
  TYPE(my_parameters), PUBLIC :: parameters
CONTAINS

  !===This function is used to define the bathymetry.
  FUNCTION bath_function(rr) RESULT(bath_func)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: x
    REAL(KIND=8), DIMENSION(SIZE(rr(1, :))) :: bath_func

    !===Define x coordinate
    x = rr(1, :)

    !===Define bathymetry
    bath_func = (parameters%h0/parameters%a**2)*(x - parameters%hL/2.d0)**2
  END FUNCTION bath_function


  FUNCTION friction(k,un) RESULT(rk)
     USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    INTEGER                                            :: k
    REAL(KIND=8), DIMENSION(mesh%np)  :: rk
    rk = - inputs%mannings*un(k+1,:)
  END FUNCTION  friction

  !===This function is used to define the initial conditions.
  !   The initial profile should be set for each system component of the system.
  !   For example:
  !     1D, shallow water equations: h, h*u
  !     2D, shallow water equations: h, h*u, h*v
  !     1D, (hyperbolic) Serre-Green-Naghdi: h, h*u, h*eta, h*w, h*beta
  !     2D, (hyperbolic) Serre-Green-Naghdi: h, h*u, h*v h*eta, h*w, h*beta
  !  For more details on how to initialize, see the paper:
  !     https://doi.org/10.1016/j.jcp.2021.110809
  FUNCTION sol_anal(k, rr, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr, 2)) :: vv, x, aux
    REAL(KIND=8) :: omega, eta

 

    !===Define x coordinate
    x = rr(1, :)

    omega = SQRT(8.d0*inputs%gravity*parameters%h0)/parameters%a
    eta = sqrt(omega**2-inputs%mannings**2)/2.d0

    aux = (parameters%h0/parameters%a**2)*(x - parameters%hL/2.d0)**2
    vv = max(parameters%h0 - aux + ((parameters%a*parameters%Bx)**2/(8.d0*inputs%gravity**2*parameters%h0))&
         *EXP(-inputs%mannings*t) &
         *(-eta*inputs%mannings*SIN(2*eta*t)+(inputs%mannings**2/4.d0 - eta**2)*COS(2*eta*t)) &
         - ((parameters%bx**2)/(4.d0*inputs%gravity))*EXP(-inputs%mannings*t) &
         -(parameters%Bx/inputs%gravity)*EXP(-inputs%mannings*t/2.d0) &
         *(eta*COS(eta*t)+(inputs%mannings/2.d0)*SIN(eta*t))*(x-parameters%hL/2),0.d0)
    SELECT CASE(k)
    CASE(1)
       RETURN
    CASE(2)
       vv =  vv * parameters%Bx * EXP(-inputs%mannings*t/2)*SIN(eta*t)
       RETURN
    CASE DEFAULT
       WRITE(*,*) 'BUG initial_profile'
       STOP
    END SELECT
    RETURN

  END FUNCTION sol_anal

  !===This is subroutine is used for initializing the problem. No need to modify
  SUBROUTINE initialize(un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(OUT) :: un
    INTEGER :: k

    !===Define constants
    inputs%gravity=9.81d0
    parameters%a =3000.d0
    parameters%h0=10.d0
    parameters%hL=10000.d0
    parameters%Bx=2.d0

    ALLOCATE (bath(mesh%np))

    !===Initialize bathymetry using bath_function
    bath = bath_function(mesh%rr)

    !===Initialization of un
    DO k = 1, inputs%syst_size
       un(k, :) = sol_anal(k, mesh%rr, inputs%time)
    END DO

    !===If we restart, we initialize time and solution from file
    IF (inputs%if_restart) THEN
       OPEN(unit=10, file='restart.'//inputs%file_name, form = 'unformatted', status = 'unknown')
       READ (10) inputs%time, un
    END IF

    !===Define max water height as pointer for everyone to read
    inputs%max_water_h = MAXVAL(un(1, :))

    !===Then define htiny. Here htiny = max_x(h(x,0))*epsilon
    inputs%htiny = inputs%epsilon_tiny * inputs%max_water_h

    !===Sanity check
    IF (inputs%max_water_h .LT. 1.d-14) THEN
       WRITE (*, *) ' --- BUG in initialize(): Forgot to initialize max_water_h !!! ', &
            inputs%max_water_h
       WRITE (*, *) ' --- Aborting! '
       STOP
    ELSE IF (MINVAL(un(1, :)) .LT. 0.d0) THEN
       WRITE (*, *) ' --- BUG in initialize(): Initial negative water depth !!! ', &
            MINVAL(un(1, :))
       WRITE (*, *) ' --- Aborting! '
       STOP
    END IF
  END SUBROUTINE initialize

END MODULE problem_setup
