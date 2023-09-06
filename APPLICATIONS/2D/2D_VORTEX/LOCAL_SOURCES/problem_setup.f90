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
     REAL(KIND=8) :: Hinfty, r0, beta, x0(k_dim), vinfty(k_dim), a
  END TYPE my_parameters
  TYPE(my_parameters), PUBLIC :: parameters
CONTAINS

  !===This function is used to define the bathymetry.
  FUNCTION bath_function(rr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: vv, x2
    REAL(KIND=8) :: a2
    INTEGER :: i
    !===Define x coordinate
    DO i = 1, SIZE(rr,2)
       x2(i)=SUM(rr(:,i)**2)
    END DO
    a2 = MAXVAL(x2)

    !===Define bathymetry
    !vv = parameters%a*parameters%Hinfty*(x2/a2)
    vv = 0.d0
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
    REAL(KIND=8),                  INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr, 2))      :: vv
    REAL(KIND=8), DIMENSION(SIZE(rr, 2))      :: x2, psi, h
    REAL(KIND=8), DIMENSION(k_dim,SIZE(rr, 2)):: xx
    REAL(KIND=8), PARAMETER :: twopi=2.d0*ACOS(-1.d0)
    INTEGER :: i

    DO i = 1, SIZE(rr,2)
       xx(:,i) = rr(:,i) - parameters%x0 - parameters%vinfty*t
       x2(i) = SUM(xx(:,i)**2)/parameters%r0**2
    END DO
    psi = (parameters%beta/twopi)*EXP(0.5d0*(1-x2))
    
!!$    h = parameters%hinfty - psi
!!$    bath = -psi**2/(2*inputs%gravity*parameters%r0**2) - h
!!$    SELECT CASE(k)
!!$    CASE(1)
!!$       vv = h
!!$       RETURN
!!$    CASE(2)
!!$       vv =  (parameters%vinfty(1) - (xx(2,:)/parameters%r0**2)*psi)*h
!!$    CASE(3)
!!$       vv =  (parameters%vinfty(2) + (xx(1,:)/parameters%r0**2)*psi)*h
!!$       RETURN
!!$    CASE DEFAULT
!!$       WRITE(*,*) 'BUG initial_profile'
!!$       STOP
!!$    END SELECT
!!$
    h = parameters%hinfty - (0.5d0/(inputs%gravity*parameters%r0**2))*psi**2 - bath_function(rr)
    SELECT CASE(k)
    CASE(1)
       vv = h
       RETURN
    CASE(2)
       vv =  (parameters%vinfty(1) - (xx(2,:)/parameters%r0**2)*psi)*h
    CASE(3)
       vv =  (parameters%vinfty(2) + (xx(1,:)/parameters%r0**2)*psi)*h
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
    inputs%gravity= 9.81d0
    parameters%a= 0.5d0
    parameters%Hinfty = 1.d0
    parameters%r0 = 1.d0
    parameters%beta = 2.d0
    parameters%vinfty(1) = 1.d0/SQRT(2.d0)
    parameters%vinfty(2) = parameters%vinfty(1)
    DO k = 1, k_dim
       parameters%x0(k) = MINVAL(mesh%rr(k,:)) + (MAXVAL(mesh%rr(k,:)) - MINVAL(mesh%rr(k,:)))/2
    END DO

    !TEST
    !parameters%x0 = 0.d0
    !TEST

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

END MODULE problem_setup
