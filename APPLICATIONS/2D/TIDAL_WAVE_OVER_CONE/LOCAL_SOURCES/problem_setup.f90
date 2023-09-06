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
     REAL(KIND=8) :: htop, rcone, scone, hcone, h0
  END TYPE my_parameters
  TYPE(my_parameters), PUBLIC :: parameters
CONTAINS

  !===This function is used to define the bathymetry.
  FUNCTION bath_function(rr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2))       :: vv
    REAL(KIND=8) :: radius
    INTEGER :: i
    vv = 0.d0
    DO i = 1, SIZE(rr,2)                                                    
       radius = SQRT((rr(1,i)-15)**2 + (rr(2,i)-13)**2)                
       IF (radius<parameters%rcone) THEN                                                 
          vv(i) = min(parameters%htop,parameters%hcone-radius/parameters%scone)                            
       ELSE                                                                      
          vv(i) = 0.d0                                                         
       END IF
    END DO
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
    REAL(KIND=8) :: vel, xs, alpha
    INTEGER :: i

    alpha = 1*parameters%h0                                                              
    xs = 15.0d0-12.96d0                                                       
    vel = SQRT(3.d0*alpha/(4*parameters%h0**3))                                          
    SELECT CASE(k)                                                            
    CASE(1)                                                                   
       !IF (t.LE.1.d-14) THEN                                                  
       !   DO i = 1, SIZE(rr,2)                                                
       !     vv(i) = MAX(parameters%h0+alpha/(COSH(vel*(rr(1,i)-xs))**2)-bath(i),0.d0)                                  
       !   END DO
       !ELSE                                                                   
       !   vv = parameters%h0                                                             
       !END IF
       DO i = 1, SIZE(rr,2)                                                
          vv(i) = MAX(parameters%h0-bath(i),0.d0)                                  
       END DO
    CASE(2)                                                                   
       !IF (t.LE.1.d-14) THEN                                                  
       !   DO i = 1, SIZE(rr,2)                                                
       !      vv(i) = alpha/(COSH(vel*(rr(1,i)-xs))**2)                      
       !      vv(i) = vv(i)*SQRT(inputs%gravity/parameters%h0)                            
       !   END DO
       !ELSE                                                                   
          vv = 0.d0                                                            
       !END IF
    CASE(3)                                                                   
       vv = 0.d0                                                              
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
    parameters%htop=0.625d0                                                              
    parameters%rcone=3.6d0                                                               
    parameters%scone=4.d0                                                                
    parameters%hcone=0.9d0                                                               
    parameters%h0=0.32  
 
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
