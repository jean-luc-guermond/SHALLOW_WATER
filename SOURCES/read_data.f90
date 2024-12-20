MODULE input_data
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE my_data
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     CHARACTER(len=200)             :: viscous_type
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
     INTEGER                        :: RKsv
     
     INTEGER                        :: syst_size
     REAL(KIND=8)                   :: Tfinal
     REAL(KIND=8)                   :: CFL
     CHARACTER(LEN=15)              :: method_type
     REAL(KIND=8)                   :: gravity
     REAL(KIND=8)                   :: mannings

     LOGICAL                        :: if_convex_limiting
     CHARACTER(LEN=30)              :: limiter_type
     LOGICAL                        :: if_relax_bounds
     REAL(KIND=8)                   :: dt, time
     
     INTEGER                        :: h_nb_Dir_bdy, ux_nb_Dir_bdy, uy_nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: h_Dir_list, ux_Dir_list, uy_Dir_list
     INTEGER                        :: nb_udotn_zero
     INTEGER                        :: nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: Dir_list
     INTEGER, DIMENSION(:), POINTER :: udotn_zero_list

     REAL(KIND=8)                   :: epsilon_tiny
     REAL(KIND=8)                   :: max_water_h
     REAL(KIND=8)                   :: htiny, hsmall

     LOGICAL                        :: if_lumped
     LOGICAL                        :: if_friction

     INTEGER                        :: nb_frame
     LOGICAL                        :: if_restart
  END type my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE 
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE character_strings
    USE space_dim
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    LOGICAL :: okay
    inputs%epsilon_tiny     = 1.d-10!==htiny
    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===") 
    READ (in_unit,*) inputs%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) inputs%file_name
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(in_unit,*) inputs%nb_dom
    ALLOCATE(inputs%list_dom(inputs%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(in_unit,*) inputs%list_dom
    CALL read_until(21, '===Type of finite element===')
    READ(in_unit,*) inputs%type_fe
    CALL read_until(in_unit, "===RK method type===")
    READ (in_unit,*) inputs%RKsv
    CALL read_until(in_unit, "===Final time===") 
    READ (in_unit,*) inputs%Tfinal
    CALL read_until(in_unit, "===CFL number===") 
    READ (in_unit,*) inputs%CFL
    CALL read_until(in_unit, "===Method type (galerkin, viscous, high)===")
    READ (in_unit,*) inputs%method_type
    SELECT CASE(inputs%method_type)
    CASE('high')
       CALL read_until(in_unit, "===Convex limiting (true/false)===")
       READ (in_unit,*) inputs%if_convex_limiting
       IF (inputs%if_convex_limiting) THEN
          CALL read_until(in_unit, "===Relax bounds?===")
          READ (in_unit,*)  inputs%if_relax_bounds
          IF (inputs%if_relax_bounds) THEN
             CALL read_until(in_unit, "===Limiter type (avg, minmod)===")
             READ (in_unit,*) inputs%limiter_type
          END IF
       END IF
    END SELECT
    CALL read_until(in_unit, "===Lumping the mass matrix? (true/false)===")
    READ (in_unit,*) inputs%if_lumped

    CALL find_string(in_unit, "===How many boundaries for u.n=0?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%nb_udotn_zero
       IF (inputs%nb_udotn_zero>0) THEN
          CALL read_until(in_unit, "===List of boundarie for u.n=0?===")
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
          READ (in_unit,*) inputs%udotn_zero_list
       ELSE
          inputs%nb_udotn_zero = 0
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
       END IF
    ELSE
       inputs%nb_udotn_zero = 0
       ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
    END IF
    CALL read_until(in_unit, "===How many Dirichlet boundaries for h?===")
    READ (in_unit,*)  inputs%h_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for h?===")
    ALLOCATE(inputs%h_Dir_list(inputs%h_nb_Dir_bdy))
    READ (in_unit,*) inputs%h_Dir_list
    CALL read_until(in_unit, "===How many Dirichlet boundaries for ux?===")
    READ (in_unit,*)  inputs%ux_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for ux?===")
    ALLOCATE(inputs%ux_Dir_list(inputs%ux_nb_Dir_bdy))
    READ (in_unit,*) inputs%ux_Dir_list
    IF (k_dim==2) THEN
       CALL read_until(in_unit, "===How many Dirichlet boundaries for uy?===")
       READ (in_unit,*)  inputs%uy_nb_Dir_bdy
       CALL read_until(in_unit, "===List of Dirichlet boundaries for uy?===")
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
       READ (in_unit,*) inputs%uy_Dir_list
    END IF

    CALL find_string(in_unit, "===Mannings Friction? (true/false)===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%if_friction
       IF (inputs%if_friction) THEN
          CALL read_until(in_unit, "===Mannings coefficient===")
          READ (in_unit,*) inputs%mannings
       END IF
    ELSE
       inputs%if_friction=.false.
       inputs%mannings = 0.d0
    END IF

    !===Postprocessing
    CALL find_string(in_unit, "===Restart?===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%if_restart
    ELSE
       inputs%if_restart=.false.
    END IF
    CALL find_string(in_unit, "===How many vtk frames?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%nb_frame
    ELSE
       inputs%nb_frame = 0
    END IF

 
 CLOSE(in_unit)
END SUBROUTINE read_my_data
END MODULE input_data


