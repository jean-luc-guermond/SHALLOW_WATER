MODULE IDP_start_shallow_water
  !USE boundary_conditions
  USE pardiso_solve
  USE Butcher_tableau
  USE space_dim
  USE mesh_handling
  USE input_data
  USE IDP_update_shallow_water
CONTAINS
  SUBROUTINE IDP_start_shallow_water_sub
    CALL read_my_data('data')
    inputs%syst_size = k_dim+1
    ERK%sv = inputs%RKsv
    CALL ERK%init
    CALL construct_mesh
    CALL IDP_construct_shallow_water_matrices
    !===Pardiso parameters
    isolve_shallow_water_pardiso = -1 !===
    CALL allocate_pardiso_parameters(1)
    pardiso_param(1)%mtype = 1    !===real and structurally symmetric matrix
    pardiso_param(1)%phase = 33   !===Direct solve
    pardiso_param(1)%parm(4) = 0  !===Direct solve
  END SUBROUTINE IDP_start_shallow_water_sub
END MODULE IDP_start_shallow_water
