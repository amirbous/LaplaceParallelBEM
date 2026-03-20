PROGRAM Electrostatic_fortran


USE IO
USE COMPUTE

IMPLICIT NONE


CHARACTER(LEN=50) :: input_problem_name


! Model parameters
TYPE :: MODEL 
    INTEGER :: N
    CHARACTER (50) :: problem_name
end type




CALL GET_INPUT_PROBLEMNAME(input_problem_name)

WRITE (*, *) "--------------------------------"
WRITE (*, *) "Electrostatic BEM Solver in Fortran"
WRITE (*, *) "Input problem name: "
WRITE (*, *) TRIM(input_problem_name)
WRITE (*, *) "--------------------------------"

END PROGRAM Electrostatic_fortran