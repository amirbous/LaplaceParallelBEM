MODULE IO

CONTAINS


SUBROUTINE GET_INPUT_PROBLEMNAME(input_problem_name)


    CHARACTER(LEN=*),  PARAMETER :: DEFAULT_PROBLEM_NAME = "default_problem"

    CHARACTER(LEN=50), intent(inout):: input_problem_name
    INTEGER :: nargs

    nargs = command_argument_count()

    IF (nargs == 0) THEN
        PRINT *, "No command line arguments provided. Using default problem name: '", DEFAULT_PROBLEM_NAME, "'"
        input_problem_name = DEFAULT_PROBLEM_NAME
    ELSE
        CALL GET_COMMAND_ARGUMENT(1, input_problem_name)
        PRINT *, "Input problem name: ", TRIM(input_problem_name)
    END IF



END SUBROUTINE GET_INPUT_PROBLEMNAME



END MODULE IO