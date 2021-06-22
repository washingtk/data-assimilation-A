#!/bin/bash

#
# program file
#
PARAMETER_FILE="global_variable.f90"
SUBROUTINE_FILE01="model_time_step.f90"
SUBROUTINE_FILE02="other_subroutines_modules.f90"
SUBROUTINE_FILE03="data_assimilation_modules.f90"
MAIN_FILE="main_lyapunov.f90"

#
# compile
#
gfortran ${PARAMETER_FILE} \
         ${SUBROUTINE_FILE01} \
         ${SUBROUTINE_FILE02} \
         ${SUBROUTINE_FILE03} \
         ${MAIN_FILE}

#
# rm mod file
#
rm *.mod
