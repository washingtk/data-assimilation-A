#!/bin/bash

#
# program file
#
PARAMETER_FILE="global_variable.f90"
SUBROUTINE_FILE01="model_time_step.f90"
SUBROUTINE_FILE02="other_subroutines_modules.f90"
SUBROUTINE_FILE03="data_assimilation_modules.f90"
# MAIN_FILE="main_kf.f90"
MAIN_FILE="main_enkf.f90"
# MAIN_FILE="main_model_ckeck.f90"
# MAIN_FILE="main_kf_beast_param.f90"
# MAIN_FILE="main_kfv02_best_param.f90"
# MAIN_FILE="main_lyapunov.f90"
# MAIN_FILE="main_3dvar_beast_param.f90"

#
# compile
#
gfortran ${PARAMETER_FILE} \
         ${SUBROUTINE_FILE01} \
         ${SUBROUTINE_FILE02} \
         ${SUBROUTINE_FILE03} \
         ${MAIN_FILE} \
         -llapack -lblas -O3

#
# rm mod file
#
rm *.mod
