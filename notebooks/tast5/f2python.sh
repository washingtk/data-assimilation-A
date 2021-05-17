#!/bin/bash

name_in_py="fortran_model"
name_in_fort="model_time_step.f90"

f2py -c -m ${name_in_py} ${name_in_fort}
