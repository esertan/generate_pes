cmake_minimum_required(VERSION 2.6)
project(makepes)
enable_language(Fortran)

set (CMAKE_VERBOSE_MAKEFILE OFF)
set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "-fopenmp -O0 -pg")

add_executable(makepes.lx main.f90 mydefs.f90 pes_init.f90 pes_calc.f90)

add_library(peslib mydefs.f90 pes_init.f90 pes_calc.f90)
