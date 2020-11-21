cmake_minimum_required(VERSION 3.17)
project(MonteCarlo C)

set(CMAKE_C_STANDARD 11)

add_library(MonteCarlo SHARED library.c library.h random.c random.h main.c)