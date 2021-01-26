cmake_minimum_required(VERSION 3.17)
project(EvolBio)

set(CMAKE_CXX_STANDARD 14)

add_executable(EvolBio main.cpp definitions.h generation.h specimen.h)