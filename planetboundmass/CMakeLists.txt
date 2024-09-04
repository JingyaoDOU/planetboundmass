cmake_minimum_required(VERSION 3.4...3.18)
project(boundcpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Include pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.10.4  # Change this to the latest version if needed
)
FetchContent_MakeAvailable(pybind11)

# Build the module
pybind11_add_module(boundcpp binding.cpp bound.cpp)
