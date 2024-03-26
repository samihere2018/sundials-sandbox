cmake_minimum_required (VERSION 3.10)
project (SUNDIALS-Sandbox)

find_package(SUNDIALS REQUIRED)

# To add an executable:
# * use the "add_executable" command, that lists the executable
#   name, followed by all source files that contribute to that
#   executable.
# * use the "target_link_libraries" command, that repeats the
#   executable name, followed by the list of libraries (from
#   SUNDIALS and elsewhere) that the executable depends on.
#
# Multiple examples of this are below, that create targets for
# each executable.  Note that each of these link against a target
# provided by SUNDIALS (information on the available targets may
# be found in the SUNDIALS documentation,
# https://sundials.readthedocs.io/en/latest/sundials/Install_link.html#using-sundials-as-a-third-party-library-in-other-cmake-projects)
# as well as the standard math library, "libm".

add_executable(testMyNVector
  testMyNVector.c
  myNVector.c)
target_link_libraries(testMyNVector PUBLIC
  SUNDIALS::core
  m)

add_executable(testMySUNLinearSolver
  testMySUNLinearSolver.c
  mySUNLinearSolver.c
  myNVector.c)
target_link_libraries(testMySUNLinearSolver PUBLIC
  SUNDIALS::core
  m)

add_executable(ark_analytic_nonlin
  ark_analytic_nonlin.c)
target_link_libraries(ark_analytic_nonlin PUBLIC
  SUNDIALS::arkode
  m)
