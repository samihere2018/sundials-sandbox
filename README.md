# sundials-sandbox
Template repostory for codes that utilize SUNDIALS.

This includes simple implementations of customized N_Vector and SUNLinearSolver modules, along with testing routines for each.  It also includes two existing SUNDIALS examples, to serve as templates for codes that utilize the ARKODE and CVODE integrators from SUNDIALS:

* Custom N_Vector module: `myNVector.h` and `myNVector.c`, with testing routine `testMyNVector.c`.

* Custom SUNLinearSolver module: `mySUNLinearSolver.h` and `mySUNLinearSolver.c`, with testing routine `testMySUNLinearSolver.c`.

* ARKODE example problem: `ark_analytic_nonlin.c`.  This uses ARKODE's "ERKStep" solver module, and the serial N_Vector module from SUNDIALS.

* CVODE example problem: `cvRoberts_dns.c`.  This uses CVODE's BDF solver, the serial N_Vector module, and the dense SUNMatrix and SUNLinearSolver modules from SUNDIALS.

This repository provides a simple `CMakeLists.txt` file to compile each of the above codes against a given SUNDIALS installation.

This setup assumes out-of-source builds, and requires that the path for the installed `SUNDIALSConfig.cmake` file be held in the CMake variable  `SUNDIALS_DIR`.  For example, from within this directory:

```
mkdir build
cd build
cmake -DSUNDIALS_DIR=/usr/local/sundials/lib/cmake/sundials ..
make
```

This has been tested using SUNDIALS v7.0.0.  Older SUNDIALS installations may be used by changing to a different branch:

* `main`: designed to work with SUNDIALS v7.0.0 (should work with any v7.x.x)

* `sundials-v6.6.1`: designed to work with SUNDIALS v6.6.1 (should work with any v6.x.x)
