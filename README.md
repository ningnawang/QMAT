# QMAT-simplify

### Intro

This is an open-source code repository for SIGGRAPH paper [Q-MAT: Computing Medial Axis Transform by Quadratic Error Minimization](https://personal.utdallas.edu/~xguo/Q-MAT.pdf).

It contains the core *simplify* function of Q-MAT.

### Requirements

- C++ (>=14)

- CGAL (5.6.1)

  NOTE: Since CGAL version 5.0, CGAL is header-only be default, which means that there is **no need to build CGAL before it can be used**. Thus, CGAL header files are included in `include/CGAL`. Also, you can see `/Eigen` and `/boost` in the same directory. These are dependencies for CGAL and are also header-only. Details can be seen at: [CGAL 5.6.1 Manual](https://doc.cgal.org/latest/Manual/thirdparty.html)
  
- CMake (3.2.6)

