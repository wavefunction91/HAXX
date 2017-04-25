https://travis-ci.org/wavefunction91/HAXX.svg?branch=master


Synopsis
========
**HAXX** (**H**amilton's Quaternion **A**lgebra for C**XX**) is a C++ software
infrastructure for the development of efficient scalar and tensorial quaternion
algorithms. HAXX can be thought of as two interdependent C++ software
libraries:

* The HAXX scalar quaternion class which handles the scalar operations
  (+,-,etc) over the quaternion numbers
* HBLAS for quaternion and mixed-type linear algebra


**HBLAS** (**H**amilton's Quaternion **B**asic **L**inear **A**lgebra
**S**ubroutines) provides a BLAS-like interface for matricies and vectors over
the quaternion numbers. As HBLAS depends solely on the HAXX scalar quaternion
infrastructure, there is no plan currently to release them separately. The
HBLAS implementaion of Level 1, 2, and 3 BLAS operations are currently based on
the original FORTRAN BLAS implementaions, with serveral extensions and
generalizations that account for the loss of scalar commutivity in the
quaternion numbers. A primary goal of HBLAS is not only to provide a
convienient and efficient interface for quaternion-quaternion linear algebra, 
but also to efficiently handle mixed-type (quaternion-real, quaternion-complex)
linear algebra through their natural embeddings into the quaternion algebra.


HAXX is actively being developed with little focus on backwards compatibility
with previous versions. The HAXX and HBLAS interfaces are constantly evolving
and can (will) change on a regular basis as new, exciting functionality is
added.

Design Goals
============
* A high-level, modern C++ API for scalar quaternion algebra (addition,
  subtraction, multiplication, division)
* Access to low level optimization and vectorization of the real arithmetic
  underlying quaternion operations
* Extension of BLAS functionality to quaternion algebra and mixed-type
  expressions (HBLAS)
* A reusable software framework to enable future scalar and tensorial
  algorithmic development using the quaternion algebra

Developers
==========
David Williams-Young (Li Research Group / University of Washington) <br />
E-Mail: dbwy at uw dot edu
