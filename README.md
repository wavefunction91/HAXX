[![Build Status](https://travis-ci.org/wavefunction91/HAXX.svg?branch=master)](https://travis-ci.org/wavefunction91/HAXX)


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
**S**ubroutines) provides a BLAS-like interface for matrices and vectors over
the quaternion numbers. As HBLAS depends solely on the HAXX scalar quaternion
infrastructure, there is no plan currently to release them separately. 
Currently, HBLAS provides an optimized (serial) software implementation of quaternion
matrix operations for AVX and AVX2 microarchitectures (see [arXiv:1903.05575](http://arxiv.org/abs/1903.05575) for details).


HAXX is currently a development code which has been hand tuned for a select few microarchitectures to demonstrate the
efficacy of such operations on modern computing platforms.
The default caching parameters shipped with HAXX are optimized for the  Intel(R) Xeon(R) CPU E5-2660 (Sandy Bridge) processor.
The API specification in HAXX is very flexible, but most of the flexibility is not directly user-facing. If there is interest
in exposing such functionality, please open a GitHub issue.

A primary goal of HBLAS is not only to provide a
convenient and efficient interface for quaternion-quaternion linear algebra, 
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
David Williams-Young (Computational Research Division / Lawrence Berkeley National Laboratory) <br />
E-Mail: dbwy at lbl dot gov
