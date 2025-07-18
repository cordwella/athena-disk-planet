#ifndef POTENTIAL_FUNCTIONS_HPP_
#define POTENTIAL_FUNCTIONS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// File added by Amelia J. Cordwell <ajc356@cam.ac.uk>
//========================================================================================
//! \file potentials.hpp
//! \brief defines class PotentialFunctionsK0
//!   Contains functions for caclulation of effective 3D potentials in 2D

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

namespace PotentialFunctions {
    // Evaluate Chebyshev polynomials
    Real chbevl(Real x, Real arra[], int n);

    // Special functions for caclulating averaged planetary potentials
    Real i0(Real x);
    Real i0e(Real x);
    Real i1(Real x);
    Real i1e(Real x);

    Real k0(Real x);
    Real k0e(Real x);
    Real k1(Real x);
    Real k1e(Real x);
} // namespace PotentialFunctions
#endif //POTENTIAL_FUNCTIONS_K0_HPP_