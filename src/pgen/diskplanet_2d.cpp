//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// File added by Amelia J. Cordwell <ajc356@cam.ac.uk>
//========================================================================================
//! \file diskplanet_alphabeta.cpp
//  \brief 2D disk with a planetary potential, beta cooling and alpha viscosity

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <memory>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"

#include "../utils/potential_functions.hpp"

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
void VelProfileCyl(const MeshBlock *pmb, const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3);

// User-defined boundary conditions for disk simulations
void DiskInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// User-defined planet source function
Real OmegaKep(const Real r);

// Holder fucntions for cell potentials
Real (*CalcPotentialCell)(const Real rad, const Real phi,
                   const Real radpl, const Real phipl);
void (*CalcPotentialDivergence)(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );

// Calculate fourth order potential
Real CalcPotentialCell4(const Real rad, const Real phi,
                        const Real radpl, const Real phipl);
void CalcPotential4Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );

Real CalcPotentialCell2(const Real rad, const Real phi,
                        const Real radpl, const Real phipl);
void CalcPotential2Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );


// Calculate bessel function type potential
Real CalcPotentialCellB(const Real rad, const Real phi,
                    const Real radpl, const Real phipl);
void CalcPotentialBDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );


Real CalcPotentialCellB_Hconst(const Real rad, const Real phi,
                    const Real radpl, const Real phipl);
void CalcPotentialB_HconstDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );

void CalcPotentialB_divHDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );

void CalcPotentialB_linear_mode_Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi );

void SourceTerms(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
void SourcePlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
void CoolingFunction(Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is, int ie);
void IsothermalReset(Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is, int ie);
void DampingZones(MeshBlock *pmb, Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is,  int ie);

void MyViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const  AthenaArray<Real> &w, const AthenaArray<Real> &bc,
                 int is, int ie, int js, int je, int ks, int ke);

void CustomMeshBeforeOutput(Mesh *pmesh, ParameterInput *pin);
void CollectBlockDataOutput(MeshBlock *pmb, ParameterInput *pin);
Real ViscCoeffCell(Real rad);

// problem parameters which are useful to make global to this file
// NB: tslope = initial temperature slope
// p0_over_rho0 = 0 value of pressure divided by zero value of temp
// written this way to allow for
Real gm0, r0, rho0, dslope, p0_over_rho0, tslope, gamma_gas, gmp, eps, scale_height_at_planet;
Real gmp_root_2_pi, smoothingB; // tiny bit quicker to evaluate the bessel function type poential
Real dfloor;
Real rout;
Real rmin, rmax;
Real gm1, igm1;

// Planet details
Real omegap, phipl0, rsm, rsm2, rsm4, tramp;
Real radpl, radpl2;
bool indirect_term;
int potential_order;

// Damping Zones
bool damping, dampdens, dampaverage, dampthermal;
// Real trampdamp; // to ramp up damping zones
Real tdamp, rdampi, rdampo; // parameters for inner damping zone to avoid reflections

// Cooling
bool strict_iso;
Real beta; // for Cooling term
Real height_slope;


// Viscosity
Real nuslope, nu_iso;
Real mdot, abs_mdot;

Real Omega0; // corotating frame's angular frequency
bool debug;

// Evaluate functions based on azimuthal averages
// This will slowdown computation significantly by requiring an additional set
// of MPI communications between nodes each timestep - however for some things
// this may be nessecary
bool evalaverage;
bool output_average;

Real mask_size;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Store azimuthal averages of v_phi, v_r, and sigma if required
  AllocateRealUserMeshDataField(10);

  for (int i=0; i<nreal_user_mesh_data_; ++i) {
    ruser_mesh_data[i].NewAthenaArray(mesh_size.nx1+2*(NGHOST));
  }
  // Setup additional output subroutine
  EnrollCustomApplyUserWorkBeforeOutput(CustomMeshBeforeOutput);

  std::stringstream msg;
  // Get parameters for gravitatonal potential of central point mass
  int orb_adv_order = pin->GetOrAddInteger("orbital_advection","OAorder",0);
  if (orb_adv_order != 0)
    std::cout << "Orbital advection is on!" << std::endl;
  else
    std::cout << "Orbital advection is off!" << std::endl;

  indirect_term = pin->GetOrAddBoolean("problem","indirect_term",false);
  std::cout << "indirect term = " << indirect_term <<std::endl;

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
  std::cout << "Omega0 = " << Omega0 <<std::endl;

  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  std::cout << "GM = " << gm0 <<std::endl;
  r0 = pin->GetOrAddReal("problem","r0",1.0);
  std::cout << "r0 = " << r0 <<std::endl;

  radpl = r0;
  radpl2 = r0 * r0;

  // Get parameters for planet
  gmp = pin->GetOrAddReal("problem","GMp",0.0);
  gmp_root_2_pi = gmp/(sqrt(2 * M_PI));
  phipl0 = pin->GetOrAddReal("problem","phipl0",M_PI);
  std::cout << "GMp = " << gmp <<std::endl;
  eps = pin->GetOrAddReal("problem","eps",0.6);
  std::cout << "eps = " << eps <<std::endl;
  omegap = sqrt((gm0+gmp)/r0/r0/r0);
  std::cout << "Omega_p = " << omegap <<std::endl;

  potential_order = pin->GetOrAddInteger("problem","potential_order", 4);
  std::cout << "Order of potential = " << potential_order <<std::endl;

  // How to calculate the planets orbit
  if (potential_order == 4) {
     CalcPotentialCell = &CalcPotentialCell4;
     CalcPotentialDivergence = &CalcPotential4Divergence;
  }else if (potential_order == 2) {
     CalcPotentialCell = &CalcPotentialCell2;
     CalcPotentialDivergence = &CalcPotential2Divergence;
  }  else if (potential_order == -1) {
     CalcPotentialCell = &CalcPotentialCellB;
     CalcPotentialDivergence = &CalcPotentialBDivergence;
  }  else if (potential_order == -2) {
     CalcPotentialCell = &CalcPotentialCellB_Hconst;
     CalcPotentialDivergence = &CalcPotentialB_HconstDivergence;
  }  else if (potential_order == -3) {
     CalcPotentialCell = &CalcPotentialCellB;
     CalcPotentialDivergence = &CalcPotentialB_divHDivergence;
  }  else if (potential_order == -4) {
     CalcPotentialCell = &CalcPotentialCellB;
     CalcPotentialDivergence = &CalcPotentialB_linear_mode_Divergence;
  } else {
    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]"
      << std::endl << "only support potential_order = 2, 4 or -1, -2, -3, -4 (bessel) at the moment!" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // smooth bessel functions by twice the cell size in the phi direction in the planet
  smoothingB = 2 * r0 * pin->GetReal("mesh", "x2max")/pin->GetInteger("mesh", "nx2");
  std::cout << "smoothingB = " << smoothingB <<std::endl;

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  std::cout << "rho0 = " << rho0 <<std::endl;
  dslope = pin->GetOrAddReal("problem","dslope",0.0);
  std::cout << "dslope = " << dslope <<std::endl;

  // Get parameters for damping zone
  damping = pin->GetOrAddBoolean("problem", "damping", true);
  dampdens = pin->GetOrAddBoolean("problem", "dampdens", true);
  dampthermal = pin->GetOrAddBoolean("problem", "dampthermal", false); // flag to also damp internal energy
  dampaverage = pin->GetOrAddBoolean("problem", "dampaverage", false);
  evalaverage = pin->GetOrAddBoolean("problem", "evalaverage", false);

  if (dampaverage){
    // Evaluations are required for average damping
    evalaverage = true;
  }
  // Whether to include f_wave in the outputs
  output_average = pin->GetOrAddBoolean("problem", "output_average", true);
  std::cout << "output_average = " << output_average << std::endl;

  tdamp = pin->GetOrAddReal("problem", "tdamp", 2.0);
  std::cout << "tdamp  = " << tdamp <<std::endl;
  std::cout << "damping  = " << damping <<std::endl;
  std::cout << "dampdens  = " << dampdens <<std::endl;
  std::cout << "dampaverage  = " << dampaverage <<std::endl;

  rdampi = pin->GetReal("problem", "rdampi");
  rdampo = pin->GetReal("problem", "rdampo");
  rmin  = pin->GetReal("mesh", "x1min");
  rmax  = pin->GetReal("mesh", "x1max");
  std::cout << "rmax = " << rmax << std::endl;

  Real rrat = pin->GetOrAddReal("mesh","x1rat",1.0);
  std::cout << "rrat = " << rrat << std::endl;

  Real rout_f = rrat*rmax;
  rout = (TWO_3RD)*(std::pow(rout_f, 3) - std::pow(rmax, 3)) /
           (std::pow(rout_f, 2) - std::pow(rmax, 2));
  std::cout << "rout = " << rout << std::endl;

  beta = pin->GetOrAddReal("problem","beta",0.0);
  std::cout << "beta = " << beta <<std::endl;

  // Get parameters for ramping up planet potential
  tramp = pin->GetOrAddReal("problem", "tramp", 10.0);

  tramp = 2.*M_PI*tramp;
  std::cout << "tramp (orbits) = " << tramp <<std::endl;

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_rho0 = pin->GetOrAddReal("problem","p0_over_rho0",0.0025);
    // default to global isothermal background
    tslope = pin->GetOrAddReal("problem","tslope", 0.0);
    gamma_gas = pin->GetReal("hydro","gamma");

    std::cout << "gamma_gas = " << gamma_gas <<std::endl;
    gm1 = (gamma_gas-1.);
    igm1 = 1./gm1;
  } else {
    p0_over_rho0=SQR(pin->GetReal("hydro","iso_sound_speed"));
    tslope = 0.0;
  }
  height_slope = 1.5 + tslope/2;
  mask_size = pin->GetOrAddReal("problem","mask_size", 1.0);
  std::cout << "mask_size = " << mask_size <<std::endl;

  std::cout << "tslope = " << tslope <<std::endl;
  std::cout << "p0_over_rho0 = " << p0_over_rho0 <<std::endl;
  dfloor=pin->GetOrAddReal("hydro","dfloor",1e-12);

  scale_height_at_planet = sqrt(p0_over_rho0)/omegap;
  mask_size = mask_size * scale_height_at_planet;

  rsm = eps * scale_height_at_planet;
  rsm2 = rsm*rsm;
  rsm4 = rsm2*rsm2;


  strict_iso = pin->GetOrAddBoolean("problem", "strict_iso", false);
  std::cout << "strict_iso = " << strict_iso << std::endl;

  // Get alpha viscosity from parameters
  Real alpha;
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  nu_iso = alpha * p0_over_rho0 / omegap;

  // temperature slope
  // slope in cs^2 = tslope
  // global isothermal causes this to cancel out
  // from nu = alpha * cs^2 / omega_k
  nuslope = 1.5 + tslope;
  std::cout << "nu_iso = " << nu_iso << std::endl;

  if (nu_iso > 0.0) {
    std::cout << "Enrolling viscosity " << nu_iso << " x R^" << nuslope <<std::endl;

    EnrollViscosityCoefficient(MyViscosity);
    // Assume initial VSS and calculate at the outer edge
    mdot = 3 * M_PI * nu_iso * ViscCoeffCell( rout ) * DenProfileCyl(rout,0.,0.);

    //if (dslope == -1.5) {
      // Use the first order correction to this value
      // Currently I've only derived it for this dslope
      // Real h_p = sqrt(p0_over_rho0);
      //Real v_r_first_order = 3 * nu_iso * (1 - 3 * p0_over_rho0 + 3 * p0_over_rho0 * p0_over_rho0) / (2 - 9 * p0_over_rho0 + 9 * p0_over_rho0 * p0_over_rho0);
      //mdot = v_r_first_order * 2 * M_PI;
    //}
  } else {
    std::cout << "No viscosity" << nuslope <<std::endl;
    mdot = 0.0;
  }
  std::cout << "nuslope = " << nuslope <<std::endl;
  std::cout << "mdot = " << mdot <<std::endl;

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiskInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiskOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  }

  // enroll user-defined source function
  EnrollUserExplicitSourceFunction(SourceTerms);
  std::cout << "enrollment success" <<std::endl;

  return;
}

// Allocate Output variables
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  std::cout << "mesh block init start" <<std::endl;

  int64_t &lx1=loc.lx1;
  int ioff = block_size.nx1*lx1;
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;

  // Allocate meshblock data for the caclulation of planet potential
  AllocateRealUserMeshBlockDataField(3 + NHYDRO);
  // Remember the momentum change due to planet potential at each location
  ruser_meshblock_data[0].NewAthenaArray(NHYDRO,block_size.nx3,
                                         block_size.nx2,block_size.nx1);

  // Save azimuthal averages of momentum and sigma into these
  for (int m=0; m<NHYDRO; ++m) { // for azimuthal averages of conserved quantities and high-cadence dts
    ruser_meshblock_data[1+m].NewAthenaArray(block_size.nx1);
  }
  ruser_meshblock_data[NHYDRO + 1].NewAthenaArray(NHYDRO,block_size.nx3,block_size.nx2,block_size.nx1);

  for (int i=is; i<=ie; i++) {
    int iglob = i+ioff;
    pmy_mesh->ruser_mesh_data[8](iglob) = pcoord->x1v(i);
  }

  // Allocate 2D output variables
  // TODO(amelia): Re-add vortensity calculation?
  AllocateUserX1OutputVariables(13);
  // Planet source terms
  SetUserX1OutputVariableName(0, "planet_src_mom1"); // i.e. dT
  SetUserX1OutputVariableName(1, "planet_src_mom2");
  SetUserX1OutputVariableName(2, "planet_src_energy");

  // Viscous flux terms
  SetUserX1OutputVariableName(3, "visc_flux_mom1"); // i.e. dG
  SetUserX1OutputVariableName(4, "visc_flux_mom2");
  SetUserX1OutputVariableName(5, "visc_flux_energy");

  // F wave and its derivative provided as an outputt
  SetUserX1OutputVariableName(6, "f_wave");
  SetUserX1OutputVariableName(7, "df_wave_dr");
  SetUserX1OutputVariableName(8, "dv_phi_dr");

  SetUserX1OutputVariableName(9, "sigma_excised");
  SetUserX1OutputVariableName(10, "sigma_excised_cell_c"); // cut off due to maximum
  SetUserX1OutputVariableName(11, "sigma_wedge_excised");
  SetUserX1OutputVariableName(12, "sigma_wedge_cell_c"); // cut off due to maximum

  std::cout << "mesh block init success" <<std::endl;

 return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::cout << "pgen start" <<std::endl;

  Real rad, phi, z;
  Real v1, v2, v3;
  AthenaArray<Real> vol(ncells1);
  AthenaArray<Real> *rumd = pmy_mesh->ruser_mesh_data;
  int64_t &lx1 = loc.lx1;
  int ioff = block_size.nx1*lx1;
  //int iglob = i+ioff-is;
 //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol);
//#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        int iglob = i+ioff;
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        phydro->u(IDN,k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(this,rad,phi,z,v1,v2,v3);
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;

        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);

          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                     + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

        rumd[6](iglob) = DenProfileCyl(rad,phi,z);
        rumd[4](iglob) = v1;
        rumd[5](iglob) = v2;

      }
    }
  }
  std::cout << "pgen success" <<std::endl;

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//========================================================================================
void MeshBlock::UserWorkInLoop(void) {

  if (evalaverage) {
    // Get sums for this mesh block
    for (int m=0; m<NHYDRO; ++m) {
      #pragma omp simd
      for (int i=is; i<=ie; ++i) {
        ruser_meshblock_data[1+m](i) = 0.0;                           // reset new array
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        #pragma omp simd
        for (int i=is; i<=ie; ++i) {
          ruser_meshblock_data[1](i) += phydro->u(IDN,k,j,i);
          // Put the sum of momenta into each of these arrays
          ruser_meshblock_data[2](i) += phydro->u(IM1,k,j,i);
          ruser_meshblock_data[3](i) += phydro->u(IM2,k,j,i);
          ruser_meshblock_data[4](i) += phydro->u(IM3,k,j,i);
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(void)
//  \brief Function called after every time step for user-defined work.
//========================================================================================
void Mesh::UserWorkInLoop(void) {
  if (evalaverage){

    MeshBlock *pmb = my_blocks(0);
    #pragma omp simd
    for (int i=0; i<mesh_size.nx1+2*(NGHOST); ++i) {
      ruser_mesh_data[0](i) = 0.0;
      ruser_mesh_data[1](i) = 0.0;
      ruser_mesh_data[2](i) = 0.0;
      ruser_mesh_data[3](i) = 0.0;
    }

    // Sum over meshblocks for each radial location
    for (int j=0; j<nblocal; ++j) {
      pmb = my_blocks(j);
      int64_t &lx1=pmb->loc.lx1;
      int ioff = pmb->block_size.nx1*lx1;
      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Global value
        int iglob = i + ioff;
        ruser_mesh_data[0](iglob) += pmb->ruser_meshblock_data[1](i); // global <Sigma>
        // Sum user meshblock data into a clobal mesh
        ruser_mesh_data[1](iglob) += pmb->ruser_meshblock_data[2](i); // global <Sigma u_r>
        ruser_mesh_data[2](iglob) += pmb->ruser_meshblock_data[3](i); // global <Sigma u_phi>
      }
    }

  #ifdef MPI_PARALLEL
    Real *dens_point = ruser_mesh_data[0].data();
    Real *u_r_point = ruser_mesh_data[1].data();
    Real *u_phi_point = ruser_mesh_data[2].data();
    // Now need to sum the <Sigma>, <Sigma u_R>, <Sigma u_phi> over meshblocks that share phi.
    // This will need MPI communication
    if (Globals::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, dens_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, u_r_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      // Sum over the different mesh datas for sigma u_r
      MPI_Reduce(MPI_IN_PLACE, u_phi_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
    } else {
      MPI_Reduce(dens_point, dens_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
      MPI_Reduce(u_r_point, u_r_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
      MPI_Reduce(u_phi_point, u_phi_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
    }

    // On node 0, broadcast each of these arrays to other values
    MPI_Bcast(dens_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_r_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_phi_point, mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  #endif

    // Create the average velocitiess
    for (int i=0; i<mesh_size.nx1+2*NGHOST; ++i) {
      ruser_mesh_data[4](i) = ruser_mesh_data[1](i)/ruser_mesh_data[0](i); // <Sigma u_R>/<Sigma>
      ruser_mesh_data[5](i) = ruser_mesh_data[2](i)/ruser_mesh_data[0](i); // <Sigma u_phi>/<Sigma>
      ruser_mesh_data[6](i) = ruser_mesh_data[0](i)/mesh_size.nx2; // <Sigma> average
    }
 }
  return;
}

// ========================================================================================
// ! \fn void CustomMeshBeforeOutput(Mesh *pmesh, ParameterInput *pin)
//  \brief Function called before every output for user-defined work.
// =======================================================================================
void CustomMeshBeforeOutput(Mesh *pmesh, ParameterInput *pin){
  if (output_average) {
    for (int i=0; i< pmesh->nblocal; ++i)
      CollectBlockDataOutput(pmesh->my_blocks(i), pin);

    // Share data accross processors
    // This is a repeat if eval average is already true
    // but this is still the most efficent approach

    MeshBlock *pmb = pmesh->my_blocks(0);
    #pragma omp simd
    for (int i=0; i< pmesh->mesh_size.nx1+2*(NGHOST); ++i) {
      pmesh->ruser_mesh_data[0](i) = 0.0;
      pmesh->ruser_mesh_data[1](i) = 0.0;
      pmesh->ruser_mesh_data[2](i) = 0.0;
      pmesh->ruser_mesh_data[3](i) = 0.0;
    }

    // Sum over meshblocks for each radial location
    for (int j=0; j< pmesh->nblocal; ++j) {
      pmb = pmesh->my_blocks(j);
      int64_t &lx1= pmb->loc.lx1;
      int ioff = pmb->block_size.nx1*lx1;

      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Global value
        int iglob = i + ioff;
        pmesh->ruser_mesh_data[0](iglob) += pmb->ruser_meshblock_data[1](i); // global <Sigma>
        // Sum user meshblock data into a clobal mesh
        pmesh->ruser_mesh_data[1](iglob) += pmb->ruser_meshblock_data[2](i); // global <Sigma u_r>
        pmesh->ruser_mesh_data[2](iglob) += pmb->ruser_meshblock_data[3](i); // global <Sigma u_phi>
        // correctly save R
        pmesh->ruser_mesh_data[8](iglob) = pmb->pcoord->x1v(i);
      }
    }

  #ifdef MPI_PARALLEL
    Real *dens_point = pmesh->ruser_mesh_data[0].data();
    Real *u_r_point = pmesh->ruser_mesh_data[1].data();
    Real *u_phi_point = pmesh->ruser_mesh_data[2].data();
    // Now need to sum the <Sigma>, <Sigma u_R>, <Sigma u_phi> over meshblocks that share phi.
    // This will need MPI communication
    if (Globals::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, dens_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, u_r_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      // Sum over the different mesh datas for sigma u_r
      MPI_Reduce(MPI_IN_PLACE, u_phi_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM, 0,
                 MPI_COMM_WORLD);
    } else {
      MPI_Reduce(dens_point, dens_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
      MPI_Reduce(u_r_point, u_r_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
      MPI_Reduce(u_phi_point, u_phi_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, MPI_SUM,
                 0, MPI_COMM_WORLD);
    }

    // On node 0, broadcast each of these arrays to other values
    MPI_Bcast(dens_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_r_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_phi_point, pmesh->mesh_size.nx1+2*(NGHOST), MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
  #endif

    // Create the average velocitiess
    for (int i=0; i< pmesh->mesh_size.nx1+2*NGHOST; ++i) {
      pmesh->ruser_mesh_data[4](i) = pmesh->ruser_mesh_data[1](i)/pmesh->ruser_mesh_data[0](i); // <Sigma u_R>/<Sigma>
      pmesh->ruser_mesh_data[5](i) = pmesh->ruser_mesh_data[2](i)/pmesh->ruser_mesh_data[0](i); // <Sigma u_phi>/<Sigma>
      pmesh->ruser_mesh_data[6](i) = pmesh->ruser_mesh_data[0](i)/pmesh->mesh_size.nx2; // <Sigma> average

    }

    // calculate d v_{phi}/dr
    for (int i=1; i< pmesh->mesh_size.nx1+NGHOST-1; ++i) {
      pmesh->ruser_mesh_data[7](i) = (
        pmesh->ruser_mesh_data[5](i+1) - pmesh->ruser_mesh_data[5](i-1))/(
        pmesh->ruser_mesh_data[8](i+1) - pmesh->ruser_mesh_data[8](i-1)); // d/dr <Sigma u_phi>/<Sigma>
    }
  }


  // Original Function parts to call
  for (int i=0; i<pmesh->nblocal; ++i)
    pmesh->my_blocks(i)->UserWorkBeforeOutput(pin);

}

void CollectBlockDataOutput(MeshBlock *pmb, ParameterInput *pin){
  // Get data on a meshblock level to add to the mesh before output

  for (int m=0; m<5; ++m) {
    #pragma omp simd
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      pmb->ruser_meshblock_data[1+m](i) = 0.0;                           // reset new array
    }
  }
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        pmb->ruser_meshblock_data[1](i) += pmb->phydro->u(IDN,k,j,i);
        // Put the sum of momenta into each of these arrays
        pmb->ruser_meshblock_data[2](i) += pmb->phydro->u(IM1,k,j,i);
        pmb->ruser_meshblock_data[3](i) += pmb->phydro->u(IM2,k,j,i);
        pmb->ruser_meshblock_data[4](i) += pmb->phydro->u(IM3,k,j,i);
      }
    }
  }
}


// ========================================================================================
// ! \fn void MeshBlock::UserWorkBeforeOutput(void)
//  \brief Function called before every output for user-defined work.
// =======================================================================================
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // Create additional output varaibles (if requested)
  int64_t &lx1=loc.lx1;
  int ioff = block_size.nx1*lx1;


  // Monitor fluxes along in x1-direction to later sum them up
  // reset variables
  for (int u=0; u<nuser_x1out_var; ++u) {
    #pragma omp simd
    for(int i=is; i<=ie; ++i) {
      user_x1out_var(u,0,0,i) = 0.0;
    }
  }

  // Calculate the effect of the planet (e.g. dT/dR)
  Real time = pmy_mesh->time;
  Real phipl = phipl0 + time*omegap;

  phipl = std::fmod(phipl, M_PI * 2);

  for(int k=ks; k<=ke; ++k) {
    for(int j=js; j<=je; ++j) {
      #pragma omp simd
      for(int i=is; i<=ie; ++i) {
        user_x1out_var(0,0,0,i) += ruser_meshblock_data[0](IM1,k,j,i); // save planet src mom1
        user_x1out_var(1,0,0,i) += ruser_meshblock_data[0](IM2,k,j,i); // save planet src mom2
        if(NON_BAROTROPIC_EOS) {
          user_x1out_var(2,0,0,i) += ruser_meshblock_data[0](IEN,k,j,i); // save planet src en
        }
      	for(int nm=0; nm<=4; ++nm){
      	    // reset data for flux dirvergences
      	    ruser_meshblock_data[NHYDRO + 1](nm, k, j, i) = 0;
            ruser_meshblock_data[0](nm, k, j, i) = 0;
      	}

        // add excised region from calcultion

        Real radius = pcoord->x1v(i);
        Real phi = pcoord->x2v(j);

        // check if outside the wedge
        // if (std::abs(radius - radpl) > mask_size or std::abs(phi - phipl) * radpl > mask_size) {
        if (std::abs(radius - radpl) > mask_size or std::abs(phi - phipl) * radpl > mask_size) {
          user_x1out_var(9,0,0,i) += phydro->u(IDN, k, j, i);
          user_x1out_var(10,0,0,i) += 1;

          if (std::abs(phi - phipl) * radpl > mask_size) {
            // wedge
            user_x1out_var(11,0,0,i)  += phydro->u(IDN, k, j, i);
            user_x1out_var(12,0,0,i) += 1;
          }
        }
      }
    }
  }

  // Calculate viscous fluxes
  // e.g. dG/dR

  // This calculates a few of the different total visocus effcts
  if (phydro->hdif.hydro_diffusion_defined) { // monitor viscous fluxes
    for(int k=ks; k<=ke; ++k) {
      for(int j=js; j<=je; ++j) {
        #pragma omp simd
        for(int i=is; i<=ie; ++i) {
          user_x1out_var(3,0,0,i) += phydro->hdif.visflx[X1DIR](IM1,k,j,i); // save visc src mom1
          user_x1out_var(4,0,0,i) += phydro->hdif.visflx[X1DIR](IM2,k,j,i); // save visc src mom2
          if(NON_BAROTROPIC_EOS) {
            user_x1out_var(5,0,0,i) += phydro->hdif.visflx[X1DIR](IEN,k,j,i); // save visc src en
          }
        }
      }
    }
  }


  if (output_average) {
    Real dphi = (2 * M_PI)/pmy_mesh->mesh_size.nx2;
    // get divergences
    // place flux divergences into an unuses slot of meshblock data
    // ensure that this block gets overriden rather than added to
    phydro->AddFluxDivergence(-1.0, ruser_meshblock_data[NHYDRO + 1]);

    for(int k=ks; k<=ke; ++k) {
      for(int j=js; j<=je; ++j) {
       #pragma omp simd
         for(int i=is; i<=ie; ++i) {
          // global i
          int iglob = i+ioff;

          // get fluxes
          Real radius = pcoord->x1v(i);

          // or really the
          Real f_wave_component    = phydro->flux[X1DIR](IM2,k,j,i) - pmy_mesh->ruser_mesh_data[5](iglob) * phydro->flux[X1DIR](IDN,k,j,i);
          // f_wave
          user_x1out_var(6,0,0,i) += radius* radius * dphi * f_wave_component;

          user_x1out_var(7,0,0,i) += radius * dphi * (
            f_wave_component + radius * (
              // divergence of phi - momentum flux
              ruser_meshblock_data[NHYDRO + 1](IM2, k, j, i)
              - pmy_mesh->ruser_mesh_data[7](iglob) * phydro->flux[X1DIR](IDN,k,j,i) // d v_phi/dr * mass flux in v_r
              - pmy_mesh->ruser_mesh_data[5](iglob) * ruser_meshblock_data[NHYDRO + 1](IDN, k, j, i)  // v_phi * divergence of mass flux
              ));

            // output d v_phi/d r
            user_x1out_var(8,0,0,i) +=  dphi * pmy_mesh->ruser_mesh_data[7](iglob);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//!\f transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (COORDINATE_SYSTEM == "cylindrical") {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    rad=fabs(pco->x1v(i)*sin(pco->x2v(j)));
    phi=pco->x3v(i);
    z=pco->x1v(i)*cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \f  computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real denmid = rho0*pow(rad/r0, dslope);
  return std::max(denmid,dfloor);
}

//----------------------------------------------------------------------------------------
//! \f  computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_rho0*pow(rad/r0, tslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! \f  computes rotational velocity in cylindrical coordinates

void VelProfileCyl(const MeshBlock *pmb, const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3) {
  Real p_over_r = PoverR(rad, phi, z);
  Real den = DenProfileCyl(rad,phi,z);
  Real vel = 0.0;

  if (NON_BAROTROPIC_EOS) {
    // Specifying a 2D disc - therefore ignore z dependence
    vel = (dslope+tslope)*p_over_r/(gm0/rad) + 1.0;
  } else {
    vel = dslope*p0_over_rho0/(gm0/rad) + 1.0;
  }

  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  if (pmb->porb->orbital_advection_defined)
    vel -= pmb->porb->OrbitalVelocity(pmb->porb, rad, phi, z);

  Real vR_visc = 0.0;

  if (pmb->phydro->hdif.nu_iso != 0) {
    vR_visc = -3 * (2 + tslope + dslope) * nu_iso * ViscCoeffCell(rad) / rad;
  }

  v1=vR_visc;
  v2=vel;
  v3=0.0;
  return;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: sets solution in ghost zones to initial values
//

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        GetCylCoord(pco,rad,phi,z,il-i,j,k);
        VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);

        prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
        prim(IM1,k,j,il-i) = v1;
        prim(IM2,k,j,il-i) = v2;
        prim(IM3,k,j,il-i) = v3;

        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }

}

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          prim(IM1,k,j,iu+i) = v1;
          prim(IM2,k,j,iu+i) = v2;
          prim(IM3,k,j,iu+i) = v3;
          if (NON_BAROTROPIC_EOS)
            prim(IPR,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
      }
    }
  }

}

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,jl-j,k);
        prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);
        prim(IM1,k,jl-j,i) = v1;
        prim(IM2,k,jl-j,i) = v2;
        prim(IM3,k,jl-j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
      }
    }
  }
}

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,ju+j,k);
        prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);
        prim(IM1,k,ju+j,i) = v1;
        prim(IM2,k,ju+j,i) = v2;
        prim(IM3,k,ju+j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
      }
    }
  }
}

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,kl-k);
        prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);
        prim(IM1,kl-k,j,i) = v1;
        prim(IM2,kl-k,j,i) = v2;
        prim(IM3,kl-k,j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
      }
    }
  }
}

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real v1(0.0), v2(0.0), v3(0.0);
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,ku+k);
        prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
        VelProfileCyl(pmb,rad,phi,z,v1,v2,v3);
        prim(IM1,ku+k,j,i) = v1;
        prim(IM2,ku+k,j,i) = v2;
        prim(IM3,ku+k,j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// \fn SourceTerms
//  \brief Adds source terms due to point mass planet, cooling and damping
//----------------------------------------------------------------------------------------

void SourceTerms(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

    // planet
    SourcePlanet(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);

    if (beta != 0.0 && NON_BAROTROPIC_EOS) {
      CoolingFunction(pmb->pcoord,dt,cons,prim,pmb->ks,pmb->ke,pmb->js,pmb->je,pmb->is,pmb->ie);
    } else if (strict_iso && NON_BAROTROPIC_EOS){
      // Reset local isothermality
      IsothermalReset(pmb->pcoord,dt,cons,prim,pmb->ks,pmb->ke,pmb->js,pmb->je,pmb->is,pmb->ie);
    }
    if (damping) {
     DampingZones(pmb, pmb->pcoord,dt,cons,prim,pmb->ks,pmb->ke,pmb->js,pmb->je,pmb->is,pmb->ie);
    }
    return;
}

//----------------------------------------------------------------------------------------
//! \fn void SourcePlanet
//  \brief Adds source terms due to point mass planet in a circular orbit
//----------------------------------------------------------------------------------------

void SourcePlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {
  // Planet coordinates
  Real phipl = phipl0 + time*omegap;

  Real ramp = 0.;
  // Ramp up potential over time
  if (time >= tramp){
    ramp = 1.0;
  } else {
    ramp = 0.5*(1.0-cos(M_PI*time/tramp));
  }

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // cell centres
        Real rad, phi, z;

        rad = pmb->pcoord->x1v(i);
        phi = pmb->pcoord->x2v(j); // central cell location in azimuth

        Real dPhidr, dPhidphi;

        (*CalcPotentialDivergence)(rad, phi, radpl, phipl,
                                   dPhidr, dPhidphi);
        if (indirect_term){
          dPhidr += gmp*cos(phi - phipl)/radpl2;
          dPhidphi -= gmp*rad*sin(phi - phipl)/radpl2;
        }

        Real den = prim(IDN,k,j,i);

        Real src1 = ramp*den*dPhidr;
        Real src2 = ramp*den*dPhidphi/rad;

        cons(IM1,k,j,i) -= src1 * dt;
        cons(IM2,k,j,i) -= src2 * dt;


        // Add to monitoring array
        pmb->ruser_meshblock_data[0](IM1,k,j,i) = -src1;
        pmb->ruser_meshblock_data[0](IM2,k,j,i) = -src2;

        if (NON_BAROTROPIC_EOS) {
          // dE/dt = F \dot v
          // Real work = -1 * dt * (src1 * cons(IM1,k,j,i)/prim(IDN,k,j,i) + src2 * cons(IM2,k,j,i)/prim(IDN,k,j,i));
          Real work = -1 * dt * (src1 * prim(IM1,k,j,i) + src2 * prim(IM2,k,j,i));
          pmb->ruser_meshblock_data[0](IEN,k,j,i) = work;

          cons(IEN,k,j,i) += work;

        }
      }
    }
  }
  return;
}

// Different potential functions
Real CalcPotentialCell4(const Real rad, const Real phi, const Real radpl,
                        const Real phipl) {
  Real pot;
  // Squared distance planet-cell and vector product r . rpl for indirect term
  Real dr2;
  dr2 = rad*rad-2.*rad*radpl*cos(phi-phipl)+radpl*radpl;
  pot = -gmp*( dr2+1.5*rsm2 )/pow( (dr2 + rsm2), 1.5);
  return pot;
}


void CalcPotential4Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi ) {

    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2;
    Real dr2prsm2 = dr2+rsm2;
    Real i32 = 1. / dr2prsm2 / sqrt(dr2prsm2);
    Real smth = i32 * (1.5 * (dr2 + 1.5 * rsm2) / dr2prsm2 - 1.); // Smoothing function

    dPhidr = gmp*2. * (rad - radpl*cos(phi-phipl)) * smth;
    dPhidphi = gmp*2.*rad*radpl*sin(phi-phipl) * smth;
}

Real CalcPotentialCell2(const Real rad, const Real phi, const Real radpl,
                        const Real phipl) {
    Real pot;
    // Squared distance planet-cell and vector product r . rpl for indirect term
    Real dr2;
    dr2 = rad*rad-2.*rad*radpl*cos(phi-phipl)+radpl*radpl;
    pot = -gmp/sqrt(dr2 + rsm2);
    return pot;
}


void CalcPotential2Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidr, Real& dPhidphi ) {

    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2;
    Real dr2prsm2 = dr2+rsm2;
    Real i32 = 1. / dr2prsm2 / sqrt(dr2prsm2);

    dPhidr   = gmp * i32 * (rad - radpl * cos(phi - phipl));
    dPhidphi = gmp * i32 * rad * radpl * sin(phi-phipl);
}


// Different planetary potential functions
Real CalcPotentialCellB(const Real rad, const Real phi, const Real radpl,
                        const Real phipl) {
  Real pot;
  Real dr2;
  dr2 = rad*rad-2.*rad*radpl*cos(phi-phipl)+radpl*radpl + smoothingB * smoothingB;

  // TODO: Move scale height at planet to actual local scale height
  // test including the derivative of the scale height?
  Real local_scale_height = scale_height_at_planet * std::pow(rad, height_slope);

  pot = - gmp_root_2_pi/local_scale_height * PotentialFunctions::k0e(
    dr2/(4 * local_scale_height * local_scale_height ));
  return pot;
}


void CalcPotentialBDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidrad, Real& dPhidphi ) {

    Real local_scale_height = scale_height_at_planet * std::pow(rad, height_slope);
    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2 + smoothingB * smoothingB;
    Real s2 = dr2/(4 * local_scale_height * local_scale_height);

    Real dPhids2 = - gmp_root_2_pi/local_scale_height * (
      PotentialFunctions::k0e(s2) - PotentialFunctions::k1e(s2));

    dPhidrad = dPhids2 * (rad - radpl * cos(phi - phipl))/(2 * local_scale_height * local_scale_height);
    dPhidphi = dPhids2 * rad * radpl  * sin(phi - phipl)/(2 * local_scale_height * local_scale_height);

}


// Different planetary potential functions
Real CalcPotentialCellB_Hconst(const Real rad, const Real phi, const Real radpl,
                        const Real phipl) {
  Real pot;
  Real dr2;
  dr2 = rad*rad-2.*rad*radpl*cos(phi-phipl)+radpl*radpl + smoothingB * smoothingB;

  pot = - gmp_root_2_pi/scale_height_at_planet * PotentialFunctions::k0e(
    dr2/(4 * scale_height_at_planet * scale_height_at_planet ));
  return pot;
}


void CalcPotentialB_HconstDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidrad, Real& dPhidphi ) {

    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2 + smoothingB * smoothingB;
    Real s2 = dr2/(4 * scale_height_at_planet * scale_height_at_planet);

    Real dPhids2 = - gmp_root_2_pi/scale_height_at_planet * (
      PotentialFunctions::k0e(s2) - PotentialFunctions::k1e(s2));

    dPhidrad = dPhids2 * (rad - radpl * cos(phi - phipl))/(2 * scale_height_at_planet * scale_height_at_planet);
    dPhidphi = dPhids2 * rad * radpl  * sin(phi - phipl)/(2 * scale_height_at_planet * scale_height_at_planet);

}


void CalcPotentialB_divHDivergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidrad, Real& dPhidphi ) {

    // Smoothing function here is just a cell length
    Real local_scale_height = scale_height_at_planet * std::pow(rad, height_slope);
    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2 + smoothingB * smoothingB;
    Real s2 = dr2/(4 * local_scale_height * local_scale_height);

    // TODO: Move scale height at planet to actual local scale height
    Real k_0 = PotentialFunctions::k0e(s2);
    Real k_1 = PotentialFunctions::k1e(s2);

    Real scale = - gmp_root_2_pi/local_scale_height;

    Real ds2dr = (
      (1 - height_slope) * rad - height_slope * radpl * radpl / rad
      - (1 - 2 * height_slope) * cos(phi - phipl))/(2 * local_scale_height * local_scale_height);

    dPhidrad   = scale * (-1 * height_slope * k_0/rad + ds2dr * (k_0 - k_1));
    dPhidphi = scale * (k_0 - k_1) * rad * radpl  * sin(phi - phipl)/(2 * local_scale_height * local_scale_height);
}


void CalcPotentialB_linear_mode_Divergence(
  const Real rad, const Real phi, const Real radpl,
  const Real phipl, Real& dPhidrad, Real& dPhidphi ) {

    // Use the additional potential components derived from the linear
    // mode decomposition process, i.e. \parital \Phi_0/\partial R + 2 \Phi_2 * \mu/R

    // Smoothing function here is just a cell length
    Real local_scale_height = scale_height_at_planet * std::pow(rad, height_slope);
    Real dr2 = rad*rad - 2. * rad * radpl * cos(phi-phipl) + radpl2 + smoothingB * smoothingB;
    Real s2 = dr2/(4 * local_scale_height * local_scale_height);

    // TODO: Move scale height at planet to actual local scale height
    Real k_0 = PotentialFunctions::k0e(s2);
    Real k_1 = PotentialFunctions::k1e(s2);

    // Real scale = - gmp_root_2_pi/( 2 * local_scale_height * local_scale_height * local_scale_height * rad);

    //dPhidrad = - scale * (
  //    (k_0 - k_1) * (height_slope * dr2 - rad * rad + rad * radpl * cos(phi-phipl))
  //    + k_0 * 6 * height_slope * local_scale_height * local_scale_height);
  //  dPhidphi = scale * (k_0 - k_1) * rad * radpl  * sin(phi - phipl);

  Real scale = - gmp_root_2_pi/local_scale_height;

  Real ds2dr = (
  (1 - height_slope) * rad - height_slope * radpl * radpl / rad
  - (1 - 2 * height_slope) * cos(phi - phipl))/(2 * local_scale_height * local_scale_height);

  dPhidrad = scale * (-1 * height_slope * k_0/rad + ds2dr * (k_0 - k_1)) - (
        height_slope * gmp_root_2_pi / (rad * 2 * local_scale_height) * k_0 );
  dPhidphi = scale * (k_0 - k_1) * rad * radpl  * sin(phi - phipl)/(2 * local_scale_height * local_scale_height);

}


//----------------------------------------------------------------------------------------
// \fn CoolingFunction
//  \brief Perform beta style cooling on the disc
//----------------------------------------------------------------------------------------
void CoolingFunction(Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is, int ie) {
  Real rad,phi,z;
  Real tc, epstht, omega;
  Real epsthc, epsthn, deltae, de;
  Real ekin, ethn;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {

      #pragma omp simd
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,k);
        omega = OmegaKep(rad);

        tc = beta/omega;
        // expected thermal energy per unit mass = Pressure/density/(gamma - 1)
        epstht = PoverR(rad,phi,z)*igm1;
        // kinetic energy per unit mass
        ekin = 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i)))/cons(IDN,k,j,i)/cons(IDN,k,j,i);
        epsthc = (cons(IEN,k,j,i)-ekin)/cons(IDN,k,j,i);
        deltae = (epsthc-epstht);


        if (strict_iso) { // reset T
          de = 0.0;
        } else {
          de = std::exp(-dt/tc)*deltae;
        }
        epsthn = epstht+de;
        ethn = epsthn*cons(IDN,k,j,i);
        cons(IEN,k,j,i) = ekin + ethn;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// \fn IsothermalReset
//  \brief Do quasi local isothermality by directly resetting the temperature
//----------------------------------------------------------------------------------------
void IsothermalReset(Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is, int ie) {
  Real rad,phi,z;
  Real ekin, ethn;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      #pragma omp simd
      for (int i=is; i<=ie; ++i) {
        GetCylCoord(pco,rad,phi,z,i,j,k);
        // Kinetic energy
        ekin = 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i)))/cons(IDN,k,j,i);

        // Expected thermal energy
        ethn = PoverR(rad,phi,z)* igm1 *cons(IDN,k,j,i);
        // reset total energy
        cons(IEN,k,j,i) = ekin + ethn;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// \fn DampingZones
//  \brief Damp inner and outer radial boundaries in density and vr
//----------------------------------------------------------------------------------------
void DampingZones(MeshBlock *pmb, Coordinates *pco, Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &prim, int ks, int ke, int js, int je, int is, int ie) {
  AthenaArray<Real> *rumd = pmb->pmy_mesh->ruser_mesh_data;
  int64_t &lx1=pmb->loc.lx1;
  int ioff = pmb->block_size.nx1*lx1;
  int nc = pmb->pmy_mesh->ncycle;
  Real tt = pmb->pmy_mesh->time;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      #pragma omp simd
      for (int i=is; i<=ie; ++i) {

        Real vrd, vphid, vzd, rhod, td, ramp;
        Real omega, tc, epstht, epsthc, expfac;
        Real deltavr, dvr, deltarho, drho, deltae, de;
        Real rad,phi,z;
        Real prev_rho, prev_mom1;
        int iglob = i + ioff;

        GetCylCoord(pco,rad,phi,z,i,j,k);

        if ((rad < rdampi) || (rad > rdampo)) {

          if (rad < rdampi) {  // Damping zones
            ramp = (rdampi-rad)/(rdampi-rmin);
          } else if (rad > rdampo) {
            ramp = (rad-rdampo)/(rmax-rdampo);
          }

          // Expected profile
          VelProfileCyl(pmb,rad,phi,z,vrd,vphid,vzd);
          rhod = DenProfileCyl(rad,phi,z);

          if (dampaverage){
            // Get average values to damp to
            vrd = rumd[4](iglob);
            rhod = rumd[6](iglob);
          }
          ramp = ramp*ramp;
          td = tdamp/OmegaKep(rad)/ramp;
          expfac = std::exp(-dt/td);

          prev_rho = cons(IDN,k,j,i);
          prev_mom1 = cons(IM1,k,j,i);
          deltavr = prev_mom1/prev_rho - vrd;
          dvr = expfac * deltavr;

          if (dampdens) {
            deltarho = cons(IDN,k,j,i) - rhod;
            drho = expfac * deltarho;
            cons(IDN,k,j,i) = rhod + drho;
          }

          // Damp v_r
          cons(IM1,k,j,i) = cons(IDN,k,j,i)*(vrd + dvr);

          if (NON_BAROTROPIC_EOS) {
            // Expected thermal energy per unit mass
            epstht = PoverR(rad,phi,z)*igm1;
            // Calculte thermal energy per unit mass
            epsthc = (cons(IEN,k,j,i)- 0.5 * (SQR(prev_mom1) + SQR(cons(IM2,k,j,i)))/prev_rho) / prev_rho;

            // New energy = kinetic + original thermal
            // or also damp thermal
            if (dampthermal) {
              deltae = epsthc - epstht;
              de = expfac*deltae;
              epsthc = epstht + de;
            }
            cons(IEN,k,j,i) = 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i)))/cons(IDN,k,j,i) + cons(IDN,k,j,i)*(epsthc);
          }
        }
      }
    }
  }

}

//----------------------------------------------------------------------------------------
// \fn OmegaKep
//  \brief Utility function for keplerian angular velocity
//----------------------------------------------------------------------------------------
Real OmegaKep(const Real r) {
  Real r3 = r*r*r;
  return std::sqrt(gm0/r3);
}



//----------------------------------------------------------------------------------------
// \fn MyViscosity
//  \brief Add radial power law viscosity
//----------------------------------------------------------------------------------------
void MyViscosity(HydroDiffusion *phdif, MeshBlock *pmb, const  AthenaArray<Real> &w, const AthenaArray<Real> &bc, int is, int ie, int js, int je, int ks, int ke) {
  Coordinates *coord = pmb->pcoord;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      #pragma omp simd
      for (int i=is; i<=ie; ++i) {
        Real rad = coord->x1v(i);
        // Pull the nu_iso from my actual work
        phdif->nu(HydroDiffusion::DiffProcess::iso,k,j,i) = nu_iso*ViscCoeffCell(rad);
      }
    }
  }
}

// simple power-law viscosity
Real ViscCoeffCell(Real rad) {
  return std::pow(rad, nuslope);
}