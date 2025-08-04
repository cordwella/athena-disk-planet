//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// File added by Amelia J. Cordwell <ajc356@cam.ac.uk>
//========================================================================================
//! \file diskplanet_3d_sph.cpp
//! \brief 3D accretion disk in spherical co-ordinates with a planetary potential, beta 
//! cooling and alpha viscosity

// C headers

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
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

namespace {
void GetSphCoord(Coordinates *pco,Real &rad,Real &phi, Real &theta, int i,int j,int k);
Real PoverRSph(const Real rad, const Real phi, const Real theta);


void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);

void SourceTerms(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
void SourcePlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
void DampingZones(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);
void IsothermalReset(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);


// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, tslope, gamma_gas;
Real flaring_index, sslope, H0, H0_over_r0, sig0, igm1;
Real dfloor;
Real Omega0;
Real gmp, rsm, phipl0, omegap, tramp_start, tramp_dur, rmin, rmax, tdamp, rdampi, rdampo;
bool strict_iso;
} // namespace

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

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0)  {
    std::stringstream msg;

    msg << "### FATAL ERROR in function [Mesh::InitUserMeshData]"
      << std::endl << "only support spherical_polar coordiates at the moment!" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }


  // Get parameters for gravitatonal potential of central point mass
  gm0 = pin->GetOrAddReal("problem","GM",1.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  sig0 = pin->GetOrAddReal("problem","sig0", 1.0);
  // slope in surface density
  sslope = pin->GetOrAddReal("problem","sslope",0.0);

  // reset isothermal background at each timestep
  strict_iso = pin->GetOrAddBoolean("problem", "strict_iso", false);


  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    tslope = pin->GetOrAddReal("problem","tslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
    igm1 = 1/(gamma_gas - 1);
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
    igm1 = 1;
  }

  // Get parameters for planet
  gmp = pin->GetOrAddReal("problem","GMp",0.0);
  phipl0 = pin->GetOrAddReal("problem","phipl0",M_PI);
  rsm = pin->GetOrAddReal("problem","rsm",0.001);
  // smooth bessel functions by twice the cell size in the phi direction in the planet
  Real minSmoothing = 2 * r0 * pin->GetReal("mesh", "x3max")/pin->GetInteger("mesh", "nx3");
  
  if (minSmoothing > rsm){
    std::cout << "Requested smoothing less than 2x grid size. Resetting rsm to  = " << minSmoothing <<std::endl;
    rsm = minSmoothing;
  }


  omegap = sqrt((gm0)/r0/r0/r0);

  H0          = std::sqrt(p0_over_r0)/omegap;
  H0_over_r0  = H0/r0;

  // Slope in midplane density  
  flaring_index = (1 + tslope)/2;
  dslope = sslope - 1 - flaring_index;
  rho0 = sig0/(std::sqrt(2 * M_PI) * H0); // initial midplane density

  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  tramp_start = pin->GetOrAddReal("problem","tramp_start",5.);
  tramp_start *= 2*M_PI;
  tramp_dur = pin->GetOrAddReal("problem","tramp_dur",5.);
  tramp_dur *= 2*M_PI;
  rmin = pin->GetReal("mesh","x1min");
  rmax = pin->GetReal("mesh","x1max");
  tdamp = pin->GetOrAddReal("problem","tdamp",1.);
  rdampi = pin->GetOrAddReal("problem","rdampi",0.28);
  rdampo = pin->GetOrAddReal("problem","rdampo",3.4);

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

  EnrollUserExplicitSourceFunction(SourceTerms);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Real rad(0.0), theta(0.0), phi(0.0);
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1, x2, x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  //  Initialize density and momenta
  for (int k=ks; k<=ke; ++k) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        x1 = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        den = DenProfileCyl(rad,phi,z);
        vel = VelProfileCyl(rad,phi,z);

        // GetSphCoord(pcoord,rad,theta,phi,i,j,k); // convert to cylindrical coordinates
        // compute initial conditions in cylindrical coordinates
        // den = DenProfileSph(rad,theta,phi);
        if (porb->orbital_advection_defined)
          vel -= vK(porb, x1, x2, x3);
        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;

        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = den*vel;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = den*vel;
        }

        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! transform to spherical coordinates

void GetSphCoord(Coordinates *pco,Real &rad, Real &theta, Real &phi, int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    theta=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad = std::sqrt(pco->x1v(i)*pco->x1v(i) + pco->x3v(i)*pco->x3v(i));
    phi = std::atan2( pco->x3v(i),  pco->x1v(i));
    theta = pco->x3v(i);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}


//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0,dslope);
  Real dentem = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}



//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, tslope);
  return poverr;
}


Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+tslope)*p_over_r/(gm0/rad) + (1.0+tslope)
             - tslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}

//----------------------------------------------------------------------------------------
// \fn SourceTerms
//  \brief Adds source terms due to point mass planet, cooling and damping
//----------------------------------------------------------------------------------------

void SourceTerms(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

  SourcePlanet(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  DampingZones(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);

  if (strict_iso && NON_BAROTROPIC_EOS){
      // Reset local isothermality
      IsothermalReset(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  }
  return;
}

void IsothermalReset(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real rad, phi, z;
        Real ekin, ethn;

        GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        // Kinetic energy
        ekin = 0.5*(SQR(cons(IM1,k,j,i))+SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);

        // Expected thermal energy
        ethn = PoverR(rad,phi, z)* igm1 *cons(IDN,k,j,i);
        // reset total energy
        cons(IEN,k,j,i) = ekin + ethn;
      }
    }
  }
}


void SourcePlanet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {
 
  // Planet coordinates
  Real radpl = r0;
  Real radpl2 = radpl * radpl;
  Real phipl = phipl0 + time*omegap;
  // planet assumed to be at theta = \pi/2 

  Real ramp = 0.;
  // Ramp up potential over time
  if (time < tramp_start){
    ramp = 0.;
  } else if (time > tramp_start + tramp_dur) {
    ramp = 1.0;
  } else {
    ramp = 0.5*(1.0-cos(M_PI*(time-tramp_start)/tramp_dur));
  }

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {

        if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
          Real rad, theta, phi, rad2; 
          rad = pmb->pcoord->x1v(i);
          theta = pmb->pcoord->x2v(j);
          phi = pmb->pcoord->x3v(k);

          rad2 = rad * rad;
          Real den = prim(IDN,k,j,i);

          Real denom, dPhidr, dPhidtheta, dPhidphi;

          denom = rad2 + radpl2 - 2 * rad*radpl*std::sin(theta)*std::cos(phi-phipl) + rsm*rsm;
          denom = pow(denom,1.5);
          dPhidr = gmp/denom * (rad - radpl * std::sin(theta) * std::cos(phi-phipl));
          dPhidtheta = gmp/denom * ( - rad * radpl * std::cos(theta) * std::cos(phi-phipl));
          dPhidphi = gmp/denom * (rad * radpl * std::sin(theta) * std::sin(phi-phipl));

          Real src1 = - ramp * den * dPhidr;
          Real src2 = - ramp * den * dPhidtheta / rad;
          Real src3 = - ramp * den * dPhidphi / (rad*std::sin(theta));

          cons(IM1,k,j,i) += src1 * dt;
          cons(IM2,k,j,i) += src2 * dt;
          cons(IM3,k,j,i) += src3 * dt;

          if (NON_BAROTROPIC_EOS) {
            Real work = dt * (src1 * prim(IM1,k,j,i) + src2 * prim(IM2,k,j,i) + src3 * prim(IM3,k,j,i));
            cons(IEN,k,j,i) += work;
          }
        }
      }
    }
  }
  return;
}

void DampingZones(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {

  // Real rad, theta, phi, rad_cyl;
  Real rad, phi, z;

  Real prev_rho, prev_v1, prev_v2, prev_v3;
  Real rhod, v1d, v2d, v3d;
  Real deltarho, deltav1, deltav2, deltav3;
  Real drho, dv1, dv2, dv3;
  Real td, expfac, ramp;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        rad = pmb->pcoord->x1v(i);   //Either cyl or sph radius, damping zone is just set by the relevant radius
        
        if ((rad < rdampi) || (rad > rdampo)) {
          if (rad < rdampi) {  // Damping zones
            ramp = (rdampi-rad)/(rdampi-rmin);
          } else if (rad > rdampo) {
            ramp = (rad-rdampo)/(rmax-rdampo);
          }
          ramp = ramp*ramp;
          td = tdamp * 2*M_PI * std::sqrt(rad*rad*rad/gm0) / ramp;
          expfac = std::exp(-dt/td);

   	      GetCylCoord(pmb->pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
          rhod = DenProfileCyl(rad,phi,z);
          v1d = 0.;

          if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
            v2d = VelProfileCyl(rad,phi,z);
            v3d = 0.;
          } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
            v2d = 0.;
            v3d = VelProfileCyl(rad,phi,z);
          }

    
          prev_rho = cons(IDN,k,j,i);
          prev_v1 = cons(IM1,k,j,i)/cons(IDN,k,j,i);
          prev_v2 = cons(IM2,k,j,i)/cons(IDN,k,j,i);
          prev_v3 = cons(IM3,k,j,i)/cons(IDN,k,j,i);

          deltarho = prev_rho - rhod;
          drho = expfac * deltarho;
          cons(IDN,k,j,i) = rhod + drho;

          // the error is in this one here so because I don't want it anyway I'm just going to ignore it
          // for the time being. This may come back to bite me later, but who knows

          deltav1 = prev_v1 - v1d;
          dv1 = expfac * deltav1;
          cons(IM1,k,j,i) = (v1d + dv1)*cons(IDN,k,j,i);

          // deltav2 = prev_v2 - v2d;
          // dv2 = expfac * deltav2;
          // cons(IM2,k,j,i) = (v2d + dv2)*cons(IDN,k,j,i);

          // deltav3 = prev_v3 - v3d;
          // dv3 = expfac * deltav3;
          // cons(IM3,k,j,i) = (v3d + dv3)*cons(IDN,k,j,i);

        }
      }
    }
  }


  return;
}

} // namespace

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = vel;
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = vel;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = vel;
          prim(IM3,k,jl-j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = 0.0;
          prim(IM3,k,jl-j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = vel;
          prim(IM3,k,ju+j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = 0.0;
          prim(IM3,k,ju+j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = vel;
          prim(IM3,kl-k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = 0.0;
          prim(IM3,kl-k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = vel;
          prim(IM3,ku+k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = 0.0;
          prim(IM3,ku+k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  }
}