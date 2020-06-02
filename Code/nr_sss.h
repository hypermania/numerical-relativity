#ifndef NR_SSS_H
#define NR_SSS_H

#include <cstdlib>
#include <iostream>
#include <vector>

#include "utility.h"
#include "numvec.h"
//#include "nummat.h"
#include "tridiag.h"
#include "particle.h"
#include "time.h"

/*
  Provides the ingredients for solving Einstein's equation
  in the spherically symmetric case, using mean field, 
  particle simulation scheme with isotropic coordinates
  and polar slicing.
*/

#define R(i, h) (((i)-1.0/2.0)*(h))

typedef struct {
  NumVec F;
  TriDiag DF;
} NonLinearTerm;

typedef struct {
  NumVec alpha;
  NumVec beta;
} Gauge;

/*
  Gives the spherical Laplacian in the r-coordinate [0, L],
  using Neumann boundary condition at r=0 and Robin boundary
  condition at r=L.
*/
TriDiag laplacian_tridiag(int N, double L);

/*
  There is a boundary term due to the Robin boundary condition,
  this function adds the boundary term to f.
*/
void add_laplacian_boundary(int N, double L, NumVec &f);

/*
  The RHS of the elliptic equation for f is nonlinear in f
  due to the non-linear dependence of density in f.
  Linearizing this term gives F(f+df) = F(f) + DF(f)df.
  This function computes F(f) and DF(f).
*/
NonLinearTerm compute_nonlinear_term(int N, double L, const ParticleProfile &particles, const NumVec &f);

/*
  Computes the gauge functions alpha(r) and beta(r) 
  given particle ensemble and (solved) metric component f.
*/
Gauge compute_gauge(int N, double L, const ParticleProfile &particles, const NumVec &f);

/*
  Given metric component function f, update the position
  and velocities of the particles in "particles" by time dt.
*/
void update_particles(int N, double L, ParticleProfile &particles, const NumVec &f, const Gauge &gauge, double dt);

/*
  Given metric component function f, update the position
  and velocities of the particles in "particles" by time dt.
*/
void update_rays(int N, double L, RayProfile &rays, const NumVec &f, const Gauge &gauge, double dt);

/*
  Evolve the metric function by dt given the gauge.
*/
void evolve_metric(int N, double L, NumVec &f, const Gauge &gauge, double dt);


#endif
