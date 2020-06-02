#ifndef PARTICLEH
#define PARTICLEH

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "utility.h"
#include "numvec.h"
//#include "nummat.h"
#include "tridiag.h"
#include "time.h"

/*
  Provides the particle ensemble class.
  A ParticleProfile struct gives all the positions and velocities
  of an ensemble of particles at some time.
*/


typedef struct {
  int M; // number of particles
  NumVec r; // radial coordinate
  NumVec u_r; // radial velocity
  NumVec u_phi; // azimuthal velocity
  NumVec m; // mass
} ParticleProfile;

typedef struct {
  int M;
  NumVec r;
} RayProfile;
  

/*
  Creates an ensemble of M static particles of mass m, 
  initially at rest, in a sphere of radius R, with uniform density.
*/
ParticleProfile uniform_static_ensemble(double M, double R, double m);

/*
  Creates an ensemble of M outward going light rays over [a, b].
  Requires M >= 2.
*/
RayProfile ray_ensemble(int M, double a, double b);

#endif
