#include "particle.h"

ParticleProfile uniform_static_ensemble(double M, double R, double m){
  ParticleProfile particles;
  particles.M = M;
  particles.r = NumVec(M, 0);
  particles.u_r = NumVec(M, 0);
  particles.u_phi = NumVec(M, 0);
  particles.m = NumVec(M, m);

  double rho = M / ((4.0/3.0) * M_PI * R * R * R);
  for(int i = 0; i < M; i++){
    double r = pow((i+1) / (rho * (4.0/3.0) * M_PI), 1.0/3.0);
    particles.r[i] = r;
  }

  return particles;
}

RayProfile ray_ensemble(int M, double a, double b){
  RayProfile rays;
  rays.M = M;
  rays.r = NumVec(M, 0);
  for(int i = 0; i < M; i++){
    rays.r[i] = a + i * (b-a) / (M-1.0);
  }

  return rays;
}
