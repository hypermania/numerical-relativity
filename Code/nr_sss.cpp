#include "nr_sss.h"

/*
  Space discretization is given by r_i = (i-1/2)h, h = L/(N-1).
  Dimension of the linear system is N-1.
*/


TriDiag laplacian_tridiag(int N, double L){
  double h = L / (N-1);
  
  TriDiag laplacian = TriDiag(N-1, std::array<double, 3>({0, 0, 0}));

  laplacian[0][1] = -R(1.5, h) * R(1.5, h) / (R(1.0, h) * R(1.0, h));
  laplacian[0][2] = R(1.5, h) * R(1.5, h) / (R(1.0, h) * R(1.0, h));
  for(int i = 2; i < N-1; i++){
    double term1 = R(i-0.5, h) * R(i-0.5, h) / (R(i, h) * R(i, h));
    double term2 = R(i+0.5, h) * R(i+0.5, h) / (R(i, h) * R(i, h));
    laplacian[i-1][0] = term1;
    laplacian[i-1][1] = -term1 - term2;
    laplacian[i-1][2] = term2;
  }
  laplacian[N-1-1][0] = R(N-1.5, h) * R(N-1.5, h) / (R(N-1, h) * R(N-1, h));
  laplacian[N-1-1][1] = -R(N-1.5, h) * R(N-1.5, h) / (R(N-1, h) * R(N-1, h)) - R(N-0.5, h) * R(N-0.5, h) / (R(N, h) * R(N, h)) * (1 - R(N-1, h)/R(N, h));

  laplacian *= 1 / (h * h);
  
  return laplacian;
}

void add_laplacian_boundary(int N, double L, NumVec &f){
  double h = L / (N-1);
  f[N-1-1] += R(N-0.5, h) * R(N-0.5, h) / (R(N-1, h) * R(N-1, h) * R(N, h) * h);
}

NonLinearTerm compute_nonlinear_term(int N, double L, const ParticleProfile &particles, const NumVec &f){
  NonLinearTerm term = {
			NumVec(N-1, 0),
			TriDiag(N-1, std::array<double, 3>({0, 0, 0}))
  };
  
  double h = L / (N-1);
  for(int i = 0; i < particles.M; i++){
    double m = particles.m[i];
    double r = particles.r[i];
    
    int pos = ceil(r/h+0.5-1);
    double f_val = f[pos];
    
    double V = particles.u_r[i] * particles.u_r[i] + (particles.u_phi[i] * particles.u_phi[i]) / (r * r);
    double W_sqr = 1 + V / (f_val * f_val * f_val * f_val);

    double F_val = (-1/(2*r*r)) * m * (1/f_val) * sqrt(W_sqr) * (1/h);
    term.F[pos] += F_val;
    double DF_val = F_val * ((1/f_val) + 2 * (1/W_sqr) * V / (f_val * f_val * f_val * f_val * f_val));
    term.DF[pos][1] += DF_val;
  }
    
  return term;
}


Gauge compute_gauge(int N, double L, const ParticleProfile &particles, const NumVec &f){
  Gauge gauge = {NumVec(N-1, 0), NumVec(N-1, 0)};
  
  double h = L / (N-1);
  double M = (f[N-2] - 1) * 2 * R(N-1, h);

  // Compute radial momentum density and stress
  NumVec S_r(N-1, 0); // Radial momentum density
  NumVec S_rr(N-1, 0); // Radial stress
  for(int i = 0; i < particles.M; i++){
    double m = particles.m[i];
    double r = particles.r[i];
    double u_r = particles.u_r[i];
    
    int pos = ceil(r/h+0.5-1);
    double f_val = f[pos];
    
    double V = u_r * u_r + (particles.u_phi[i] * particles.u_phi[i]) / (r * r);
    double W_sqr = 1 + V / (f_val * f_val * f_val * f_val);

    double nW = 1 / (4 * M_PI * (f_val*f_val*f_val*f_val*f_val*f_val) * (r*r) * h);
    S_r[pos] += m * nW * u_r;
    S_rr[pos] += m * nW / sqrt(W_sqr) * u_r * u_r;
  }

  double alpha_prefactor = 1 - (M / (2 * L)) * (M / (2 * L));
  double alpha_accum = 0;
  double beta_accum = 0;
  // Integrate from r_max=L to r, obtain alpha and beta along the way

  double f_N = (R(N-1, h) * f[N-2] + h) / R(N, h);
  double drlnA = (f_N - f[N-2-1]) / (h * f[N-2]);
  double alpha_integrand = 0.5 * R(N-1, h) * (drlnA * drlnA + 8 * M_PI * S_rr[N-2]) / (1 + R(N-1, h) * drlnA);
  alpha_accum -= alpha_integrand * h;
  gauge.alpha[N-2] = alpha_prefactor * exp(alpha_accum);
  
  double beta_integrand = gauge.alpha[N-2] * 4 * M_PI * S_r[N-2] / (1 + R(N-1, h) * drlnA);
  beta_accum -= beta_integrand * h;
  gauge.beta[N-2] = R(N-1, h) * beta_accum;

  for(int i = N-2; i >= 2; i--){ // i is the convention on the writeup (i = 1,...,N-1)
    double drlnA = (f[i+1-1] - f[i-1-1]) / (h * f[i-1]);
    double alpha_integrand = 0.5 * R(i, h) * (drlnA * drlnA + 8 * M_PI * S_rr[i-1]) / (1 + R(i, h) * drlnA);
    alpha_accum -= alpha_integrand * h;
    gauge.alpha[i-1] = alpha_prefactor * exp(alpha_accum);

    double beta_integrand = gauge.alpha[i-1] * 4 * M_PI * S_r[i-1] / (1 + R(i, h) * drlnA);
    beta_accum -= beta_integrand * h;
    gauge.beta[i-1] = R(i, h) * beta_accum;
  }

  double f_0 = f[0];
  drlnA = (f[1] - f_0) / (h * f[0]);
  alpha_integrand = 0.5 * R(1, h) * (drlnA * drlnA + 8 * M_PI * S_rr[0]) / (1 + R(1, h) * drlnA);
  alpha_accum -= alpha_integrand * h;
  gauge.alpha[0] = alpha_prefactor * exp(alpha_accum);
  
  beta_integrand = gauge.alpha[0] * 4 * M_PI * S_r[0] / (1 + R(1, h) * drlnA);
  beta_accum -= beta_integrand * h;
  gauge.beta[0] = R(1, h) * beta_accum;
  
  return gauge;
}

void update_particles(int N, double L, ParticleProfile &particles, const NumVec &f, const Gauge &gauge, double dt){
  double h = L / (N-1);
  for(int i = 0; i < particles.M; i++){
    double r = particles.r[i];
    double u_r = particles.u_r[i];
    double u_phi = particles.u_phi[i];
    
    int pos = ceil(r/h+0.5-1);
    double f_val = f[pos];
    
    double A_sqr = f_val * f_val * f_val * f_val;
    double W = sqrt(1 + (u_r * u_r + (u_phi * u_phi) / (r * r)) / A_sqr);

    double alpha = gauge.alpha[pos];
    double beta = gauge.beta[pos];
    double u0 = W / alpha;


    double dralpha, drbeta, drlnA;
    if(pos == 0){
      dralpha = (gauge.alpha[1] - gauge.alpha[0]) / (2 * h);
      drbeta = (gauge.beta[1] - gauge.beta[0]) / (2 * h);
      drlnA = (f[1] - f[0]) / (h * f[0]);
    } else if(pos == N-2) {
      // This should never happen, all particles should be within [0, L)
      dralpha = 0;
      drbeta = 0;
      drlnA = ((R(N-1, h) * f[N-2] + h) / R(N, h) - f[N-2-1]) / (h * f[N-2]);
    } else {
      dralpha = (gauge.alpha[pos+1] - gauge.alpha[pos-1]) / (2 * h);
      drbeta = (gauge.beta[pos+1] - gauge.beta[pos-1]) / (2 * h);
      drlnA = (f[pos+1] - f[pos-1]) / (h * f[pos]);
    }
    
    double dr_dt = u_r / (u0 * A_sqr) - beta;
    double du_r_dt = -W * dralpha + u_r * drbeta + u_r * u_r * drlnA / (u0 * A_sqr) + (u_phi * u_phi / (u0 * R(pos+1, h) * R(pos+1, h) * A_sqr)) * (1 / R(pos+1, h) + drlnA);

    particles.r[i] += dr_dt * dt;
    particles.u_r[i] += du_r_dt * dt;
  }
}

void update_rays(int N, double L, RayProfile &rays, const NumVec &f, const Gauge &gauge, double dt){
  double h = L / (N-1);
  
  for(int i = 0; i < rays.M; i++){
    double r = rays.r[i];
    int pos = ceil(r/h+0.5-1);
    double A = f[pos] * f[pos];
    double alpha = gauge.alpha[pos];
    double beta = gauge.beta[pos];

    double dr_dt = alpha / A - beta;
    rays.r[i] += dr_dt * dt;
  }
}

void evolve_metric(int N, double L, NumVec &f, const Gauge &gauge, double dt){
  double h = L / (N-1);
  NumVec df(N-1, 0);
  df[0] = dt * gauge.beta[0] * ((f[1] - f[0]) / (2 * h) + f[0] / (2 * R(1, h)));
  for(int i = 1; i <= N-3; i++){
    df[i] = dt * gauge.beta[i] * ((f[i+1] - f[i-1]) / (2 * h) + f[i] / (2 * R(i+1, h)));
  }
  df[N-2] = dt * gauge.beta[N-2] * (((R(N-1, h) * f[N-2] + h) / R(N, h) - f[N-2-1]) / (2 * h) + f[N-2] / (2 * R(N-1, h)));

  f += df;
}
