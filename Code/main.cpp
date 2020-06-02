#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include "math.h"
#include "time.h"
#include "stdio.h"

#include "utility.h"
#include "numvec.h"
//#include "nummat.h"
#include "tridiag.h"
#include "linalg.h"
#include "nr_sss.h"
#include "particle.h"


int main(int argc, char **argv){

  // Set the discretization scheme: N-1 points over [0, L]
  // Dimension of the linear system is N-1
  int N = 10000;
  double L = 150;
  
  // Initial particle ensemble config
  int nParticle = 200000; // Total number of particles composing the dust star
  double rStar = 10; // Radius of the dust star
  double mParticle = 0.00001; // Rest mass of each particle
  ParticleProfile particles = uniform_static_ensemble(nParticle, rStar, mParticle);

  // Setup particle tracers
  std::vector<unsigned int> idxParticleTracers(5, 0);
  for(int i = 0; i < 5; i++){
    idxParticleTracers[i] = floor(0.2 * (i+1) * nParticle) - 1;
  }

  // Setup Lagrangian matter tracer
  std::vector<double> restMassFractions(5, 0);
  for(int i = 0; i < 5; i++){
    restMassFractions[i] = 0.2 * (i+1);
  }
  
  // Ray profile for tracking the event horizon
  // Note that the rays cannot go out of boundary [0, L], otherwise an error will occur
  //RayProfile rays = ray_ensemble(100, 0.5, 1.5);
  std::vector<int> rayStartingStep(50, 0);
  for(int i = 0; i < 50; i++){
    rayStartingStep[i] = 28750 + 50 * i;
  }
  double rayLowerRadius = 0.05;
  double rayUpperRadius = 1.0;
  int rayNum = 40;
  std::vector<RayProfile> rays;
  int nextRaySet = 0;
  
  // Precompute negative Laplacian matrix + lambda for future use
  TriDiag laplacian = laplacian_tridiag(N, L);  // Compute the Laplacian
  double lambda = 2;  // add lambda for stability
  TriDiag T = (-1) * laplacian + lambda * identity_tridiag(N-1);
  
  // Solve for f at t=0
  NumVec f(N-1, 1.0);  // Initial guess for f
  while(true){
    NonLinearTerm term = compute_nonlinear_term(N, L, particles, f);
    
    TriDiag A = T + term.DF;
    NumVec X = laplacian * f;
    add_laplacian_boundary(N, L, X);
    X -= term.F;
    
    solve_tridiag_inplace(A, X);
    f += X;
    std::cout << "||df|| = " << max_norm(X) << std::endl;    
    if(max_norm(X) < 1e-6){ // Stop when the solution is barely changing
      break;
    }
  }

  
  // Measure how much time it takes to run the program
  clock_t t;
  t = clock();

  // Iteration parameters
  double dt = 0.004; // The time step
  int nSteps = 40000; // How many steps to run
  int savePeriod = 4000; // Save file every "savePeriod" steps
  
  for(int n = 0; n < nSteps; n++){
    // Compute gauge functions
    Gauge gauge = compute_gauge(N, L, particles, f);

    // Print some quantities on screen
    // The quantities are: number of steps, t, central alpha, central f
    printf("(n=%.5d) t=%.4f, alpha_c=%.4e, f_c=%.4f\n", n, n*dt, gauge.alpha[0], f[0]);

    
    // Update the matter distribution
    update_particles(N, L, particles, f, gauge, dt);

    // Update the spectator light rays
    if(nextRaySet < rayStartingStep.size() && rayStartingStep[nextRaySet] == n){
      rays.push_back(ray_ensemble(rayNum, rayLowerRadius, rayUpperRadius));
      nextRaySet++;
    }
    for(long unsigned int i = 0; i < rays.size(); i++){
      update_rays(N, L, rays[i], f, gauge, dt);
    }
    
    // Solve for the metric given the new matter distribution
    evolve_metric(N, L, f, gauge, dt); // Optional update on f to speed up convergence
    while(true){
      NonLinearTerm term = compute_nonlinear_term(N, L, particles, f);

      TriDiag A = T + term.DF;
      NumVec X = laplacian * f;
      add_laplacian_boundary(N, L, X);
      X -= term.F;
      
      solve_tridiag_inplace(A, X);
      f += X;
      if(max_norm(X) < 1e-6){
	break;
      }
    }

    
    // Write f and alpha to file periodically    
    char filename[128];
    std::ofstream file;

    if(n % savePeriod == 0){
      sprintf(filename, "solution/%04.4f.tsv", n * dt);
      file.open(filename);
      for(int i = 0; i < N-1; i++){
	file << std::setprecision(10) << R(i+1, L/(N-1)) << "\t" << f[i] << "\n";
      }
      file.close();
      
      sprintf(filename, "alpha/%04.4f.tsv", n * dt);
      file.open(filename);
      for(int i = 0; i < N-1; i++){
	file << std::setprecision(10) << R(i+1, L/(N-1)) << "\t" << gauge.alpha[i] << "\n";
      }
      file.close();
    }

    
    // Write particles to file
    sprintf(filename, "particles/particles.tsv");
    file.open(filename, std::ofstream::app);
    file << n * dt;
    for(long unsigned int i = 0; i < idxParticleTracers.size(); i++){
      double r = particles.r[idxParticleTracers[i]];
      int pos = ceil(r/(L/(N-1))+0.5-1);
      int f_val = f[pos];
      file << std::setprecision(10) << "\t" << (f_val * f_val * r);
    }
    file << "\n";
    file.close();

    // Write rays to file
    for(long unsigned int j = 0; j < rays.size(); j++){
      sprintf(filename, "rays/rays%d.tsv", (int)j);
      file.open(filename, std::ofstream::app);
      file << n * dt;
      for(int i = 0; i < rays[j].M; i++){
	double r = rays[j].r[i];
	int pos = ceil(r/(L/(N-1))+0.5-1);
	int f_val = f[pos];
	file << std::setprecision(10) << "\t" << (f_val * f_val * r);
      }
      file << "\n";
      file.close();
    }

    // Write Lagrangian matter tracer to file
    sprintf(filename, "mat_tracers/mat_tracers.tsv");
    file.open(filename, std::ofstream::app);
    file << n * dt;
    NumVec r_sorted(particles.r);
    std::sort(r_sorted.begin(), r_sorted.end(), std::less<double>());
    for(long unsigned int i = 0; i < restMassFractions.size(); i++){
      int idx = floor(restMassFractions[i] * particles.M - 1);
      double r = r_sorted[idx];
      int pos = ceil(r/(L/(N-1))+0.5-1);
      int f_val = f[pos];
      file << std::setprecision(10) << "\t" << (f_val * f_val * r);
    }
    file << "\n";
    file.close();
    
  }

  // Print the time used for evolving the system
  t = clock() - t;
  std::cout << "Time used = " << (float)t/(CLOCKS_PER_SEC) << std::endl;

  
}
