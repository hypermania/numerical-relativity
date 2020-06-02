#include "linalg.h"

/*
NumVec solve(const NumMat &A, const NumVec &v){

  if(A.size() != A[0].size() || A.size() != v.size()){
    error("solve: incompatible matrix/vector dimensions");
  }
  
  int n;
  n = A.size();
  
  NumMat U(A);
  for(int i = 0; i < n; i++){
    U[i].push_back(v[i]);
  }

  for(int c = 0; c < n; c++){
    // Possible row swap
    if(U[c][c] == 0){
      // Find row to swap
      int e;
      for(e = c + 1; e < n; e++){
	if(U[e][c] != 0){
	  break;
	}
      }
      if(e == n){
	error("solve: singular matrix");
      }
      std::swap(U[c], U[e]);
    }
    // Row addition
    for(int r = c + 1; r < n; r++){
      U[r] = U[r] - (U[r][c] / U[c][c]) * U[c];
    }
  }

  //  std::cout << U << std::endl;
  
  NumVec x(n, 0);
  for(int i = n - 1; i >= 0; i--){
    x[i] = U[i][n];
    for(int j = i + 1; j < n; j++){
      x[i] -= x[j] * U[i][j];
    }
    x[i] /= U[i][i];
  }
  
  return x;
}
*/

NumVec solve_tridiag(const TriDiag &A, const NumVec &v){
  
  if(A.size() != v.size() || A.size() <= 2){
    error("solve: incompatible matrix/vector dimensions");
  }

  //clock_t t;
  //t = clock();
  
  int n = A.size();
  
  TriDiag U = TriDiag(A);
  NumVec x = NumVec(v);

  /*
  t = clock() - t;
  std::cout << "copy = " << (float)t/(CLOCKS_PER_SEC) << std::endl;
  t = clock();
  */
  
  for(int c = 0; c < n - 1; c++){
    /*
    // Possible row swap
    if(U[c][1] == 0){
      // Find row to swap
      for(int e = c + 1; e < n; e++){
	if(U[e][c] != 0){
	  break;
	}
      }
      if(e == n){
	error("solve: singular matrix");
      }
      std::swap(U[c], U[e]);
    }
    */
    // Row addition
    int r = c + 1;
    double scale = U[r][0] / U[c][1];
    U[r][0] = U[r][0] - scale * U[c][1];
    U[r][1] = U[r][1] - scale * U[c][2];
    x[r] = x[r] - scale * x[c];
  }

  /*
  t = clock() - t;
  std::cout << "Gaussian elim = " << (float)t/(CLOCKS_PER_SEC) << std::endl;
  t = clock();
  */
  
  x[n - 1] /= U[n - 1][1];
  for(int i = n - 2; i >= 0; i--){
    x[i] -= x[i + 1] * U[i][2];
    x[i] /= U[i][1];
  }

  /*
  t = clock() - t;
  std::cout << "backward sub = " << (float)t/(CLOCKS_PER_SEC) << std::endl;
  t = clock();
  */
  
  return x;
}

NumVec &solve_tridiag_inplace(TriDiag &U, NumVec &x){
  if(U.size() != x.size() || U.size() <= 2){
    error("solve: incompatible matrix/vector dimensions");
  }
  
  int n = U.size();
  
  for(int c = 0; c < n - 1; c++){
    // Row addition
    int r = c + 1;
    double scale = U[r][0] / U[c][1];
    U[r][0] = U[r][0] - scale * U[c][1];
    U[r][1] = U[r][1] - scale * U[c][2];
    x[r] = x[r] - scale * x[c];
  }
  
  x[n - 1] /= U[n - 1][1];
  for(int i = n - 2; i >= 0; i--){
    x[i] -= x[i + 1] * U[i][2];
    x[i] /= U[i][1];
  }
  return x;
}

/*
NumMat identity_matrix(int n){
  NumMat I = NumMat(n, NumVec(n, 0));
  for(int i = 0; i < n; i++){
    I[i][i] = 1;
  }
  return I;
}
*/

TriDiag identity_tridiag(int n){
  TriDiag I = TriDiag(n, std::array<double, 3>());
  for(int i = 0; i < n; i++){
    I[i][0] = 0;
    I[i][1] = 1;
    I[i][2] = 0;
  }
  return I;
}

