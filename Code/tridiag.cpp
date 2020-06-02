#include "tridiag.h"

TriDiag operator+(const TriDiag& A, const TriDiag& B){
  
  unsigned int n = A.size();
  TriDiag C = TriDiag(A);
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      C[i][j] += B[i][j];
    }
  }
  return C;
}

TriDiag operator-(const TriDiag& A, const TriDiag& B){
  
  unsigned int n = A.size();
  TriDiag C = TriDiag(A);
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      C[i][j] -= B[i][j];
    }
  }
  return C;
}

TriDiag operator*(double a, const TriDiag& B){

  unsigned int n = B.size();
  TriDiag C = TriDiag(B);
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      C[i][j] *= a;
    }
  }
  return C;
}

NumVec operator*(const TriDiag& A, const NumVec& v){

  unsigned int n = A.size();
  NumVec w = NumVec(n);
  // w[i] = A[i][i-1]v[i-1] + A[i][i]v[i] + A[i][i+1]v[i+1]
  w[0] = A[0][1] * v[0] + A[0][2] * v[1];
  for(unsigned int i = 1; i < n - 1; i++){
    w[i] = A[i][0] * v[i - 1] + A[i][1] * v[i] + A[i][2] * v[i + 1];
  }
  w[n - 1] = A[n - 1][0] * v[n - 2] + A[n - 1][1] * v[n - 1];
  return w;
}

TriDiag &operator+=(TriDiag &A, const TriDiag &B){
  unsigned int n = A.size();
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      A[i][j] += B[i][j];
    }
  }
  return A;
}

TriDiag &operator-=(TriDiag &A, const TriDiag &B){
  unsigned int n = A.size();
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      A[i][j] -= B[i][j];
    }
  }
  return A;
}

TriDiag &operator*=(TriDiag &A, double a){
  unsigned int n = A.size();
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < 3; j++){
      A[i][j] *= a;
    }
  }
  return A;
}


std::ostream& operator<<(std::ostream& os, const TriDiag& A){

  unsigned int n = A.size();
  for(unsigned int i = 0; i < n; i++){
    for(unsigned int j = 0; j < n; j++){
      if(j == i - 1){
	os << A[i][0] << "\t";
      } else if(j == i){
	os << A[i][1] << "\t";
      } else if(j == i + 1){
	os << A[i][2] << "\t";
      } else {
	os << "0" << "\t";
      }
    }
    os << "\n";
  }

  return os;
}

/*
TriDiag nummat_to_tridiag(const NumMat &A){
  int n = A.size();
  TriDiag T = TriDiag(n);

  T[0][1] = A[0][0];
  T[0][2] = A[0][1];
  for(int i = 1; i < n - 1; i++){
    T[i][0] = A[i][i - 1];
    T[i][1] = A[i][i];
    T[i][2] = A[i][i + 1];
  }
  T[n - 1][0] = A[n - 1][n - 2];
  T[n - 1][1] = A[n - 1][n - 1];
  return T;
}
*/
