#include "numvec.h"


NumVec operator+(const NumVec& A, const NumVec& B)
{ // A+B
  if( A.size() != B.size() )
    error("NumVec +: incompatible sizes");
  NumVec C(A);
  unsigned int i;
  for(i=0; i<A.size(); i++)
    C[i] += B[i];
  return C;
}

NumVec operator*( double a, const NumVec& B)
{ // a*B
  NumVec C(B);
  unsigned int i;
  for(i=0; i<B.size(); i++ )
    C[i] *= a;
  return C;
}

NumVec operator-(const NumVec& A, const NumVec& B)
{ // A-B
  if( A.size() != B.size() )
    error("NumVec -: incompatible sizes");
  NumVec C(A);
  unsigned int i;
  for(i=0; i<A.size(); i++)
    C[i] -= B[i];
  return C;

}

double operator,(const NumVec& A, const NumVec& B)
{ // (A,B)
  if( A.size() != B.size() )
    error("NumVec ,: incompatible sizes");
  double sum = 0.0;
  unsigned int i;
  for( i=0; i< A.size(); i++ )
    sum += A[i]*B[i];

  return sum;
}

NumVec &operator+=(NumVec &A, const NumVec &B){
  if(A.size() != B.size()){
    error("NumVec ,: incompatible sizes");
  }
  unsigned int size = A.size();
  for(unsigned int i = 0; i < size; i++){
    A[i] += B[i];
  }
  return A;
}

NumVec &operator-=(NumVec &A, const NumVec &B){
  if(A.size() != B.size()){
    error("NumVec ,: incompatible sizes");
  }
  unsigned int size = A.size();
  for(unsigned int i = 0; i < size; i++){
    A[i] -= B[i];
  }
  return A;
}

NumVec &operator*=(NumVec &A, double a){
  unsigned int size = A.size();
  for(unsigned int i = 0; i < size; i++){
    A[i] *= a;
  }
  return A;
}

double max_norm(const NumVec &A){
  double max = 0;
  unsigned int size = A.size();
  for(unsigned int i = 0; i < size; i++){
    if(std::abs(A[i]) > max){
      max = std::abs(A[i]);
    }
  }
  return max;
}

std::ostream& operator<<(std::ostream& os, const NumVec& A){
 // output for small numeric vectors
  unsigned int i;
   for(i=0; i<A.size(); i++){
     os << A[i] << "\t";
   }
 os << "\n";
 return os;
}
