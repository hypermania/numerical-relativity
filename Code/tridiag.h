#ifndef TRIDIAGH
#define TRIDIAGH

#include <cstdlib>
#include <iostream>
#include <vector>
#include <array>

#include "utility.h"
#include "numvec.h"
//#include "nummat.h"


//a[i - 1], a[i], a[i + 1]
typedef std::vector<std::array<double, 3>> TriDiag;

TriDiag operator+(const TriDiag &A, const TriDiag &B);
TriDiag operator-(const TriDiag &A, const TriDiag &B);
TriDiag operator*(double a, const TriDiag &B);
NumVec operator*(const TriDiag &A, const NumVec &v);


// In-place variants of the above operators
TriDiag &operator+=(TriDiag &A, const TriDiag &B);
TriDiag &operator-=(TriDiag &A, const TriDiag &B);
TriDiag &operator*=(TriDiag &A, double a);
//NumVec &operator*=(NumVec &v, const TriDiag &A);


// Output stream
std::ostream &operator<<(std::ostream &os, const TriDiag &A);

// Utility function
//TriDiag nummat_to_tridiag(const NumMat &A);

#endif
