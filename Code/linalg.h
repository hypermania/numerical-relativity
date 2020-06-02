#ifndef LINALGH
#define LINALGH

#include <cstdlib>
#include <iostream>
#include <vector>

#include "utility.h"
#include "numvec.h"
//#include "nummat.h"
#include "tridiag.h"

// Solves Ax=v, assuming A is a non-singular nxn matrix
//NumVec solve(const NumMat &A, const NumVec &v);
NumVec solve_tridiag(const TriDiag &A, const NumVec &v);
NumVec &solve_tridiag_inplace(TriDiag &U, NumVec &x);

//NumMat identity_matrix(int n);
TriDiag identity_tridiag(int n);

#endif
