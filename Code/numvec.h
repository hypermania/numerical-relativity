#ifndef NUMVECH
#define NUMVECH
// interface for numerical vectors

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "utility.h"

typedef std::vector<double> NumVec;


NumVec operator+(const NumVec& A, const NumVec& B); // A+B
NumVec operator-(const NumVec& A, const NumVec& B); // A-B
NumVec operator*(double a, const NumVec& B); // a*B
double operator,(const NumVec& A, const NumVec& B); // (A,B)

// In-place variants of the above operators
NumVec &operator+=(NumVec &A, const NumVec &B);
NumVec &operator-=(NumVec &A, const NumVec &B);
NumVec &operator*=(NumVec &A, double a);

double max_norm(const NumVec &A);


std::ostream& operator<<(std::ostream& os, const NumVec& A);
#endif






