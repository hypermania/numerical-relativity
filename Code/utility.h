// utility.h

#ifndef UTILTYH
#define UTILTYH
#include <cstdlib>
#include <iostream>
#include <string>

void error( std::string s );
void warning( std::string s, std::ostream& out=std::cerr );
inline double sq( double t )
{return t*t;}

#endif
