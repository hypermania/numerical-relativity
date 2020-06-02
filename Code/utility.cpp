// utility.cp

#include "utility.h"

void error( std::string s )
{// write message and die
  std::cerr << "*** Error ***\n";
  std::cerr << s << std::endl;
  exit(1);
}

void warning( std::string s, std::ostream& out )
{// write message but keep cranking
  out << "**warning: " << s << std::endl;
}
