#include <iostream>
#include "read.h"

extern "C"{
 void froutine_(void) ;
 void finit_(void) ;
 void readpt_(void) ;
}


int main()
{
	std::cout << "hello world\n";
  //readpt_() ;
  //BaryonRich br("Au30mix_i1_Bps.dat");
  Fireball("Au30mix_i1_Fps.dat");
	return 0 ;
}
