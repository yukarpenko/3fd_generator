#include <iostream>


extern "C"{
 void froutine_(void) ;
 void finit_(void) ;
 void readpt_(void) ;
}


int main()
{
	std::cout << "hello world\n" ;
  readpt_() ;
	return 0 ;
}
