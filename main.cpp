#include <iostream>


extern "C"{
 void froutine_(void) ;
 void finit_(void) ;
}


int main()
{
	std::cout << "hello world\n" ;
  finit_() ;
	froutine_() ;
	return 0 ;
}
