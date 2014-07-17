#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE testModule

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Unit tests
#include "integratorGaussLegendreTest.h"



//==============================================================================
bool
init_function()
{
// do your own initialization here
// if it successful return true
// But, you CAN'T use testing tools here
return true;
}

//==============================================================================
int
main( int argc, char* argv[] )
{
  int res=::boost::unit_test::unit_test_main( &init_function, argc, argv );
  
  return res;
}



