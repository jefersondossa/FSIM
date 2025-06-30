
#include <iostream>
#include <math.h>

#include "DenseEigen.h"

#include "Tensor.h"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
using namespace Catch::literals;

#define fTolerance  1.e-10

void CheckInvariants(Tensor &tensor)
{
    // std::cout << tensor.MatrixForm() << '\n';
    
    //Check I2
    double I2 = tensor.I2();
    REQUIRE(fabs(I2) > fTolerance);

    //Check I3
    double I3 = tensor.I3();
    REQUIRE(fabs(I3) > fTolerance);

    //Check J2
    double J2 = tensor.J2();
    REQUIRE(fabs(J2) > fTolerance);

    //Check J3
    double J3 = tensor.J3();
    REQUIRE(fabs(J3) > fTolerance);

}


TEST_CASE("tensor_test","[tensor]")
{

    SECTION("Check Invariants"){
        int nTests = 10;
        for (int i = 0; i < nTests; i++){
            //Create a tensor with random seed
            VecDouble testTensor = VecDouble::Random(6);
            Tensor t(testTensor);        
            CheckInvariants(t);
        }
    }
    

    

    
    
}
