
#include <iostream>
#include <math.h>

#include "DenseEigen.h"

#include "Tensor3D.h"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
using namespace Catch::literals;

#define fTolerance  1.e-10

void CheckInvariants(Tensor3D &Tensor3D)
{
    // std::cout << Tensor3D.MatrixForm() << std::endl;
    
    //Check I2
    REAL I2 = Tensor3D.I2();
    REQUIRE(fabs(I2) > fTolerance);

    //Check I3
    REAL I3 = Tensor3D.I3();
    REQUIRE(fabs(I3) > fTolerance);

    //Check J2
    REAL J2 = Tensor3D.J2();
    REQUIRE(fabs(J2) > fTolerance);

    //Check J3
    REAL J3 = Tensor3D.J3();
    REQUIRE(fabs(J3) > fTolerance);

}


TEST_CASE("Tensor3D_test","[Tensor3D]")
{

    SECTION("Check Invariants"){
        int nTests = 10;
        for (int i = 0; i < nTests; i++){
            //Create a Tensor3D with random seed
            VecDouble testTensor3D = VecDouble::Random(6);
            Tensor3D t(testTensor3D);        
            CheckInvariants(t);
        }
    }
    

    

    
    
}
