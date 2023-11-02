#include <stdio.h>
#include "IntRuleHexahedron.h"


static int ComputingSymmetricCubatureRule(int order, MatrixDouble &Points, VecDouble &Weights);
void TransformBarycentricCoordInCartesianCoord(long double baryvec[], long double weightvec[], MatrixDouble &Points, VecDouble &Weight);

IntRuleHexahedron::IntRuleHexahedron() {
    PanicButton(); // this integration rule was not implemented yet
}

IntRuleHexahedron::IntRuleHexahedron(int order) : IntRule(order) {
    SetOrder(order);
}

void IntRuleHexahedron::SetOrder(int order) {

    if (order < 0 || order > MaxOrder()) {
        PanicButton();
    }

}
