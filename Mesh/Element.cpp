#include "Element.h"

void Element::setIntegPointWeightFunction() {
    PanicButton(); //It produces wrong results for constant weight functions
    // VecDouble xsi(DIM);
    // ShapeFunction shapeQuad(DIM,DEG);
    // VecDouble phi_(fMesh->NElNodes());
    
    // IntegQuadrature nQuad(DIM,DEG);
    // for(int i = 0; i < nQuad.getNumberOfIntegrationPoints(); i++) {
    //     intPointWeightFunctionPrev[i] = intPointWeightFunction[i];
    //     intPointWeightFunction[i] = 0.;
    // }

    // int index=0;

    // for(int it = 0; it < nQuad.getNumberOfIntegrationPoints(); it++){
        
    //    xsi[0] = nQuad.PointList(index,0);
    //    xsi[1] = nQuad.PointList(index,1);
            
    //    //Computes the velocity shape functions
    //    shapeQuad.Shape(xsi,phi_);

    //    for (int j=0; j<fMesh->NElNodes(); j++){
    //        intPointWeightFunction[index] += phi_[j] * fMesh->NodeVec()[fConnect[j]] -> getWeightFunction();
    //    };
    //    // intPointWeightFunction(index) = 1.;
    //    index++;
    // }; 

    // index = 0;
     return;
};