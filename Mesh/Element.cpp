#include "Element.h"
#include <cassert>


void Element::setConnectivity(VecInt connect){
    for (int iconnect = 0; iconnect < connect.size(); iconnect++){
        fConnect[iconnect] = fMesh->ConnectVec()[connect[iconnect]];
    } 
}
