#include "IntRule0d.h"
#include <iostream> 
#include <vector>
#include <math.h>
#include <cmath>
using namespace std;


IntRule0d::IntRule0d(){

}

IntRule0d::IntRule0d(int order) {

#ifdef DEBUG_BUILD 
    if (order != 0) {
        PanicButton();
    }
#endif

    SetOrder(order);

}


void IntRule0d::SetOrder(int order) {
    fOrder = order;  

#ifdef DEBUG_BUILD    
    if (order != 0) {
        fOrder = 0;
    }
#endif
    int npoints = 1;
    
    fPoints.resize(npoints, 0);
    fWeights.resize(npoints);
    fWeights[0] = 1.;
    
}
