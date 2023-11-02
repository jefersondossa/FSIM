#include "IntRule0d.h"
#include <iostream> 
#include <vector>
#include <math.h>
#include <cmath>
using namespace std;


IntRule0d::IntRule0d(){

}

IntRule0d::IntRule0d(int order) {
    if (order != 0) {
        PanicButton();
    }

    SetOrder(order);

}


void IntRule0d::SetOrder(int order) {
    fOrder = order;  
    
    if (order != 0) {
        fOrder = 0;
        // PanicButton();//Apenas comentar esse debug stop e adicionar a linha de cima
    }

    int npoints = 1;
    
    fPoints.resize(npoints, 0);
    fWeights.resize(npoints);
    fWeights[0] = 1.;
    
}
