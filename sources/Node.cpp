#include "Node.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------CLEAR VARIABLES--------------------------------
//------------------------------------------------------------------------------
template<int DIM, int DEG>
void Node<DIM,DEG>::clearVariables(){
    pressure_ = 0.;   
    elemCorresp = 0;    
    weightFunction_ = 0.;

    for (int i = 0; i < DIM; ++i){
        if (constrainType[i] != 1){
            velocity_[i] = 0.;   
            previousVelocity_[i] = 0.;   
            acceleration_[i] = 0.;
            previousAcceleration_[i] = 0.;
        }
        lagMultiplier_[i] = 0.;
        xsiCorresp[i] = 0.;
    }
    
    return;
}; 

template class Node<2,1>;
template class Node<2,2>;
template class Node<2,3>;
template class Node<3,1>;
template class Node<3,2>;
template class Node<3,3>;