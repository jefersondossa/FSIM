//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef INTEG_QUADRATURE_H
#define INTEG_QUADRATURE_H

#include "QuadraticShapeFunction.hpp"

/// Defines the domain integration Hammer quadrature

template<int DIM>
class IntegQuadrature{
public:
    /// Returns the index of the first integration point
    /// @return first integration point index
    double* begin() {
        return std::begin(pointWeight);
    }

    /// Returns the index of the last integration point
    /// @return last integration point index
    double* end() {
        return std::end(pointWeight);
    }

    int getNumberOfIntegrationPoints(){
        return (*(&pointWeight+1) - pointWeight);
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j); 

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i);

    IntegQuadrature(){
        pointCoord[0][0] = 1. / 3.;
        pointCoord[0][1] = 1. / 3.;
            
        pointCoord[1][0] = (9. + 2. * sqrt(15.)) / 21.;
        pointCoord[1][1] = (6. - sqrt(15.)) / 21.;
          
        pointCoord[2][0] = (6. - sqrt(15.)) / 21.;
        pointCoord[2][1] = (9. + 2. * sqrt(15.)) / 21.;
          
        pointCoord[3][0] = (6. - sqrt(15.)) / 21.;
        pointCoord[3][1] = (6. - sqrt(15.)) / 21.;
          
        pointCoord[4][0] = (6. + sqrt(15.)) / 21.;
        pointCoord[4][1] = (6. + sqrt(15.)) / 21.;
          
        pointCoord[5][0] = (9. - 2. * sqrt(15.)) / 21.;
        pointCoord[5][1] = (6. + sqrt(15.)) / 21.;
          
        pointCoord[6][0] = (6. + sqrt(15.)) / 21.;
        pointCoord[6][1] = (9. - 2. * sqrt(15.)) / 21.;

        pointWeight[0] = 0.11250;
        pointWeight[1] = (155. - sqrt(15.)) / 2400.;
        pointWeight[2] = (155. - sqrt(15.)) / 2400.;
        pointWeight[3] = (155. - sqrt(15.)) / 2400.;
        pointWeight[4] = (155. + sqrt(15.)) / 2400.;
        pointWeight[5] = (155. + sqrt(15.)) / 2400.;
        pointWeight[6] = (155. + sqrt(15.)) / 2400.;
    }


private:
    //List of integration points coordinates
    double pointCoord[7][2];

    //List of integration points weights
    double pointWeight[7];

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::PointList(int i, int j){

    return pointCoord[i][j];
};

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::WeightList(int i){
    
    return pointWeight[i];
};

#endif
