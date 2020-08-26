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

#ifndef INTEG_QUADRATURE11_H
#define INTEG_QUADRATURE11_H

#include "QuadraticShapeFunction.hpp"

/// Defines a special quadrature rule (equal to the normal quadradure rule) used for the integration of problems with coincident meshes

template<int DIM>
class IntegQuadratureSpecial{
public:
    int getNumberOfIntegrationPoints(){
        return (*(&pointWeight+1) - pointWeight);
    }

    //Returns the index of the first integration point
    double* begin() {
        return std::begin(pointWeight);
    }

    //Returns the index of the last integration point
    double* end() {
        return std::end(pointWeight);
    }

    //Returns the integration point coordinate
    double PointList(int i, int j); 

    //Retuns the integration point weight
    double WeightList(int i);
  
    //Interpolate variables
    double interpolateQuadraticVariable(double nValues[], int point);

    IntegQuadratureSpecial(){
        pointCoord[0][0] = 1. / 3.;
        pointCoord[0][1] = 1. / 3.;
            
        pointCoord[1][0] = (9. + 2. * std::sqrt(15.)) / 21.;
        pointCoord[1][1] = (6. - std::sqrt(15.)) / 21.;
          
        pointCoord[2][0] = (6. - std::sqrt(15.)) / 21.;
        pointCoord[2][1] = (9. + 2. * std::sqrt(15.)) / 21.;
          
        pointCoord[3][0] = (6. - std::sqrt(15.)) / 21.;
        pointCoord[3][1] = (6. - std::sqrt(15.)) / 21.;
          
        pointCoord[4][0] = (6. + std::sqrt(15.)) / 21.;
        pointCoord[4][1] = (6. + std::sqrt(15.)) / 21.;
          
        pointCoord[5][0] = (9. - 2. * std::sqrt(15.)) / 21.;
        pointCoord[5][1] = (6. + std::sqrt(15.)) / 21.;
          
        pointCoord[6][0] = (6. + std::sqrt(15.)) / 21.;
        pointCoord[6][1] = (9. - 2. * std::sqrt(15.)) / 21.;

        pointWeight[0] = 0.11250;
        pointWeight[1] = (155. - std::sqrt(15.)) / 2400.;
        pointWeight[2] = (155. - std::sqrt(15.)) / 2400.;
        pointWeight[3] = (155. - std::sqrt(15.)) / 2400.;
        pointWeight[4] = (155. + std::sqrt(15.)) / 2400.;
        pointWeight[5] = (155. + std::sqrt(15.)) / 2400.;
        pointWeight[6] = (155. + std::sqrt(15.)) / 2400.;
    }

private:
    //List of integration points coordinates
    double pointCoord[7][2];

    //List of integration points weights
    double pointWeight[7];

    //Defines shape functions
    QuadShapeFunction<DIM> shapeQuad;

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadratureSpecial<2>::PointList(int i, int j){
    
    return pointCoord[i][j];
};


//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadratureSpecial<2>::WeightList(int i){

    return pointWeight[i];
};


//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double IntegQuadratureSpecial<2>::interpolateQuadraticVariable(double nValues[], int point){
    
    double xsi[2];

    double int_value = 0.;

    xsi[0] = PointList(point,0);
    xsi[1] = PointList(point,1);

    double phi_[6] = {};
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 6; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};


#endif
