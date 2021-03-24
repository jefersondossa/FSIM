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
#ifndef PART_QUADRATURE_H
#define PART_QUADRATURE_H

#include "QuadraticShapeFunction.hpp"

/// Defines the domain integration sub-element partitioned Hammer quadrature

template<int DIM>
class IntegQuadratureSpecial{

public:
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
        pointCoord[0][0] =   0.166666666666667   ;
        pointCoord[1][0] =   0.050643253661728   ;
        pointCoord[2][0] =   0.398713492676543   ;
        pointCoord[3][0] =   0.050643253661728   ;
        pointCoord[4][0] =   0.235071032052557   ;
        pointCoord[5][0] =   0.235071032052557   ;
        pointCoord[6][0] =   0.029857935894885   ;
        pointCoord[7][0] =   0.333333333333333   ;
        pointCoord[8][0] =   0.101286507323457   ;
        pointCoord[9][0] =   0.449356746338272   ;
        pointCoord[10][0] =   0.449356746338272   ;
        pointCoord[11][0] =   0.264928967947442   ;
        pointCoord[12][0] =   0.470142064105115   ;
        pointCoord[13][0] =   0.264928967947442   ;
        pointCoord[14][0] =   0.666666666666667   ;
        pointCoord[15][0] =   0.550643253661728   ;
        pointCoord[16][0] =   0.898713492676543   ;
        pointCoord[17][0] =   0.550643253661728   ;
        pointCoord[18][0] =   0.735071032052557   ;
        pointCoord[19][0] =   0.735071032052558   ;
        pointCoord[20][0] =   0.529857935894885   ;
        pointCoord[21][0] =   0.166666666666667   ;
        pointCoord[22][0] =   0.050643253661728   ;
        pointCoord[23][0] =   0.398713492676543   ;
        pointCoord[24][0] =   0.050643253661728   ;
        pointCoord[25][0] =   0.235071032052557   ;
        pointCoord[26][0] =   0.235071032052557   ;
        pointCoord[27][0] =   0.029857935894885   ;

        pointCoord[0][1] =   0.166666666666667   ;
        pointCoord[1][1] =   0.050643253661729   ;
        pointCoord[2][1] =   0.050643253661729   ;
        pointCoord[3][1] =   0.398713492676544   ;
        pointCoord[4][1] =   0.029857935894885   ;
        pointCoord[5][1] =   0.235071032052557   ;
        pointCoord[6][1] =   0.235071032052557   ;
        pointCoord[7][1] =   0.333333333333333   ;
        pointCoord[8][1] =   0.449356746338272   ;
        pointCoord[9][1] =   0.101286507323457   ;
        pointCoord[10][1] =   0.449356746338272   ;
        pointCoord[11][1] =   0.264928967947442   ;
        pointCoord[12][1] =   0.264928967947442   ;
        pointCoord[13][1] =   0.470142064105115   ;
        pointCoord[14][1] =   0.166666666666667   ;
        pointCoord[15][1] =   0.050643253661729   ;
        pointCoord[16][1] =   0.050643253661729   ;
        pointCoord[17][1] =   0.398713492676544   ;
        pointCoord[18][1] =   0.029857935894885   ;
        pointCoord[19][1] =   0.235071032052557   ;
        pointCoord[20][1] =   0.235071032052557   ;
        pointCoord[21][1] =   0.666666666666667   ;
        pointCoord[22][1] =   0.550643253661729   ;
        pointCoord[23][1] =   0.550643253661729   ;
        pointCoord[24][1] =   0.898713492676544   ;
        pointCoord[25][1] =   0.529857935894885   ;
        pointCoord[26][1] =   0.735071032052558   ;
        pointCoord[27][1] =   0.735071032052557   ;


        pointWeight[0] =    0.0281250000000000  ;
        pointWeight[1] =    0.0157423975681034  ;
        pointWeight[2] =    0.0157423975681034  ;
        pointWeight[3] =    0.0157423975681034  ;
        pointWeight[4] =    0.0165492690985632  ;
        pointWeight[5] =    0.0165492690985632  ;
        pointWeight[6] =    0.0165492690985632  ;
        pointWeight[7] =    0.0281250000000000  ;
        pointWeight[8] =    0.0157423975681034  ;
        pointWeight[9] =    0.0157423975681034  ;
        pointWeight[10] =   0.0157423975681034  ;
        pointWeight[11] =   0.0165492690985632  ;
        pointWeight[12] =   0.0165492690985632  ;
        pointWeight[13] =   0.0165492690985632  ;
        pointWeight[14] =   0.0281250000000000  ;
        pointWeight[15] =   0.0157423975681034  ;
        pointWeight[16] =   0.0157423975681034  ;
        pointWeight[17] =   0.0157423975681034  ;
        pointWeight[18] =   0.0165492690985632  ;
        pointWeight[19] =   0.0165492690985632  ;
        pointWeight[20] =   0.0165492690985632  ;
        pointWeight[21] =   0.0281250000000000  ;
        pointWeight[22] =   0.0157423975681034  ;
        pointWeight[23] =   0.0157423975681034  ;
        pointWeight[24] =   0.0157423975681034  ;
        pointWeight[25] =   0.0165492690985632  ;
        pointWeight[26] =   0.0165492690985632  ;
        pointWeight[27] =   0.0165492690985632  ;

    
    }


private:
    //List of integration points coordinates
    double pointCoord[28][2];

    //List of integration points weights
    double pointWeight[28];

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






