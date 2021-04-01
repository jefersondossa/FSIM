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

template<int DIM, int DEG>
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

    int getNumberOfIntegrationPoints() const {
        return numIntegPoints;
    }

    /// Returns the integration point coordinate
    /// @param int integration point index @param int adimensional direction
    /// @return integration point adimensional coordinates
    double PointList(int i, int j){return pointCoord[i][j];};

    /// Retuns the integration point weight
    /// @param int integration point index @return integration point weight
    double WeightList(int i) {return pointWeight[i];};

    /// Sets the domain quadrature and store the values of the integration points and weights 
    void setQuadrature();

    /// Constructor of the domain integration quadrature
    IntegQuadrature(){
        setQuadrature();
    }

private:
    const static int numIntegPoints = -5*DIM-8*DEG+6*DIM*DEG+9;

    ///List of integration points coordinates
    double pointCoord[numIntegPoints][DIM];

    ///List of integration points weights
    double pointWeight[numIntegPoints];

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
void IntegQuadrature<2,1>::setQuadrature(){


    pointCoord[0][0] = 0.166666666666666666666666666666667e+0;
    pointCoord[0][1] = 0.666666666666666666666666666666667e+0;

    pointCoord[1][0] = 0.666666666666666666666666666666667e+0;
    pointCoord[1][1] = 0.166666666666666666666666666666667e+0;

    pointCoord[2][0] = 0.166666666666666666666666666666667e+0;
    pointCoord[2][1] = 0.166666666666666666666666666666667e+0;


    pointWeight[0]  =  0.166666666666666666666666666666667e+0;
    pointWeight[1]  =  0.166666666666666666666666666666667e+0;
    pointWeight[2]  =  0.166666666666666666666666666666667e+0;

    return;
}


template<>
void IntegQuadrature<2,2>::setQuadrature(){

    double g1 = (6. - std::sqrt(15.)) / 21.;
    double g2 = (6. + std::sqrt(15.)) / 21.;

    double w1 = (155. - std::sqrt(15.)) / 2400.;
    double w2 = (155. + std::sqrt(15.)) / 2400.;

    pointCoord[0][0] = 1. / 3.;
    pointCoord[0][1] = 1. / 3.;
        
    pointCoord[1][0] = 1. - 2.*g1;
    pointCoord[1][1] = g1;
      
    pointCoord[2][0] = g1;
    pointCoord[2][1] = 1. - 2.*g1;
      
    pointCoord[3][0] = g1;
    pointCoord[3][1] = g1;
      
    pointCoord[4][0] = g2;
    pointCoord[4][1] = g2;
      
    pointCoord[5][0] = 1. - 2.*g2;
    pointCoord[5][1] = g2;
      
    pointCoord[6][0] = g2;
    pointCoord[6][1] = 1. - 2.*g2;

    pointWeight[0] = 0.11250;
    pointWeight[1] = w1;
    pointWeight[2] = w1;
    pointWeight[3] = w1;
    pointWeight[4] = w2;
    pointWeight[5] = w2;
    pointWeight[6] = w2;

    // double g1 = (8. - std::sqrt(10.) + std::sqrt(38. - 44.*std::sqrt(2./5.))) / 18.;
    // double g2 = (8. - std::sqrt(10.) - std::sqrt(38. - 44.*std::sqrt(2./5.))) / 18.;
        
    // double w1 = (620. + std::sqrt(213125. - 53320.*std::sqrt(10))) / 3720.;
    // double w2 = (620. - std::sqrt(213125. - 53320.*std::sqrt(10))) / 3720.;

    // pointCoord[0][0] = g1;
    // pointCoord[0][1] = g1;
      
    // pointCoord[1][0] = g1;
    // pointCoord[1][1] = 1. - 2.*g1;
      
    // pointCoord[2][0] = 1. - 2.*g1;
    // pointCoord[2][1] = g1;
      
    // pointCoord[3][0] = g2;
    // pointCoord[3][1] = g2;
      
    // pointCoord[4][0] = g2;
    // pointCoord[4][1] = 1. - 2.*g2;
      
    // pointCoord[5][0] = 1. - 2.*g2;
    // pointCoord[5][1] = g2;

    // pointWeight[0] = w1;
    // pointWeight[1] = w1;
    // pointWeight[2] = w1;
    // pointWeight[3] = w2;
    // pointWeight[4] = w2;
    // pointWeight[5] = w2;
    
};


template<>
void IntegQuadrature<3,1>::setQuadrature(){


    pointCoord[0][0] = 0.138196601125010515179541316563436e+0;
    pointCoord[0][1] = 0.138196601125010515179541316563436e+0;
    pointCoord[0][2] = 0.585410196624968454461376050309691e+0;

    pointCoord[1][0] = 0.138196601125010515179541316563436e+0;
    pointCoord[1][1] = 0.585410196624968454461376050309691e+0;
    pointCoord[1][2] = 0.138196601125010515179541316563436e+0;

    pointCoord[2][0] = 0.585410196624968454461376050309691e+0;
    pointCoord[2][1] = 0.138196601125010515179541316563436e+0;
    pointCoord[2][2] = 0.138196601125010515179541316563436e+0;

    pointCoord[3][0] = 0.138196601125010515179541316563436e+0;
    pointCoord[3][1] = 0.138196601125010515179541316563436e+0;
    pointCoord[3][2] = 0.138196601125010515179541316563436e+0;


    pointWeight[0] =  0.416666666666666666666666666666667e-1;
    pointWeight[1] =  0.416666666666666666666666666666667e-1;
    pointWeight[2] =  0.416666666666666666666666666666667e-1;
    pointWeight[3] =  0.416666666666666666666666666666667e-1;


    return;
}

template<>
void IntegQuadrature<3,2>::setQuadrature(){

    pointCoord[0][0] = 0.310885919263300609797345733763458e+0;
    pointCoord[0][1] = 0.310885919263300609797345733763458e+0;
    pointCoord[0][2] = 0.673422422100981706079627987096265e-1;

    pointCoord[1][0] = 0.310885919263300609797345733763458e+0;
    pointCoord[1][1] = 0.673422422100981706079627987096265e-1;
    pointCoord[1][2] = 0.310885919263300609797345733763458e+0;

    pointCoord[2][0] = 0.673422422100981706079627987096265e-1;
    pointCoord[2][1] = 0.310885919263300609797345733763458e+0;
    pointCoord[2][2] = 0.310885919263300609797345733763458e+0;

    pointCoord[3][0] = 0.310885919263300609797345733763458e+0;
    pointCoord[3][1] = 0.310885919263300609797345733763458e+0;
    pointCoord[3][2] = 0.310885919263300609797345733763458e+0;

    pointCoord[4][0] = 0.927352503108912264023239137370306e-1;
    pointCoord[4][1] = 0.927352503108912264023239137370306e-1;
    pointCoord[4][2] = 0.721794249067326320793028258788908e+0;

    pointCoord[5][0] = 0.927352503108912264023239137370306e-1;
    pointCoord[5][1] = 0.721794249067326320793028258788908e+0;
    pointCoord[5][2] = 0.927352503108912264023239137370306e-1;

    pointCoord[6][0] = 0.721794249067326320793028258788908e+0;
    pointCoord[6][1] = 0.927352503108912264023239137370306e-1;
    pointCoord[6][2] = 0.927352503108912264023239137370306e-1;

    pointCoord[7][0] = 0.927352503108912264023239137370306e-1;
    pointCoord[7][1] = 0.927352503108912264023239137370306e-1;
    pointCoord[7][2] = 0.927352503108912264023239137370306e-1;

    pointCoord[8][0] = 0.455037041256496494918805262793394e-1;
    pointCoord[8][1] = 0.454496295874350350508119473720661e+0;
    pointCoord[8][2] = 0.454496295874350350508119473720661e+0;

    pointCoord[9][0] = 0.454496295874350350508119473720661e+0;
    pointCoord[9][1] = 0.455037041256496494918805262793394e-1;
    pointCoord[9][2] = 0.454496295874350350508119473720661e+0;

    pointCoord[10][0] = 0.455037041256496494918805262793394e-1;
    pointCoord[10][1] = 0.455037041256496494918805262793394e-1;
    pointCoord[10][2] = 0.454496295874350350508119473720661e+0;

    pointCoord[11][0] = 0.455037041256496494918805262793394e-1;
    pointCoord[11][1] = 0.454496295874350350508119473720661e+0;
    pointCoord[11][2] = 0.455037041256496494918805262793394e-1;

    pointCoord[12][0] = 0.454496295874350350508119473720661e+0;
    pointCoord[12][1] = 0.455037041256496494918805262793394e-1;
    pointCoord[12][2] = 0.455037041256496494918805262793394e-1;

    pointCoord[13][0] = 0.454496295874350350508119473720661e+0;
    pointCoord[13][1] = 0.454496295874350350508119473720661e+0;
    pointCoord[13][2] = 0.455037041256496494918805262793394e-1;

    pointWeight[0] = 0.187813209530026417998642753888811e-1;
    pointWeight[1] = 0.187813209530026417998642753888811e-1;
    pointWeight[2] = 0.187813209530026417998642753888811e-1;
    pointWeight[3] = 0.187813209530026417998642753888811e-1;
    pointWeight[4] = 0.122488405193936582572850342477213e-1;
    pointWeight[5] = 0.122488405193936582572850342477213e-1;
    pointWeight[6] = 0.122488405193936582572850342477213e-1;
    pointWeight[7] = 0.122488405193936582572850342477213e-1;
    pointWeight[8] = 0.709100346284691107301157135337624e-2;
    pointWeight[9] = 0.709100346284691107301157135337624e-2;
    pointWeight[10] = 0.709100346284691107301157135337624e-2;
    pointWeight[11] = 0.709100346284691107301157135337624e-2;
    pointWeight[12] = 0.709100346284691107301157135337624e-2;
    pointWeight[13] = 0.709100346284691107301157135337624e-2;
 
};


#endif
