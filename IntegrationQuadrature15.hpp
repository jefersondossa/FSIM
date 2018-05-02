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
#include "LinearShapeFunction.hpp"

template<int DIM>
class IntegQuadrature{
    //Defines the type "PointCoord" which stores the 
    //integration points coordinates
    typedef ublas::bounded_matrix<double, 8*DIM-9,DIM>  PointCoord;

    //Defines the type "PointWeight" which stores the
    //integration points weights
    typedef ublas::bounded_vector<double, 8*DIM-9>      PointWeight;

    //Integration point logical vector
    typedef ublas::bounded_vector<bool, 8*DIM-9>        PointLogical;

    //Defines the numerical integration iterator
    typedef typename PointWeight::iterator              QuadratureListIt;

    //Defines vector of nodal values
    typedef ublas::bounded_vector<double, 4*DIM-2>      NodalValuesQuad;
    typedef ublas::bounded_vector<double, DIM+1>        NodalValuesLin;

public:
    //Returns the index of the first integration point
    QuadratureListIt begin() {
        return pointWeight.begin();
    }

    //Returns the index of the last integration point
    QuadratureListIt end() {
        return pointWeight.end();
    }

    //Returns the integration point coordinate
    double PointList(int i, int j); 

    //Retuns the integration point weight
    double WeightList(int i);
  
    //Interpolate variables
    double interpolateQuadraticVariable(NodalValuesQuad nValues, int point);
    double interpolateLinearVariable(NodalValuesLin nValues, int point);

private:
    //List of integration points coordinates
    PointCoord pointCoord;

    //List of integration points weights
    PointWeight pointWeight;

    //Defines shape functions
    QuadShapeFunction<DIM> shapeQuad;
    LinShapeFunction<DIM>  shapeLin;

    //Values of velocity shape functins
    typename QuadShapeFunction<DIM>::Values      phi_;     
    //Values of pressure shape functins
    typename LinShapeFunction<DIM>::Values       phip_; 
};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::PointList(int i, int j){
    
    pointCoord(0,0) = 1./3.;
    pointCoord(0,1) = 1./3.;
        
    pointCoord(1,0) = 0.797426985353087;
    pointCoord(1,1) = 0.101286507323456;
      
    pointCoord(2,0) = 0.101286507323456;
    pointCoord(2,1) = 0.797426985353087;
      
    pointCoord(3,0) = 0.101286507323456;
    pointCoord(3,1) = 0.101286507323456;
      
    pointCoord(4,0) = 0.470142064105115;
    pointCoord(4,1) = 0.470142064105115;
      
    pointCoord(5,0) = 0.059715871789770;
    pointCoord(5,1) = 0.470142064105115;
      
    pointCoord(6,0) = 0.470142064105115;
    pointCoord(6,1) = 0.059715871789770;
  
    return pointCoord(i,j);
};

template<>
double IntegQuadrature<3>::PointList(int i, int j){
    
    pointCoord(0,0) = 1. / 4.;
    pointCoord(0,1) = 1. / 4.;
    pointCoord(0,2) = 1. / 4.;

    pointCoord(1,0) = 0.;
    pointCoord(1,1) = 1. / 3.;
    pointCoord(1,2) = 1. / 3.;

    pointCoord(2,0) = 1. / 3.;
    pointCoord(2,1) = 1. / 3.;
    pointCoord(2,2) = 1. / 3.;

    pointCoord(3,0) = 1. / 3.;
    pointCoord(3,1) = 1. / 3.;
    pointCoord(3,2) = 0.;

    pointCoord(4,0) = 1. / 3.;
    pointCoord(4,1) = 0.;
    pointCoord(4,2) = 1. / 3.;

    pointCoord(5,0) = 8. / 11.;
    pointCoord(5,1) = 1. / 11.;
    pointCoord(5,2) = 1. / 11.;

    pointCoord(6,0) = 1. / 11.;
    pointCoord(6,1) = 1. / 11.;
    pointCoord(6,2) = 1. / 11.;

    pointCoord(7,0) = 1. / 11.;
    pointCoord(7,1) = 1. / 11.;
    pointCoord(7,2) = 8. / 11.;

    pointCoord(8,0) = 1. / 11.;
    pointCoord(8,1) = 8. / 11.;
    pointCoord(8,2) = 1. / 11.;

    pointCoord(9,0) = 0.4334498464263357;
    pointCoord(9,1) = 0.0665501535736643;
    pointCoord(9,2) = 0.0665501535736643;

    pointCoord(10,0) = 0.0665501535736643;
    pointCoord(10,1) = 0.4334498464263357;
    pointCoord(10,2) = 0.0665501535736643;

    pointCoord(11,0) = 0.0665501535736643;
    pointCoord(11,1) = 0.0665501535736643;
    pointCoord(11,2) = 0.4334498464263357;

    pointCoord(12,0) = 0.0665501535736643;
    pointCoord(12,1) = 0.4334498464263357;
    pointCoord(12,2) = 0.4334498464263357;

    pointCoord(13,0) = 0.4334498464263357;
    pointCoord(13,1) = 0.0665501535736643;
    pointCoord(13,2) = 0.4334498464263357;

    pointCoord(14,0) = 0.4334498464263357;
    pointCoord(14,1) = 0.4334498464263357;
    pointCoord(14,2) = 0.0665501535736643;

    return pointCoord(i,j);
}

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::WeightList(int i){
    
    pointWeight(0) = 0.11250;
    pointWeight(1) = 0.125939180544827 / 2.;
    pointWeight(2) = 0.125939180544827 / 2.;
    pointWeight(3) = 0.125939180544827 / 2.;
    pointWeight(4) = 0.132394152788506 / 2.;
    pointWeight(5) = 0.132394152788506 / 2.;
    pointWeight(6) = 0.132394152788506 / 2.; 

    return pointWeight(i);
};

template<>
double IntegQuadrature<3>::WeightList(int i){
    
    pointWeight(0) = 0.1817020685825351 / 6.;
    pointWeight(1) = 0.0361607142857143 / 6.;
    pointWeight(2) = 0.0361607142857143 / 6.;
    pointWeight(3) = 0.0361607142857143 / 6.;
    pointWeight(4) = 0.0361607142857143 / 6.;
    pointWeight(5) = 0.0698714945161738 / 6.;
    pointWeight(6) = 0.0698714945161738 / 6.; 
    pointWeight(7) = 0.0698714945161738 / 6.;
    pointWeight(8) = 0.0698714945161738 / 6.;
    pointWeight(9) = 0.0656948493683187 / 6.;
    pointWeight(10) = 0.0656948493683187 / 6.; 
    pointWeight(11) = 0.0656948493683187 / 6.;
    pointWeight(12) = 0.0656948493683187 / 6.;
    pointWeight(13) = 0.0656948493683187 / 6.;
    pointWeight(14) = 0.0656948493683187 / 6.;

    return pointWeight(i);
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 2> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 6; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};

template<>
double IntegQuadrature<3>::interpolateQuadraticVariable(NodalValuesQuad nValues,
                                                        int point){
    
    ublas::bounded_vector<double, 3> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    xsi(2) = PointList(point,2);
    
    shapeQuad.evaluate(xsi,phi_);
    
    for (int i = 0; i < 10; i++){
        int_value += nValues(i) * phi_(i);
    };

    return int_value;
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::interpolateLinearVariable(NodalValuesLin nValues, 
                                                     int point){
    
    ublas::bounded_vector<double, 2> xsi;

    double int_value = 0.;

    xsi(0) = PointList(point,0);
    xsi(1) = PointList(point,1);
    
    shapeLin.evaluate(xsi,phip_);
    
    for (int i = 0; i < 3; i++){
        int_value += nValues(i) * phip_(i);
    };

    return int_value;
};

#endif
