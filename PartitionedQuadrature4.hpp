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
#include "LinearShapeFunction.hpp"


template<int DIM>
class PartQuadrature{

public:
public:
    //Defines the type "PointCoord" which stores the 
    //integration points coordinates
    typedef ublas::bounded_matrix<double, 28,DIM>  PointCoord;

    //Defines the type "PointWeight" which stores the
    //integration points weights
    typedef ublas::bounded_vector<double, 28>      PointWeight;

    //Integration point logical vector
    typedef ublas::bounded_vector<bool, 28>        PointLogical;

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
double PartQuadrature<2>::PointList(int i, int j){

pointCoord(0	,0) = 	0.333333333333333 	;
pointCoord(1	,0) = 	0.101286507323457 	;
pointCoord(2	,0) = 	0.449356746338272 	;
pointCoord(3	,0) = 	0.449356746338272 	;
pointCoord(4	,0) = 	0.264928967947442 	;
pointCoord(5	,0) = 	0.470142064105115 	;
pointCoord(6	,0) = 	0.264928967947442 	;
pointCoord(7	,0) = 	0.166666666666667 	;
pointCoord(8	,0) = 	0.050643253661728 	;
pointCoord(9	,0) = 	0.398713492676543 	;
pointCoord(10	,0) = 	0.050643253661728 	;
pointCoord(11	,0) = 	0.235071032052557 	;
pointCoord(12	,0) = 	0.235071032052557 	;
pointCoord(13	,0) = 	0.029857935894885 	;
pointCoord(14	,0) = 	0.666666666666667 	;
pointCoord(15	,0) = 	0.550643253661728 	;
pointCoord(16	,0) = 	0.898713492676543 	;
pointCoord(17	,0) = 	0.550643253661728 	;
pointCoord(18	,0) = 	0.735071032052557 	;
pointCoord(19	,0) = 	0.735071032052558 	;
pointCoord(20	,0) = 	0.529857935894885 	;
pointCoord(21	,0) = 	0.166666666666667 	;
pointCoord(22	,0) = 	0.050643253661728 	;
pointCoord(23	,0) = 	0.398713492676543 	;
pointCoord(24	,0) = 	0.050643253661728 	;
pointCoord(25	,0) = 	0.235071032052557 	;
pointCoord(26	,0) = 	0.235071032052557 	;
pointCoord(27	,0) = 	0.029857935894885 	;

pointCoord(0	,1) = 	0.166666666666667 	;
pointCoord(1	,1) = 	0.050643253661729 	;
pointCoord(2	,1) = 	0.050643253661729 	;
pointCoord(3	,1) = 	0.398713492676544 	;
pointCoord(4	,1) = 	0.029857935894885 	;
pointCoord(5	,1) = 	0.235071032052557 	;
pointCoord(6	,1) = 	0.235071032052557 	;
pointCoord(7	,1) = 	0.333333333333333 	;
pointCoord(8	,1) = 	0.101286507323457 	;
pointCoord(9	,1) = 	0.449356746338272 	;
pointCoord(10	,1) = 	0.449356746338272 	;
pointCoord(11	,1) = 	0.264928967947442 	;
pointCoord(12	,1) = 	0.470142064105115 	;
pointCoord(13	,1) = 	0.264928967947442 	;
pointCoord(14	,1) = 	0.166666666666667 	;
pointCoord(15	,1) = 	0.050643253661729 	;
pointCoord(16	,1) = 	0.050643253661729 	;
pointCoord(17	,1) = 	0.398713492676544 	;
pointCoord(18	,1) = 	0.029857935894885 	;
pointCoord(19	,1) = 	0.235071032052557 	;
pointCoord(20	,1) = 	0.235071032052557 	;
pointCoord(21	,1) = 	0.666666666666667 	;
pointCoord(22	,1) = 	0.550643253661729 	;
pointCoord(23	,1) = 	0.550643253661729 	;
pointCoord(24	,1) = 	0.898713492676544 	;
pointCoord(25	,1) = 	0.529857935894885 	;
pointCoord(26	,1) = 	0.735071032052558 	;
pointCoord(27	,1) = 	0.735071032052557 	;

    return pointCoord(i,j);
};

template<>
double PartQuadrature<3>::PointList(int i, int j){
    
    const double a = (1. + sqrt(5. / 14.)) / 4.;
    const double b = (1. - sqrt(5. / 14.)) / 4.;

    pointCoord(0,0) = 1. / 4.;
    pointCoord(0,1) = 1. / 4.;
    pointCoord(0,2) = 1. / 4.;

    pointCoord(1,0) = 11. / 14.;
    pointCoord(1,1) = 1. / 14.;
    pointCoord(1,2) = 1. / 14.;

    pointCoord(2,0) = 1. / 14.;
    pointCoord(2,1) = 11. / 14.;
    pointCoord(2,2) = 1. / 14.;

    pointCoord(3,0) = 1. / 14.;
    pointCoord(3,1) = 1. / 14.;
    pointCoord(3,2) = 11. / 14.;

    pointCoord(4,0) = 1. / 14.;
    pointCoord(4,1) = 1. / 14.;
    pointCoord(4,2) = 1. / 14.;

    pointCoord(5,0) = a;
    pointCoord(5,1) = a;
    pointCoord(5,2) = b;

    pointCoord(6,0) = a;
    pointCoord(6,1) = b;
    pointCoord(6,2) = a;

    pointCoord(7,0) = a;
    pointCoord(7,1) = b;
    pointCoord(7,2) = b;

    pointCoord(8,0) = b;
    pointCoord(8,1) = a;
    pointCoord(8,2) = a;

    pointCoord(9,0) = b;
    pointCoord(9,1) = a;
    pointCoord(9,2) = b;

    pointCoord(10,0) = b;
    pointCoord(10,1) = b;
    pointCoord(10,2) = a;

    return pointCoord(i,j);
}

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::WeightList(int i){
    
pointWeight(0) = 	0.0281250000000000 	;
pointWeight(1) = 	0.0157423975681034 	;
pointWeight(2) = 	0.0157423975681034 	;
pointWeight(3) = 	0.0157423975681034 	;
pointWeight(4) = 	0.0165492690985632 	;
pointWeight(5) = 	0.0165492690985632 	;
pointWeight(6) = 	0.0165492690985632 	;
pointWeight(7) = 	0.0281250000000000 	;
pointWeight(8) = 	0.0157423975681034 	;
pointWeight(9) = 	0.0157423975681034 	;
pointWeight(10) = 	0.0157423975681034 	;
pointWeight(11) = 	0.0165492690985632 	;
pointWeight(12) = 	0.0165492690985632 	;
pointWeight(13) = 	0.0165492690985632 	;
pointWeight(14) = 	0.0281250000000000 	;
pointWeight(15) = 	0.0157423975681034 	;
pointWeight(16) = 	0.0157423975681034 	;
pointWeight(17) = 	0.0157423975681034 	;
pointWeight(18) = 	0.0165492690985632 	;
pointWeight(19) = 	0.0165492690985632 	;
pointWeight(20) = 	0.0165492690985632 	;
pointWeight(21) = 	0.0281250000000000 	;
pointWeight(22) = 	0.0157423975681034 	;
pointWeight(23) = 	0.0157423975681034 	;
pointWeight(24) = 	0.0157423975681034 	;
pointWeight(25) = 	0.0165492690985632 	;
pointWeight(26) = 	0.0165492690985632 	;
pointWeight(27) = 	0.0165492690985632 	;



    return pointWeight(i);
};

template<>
double PartQuadrature<3>::WeightList(int i){
    
    pointWeight(0) = -74. / 5625.;
    pointWeight(1) = 343. / 45000.;
    pointWeight(2) = 343. / 45000.;
    pointWeight(3) = 343. / 45000.;
    pointWeight(4) = 343. / 45000.;
    pointWeight(5) = 56. / 2250.;
    pointWeight(6) = 56. / 2250.; 
    pointWeight(7) = 56. / 2250.;
    pointWeight(8) = 56. / 2250.; 
    pointWeight(9) = 56. / 2250.;
    pointWeight(10) = 56. / 2250.; 

    return pointWeight(i);
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double PartQuadrature<2>::interpolateQuadraticVariable(
                                                      NodalValuesQuad nValues,
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
double PartQuadrature<3>::interpolateQuadraticVariable(
                                                       NodalValuesQuad nValues,
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
double PartQuadrature<2>::interpolateLinearVariable(
                                                     NodalValuesLin nValues, 
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
