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

#ifndef INTEG_QUADRATURE37_H
#define INTEG_QUADRATURE37_H

#include "QuadraticShapeFunction.hpp"
#include "LinearShapeFunction.hpp"

template<int DIM>
class IntegQuadratureSpecial{
public:  
    //Defines the type "PointCoord" which stores the 
    //integration points coordinates
    typedef ublas::bounded_matrix<double, (89 - 26 * DIM),DIM>  PointCoord;

    //Defines the type "PointWeight" which stores the
    //integration points weights
    typedef ublas::bounded_vector<double, (89 - 26 * DIM)>      PointWeight;

    //Integration point logical vector
    typedef ublas::bounded_vector<bool, (89 - 26 * DIM)>        PointLogical;

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
double IntegQuadratureSpecial<2>::PointList(int i, int j){
    
    // pointCoord(0,0) = 1./3.;
    // pointCoord(0,1) = 1./3.;
        
    // pointCoord(1,0) = 0.797426985353087;
    // pointCoord(1,1) = 0.101286507323456;
      
    // pointCoord(2,0) = 0.101286507323456;
    // pointCoord(2,1) = 0.797426985353087;
      
    // pointCoord(3,0) = 0.101286507323456;
    // pointCoord(3,1) = 0.101286507323456;
      
    // pointCoord(4,0) = 0.470142064105115;
    // pointCoord(4,1) = 0.470142064105115;
      
    // pointCoord(5,0) = 0.059715871789770;
    // pointCoord(5,1) = 0.470142064105115;
      
    // pointCoord(6,0) = 0.470142064105115;
    // pointCoord(6,1) = 0.059715871789770;

    
    pointCoord(0,0) = 1. / 3.;
    pointCoord(0,1) = 1. / 3.;
        
    pointCoord(1,0) = 0.950275662924105565450352089520;
    pointCoord(1,1) = 0.024862168537947217274823955239;
      
    pointCoord(2,0) = 0.024862168537947217274823955239;
    pointCoord(2,1) = 0.950275662924105565450352089520;
      
    pointCoord(3,0) = 0.024862168537947217274823955239;
    pointCoord(3,1) = 0.024862168537947217274823955239;
      
    pointCoord(4,0) = 0.171614914923835347556304795551;
    pointCoord(4,1) = 0.414192542538082326221847602214;
      
    pointCoord(5,0) = 0.414192542538082326221847602214;
    pointCoord(5,1) = 0.171614914923835347556304795551;
      
    pointCoord(6,0) = 0.414192542538082326221847602214;
    pointCoord(6,1) = 0.414192542538082326221847602214;

    pointCoord(7,0) = 0.539412243677190440263092985511;
    pointCoord(7,1) = 0.230293878161404779868453507244;

    pointCoord(8,0) = 0.230293878161404779868453507244;
    pointCoord(8,1) = 0.539412243677190440263092985511;

    pointCoord(9,0) = 0.230293878161404779868453507244;
    pointCoord(9,1) = 0.230293878161404779868453507244;

    pointCoord(10,0) = 0.772160036676532561750285570113;
    pointCoord(10,1) = 0.113919981661733719124857214943;

    pointCoord(11,0) = 0.113919981661733719124857214943;
    pointCoord(11,1) = 0.772160036676532561750285570113;

    pointCoord(12,0) = 0.113919981661733719124857214943;
    pointCoord(12,1) = 0.113919981661733719124857214943;

    pointCoord(13,0) = 0.009085399949835353883572964740;
    pointCoord(13,1) = 0.495457300025082323058213517632;

    pointCoord(14,0) = 0.495457300025082323058213517632;
    pointCoord(14,1) = 0.009085399949835353883572964740;

    pointCoord(15,0) = 0.495457300025082323058213517632;
    pointCoord(15,1) = 0.495457300025082323058213517632;

    pointCoord(16,0) = 0.062277290305886993497083640527;
    pointCoord(16,1) = 0.468861354847056503251458179727;

    pointCoord(17,0) = 0.468861354847056503251458179727;
    pointCoord(17,1) = 0.062277290305886993497083640527;

    pointCoord(18,0) = 0.468861354847056503251458179727;
    pointCoord(18,1) = 0.468861354847056503251458179727;

    pointCoord(19,0) = 0.022076289653624405142446876931;
    pointCoord(19,1) = 0.851306504174348550389457672223;

    pointCoord(20,0) = 0.022076289653624405142446876931;
    pointCoord(20,1) = 0.126617206172027096933163647918;

    pointCoord(21,0) = 0.851306504174348550389457672223;
    pointCoord(21,1) = 0.022076289653624405142446876931;

    pointCoord(22,0) = 0.851306504174348550389457672223;
    pointCoord(22,1) = 0.126617206172027096933163647918;

    pointCoord(23,0) = 0.126617206172027096933163647918;
    pointCoord(23,1) = 0.022076289653624405142446876931;

    pointCoord(24,0) = 0.126617206172027096933163647918;
    pointCoord(24,1) = 0.851306504174348550389457672223;

    pointCoord(25,0) = 0.018620522802520968955913511549;
    pointCoord(25,1) = 0.689441970728591295496647976487;

    pointCoord(26,0) = 0.018620522802520968955913511549;
    pointCoord(26,1) = 0.291937506468887771754472382212;

    pointCoord(27,0) = 0.689441970728591295496647976487;
    pointCoord(27,1) = 0.018620522802520968955913511549;

    pointCoord(28,0) = 0.689441970728591295496647976487;
    pointCoord(28,1) = 0.291937506468887771754472382212;

    pointCoord(29,0) = 0.291937506468887771754472382212;
    pointCoord(29,1) = 0.018620522802520968955913511549;

    pointCoord(30,0) = 0.291937506468887771754472382212;
    pointCoord(30,1) = 0.689441970728591295496647976487;

    pointCoord(31,0) = 0.096506481292159228736516560903;
    pointCoord(31,1) = 0.635867859433872768286976979827;

    pointCoord(32,0) = 0.096506481292159228736516560903;
    pointCoord(32,1) = 0.267625659273967961282458816185;

    pointCoord(33,0) = 0.635867859433872768286976979827;
    pointCoord(33,1) = 0.096506481292159228736516560903;

    pointCoord(34,0) = 0.635867859433872768286976979827;
    pointCoord(34,1) = 0.267625659273967961282458816185;

    pointCoord(35,0) = 0.267625659273967961282458816185;
    pointCoord(35,1) = 0.096506481292159228736516560903;

    pointCoord(36,0) = 0.267625659273967961282458816185;
    pointCoord(36,1) = 0.635867859433872768286976979827;

    return pointCoord(i,j);
};

template<>
double IntegQuadratureSpecial<3>::PointList(int i, int j){
    
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
double IntegQuadratureSpecial<2>::WeightList(int i){
    
    // pointWeight(0) = 0.11250;
    // pointWeight(1) = 0.125939180544827 / 2.;
    // pointWeight(2) = 0.125939180544827 / 2.;
    // pointWeight(3) = 0.125939180544827 / 2.;
    // pointWeight(4) = 0.132394152788506 / 2.;
    // pointWeight(5) = 0.132394152788506 / 2.;
    // pointWeight(6) = 0.132394152788506 / 2.; 

    pointWeight(0) = 0.051739766065744133555179145422 / 2.;
    pointWeight(1) = 0.008007799555564801597804123460 / 2.;
    pointWeight(2) = 0.008007799555564801597804123460 / 2.;
    pointWeight(3) = 0.008007799555564801597804123460 / 2.;
    pointWeight(4) = 0.046868898981821644823226732071 / 2.;
    pointWeight(5) = 0.046868898981821644823226732071 / 2.;
    pointWeight(6) = 0.046868898981821644823226732071 / 2.; 
    pointWeight(7) = 0.046590940183976487960361770070 / 2.; 
    pointWeight(8) = 0.046590940183976487960361770070 / 2.; 
    pointWeight(9) = 0.046590940183976487960361770070 / 2.; 
    pointWeight(10) = 0.031016943313796381407646220131 / 2.; 
    pointWeight(11) = 0.031016943313796381407646220131 / 2.;
    pointWeight(12) = 0.031016943313796381407646220131 / 2.;
    pointWeight(13) = 0.010791612736631273623178240136 / 2.;
    pointWeight(14) = 0.010791612736631273623178240136 / 2.;
    pointWeight(15) = 0.010791612736631273623178240136 / 2.;
    pointWeight(16) = 0.032195534242431618819414482205 / 2.;
    pointWeight(17) = 0.032195534242431618819414482205 / 2.;
    pointWeight(18) = 0.032195534242431618819414482205 / 2.;
    pointWeight(19) = 0.015445834210701583817692900053 / 2.;
    pointWeight(20) = 0.015445834210701583817692900053 / 2.;
    pointWeight(21) = 0.015445834210701583817692900053 / 2.;
    pointWeight(22) = 0.015445834210701583817692900053 / 2.;
    pointWeight(23) = 0.015445834210701583817692900053 / 2.;
    pointWeight(24) = 0.015445834210701583817692900053 / 2.;
    pointWeight(25) = 0.017822989923178661888748319485 / 2.;
    pointWeight(26) = 0.017822989923178661888748319485 / 2.;
    pointWeight(27) = 0.017822989923178661888748319485 / 2.;
    pointWeight(28) = 0.017822989923178661888748319485 / 2.;
    pointWeight(29) = 0.017822989923178661888748319485 / 2.;
    pointWeight(30) = 0.017822989923178661888748319485 / 2.;
    pointWeight(31) = 0.037038683681384627918546472190 / 2.;
    pointWeight(32) = 0.037038683681384627918546472190 / 2.;
    pointWeight(33) = 0.037038683681384627918546472190 / 2.;
    pointWeight(34) = 0.037038683681384627918546472190 / 2.;
    pointWeight(35) = 0.037038683681384627918546472190 / 2.;
    pointWeight(36) = 0.037038683681384627918546472190 / 2.;

    return pointWeight(i);
};

template<>
double IntegQuadratureSpecial<3>::WeightList(int i){
    
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
double IntegQuadratureSpecial<2>::interpolateQuadraticVariable(
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
double IntegQuadratureSpecial<3>::interpolateQuadraticVariable(
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
double IntegQuadratureSpecial<2>::interpolateLinearVariable(
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
