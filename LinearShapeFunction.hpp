//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------LINEAR SHAPE FUNCTION-----------------------------
//------------------------------------------------------------------------------

#ifndef LINSHAPEFUNCTION_H
#define LINSHAPEFUNCTION_H

using namespace boost::numeric;

template<int DIM>
class LinShapeFunction {
public:
    
    //Defines type "Coords" to allocate the integration points coordinates
    typedef ublas::bounded_vector<double, DIM>     Coords;
    
    //Number of interpolation notes
    static const int numIntpNodes = DIM+1;
    
    //Defines the type "Values" which stores the shape function values
    typedef ublas::bounded_vector<double, numIntpNodes>           Values;

    //Defines the type "ValuesDeriv" which stores the shape function derivatives
    typedef ublas::bounded_matrix<double, DIM, numIntpNodes>      ValueDeriv;

public:
    //Evaluates the shape function value
    void evaluate(const Coords& xi, Values& phi) const;

    //Evaluates the values of the shape funtion derivatives    
    void evaluateGradient(const Coords& xi, ValueDeriv& dphi) const;   

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines linear shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void LinShapeFunction<2>::evaluate(const Coords& xi, Values& phi) const {    

    const double xsi1 = xi(0);
    const double xsi2 = xi(1);
    const double xsi3 = 1. - xsi1 - xsi2;
    
    phi(0) = xsi3;
    phi(1) = xsi1;
    phi(2) = xsi2;
    
    // element conectivity
    //     2
    //     01
    
    return;
}

template<>
void LinShapeFunction<3>::evaluate(const Coords& xi, Values& phi) const {

    const double xsi1 = xi(0);
    const double xsi2 = xi(1);
    const double xsi3 = xi(2);
    const double xsi4 = 1. - xsi1 - xsi2 - xsi3;
    
    phi(0) = xsi1;
    phi(1) = xsi2;
    phi(2) = xsi3;
    phi(3) = xsi4;
    
    // element conectivity
    //layer 1
    //     2
    //     1 3
    //layer 2
    //     0
    
    return;
}

//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void LinShapeFunction<2>::evaluateGradient(const Coords& xi, \
                                           ValueDeriv& dphi) const {
    // const double xsi1 = xi(0);
    // const double xsi2 = xi(1);
    // const double xsi3 = 1. - xsi1 - xsi2;
    
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;
    
    // element conectivity
    //     2
    //     01
     
    return;
}

template<>
void LinShapeFunction<3>::evaluateGradient(const Coords& xi, \
                                           ValueDeriv& dphi) const {

    // const double xsi1 = xi(0);
    // const double xsi2 = xi(1);
    // const double xsi3 = xi(2);
    // const double xsi4 = 1. - xsi1 - xsi2 - xsi3;
    
    dphi(0,0) = 1.;
    dphi(1,0) = 0.;
    dphi(2,0) = 0.;
    
    dphi(0,1) = 0.;
    dphi(1,1) = 1.;
    dphi(2,1) = 0.;
    
    dphi(0,2) = 0.;
    dphi(1,2) = 0.;
    dphi(2,2) = 1.;
    
    dphi(0,3) = -1.;
    dphi(1,3) = -1.;
    dphi(2,3) = -1.;
    
        // element conectivity
    //layer 1
    //     2
    //     1 3
    //layer 2
    //     0
    
    return;
}

#endif
