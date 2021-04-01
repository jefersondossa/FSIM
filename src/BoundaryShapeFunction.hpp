//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2020 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------BOUNDARY SHAPE FUNCTIONS----------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_SHAPEFUNCTION_H
#define BOUND_SHAPEFUNCTION_H

/// Defines the fluid boundary shape functions
template<int DIM, int DEG>
class BoundShapeFunction {
public:
    /// Conputes the boundary shape function value
    /// @param double Adimensional coordinates 
    /// @param double* boundary shape function
    /// @param double* boundary shape function derivatives
    void getShapeFunction(double* Xsi, double *phi_, double **dphi_);
};

//------------------------------------------------------------------------------
//--------------COMPUTE THE SHAPE FUNCTION VALUES AND DERIVATIVES---------------
//------------------------------------------------------------------------------
template<>
void BoundShapeFunction<2,1>::getShapeFunction(double* Xsi, double *phi_, double **dphi_){
        
    double aux;
    double Nnos = 2;
    double AdimCoord[2];
    AdimCoord[0] = -1.;
    AdimCoord[1] =  1.;

    
    for (int j=0; j<Nnos; j++) {

        phi_[j] = 1.0;
        dphi_[j][0] = 0.0;

        for (int i = 0; i < Nnos; i++) {

            aux = 1.0;

            if(i != j){
                
                phi_[j] = phi_[j] * (Xsi[0] - AdimCoord[i]) / 
                    (AdimCoord[j] - AdimCoord[i]);
                
                for (int k=0; k<Nnos; k++) {
                    if ((i != k) && (j != k)) aux = aux*(Xsi[0] - AdimCoord[k]);
                };
                dphi_[j][0] += aux;
            };
        };
    };

    for (int i=0; i<Nnos; i++) {
        for (int k=0; k<Nnos; k++) {
            if (i != k) dphi_[i][0] = dphi_[i][0]/(AdimCoord[i]-AdimCoord[k]);
        };
    };    
    
    return;
};

template<>
void BoundShapeFunction<2,2>::getShapeFunction(double* Xsi, double *phi_, double **dphi_){
        
    double aux;
    double Nnos = 3;
    double AdimCoord[3];
    AdimCoord[0] = -1.;
    AdimCoord[1] =  0.;
    AdimCoord[2] =  1.;

    
    for (int j=0; j<Nnos; j++) {

        phi_[j] = 1.0;
        dphi_[j][0] = 0.0;

        for (int i = 0; i < Nnos; i++) {

            aux = 1.0;

            if(i != j){
                
                phi_[j] = phi_[j] * (Xsi[0] - AdimCoord[i]) / 
                    (AdimCoord[j] - AdimCoord[i]);
                
                for (int k=0; k<Nnos; k++) {
                    if ((i != k) && (j != k)) aux = aux*(Xsi[0] - AdimCoord[k]);
                };
                dphi_[j][0] += aux;
            };
        };
    };

    for (int i=0; i<Nnos; i++) {
        for (int k=0; k<Nnos; k++) {
            if (i != k) dphi_[i][0] = dphi_[i][0]/(AdimCoord[i]-AdimCoord[k]);
        };
    };    
    
    return;
};

template<>
void BoundShapeFunction<3,1>::getShapeFunction(double* Xsi, double *phi_, double **dphi_){
    
    double xsi1 = Xsi[0];
    double xsi2 = Xsi[1];
    double xsi3 = 1. - xsi1 - xsi2;
    
    phi_[0] = xsi3;
    phi_[1] = xsi1;
    phi_[2] = xsi2;

    dphi_[0][0] = -1.;
    dphi_[0][1] = -1.;
    
    dphi_[1][0] = 1.;
    dphi_[1][1] = 0.;
    
    dphi_[2][0] = 0.;
    dphi_[2][1] = 1.;

    return;
}

template<>
void BoundShapeFunction<3,2>::getShapeFunction(double* Xsi, double *phi_, double **dphi_){

    double xsi1 = Xsi[0];
    double xsi2 = Xsi[1];
    double xsi3 = 1. - xsi1 - xsi2;

    phi_[0] = xsi3 * (2.0 * xsi3 - 1.0);
    phi_[1] = xsi1 * (2.0 * xsi1 - 1.0);
    phi_[2] = xsi2 * (2.0 * xsi2 - 1.0);
    phi_[3] = 4.0 * xsi3 * xsi1;
    phi_[4] = 4.0 * xsi1 * xsi2;
    phi_[5] = 4.0 * xsi2 * xsi3;


    dphi_[1][0] = 4. * xsi1 - 1.;
    dphi_[1][1] = 0.;

    dphi_[2][0] = 0.;
    dphi_[2][1] = 4. * xsi2 - 1.;

    dphi_[0][0] = -4. * xsi3 + 1.;
    dphi_[0][1] = -4. * xsi3 + 1.;

    dphi_[4][0] = 4. * xsi2;
    dphi_[4][1] = 4. * xsi1;

    dphi_[5][0] = -4. * xsi2;
    dphi_[5][1] = 4. * (xsi3 - xsi2);

    dphi_[3][0] = 4. * (xsi3 - xsi1);
    dphi_[3][1] = -4. * xsi1;

    return;
}

#endif
