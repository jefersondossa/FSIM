//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------------FLUID BOUNDARY--------------------------------
//------------------------------------------------------------------------------

#ifndef BOUNDARY_H
#define BOUNDARY_H

template<int DIM>
class Boundary{

public: 
    typedef ublas::bounded_vector<int, 3*(DIM-1)>         BoundConnect;

    //Defines the boundary element
    Boundary(BoundConnect& connec, int index,  \
             int *constrain, double *values, int gr){
        connectB_ = connec;
        index_ = index;
        group_ = gr;

        constrainType[0] = constrain[0];
        constrainType[1] = constrain[1];
        constrainType[2] = constrain[2];
        
        constrainValue[0] = values[0];
        constrainValue[1] = values[1];
        constrainValue[2] = values[2];

        element_ = 0;
        elementSide_ = 0;
    };     

    //Returns the constrain type
    int getConstrain(int dir){return constrainType[dir];}
    //Returns the constrain value
    double getConstrainValue(int dir){return constrainValue[dir];}
    //Returns the boundary element connectivity
    BoundConnect getBoundaryConnectivity(){return connectB_;}

    void setBoundaryGroup(int gr){group_ = gr;}
    int getBoundaryGroup(){return group_;}

    void setElement(int el){element_ = el;}
    void setElementSide(int el){elementSide_ = el;}

    int getElement(){return element_;}
    int getElementSide(){return elementSide_;}

private:
    BoundConnect connectB_;         //Boundary element connectivity
    int          index_;            //Boundary element index
    int          constrainType[3];  //Element type of constrain
    double       constrainValue[3]; //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group

};






















#endif
